## Steps 3-4: kappa bounds, support functions, per-group precision floors.

#' Proposition 2 level: P(Gamma(d) > r) / P(Gamma(d) <= r) = epsilon.
#' @noRd
.gamma_level_for_budget <- function(d, epsilon) {
  if (!is.finite(epsilon) || epsilon <= 0 || epsilon >= 1) {
    stop("'epsilon' must be in (0, 1).", call. = FALSE)
  }
  d <- as.numeric(d)
  if (d < 1) stop("'d' must be at least 1.", call. = FALSE)
  uniroot(
    function(r) {
      pgamma(r, d, lower.tail = FALSE) / pgamma(r, d, lower.tail = TRUE) - epsilon
    },
    c(1e-8, 200 * d + 200),
    tol = 1e-10
  )$root
}

#' @noRd
.group_floor_kappa_crude <- function(j, Lambda_b, Gbar_list, p_re, J) {
  o <- .group_floor_idx(j, p_re)
  n_dim <- J * p_re
  L <- Lambda_b[-o, -o, drop = FALSE]
  idx_other <- setdiff(seq_len(n_dim), o)
  Gb <- matrix(0, length(idx_other), length(idx_other))
  for (kj in setdiff(seq_len(J), j)) {
    rg <- .group_floor_idx(kj, p_re)
    loc <- match(rg, idx_other)
    Gb[loc, loc] <- Gbar_list[[kj]]
  }
  0.5 * as.numeric(determinant(diag(nrow(Gb)) + solve(L, Gb), logarithm = TRUE)$modulus)
}

#' @noRd
.group_floor_kappa_laplace <- function(j, engine, cv, r_level) {
  o <- .group_floor_idx(j, engine$p_re)
  logdet_mj <- function(b) {
    as.numeric(determinant(engine$He_f(b)[-o, -o, drop = FALSE], logarithm = TRUE)$modulus)
  }
  sup <- .group_floor_supp(engine, cv, r_level)
  0.5 * abs(logdet_mj(sup$beta) - logdet_mj(engine$b_dag))
}

#' @noRd
.group_floor_pathb <- function(engine, t, cv, start) {
  rhs <- function(x) t * cv
  engine$newton(start, rhs)
}

#' @noRd
.group_floor_supp <- function(engine, cv, r_level, start = engine$b_dag) {
  pathb <- function(t) .group_floor_pathb(engine, t, cv, start)
  hi <- 1
  bh <- pathb(hi)
  while (engine$Xi(bh) < r_level && hi < 1e12) {
    hi <- hi * 4
    bh <- pathb(hi)
  }
  tt <- uniroot(
    function(t) engine$Xi(pathb(t)) - r_level,
    c(0, hi),
    tol = 1e-10
  )$root
  b <- pathb(tt)
  list(val = drop(crossprod(cv, b)), beta = b, t = tt)
}

#' @noRd
.group_floor_omega <- function(Gamma, P_b) {
  ePb <- eigen(P_b, symmetric = TRUE)
  Pmh <- ePb$vectors %*% diag(1 / sqrt(ePb$values)) %*% t(ePb$vectors)
  1 / (1 + min(eigen(Pmh %*% Gamma %*% Pmh, symmetric = TRUE)$values))
}

#' @noRd
.group_floor_one_group <- function(j,
                                   engine,
                                   group_data,
                                   family_hook,
                                   r_level) {
  Z <- group_data$Z_list[[j]]
  n_i <- nrow(Z)
  if (length(r_level) == 1L) {
    r_level <- rep(r_level, n_i)
  } else if (length(r_level) != n_i) {
    stop("length(r_level) must be 1 or the number of rows in group ", j, ".",
         call. = FALSE)
  }

  lo <- numeric(n_i)
  hi <- numeric(n_i)

  for (i in seq_len(n_i)) {
    cv <- engine$emb(j, Z[i, ])
    hi[i] <- .group_floor_supp(engine, cv, r_level[i])$val
    lo[i] <- -.group_floor_supp(engine, engine$emb(j, -Z[i, ]), r_level[i])$val
  }

  es <- pmax(abs(lo), abs(hi))

  wbar <- vapply(seq_len(n_i), function(i) {
    family_hook$weight_min_on_interval(
      lo = lo[i],
      hi = hi[i],
      size = group_data$size_list[[j]][i],
      wt = group_data$wt_list[[j]][i]
    )
  }, numeric(1L))

  if (isTRUE(family_hook$gaussian)) {
    wbar <- wbar / group_data$dispersion
  }

  Gamma <- crossprod(Z, wbar * Z)
  list(
    lo = lo,
    hi = hi,
    es = es,
    wbar = wbar,
    Gamma = Gamma,
    omega = .group_floor_omega(Gamma, engine$Lambda_b[.group_floor_idx(j, engine$p_re),
                                                     .group_floor_idx(j, engine$p_re)])
  )
}

#' @noRd
.group_floor_build_gbar <- function(group_data, family_hook) {
  lapply(seq_len(group_data$J), function(j) {
    Z <- group_data$Z_list[[j]]
    wmax <- family_hook$w_max(group_data$size_list[[j]], group_data$wt_list[[j]])
    if (isTRUE(family_hook$gaussian)) {
      wmax <- wmax / group_data$dispersion
    }
    crossprod(Z, wmax * Z)
  })
}

#' @noRd
.group_floor_levels <- function(epsilon,
                                union_over,
                                J,
                                p_re,
                                n_obs_total,
                                kappa,
                                gaussian = FALSE,
                                kappa_method = "crude") {
  r_joint <- .gamma_level_for_budget(J * p_re, epsilon)

  if (identical(union_over, "groups")) {
    r_base <- .gamma_level_for_budget(p_re, epsilon / J)
    R <- if (gaussian || identical(kappa_method, "none")) {
      rep(r_base, J)
    } else {
      r_base + 2 * kappa
    }
    scheme <- "groups"
  } else if (identical(union_over, "functionals")) {
    r_base <- .gamma_level_for_budget(1L, epsilon / n_obs_total)
    R <- rep(r_base, J)
    scheme <- "functionals"
  } else {
    stop("'union_over' must be 'groups' or 'functionals'.", call. = FALSE)
  }

  list(
    epsilon = epsilon,
    scheme = scheme,
    r_joint = r_joint,
    r_group = if (identical(union_over, "groups")) r_base else NA_real_,
    r_functional = if (identical(union_over, "functionals")) r_base else NA_real_,
    kappa = kappa,
    R = R,
    s = sqrt(2 * R)
  )
}
