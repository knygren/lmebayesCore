## Step 2: posterior engine (mode, Xi, gradient, Hessian).

#' @noRd
.group_floor_prepare_group_data <- function(design,
                                            family_hook,
                                            group_levels,
                                            offset = NULL,
                                            weights = NULL,
                                            dispersion = 1) {
  .lmerb_validate_design(design)
  if (is.null(design$y) || is.null(design$D)) {
    stop("'design' must contain 'y' and 'D'.", call. = FALSE)
  }

  group_info <- normalize_group(design$group, nrow(design$D))
  idx_map <- match(group_levels, as.character(group_info$ids))
  if (anyNA(idx_map)) {
    stop("Some group_levels are not present in design$group.", call. = FALSE)
  }
  row_idx <- group_info$rows[idx_map]
  J <- length(group_levels)
  p_re <- length(design$groupef.names)

  n_obs <- vapply(row_idx, length, integer(1L))
  l2 <- nrow(design$D)

  if (is.null(offset)) {
    offset <- if (!is.null(design$offset)) design$offset else rep(0, l2)
  }
  offset <- as.numeric(offset)
  if (length(offset) == 1L) offset <- rep(offset, l2)

  if (is.null(weights)) {
    weights <- if (!is.null(design$weights)) design$weights else rep(1, l2)
  }
  weights <- as.numeric(weights)
  if (length(weights) == 1L) weights <- rep(weights, l2)

  if (identical(family_hook$name, "binomial_logit")) {
    size <- weights
  } else {
    size <- rep(1, l2)
  }

  Z_list <- vector("list", J)
  y_list <- vector("list", J)
  eta_offset <- vector("list", J)
  size_list <- vector("list", J)
  wt_list <- vector("list", J)

  for (j in seq_len(J)) {
    rows <- row_idx[[j]]
    Z_list[[j]] <- design$D[rows, , drop = FALSE]
    y_list[[j]] <- design$y[rows]
    eta_offset[[j]] <- offset[rows]
    size_list[[j]] <- size[rows]
    wt_list[[j]] <- if (identical(family_hook$name, "binomial_logit")) {
      rep(1, length(rows))
    } else {
      weights[rows]
    }
  }

  list(
    group_levels = group_levels,
    row_idx = row_idx,
    Z_list = Z_list,
    y_list = y_list,
    eta_offset = eta_offset,
    size_list = size_list,
    wt_list = wt_list,
    n_obs = n_obs,
    dispersion = dispersion,
    J = J,
    p_re = p_re
  )
}

#' @noRd
.group_floor_engine <- function(prior_int,
                              group_data,
                              family_hook,
                              b_dag = NULL) {
  Lambda_b <- prior_int$Lambda_beta
  mu_b <- prior_int$mu_beta
  J <- group_data$J
  p_re <- group_data$p_re
  n_dim <- prior_int$n_dim

  eta_j <- function(j, bj) {
    drop(group_data$Z_list[[j]] %*% bj + group_data$eta_offset[[j]])
  }

  negll_j <- function(j, bj) {
    family_hook$negll_j(
      y = group_data$y_list[[j]],
      eta = eta_j(j, bj),
      size = group_data$size_list[[j]],
      wt = group_data$wt_list[[j]],
      dispersion = group_data$dispersion
    )
  }

  f <- function(b) {
    ll <- sum(vapply(seq_len(J), function(j) negll_j(j, b[.group_floor_idx(j, p_re)]),
                     numeric(1L)))
    d <- b - mu_b
    0.5 * drop(crossprod(d, Lambda_b %*% d)) + ll
  }

  gr_f <- function(b) {
    g <- as.vector(Lambda_b %*% (b - mu_b))
    for (j in seq_len(J)) {
      o <- .group_floor_idx(j, p_re)
      bj <- b[o]
      Z <- group_data$Z_list[[j]]
      e <- eta_j(j, bj)
      y <- group_data$y_list[[j]]
      sz <- group_data$size_list[[j]]
      wt <- group_data$wt_list[[j]]
      if (identical(family_hook$name, "binomial_logit")) {
        g[o] <- g[o] - drop(crossprod(Z, y - sz * plogis(e)))
      } else if (identical(family_hook$name, "poisson_log")) {
        g[o] <- g[o] - drop(crossprod(Z, y - sz * wt * exp(e)))
      } else if (identical(family_hook$name, "gaussian")) {
        g[o] <- g[o] - drop(crossprod(Z, wt * (y - e) / group_data$dispersion))
      }
    }
    g
  }

  Gmat <- function(j, bj) {
    family_hook$Gmat_j(
      Z = group_data$Z_list[[j]],
      eta = eta_j(j, bj),
      size = group_data$size_list[[j]],
      wt = group_data$wt_list[[j]],
      dispersion = group_data$dispersion
    )
  }

  He_f <- function(b) {
    H <- Lambda_b
    for (j in seq_len(J)) {
      o <- .group_floor_idx(j, p_re)
      H[o, o] <- H[o, o] + Gmat(j, b[o])
    }
    0.5 * (H + t(H))
  }

  newton <- function(b0, rhs = function(x) numeric(n_dim)) {
    b <- b0
    for (k in seq_len(300L)) {
      s <- solve(He_f(b), gr_f(b) - rhs(b))
      st <- 1
      while (st > 1e-10 && max(abs(s * st)) > 5) st <- st / 2
      b <- b - st * s
      if (max(abs(s)) < 1e-10) break
    }
    b
  }

  if (is.null(b_dag)) {
    b_dag <- newton(mu_b)
  }
  f_dag <- f(b_dag)
  Xi <- function(b) f(b) - f_dag

  emb <- function(j, u) {
    v <- numeric(n_dim)
    v[.group_floor_idx(j, p_re)] <- u
    v
  }

  list(
    Lambda_b = Lambda_b,
    mu_b = mu_b,
    f = f,
    gr_f = gr_f,
    He_f = He_f,
    Gmat = Gmat,
    Xi = Xi,
    b_dag = b_dag,
    f_dag = f_dag,
    newton = newton,
    emb = emb,
    eta_j = eta_j,
    negll_j = negll_j,
    n_dim = n_dim,
    J = J,
    p_re = p_re
  )
}

#' @noRd
.group_floor_mode_from_icm <- function(design,
                                       family,
                                       measurement_prior_list,
                                       group_levels,
                                       tol = 1e-10,
                                       maxit = 200L) {
  pm <- glmerb_posterior_mode(
    design = design,
    family = family,
    measurement_prior_list = measurement_prior_list,
    tol = tol,
    maxit = maxit
  )
  b_mat <- pm$b_mean
  rn <- rownames(b_mat)
  if (!is.null(rn)) {
    b_mat <- b_mat[group_levels, , drop = FALSE]
  }
  list(
    beta = .group_floor_stack_beta(b_mat),
    icm = pm
  )
}
