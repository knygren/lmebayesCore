## Chapter C05 restricted Gibbs: EM fixed point gamma* and Jacobian.

#' Conditional-mean E-step for the population mean map.
#' @noRd
.c05_estep <- function(design,
                                    fixef,
                                    p11,
                                    measurement_prior_list,
                                    family,
                                    estep = c("exact", "aghq", "mc"),
                                    n = 10000L,
                                    mc_seed = NULL) {
  estep <- match.arg(estep)
  if (identical(estep, "aghq")) {
    stop("estep = \"aghq\" is not implemented yet; use \"mc\".", call. = FALSE)
  }
  if (!identical(family$family, "gaussian") && identical(estep, "exact")) {
    stop(
      "estep = \"exact\" is only available for gaussian(); ",
      "use \"mc\" or \"aghq\".",
      call. = FALSE
    )
  }

  group_levels <- p11$group_levels
  J <- p11$J
  p_re <- p11$p_re
  re_names <- p11$re_names
  g_chr <- as.character(design$group)
  gamma <- .c05_gamma_from_fixef(fixef, p11)
  sigma2 <- measurement_prior_list$group.dispersion
  P_b <- p11$P_b
  Sigma_b <- measurement_prior_list$group.Sigma
  is_gaussian <- identical(family$family, "gaussian")

  b_mean <- matrix(
    0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )
  b_mc_se <- NULL
  V_list <- stats::setNames(vector("list", J), group_levels)

  if (is_gaussian && identical(estep, "exact")) {
    if (is.null(sigma2)) {
      stop(
        "'measurement_prior_list$group.dispersion' is required for gaussian().",
        call. = FALSE
      )
    }
    sigma2 <- as.numeric(sigma2)
    if (!(length(sigma2) %in% c(1L, J))) {
      stop(
        "'measurement_prior_list$group.dispersion' must have length 1 or J.",
        call. = FALSE
      )
    }

    if (!is.null(p11$lmerb_system)) {
      b_out <- .lmerb_posterior_b_given_gamma(p11$lmerb_system, design, fixef)
      b_mean <- b_out$b
      for (lev in group_levels) {
        V_list[[lev]] <- solve(p11$lmerb_system$post_P_j_list[[lev]])
      }
    } else {
      for (jj in seq_len(J)) {
        lev <- group_levels[jj]
        rows <- which(g_chr == lev)
        Z_j <- design$D[rows, , drop = FALSE]
        y_j <- design$y[rows]
        sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
        H_j <- p11$H_list[[lev]]
        mu_j <- as.vector(H_j %*% gamma)

        post_P_j <- crossprod(Z_j) / sigma2_j + P_b
        post_v_j <- crossprod(Z_j, y_j) / sigma2_j + P_b %*% mu_j
        b_mean[jj, ] <- solve(post_P_j, post_v_j)
        V_list[[lev]] <- solve(post_P_j)
      }
    }
  } else if (identical(estep, "mc")) {
    if (!is.null(mc_seed)) set.seed(mc_seed)
    mu_all <- as.matrix(build_mu_all(design, fixef, group_levels = group_levels)$mu_all)
    if (is_gaussian) {
      if (is.null(sigma2)) {
        stop(
          "'measurement_prior_list$group.dispersion' is required for gaussian().",
          call. = FALSE
        )
      }
      sigma2 <- as.numeric(sigma2)
    }

    for (jj in seq_len(J)) {
      lev <- group_levels[jj]
      rows <- which(g_chr == lev)
      y_j <- design$y[rows]
      Z_j <- design$D[rows, , drop = FALSE]
      mu_j <- mu_all[, jj]

      pf_j <- if (is_gaussian) {
        sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
        glmbayesCore::dNormal(
          mu = mu_j, Sigma = Sigma_b, dispersion = sigma2_j
        )
      } else {
        glmbayesCore::dNormal(mu = mu_j, Sigma = Sigma_b)
      }
      draws <- matrix(0, n, p_re)
      for (m in seq_len(n)) {
        fit_j <- glmbayesCore::rglmb(
          n = 1L,
          y = y_j,
          x = Z_j,
          family = family,
          pfamily = pf_j,
          verbose = FALSE
        )
        draws[m, ] <- .minorization_rglmb_draw(fit_j)
      }
      b_mean[jj, ] <- colMeans(draws)
      V_list[[lev]] <- stats::cov(draws)
      if (is.null(b_mc_se)) {
        b_mc_se <- matrix(
          0, nrow = J, ncol = p_re,
          dimnames = list(group_levels, re_names)
        )
      }
      b_mc_se[jj, ] <- apply(draws, 2L, stats::sd) / sqrt(n)
    }
  }

  list(b_mean = b_mean, b_mc_se = b_mc_se, V_list = V_list, estep = estep, n = n)
}

#' One C05 mean-map update M(gamma) given conditional means b_j.
#' @noRd
.c05_mean_map <- function(b_mean, p11) {
  if (!is.null(p11$lmerb_system) && is.null(b_mean)) {
    return(.c05_mean_map_lmerb(b_mean, p11))
  }

  rhs <- p11$Lambda_gamma %*% p11$mu_0
  for (j in seq_len(p11$J)) {
    b_j <- b_mean[j, ]
    H_j <- p11$H_list[[p11$group_levels[[j]]]]
    rhs <- rhs + t(H_j) %*% p11$P_b %*% b_j
  }
  gamma <- solve(p11$P11, rhs)
  .c05_fixef_from_gamma(gamma, p11)
}

#' Block~2 update using the exact Gaussian Schur system (\code{lmerb}).
#' @noRd
.c05_mean_map_lmerb <- function(b_mean, p11) {
  system <- p11$lmerb_system
  gamma <- solve(system$M, system$v)
  .c05_fixef_from_gamma(gamma, p11)
}

#' Conditional-mode Block~1 update for ICM initialization.
#' @noRd
.c05_block1_modes <- function(design,
                              fixef,
                              p11,
                              measurement_prior_list,
                              family) {
  group_levels <- p11$group_levels
  J <- p11$J
  p_re <- p11$p_re
  re_names <- p11$re_names
  g_chr <- as.character(design$group)
  gamma <- .c05_gamma_from_fixef(fixef, p11)
  sigma2 <- measurement_prior_list$group.dispersion
  P_b <- p11$P_b
  Sigma_b <- measurement_prior_list$group.Sigma
  is_gaussian <- identical(family$family, "gaussian")

  b_mode <- matrix(
    0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )

  if (is_gaussian) {
    sigma2 <- as.numeric(sigma2)
    for (jj in seq_len(J)) {
      lev <- group_levels[jj]
      rows <- which(g_chr == lev)
      Z_j <- design$D[rows, , drop = FALSE]
      y_j <- design$y[rows]
      sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
      H_j <- p11$H_list[[lev]]
      mu_j <- as.vector(H_j %*% gamma)
      post_P_j <- crossprod(Z_j) / sigma2_j + P_b
      post_v_j <- crossprod(Z_j, y_j) / sigma2_j + P_b %*% mu_j
      b_mode[jj, ] <- solve(post_P_j, post_v_j)
    }
  } else {
    mu_all <- as.matrix(
      build_mu_all(design, fixef, group_levels = group_levels)$mu_all
    )
    for (jj in seq_len(J)) {
      lev <- group_levels[jj]
      rows <- which(g_chr == lev)
      y_j <- design$y[rows]
      Z_j <- design$D[rows, , drop = FALSE]
      mu_j <- mu_all[, jj]
      pf_j <- glmbayesCore::dNormal(mu = mu_j, Sigma = Sigma_b)
      fit_j <- glmbayesCore::rglmb(
        n = 1L,
        y = y_j,
        x = Z_j,
        family = family,
        pfamily = pf_j,
        verbose = FALSE
      )
      b_mode[jj, ] <- .minorization_rglmb_mode(fit_j)
    }
  }

  b_mode
}

#' ICM start for \code{population_mode()} (conditional modes + C05 mean map).
#' @noRd
.c05_icm_init <- function(design,
                          fixef_start,
                          p11,
                          measurement_prior_list,
                          family,
                          tol = 1e-8,
                          maxit = 200L) {
  if (identical(family$family, "gaussian") &&
      !is.null(p11$lmerb_system)) {
    pm <- lmerb_posterior_mean(design, measurement_prior_list, tol, maxit)
    return(list(
      fixef = pm$fixef,
      b_mode = pm$b_mean,
      converged = pm$converged,
      iterations = pm$iterations,
      delta = pm$delta
    ))
  }

  fixef <- fixef_start
  converged <- FALSE
  delta <- NA_real_

  for (iter in seq_len(maxit)) {
    b_mode <- .c05_block1_modes(
      design = design,
      fixef = fixef,
      p11 = p11,
      measurement_prior_list = measurement_prior_list,
      family = family
    )
    fixef_new <- .c05_mean_map(b_mode, p11)
    delta <- .c05_gamma_delta(fixef, fixef_new, p11)
    fixef <- fixef_new
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "ICM initialization did not converge in ", maxit,
      " iterations (final delta = ", signif(delta, 3L), ").",
      call. = FALSE
    )
  }

  list(
    fixef = fixef,
    b_mode = b_mode,
    converged = converged,
    iterations = iter,
    delta = delta
  )
}

#' Mahalanobis change in \eqn{\gamma} under \eqn{P_{11}}.
#' @noRd
.c05_gamma_delta <- function(fixef_old, fixef_new, p11) {
  d_gamma <- .c05_gamma_from_fixef(fixef_new, p11) -
    .c05_gamma_from_fixef(fixef_old, p11)
  sqrt(as.numeric(crossprod(d_gamma, p11$P11 %*% d_gamma)))
}

#' MC noise floor for \eqn{\|\Delta\gamma\|_{P_{11}}} from \code{b_mc_se}.
#'
#' Propagates independent elementwise MC standard errors in \code{b_mean}
#' through the C05 mean map and returns one MC standard deviation of the
#' resulting \eqn{P_{11}}-norm update (i.e. stop when \code{delta} is not
#' statistically distinguishable from MC noise at the current \code{n}).
#' @noRd
.c05_mc_delta_floor <- function(b_mc_se, p11) {
  if (is.null(b_mc_se)) {
    return(0)
  }
  se2 <- as.numeric(b_mc_se)
  if (!any(is.finite(se2) & se2 > 0)) {
    return(0)
  }

  P11_inv <- chol2inv(p11$chol_P11)
  Cov_gamma <- matrix(0, p11$q, p11$q)
  for (j in seq_len(p11$J)) {
    lev <- p11$group_levels[[j]]
    H_j <- p11$H_list[[lev]]
    A_j <- P11_inv %*% t(H_j) %*% p11$P_b
    se2_j <- b_mc_se[j, ]^2
    Cov_gamma <- Cov_gamma + A_j %*% diag(se2_j, p11$p_re) %*% t(A_j)
  }

  sq <- sum(diag(p11$P11 %*% Cov_gamma))
  sqrt(max(sq, 0))
}

#' Effective EM tolerance when the E-step is Monte Carlo.
#' @noRd
.c05_em_tol <- function(tol, estep, b_mc_se, p11) {
  if (!identical(estep, "mc")) {
    return(tol)
  }
  max(tol, .c05_mc_delta_floor(b_mc_se, p11))
}

#' Jacobian J, spectrum, and closure objects at the current E-step.
#' @noRd
.c05_jacobian <- function(p11, estep_out) {
  q <- p11$q
  P11_inv <- chol2inv(p11$chol_P11)
  J_mat <- matrix(0, q, q)

  for (j in seq_len(p11$J)) {
    lev <- p11$group_levels[[j]]
    H_j <- p11$H_list[[lev]]
    V_j <- estep_out$V_list[[lev]]
    J_mat <- J_mat + P11_inv %*% t(H_j) %*% p11$P_b %*% V_j %*% p11$P_b %*% H_j
  }

  spec <- .c05_coupling_spectrum(p11, J_mat, strict = FALSE)
  kappa <- spec$kappa_max
  rho <- spec$rho

  Iq <- diag(q)
  Sigma_pi <- tryCatch(
    solve(p11$P11 %*% (Iq - J_mat)),
    error = function(e) NULL
  )

  list(
    tilde_J = J_mat,
    kappa_spectrum = spec$kappa,
    weights = spec$weights,
    kappa = kappa,
    rho = rho,
    Sigma_pi = Sigma_pi
  )
}

#' EM fixed point for the C05 population mode gamma*.
#'
#' Finds the fixed point of the conditional-mean map (Chapter C05 Stage 1).
#' Uses conditional means for Block~1, not modes.
#'
#' @param design A \code{\link{model_setup}} list.
#' @param pfamily_list Block~2 prior list from \code{\link{pfamily_list}()}.
#' @param family A \code{\link[stats]{family}} object.
#' @param dispprior_list Optional Block~1 dispersion prior for \code{gaussian()}.
#' @param estep E-step tier: \code{"exact"} (default Gaussian closed form),
#'   \code{"mc"} (simulated conditional means; valid for Gaussian and
#'   non-Gaussian), or \code{"aghq"} (not yet implemented).
#' @param acceleration \code{"none"} or \code{"squarem"} (not yet implemented).
#' @param n Number of \eqn{\beta_j} draws per group when \code{estep = "mc"}.
#'   The sample mean \eqn{\bar b_j} has elementwise MC standard error
#'   \eqn{\mathrm{sd}(\beta_j \mid \gamma, y)/\sqrt{n}}; relative error on
#'   means is therefore \eqn{O(1/\sqrt{n})} (e.g. \code{n = 10000} targets
#'   \eqn{\approx 1\%} of the posterior scale).
#' @param mc_seed Optional seed for the MC E-step (first iteration only).
#' @param icm_init If \code{TRUE}, run iterated conditional modes (Block~1
#'   modes via \code{\link[glmbayesCore]{rglmb}}, Block~2 via the C05 mean map)
#'   before EM. Gaussian Block~1 uses exact conditional modes.
#' @param icm_tol ICM convergence tolerance on \eqn{\|\Delta\gamma\|_{P_{11}}}.
#' @param icm_maxit Maximum ICM iterations.
#' @param tol Convergence tolerance on the Mahalanobis change in \code{fixef}
#'   under the C05 metric \eqn{\|\Delta\gamma\|_{P_{11}}}. When
#'   \code{estep = "mc"}, the effective tolerance is
#'   \code{max(tol, mc_delta_floor)}, where \code{mc_delta_floor} propagates
#'   the terminal \code{b_mc_se} (MC standard errors at draw count \code{n})
#'   through the mean map: EM stops once successive iterates differ by less
#'   than one MC standard deviation in the \eqn{P_{11}} metric.
#' @param maxit Maximum EM iterations.
#' @return A list with \code{fixef} (\eqn{\gamma^\star}), \code{gamma_star},
#'   \code{b_mean}, \code{b_mc_se} (terminal MC standard errors when
#'   \code{estep = "mc"}), \code{V_list}, \code{n}, model context
#'   (\code{design}, \code{family}, \code{measurement_prior_list}, \code{p11}),
#'   refresh objects (\code{P11}, \code{Sigma_star}), \code{tilde_J},
#'   \code{kappa_spectrum}, \code{weights}, \code{kappa}, \code{rho},
#'   \code{Sigma_pi}, \code{eps_star_closure},
#'   \code{icm} initialization diagnostics, EM diagnostics (\code{tol_eff},
#'   \code{mc_delta_floor} when \code{estep = "mc"}), and related fields.
#' @export
population_mode <- function(design,
                                             pfamily_list,
                                             family = gaussian(),
                                             dispprior_list = NULL,
                                             estep = c("exact", "aghq", "mc"),
                                             acceleration = c("none", "squarem"),
                                             n = 10000L,
                                             mc_seed = NULL,
                                             icm_init = TRUE,
                                             icm_tol = 1e-8,
                                             icm_maxit = 200L,
                                             tol = 1e-10,
                                             maxit = 200L) {
  estep <- match.arg(estep)
  acceleration <- match.arg(acceleration)
  if (!identical(acceleration, "none")) {
    stop("acceleration = \"squarem\" is not implemented yet.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!(n >= 1L)) {
    stop("'n' must be at least 1.", call. = FALSE)
  }

  prep <- .c05_validate(
    design, pfamily_list, family, dispprior_list
  )
  mpl <- prep$measurement_prior_list
  p11 <- .c05_p11(
    design,
    prep$prior_pack,
    prep$group_levels,
    measurement_prior_list = mpl,
    family = family
  )

  fixef <- lapply(mpl$pop.prior_list, `[[`, "mu")
  names(fixef) <- prep$re_names

  icm <- NULL
  if (isTRUE(icm_init)) {
    icm <- .c05_icm_init(
      design = design,
      fixef_start = fixef,
      p11 = p11,
      measurement_prior_list = mpl,
      family = family,
      tol = icm_tol,
      maxit = icm_maxit
    )
    fixef <- icm$fixef
  }

  converged <- FALSE
  delta <- NA_real_
  em_iterations <- 0L
  tol_eff <- tol
  mc_delta_floor <- NA_real_

  if (identical(family$family, "gaussian") && identical(estep, "exact") &&
      !is.null(p11$lmerb_system)) {
    fixef <- .c05_mean_map_lmerb(NULL, p11)
    estep_out <- .c05_estep(
      design = design,
      fixef = fixef,
      p11 = p11,
      measurement_prior_list = mpl,
      family = family,
      estep = estep,
      n = n
    )
    em_iterations <- 1L
    converged <- TRUE
    delta <- .c05_gamma_delta(fixef, fixef, p11)
  } else {
    for (iter in seq_len(maxit)) {
      em_iterations <- iter
      estep_out <- .c05_estep(
        design = design,
        fixef = fixef,
        p11 = p11,
        measurement_prior_list = mpl,
        family = family,
        estep = estep,
        n = n,
        mc_seed = if (iter == 1L) mc_seed else NULL
      )

      fixef_new <- .c05_mean_map(estep_out$b_mean, p11)

      delta <- .c05_gamma_delta(fixef, fixef_new, p11)
      mc_delta_floor <- if (identical(estep, "mc")) {
        .c05_mc_delta_floor(estep_out$b_mc_se, p11)
      } else {
        NA_real_
      }
      tol_eff <- .c05_em_tol(tol, estep, estep_out$b_mc_se, p11)

      fixef <- fixef_new
      if (delta < tol_eff) {
        converged <- TRUE
        break
      }
    }
  }

  if (!converged) {
    msg <- paste0(
      "population_mode() did not converge in ", maxit,
      " iterations (final delta = ", signif(delta, 3L)
    )
    if (identical(estep, "mc") && is.finite(mc_delta_floor)) {
      msg <- paste0(
        msg,
        ", effective tol = ", signif(tol_eff, 3L),
        " from MC n = ", n,
        ", mc_delta_floor = ", signif(mc_delta_floor, 3L), ")"
      )
    } else {
      msg <- paste0(msg, ").")
    }
    warning(msg, call. = FALSE)
  }

  jac <- .c05_jacobian(p11, estep_out)
  eps_star_closure <- .c05_epsilon_closure(jac$tilde_J)

  stationarity <- {
    fixef_check <- .c05_mean_map(estep_out$b_mean, p11)
    .c05_gamma_delta(fixef, fixef_check, p11)
  }

  list(
    fixef = fixef,
    gamma_star = .c05_gamma_from_fixef(fixef, p11),
    b_mean = estep_out$b_mean,
    b_mc_se = estep_out$b_mc_se,
    V_list = estep_out$V_list,
    icm = icm,
    design = design,
    family = family,
    measurement_prior_list = mpl,
    group_levels = prep$group_levels,
    p11 = p11,
    P11 = p11$P11,
    P11_RE = p11$P11_RE,
    Sigma_star = p11$Sigma_star,
    chol_P11 = p11$chol_P11,
    tilde_J = jac$tilde_J,
    kappa_spectrum = jac$kappa_spectrum,
    weights = jac$weights,
    kappa = jac$kappa,
    rho = jac$rho,
    Sigma_pi = jac$Sigma_pi,
    eps_star_closure = eps_star_closure,
    stationarity = stationarity,
    converged = converged,
    iterations = em_iterations,
    delta = delta,
    tol_eff = tol_eff,
    mc_delta_floor = mc_delta_floor,
    estep = estep_out$estep,
    n = if (identical(estep_out$estep, "mc")) n else NA_integer_,
    q = p11$q,
    call = match.call()
  )
}
