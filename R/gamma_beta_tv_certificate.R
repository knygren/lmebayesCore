## Sharpest displayed gamma-beta TV certificate (Rosenthal route).

#' @noRd
.c05_rosenthal_expanded_box <- function(drift, eps, k, alpha, display_mode) {
  km <- drift$kappa_max_lb
  q <- drift$q
  kap_sum <- sum(drift$kappa_lb)
  b_inner <- 1 - km^2 + q / 2 + 0.5 * kap_sum
  drift_num <- 1 + 2 * b_inner + km^2 * drift$V_gamma_0
  drift_den <- 1 + 2 * b_inner / (1 - km^2)
  U_num <- 1 + 2 * b_inner + km^2 * drift$V_sup

  sprintf(
    paste0(
      "||P_gamma^k(gamma*,.) - pi_{gamma|B~}||_TV <= ",
      "(1 - eps)^floor(alpha*k) + (U/alpha) * (num/den) * alpha^k\n",
      "  eps = %.6g; alpha = %.6g; k = %d; display = %s\n",
      "  kappa_max^LB = %.6g; sum kappa_i^LB = %.6g; q = %d\n",
      "  b = 1 - kappa_max^2 + q/2 + (1/2) sum kappa_i = %.6g\n",
      "  U = 1 + 2b + kappa_max^2 * V_sup = %.6g\n",
      "  num = 1 + 2b + kappa_max^2 * V(gamma_0) = %.6g\n",
      "  den = 1 + 2b / (1 - kappa_max^2) = %.6g"
    ),
    eps, alpha, k, display_mode,
    km, kap_sum, q,
    b_inner,
    U_num,
    drift_num,
    drift_den
  )
}

#' Sharpest displayed gamma-beta total-variation certificate.
#'
#' Assembles the Rosenthal bound from
#' \code{inst/GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md} (sharpest displayed
#' box, section 3.1) using marginal-mode \eqn{\widetilde B(\delta_2)},
#' \eqn{\varepsilon(\gamma^\star)}, and floor spectrum \eqn{\kappa_i^{\mathrm{LB}}}.
#' Does not modify \code{\link{certificate}} (restricted gamma-only route).
#'
#' @inheritParams population_mode
#' @param em_tol EM convergence tolerance passed to \code{\link{population_mode}}
#'   as \code{tol}.
#' @param delta_2 Tail budget for \eqn{r_{\mathrm{Gauss}}(n,\delta_2)} and the
#'   asymptotic full-\eqn{\pi_\gamma} correction when
#'   \code{include_full_pi_gamma = TRUE}.
#' @param k If set, evaluate the inner Rosenthal bound at this sweep count.
#' @param tol If set, invert the inner bound for a sweep count (must exceed any
#'   tail terms you require).
#' @param alpha Rosenthal tuning in \eqn{(\lambda^{\mathrm{LB}},1)}; \code{NULL}
#'   optimizes at \code{k}.
#' @param display_mode \code{"sharp"} (default) or \code{"general"}.
#' @param include_full_pi_gamma If \code{TRUE}, add \code{delta_2} for the step
#'   from \eqn{\pi_{\gamma\mid\widetilde B}} to full \eqn{\pi_\gamma}.
#' @param kappa_method Passed to \code{\link{beta_marginal_safe_set}}.
#' @param mode_method Ignored (always marginal Newton for this certificate).
#' @param delta Tail budget for \code{display_mode = "general"} minorization via
#'   \code{\link{epsilon}} (not used in the default sharp display).
#' @param verbose If \code{TRUE}, pass \code{verbose} to
#'   \code{\link{beta_marginal_safe_set}}.
#' @return An object of class \code{"gamma_beta_tv_certificate"} with
#'   \code{mode}, \code{beta_set}, \code{epsilon}, \code{floor_spectrum},
#'   \code{rosenthal}, \code{delta_2}, \code{full_pi_gamma}, \code{certified},
#'   and \code{call}.
#' @seealso \code{\link{certificate}}, \code{\link{beta_marginal_safe_set}},
#'   \code{\link{rosenthal_tv_bound}}
#' @export
gamma_beta_tv_certificate <- function(design,
                                      pfamily_list,
                                      family = gaussian(),
                                      delta_2 = 0.01,
                                      dispprior_list = NULL,
                                      k = NULL,
                                      tol = NULL,
                                      alpha = NULL,
                                      display_mode = c("sharp", "general"),
                                      include_full_pi_gamma = TRUE,
                                      estep = c("exact", "aghq", "mc"),
                                      kappa_method = c("laplace", "crude", "none"),
                                      mode_method = "marginal_newton",
                                      acceleration = c("none", "squarem"),
                                      n = 10000L,
                                      mc_seed = NULL,
                                      icm_init = TRUE,
                                      icm_tol = 1e-8,
                                      icm_maxit = 200L,
                                      em_tol = 1e-10,
                                      maxit = 200L,
                                      delta = NULL,
                                      verbose = FALSE) {
  display_mode <- match.arg(display_mode)
  estep <- match.arg(estep)
  kappa_method <- match.arg(kappa_method)
  acceleration <- match.arg(acceleration)

  mode <- population_mode(
    design = design,
    pfamily_list = pfamily_list,
    family = family,
    dispprior_list = dispprior_list,
    estep = estep,
    acceleration = acceleration,
    n = n,
    mc_seed = mc_seed,
    icm_init = icm_init,
    icm_tol = icm_tol,
    icm_maxit = icm_maxit,
    tol = em_tol,
    maxit = maxit
  )

  beta_set <- beta_marginal_safe_set(
    design = design,
    pfamily_list = pfamily_list,
    family = family,
    delta_2 = delta_2,
    dispprior_list = dispprior_list,
    kappa_method = kappa_method,
    verbose = verbose
  )

  use_closure <- identical(family$family, "gaussian") &&
    identical(estep, "exact")
  eps_obj <- if (use_closure) {
    epsilon_star(mode, method = "closure")
  } else {
    epsilon_optimize(mode, n = n, mc_seed = mc_seed)
  }

  floor_spec <- floor_coupling_spectrum(mode, beta_set)

  eps_use <- if (identical(display_mode, "sharp")) {
    eps_obj$eps_star
  } else {
    if (is.null(delta)) {
      stop("'delta' is required when display_mode = \"general\".", call. = FALSE)
    }
    epsilon(eps_obj$eps_star, delta = delta, mode = mode)$eps
  }

  d_use <- NULL
  if (identical(display_mode, "general")) {
    d_use <- epsilon(eps_obj$eps_star, delta = delta, mode = mode)$d
  }

  rosenthal <- NULL
  sweeps <- NULL

  if (!is.null(k)) {
    rosenthal <- rosenthal_tv_bound(
      k = k,
      eps = eps_use,
      spectrum = floor_spec,
      alpha = alpha,
      mode = mode,
      display_mode = display_mode,
      d = d_use
    )
  } else if (!is.null(tol)) {
    sweeps <- .c05_rosenthal_sweeps_for_tol(
      tol = tol,
      eps = eps_use,
      spectrum = floor_spec,
      display_mode = display_mode,
      mode = mode,
      d = d_use,
      alpha = alpha
    )
    sweeps$tol <- tol
    rosenthal <- sweeps$rosenthal
    k <- sweeps$k
  }

  inner_bound <- if (!is.null(rosenthal)) rosenthal$bound else NULL
  full_bound <- if (isTRUE(include_full_pi_gamma) && !is.null(inner_bound)) {
    inner_bound + delta_2
  } else {
    inner_bound
  }

  certified_kappa <- isTRUE(beta_set$certified$kappa)
  certified <- list(
    gamma_em = isTRUE(mode$converged),
    epsilon = isTRUE(eps_obj$certified),
    kappa_lb = certified_kappa,
    delta_2 = "asymptotic_laplace",
    sharpest_display = if (identical(display_mode, "sharp")) "limit" else "general"
  )

  structure(
    list(
      mode = mode,
      beta_set = beta_set,
      epsilon = eps_obj,
      floor_spectrum = floor_spec,
      rosenthal = rosenthal,
      delta_2 = delta_2,
      delta = delta,
      display_mode = display_mode,
      k = k,
      sweeps = sweeps,
      full_pi_gamma = list(
        include = isTRUE(include_full_pi_gamma),
        inner_bound = inner_bound,
        full_bound = full_bound,
        tail_mass_label = "asymptotic Laplace at r_Gauss"
      ),
      certified = certified,
      eps_use = eps_use,
      call = match.call()
    ),
    class = "gamma_beta_tv_certificate"
  )
}

#' @export
print.gamma_beta_tv_certificate <- function(x, digits = 4, ...) {
  cat("Gamma-beta TV certificate (Rosenthal / sharpest displayed route)\n\n")
  cat("  gamma* EM: ", x$mode$iterations,
      if (isTRUE(x$certified$gamma_em)) " (converged)" else " (not converged)",
      "\n", sep = "")
  bd <- x$beta_set$mode$beta_dagger
  cat("  beta_dagger (marginal Newton): ",
      paste(signif(bd, digits), collapse = ", "), "\n", sep = "")
  lv <- x$beta_set$level
  cat("  delta_2: ", x$delta_2,
      "  r_Gauss: ", signif(lv$r_gauss, digits),
      "  min omega: ",
      signif(min(vapply(x$beta_set$per_group, function(p) p$omega, 0)), digits),
      "\n", sep = "")

  fs <- x$floor_spectrum
  cat("  kappa_max^LB: ", signif(fs$kappa_max_lb, digits),
      "  lambda^LB: ", signif(fs$lambda_lb, digits),
      "  q: ", fs$q, "\n", sep = "")
  cat("  kappa_i^LB: ",
      paste(signif(fs$kappa_lb, digits), collapse = ", "), "\n", sep = "")

  cat("\n  eps(gamma*): ", signif(x$epsilon$eps_star, digits),
      "  eps (bound): ", signif(x$eps_use, digits),
      "  display: ", x$display_mode, "\n", sep = "")

  cat("\n  certified: gamma_em=", x$certified$gamma_em,
      " eps=", x$certified$epsilon,
      " kappa_lb=", if (x$certified$kappa_lb) "TRUE" else "diagnostic",
      " delta_2=", x$certified$delta_2,
      " sharpest=", x$certified$sharpest_display, "\n", sep = "")

  if (!is.null(x$rosenthal)) {
    ros <- x$rosenthal
    cat("\n--- Expanded Rosenthal box (inner: pi_{gamma|B~}) ---\n")
    cat(.c05_rosenthal_expanded_box(
      drift = ros$drift_block,
      eps = ros$eps,
      k = ros$k,
      alpha = ros$alpha,
      display_mode = x$display_mode
    ), "\n", sep = "")
    cat("\n  inner bound: ", signif(ros$bound, digits),
        "  (minorization: ", signif(ros$minorization, digits),
        " + drift: ", signif(ros$drift, digits), ")\n", sep = "")
    if (isTRUE(x$full_pi_gamma$include)) {
      cat("  full pi_gamma bound (inner + delta_2): ",
          signif(x$full_pi_gamma$full_bound, digits),
          "  [", x$full_pi_gamma$tail_mass_label, "]\n", sep = "")
    }
  } else if (!is.null(x$sweeps)) {
    cat("\n  sweeps for tol=", signif(x$sweeps$tol, digits),
        ": k=", x$sweeps$k,
        "  inner bound=", signif(x$sweeps$bound, digits), "\n", sep = "")
  } else {
    cat("\n  (Set 'k' or 'tol' to evaluate the Rosenthal bound.)\n")
  }

  invisible(x)
}

#' @export
format.gamma_beta_tv_certificate <- function(x, ...) {
  if (is.null(x$rosenthal)) {
    return("gamma_beta_tv_certificate (bound not evaluated; set k or tol)")
  }
  .c05_rosenthal_expanded_box(
    drift = x$rosenthal$drift_block,
    eps = x$rosenthal$eps,
    k = x$rosenthal$k,
    alpha = x$rosenthal$alpha,
    display_mode = x$display_mode
  )
}
