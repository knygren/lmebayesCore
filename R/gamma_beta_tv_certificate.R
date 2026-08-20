## Sharpest displayed gamma-beta TV certificate (orchestrator).

#' @noRd
.c05_rosenthal_expanded_box <- function(drift, eps, k, alpha, display_mode) {
  km <- drift$kappa_max_lb
  q <- drift$q
  kap_sum <- sum(drift$kappa_lb)
  b_inner <- drift$b_drift
  drift_num <- drift$drift_numerator
  drift_den <- drift$drift_denominator
  U_num <- drift$U

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
#' Orchestrates categories 1--8: \code{\link{population_mode}},
#' \code{\link{epsilon_star}}, \code{\link{beta_marginal_mode}},
#' \code{\link{beta_marginal_safe_set}}, \code{\link{group_precision_floor}},
#' \code{\link{floor_coupling_eigenvalues}}, and
#' \code{\link{optimal_rosenthal_tv_bound}}.
#'
#' @inheritParams population_mode
#' @inheritParams optimal_rosenthal_tv_bound
#' @param em_tol EM convergence tolerance passed to \code{\link{population_mode}}
#'   as \code{tol}.
#' @param delta Tail budget for \code{display_mode = "general"} minorization via
#'   \code{\link{epsilon}} (not used in the default sharp display).
#' @param tol Alias for \code{inner_tol} when inverting for sweep count.
#' @param include_full_pi_gamma If \code{TRUE}, report full \eqn{\pi_\gamma}
#'   bound as inner bound plus \code{delta_2}.
#' @param kappa_method Passed to \code{\link{group_precision_floor}}.
#' @param verbose If \code{TRUE}, pass \code{verbose} to
#'   \code{\link{beta_marginal_safe_set}}.
#' @return An object of class \code{"gamma_beta_tv_certificate"}.
#' @seealso \code{\link{certificate}}, \code{\link{optimal_rosenthal_tv_bound}}
#' @export
gamma_beta_tv_certificate <- function(design,
                                      pfamily_list,
                                      family = gaussian(),
                                      delta_2 = 0.01,
                                      dispprior_list = NULL,
                                      k = NULL,
                                      tol = NULL,
                                      inner_tol = NULL,
                                      total_tol = NULL,
                                      alpha = NULL,
                                      display_mode = c("sharp", "general"),
                                      include_full_pi_gamma = TRUE,
                                      estep = c("exact", "aghq", "mc"),
                                      kappa_method = c("laplace", "crude", "none"),
                                      acceleration = c("none", "squarem"),
                                      n = 10000L,
                                      mc_seed = NULL,
                                      icm_init = TRUE,
                                      icm_tol = 1e-8,
                                      icm_maxit = 200L,
                                      em_tol = 1e-10,
                                      maxit = 200L,
                                      delta = NULL,
                                      optimize_alpha = TRUE,
                                      verbose = FALSE) {
  display_mode <- match.arg(display_mode)
  estep <- match.arg(estep)
  kappa_method <- match.arg(kappa_method)
  acceleration <- match.arg(acceleration)

  if (is.null(inner_tol) && !is.null(tol)) {
    inner_tol <- tol
  }

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

  beta_mode <- beta_marginal_mode(
    design = design,
    pfamily_list = pfamily_list,
    family = family,
    dispprior_list = dispprior_list
  )

  beta_set <- beta_marginal_safe_set(
    delta_2 = delta_2,
    beta_mode = beta_mode,
    verbose = verbose
  )

  floor_obj <- group_precision_floor(
    beta_mode = beta_mode,
    beta_set = beta_set,
    kappa_method = kappa_method
  )

  use_closure <- identical(family$family, "gaussian") &&
    identical(estep, "exact")
  eps_obj <- if (use_closure) {
    epsilon_star(mode, method = "closure")
  } else {
    epsilon_optimize(mode, n = n, mc_seed = mc_seed)
  }

  eigenvalues <- floor_coupling_eigenvalues(mode, floor_obj)

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

  optimal <- NULL
  if (!is.null(k) || !is.null(inner_tol) || !is.null(total_tol)) {
    optimal <- optimal_rosenthal_tv_bound(
      mode = mode,
      eigenvalues = eigenvalues,
      eps = eps_use,
      k = k,
      alpha = alpha,
      inner_tol = inner_tol,
      total_tol = total_tol,
      delta_2 = if (isTRUE(include_full_pi_gamma)) delta_2 else 0,
      include_full_pi_gamma = include_full_pi_gamma,
      optimize_alpha = optimize_alpha,
      display_mode = display_mode,
      d = d_use
    )
  }

  rosenthal <- if (!is.null(optimal)) optimal$rosenthal else NULL

  certified <- list(
    gamma_em = isTRUE(mode$converged),
    epsilon = isTRUE(eps_obj$certified),
    kappa_lb = isTRUE(floor_obj$certified),
    delta_2 = "asymptotic_laplace",
    sharpest_display = if (identical(display_mode, "sharp")) "limit" else "general"
  )

  structure(
    list(
      mode = mode,
      beta_mode = beta_mode,
      beta_set = beta_set,
      floor = floor_obj,
      epsilon = eps_obj,
      eigenvalues = eigenvalues,
      optimal = optimal,
      rosenthal = rosenthal,
      delta_2 = delta_2,
      delta = delta,
      display_mode = display_mode,
      full_pi_gamma = if (!is.null(optimal)) {
        list(
          include = isTRUE(include_full_pi_gamma),
          inner_bound = optimal$inner_bound,
          full_bound = optimal$full_bound,
          tail_mass_label = "asymptotic Laplace at r_Gauss"
        )
      } else {
        NULL
      },
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
  cat("  beta_dagger ||beta||: ",
      signif(sqrt(sum(x$beta_mode$beta^2)), digits), "\n", sep = "")
  lv <- x$beta_set$level
  cat("  delta_2: ", x$delta_2,
      "  r_Gauss: ", signif(lv$r_gauss, digits),
      "  min omega: ",
      signif(min(vapply(x$floor$per_group, function(p) p$omega, 0)), digits),
      "\n", sep = "")

  ev <- x$eigenvalues
  drift <- if (!is.null(x$optimal)) {
    x$optimal$drift
  } else {
    rosenthal_drift_constants(ev, display_mode = x$display_mode)
  }
  cat("  kappa_max^LB: ", signif(ev$kappa_max_lb, digits),
      "  lambda^LB: ", signif(drift$lambda_lb, digits),
      "  q: ", ev$q, "\n", sep = "")
  cat("  kappa_i^LB: ",
      paste(signif(ev$kappa_lb, digits), collapse = ", "), "\n", sep = "")

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
      drift = ros$drift_constants,
      eps = ros$eps,
      k = ros$k,
      alpha = ros$alpha,
      display_mode = x$display_mode
    ), "\n", sep = "")
    cat("\n  inner bound: ", signif(ros$bound, digits),
        "  (minorization: ", signif(ros$minorization, digits),
        " + drift: ", signif(ros$drift, digits), ")\n", sep = "")
    if (!is.null(x$full_pi_gamma) && isTRUE(x$full_pi_gamma$include)) {
      cat("  full pi_gamma bound (inner + delta_2): ",
          signif(x$full_pi_gamma$full_bound, digits),
          "  [", x$full_pi_gamma$tail_mass_label, "]\n", sep = "")
    }
  } else {
    cat("\n  (Set 'k', 'inner_tol', or 'total_tol' to evaluate the bound.)\n")
  }

  invisible(x)
}

#' @export
format.gamma_beta_tv_certificate <- function(x, ...) {
  if (is.null(x$rosenthal)) {
    return("gamma_beta_tv_certificate (bound not evaluated; set k or tol)")
  }
  .c05_rosenthal_expanded_box(
    drift = x$rosenthal$drift_constants,
    eps = x$rosenthal$eps,
    k = x$rosenthal$k,
    alpha = x$rosenthal$alpha,
    display_mode = x$display_mode
  )
}
