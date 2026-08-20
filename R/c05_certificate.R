## End-to-end restricted Gibbs minorization certificate.

#' Restricted two-block Gibbs minorization certificate (Chapter C05).
#'
#' Computes the EM fixed point \eqn{\gamma^\star}, the mode profile
#' \eqn{\varepsilon(\gamma^\star)}, and the spectrum-calibrated restricted-chain
#' multiplier \eqn{\varepsilon = \varepsilon_d Q(\widetilde C_d)}.
#'
#' @details
#' The certificate applies to the \emph{restricted} Gibbs chain on
#' \eqn{\widetilde C_d}, not necessarily the production sampler; see
#' \code{inst/CHAPTER_C05_IMPLEMENTATION.md} section 8. Set sizing inverts
#' the weighted-\eqn{\chi^2} escape tail on the coupling spectrum at
#' \eqn{\gamma^\star} (section 4A.3). Outside Gaussian closure the exact
#' escape mass may differ from \code{delta}.
#'
#' @inheritParams population_mode
#' @param em_tol EM convergence tolerance passed to \code{\link{population_mode}}
#'   as \code{tol} (distinct from sweep-count \code{tol} below).
#' @param delta Tail / escape probability budget in \code{(0, 1)}.
#' @param tol If supplied, compute a sweep count \code{n} such that
#'   \code{(1 - eps)^n + delta <= tol} (Theorem 2).
#' @return An object of class \code{"minorization_certificate"}.
#' @seealso \code{\link{population_mode}},
#'   \code{\link{epsilon}}, \code{\link{deficiency_spectrum}},
#'   \code{\link{gamma_beta_tv_certificate}}
#' @export
certificate <- function(design,
                                         pfamily_list,
                                         family = gaussian(),
                                         delta = 0.01,
                                         dispprior_list = NULL,
                                         estep = c("exact", "aghq", "mc"),
                                         acceleration = c("none", "squarem"),
                                         icm_init = TRUE,
                                         icm_tol = 1e-8,
                                         icm_maxit = 200L,
                                         n = 10000L,
                                         mc_seed = NULL,
                                         tol = NULL,
                                         em_tol = 1e-10,
                                         maxit = 200L) {
  estep <- match.arg(estep)
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

  use_closure <- identical(family$family, "gaussian") &&
    identical(estep, "exact")
  eps_star <- if (use_closure) {
    epsilon_star(mode, method = "closure")
  } else {
    epsilon_optimize(
      mode,
      n = n,
      mc_seed = mc_seed
    )
  }

  constants <- epsilon(
    eps_star = eps_star$eps_star,
    delta = delta,
    mode = mode
  )

  sweeps <- NULL
  if (!is.null(tol)) {
    if (!(tol > delta)) {
      stop(
        "'tol' must exceed 'delta' (Theorem 2 bound never drops below delta).",
        call. = FALSE
      )
    }
    sweeps <- list(
      tol = tol,
      n = ceiling(log(tol - delta) / log1p(-constants$eps))
    )
  }

  certified <- isTRUE(eps_star$certified) && isTRUE(mode$converged)

  structure(
    list(
      mode = mode,
      epsilon = eps_star,
      constants = constants,
      sweeps = sweeps,
      delta = delta,
      certified = certified,
      route = "spectrum_calibrated",
      call = match.call()
    ),
    class = "minorization_certificate"
  )
}

#' @export
print.minorization_certificate <- function(x, digits = 4, ...) {
  cat("Restricted Gibbs minorization certificate (spectrum-calibrated route)\n\n")
  cat("  gamma* EM iterations: ", x$mode$iterations,
      if (isTRUE(x$mode$converged)) " (converged)" else " (not converged)",
      "\n", sep = "")
  cat("  stationarity ||P11(M(gamma*)-gamma*)||: ",
      signif(x$mode$stationarity, digits), "\n", sep = "")
  cat("  kappa (EM rate): ", signif(x$mode$kappa, digits),
      "  q: ", x$mode$q, "\n", sep = "")
  if (!is.null(x$constants$spectrum$kappa_max)) {
    cat("  kappa_max (spectrum): ",
        signif(x$constants$spectrum$kappa_max, digits), "\n", sep = "")
  }
  cat("  eps(gamma*): ", signif(x$epsilon$eps_star, digits),
      " [", x$epsilon$method, "]\n", sep = "")
  cat("  delta: ", signif(x$delta, digits),
      "  r: ", signif(x$constants$r, digits),
      "  d: ", signif(x$constants$d, digits), "\n", sep = "")
  if (!is.null(x$constants$spectrum$kappa)) {
    cat("  kappa (A): ", paste(signif(x$constants$spectrum$kappa, digits),
                            collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$constants$spectrum$weights)) {
    cat("  weights w=k/(1-k) (B): ",
        paste(signif(x$constants$spectrum$weights, digits), collapse = ", "),
        "\n", sep = "")
  }
  if (!is.null(x$constants$cal$escape_mass)) {
    cat("  escape mass (ref): ",
        signif(x$constants$cal$escape_mass, digits), "\n", sep = "")
  }
  cat("  eps_d: ", signif(x$constants$eps_d, digits),
      "  Q(C) lb: ", signif(x$constants$Q_mass_lb, digits), "\n", sep = "")
  cat("  eps (restricted multiplier): ", signif(x$constants$eps, digits), "\n",
      sep = "")
  if (!is.null(x$sweeps)) {
    cat("  n sweeps (tol = ", signif(x$sweeps$tol, digits), "): ",
        x$sweeps$n, "\n", sep = "")
  }
  cat("  certified: ", x$certified, "\n", sep = "")
  cat("\nBound is for the restricted chain on C_tilde_d; see ",
      "inst/CHAPTER_C05_IMPLEMENTATION.md section 8.\n", sep = "")
  invisible(x)
}
