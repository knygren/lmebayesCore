#' Marginal-mode beta safe set for Rosenthal / TV certificates
#'
#' Build the beta-safe set \eqn{\widetilde B(\delta_2) = \{\Xi(\beta)
#' \le r_{\mathrm{Gauss}}(n,\delta_2)\}} at fixed \eqn{\beta^\dagger} from
#' \code{\link{beta_marginal_mode}}. Does not compute data-precision floors;
#' use \code{\link{group_precision_floor}} for \eqn{\Gamma_j^{\mathrm{LB}}}.
#'
#' @param design A \code{\link{model_setup}} list (used if \code{beta_mode} is
#'   \code{NULL}).
#' @param pfamily_list Block~2 prior list (used if \code{beta_mode} is \code{NULL}).
#' @param family A \code{\link[stats]{family}} object (used if \code{beta_mode} is
#'   \code{NULL}).
#' @param delta_2 Tail mass budget in \eqn{(0,1)} for the Laplace reference at
#'   \eqn{r_{\mathrm{Gauss}}}.
#' @param beta_mode Optional result of \code{\link{beta_marginal_mode}}; if
#'   \code{NULL}, marginal mode is computed from \code{design}.
#' @param dispprior_list,offset,weights Passed to \code{\link{beta_marginal_mode}}
#'   when \code{beta_mode} is \code{NULL}.
#' @param tol,maxit Passed to \code{\link{beta_marginal_mode}} when
#'   \code{beta_mode} is \code{NULL}.
#' @param verbose If \code{TRUE}, emit \code{\link{message}} notes about the
#'   calibrated level.
#' @return An object of class \code{"beta_marginal_safe_set"} with \code{level},
#'   \code{beta_mode}, \code{delta_2}, \code{P_b}, and certification metadata.
#' @seealso \code{\link{beta_marginal_mode}}, \code{\link{group_precision_floor}},
#'   \code{\link{gamma_beta_tv_certificate}}
#' @export
beta_marginal_safe_set <- function(design = NULL,
                                   pfamily_list = NULL,
                                   family = NULL,
                                   delta_2 = 0.01,
                                   beta_mode = NULL,
                                   dispprior_list = NULL,
                                   offset = NULL,
                                   weights = NULL,
                                   tol = 1e-10,
                                   maxit = 200L,
                                   verbose = FALSE) {
  if (is.null(beta_mode)) {
    if (is.null(design) || is.null(pfamily_list) || is.null(family)) {
      stop(
        "'design', 'pfamily_list', and 'family' are required when ",
        "'beta_mode' is NULL.",
        call. = FALSE
      )
    }
    beta_mode <- beta_marginal_mode(
      design = design,
      pfamily_list = pfamily_list,
      family = family,
      dispprior_list = dispprior_list,
      offset = offset,
      weights = weights,
      tol = tol,
      maxit = maxit
    )
  } else if (!inherits(beta_mode, "beta_marginal_mode")) {
    stop("'beta_mode' must be from beta_marginal_mode().", call. = FALSE)
  }

  prep <- beta_mode$prep
  if (is.null(prep)) {
    stop("'beta_mode' is missing internal 'prep'; recompute with beta_marginal_mode().",
         call. = FALSE)
  }

  level <- .group_floor_resolve_levels(
    epsilon = delta_2,
    level_method = "r_gauss_joint",
    union_over = "groups",
    J = prep$J,
    p_re = prep$p_re,
    n_obs_total = prep$n_obs_total,
    kappa = rep(0, prep$J),
    gaussian = prep$family_hook$gaussian,
    kappa_method = "none",
    inflate_kappa = FALSE
  )

  if (isTRUE(verbose)) {
    message(
      "beta_marginal_safe_set: r_Gauss = ", signif(level$r_gauss, 4),
      " at delta_2 = ", delta_2
    )
  }

  structure(
    list(
      level = level,
      beta_mode = beta_mode,
      delta_2 = delta_2,
      P_b = beta_mode$P_b,
      certified = list(delta_2 = "asymptotic_laplace"),
      family = beta_mode$family,
      call = match.call()
    ),
    class = "beta_marginal_safe_set"
  )
}

#' @export
print.beta_marginal_safe_set <- function(x, digits = 4, ...) {
  lv <- x$level
  cat("Marginal-mode beta safe set (r_Gauss joint level)\n")
  cat("  family:", x$family$family, "(", x$family$link, ")\n", sep = "")
  cat("  delta_2:", lv$epsilon,
      "  r_Gauss:", signif(lv$r_gauss, digits),
      "  s:", signif(lv$s[1L], digits), "\n", sep = "")
  cat("  laplace_tail_mass (design):", lv$laplace_tail_mass, "\n")
  invisible(x)
}
