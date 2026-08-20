#' Marginal-mode beta posterior mode (gamma integrated out)
#'
#' Finds the Newton mode \eqn{\beta^\dagger} of the gamma-integrated marginal
#' posterior on stacked group coefficients. Reusable input for
#' \code{\link{beta_marginal_safe_set}} and \code{\link{group_precision_floor}}.
#'
#' @param design A \code{\link{model_setup}} list.
#' @param pfamily_list Block~2 prior list from \code{\link{pfamily_list}()}.
#' @param family A \code{\link[stats]{family}} object.
#' @param dispprior_list Optional Block~1 dispersion prior for \code{gaussian()}.
#' @param offset,weights Observation \code{offset} and \code{weights}; default
#'   from \code{design}.
#' @param tol,maxit Newton controls (unused beyond initial solve).
#' @return An object of class \code{"beta_marginal_mode"} with \code{beta},
#'   \code{hessian}, \code{f_mode}, \code{engine}, \code{prior_int},
#'   \code{group_data}, \code{P_b}, and \code{family}.
#' @seealso \code{\link{beta_marginal_safe_set}}, \code{\link{group_precision_floor}}
#' @export
beta_marginal_mode <- function(design,
                               pfamily_list,
                               family,
                               dispprior_list = NULL,
                               offset = NULL,
                               weights = NULL,
                               tol = 1e-10,
                               maxit = 200L) {
  prep <- .group_floor_prepare(
    design = design,
    pfamily_list = pfamily_list,
    family = family,
    dispprior_list = dispprior_list,
    offset = offset,
    weights = weights,
    fn_name = "beta_marginal_mode"
  )

  mode_res <- .c05_beta_marginal_mode(
    prior_int = prep$prior_int,
    group_data = prep$group_data,
    family_hook = prep$family_hook
  )

  structure(
    list(
      beta = mode_res$beta,
      beta_dagger = mode_res$beta,
      hessian = mode_res$hessian,
      f_mode = mode_res$f_mode,
      engine = mode_res$engine,
      prior_int = prep$prior_int,
      group_data = prep$group_data,
      prep = prep,
      P_b = prep$prior_int$P_b,
      method = "marginal_newton",
      family = prep$family_hook$family,
      call = match.call()
    ),
    class = "beta_marginal_mode"
  )
}

#' @export
print.beta_marginal_mode <- function(x, digits = 4, ...) {
  cat("Marginal-mode beta (Newton on gamma-integrated posterior)\n")
  cat("  family:", x$family$family, "(", x$family$link, ")\n", sep = "")
  cat("  f_mode:", signif(x$f_mode, digits), "\n")
  cat("  ||beta||:", signif(sqrt(sum(x$beta^2)), digits), "\n")
  invisible(x)
}
