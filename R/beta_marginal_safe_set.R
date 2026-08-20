#' Marginal-mode beta safe set for Rosenthal / TV certificates
#'
#' Build the certified beta-safe set \eqn{\widetilde B(\delta_2) = \{\Xi(\beta)
#' \le r_{\mathrm{Gauss}}(n,\delta_2)\}} using the gamma-integrated marginal
#' posterior, Newton mode \eqn{\beta^\dagger}, and joint Gaussian-reference
#' level calibration. See \code{inst/BETA_MARGINAL_MODE_LEVELSET.md}.
#'
#' @inheritParams group_precision_floor
#' @param delta_2 Tail mass budget in \eqn{(0,1)} for the Laplace reference at
#'   \eqn{r_{\mathrm{Gauss}}}.
#' @param kappa_method Passed through for metadata only; joint \eqn{r} does not
#'   inflate by \eqn{2\kappa_j} on this route.
#' @return An object of class \code{"beta_marginal_safe_set"} with
#'   \code{Gamma_lb}, \code{P_b}, \code{level}, marginal \code{mode}
#'   (\code{beta_dagger}, \code{hessian}), \code{per_group}, \code{engine},
#'   and certification flags.
#' @seealso \code{\link{group_precision_floor}},
#'   \code{\link{gamma_beta_tv_certificate}}
#' @export
beta_marginal_safe_set <- function(design,
                                  pfamily_list,
                                  family,
                                  delta_2 = 0.01,
                                  dispprior_list = NULL,
                                  offset = NULL,
                                  weights = NULL,
                                  kappa_method = c("laplace", "crude", "none"),
                                  tol = 1e-10,
                                  maxit = 200L,
                                  verbose = FALSE) {
  kappa_method <- match.arg(kappa_method)

  prep <- .group_floor_prepare(
    design = design,
    pfamily_list = pfamily_list,
    family = family,
    dispprior_list = dispprior_list,
    offset = offset,
    weights = weights,
    fn_name = "beta_marginal_safe_set"
  )

  mode_info <- .group_floor_resolve_mode(
    prep = prep,
    design = design,
    mode_method = "marginal_newton",
    tol = tol,
    maxit = maxit
  )

  engine <- .group_floor_engine(
    prior_int = prep$prior_int,
    group_data = prep$group_data,
    family_hook = prep$family_hook,
    b_dag = mode_info$beta
  )

  if (is.null(mode_info$hessian)) {
    mode_info$hessian <- engine$He_f(engine$b_dag)
  }

  kappa <- .group_floor_compute_kappa(
    prep = prep,
    engine = engine,
    kappa_method = kappa_method,
    epsilon = delta_2,
    level_method = "r_gauss_joint"
  )

  level <- .group_floor_resolve_levels(
    epsilon = delta_2,
    level_method = "r_gauss_joint",
    union_over = "groups",
    J = prep$J,
    p_re = prep$p_re,
    n_obs_total = prep$n_obs_total,
    kappa = kappa,
    gaussian = prep$family_hook$gaussian,
    kappa_method = kappa_method,
    inflate_kappa = FALSE
  )

  floors <- .group_floor_assemble_floors(
    prep = prep,
    engine = engine,
    level = level,
    union_over = "groups"
  )

  certified_kappa <- identical(kappa_method, "crude") ||
    (prep$family_hook$gaussian && identical(kappa_method, "none"))

  if (isTRUE(verbose) && !is.null(mode_info$icm)) {
    message("beta_marginal_safe_set uses marginal Newton mode, not ICM.")
  }

  structure(
    list(
      Gamma_lb = floors$Gamma_lb,
      P_b = prep$prior_int$P_b,
      level = level,
      mode = list(
        beta_dagger = engine$b_dag,
        f_mode = engine$f_dag,
        hessian = mode_info$hessian,
        method = mode_info$method
      ),
      per_group = floors$per_group,
      engine = engine,
      delta_2 = delta_2,
      certified = list(
        kappa = certified_kappa,
        delta_2 = "asymptotic_laplace"
      ),
      kappa_method = kappa_method,
      family = prep$family_hook$family,
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
  cat("  laplace_tail_mass (design):", lv$laplace_tail_mass,
      "  kappa_method:", x$kappa_method, "\n")
  tab <- do.call(rbind, lapply(names(x$per_group), function(g) {
    pg <- x$per_group[[g]]
    data.frame(
      group = g,
      eta_star_max = max(pg$es),
      wbar_min = min(pg$wbar),
      omega = pg$omega,
      stringsAsFactors = FALSE
    )
  }))
  print(tab, row.names = FALSE, digits = digits)
  invisible(x)
}
