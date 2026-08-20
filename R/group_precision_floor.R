#' Certified group-level data precision lower bounds
#'
#' Compute per-group lower-bound data precision matrices \eqn{\Gamma_j} for a
#' two-block GLMM by integrating population parameters out of the prior,
#' finding the joint posterior mode, certifying per-group profile level sets at
#' total escape budget \code{epsilon}, and minimizing observation weights on
#' those sets (exactly for log-concave / monotone canonical GLM weights).
#'
#' @details
#' Theoretical basis: \code{inst/LAPLACE_PROFILE_GROUP_MARGINAL_BOUND.md}
#' (Assumption (L): Gaussian in \eqn{\beta_{-j}} given \eqn{\beta_j};
#' approximate marginal profile with a log-determinant correction; Prop 2 on
#' the profile level only; \eqn{\kappa} prices log-det variation).
#' For joint \eqn{r_{\mathrm{Gauss}}} level sets see
#' \code{\link{beta_marginal_safe_set}} and
#' \code{inst/BETA_MARGINAL_MODE_LEVELSET.md}.
#'
#' @param design A \code{\link{model_setup}} list with \code{y}, \code{D},
#'   \code{group}, \code{W}, and \code{groupef.names}.
#' @param pfamily_list Block~2 prior list from \code{\link{pfamily_list}()}.
#' @param family A \code{\link[stats]{family}} object.
#' @param epsilon Total escape probability budget in \eqn{(0,1)}.
#' @param dispprior_list Optional Block~1 dispersion prior for \code{gaussian()}.
#' @param offset,weights Observation \code{offset} and \code{weights}; default
#'   from \code{design}.
#' @param kappa_method \code{"crude"} (certified), \code{"laplace"} (diagnostic),
#'   or \code{"none"}.
#' @param union_over \code{"groups"} or \code{"functionals"} (Prop 2 schemes).
#' @param mode_method \code{"icm"} (default, legacy) or \code{"marginal_newton"}.
#' @param level_method \code{"prop2_groups"} (default), \code{"prop2_functionals"},
#'   or \code{"r_gauss_joint"}.
#' @param tol,maxit Controls for ICM or marginal Newton mode finding.
#' @param verbose If \code{TRUE}, emit \code{\link{message}} notes.
#' @return An object of class \code{"group_precision_floor"}: a list with
#'   \code{Gamma_lb}, \code{P_b}, \code{level}, \code{mode}, \code{per_group},
#'   and \code{certified}.
#' @seealso \code{\link{beta_marginal_safe_set}}, \code{\link{glmerb_posterior_mode}}
#' @export
group_precision_floor <- function(design,
                                pfamily_list,
                                family,
                                epsilon = 0.01,
                                dispprior_list = NULL,
                                offset = NULL,
                                weights = NULL,
                                kappa_method = c("crude", "laplace", "none"),
                                union_over = c("groups", "functionals"),
                                mode_method = c("icm", "marginal_newton"),
                                level_method = c("prop2_groups", "prop2_functionals", "r_gauss_joint"),
                                tol = 1e-10,
                                maxit = 200L,
                                verbose = FALSE) {
  kappa_method <- match.arg(kappa_method)
  union_over <- match.arg(union_over)
  mode_method <- match.arg(mode_method)
  level_method <- match.arg(level_method)

  prep <- .group_floor_prepare(
    design = design,
    pfamily_list = pfamily_list,
    family = family,
    dispprior_list = dispprior_list,
    offset = offset,
    weights = weights,
    fn_name = "group_precision_floor"
  )

  mode_info <- .group_floor_resolve_mode(
    prep = prep,
    design = design,
    mode_method = mode_method,
    tol = tol,
    maxit = maxit
  )

  engine <- .group_floor_engine(
    prior_int = prep$prior_int,
    group_data = prep$group_data,
    family_hook = prep$family_hook,
    b_dag = mode_info$beta
  )

  kappa <- .group_floor_compute_kappa(
    prep = prep,
    engine = engine,
    kappa_method = kappa_method,
    epsilon = epsilon,
    level_method = level_method
  )

  level <- .group_floor_resolve_levels(
    epsilon = epsilon,
    level_method = level_method,
    union_over = union_over,
    J = prep$J,
    p_re = prep$p_re,
    n_obs_total = prep$n_obs_total,
    kappa = kappa,
    gaussian = prep$family_hook$gaussian,
    kappa_method = kappa_method,
    inflate_kappa = !identical(level_method, "r_gauss_joint")
  )

  if (identical(level_method, "prop2_groups") &&
      identical(kappa_method, "crude") &&
      any(level$R > level$r_joint + 1e-8)) {
    if (isTRUE(verbose)) {
      message(
        "Crude kappa inflates R_j above the joint level r(n): ",
        "max(R)/r_joint = ",
        signif(max(level$R) / level$r_joint, 4),
        ". The per-group certificate may be no tighter than the joint scheme."
      )
    }
  }

  uo <- if (identical(level_method, "prop2_functionals")) {
    "functionals"
  } else {
    "groups"
  }

  floors <- .group_floor_assemble_floors(
    prep = prep,
    engine = engine,
    level = level,
    union_over = uo
  )

  certified <- identical(kappa_method, "crude") ||
    (prep$family_hook$gaussian && identical(kappa_method, "none"))

  structure(
    list(
      Gamma_lb = floors$Gamma_lb,
      P_b = prep$prior_int$P_b,
      level = level,
      mode = list(
        beta = engine$b_dag,
        f_mode = engine$f_dag,
        hessian = mode_info$hessian,
        method = mode_info$method,
        icm = mode_info$icm
      ),
      per_group = floors$per_group,
      certified = certified,
      kappa_method = kappa_method,
      mode_method = mode_method,
      level_method = level_method,
      union_over = uo,
      family = prep$family_hook$family,
      call = match.call()
    ),
    class = "group_precision_floor"
  )
}

#' @export
print.group_precision_floor <- function(x, digits = 4, ...) {
  lv <- x$level
  cat("Group-level precision floor certificate\n")
  cat("  family:", x$family$family, "(", x$family$link, ")\n", sep = "")
  cat("  epsilon:", lv$epsilon, " scheme:", lv$scheme,
      " mode:", x$mode_method, " level:", x$level_method,
      " kappa_method:", x$kappa_method,
      " certified:", x$certified, "\n")
  if (!is.na(lv$r_group)) {
    cat("  r_group:", signif(lv$r_group, digits),
        " r_joint:", signif(lv$r_joint, digits), "\n", sep = "")
  }
  if (!is.null(lv$r_gauss)) {
    cat("  r_Gauss:", signif(lv$r_gauss, digits), "\n", sep = "")
  }
  tab <- do.call(rbind, lapply(names(x$per_group), function(g) {
    pg <- x$per_group[[g]]
    data.frame(
      group = g,
      R = lv$R[match(g, names(x$per_group))],
      s = lv$s[match(g, names(x$per_group))],
      eta_star_max = max(pg$es),
      wbar_min = min(pg$wbar),
      omega = pg$omega,
      stringsAsFactors = FALSE
    )
  }))
  print(tab, row.names = FALSE, digits = digits)
  invisible(x)
}
