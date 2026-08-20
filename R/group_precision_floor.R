#' Certified group-level data precision lower bounds
#'
#' Compute per-group lower-bound data precision matrices \eqn{\Gamma_j} on a
#' certified beta level set. Typically used after
#' \code{\link{beta_marginal_mode}} and \code{\link{beta_marginal_safe_set}}.
#'
#' @param design Used when \code{beta_mode} and \code{beta_set} are both
#'   \code{NULL} (legacy standalone path).
#' @param pfamily_list Block~2 prior list.
#' @param family A \code{\link[stats]{family}} object.
#' @param beta_mode Result of \code{\link{beta_marginal_mode}} (preferred).
#' @param beta_set Result of \code{\link{beta_marginal_safe_set}} supplying the
#'   joint level \eqn{r_{\mathrm{Gauss}}}; if \code{NULL}, \code{epsilon} and
#'   \code{level_method} define the set.
#' @param epsilon Escape budget when building levels without \code{beta_set}.
#' @param kappa_method \code{"crude"} (certified), \code{"laplace"}, or
#'   \code{"none"}.
#' @param union_over \code{"groups"} or \code{"functionals"} (Prop 2 schemes).
#' @param mode_method Default \code{"marginal_newton"}; ignored when
#'   \code{beta_mode} is supplied.
#' @param level_method Default \code{"r_gauss_joint"} when \code{beta_set} is
#'   supplied; otherwise \code{"prop2_groups"} for legacy.
#' @inheritParams beta_marginal_mode
#' @param verbose If \code{TRUE}, emit \code{\link{message}} notes (legacy Prop
#'   2 path).
#' @return An object of class \code{"group_precision_floor"}.
#' @seealso \code{\link{beta_marginal_safe_set}}, \code{\link{floor_coupling_eigenvalues}}
#' @export
group_precision_floor <- function(design = NULL,
                                pfamily_list = NULL,
                                family = NULL,
                                beta_mode = NULL,
                                beta_set = NULL,
                                epsilon = 0.01,
                                dispprior_list = NULL,
                                offset = NULL,
                                weights = NULL,
                                kappa_method = c("crude", "laplace", "none"),
                                union_over = c("groups", "functionals"),
                                mode_method = c("marginal_newton", "icm"),
                                level_method = c("prop2_groups", "prop2_functionals", "r_gauss_joint"),
                                tol = 1e-10,
                                maxit = 200L,
                                verbose = FALSE) {
  kappa_method <- match.arg(kappa_method)
  union_over <- match.arg(union_over)
  mode_method <- match.arg(mode_method)
  level_method <- match.arg(level_method)

  if (!is.null(beta_set) && !inherits(beta_set, "beta_marginal_safe_set")) {
    stop("'beta_set' must be from beta_marginal_safe_set().", call. = FALSE)
  }
  if (!is.null(beta_mode) && !inherits(beta_mode, "beta_marginal_mode")) {
    stop("'beta_mode' must be from beta_marginal_mode().", call. = FALSE)
  }

  if (!is.null(beta_set)) {
    if (is.null(beta_mode)) {
      beta_mode <- beta_set$beta_mode
    }
    level <- beta_set$level
    level_method <- "r_gauss_joint"
    epsilon <- beta_set$delta_2
  }

  if (!is.null(beta_mode)) {
    prep <- beta_mode$prep
    if (is.null(prep)) {
      stop("'beta_mode' is missing 'prep'; recompute with beta_marginal_mode().",
           call. = FALSE)
    }
    engine <- beta_mode$engine
    mode_info <- list(
      beta = beta_mode$beta,
      hessian = beta_mode$hessian,
      method = beta_mode$method,
      icm = NULL
    )
  } else {
    if (is.null(design) || is.null(pfamily_list) || is.null(family)) {
      stop(
        "'design', 'pfamily_list', and 'family' are required when ",
        "'beta_mode' is NULL.",
        call. = FALSE
      )
    }
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
  }

  if (is.null(beta_set)) {
    if (identical(level_method, "prop2_groups") ||
        identical(level_method, "prop2_functionals")) {
      # legacy default path
    } else if (identical(level_method, "r_gauss_joint")) {
      # ok with epsilon as delta_2
    }
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
  } else {
    kappa <- rep(0, prep$J)
  }

  if (is.null(beta_set) &&
      identical(level_method, "prop2_groups") &&
      identical(kappa_method, "crude") &&
      any(level$R > level$r_joint + 1e-8)) {
    if (isTRUE(verbose)) {
      message(
        "Crude kappa inflates R_j above the joint level r(n): ",
        "max(R)/r_joint = ",
        signif(max(level$R) / level$r_joint, 4),
        "."
      )
    }
  }

  uo <- if (identical(level$scheme, "functionals")) {
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
      beta_mode = beta_mode,
      beta_set = beta_set,
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
      mode_method = if (!is.null(beta_mode)) "marginal_newton" else mode_method,
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
