## Shared preparation for group_precision_floor() and beta_marginal_safe_set().

#' @noRd
.group_floor_prepare <- function(design,
                                 pfamily_list,
                                 family,
                                 dispprior_list = NULL,
                                 offset = NULL,
                                 weights = NULL,
                                 fn_name = "group_precision_floor") {
  .lmerb_validate_design(design)
  if (is.null(design$y) || is.null(design$D)) {
    stop("'design' must contain 'y' and 'D'.", call. = FALSE)
  }

  family_hook <- .group_floor_family(family)
  group_levels <- levels(design$group)
  J <- length(group_levels)
  re_names <- design$groupef.names
  p_re <- length(re_names)

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, re_names, J = J
  )

  prior_pack <- .priors_from_pfamily_list(
    pfamily_list = pfamily_list,
    group.dispersion = .lmebayes_dispprior_list_as_group_dispersion(dispprior_list),
    design = design,
    family = family_hook$family,
    fn_name = fn_name
  )

  if (isTRUE(prior_pack$any_non_normal)) {
    stop(
      fn_name, "() requires fixed dNormal tau^2 in 'pfamily_list' ",
      "(dIndependent_Normal_Gamma is not supported in v1).",
      call. = FALSE
    )
  }

  measurement_prior_list <- list(
    group.Sigma = prior_pack$group.Sigma,
    pop.prior_list = prior_pack$pop.prior_list
  )
  if (!is.null(prior_pack$group.dispersion)) {
    measurement_prior_list$group.dispersion <- prior_pack$group.dispersion
  }

  dispersion <- if (identical(family_hook$name, "gaussian")) {
    as.numeric(prior_pack$group.dispersion)[1L]
  } else {
    1
  }

  prior_int <- .group_floor_integrate_gamma(
    design = design,
    prior_pack = prior_pack,
    group_levels = group_levels
  )

  group_data <- .group_floor_prepare_group_data(
    design = design,
    family_hook = family_hook,
    group_levels = group_levels,
    offset = offset,
    weights = weights,
    dispersion = dispersion
  )

  list(
    family_hook = family_hook,
    group_levels = group_levels,
    J = J,
    p_re = p_re,
    re_names = re_names,
    prior_int = prior_int,
    group_data = group_data,
    measurement_prior_list = measurement_prior_list,
    dispersion = dispersion,
    n_dim = prior_int$n_dim,
    n_obs_total = sum(group_data$n_obs)
  )
}

#' @noRd
.group_floor_resolve_mode <- function(prep,
                                      design,
                                      mode_method = c("icm", "marginal_newton"),
                                      tol = 1e-10,
                                      maxit = 200L) {
  mode_method <- match.arg(mode_method)
  if (identical(mode_method, "icm")) {
    mode_res <- .group_floor_mode_from_icm(
      design = design,
      family = prep$family_hook$family,
      measurement_prior_list = prep$measurement_prior_list,
      group_levels = prep$group_levels,
      tol = tol,
      maxit = maxit
    )
    list(
      beta = mode_res$beta,
      method = "icm",
      icm = mode_res$icm,
      hessian = NULL
    )
  } else {
    mode_res <- .c05_beta_marginal_mode(
      prior_int = prep$prior_int,
      group_data = prep$group_data,
      family_hook = prep$family_hook
    )
    list(
      beta = mode_res$beta,
      method = "marginal_newton",
      icm = NULL,
      hessian = mode_res$hessian,
      f_mode = mode_res$f_mode
    )
  }
}

#' @noRd
.group_floor_resolve_levels <- function(epsilon,
                                        level_method = c("prop2_groups", "prop2_functionals", "r_gauss_joint"),
                                        union_over = c("groups", "functionals"),
                                        J,
                                        p_re,
                                        n_obs_total,
                                        kappa,
                                        gaussian = FALSE,
                                        kappa_method = "crude",
                                        inflate_kappa = TRUE) {
  level_method <- match.arg(level_method)
  union_over <- match.arg(union_over)

  if (identical(level_method, "r_gauss_joint")) {
    rg <- .c05_beta_r_gauss_level(J * p_re, epsilon)
    R <- rep(rg$r, J)
    return(list(
      epsilon = epsilon,
      scheme = "r_gauss_joint",
      level_method = level_method,
      r_joint = rg$r,
      r_group = rg$r,
      r_functional = NA_real_,
      r_gauss = rg$r,
      r2 = rg$r2,
      laplace_tail_mass = rg$laplace_tail_mass,
      kappa = kappa,
      R = R,
      s = sqrt(2 * R)
    ))
  }

  if (identical(level_method, "prop2_functionals")) {
    union_over <- "functionals"
  } else {
    union_over <- "groups"
  }

  lv <- .group_floor_levels(
    epsilon = epsilon,
    union_over = union_over,
    J = J,
    p_re = p_re,
    n_obs_total = n_obs_total,
    kappa = if (isTRUE(inflate_kappa)) kappa else rep(0, J),
    gaussian = gaussian,
    kappa_method = if (isTRUE(inflate_kappa)) kappa_method else "none"
  )
  lv$level_method <- level_method
  lv$r_gauss <- .c05_beta_r_gauss_level(J * p_re, epsilon)$r
  lv$laplace_tail_mass <- epsilon
  lv
}

#' @noRd
.group_floor_compute_kappa <- function(prep,
                                       engine,
                                       kappa_method,
                                       epsilon,
                                       level_method,
                                       r_for_laplace = NULL) {
  J <- prep$J
  p_re <- prep$p_re
  kappa <- rep(0, J)
  if (prep$family_hook$gaussian || identical(kappa_method, "none")) {
    return(kappa)
  }
  if (identical(kappa_method, "crude")) {
    Gbar <- .group_floor_build_gbar(prep$group_data, prep$family_hook)
    return(vapply(seq_len(J), function(j) {
      .group_floor_kappa_crude(j, prep$prior_int$Lambda_beta, Gbar, p_re, J)
    }, numeric(1L)))
  }
  if (identical(kappa_method, "laplace")) {
    if (is.null(r_for_laplace)) {
      r_for_laplace <- if (identical(level_method, "r_gauss_joint")) {
        .c05_beta_r_gauss_level(J * p_re, epsilon)$r
      } else {
        .gamma_level_for_budget(p_re, epsilon / J)
      }
    }
    return(vapply(seq_len(J), function(j) {
      cv <- engine$emb(j, prep$group_data$Z_list[[j]][1, ])
      .group_floor_kappa_laplace(j, engine, cv, r_for_laplace)
    }, numeric(1L)))
  }
  kappa
}

#' @noRd
.group_floor_assemble_floors <- function(prep,
                                         engine,
                                         level,
                                         union_over = c("groups", "functionals")) {
  union_over <- match.arg(union_over)
  J <- prep$J
  group_levels <- prep$group_levels
  per_group <- vector("list", J)
  Gamma_lb <- vector("list", J)
  names(Gamma_lb) <- group_levels
  names(per_group) <- group_levels

  for (j in seq_len(J)) {
    r_j <- if (identical(union_over, "functionals")) {
      rep(level$r_functional, prep$group_data$n_obs[j])
    } else {
      level$R[j]
    }
    fg <- .group_floor_one_group(
      j = j,
      engine = engine,
      group_data = prep$group_data,
      family_hook = prep$family_hook,
      r_level = r_j
    )
    Gamma_lb[[j]] <- fg$Gamma
    per_group[[j]] <- list(
      group = group_levels[j],
      Gamma = fg$Gamma,
      eta_range = cbind(lo = fg$lo, hi = fg$hi),
      es = fg$es,
      wbar = fg$wbar,
      omega = fg$omega
    )
  }

  list(Gamma_lb = Gamma_lb, per_group = per_group)
}
