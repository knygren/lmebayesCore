## Internal helpers for matrix-level mixed-model samplers (rlmerb / rglmerb).
## lmebayes calls these via glmbayesCore::: from lmerb() / glmerb() only.

#' @noRd
.lmebayes_resolve_dispersion_ranef <- function(
    dispersion_ranef,
    family,
    design = NULL,
    fn_name = "lmerb"
) {
  has_dispersion <- family$family %in%
    c("gaussian", "Gamma", "quasipoisson", "quasibinomial")

  if (!has_dispersion) {
    if (!is.null(dispersion_ranef)) {
      stop(
        "'dispersion_ranef' must be NULL for family = ", family$family,
        "() (no observation-level dispersion).",
        call. = FALSE
      )
    }
    return(list(
      mode                   = "none",
      dispersion_fix         = NULL,
      dispersion_prior_list  = NULL,
      dispersion_pfamily     = NULL
    ))
  }

  if (inherits(dispersion_ranef, "pfamily")) {
    if (!identical(dispersion_ranef$pfamily, "dGamma")) {
      stop(
        fn_name, "(): 'dispersion_ranef' pfamily must be dGamma(); got ",
        dispersion_ranef$pfamily, ". RE priors belong in 'pfamily_list'.",
        call. = FALSE
      )
    }
    pl <- dispersion_ranef$prior_list
    if (!isTRUE(pl$Inv_Dispersion)) {
      stop(
        fn_name, "(): dGamma() observation-dispersion prior requires ",
        "Inv_Dispersion = TRUE.",
        call. = FALSE
      )
    }
    if (is.null(design) || is.null(design$residual_var)) {
      stop(
        fn_name, "(): a model_setup with residual_var is required for ",
        "dGamma() dispersion_ranef (plug-in sigma^2).",
        call. = FALSE
      )
    }
    return(list(
      mode                  = "gamma",
      dispersion_fix        = as.numeric(design$residual_var),
      dispersion_prior_list = pl,
      dispersion_pfamily    = dispersion_ranef
    ))
  }

  if (is.null(dispersion_ranef) || !is.numeric(dispersion_ranef) ||
      length(dispersion_ranef) != 1L || !is.finite(dispersion_ranef) ||
      dispersion_ranef <= 0) {
    stop(
      "'dispersion_ranef' must be a single positive number or a dGamma() ",
      "pfamily for family = ", family$family, "().",
      call. = FALSE
    )
  }
  list(
    mode                  = "fixed",
    dispersion_fix        = as.numeric(dispersion_ranef),
    dispersion_prior_list = NULL,
    dispersion_pfamily    = NULL
  )
}

#' @noRd
.lmebayes_validate_dispersion_ranef <- function(
    dispersion_ranef,
    family,
    fn_name = "lmerb"
) {
  resolved <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = NULL,
    fn_name          = fn_name
  )
  resolved$dispersion_fix
}

#' @noRd
.lmebayes_priors_from_pfamily_list <- function(pfamily_list,
                                               dispersion_ranef,
                                               design,
                                               family,
                                               fn_name = "lmerb") {

  re_names <- design$re_coef_names
  p_re     <- length(re_names)

  ## --- dispersion_ranef (Block 1 measurement dispersion) -------------------
  disp_res <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = design,
    fn_name          = fn_name
  )
  dispersion_ranef <- disp_res$dispersion_fix

  ## --- pfamily_list ---------------------------------------------------------
  if (!is.list(pfamily_list) || length(pfamily_list) != p_re) {
    stop(
      "'pfamily_list' must be a list with one pfamily per random-effect ",
      "component (", p_re, " expected: ", paste(re_names, collapse = ", "),
      "). Build it with pfamily_list(Prior_Setup_lmebayes(...)).",
      call. = FALSE
    )
  }
  if (is.null(names(pfamily_list)) || !setequal(names(pfamily_list), re_names)) {
    stop(
      "Names of 'pfamily_list' must match the random-effect coefficient ",
      "names: ", paste(re_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  pfamily_list <- pfamily_list[re_names]

  prior_list <- stats::setNames(vector("list", p_re), re_names)
  tau2   <- stats::setNames(numeric(p_re), re_names)
  ptypes <- stats::setNames(character(p_re), re_names)

  for (k in re_names) {
    pf <- pfamily_list[[k]]
    if (!inherits(pf, "pfamily")) {
      stop("pfamily_list[[\"", k, "\"]] must be a pfamily object.",
           call. = FALSE)
    }
    if (!pf$pfamily %in% c("dNormal", "dIndependent_Normal_Gamma")) {
      stop(
        fn_name, "() supports only dNormal and dIndependent_Normal_Gamma ",
        "pfamilies in 'pfamily_list'; component \"", k, "\" is ",
        pf$pfamily, ".",
        call. = FALSE
      )
    }
    ptypes[[k]] <- pf$pfamily

    par_names <- colnames(design$X_hyper[[k]])
    q_k <- length(par_names)

    mu_k <- as.numeric(pf$prior_list$mu)
    if (length(mu_k) != q_k) {
      stop(
        "pfamily_list[[\"", k, "\"]]$prior_list$mu has length ",
        length(mu_k), " but the hyper design has ", q_k, " column(s): ",
        paste(par_names, collapse = ", "), ".",
        call. = FALSE
      )
    }
    mu_nms <- rownames(pf$prior_list$mu)
    if (!is.null(mu_nms) && all(nzchar(mu_nms))) {
      if (!setequal(mu_nms, par_names)) {
        stop(
          "Parameter names of pfamily_list[[\"", k, "\"]] (",
          paste(mu_nms, collapse = ", "), ") do not match the hyper design ",
          "columns (", paste(par_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      ord <- match(par_names, mu_nms)
      mu_k <- mu_k[ord]
      Sigma_k <- as.matrix(pf$prior_list$Sigma)[ord, ord, drop = FALSE]
    } else {
      Sigma_k <- as.matrix(pf$prior_list$Sigma)
    }
    names(mu_k) <- par_names
    dimnames(Sigma_k) <- list(par_names, par_names)

    ## Keep the pfamily object itself aligned with the hyper-design column
    ## order: it is passed straight to the v2 sampler as the Block 2 source
    ## of truth, so its mu/Sigma must match x_hyper[[k]].
    pfamily_list[[k]]$prior_list$mu <-
      matrix(mu_k, ncol = 1L, dimnames = list(par_names, NULL))
    pfamily_list[[k]]$prior_list$Sigma <- Sigma_k

    if (identical(pf$pfamily, "dNormal")) {
      d_k <- pf$prior_list$dispersion
      if (isTRUE(pf$prior_list$ddef)) {
        warning(
          fn_name, ": pfamily_list[[\"", k, "\"]] uses the default ",
          "dispersion = 1 (none was supplied to dNormal()); the Block 1 ",
          "random-effect variance tau^2 for \"", k, "\" is therefore 1.",
          call. = FALSE
        )
      }
    } else {
      ## ING: disp_lower/disp_upper fix the truncation window and lambda*
      ## calibration; plug-in tau^2 for Sigma_ranef is rate/shape = 1/E[1/tau^2].
      d_k <- pf$prior_list$disp_lower
      if (is.null(d_k) || !is.numeric(d_k) || length(d_k) != 1L ||
          !is.finite(d_k) || d_k <= 0) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a positive scalar ",
          "'disp_lower' (lower dispersion truncation) for lambda* calibration.",
          call. = FALSE
        )
      }
      u_k <- pf$prior_list$disp_upper
      if (is.null(u_k) || !is.numeric(u_k) || length(u_k) != 1L ||
          !is.finite(u_k) || u_k <= as.numeric(d_k)) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a finite scalar ",
          "'disp_upper' > 'disp_lower' (upper dispersion truncation), so ",
          "the tau^2 truncation window is fixed across Gibbs sweeps. ",
          "pfamily_list(Prior_Setup_lmebayes(...)) sets both bounds to the ",
          "0.01/0.99 prior dispersion quantiles by default.",
          call. = FALSE
        )
      }
    }

    tau2_k <- .two_block_tau2_ref_from_pfamily(pf)
    tau2[[k]] <- tau2_k
    prior_list[[k]] <- list(
      mu_fixef         = mu_k,
      Sigma_fixef      = Sigma_k,
      dispersion_fixef = tau2_k
    )
  }

  Sigma_ranef <- diag(unname(tau2), nrow = p_re, ncol = p_re)
  dimnames(Sigma_ranef) <- list(re_names, re_names)

  list(
    pfamily_list          = pfamily_list,
    dispersion_ranef      = dispersion_ranef,
    dispersion_mode       = disp_res$mode,
    dispersion_pfamily    = disp_res$dispersion_pfamily,
    dispersion_prior_list = disp_res$dispersion_prior_list,
    Sigma_ranef           = Sigma_ranef,
    prior_list            = prior_list,
    ptypes         = ptypes,
    any_non_normal = any(ptypes != "dNormal")
  )
}

#' @noRd
.lmebayes_run_lmm_engine <- function(
    n,
    design,
    prior,
    disp_info,
    fixef_start   = NULL,
    m_convergence = NULL,
    tv_tol        = 0.01,
    progbar       = TRUE,
    verbose       = FALSE,
    gap_tol             = 0.0196,
    mode_gap_max        = 1.0,
    diag_sweeps         = FALSE
) {
  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  P            <- solve(prior$Sigma_ranef)
  common_args  <- list(
    n             = n,
    y             = design$y,
    x             = design$Z,
    block         = design$groups,
    x_hyper       = design$X_hyper,
    P             = P,
    pfamily_list  = prior$pfamily_list,
    start         = fixef_start,
    m_convergence = m_convergence,
    tv_tol        = tv_tol,
    re_coef_names = re_names,
    group_levels  = group_levels,
    group_name    = design$group_name,
    progbar       = progbar,
    verbose       = verbose
  )
  if (identical(disp_info$mode, "gamma")) {
    do.call(
      rLMMindepNormalGamma_reg,
      c(
        common_args,
        list(
          prior_list     = disp_info$dispersion_prior_list,
          dispersion_fix = disp_info$dispersion_fix
        )
      )
    )
  } else if (isTRUE(prior$any_non_normal)) {
    do.call(
      rLMMNormal_reg_estimated_vcov,
      c(
        common_args,
        list(
          prior_list    = list(dispersion = disp_info$dispersion_fix),
          gap_tol       = gap_tol,
          mode_gap_max  = mode_gap_max,
          diag_sweeps   = diag_sweeps,
          stage_verbose = verbose
        )
      )
    )
  } else {
    do.call(
      rLMMNormal_reg,
      c(
        common_args,
        list(prior_list = list(dispersion = disp_info$dispersion_fix))
      )
    )
  }
}

#' @noRd
.lmebayes_block1_prior_list <- function(
    measurement_prior_list,
    dispersion_ranef = NULL
) {
  if (is.null(measurement_prior_list$Sigma_ranef)) {
    stop("measurement_prior_list must contain 'Sigma_ranef'.", call. = FALSE)
  }
  P <- solve(measurement_prior_list$Sigma_ranef)
  dispersion <- if (!is.null(dispersion_ranef)) {
    dispersion_ranef
  } else {
    measurement_prior_list$dispersion_ranef
  }
  if (is.null(dispersion)) {
    list(P = P, ddef = TRUE)
  } else {
    list(P = P, dispersion = dispersion, ddef = FALSE)
  }
}

#' @noRd
.lmebayes_add_fixef_summaries <- function(x) {
  if (!is.null(x$fixef)) {
    x$fixef.means <- lapply(x$fixef, colMeans)
  }
  if (!is.null(x$fixef.dispersion)) {
    x$fixef.dispersion.mean <- colMeans(x$fixef.dispersion)
  }
  if (!is.null(x$fixef.iters) && !is.null(x$m_convergence)) {
    x$fixef.iters.mean <- colMeans(x$fixef.iters) / x$m_convergence
  }
  if (!is.null(x$ranef.iters) && !is.null(x$m_convergence)) {
    x$ranef.iters.mean <- mean(x$ranef.iters) / x$m_convergence
  }
  x
}

#' @noRd
.lmebayes_block2_icm_labels <- function(prior, family = gaussian()) {
  any_ing <- isTRUE(prior$any_non_normal)
  is_gauss <- is.null(family) || identical(family$family, "gaussian")
  ref_label <- "prior mean"
  if (any_ing) {
    icm_label   <- "gamma @ lmer tau2"
    icm_verbose <- "Block 2 start at lmer tau^2 plug-in"
    conv_label  <- "Plug-in fixed point"
  } else if (is_gauss) {
    icm_label   <- "ICM mean"
    icm_verbose <- "ICM posterior mean"
    conv_label  <- "ICM"
  } else {
    icm_label   <- "ICM mode"
    icm_verbose <- "ICM posterior mode"
    conv_label  <- "ICM"
  }
  list(
    ref_label   = ref_label,
    icm_label   = icm_label,
    icm_verbose = icm_verbose,
    conv_label  = conv_label
  )
}

#' @noRd
.lmebayes_print_icm_fixef_table <- function(
    prior_list,
    re_names,
    fixef_icm,
    icm_info,
    ref_label,
    icm_label,
    conv_label = "ICM",
    header,
    verbose
) {
  if (!isTRUE(verbose) || is.null(fixef_icm)) {
    return(invisible(NULL))
  }
  fixef_ref <- lapply(prior_list, `[[`, "mu_fixef")
  names(fixef_ref) <- re_names
  hdr <- sprintf("  %-18s  %-30s  %14s  %18s",
                 "RE component", "parameter", ref_label, icm_label)
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat(header, "\n")
  cat(hdr, "\n")
  cat(sep, "\n")
  for (k in re_names) {
    nms_k  <- names(fixef_ref[[k]])
    ref_v  <- fixef_ref[[k]]
    icm_v  <- fixef_icm[[k]]
    for (nm in nms_k) {
      cat(sprintf("  %-18s  %-30s  %14.4f  %18.4f\n",
                  k, nm, ref_v[[nm]], icm_v[[nm]]))
    }
  }
  if (!is.null(icm_info)) {
    cat(sprintf("  (%s converged: %s, %d iter, delta = %.2e)\n\n",
                conv_label,
                icm_info$converged, icm_info$iterations, icm_info$delta))
  } else {
    cat("\n")
  }
  invisible(NULL)
}

#' @noRd
.lmebayes_print_ranef_mode_reference <- function(
    ranef_mode,
    re_names,
    group_levels,
    verbose
) {
  invisible(NULL)
}

#' @noRd
.lmebayes_print_fixef_init <- function(
    fixef_init,
    re_names,
    verbose,
    header = "--- main-stage fixef.init (pilot colMeans) ---"
) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat(header, "\n")
  for (k in re_names) {
    for (nm in names(fixef_init[[k]])) {
      cat(sprintf("  %-18s  %-30s  %12.4f\n",
                  k, nm, fixef_init[[k]][[nm]]))
    }
  }
  cat("\n")
  invisible(NULL)
}
