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
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
      stop(
        fn_name, "(): dGamma() dispersion_ranef prior_list requires positive ",
        "'shape' and 'rate'.",
        call. = FALSE
      )
    }
    return(list(
      mode                  = "gamma",
      dispersion_fix        = shape / rate,
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
      ## calibration only; ICM plug-in tau^2 comes from the pfamily spec.
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

    tau2_k <- .two_block_tau2_plug_in_from_pfamily(pf)
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

#' Resolve Block~1 \eqn{\sigma^2} prior weight into observation-scale \code{n_prior}.
#'
#' \eqn{n_{\mathrm{prior}} = w/(1-w)\times n} with \eqn{w =} \code{pwt_measurement}.
#' Independent of Block~2 fixef \code{pwt} and Block~2 \eqn{\tau^2}
#' \code{pwt_dispersion}.
#' @noRd
.lmebayes_resolve_measurement_disp_prior <- function(
    pwt_measurement,
    n_prior_measurement,
    n_obs
) {
  if (!is.null(pwt_measurement) && !is.null(n_prior_measurement)) {
    stop(
      "Supply at most one of 'pwt_measurement' and 'n_prior_measurement'.",
      call. = FALSE
    )
  }
  if (!is.numeric(n_obs) || length(n_obs) != 1L || !is.finite(n_obs) ||
      n_obs <= 0) {
    stop("'n_obs' must be a positive finite scalar.", call. = FALSE)
  }

  if (!is.null(pwt_measurement)) {
    if (!is.numeric(pwt_measurement) || length(pwt_measurement) != 1L ||
        is.na(pwt_measurement) || pwt_measurement <= 0 || pwt_measurement >= 1) {
      stop(
        "'pwt_measurement' must be a scalar in (0, 1).",
        call. = FALSE
      )
    }
    w <- as.numeric(pwt_measurement)
    n_prior <- w / (1 - w) * n_obs
    src <- "user-supplied (pwt_measurement)"
  } else if (!is.null(n_prior_measurement)) {
    if (!is.numeric(n_prior_measurement) || length(n_prior_measurement) != 1L ||
        is.na(n_prior_measurement) || n_prior_measurement <= 0 ||
        !is.finite(n_prior_measurement)) {
      stop(
        "'n_prior_measurement' must be a positive finite scalar.",
        call. = FALSE
      )
    }
    n_prior <- as.numeric(n_prior_measurement)
    w <- n_prior / (n_prior + n_obs)
    src <- "user-supplied (n_prior_measurement)"
  } else {
    w <- 0.01
    n_prior <- w / (1 - w) * n_obs
    src <- "default (pwt_measurement = 0.01)"
  }

  if (n_prior > n_obs) {
    stop(
      "Measurement dispersion prior requires n_prior <= n (equivalently ",
      "pwt_measurement <= 0.5); got n_prior = ", signif(n_prior, 4),
      ", n = ", n_obs, ".",
      call. = FALSE
    )
  }

  list(
    pwt_measurement     = w,
    n_prior_measurement = n_prior,
    source              = src
  )
}

#' Central 98% prior-mass \eqn{\sigma^2}/\eqn{\tau^2} window from calibrated precision prior
#'
#' Precision \eqn{1/\sigma^2 \sim \mathrm{Gamma}(\code{shape}, \code{rate})};
#' bounds are 0.01/0.99 quantiles inverted to the variance scale.
#' @noRd
.lmebayes_ing_prior_quantile_window <- function(shape, rate) {
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      "ING prior quantile window requires positive finite shape and rate.",
      call. = FALSE
    )
  }
  list(
    disp_lower = 1 / stats::qgamma(0.99, shape = shape, rate = rate),
    disp_upper = 1 / stats::qgamma(0.01, shape = shape, rate = rate)
  )
}

#' Prospective \code{dGamma()} measurement \eqn{\sigma^2} calibration from setup
#'
#' Mean-matched inverse-Gamma hyperparameters for Block~1 ING (same algebra as
#' \code{ing_prior} for \eqn{\tau^2_k}, with \eqn{\hat\sigma^2} =
#' \code{dispersion_ranef}, \eqn{p = p_{\mathrm{re}}}, and
#' \eqn{n_{\mathrm{prior}} = \mathrm{pwt\_measurement}/(1-\mathrm{pwt\_measurement})\times n} on the total
#' observation count).  Truncation bounds are the central 98% prior-mass
#' interval from the same \code{shape}/\code{rate}.
#' @noRd
.lmebayes_calibrate_ing_prior_measurement <- function(
    design,
    dispersion_ranef,
    n_prior
) {
  p_re <- length(design$re_coef_names)
  n    <- length(design$y)
  if (p_re < 1L) {
    stop(
      "Measurement dispersion calibration requires at least one random coefficient.",
      call. = FALSE
    )
  }

  if (!is.numeric(n_prior) || length(n_prior) != 1L || !is.finite(n_prior) ||
      n_prior <= 0) {
    stop(
      "'n_prior' must be a positive finite scalar for measurement dispersion calibration.",
      call. = FALSE
    )
  }
  if (n_prior > n) {
    stop(
      "Measurement dispersion prior requires n_prior <= n; got n_prior = ",
      signif(n_prior, 4), ", n = ", n, ".",
      call. = FALSE
    )
  }

  shape <- (n_prior + 1) / 2 + p_re / 2
  rate  <- dispersion_ranef * (n_prior + p_re - 1) / 2
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      "Measurement dispersion ING calibration produced non-positive shape/rate.",
      call. = FALSE
    )
  }

  win <- .lmebayes_ing_prior_quantile_window(shape, rate)

  list(
    sigma2_hat  = dispersion_ranef,
    shape       = shape,
    rate        = rate,
    disp_lower  = win$disp_lower,
    disp_upper  = win$disp_upper,
    n_prior     = n_prior,
    n_effective = n,
    p_re        = p_re
  )
}

#' Limiting-posterior \eqn{\sigma^2}/\eqn{\tau^2} truncation window (lmebayes default)
#'
#' Central 98% mass of \code{Gamma((J+1)/2, d_hat*(J-1)/2)} inverted to the
#' variance scale; see \code{inst/ING_TRUNCATION_WINDOW.md} in \pkg{lmebayes}.
#' @noRd
.lmebayes_ing_limiting_posterior_window <- function(d_hat, J) {
  if (!is.numeric(d_hat) || length(d_hat) != 1L || !is.finite(d_hat) ||
      d_hat <= 0) {
    stop(
      "'d_hat' must be a positive finite scalar (classical variance plug-in).",
      call. = FALSE
    )
  }
  J <- as.integer(J[1L])
  if (!is.finite(J) || J < 1L) {
    stop("'J' must be a positive integer (number of groups).", call. = FALSE)
  }
  a_inf <- (J + 1) / 2
  b_inf <- as.numeric(d_hat) * (J - 1) / 2
  if (b_inf <= 0) {
    stop(
      "Limiting-posterior ING window requires J >= 2 (got J = ", J, ").",
      call. = FALSE
    )
  }
  list(
    disp_lower = 1 / stats::qgamma(0.99, shape = a_inf, rate = b_inf),
    disp_upper = 1 / stats::qgamma(0.01, shape = a_inf, rate = b_inf)
  )
}

#' Build shared ING Block~1 measurement \code{prior_list} for lmebayes glue
#' @noRd
.lmebayes_ing_measurement_prior_list <- function(prior, disp_info, design) {
  re_names <- design$re_coef_names
  p_re     <- length(re_names)
  pl       <- disp_info$dispersion_prior_list
  if (is.null(pl$shape) || is.null(pl$rate)) {
    stop(
      "dGamma() dispersion_ranef prior_list must contain 'shape' and 'rate'.",
      call. = FALSE
    )
  }
  mu <- matrix(
    0,
    nrow = p_re,
    ncol = 1L,
    dimnames = list(re_names, NULL)
  )
  Sigma <- as.matrix(prior$Sigma_ranef)
  if (nrow(Sigma) != p_re || ncol(Sigma) != p_re) {
    stop(
      "prior$Sigma_ranef must be ", p_re, " x ", p_re, ".",
      call. = FALSE
    )
  }
  out <- list(
    mu            = mu,
    Sigma         = Sigma,
    shape         = pl$shape,
    rate          = pl$rate,
    max_disp_perc = if (!is.null(pl$max_disp_perc)) pl$max_disp_perc else 0.99
  )
  if (!is.null(pl$disp_lower)) out$disp_lower <- pl$disp_lower
  if (!is.null(pl$disp_upper)) out$disp_upper <- pl$disp_upper
  if (is.null(out$disp_lower) || is.null(out$disp_upper)) {
    win <- .lmebayes_ing_prior_quantile_window(
      shape = as.numeric(pl$shape[1L]),
      rate  = as.numeric(pl$rate[1L])
    )
    if (is.null(out$disp_lower)) out$disp_lower <- win$disp_lower
    if (is.null(out$disp_upper)) out$disp_upper <- win$disp_upper
  }
  if (out$disp_upper <= out$disp_lower) {
    stop(
      "dGamma() measurement prior: implied disp_upper must exceed disp_lower.",
      call. = FALSE
    )
  }
  out
}

#' Shared matrix-level arguments for LMM reg routes
#' @noRd
.lmebayes_matrix_args_lmm <- function(
    n,
    design,
    prior,
    disp_info,
    tv_tol        = 0.01,
    progbar       = TRUE,
    verbose       = FALSE,
    gap_tol       = 0.0196,
    mode_gap_max  = 1.0,
    diag_sweeps   = FALSE
) {
  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  P            <- solve(prior$Sigma_ranef)

  args <- list(
    n             = n,
    y             = design$y,
    x             = design$Z,
    block         = design$groups,
    x_hyper       = design$X_hyper,
    P             = P,
    pfamily_list  = prior$pfamily_list,
    tv_tol        = tv_tol,
    re_coef_names = re_names,
    group_levels  = group_levels,
    group_name    = design$group_name,
    progbar       = progbar,
    verbose       = verbose
  )

  if (identical(disp_info$mode, "gamma")) {
    args$prior_list <- .lmebayes_ing_measurement_prior_list(
      prior     = prior,
      disp_info = disp_info,
      design    = design
    )
  } else {
    args$prior_list <- list(dispersion = disp_info$dispersion_fix)
  }

  if (isTRUE(prior$any_non_normal)) {
    args$gap_tol       <- gap_tol
    args$mode_gap_max  <- mode_gap_max
    args$diag_sweeps   <- diag_sweeps
    args$stage_verbose <- verbose
  }

  args
}

#' Shared matrix-level arguments for GLMM reg routes
#' @noRd
.lmebayes_matrix_args_glmm <- function(
    n,
    design,
    prior,
    family,
    gap_tol       = 0.0196,
    tv_tol        = 0.01,
    mode_gap_max  = 1.0,
    verbose       = FALSE,
    progbar       = FALSE,
    collect_block1 = TRUE
) {
  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  block1_prior <- .lmebayes_block1_prior_list(prior, dispersion_ranef = NULL)

  list(
    n               = n,
    y               = design$y,
    x               = design$Z,
    block           = design$groups,
    x_hyper         = design$X_hyper,
    prior_list      = block1_prior,
    pfamily_list    = prior$pfamily_list,
    family          = family,
    re_coef_names   = re_names,
    group_levels    = group_levels,
    group_name      = design$group_name,
    gap_tol         = gap_tol,
    tv_tol          = tv_tol,
    mode_gap_max    = mode_gap_max,
    verbose         = verbose,
    progbar         = progbar,
    stage_verbose   = verbose,
    collect_block1  = collect_block1
  )
}

#' @noRd
.lmebayes_run_lmm_engine <- function(
    n,
    design,
    prior,
    disp_info,
    tv_tol        = 0.01,
    progbar       = TRUE,
    verbose       = FALSE,
    gap_tol             = 0.0196,
    mode_gap_max        = 1.0,
    diag_sweeps         = FALSE
) {
  route_key <- .lmebayes_reg_route_key(
    family         = gaussian(),
    disp_mode      = disp_info$mode,
    any_non_normal = prior$any_non_normal
  )
  route <- .lmebayes_reg_route_fn(route_key)
  args  <- .lmebayes_matrix_args_lmm(
    n             = n,
    design        = design,
    prior         = prior,
    disp_info     = disp_info,
    tv_tol        = tv_tol,
    progbar       = progbar,
    verbose       = verbose,
    gap_tol       = gap_tol,
    mode_gap_max  = mode_gap_max,
    diag_sweeps   = diag_sweeps
  )
  out <- do.call(route$export_fn, args)
  .lmebayes_attach_sigma2(out, disp_info)
}

#' @noRd
.lmebayes_run_glmm_engine <- function(
    n,
    design,
    prior,
    family,
    gap_tol       = 0.0196,
    tv_tol        = 0.01,
    mode_gap_max  = 1.0,
    verbose       = FALSE,
    progbar       = FALSE,
    collect_block1 = TRUE
) {
  route_key <- .lmebayes_reg_route_key(
    family         = family,
    disp_mode      = "none",
    any_non_normal = prior$any_non_normal
  )
  route <- .lmebayes_reg_route_fn(route_key)
  args  <- .lmebayes_matrix_args_glmm(
    n              = n,
    design         = design,
    prior          = prior,
    family         = family,
    gap_tol        = gap_tol,
    tv_tol         = tv_tol,
    mode_gap_max   = mode_gap_max,
    verbose        = verbose,
    progbar        = progbar,
    collect_block1 = collect_block1
  )
  out <- do.call(route$export_fn, args)
  disp_none <- list(mode = "none")
  .lmebayes_attach_sigma2(out, disp_none)
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

#' Attach \code{sigma2} / \code{sigma2.mean} from dispersion mode and sampler draws.
#'
#' Fixed measurement dispersion returns a scalar; \code{dGamma()} returns the
#' length-\code{n} vector from the final inner sweep (\code{dispersion_ranef});
#' families without observation-level dispersion get \code{NULL}.
#' @noRd
.lmebayes_attach_sigma2 <- function(out, disp_info) {
  mode <- disp_info$mode
  if (identical(mode, "none")) {
    out$sigma2 <- NULL
    out$sigma2.mean <- NULL
    return(out)
  }
  if (identical(mode, "fixed")) {
    val <- as.numeric(disp_info$dispersion_fix)
    out$sigma2 <- val
    out$sigma2.mean <- val
    return(out)
  }
  if (identical(mode, "gamma")) {
    dr <- out$dispersion_ranef
    if (is.null(dr)) {
      stop(
        "Internal error: dGamma measurement dispersion requires ",
        "'dispersion_ranef' draws on the sampler output.",
        call. = FALSE
      )
    }
    out$sigma2 <- as.numeric(dr)
    out$sigma2.mean <- mean(out$sigma2)
    return(out)
  }
  stop("Unknown dispersion mode: ", mode, call. = FALSE)
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

#' Per-group full column-rank flag for Block~1 \code{Z_j} (same rule as
#' \code{model_setup()$re_rank}).
#' @noRd
.lmebayes_re_rank_from_Z <- function(Z, groups, group_levels = NULL) {
  Z <- as.matrix(Z)
  g_chr <- as.character(groups)
  levs <- if (is.null(group_levels)) {
    unique(g_chr)
  } else {
    as.character(group_levels)
  }
  p_re <- ncol(Z)
  stats::setNames(
    vapply(
      levs,
      function(lev) {
        rows <- which(g_chr == lev)
        Z_j  <- Z[rows, , drop = FALSE]
        nrow(Z_j) >= p_re &&
          Matrix::rankMatrix(Z_j, method = "qr")[1L] == p_re
      },
      logical(1L)
    ),
    levs
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
