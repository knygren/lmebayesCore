## rglmerb(..., simulate = FALSE) point estimates (glmerb only).

#' Point estimates for Bayesian GLMMs (\code{rglmerb(..., simulate = FALSE)})
#'
#' Computes the joint posterior mode (ICM) at fixed variance-component
#' plug-ins, or the exact Gaussian mean when \code{family = gaussian()},
#' without MCMC draws. Called from \code{\link{rglmerb}} when
#' \code{simulate = FALSE}; not exported.
#'
#' @inheritParams rglmerb
#' @param tol,maxit ICM controls passed to \code{\link{glmerb_posterior_mode}}
#'   (unused for the closed-form Gaussian mean).
#' @param offset_missing,weights_missing See \code{\link{rlmerb_point}}.
#' @return An object of class \code{c("rglmerb", "list")} with mode fields,
#'   full \code{prior}, thin \code{Prior}, and \code{NULL} draw slots.
#' @seealso \code{\link{rglmerb}}, \code{\link{glmerb_posterior_mode}},
#'   \code{\link{rlmerb_point}}
#' @keywords internal
rglmerb_point <- function(
    design,
    pfamily_list,
    family = stats::poisson(),
    dispprior_list = NULL,
    offset = NULL,
    weights = 1,
    verbose = TRUE,
    print_icm_table = TRUE,
    tol = 1e-10,
    maxit = 200L,
    offset_missing = NULL,
    weights_missing = NULL
) {
  cl <- match.call()
  if (is.null(offset_missing)) offset_missing <- missing(offset)
  if (is.null(weights_missing)) weights_missing <- missing(weights)

  if (!inherits(family, "family") || is.null(family$family)) {
    stop("'family' must be a family object.", call. = FALSE)
  }
  prep <- .rlmerb_prepare_prior(
    design          = design,
    pfamily_list    = pfamily_list,
    dispprior_list  = dispprior_list,
    family          = family,
    fn_name         = "rglmerb",
    offset          = offset,
    weights         = weights,
    offset_missing  = offset_missing,
    weights_missing = weights_missing
  )
  prior <- prep$prior
  disp_info <- prep$disp_info
  re_names <- prep$re_names

  icm <- .rlmerb_icm_at_fixed_vc(
    design = design,
    prior  = prior,
    family = family,
    tol    = tol,
    maxit  = maxit
  )

  if (isTRUE(print_icm_table)) {
    icm_lbl <- .lmebayes_block2_icm_labels(prior, family)
    .lmebayes_print_icm_fixef_table(
      prior_list = prior$pop.prior_list,
      re_names   = re_names,
      fixef_icm  = icm$fixef,
      icm_info   = icm$icm_info,
      ref_label  = icm_lbl$ref_label,
      icm_label  = icm$icm_label,
      conv_label = icm_lbl$conv_label,
      header     = "--- rglmerb: population effects ---",
      verbose    = verbose
    )
  }

  n_obs <- length(design$y)
  off <- prep$ow$offset
  wts <- prep$ow$weights
  offset2 <- if (is.null(off)) {
    rep.int(0, n_obs)
  } else if (length(off) == 1L) {
    rep.int(as.numeric(off), n_obs)
  } else {
    as.numeric(off)
  }
  prior.weights <- if (length(wts) == 1L) {
    rep.int(as.numeric(wts), n_obs)
  } else {
    as.numeric(wts)
  }

  out <- list(
    call                  = cl,
    family                = family,
    popef.mode            = icm$fixef,
    popef.init            = icm$fixef_init,
    groupef.mode          = icm$b_mean,
    popef.means           = NULL,
    popef                 = NULL,
    groupef               = NULL,
    popef.dispersion      = NULL,
    popef.dispersion.mean = NULL,
    popef.iters           = NULL,
    popef.iters.mean      = NULL,
    groupef.iters         = NULL,
    groupef.iters.mean    = NULL,
    group.dispersion      = prior$group.dispersion,
    group.dispersion.mean = prior$group.dispersion,
    group.dispersion.mode = prior$group.dispersion,
    group.dispersion.iters = NULL,
    group.dispersion.iters.mean = NULL,
    tau2.mode             = icm$tau2,
    joint_mode            = icm$joint_mode,
    icm_info              = icm$icm_info,
    m_convergence         = NULL,
    convergence           = NULL,
    convergence_info      = NULL,
    pilot                 = NULL,
    sweep_history         = NULL,
    prior                 = prior,
    Prior                 = .rlmerb_thin_Prior(prior, disp_info, prep$block1_prior),
    design                = design,
    offset                = off,
    offset2               = offset2,
    prior.weights         = prior.weights
  )
  class(out) <- c("rglmerb", "list")
  out
}
