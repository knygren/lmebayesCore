#' The Bayesian Generalized Linear Mixed-Effects Model Distribution
#'
#' \code{rglmerb} generates posterior draws for Bayesian generalized linear
#' mixed models from \code{model_setup} design objects and the same prior
#' arguments as \code{\link{rGLMM_reg}} / \code{\link{rLMMNormal_reg}}
#' (\code{pfamily_list}, \code{dispprior_list}). Routes by response family:
#' \itemize{
#'   \item \code{family = gaussian()} delegates to
#'     \code{\link{rLMMNormal_reg_known_vcov}},
#'     \code{\link{rLMMNormal_reg_estimated_vcov}},
#'     \code{\link{rLMMindepNormalGamma_reg_known_vcov}}, or
#'     \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}} according to
#'     \code{dispprior_list} and population \code{pfamily_list}.
#'   \item Non-Gaussian families delegate via \code{\link{rGLMM_reg}}
#'     (pilot stage always unless \code{n_pilot = 0L}; routes differ in
#'     eigenvalue-bound complexity). Inner sweeps use
#'     \code{\link{rGLMM_sweep}}.
#' }
#'
#' For formula interfaces, \code{glmerb()} in the lmebayes package wraps this
#' function. Design matrices are built with \code{\link{model_setup}} and
#' priors with \code{\link{Prior_Setup_GLMM}}.
#'
#' @param n Integer. Number of independent chains in the main stage.
#' @param design A \code{model_setup} object (from \code{\link{model_setup}}).
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects, one per
#'   random-effect coefficient in \code{design$groupef.names} (same as
#'   \code{\link{rGLMM_reg}} / \code{\link{rLMMNormal_reg}}).
#' @param family A \code{\link[stats]{family}} object. Default \code{poisson()}.
#' @param dispprior_list Block~1 observation-dispersion prior, in the same
#'   shapes as the \code{rLMM*} / \code{rGLMM*} engines. Required for
#'   \code{family = gaussian()} (\code{list(dispersion = )} or
#'   \code{dGamma()} / \code{dGamma_list()}); must be \code{NULL} (default)
#'   for \code{poisson()} and \code{binomial()}.
#' @param offset,weights Observation offset and prior weights, as in
#'   \code{\link{rGLMM_reg}} / \code{\link{rLMMNormal_reg}}
#'   (\code{offset = NULL}, \code{weights = 1}). When omitted, inherit
#'   \code{design$offset} / \code{design$weights} if present. Explicit values
#'   always override the design. Echoed on the fit; not yet used in Gibbs
#'   sweeps.
#' @param gap_tol Legacy mode--mean gap for deriving the pilot chain count when
#'   \code{tv_tol} is \code{NULL}. Ignored for Gaussian without ING population
#'   components.
#' @param tv_tol Total variation tolerance for convergence calibration.
#'   Inner Gibbs sweeps and pilot chain counts are derived internally.
#' @param mode_gap_max Pilot inner-sweep calibration (non-Gaussian and
#'   Gaussian+ING only).
#' @param collect_block1 Collect \code{groupef} draws from main chains
#'   (non-Gaussian only).
#' @param verbose Print stage headers and diagnostics.
#' @param progbar Progress bars when \code{verbose} is \code{FALSE}.
#' @param sim_method Simulation method for \code{family = gaussian()}:
#'   \code{"DEFAULT"} or \code{"TWO_BLOCK_GIBBS"}; see \code{\link{rlmerb}}.
#'   Ignored (two-block Gibbs only) for non-Gaussian families.
#' @return Object of class \code{c("rglmerb", "list")} with \code{popef.*},
#'   \code{groupef}/\code{groupef.mode}, and (Gaussian only)
#'   \code{group.dispersion}/\code{group.dispersion.mean}, plus
#'   \code{Prior}, \code{design}, and \code{family}. When a pilot stage runs,
#'   nested \code{pilot} (\code{n}, \code{chisq}, \code{draws}) is included.
#' @seealso \code{\link{rlmerb}}, \code{\link{rLMMNormal_reg_known_vcov}},
#'   \code{\link{rLMMNormal_reg_estimated_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg_known_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}}, \code{\link{rGLMM_reg}},
#'   \code{\link{Prior_Setup_GLMM}}
#' @name rglmerb
NULL

#' @rdname rglmerb
#' @export
rglmerb <- function(
    n,
    design,
    pfamily_list,
    family              = poisson(),
    dispprior_list      = NULL,
    offset              = NULL,
    weights             = 1,
    gap_tol             = 0.0196,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    collect_block1      = TRUE,
    verbose             = TRUE,
    progbar             = FALSE,
    sim_method          = "DEFAULT"
) {
  cl <- match.call()

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  sim_method <- .rLMM_validate_sim_method(sim_method, fn_name = "rglmerb")

  if (!inherits(design, "model_setup")) {
    stop("'design' must be a model_setup object.", call. = FALSE)
  }

  if (!inherits(family, "family") || is.null(family$family)) {
    stop("'family' must be a family object.", call. = FALSE)
  }

  if (missing(pfamily_list) || is.null(pfamily_list)) {
    stop(
      "'pfamily_list' is required for rglmerb() (named Block~2 pfamily list, ",
      "same as rGLMM_reg() / rLMMNormal_reg()).",
      call. = FALSE
    )
  }

  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }

  is_gaussian <- identical(family$family, "gaussian")

  group.dispersion <- .lmebayes_dispprior_list_as_group_dispersion(dispprior_list)
  prior <- priors_from_pfamily_list(
    pfamily_list     = pfamily_list,
    group.dispersion = group.dispersion,
    design           = design,
    family           = family,
    fn_name          = "rglmerb"
  )
  disp_info <- list(
    mode                  = prior$dispersion_mode,
    dispersion_fix        = prior$group.dispersion,
    dispersion_prior_list = prior$dispersion_prior_list,
    dispersion_pfamily    = prior$dispersion_pfamily,
    window_diagnostics    = prior$window_diagnostics
  )

  re_names     <- design$groupef.names
  group_levels <- levels(design$group)

  if (is_gaussian) {
    block1_prior <- .lmebayes_block1_prior_list(
      prior,
      group.dispersion = disp_info$dispersion_fix
    )

    out <- .lmebayes_run_lmm_engine(
      n             = n,
      design        = design,
      prior         = prior,
      disp_info     = disp_info,
      tv_tol        = tv_tol,
      progbar       = progbar,
      verbose       = verbose,
      gap_tol       = gap_tol,
      mode_gap_max  = mode_gap_max,
      sim_method    = sim_method,
      offset        = offset,
      weights       = weights,
      offset_missing = missing(offset),
      weights_missing = missing(weights)
    )

    icm_lbl <- .lmebayes_block2_icm_labels(prior, family)
    .lmebayes_print_icm_fixef_table(
      prior_list = prior$pop.prior_list,
      re_names   = re_names,
      fixef_icm  = out$popef.mode,
      icm_info   = out$convergence_info$icm_info,
      ref_label  = icm_lbl$ref_label,
      icm_label  = icm_lbl$icm_label,
      conv_label = icm_lbl$conv_label,
      header     = "--- rglmerb: population effects ---",
      verbose    = verbose
    )

    out <- .lmebayes_add_popef_summaries(out)
    out$call        <- cl
    out$convergence <- out$convergence_info
    out$Prior       <- list(
      block1_prior          = block1_prior,
      pfamily_list          = prior$pfamily_list,
      group.dispersion      = disp_info$dispersion_fix,
      dispersion_mode       = disp_info$mode,
      dispersion_pfamily    = disp_info$dispersion_pfamily,
      dispersion_prior_list = disp_info$dispersion_prior_list
    )
    out$design      <- design
    out$family      <- family

    if (!is.null(out$pilot$n) && out$pilot$n > 0L) {
      .lmebayes_print_fixef_init(
        out$popef.init,
        re_names,
        verbose,
        header = "--- rglmerb: main-stage popef.init (pilot colMeans) ---"
      )
    }

    class(out)      <- c("rglmerb", "list")
    return(out)
  }

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  block1_prior <- .lmebayes_block1_prior_list(prior, group.dispersion = NULL)

  out <- .lmebayes_run_glmm_engine(
    n              = n,
    design         = design,
    prior          = prior,
    family         = family,
    gap_tol        = gap_tol,
    tv_tol         = tv_tol,
    mode_gap_max   = mode_gap_max,
    verbose        = verbose,
    progbar        = progbar,
    collect_block1 = collect_block1,
    offset         = offset,
    weights        = weights,
    offset_missing = missing(offset),
    weights_missing = missing(weights)
  )
  out$call <- cl

  icm_lbl <- .lmebayes_block2_icm_labels(prior, family)
  .lmebayes_print_icm_fixef_table(
    prior_list = prior$pop.prior_list,
    re_names   = re_names,
    fixef_icm  = out$popef.mode,
    icm_info   = out$convergence_info$icm_info,
    ref_label  = icm_lbl$ref_label,
    icm_label  = icm_lbl$icm_label,
    conv_label = icm_lbl$conv_label,
    header     = "--- rglmerb: population effects ---",
    verbose    = verbose
  )

  .lmebayes_print_ranef_mode_reference(
    out$groupef.mode, re_names, group_levels, verbose
  )

  if (!is.null(out$pilot$n) && out$pilot$n > 0L) {
    .lmebayes_print_fixef_init(
      out$popef.init,
      re_names,
      verbose,
      header = "--- rglmerb: main-stage popef.init (pilot colMeans) ---"
    )
  }

  out <- .lmebayes_add_popef_summaries(out)
  out$call        <- cl
  out$convergence <- out$convergence_info
  out$Prior       <- list(
    block1_prior          = block1_prior,
    pfamily_list          = prior$pfamily_list,
    group.dispersion      = disp_info$dispersion_fix,
    dispersion_mode       = disp_info$mode,
    dispersion_pfamily    = disp_info$dispersion_pfamily,
    dispersion_prior_list = disp_info$dispersion_prior_list
  )
  out$design      <- design
  out$family      <- family

  class(out) <- c("rglmerb", "list")
  out
}
