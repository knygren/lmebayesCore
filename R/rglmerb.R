#' Bayesian generalized linear mixed-effects model sampler
#'
#' Matrix-level sampler for \code{model_setup} design objects and normalized
#' prior containers. Routes by response family:
#' \itemize{
#'   \item \code{family = gaussian()} delegates to
#'     \code{\link{rLMMNormal_reg_known_vcov}},
#'     \code{\link{rLMMNormal_reg_estimated_vcov}},
#'     \code{\link{rLMMindepNormalGamma_reg_known_vcov}}, or
#'     \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}} according to
#'     \code{dispersion_ranef} and Block~2 \code{pfamily_list}.
#'   \item Non-Gaussian families delegate via \code{\link{rGLMM_reg}} routes
#'     (through internal \code{.lmebayes_run_glmm_engine()} and
#'     \code{REG_ROUTE_TABLE}; pilot stage always unless \code{n_pilot = 0L};
#'     routes differ in eigenvalue-bound complexity). Inner sweeps use
#'     \code{\link{rGLMM_sweep}}.
#' }
#'
#' For formula-level fitting, \code{glmerb()} in the lmebayes package wraps this sampler.
#' Design matrices are built with \code{\link{model_setup}} and priors with
#' \code{\link{Prior_Setup_lmebayes}}.
#'
#' @param n Integer. Number of independent chains in the main stage.
#' @param design A \code{model_setup} object (from \code{\link{model_setup}}).
#' @param prior A normalized prior container with \code{Sigma_ranef},
#'   \code{prior_list}, and \code{pfamily_list}.
#' @param family A \code{\link[stats]{family}} object. Default \code{poisson()}.
#' @param dispersion_ranef Observation-level measurement dispersion \eqn{\sigma^2}:
#'   required positive scalar for \code{family = gaussian()}, or a
#'   \code{\link[glmbayesCore]{dGamma}()} pfamily with \code{Inv_Dispersion = TRUE}; must be
#'   \code{NULL} (default) for \code{poisson()} and \code{binomial()}.
#' @param gap_tol Legacy mode--mean gap for deriving the pilot chain count when
#'   \code{tv_tol} is \code{NULL}. Ignored for Gaussian without ING Block~2
#'   components.
#' @param tv_tol Total variation tolerance for convergence calibration.
#'   Inner Gibbs sweeps and pilot chain counts are derived internally.
#' @param mode_gap_max Pilot inner-sweep calibration (non-Gaussian and
#'   Gaussian+ING only).
#' @param collect_block1 Collect Block~1 \code{coefficients} from main chains
#'   (non-Gaussian only).
#' @param verbose Print stage headers and diagnostics.
#' @param progbar Progress bars when \code{verbose} is \code{FALSE}.
#' @param sim_method Sampling engine for \code{family = gaussian()}:
#'   \code{"DEFAULT"} or \code{"TWO_BLOCK_GIBBS"}; see \code{\link{rlmerb}}.
#'   Ignored (two-block Gibbs is the only engine) for non-Gaussian families.
#' @return Object of class \code{c("rglmerb", "list")} with Block~2 fields in
#'   the \code{fixef.*} namespace, plus \code{ranef.mode}, \code{sigma2}
#'   (Gaussian only: scalar or length-\code{n} vector as for \code{\link{rlmerb}}),
#'   \code{sigma2.mean}, \code{Prior}, \code{design}, and \code{family}. When a
#'   pilot stage runs (ING Block~2 and/or dGamma measurement dispersion),
#'   \code{n_pilot}, \code{pilot}, and \code{pilot_chisq} are included;
#'   otherwise \code{n_pilot} is \code{0L}.
#' @seealso \code{\link{rlmerb}}, \code{\link{rLMMNormal_reg_known_vcov}},
#'   \code{\link{rLMMNormal_reg_estimated_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg_known_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}}, \code{\link{rGLMM_reg}},
#'   \code{\link{Prior_Setup_lmebayes}}
#' @name rglmerb
#' @title The Bayesian Generalized Linear Mixed-Effects Model Distribution
NULL

#' @rdname rglmerb
#' @export
rglmerb <- function(
    n,
    design,
    prior,
    family              = poisson(),
    dispersion_ranef    = NULL,
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

  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }

  is_gaussian <- identical(family$family, "gaussian")

  disp_info <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = design,
    fn_name          = "rglmerb"
  )

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)

  if (is_gaussian) {
    block1_prior <- .lmebayes_block1_prior_list(
      prior,
      dispersion_ranef = disp_info$dispersion_fix
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
      sim_method    = sim_method
    )

    icm_lbl <- .lmebayes_block2_icm_labels(prior, family)
    .lmebayes_print_icm_fixef_table(
      prior_list = prior$prior_list,
      re_names   = re_names,
      fixef_icm  = out$fixef.mode,
      icm_info   = out$icm_info,
      ref_label  = icm_lbl$ref_label,
      icm_label  = icm_lbl$icm_label,
      conv_label = icm_lbl$conv_label,
      header     = "--- rglmerb: Block 2 fixed effects ---",
      verbose    = verbose
    )

    out <- .lmebayes_add_fixef_summaries(out)
    out$call        <- cl
    out$convergence <- out$convergence_info
    out$Prior       <- list(
      block1_prior          = block1_prior,
      pfamily_list          = prior$pfamily_list,
      dispersion_ranef      = disp_info$dispersion_fix,
      dispersion_mode       = disp_info$mode,
      dispersion_pfamily    = disp_info$dispersion_pfamily,
      dispersion_prior_list = disp_info$dispersion_prior_list
    )
    out$design      <- design
    out$family      <- family

    if (!is.null(out$n_pilot) && out$n_pilot > 0L) {
      .lmebayes_print_fixef_init(
        out$fixef.init,
        re_names,
        verbose,
        header = "--- rglmerb: main-stage fixef.init (pilot colMeans) ---"
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

  block1_prior <- .lmebayes_block1_prior_list(prior, dispersion_ranef = NULL)

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
    collect_block1 = collect_block1
  )
  out$call <- cl

  icm_lbl <- .lmebayes_block2_icm_labels(prior, family)
  .lmebayes_print_icm_fixef_table(
    prior_list = prior$prior_list,
    re_names   = re_names,
    fixef_icm  = out$fixef.mode,
    icm_info   = out$icm_info,
    ref_label  = icm_lbl$ref_label,
    icm_label  = icm_lbl$icm_label,
    conv_label = icm_lbl$conv_label,
    header     = "--- rglmerb: Block 2 fixed effects ---",
    verbose    = verbose
  )

  .lmebayes_print_ranef_mode_reference(
    out$ranef.mode, re_names, group_levels, verbose
  )

  if (!is.null(out$n_pilot) && out$n_pilot > 0L) {
    .lmebayes_print_fixef_init(
      out$fixef.init,
      re_names,
      verbose,
      header = "--- rglmerb: main-stage fixef.init (pilot colMeans) ---"
    )
  }

  out <- .lmebayes_add_fixef_summaries(out)
  out$call        <- cl
  out$convergence <- out$convergence_info
  out$Prior       <- list(
    block1_prior          = block1_prior,
    pfamily_list          = prior$pfamily_list,
    dispersion_ranef      = disp_info$dispersion_fix,
    dispersion_mode       = disp_info$mode,
    dispersion_pfamily    = disp_info$dispersion_pfamily,
    dispersion_prior_list = disp_info$dispersion_prior_list
  )
  out$design      <- design
  out$family      <- family

  class(out) <- c("rglmerb", "list")
  out
}
