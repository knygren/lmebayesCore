#' The Bayesian Linear Mixed-Effects Model Distribution
#'
#' \code{rlmerb} generates posterior draws for Bayesian linear mixed models,
#' parallel to \code{\link[glmbayesCore]{rlmb}} and \code{\link{rglmerb}}. It
#' takes a \code{model_setup} design plus the same Block~2 / Block~1 prior
#' arguments as \code{\link{rLMMNormal_reg}} /
#' \code{\link{rLMMindepNormalGamma_reg}} (\code{pfamily_list},
#' \code{dispprior_list}), computes the ICM posterior mean internally, and
#' delegates sampling to
#' \code{\link{rLMMNormal_reg_known_vcov}},
#' \code{\link{rLMMNormal_reg_estimated_vcov}},
#' \code{\link{rLMMindepNormalGamma_reg_known_vcov}}, or
#' \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}} according to
#' \code{dispprior_list} and population \code{pfamily_list}.
#'
#' When \code{simulate = FALSE}, no MCMC draws are produced: the call
#' delegates to the internal \code{\link{rlmerb_point}} and returns the
#' exact joint Gaussian posterior mean at fixed variance-component plug-ins
#' (draw slots are \code{NULL}).
#'
#' For formula interfaces, \code{lmerb()} in the lmebayes package wraps this
#' function.
#'
#' @param n Integer. Number of stored draws (each draw is one full pass through
#'   \code{m_convergence} inner Gibbs sweeps). Required when
#'   \code{simulate = TRUE}; ignored when \code{simulate = FALSE}.
#' @param design A \code{model_setup} object (from \code{\link{model_setup}})
#'   supplying \code{y}, \code{D}, \code{group}, \code{W},
#'   \code{group_name}, and \code{groupef.names}.
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects, one per
#'   random-effect coefficient in \code{design$groupef.names} (same as
#'   \code{\link{rLMMNormal_reg}}).
#' @param dispprior_list Block~1 observation-dispersion prior, in the same
#'   shapes accepted by the \code{rLMM*} engines: \code{list(dispersion = )}
#'   for fixed \eqn{\sigma^2} (scalar or per-group vector); a
#'   \code{\link[glmbayesCore]{dGamma}()} pfamily (pooled ING); or a named
#'   list of \code{dGamma()} pfamilies (per-group ING). A bare positive
#'   numeric scalar/vector is also accepted (treated as fixed dispersion).
#' @param offset,weights Observation offset and prior weights, as in
#'   \code{\link{rLMMNormal_reg}} (\code{offset = NULL}, \code{weights = 1}).
#'   When omitted, inherit \code{design$offset} / \code{design$weights} from
#'   \code{\link{model_setup}} if present. Explicit values always override
#'   the design. Echoed on the fit as \code{prior.weights} / \code{offset} /
#'   \code{offset2}; not yet used in Gibbs sweeps.
#' @param tv_tol Single numeric in \code{(0, 1)}. Total variation tolerance
#'   used for convergence calibration. Default \code{0.01}.
#'   Inner Gibbs sweeps per stored draw are derived from Theorem~3.
#'   Ignored when \code{simulate = FALSE}.
#' @param gap_tol Legacy mode--mean gap tolerance for the pilot stage when
#'   any population component uses \code{dIndependent_Normal_Gamma} and
#'   \code{tv_tol} is \code{NULL}. Ignored for all-\code{dNormal} models.
#'   Ignored when \code{simulate = FALSE}.
#' @param mode_gap_max Pilot inner-sweep calibration for ING population models
#'   (default \code{1.0}). Ignored for all-\code{dNormal} models.
#'   Ignored when \code{simulate = FALSE}.
#' @param progbar Logical. Show a text progress bar during sampling.
#'   Default \code{TRUE}. Ignored when \code{simulate = FALSE}.
#' @param verbose Logical. Print the reference-vs-ICM table and the convergence
#'   calibration line. Default \code{TRUE}.
#' @param print_icm_table Logical. When \code{FALSE}, skip the reference-vs-ICM
#'   table. The convergence calibration line from the Core sampler still follows
#'   \code{verbose}. Default \code{TRUE}.
#' @param diag_sweeps Diagnostic flag for ING population models with a pilot stage.
#'   When \code{TRUE}, auto-print one combined population chain-mean table per
#'   stage when each stage finishes; \code{sweep_history} is always stored on
#'   the fit. Default \code{FALSE}. Ignored when \code{simulate = FALSE}.
#' @param sim_method Simulation method: \code{"DEFAULT"} or
#'   \code{"TWO_BLOCK_GIBBS"}. Only affects the fixed-dispersion,
#'   known-variance-components route (scalar or per-group fixed
#'   \code{dispprior_list}, all population \code{dNormal()}): \code{"DEFAULT"}
#'   draws directly from the exact multivariate-normal posterior via
#'   \code{\link{rLMMNormal_reg_known_vcov_iid}} (no Gibbs sweeps, no
#'   burn-in); \code{"TWO_BLOCK_GIBBS"} forces two-block Gibbs sampling
#'   (\code{\link{rLMMNormal_reg_known_vcov_two_bg}}) instead. Every other
#'   route (\code{dGamma()}/\code{dIndependent_Normal_Gamma} components,
#'   variance components estimated) only supports two-block Gibbs, so both
#'   values behave identically there. Ignored when \code{simulate = FALSE}.
#' @param simulate Logical (default \code{TRUE}). When \code{TRUE}, draw
#'   posterior samples. When \code{FALSE}, return point estimates only via
#'   \code{\link{rlmerb_point}} (no MCMC).
#' @return An object of class \code{c("rlmerb", "list")} with population fields
#'   in the \code{popef.*} namespace, group draws in \code{groupef}/
#'   \code{groupef.mode}, observation residual variance in
#'   \code{group.dispersion}/\code{group.dispersion.mean} (fixed scalar or
#'   vector, or draws when a Gamma measurement prior was used),
#'   \code{m_convergence}, \code{convergence},
#'   \code{convergence_info$sim_method_used} (\code{"DEFAULT"} or
#'   \code{"TWO_BLOCK_GIBBS"}), full packed \code{prior}
#'   (from \code{\link{priors_from_pfamily_list}}), thin \code{Prior}, and
#'   \code{design}. When \code{simulate = FALSE}, draw / pilot / convergence
#'   slots are \code{NULL}.
#' @seealso \code{\link{rglmerb}}, \code{\link{rlmerb_point}},
#'   \code{\link{rLMMNormal_reg_known_vcov}},
#'   \code{\link{rLMMNormal_reg_estimated_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg_known_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}},
#'   \code{\link{Prior_Setup_GLMM}},
#'   \code{\link[glmbayesCore]{rlmb}}
#' @export
rlmerb <- function(
    n,
    design,
    pfamily_list,
    dispprior_list,
    offset              = NULL,
    weights             = 1,
    tv_tol        = 0.01,
    progbar         = TRUE,
    verbose         = TRUE,
    print_icm_table = TRUE,
    gap_tol             = 0.0196,
    mode_gap_max        = 1.0,
    diag_sweeps         = FALSE,
    sim_method          = "DEFAULT",
    simulate            = TRUE
) {
  cl <- match.call()

  if (!isTRUE(simulate)) {
    out <- rlmerb_point(
      design            = design,
      pfamily_list      = pfamily_list,
      dispprior_list    = dispprior_list,
      offset            = offset,
      weights           = weights,
      verbose           = verbose,
      print_icm_table   = print_icm_table,
      offset_missing    = missing(offset),
      weights_missing   = missing(weights)
    )
    out$call <- cl
    return(out)
  }

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  sim_method <- .rLMM_validate_sim_method(sim_method, fn_name = "rlmerb")

  if (missing(pfamily_list) || is.null(pfamily_list)) {
    stop(
      "'pfamily_list' is required for rlmerb() (named Block~2 pfamily list, ",
      "same as rLMMNormal_reg()).",
      call. = FALSE
    )
  }
  if (missing(dispprior_list)) {
    stop(
      "'dispprior_list' is required for rlmerb(): list(dispersion = ...) for ",
      "fixed sigma^2, or a dGamma()/dGamma_list() measurement prior ",
      "(same as rLMM* engines).",
      call. = FALSE
    )
  }

  prep <- .rlmerb_prepare_prior(
    design          = design,
    pfamily_list    = pfamily_list,
    dispprior_list  = dispprior_list,
    family          = stats::gaussian(),
    fn_name         = "rlmerb",
    offset          = offset,
    weights         = weights,
    offset_missing  = missing(offset),
    weights_missing = missing(weights)
  )
  prior <- prep$prior
  disp_info <- prep$disp_info
  re_names <- prep$re_names
  block1_prior <- prep$block1_prior

  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  out <- .lmebayes_run_lmm_engine(
    n               = n,
    design          = design,
    prior           = prior,
    disp_info       = disp_info,
    tv_tol          = tv_tol,
    progbar         = progbar,
    verbose         = verbose,
    gap_tol             = gap_tol,
    mode_gap_max        = mode_gap_max,
    diag_sweeps         = diag_sweeps,
    sim_method          = sim_method,
    offset              = offset,
    weights             = weights,
    offset_missing      = missing(offset),
    weights_missing     = missing(weights)
  )

  if (isTRUE(print_icm_table)) {
    icm_lbl <- .lmebayes_block2_icm_labels(prior, stats::gaussian())
    .lmebayes_print_icm_fixef_table(
      prior_list = prior$pop.prior_list,
      re_names   = re_names,
      fixef_icm  = out$popef.mode,
      icm_info   = out$convergence_info$icm_info,
      ref_label  = icm_lbl$ref_label,
      icm_label  = icm_lbl$icm_label,
      conv_label = icm_lbl$conv_label,
      header     = "--- rlmerb: population effects ---",
      verbose    = verbose
    )
  }

  out <- .lmebayes_add_popef_summaries(out)
  out$call       <- cl
  out$convergence <- out$convergence_info
  out$prior      <- prior
  out$Prior      <- .rlmerb_thin_Prior(prior, disp_info, block1_prior)
  out$design     <- design

  if (!is.null(out$pilot$n) && out$pilot$n > 0L) {
    .lmebayes_print_fixef_init(
      out$popef.init,
      re_names,
      verbose,
      header = "--- rlmerb: main-stage popef.init (pilot colMeans) ---"
    )
  }

  class(out) <- c("rlmerb", "list")
  out
}
