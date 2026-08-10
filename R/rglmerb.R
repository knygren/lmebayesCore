#' The Bayesian Generalized Linear Mixed-Effects Model Distribution
#'
#' \code{rglmerb} draws from the posterior of a Bayesian generalized linear
#' mixed model --- Poisson, binomial, Gamma, or Gaussian --- with a single
#' grouping factor and uncorrelated random effects, given a design from
#' \code{\link{model_setup}} and priors for the population level
#' (\code{pfamily_list}) and, for Gaussian responses only, the residual
#' variance (\code{dispprior_list}).
#'
#' The returned sample is \code{n} mutually independent draws, one per
#' replicate chain, taken at that chain's final sweep. It is analyzed like an
#' iid sample rather than like an MCMC trace.
#'
#' This is the matrix-level interface; for a formula interface use
#' \code{glmerb()} in \pkg{lmebayes}. \code{\link{rlmerb}} is the
#' Gaussian-only counterpart.
#'
#' @details
#' ## The model being estimated
#'
#' \code{rglmerb} estimates the same two-stage structure as
#' \code{\link{rlmerb}}, with a link function between the linear predictor
#' and the response mean. For group \eqn{j} with link \eqn{g}:
#' \deqn{y_j \mid \beta_j \sim \mathrm{family}\big(g^{-1}(D_j \beta_j)\big)
#'       \qquad \textrm{(within group)}}
#' \deqn{\beta_j \mid \gamma, \Psi \sim N(\mathcal{W}_j \gamma,\; \Psi)
#'       \qquad \textrm{(across groups)}}
#'
#' The across-group stage is unchanged and remains Gaussian: group
#' coefficients scatter around a population expectation
#' \eqn{\mathcal{W}_j\gamma} with diagonal covariance
#' \eqn{\Psi = \mathrm{diag}(\tau^2_k)}. What changes is the within-group
#' stage and, with it, the interpretation of \eqn{\beta_j}, which now lives
#' on the link scale --- log rates for \code{poisson()}, log-odds for
#' \code{binomial()}. Summaries that are meaningful on the response scale
#' need the inverse link applied to the draws.
#'
#' Supported families are \code{gaussian()}, \code{poisson()},
#' \code{binomial()}, and \code{Gamma()} with their standard links. The
#' restriction is log-concavity of the likelihood, which the sampler
#' requires \insertCite{Nygren2006}{lmebayesCore}.
#'
#' ## What you supply
#'
#' \describe{
#'   \item{\code{design}}{A \code{\link{model_setup}} object built with the
#'     matching \code{family}. For non-Gaussian responses its estimability
#'     report matters more than usual: a group with no variation in the
#'     response (all-zero counts, a single binomial outcome level) has no
#'     finite maximum likelihood estimate, and is reported separately from
#'     the rank check. See
#'     \code{vignette("Chapter-B08", package = "lmebayesCore")}.}
#'   \item{\code{pfamily_list}}{The population-level prior, one component per
#'     random-effect coefficient, declaring whether each between-group
#'     variance \eqn{\tau^2_k} is known (\code{dNormal}) or estimated
#'     (\code{dIndependent_Normal_Gamma}).}
#'   \item{\code{dispprior_list}}{Only meaningful for \code{gaussian()},
#'     where it gives the residual variance. \code{poisson()} and
#'     \code{binomial()} have no free dispersion parameter, so it must be
#'     \code{NULL}.}
#' }
#'
#' ## What the returned draws are
#'
#' As in \code{\link{rlmerb}}, \code{rglmerb} runs \code{n} independent
#' replicate chains and stores one draw from each, taken at the final sweep,
#' so the \code{n} draws are mutually independent and ordinary Monte Carlo
#' standard errors apply. There is nothing to thin and no autocorrelation to
#' diagnose.
#'
#' The strength of the accompanying guarantee is weaker here than in the
#' Gaussian case, and the difference is worth understanding rather than
#' glossing over. The calibration behind \code{m_convergence} is derived for
#' Gaussian targets. A Poisson or binomial posterior is not Gaussian, so the
#' engine applies the same calibration at curvature evaluated near the
#' posterior mode. The sweep count is therefore a well-motivated
#' \strong{budget} rather than a certificate, and \code{tv_tol} should be
#' read in that spirit. Only \code{family = gaussian()} with fixed
#' \eqn{\sigma^2} and all-\code{dNormal} components yields exact draws.
#'
#' A second difference follows from the same fact. Because the posterior mode
#' is not the posterior mean when the target is skewed, chains started at the
#' mode begin at an unknown distance from the target. \code{rglmerb}
#' therefore runs a short \strong{pilot stage} first and uses its cross-chain
#' mean as the starting point for the main chains, which shortens the main
#' stage. Pilot draws are diagnostic only and are discarded --- they are
#' never part of the returned sample.
#' \code{vignette("Chapter-B03", package = "lmebayesCore")} covers the
#' trade-off between pilot size and main-stage length.
#'
#' ## Point estimates without sampling
#'
#' \code{simulate = FALSE} returns the joint posterior mode by iterated
#' conditional modes (or the exact mean for \code{gaussian()}), at fixed
#' variance-component plug-ins. No sweeps run and no RNG is consumed, so the
#' result is deterministic --- useful while iterating on model
#' specification. All draw, pilot, and convergence slots are \code{NULL}.
#'
#' @references
#'   \insertAllCited{}
#' @param n Integer. Number of independent chains in the main stage.
#'   Required when \code{simulate = TRUE}; ignored when \code{simulate = FALSE}.
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
#'   components. Ignored when \code{simulate = FALSE}.
#' @param tv_tol Total variation tolerance for convergence calibration.
#'   Inner Gibbs sweeps and pilot chain counts are derived internally.
#'   Ignored when \code{simulate = FALSE}.
#' @param mode_gap_max Pilot inner-sweep calibration (non-Gaussian and
#'   Gaussian+ING only). Ignored when \code{simulate = FALSE}.
#' @param collect_block1 Collect \code{groupef} draws from main chains
#'   (non-Gaussian only). Ignored when \code{simulate = FALSE}.
#' @param verbose Print stage headers and diagnostics.
#' @param progbar Progress bars when \code{verbose} is \code{FALSE}.
#'   Ignored when \code{simulate = FALSE}.
#' @param sim_method Simulation method for \code{family = gaussian()}:
#'   \code{"DEFAULT"} or \code{"TWO_BLOCK_GIBBS"}; see \code{\link{rlmerb}}.
#'   Ignored (two-block Gibbs only) for non-Gaussian families.
#'   Ignored when \code{simulate = FALSE}.
#' @param simulate Logical (default \code{TRUE}). When \code{TRUE}, draw
#'   posterior samples. When \code{FALSE}, return point estimates only via
#'   \code{\link{rglmerb_point}} (no MCMC).
#' @return An object of class \code{c("rglmerb", "list")}, with the same
#'   layout as \code{\link{rlmerb}}. Each draw slot holds one row per
#'   replicate chain, taken at that chain's final sweep.
#'
#'   \describe{
#'     \item{\code{popef}}{Named list of population draws \eqn{\gamma}, one
#'       \code{n} \eqn{\times q_k} matrix per random-effect coefficient, on
#'       the link scale. \code{popef.means} and \code{popef.mode} are the
#'       draw mean and the posterior mode.}
#'     \item{\code{groupef}}{Long-format data frame of group coefficients
#'       \eqn{\beta_j} (\code{n} \eqn{\times J} rows: \code{draw}, the
#'       grouping factor, then one column per coefficient), also on the link
#'       scale. \code{groupef.mode} is the \eqn{J \times p_{re}} mode matrix.
#'       Collected only when \code{collect_block1 = TRUE}.}
#'     \item{\code{popef.dispersion}}{Draws of the between-group variances
#'       \eqn{\tau^2_k}; constant for \code{dNormal} components, varying for
#'       estimated ones. \code{popef.dispersion.mean} is the column mean.}
#'     \item{\code{group.dispersion}}{Residual variance, Gaussian responses
#'       only; absent for \code{poisson()} and \code{binomial()}.}
#'     \item{\code{m_convergence}, \code{convergence_info}}{Sweeps run per
#'       chain and the calibration behind that count (\code{method},
#'       \code{lambda_star}, \code{draw_engine},
#'       \code{sim_method_used}, \code{icm_info}).}
#'     \item{\code{pilot}}{Present when a pilot stage ran, holding its size
#'       \code{n}, the \code{chisq} calibration constant, and its
#'       \code{draws}. Diagnostic only --- these are not posterior draws.}
#'   }
#'
#'   Inputs are echoed back in \code{design}, \code{family}, \code{n},
#'   \code{call}, \code{pfamily_list}, \code{dispprior_list},
#'   \code{prior.weights} and \code{offset}, alongside the packed
#'   \code{prior} and thin \code{Prior} summaries.
#'
#'   When \code{simulate = FALSE}, every draw, pilot and convergence slot is
#'   \code{NULL} and only the mode fields are populated.
#' @seealso \code{\link{rlmerb}}, \code{\link{rglmerb_point}},
#'   \code{\link{rLMMNormal_reg_known_vcov}},
#'   \code{\link{rLMMNormal_reg_estimated_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg_known_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}}, \code{\link{rGLMM_reg}},
#'   \code{\link{Prior_Setup_GLMM}}
#' @example inst/examples/Ex_rglmerb.R
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
    sim_method          = "DEFAULT",
    simulate            = TRUE
) {
  cl <- match.call()

  if (!isTRUE(simulate)) {
    out <- rglmerb_point(
      design            = design,
      pfamily_list      = pfamily_list,
      family            = family,
      dispprior_list    = dispprior_list,
      offset            = offset,
      weights           = weights,
      verbose           = verbose,
      print_icm_table   = TRUE,
      offset_missing    = missing(offset),
      weights_missing   = missing(weights)
    )
    out$call <- cl
    return(out)
  }

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  sim_method <- .rLMM_validate_sim_method(sim_method, fn_name = "rglmerb")

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

  prep <- .rlmerb_prepare_prior(
    design          = design,
    pfamily_list    = pfamily_list,
    dispprior_list  = dispprior_list,
    family          = family,
    fn_name         = "rglmerb",
    offset          = offset,
    weights         = weights,
    offset_missing  = missing(offset),
    weights_missing = missing(weights)
  )
  prior <- prep$prior
  disp_info <- prep$disp_info
  re_names <- prep$re_names
  block1_prior <- prep$block1_prior
  group_levels <- levels(design$group)

  if (is_gaussian) {
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
    out$prior       <- prior
    out$Prior       <- .rlmerb_thin_Prior(prior, disp_info, block1_prior)
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

  ## Non-Gaussian: rebuild block1 with NULL dispersion (prepare used gaussian rule).
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
  out$prior       <- prior
  out$Prior       <- .rlmerb_thin_Prior(prior, disp_info, block1_prior)
  out$design      <- design
  out$family      <- family

  class(out) <- c("rglmerb", "list")
  out
}
