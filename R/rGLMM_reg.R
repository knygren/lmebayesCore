#' Simulation Functions for Generalized Linear Mixed Models
#'
#' @description
#' Simulation functions for generating posterior draws from Bayesian
#' generalized linear mixed models (and optionally Gaussian models). Parallel
#' to the glmbayes \code{\link[glmbayes]{simfuncs}} for GLMs; typically called
#' from \code{\link{rglmerb}} (and from \code{glmerb()} in \strong{lmebayes}).
#' Each stored draw runs \code{m_convergence} inner Gibbs sweeps via
#' \code{\link{rGLMM_sweep}}. Gaussian models with observation dispersion
#' usually use the \code{\link{rLMM_reg}} routes via \code{\link{rlmerb}}.
#'
#' @section Simulation routes:
#' Both routes run a \strong{pilot stage for non-Gaussian} families (local-Gaussian
#' rate calibration and chain-mean initialisation; skip only with
#' \code{n_pilot = 0L}). The route split is \strong{not} whether a pilot runs,
#' but how \strong{eigenvalue bounds} are set for Theorem~3 and post-pilot
#' upper-bound calibration:
#' \describe{
#'   \item{\code{rGLMM_reg_known_vcov}}{
#'     All population components \code{dNormal} (known \eqn{\tau^2_k}): standard
#'     fixed-dispersion rate at the mode; post-pilot eigenvalue upper bound from
#'     pilot \code{groupef} without ING \code{disp_lower} conservatism.}
#'   \item{\code{rGLMM_reg_estimated_vcov}}{
#'     At least one ING population component: conservative \code{disp_lower}
#'     plug-in in \code{\link{two_block_rate_from_pfamily_list}} (upper bound
#'     over truncated \eqn{\tau^2}); pilot updates \eqn{\tau^2} starts from
#'     dispersion draws.}
#' }
#'
#' @section Dispatcher:
#' \code{\link{rGLMM_reg}} validates population \code{pfamily_list} and delegates
#' to the appropriate route.
#'
#' @param n Number of stored main-stage draws. If \code{length(n) > 1}, the
#'   length is used.
#' @param y Response vector of length \code{l2} (\code{= nrow(D)}).
#' @param D Level-1 design matrix (\code{l2 x p_re}). Must have unique,
#'   non-empty \code{colnames(D)}: these are the random-effect coefficient
#'   names used to key \code{W} and \code{pfamily_list} (there is no
#'   separate \code{groupef.names} argument to override them).
#' @param group Grouping factor of length \code{l2} (must be a \code{factor};
#'   \code{levels(group)} fixes the row order of \code{groupef} draws -- there
#'   is no separate \code{group_levels} argument. To use a level order/superset
#'   not present in the observed data, construct \code{group} as
#'   \code{factor(observed_group, levels = full_superset)} yourself. The name
#'   used for the grouping column in \code{groupef} (\code{group_name})
#'   is resolved from \code{attr(group, "group_name")} if set, otherwise from
#'   \code{group}'s own variable name via \code{substitute()} -- this only
#'   works when \code{group} is passed as a bare variable (e.g.
#'   \code{group = school_id}); otherwise attach the name yourself via
#'   \code{attr(group, "group_name") <- "school_id"}.
#' @param W Named list of group-level design matrices (\code{J x q_k}),
#'   one per column of \code{D}.
#' @param pfamily_list Named list of population \code{pfamily} objects. One
#'   \eqn{\tau^2_k} plug-in per component (fixed \code{dispersion} for
#'   \code{dNormal}, prior mean \eqn{rate/(shape - 1)} for
#'   \code{dIndependent_Normal_Gamma}) is assembled into the diagonal group
#'   random-effect prior precision.
#' @param dispprior_list Observation-dispersion prior container:
#'   \code{dispersion} (required for \code{gaussian()}), optional \code{ddef}.
#'   For non-Gaussian families this is often an empty list. The group
#'   random-effect prior precision (formerly a separate \code{P}/\code{Sigma}
#'   field) is always derived internally from \code{pfamily_list};
#'   \code{dispprior_list} must not contain \code{P} or \code{Sigma}.
#' @param icm_tol,icm_maxit ICM convergence controls for the internal population start.
#' @param offset,weights Observation offset and prior weights (glmbayes-style
#'   formals: \code{offset = NULL}, \code{weights = 1}). Normalized to length
#'   \code{length(y)} and echoed on the return as \code{offset}/
#'   \code{offset2}/\code{prior.weights}. \strong{Not yet used} by the
#'   mixed-model sampling path (ICM / sweeps still assume unit weights and
#'   zero offset).
#' @param family Likelihood family (length-\code{l2} recycling rules follow
#'   the group-level likelihood).
#' @param gap_tol Legacy mode--mean gap for pilot chain count when \code{tv_tol}
#'   is \code{NULL}.
#' @param tv_tol Total-variation tolerance for Theorem~3 calibration.
#'   \code{n_pilot}, \code{m_convergence_pilot}, and main-stage
#'   \code{m_convergence} are derived internally.
#' @param mode_gap_max Pilot inner-sweep calibration when a pilot stage runs.
#' @param Gridtype,n_envopt,use_parallel,use_opencl Reserved (not yet forwarded).
#' @param verbose Print stage headers and convergence calibration lines.
#' @param progbar Progress bars during sampling.
#' @param stage_verbose Print pilot chi-squared and post-pilot UB diagnostics.
#' @param rate_calibration Optional rate object for \code{stage_verbose}.
#' @param collect_block1 Collect \code{groupef} draws from main chains.
#'
#' @return An object of class \code{c("<route>", "rGLMM_reg", "list")}, where
#'   \code{<route>} is \code{"rGLMM_reg_known_vcov"} or
#'   \code{"rGLMM_reg_estimated_vcov"}. Components use package
#'   \strong{group}/\strong{population} names (see \file{inst/notation.md}),
#'   in glm/glmbayes-style order:
#'   \describe{
#'     \item{\code{groupef}}{Draws of non-centered group coefficients
#'       \eqn{\beta_j} (when \code{collect_block1} is \code{TRUE}): grouping
#'       column plus one column per \code{colnames(D)}.}
#'     \item{\code{groupef.mode}}{\eqn{J \times p_{re}} matrix of ICM group
#'       coefficients \eqn{\hat\beta_j}. \strong{Not} \code{lme4}'s mean-zero
#'       \eqn{u_j}.}
#'     \item{\code{groupef.iters}}{Optional group-level envelope iteration counts.}
#'     \item{\code{popef}}{Named list of \code{n x q_k} matrices of population
#'       coefficient draws \eqn{\gamma_k}.}
#'     \item{\code{popef.mode}, \code{popef.init}}{ICM population point estimates
#'       and main-stage starts.}
#'     \item{\code{popef.dispersion}, \code{popef.iters}}{Optional population
#'       per-draw diagnostics when produced by the sampler.}
#'     \item{\code{pfamily_list}, \code{dispprior_list}}{Population and
#'       observation-dispersion priors that were used.}
#'     \item{\code{prior.weights}, \code{offset}, \code{offset2}}{Normalized
#'       copies of the \code{weights}/\code{offset} arguments (glmbayes
#'       naming). Not yet consumed by sampling.}
#'     \item{\code{any_non_normal}}{Whether any population component is not
#'       \code{dNormal}.}
#'     \item{\code{family}, \code{design}, \code{n}}{Likelihood family, matrix
#'       inputs (including echoed \code{weights}/\code{offset}), and chain
#'       count.}
#'     \item{\code{call}}{Matched call.}
#'     \item{\code{m_convergence}}{Inner Gibbs sweeps per stored main-stage
#'       draw.}
#'     \item{\code{convergence_info}}{Theorem~3 / UB calibration details,
#'       including \code{draw_engine} and \code{icm_info}.}
#'     \item{\code{pilot}}{When a pilot ran: list with \code{n},
#'       \code{m_convergence}, \code{chisq}, and \code{draws}; otherwise
#'       \code{NULL}.}
#'     \item{\code{sweep_history}}{Main-stage sweep history when collected.}
#'   }
#'
#' @param ... further arguments passed to or from other methods.
#' @family simfuncs
#' @seealso \code{\link{rGLMM_sweep}}, \code{\link{rLMM_reg}},
#'   \code{\link{rglmerb}}, \code{\link{print_groupef}}
#' @example inst/examples/Ex_rGLMM_reg.R
#' @name rGLMM_reg
#' @order 1
NULL

#' Shared matrix-level validation for GLMM replicate-chain engines
#'
#' \code{groupef.names} and \code{group_levels} are no longer separate
#' arguments: they are always \code{colnames(D)} and \code{levels(group)}
#' respectively. \code{group_name} must already be resolved by the caller
#' (see \code{\link{.lmebayes_resolve_group_name}}); this function only
#' sanity-checks it. \code{dispprior_list} must not contain \code{P}/
#' \code{Sigma}: the group random-effect prior precision is derived
#' internally from \code{pfamily_list} and injected here.
#' @noRd
.rGLMM_validate_matrix_inputs <- function(
    n,
    y,
    D,
    group,
    W,
    tv_tol,
    group_name,
    family,
    mode_gap_max,
    gap_tol,
    pfamily_list,
    dispprior_list
) {
  family <- .two_block_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  gap_tol <- .two_block_validate_gap_tol(gap_tol)

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be a single positive finite number.",
           call. = FALSE)
    }
  }

  y <- as.vector(y)
  D <- as.matrix(D)
  l2 <- nrow(D)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(D).", call. = FALSE)
  }

  re_names <- colnames(D)
  if (is.null(re_names) || length(re_names) != ncol(D) || anyNA(re_names) ||
      any(!nzchar(re_names)) || anyDuplicated(re_names)) {
    stop(
      "'D' must have unique, non-empty column names (colnames(D)); ",
      "there is no 'groupef.names' argument to override this.",
      call. = FALSE
    )
  }

  if (!is.factor(group)) {
    stop(
      "'group' must be a factor (wrap with factor(group, levels = ...) ",
      "to control level order or supply a fixed superset of levels); ",
      "there is no 'group_levels' argument to override this.",
      call. = FALSE
    )
  }
  group_levels <- levels(group)
  if (length(group_levels) < 1L) {
    stop("'group' must have at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    stop(
      "'group_name' could not be derived and was not supplied; ",
      "pass 'group_name' explicitly.",
      call. = FALSE
    )
  }

  if (!is.list(W) || is.data.frame(W)) {
    stop("'W' must be a list of design matrices.", call. = FALSE)
  }
  if (length(W) != length(re_names)) {
    stop("length(W) must equal ncol(D) = ", length(re_names), ".",
         call. = FALSE)
  }
  if (!setequal(names(W), re_names)) {
    stop(
      "names(W) must match colnames(D): ",
      paste(re_names, collapse = ", "), ".", call. = FALSE
    )
  }
  W <- W[re_names]

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, re_names, J = length(group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)

  if (!is.null(tv_tol)) {
    if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
        !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
      stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
    }
  }

  if (!is.null(dispprior_list$P) || !is.null(dispprior_list$Sigma)) {
    stop(
      "'dispprior_list' must not contain 'P'/'Sigma'; the group random-effect ",
      "prior precision is derived internally from 'pfamily_list'.",
      call. = FALSE
    )
  }
  dispprior_list$P <- .rLMM_P_from_pfamily_list(pfamily_list, re_names)

  .two_block_validate_block1_prior(dispprior_list, family = family)

  list(
    n              = n,
    y              = y,
    D              = D,
    group          = group,
    W              = W,
    re_names       = re_names,
    group_levels   = group_levels,
    group_name     = group_name,
    family         = family,
    is_gaussian    = is_gaussian,
    dispprior_list = dispprior_list,
    pfamily_list   = pfamily_list,
    pf_summary     = pf_summary,
    ptypes         = pf_summary$ptypes,
    any_non_normal = pf_summary$any_non_normal,
    tv_tol         = tv_tol,
    gap_tol        = gap_tol
  )
}

#' Main GLMM sampling pipeline entry (ICM, optional pilot, main sweep)
#' @noRd
.rGLMM_reg_run <- function(
    inp,
    icm_tol,
    icm_maxit,
    mode_gap_max,
    verbose,
    progbar,
    stage_verbose,
    rate_calibration,
    collect_block1,
    engine_label,
    result_class,
    cl
) {
  .rGLMM_reg_run_pipeline(
    inp                = inp,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    mode_gap_max       = mode_gap_max,
    verbose            = verbose,
    progbar            = progbar,
    stage_verbose      = stage_verbose,
    rate_calibration   = rate_calibration,
    collect_block1     = collect_block1,
    engine_label       = engine_label,
    result_class       = result_class,
    cl                 = cl
  )
}

#' Full GLMM sampling pipeline (same as \code{.rGLMM_reg_run}; reserved for
#' future pilot-policy split)
#' @noRd
.rGLMM_reg_run_with_pilot <- function(
    inp,
    icm_tol,
    icm_maxit,
    mode_gap_max,
    verbose,
    progbar,
    stage_verbose,
    rate_calibration,
    collect_block1,
    engine_label,
    result_class,
    cl
) {
  .rGLMM_reg_run_pipeline(
    inp                = inp,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    mode_gap_max       = mode_gap_max,
    verbose            = verbose,
    progbar            = progbar,
    stage_verbose      = stage_verbose,
    rate_calibration   = rate_calibration,
    collect_block1     = collect_block1,
    engine_label       = engine_label,
    result_class       = result_class,
    cl                 = cl
  )
}

#' Shared GLMM replicate-chain pipeline (sweep-outer driver)
#' @noRd
.rGLMM_reg_run_pipeline <- function(
    inp,
    icm_tol,
    icm_maxit,
    mode_gap_max,
    verbose,
    progbar,
    stage_verbose,
    rate_calibration,
    collect_block1,
    engine_label,
    result_class,
    cl
) {
  n              <- inp$n
  y              <- inp$y
  D              <- inp$D
  group          <- inp$group
  W              <- inp$W
  re_names       <- inp$re_names
  group_levels   <- inp$group_levels
  group_name     <- inp$group_name
  family         <- inp$family
  is_gaussian    <- inp$is_gaussian
  dispprior_list     <- inp$dispprior_list
  pfamily_list   <- inp$pfamily_list
  pf_summary     <- inp$pf_summary
  ptypes         <- inp$ptypes
  any_non_normal <- inp$any_non_normal
  tv_tol         <- inp$tv_tol
  gap_tol        <- inp$gap_tol
  n_pilot_arg         <- NULL
  m_convergence_pilot <- NULL

  will_pilot <- .two_block_pilot_will_run(
    is_gaussian, n_pilot_arg, gap_tol, tv_tol,
    any_non_normal = any_non_normal
  )
  run_pilot <- will_pilot
  run_ub    <- will_pilot && !is.null(tv_tol)

  if (run_pilot && is.null(m_convergence_pilot)) {
    m_convergence_pilot <- if (!is.null(tv_tol)) {
      NULL
    } else {
      10L
    }
  }

  icm_info <- NULL
  design_icm <- list(
    y             = y,
    D             = D,
    group         = factor(group, levels = group_levels),
    W             = W,
    groupef.names = re_names,
    group_name    = group_name
  )
  icm <- .two_block_icm_at_start(
    design       = design_icm,
    prior_list   = dispprior_list,
    pfamily_list = pfamily_list,
    re_names     = re_names,
    family       = family,
    tol          = icm_tol,
    maxit        = icm_maxit
  )
  fixef_mode <- icm$start
  b_start    <- icm$b_start
  ranef_mode <- b_start
  icm_info   <- icm$icm
  if (isTRUE(verbose)) {
    if (isTRUE(any_non_normal)) {
      icm_what <- "population start at lmer tau^2 plug-in"
    } else if (is_gaussian) {
      icm_what <- "ICM posterior mean"
    } else {
      icm_what <- "ICM posterior mode"
    }
    cat(sprintf(
      "  %s: %s (converged: %s, %d iter, delta = %.2e)\n\n",
      engine_label,
      icm_what,
      icm_info$converged,
      icm_info$iterations,
      icm_info$delta
    ))
  }

  design <- list(
    y             = y,
    D             = D,
    group         = factor(group, levels = group_levels),
    W             = W,
    groupef.names = re_names,
    group_name    = group_name
  )

  fixef_mode_ref <- fixef_mode
  b_mode_ref     <- b_start
  progbar_use    <- isTRUE(progbar) || isTRUE(verbose) || isTRUE(stage_verbose)

  rate <- .rGLMM_rate_at_mode(
    design       = design,
    dispprior_list   = dispprior_list,
    pfamily_list = pfamily_list,
    family       = family,
    b_mode       = b_start,
    group_levels = group_levels,
    is_gaussian  = is_gaussian
  )

  m_min <- NULL
  if (!is.null(tv_tol)) {
    m_min <- .two_block_cap_inner_sweeps(
      two_block_l_for_tv(rate, tv_tol, method = "theorem3") + 1L
    )
  }

  p_dim            <- sum(vapply(fixef_mode, length, integer(1L)))
  D_max            <- if (!is.null(mode_gap_max)) sqrt(p_dim) * mode_gap_max else 0
  m_pilot_from_gap <- NULL

  if (run_pilot && is.null(m_convergence_pilot) && !is.null(tv_tol)) {
    erf1_inv_tv <- stats::qnorm((tv_tol + 1) / 2) / sqrt(2)
    c_tol       <- erf1_inv_tv * 2 * sqrt(2)
    m_pilot_from_gap <- .two_block_m_pilot_from_gap(rate, D_max, c_tol, m_min)
    m_convergence_pilot <- m_pilot_from_gap
  }

  pilot_plan <- .two_block_resolve_pilot_plan(
    is_gaussian         = is_gaussian,
    n                   = n,
    n_pilot_arg         = n_pilot_arg,
    gap_tol             = gap_tol,
    tv_tol              = tv_tol,
    m_convergence_user  = NULL,
    m_convergence_pilot = m_convergence_pilot,
    rate                = rate,
    p_dim               = p_dim,
    m_min               = m_min,
    any_non_normal      = any_non_normal
  )
  n_pilot          <- pilot_plan$n_pilot
  m_convergence    <- pilot_plan$m_convergence
  pilot_cost_opt   <- pilot_plan$pilot_cost_opt
  run_pilot        <- n_pilot > 0L
  run_ub           <- run_pilot && !is.null(tv_tol)

  if (is.null(m_min) && !run_pilot) {
    m_convergence <- 10L
  }

  if (is.null(rate_calibration) && !is.null(tv_tol)) {
    rate_calibration <- list(
      lambda_star = rate$lambda_star,
      eigenvalues = rate$eigenvalues,
      m_min       = m_min
    )
  }

  calib_label <- if (is_gaussian) {
    "exact (Gaussian posterior)"
  } else {
    sprintf("approximate (local-Gaussian at mode, %s)", family$family)
  }
  if (isTRUE(any_non_normal)) {
    calib_label <- paste0(
      calib_label,
      "; conservative: non-dNormal RE dispersion (disp_lower plug-in)"
    )
  }

  if (isTRUE(verbose) && !is.null(tv_tol)) {
    cat(sprintf(
      paste0(
        "--- %s: convergence calibration [%s]:\n",
        "    lambda* = %.4f, tv_tol = %g => m_min = %d (mode start), ",
        "main m_convergence = %d ---\n\n"
      ),
      engine_label, calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
    ))
    if (run_pilot && !is.null(mode_gap_max) && !is.null(m_pilot_from_gap)) {
      cat(sprintf(
        paste0(
          "--- %s: pilot sweep calibration [mode_gap_max = %g SD/dim, p = %d, ",
          "D_max = %.4f]:\n    m_min = %d, lambda* = %.4f => ",
          "m_convergence_pilot = %d ---\n\n"
        ),
        engine_label, mode_gap_max, p_dim, D_max, m_min,
        rate$lambda_star, m_convergence_pilot
      ))
    }
    if (run_pilot) {
      .two_block_print_pilot_plan(
        pilot_plan          = pilot_plan,
        n                   = n,
        m_convergence_pilot = m_convergence_pilot,
        rate                = rate,
        tv_tol              = tv_tol,
        p                   = p_dim,
        verbose             = verbose
      )
    }
  }

  method_label <- if (is_gaussian) "exact" else "local_gaussian_mode"
  if (isTRUE(any_non_normal)) {
    method_label <- paste0(method_label, "+disp_lower_bound")
  }

  convergence_info <- list(
    method              = method_label,
    tv_tol              = tv_tol,
    gap_tol             = gap_tol,
    n_pilot             = n_pilot,
    n_pilot_source      = pilot_plan$n_pilot_source,
    n_pilot_gap_tol     = pilot_plan$n_pilot_gap_tol,
    lambda_star         = rate$lambda_star,
    eigenvalues         = rate$eigenvalues,
    m_min               = m_min,
    m_certificate       = pilot_plan$m_certificate,
    m_convergence       = m_convergence,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    mode_gap_max        = if (run_pilot) mode_gap_max else NULL,
    m_pilot_from_gap    = if (run_pilot) m_pilot_from_gap else NULL,
    pilot_cost_opt      = pilot_cost_opt,
    draw_engine         = "rGLMM_sweep"
  )

  m_convergence_used <- m_convergence
  fixef_init         <- fixef_mode
  pilot_res          <- NULL
  pilot_chisq        <- NULL
  pilot_ub           <- NULL
  tau2_start_main    <- .two_block_tau2_start_from_pfamily(pfamily_list, re_names)

  if (run_pilot) {
    if (isTRUE(verbose)) {
      cat(sprintf(
        "--- %s [sweep-outer]: pilot stage (%d chains; m_convergence_pilot = %d) ---\n\n",
        engine_label, n_pilot, m_convergence_pilot
      ))
    }

    pilot_raw <- rGLMM_sweep(
      n_chains       = n_pilot,
      start_fixef    = fixef_mode,
      inner_sweeps   = m_convergence_pilot,
      design         = design,
      block1_prior   = dispprior_list,
      pfamily_list   = pfamily_list,
      family         = family,
      re_names       = re_names,
      group_levels   = group_levels,
      collect_block1 = collect_block1,
      progbar        = progbar_use,
      stage_label    = "pilot",
      fixef_mode     = fixef_mode_ref,
      b_mode         = b_mode_ref,
      b_start        = b_mode_ref,
      ptypes         = ptypes
    )

    pilot_chisq <- .two_block_pilot_chisq_test(
      fixef_draws = pilot_raw$fixef_draws,
      re_names    = re_names,
      fixef_mode  = fixef_mode,
      n_pilot     = n_pilot
    )

    if (isTRUE(stage_verbose) || isTRUE(verbose)) {
      cat(sprintf(
        "--- %s: pilot vs mode chi-squared test: p = %.4g (df = %d, n_pilot = %d) ---\n\n",
        engine_label,
        pilot_chisq$p_value, pilot_chisq$df, pilot_chisq$n_pilot
      ))
    }

    fixef_init <- .two_block_fixef_colmeans(
      pilot_raw$fixef_draws, re_names, fixef_mode
    )

    if (isTRUE(any_non_normal)) {
      tau2_start_main <- .two_block_tau2_start_from_dispersion_draws(
        pilot_raw$dispersion_fixef_draws, re_names
      )
    }

    if (run_ub) {
      pilot_ub <- .two_block_pilot_ub_from_coefficients(
        pilot_coefficients = pilot_raw$coefficients,
        n_pilot            = n_pilot,
        re_names           = re_names,
        group_levels       = group_levels,
        group_name         = group_name,
        x                  = D,
        group              = group,
        x_hyper            = W,
        prior_list         = dispprior_list,
        pfamily_list       = pfamily_list,
        family             = family,
        tv_tol             = tv_tol
      )
      if (pilot_ub$m_min_upper > m_convergence_used) {
        m_convergence_used <- pilot_ub$m_min_upper
      }
      convergence_info$lambda_star_upper <- pilot_ub$rate_upper$lambda_star
      convergence_info$eigenvalues_upper <- pilot_ub$max_eigenvalues
      convergence_info$m_min_upper       <- pilot_ub$m_min_upper
      convergence_info$i_max_rate        <- pilot_ub$i_max_rate
      convergence_info$lambda_star_vec   <- pilot_ub$lambda_star_vec
      convergence_info$m_convergence     <- m_convergence_used
    }

    if (isTRUE(stage_verbose) && run_ub) {
      .two_block_print_pilot_stage_diagnostics(
        n_pilot            = n_pilot,
        n_main             = n,
        pilot_ub           = pilot_ub,
        rate_calibration   = rate_calibration,
        m_convergence_used = m_convergence_used
      )
    } else if (isTRUE(verbose)) {
      cat(sprintf(
        "--- %s [sweep-outer]: pilot complete; main stage (%d chains; m_convergence = %d) ---\n\n",
        engine_label, n, m_convergence_used
      ))
    }

    pilot_res <- .rGLMM_format_v6_out(
      v6_out       = pilot_raw,
      n            = n_pilot,
      re_names     = re_names,
      group_levels = group_levels,
      fixef_mode   = fixef_mode,
      fixef_init   = fixef_mode
    )
  } else if (isTRUE(verbose)) {
    cat(sprintf(
      "--- %s [sweep-outer]: main stage (%d chains; m_convergence = %d) ---\n\n",
      engine_label, n, m_convergence_used
    ))
  }

  main_raw <- rGLMM_sweep(
    n_chains       = n,
    start_fixef    = fixef_init,
    inner_sweeps   = m_convergence_used,
    design         = design,
    block1_prior   = dispprior_list,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    stage_label    = "main",
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref,
    b_start        = b_mode_ref,
    ptypes         = ptypes,
    tau2_start     = tau2_start_main
  )

  draw_engine_args <- list(
    n_chains       = n,
    start_fixef    = fixef_init,
    inner_sweeps   = m_convergence_used,
    design         = design,
    block1_prior   = dispprior_list,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    stage_label    = "main",
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref,
    b_start        = b_mode_ref,
    ptypes         = ptypes,
    tau2_start     = tau2_start_main
  )

  main_res <- .rGLMM_format_v6_out(
    v6_out       = main_raw,
    n            = n,
    re_names     = re_names,
    group_levels = group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_init
  )

  .lmebayes_assemble_reg_result(
    staged              = main_res,
    call                = cl,
    m_convergence       = m_convergence_used,
    convergence_info    = convergence_info,
    pfamily_list        = pfamily_list,
    dispprior_list      = dispprior_list,
    family              = family,
    groupef.mode        = ranef_mode,
    any_non_normal      = pf_summary$any_non_normal,
    design              = design,
    result_class        = result_class,
    parent_class        = "rGLMM_reg",
    draw_engine         = "rGLMM_sweep",
    icm_info            = icm_info,
    pilot_draws         = if (run_pilot) pilot_res else NULL,
    n_pilot             = if (run_pilot) n_pilot else NULL,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    pilot_chisq         = if (run_pilot) pilot_chisq else NULL,
    pilot_ub            = if (run_ub) pilot_ub else NULL,
    tv_tol              = if (run_ub) tv_tol else NULL,
    offset              = inp$offset,
    weights             = if (is.null(inp$weights)) 1 else inp$weights
  )
}

#' Local-Gaussian rate at the ICM mode
#' @noRd
.rGLMM_rate_at_mode <- function(
    design,
    dispprior_list,
    pfamily_list,
    family,
    b_mode,
    group_levels,
    is_gaussian
) {
  if (is_gaussian) {
    two_block_rate_from_pfamily_list(
      x                 = design$D,
      group             = design$group,
      x_hyper           = design$W,
      prior_list_block1 = dispprior_list,
      pfamily_list      = pfamily_list,
      family            = gaussian(),
      group_levels      = group_levels
    )
  } else {
    mode_w <- two_block_mode_weights(
      x            = design$D,
      group        = design$group,
      b_mode       = b_mode,
      family       = family,
      group_levels = group_levels
    )
    two_block_rate_from_pfamily_list(
      x                 = design$D,
      group             = design$group,
      x_hyper           = design$W,
      prior_list_block1 = dispprior_list,
      pfamily_list      = pfamily_list,
      weights           = mode_w$weights,
      family            = family,
      group_levels      = group_levels
    )
  }
}

#' Format v6 batch output for staged \code{fixef.*} naming
#' @noRd
.rGLMM_format_v6_out <- function(
    v6_out,
    n,
    re_names,
    group_levels,
    fixef_mode,
    fixef_init
) {
  x <- list(
    fixef_draws            = v6_out$fixef_draws,
    coefficients           = v6_out$coefficients,
    dispersion_fixef_draws = v6_out$dispersion_fixef_draws,
    iters_fixef_draws      = v6_out$iters_fixef_draws,
    iters_ranef_draws      = v6_out$iters_ranef_draws,
    sweep_history          = v6_out$sweep_history,
    n                      = n
  )
  .two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}

#' @describeIn rGLMM_reg All population components \code{dNormal} (known
#'   \eqn{\tau^2_k}). Non-Gaussian: pilot always (unless \code{n_pilot = 0L});
#'   standard fixed-dispersion eigenvalue rate bounds (no ING
#'   \code{disp_lower} path).
#' @export
rGLMM_reg_known_vcov <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    icm_tol             = 1e-10,
    icm_maxit           = 200L,
    offset              = NULL,
    weights             = 1,
    family              = gaussian(),
    gap_tol             = 0.0196,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    Gridtype            = 2,
    n_envopt            = NULL,
    use_parallel        = TRUE,
    use_opencl          = FALSE,
    verbose             = FALSE,
    progbar             = FALSE,
    stage_verbose       = FALSE,
    rate_calibration    = NULL,
    collect_block1      = TRUE
) {
  cl <- match.call()
  fn_name <- "rGLMM_reg_known_vcov"

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = fn_name
  )

  inp <- .rGLMM_validate_matrix_inputs(
    n, y, D, group, W, tv_tol,
    group_name, family, mode_gap_max,
    gap_tol, pfamily_list, dispprior_list
  )
  if (!inp$pf_summary$all_dNormal) {
    stop(
      fn_name, "(): all population components must be dNormal(); ",
      "use rGLMM_reg_estimated_vcov() or rGLMM_reg().",
      call. = FALSE
    )
  }
  inp$offset <- offset
  inp$weights <- weights

  .rGLMM_reg_run(
    inp                = inp,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    mode_gap_max       = mode_gap_max,
    verbose            = verbose,
    progbar            = progbar,
    stage_verbose      = stage_verbose,
    rate_calibration   = rate_calibration,
    collect_block1     = collect_block1,
    engine_label       = fn_name,
    result_class       = "rGLMM_reg_known_vcov",
    cl                 = cl
  )
}

#' @describeIn rGLMM_reg ING population components (estimated \eqn{\tau^2_k}).
#'   Non-Gaussian: pilot always (unless \code{n_pilot = 0L}); conservative
#'   \code{disp_lower} eigenvalue bounds and pilot-updated \eqn{\tau^2} starts.
#' @export
rGLMM_reg_estimated_vcov <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    icm_tol             = 1e-10,
    icm_maxit           = 200L,
    offset              = NULL,
    weights             = 1,
    family              = gaussian(),
    gap_tol             = 0.0196,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    Gridtype            = 2,
    n_envopt            = NULL,
    use_parallel        = TRUE,
    use_opencl          = FALSE,
    verbose             = FALSE,
    progbar             = FALSE,
    stage_verbose       = FALSE,
    rate_calibration    = NULL,
    collect_block1      = TRUE
) {
  cl <- match.call()
  fn_name <- "rGLMM_reg_estimated_vcov"

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = fn_name
  )

  inp <- .rGLMM_validate_matrix_inputs(
    n, y, D, group, W, tv_tol,
    group_name, family, mode_gap_max,
    gap_tol, pfamily_list, dispprior_list
  )
  if (inp$pf_summary$all_dNormal) {
    stop(
      fn_name, "(): at least one population component must not be dNormal(); ",
      "use rGLMM_reg_known_vcov() or rGLMM_reg().",
      call. = FALSE
    )
  }
  inp$offset <- offset
  inp$weights <- weights

  .rGLMM_reg_run_with_pilot(
    inp                = inp,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    mode_gap_max       = mode_gap_max,
    verbose            = verbose,
    progbar            = progbar,
    stage_verbose      = stage_verbose,
    rate_calibration   = rate_calibration,
    collect_block1     = collect_block1,
    engine_label       = fn_name,
    result_class       = "rGLMM_reg_estimated_vcov",
    cl                 = cl
  )
}

#' @describeIn rGLMM_reg Route by population \code{pfamily_list} to known or
#'   estimated \eqn{\tau^2} simulation routes.
#' @export
rGLMM_reg <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    icm_tol             = 1e-10,
    icm_maxit           = 200L,
    offset              = NULL,
    weights             = 1,
    family              = gaussian(),
    gap_tol             = 0.0196,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    Gridtype            = 2,
    n_envopt            = NULL,
    use_parallel        = TRUE,
    use_opencl          = FALSE,
    verbose             = FALSE,
    progbar             = FALSE,
    stage_verbose       = FALSE,
    rate_calibration    = NULL,
    collect_block1      = TRUE
) {
  cl <- match.call()

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = "rGLMM_reg"
  )

  inp <- .rGLMM_validate_matrix_inputs(
    n, y, D, group, W, tv_tol,
    group_name, family, mode_gap_max,
    gap_tol, pfamily_list, dispprior_list
  )

  route_fn <- if (inp$pf_summary$all_dNormal) {
    rGLMM_reg_known_vcov
  } else {
    rGLMM_reg_estimated_vcov
  }
  mc <- match.call(expand.dots = FALSE)
  mc[[1L]] <- route_fn
  out <- eval(mc, parent.frame())
  out$call <- cl
  out
}
