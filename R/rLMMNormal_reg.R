#' Shared matrix-level validation for LMM replicate-chain engines
#' @noRd
.rLMM_validate_matrix_inputs <- function(
    n,
    y,
    x,
    x_hyper,
    tv_tol,
    m_convergence,
    re_coef_names,
    group_levels,
    group_name,
    block
) {
  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }

  if (!is.null(m_convergence)) {
    if (!is.numeric(m_convergence) || length(m_convergence) != 1L ||
        !is.finite(m_convergence) || m_convergence < 1) {
      stop("'m_convergence' must be NULL or a single integer >= 1.",
           call. = FALSE)
    }
    m_convergence <- as.integer(m_convergence)
  }

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  if (is.null(re_coef_names) || length(re_coef_names) != ncol(x)) {
    re_coef_names <- if (ncol(x) >= 1L) {
      cn <- colnames(x)
      if (is.null(cn) || length(cn) != ncol(x)) {
        paste0("RE", seq_len(ncol(x)))
      } else {
        cn
      }
    } else {
      stop("'x' must have at least one column.", call. = FALSE)
    }
  }
  colnames(x) <- re_coef_names
  re_names <- re_coef_names

  group_levels <- as.character(group_levels)
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    group_name <- tryCatch(
      deparse(substitute(block))[1L],
      error = function(e) "group"
    )
    if (!nzchar(group_name)) group_name <- "group"
  }

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop("length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
         call. = FALSE)
  }
  if (!setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }

  list(
    n             = n,
    y             = y,
    x             = x,
    x_hyper       = x_hyper,
    tv_tol        = tv_tol,
    m_convergence = m_convergence,
    re_names      = re_names,
    group_levels  = group_levels,
    group_name    = group_name
  )
}

#' @noRd
.rLMM_validate_P <- function(P, p_re, fn_name = "rLMMNormal_reg") {
  P <- as.matrix(P)
  if (!is.matrix(P) || nrow(P) != p_re || ncol(P) != p_re) {
    stop(
      fn_name, "(): 'P' must be a ", p_re, " x ", p_re, " matrix.",
      call. = FALSE
    )
  }
  if (!isSymmetric(P)) {
    stop(fn_name, "(): 'P' must be symmetric.", call. = FALSE)
  }
  P
}

#' @noRd
.rLMM_validate_fixed_dispersion_prior_list <- function(
    prior_list,
    fn_name = "rLMMNormal_reg"
) {
  if (!is.list(prior_list) || is.null(prior_list$dispersion)) {
    stop(
      fn_name, "(): 'prior_list' must contain 'dispersion' (fixed sigma^2).",
      call. = FALSE
    )
  }
  d <- prior_list$dispersion
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d <= 0) {
    stop(
      fn_name, "(): 'prior_list$dispersion' must be a single positive number.",
      call. = FALSE
    )
  }
  extra <- setdiff(
    names(prior_list),
    c("dispersion", "ddef")
  )
  if (length(extra)) {
    stop(
      fn_name, "(): 'prior_list' must contain fixed dispersion only; ",
      "unexpected fields: ", paste(extra, collapse = ", "), ".",
      call. = FALSE
    )
  }
  as.numeric(d)
}

#' @noRd
.rLMM_validate_dGamma_dispersion_prior_list <- function(
    prior_list,
    fn_name = "rLMMindepNormalGamma_reg"
) {
  if (!is.list(prior_list)) {
    stop(fn_name, "(): 'prior_list' must be a list.", call. = FALSE)
  }
  req <- c("shape", "rate", "beta", "Inv_Dispersion")
  miss <- req[!req %in% names(prior_list)]
  if (length(miss)) {
    stop(
      fn_name, "(): 'prior_list' must contain ",
      paste(req, collapse = ", "), " (from dGamma()).",
      call. = FALSE
    )
  }
  if (!isTRUE(prior_list$Inv_Dispersion)) {
    stop(
      fn_name, "(): dGamma() observation-dispersion prior requires ",
      "Inv_Dispersion = TRUE.",
      call. = FALSE
    )
  }
  shape <- prior_list$shape
  rate  <- prior_list$rate
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) ||
      shape <= 0) {
    stop(fn_name, "(): 'prior_list$shape' must be a positive scalar.",
         call. = FALSE)
  }
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) ||
      rate <= 0) {
    stop(fn_name, "(): 'prior_list$rate' must be a positive scalar.",
         call. = FALSE)
  }
  beta <- as.matrix(prior_list$beta)
  if (nrow(beta) != 1L || ncol(beta) != 1L) {
    stop(
      fn_name, "(): 'prior_list$beta' must be a 1 x 1 matrix for ",
      "observation-level dispersion.",
      call. = FALSE
    )
  }
  if (!is.null(prior_list$disp_lower) && !is.null(prior_list$disp_upper)) {
    if (prior_list$disp_upper <= prior_list$disp_lower) {
      stop(
        fn_name, "(): 'prior_list$disp_upper' must exceed 'disp_lower'.",
        call. = FALSE
      )
    }
  }
  prior_list
}

#' Observation-level linear predictor from group random effects
#' @noRd
.rLMM_observation_mu <- function(x, block, b_mat, group_levels) {
  g_chr <- as.character(block)
  g_idx <- match(g_chr, group_levels)
  if (anyNA(g_idx)) {
    stop("block levels not found in 'group_levels'.", call. = FALSE)
  }
  b_obs <- b_mat[g_idx, , drop = FALSE]
  rowSums(x * b_obs)
}

#' Group-level random-effect matrix from \code{coefficients} output
#' @noRd
.rLMM_b_matrix_from_coefficients <- function(
    coef_df,
    re_names,
    group_levels,
    group_name
) {
  b_mat <- matrix(
    NA_real_,
    nrow = length(group_levels),
    ncol = length(re_names),
    dimnames = list(group_levels, re_names)
  )
  for (lev in group_levels) {
    rows <- coef_df[[group_name]] == lev
    if (!any(rows)) {
      stop("coefficients missing group level: ", lev, call. = FALSE)
    }
    hit <- which(rows)[1L]
    b_mat[lev, ] <- as.numeric(coef_df[hit, re_names, drop = TRUE])
  }
  b_mat
}

#' ICM start for matrix-level LMM engines
#'
#' Iterated conditional modes for Block~2 hyperparameters at fixed Block~1
#' dispersion; used when \code{start = NULL} in \code{rLMMNormal_reg}.
#' @noRd
.rLMM_icm_at_start <- function(
    y,
    x,
    block,
    x_hyper,
    prior_list_block1,
    pfamily_list,
    re_names,
    group_levels,
    group_name,
    icm_tol,
    icm_maxit,
    verbose,
    engine_label
) {
  design_icm <- list(
    y             = y,
    Z             = x,
    groups        = factor(block, levels = group_levels),
    X_hyper       = x_hyper,
    re_coef_names = re_names,
    group_name    = group_name
  )
  icm <- .two_block_icm_at_start(
    design       = design_icm,
    prior_list   = prior_list_block1,
    pfamily_list = pfamily_list,
    re_names     = re_names,
    family       = gaussian(),
    tol          = icm_tol,
    maxit        = icm_maxit
  )
  if (isTRUE(verbose)) {
    cat(sprintf(
      "  %s: ICM posterior mean (converged: %s, %d iter, delta = %.2e)\n\n",
      engine_label,
      icm$icm$converged,
      icm$icm$iterations,
      icm$icm$delta
    ))
  }
  icm
}

#' Calibrate inner Gibbs sweeps for matrix-level LMM engines
#'
#' Computes the two-block convergence rate at the chain start, then sets
#' \code{m_convergence} to at least \code{l_for_tv(tv_tol) + 1} (Theorem 3).
#' @noRd
.rLMM_calibrate_m_convergence <- function(
    x,
    block,
    x_hyper,
    prior_list_block1,
    pfamily_list,
    group_levels,
    m_convergence,
    tv_tol,
    any_non_normal,
    engine_label,
    verbose
) {
  rate <- two_block_rate_v2(
    x                 = x,
    block             = block,
    x_hyper           = x_hyper,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    family            = gaussian(),
    group_levels      = group_levels
  )
  m_min <- two_block_l_for_tv(rate, tv_tol, method = "theorem3") + 1L
  if (is.null(m_convergence)) {
    m_convergence <- m_min
  } else if (m_convergence < m_min) {
    warning(
      engine_label, ": m_convergence = ", m_convergence,
      " is below the derived minimum m_min = ", m_min,
      " for tv_tol = ", tv_tol, "; using m_min instead.",
      call. = FALSE
    )
    m_convergence <- m_min
  }
  calib_label <- if (isTRUE(any_non_normal)) {
    "conservative: non-dNormal RE dispersion (disp_lower plug-in)"
  } else {
    "exact"
  }
  if (isTRUE(verbose)) {
    cat(sprintf(
      paste0(
        "--- %s: convergence calibration [%s]: lambda* = %.4f, ",
        "tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n"
      ),
      engine_label, calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
    ))
  }
  list(
    m_convergence     = m_convergence,
    convergence_info  = list(
      method        = if (isTRUE(any_non_normal)) "disp_lower_bound" else "exact",
      tv_tol        = tv_tol,
      lambda_star   = rate$lambda_star,
      eigenvalues   = rate$eigenvalues,
      m_min         = m_min,
      m_convergence = m_convergence
    )
  )
}

#' Shared sampling pipeline for \code{rLMMNormal_reg_*_vcov} routes
#' @noRd
.rLMMNormal_reg_run <- function(
    inp,
    block,
    P,
    dispersion,
    pfamily_list,
    pf_summary,
    start,
    icm_tol,
    icm_maxit,
    progbar,
    verbose,
    engine_label,
    any_non_normal,
    draw_engine,
    result_class,
    cl
) {
  prior_list_block1 <- list(
    P          = P,
    dispersion = dispersion,
    ddef       = FALSE
  )
  .two_block_validate_block1_prior(prior_list_block1, family = gaussian())

  icm_info   <- NULL
  ranef_mode <- NULL

  if (is.null(start)) {
    icm <- .rLMM_icm_at_start(
      y                 = inp$y,
      x                 = inp$x,
      block             = block,
      x_hyper           = inp$x_hyper,
      prior_list_block1 = prior_list_block1,
      pfamily_list      = pfamily_list,
      re_names          = inp$re_names,
      group_levels      = inp$group_levels,
      group_name        = inp$group_name,
      icm_tol           = icm_tol,
      icm_maxit         = icm_maxit,
      verbose           = verbose,
      engine_label      = engine_label
    )
    start      <- icm$start
    ranef_mode <- icm$b_start
    icm_info   <- icm$icm
  } else {
    if (!is.list(start) || is.null(names(start))) {
      stop("'start' must be a named list or NULL.", call. = FALSE)
    }
    if (!setequal(names(start), inp$re_names)) {
      stop("names(start) must match re_coef_names.", call. = FALSE)
    }
    start <- start[inp$re_names]
  }
  fixef_mode <- start

  calib <- .rLMM_calibrate_m_convergence(
    x                 = inp$x,
    block             = block,
    x_hyper           = inp$x_hyper,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    group_levels      = inp$group_levels,
    m_convergence     = inp$m_convergence,
    tv_tol            = inp$tv_tol,
    any_non_normal    = any_non_normal,
    engine_label      = engine_label,
    verbose           = verbose
  )
  m_convergence    <- calib$m_convergence
  convergence_info <- calib$convergence_info
  convergence_info$draw_engine <- draw_engine

  out <- two_block_rNormal_reg_v2(
    n                 = inp$n,
    y                 = inp$y,
    x                 = inp$x,
    block             = block,
    x_hyper           = inp$x_hyper,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    fixef_start       = start,
    re_coef_names     = inp$re_names,
    group_levels      = inp$group_levels,
    group_name        = inp$group_name,
    family            = gaussian(),
    m_convergence     = m_convergence,
    progbar           = progbar
  )

  staged <- .rLMM_format_v2_out(
    v2_out       = out,
    n            = inp$n,
    re_names     = inp$re_names,
    group_levels = inp$group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_mode
  )

  staged$call             <- cl
  staged$m_convergence    <- m_convergence
  staged$convergence_info <- convergence_info
  staged$draw_engine      <- draw_engine
  staged$pfamily_list     <- pfamily_list
  staged$prior_list       <- prior_list_block1
  staged$family           <- gaussian()
  staged$ranef.mode       <- ranef_mode
  staged$icm_info         <- icm_info
  staged$ptypes           <- pf_summary$ptypes
  staged$any_non_normal   <- any_non_normal

  class(staged) <- c(result_class, "rLMMNormal_reg", "list")
  staged
}

#' Replicate-chain Gibbs sampling for Bayesian LMMs with fixed observation dispersion
#'
#' Matrix-level LMM sampler with fixed observation-level \eqn{\sigma^2} and
#' random-effect prior precision \code{P}.  Validates inputs and dispatches to
#' \code{\link{rLMMNormal_reg_known_vcov}} or
#' \code{\link{rLMMNormal_reg_estimated_vcov}} according to Block~2
#' \code{pfamily_list}.
#'
#' @details
#' Steps 1--3 (validation and Block~2 classification) run here; sampling runs
#' in the route-specific functions. See their documentation for the full
#' pipeline (ICM start, convergence calibration, \code{two_block_rNormal_reg_v2},
#' output formatting).
#'
#' Each stored draw runs \code{m_convergence} inner Gibbs sweeps of the
#' two-block sampler: Block~1 updates group random effects \eqn{b} given
#' Block~2 hyperparameters \eqn{\gamma}; Block~2 updates \eqn{\gamma} (and
#' optional ING \eqn{\tau^2_k}) given \eqn{b}.
#'
#' @param n Number of stored draws. If \code{length(n) > 1}, the length is
#'   taken to be the number required.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Level-1 design matrix \code{Z} (\code{l2 x p_re}).
#' @param block Grouping factor or block partition of length \code{l2}.
#' @param x_hyper Named list of group-level design matrices \code{X_k}
#'   (\code{J x q_k}), one per column of \code{x}.
#' @param P Random-effect prior precision matrix (\code{p_re x p_re}).
#' @param prior_list List containing fixed observation dispersion only:
#'   \code{list(dispersion = sigma2)}.
#' @param pfamily_list Named list of \code{pfamily} objects for Block~2.
#' @param start Named list of Block~2 hyper-parameter vectors at which each
#'   replicate chain is initialised (typically the ICM posterior mean). When
#'   \code{NULL} (default), the ICM posterior mean is computed internally via
#'   \code{\link{lmerb_posterior_mean}}.
#' @param icm_tol,icm_maxit Convergence controls passed to
#'   \code{lmerb_posterior_mean()} when \code{start = NULL}.
#' @param m_convergence Inner Gibbs steps per stored draw. When \code{NULL},
#'   derived from Theorem~3 at \code{start} using \code{tv_tol}.
#' @param tv_tol Total variation tolerance in \code{(0, 1)} for convergence
#'   calibration.
#' @param re_coef_names Character vector naming columns of \code{x}.
#' @param group_levels Character vector defining row order of Block~1 draws.
#' @param group_name Name for the grouping column in \code{coefficients}.
#' @param progbar Logical; show a text progress bar during sampling.
#' @param verbose Print the convergence calibration line.
#' @return Object of class \code{c("rLMMNormal_reg", "list")} (plus a route-specific
#'   prefix class).  Includes \code{ptypes} and \code{any_non_normal}.
#' @family simfuncs
#' @seealso \code{\link{rLMMNormal_reg_known_vcov}},
#'   \code{\link{rLMMNormal_reg_estimated_vcov}},
#'   \code{\link{rLMMindepNormalGamma_reg}}, \code{\link{rGLMM}}
#' @export
rLMMNormal_reg <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    start           = NULL,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    m_convergence   = NULL,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()

  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol, m_convergence,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names))
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(prior_list)
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)

  route_fn <- if (pf_summary$all_dNormal) {
    rLMMNormal_reg_known_vcov
  } else {
    rLMMNormal_reg_estimated_vcov
  }
  mc <- match.call(expand.dots = FALSE)
  mc[[1L]] <- route_fn
  out <- eval(mc, parent.frame())
  out$call <- cl
  out
}

#' LMM replicate chains with known random-effect covariance (all \code{dNormal})
#'
#' Block~2 route used when every component of \code{pfamily_list} is
#' \code{dNormal}: component scales \eqn{\tau^2_k} are fixed, so the Block~1
#' prior covariance \code{Sigma_ranef} is known for sampling.  Convergence
#' calibration uses the exact two-block rate (Theorem~3).
#'
#' @inheritParams rLMMNormal_reg
#' @return Object of class
#'   \code{c("rLMMNormal_reg_known_vcov", "rLMMNormal_reg", "list")}.
#' @family simfuncs
#' @seealso \code{\link{rLMMNormal_reg}}, \code{\link{rLMMNormal_reg_estimated_vcov}}
#' @export
rLMMNormal_reg_known_vcov <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    start           = NULL,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    m_convergence   = NULL,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMNormal_reg_known_vcov"

  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol, m_convergence,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names), fn_name = fn_name)
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(
    prior_list, fn_name = fn_name
  )
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (!pf_summary$all_dNormal) {
    stop(
      fn_name, "(): all Block~2 components must be dNormal(); ",
      "use rLMMNormal_reg_estimated_vcov() or rLMMNormal_reg().",
      call. = FALSE
    )
  }

  .rLMMNormal_reg_run(
    inp              = inp,
    block            = block,
    P                = P,
    dispersion       = dispersion,
    pfamily_list     = pfamily_list,
    pf_summary       = pf_summary,
    start            = start,
    icm_tol          = icm_tol,
    icm_maxit        = icm_maxit,
    progbar          = progbar,
    verbose          = verbose,
    engine_label     = fn_name,
    any_non_normal   = FALSE,
    draw_engine      = "two_block_rNormal_reg_v2_known_vcov",
    result_class     = "rLMMNormal_reg_known_vcov",
    cl               = cl
  )
}

#' LMM replicate chains with estimated random-effect covariance scales
#'
#' Block~2 route used when at least one \code{pfamily_list} component is not
#' \code{dNormal} (typically \code{dIndependent_Normal_Gamma}): component
#' scales \eqn{\tau^2_k} are updated during Gibbs sampling, so the RE
#' covariance is estimated.  Convergence calibration uses a conservative
#' \code{disp_lower} plug-in rate bound.
#'
#' @param gap_tol Legacy mode--mean gap tolerance used to derive the pilot
#'   chain count when \code{tv_tol} is \code{NULL}. Default \code{0.0196}.
#'   Set \code{NULL} to skip the pilot unless \code{tv_tol} is set.
#' @param mode_gap_max When the pilot runs and \code{tv_tol} is set, used to
#'   derive pilot inner sweeps from the mode--mean gap bound.
#' @param diag_sweeps When \code{TRUE}, print one combined Block~2 chain-mean
#'   table per stage when that stage finishes. Sweep history is always stored
#'   on \code{$sweep_history} (and \code{$pilot$sweep_history} when a pilot
#'   runs).
#' @param stage_verbose Print pilot chi-squared and post-pilot UB diagnostics
#'   between stages (in addition to \code{verbose} stage headers).
#' @inheritParams rLMMNormal_reg
#' @return Object of class
#'   \code{c("rLMMNormal_reg_estimated_vcov", "rLMMNormal_reg", "list")}.
#'   Includes optional \code{pilot}, \code{pilot_chisq}, and
#'   \code{sweep_history} when a pilot stage runs (sweep-outer engine).
#' @family simfuncs
#' @seealso \code{\link{rLMMNormal_reg}}, \code{\link{rLMMNormal_reg_known_vcov}},
#'   \code{\link{rGLMM_sweep}}
#' @export
rLMMNormal_reg_estimated_vcov <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    start           = NULL,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    m_convergence   = NULL,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE,
    gap_tol         = 0.0196,
    mode_gap_max    = 1.0,
    diag_sweeps     = FALSE,
    stage_verbose   = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMNormal_reg_estimated_vcov"

  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol, m_convergence,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names), fn_name = fn_name)
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(
    prior_list, fn_name = fn_name
  )
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (pf_summary$all_dNormal) {
    stop(
      fn_name, "(): at least one Block~2 component must not be dNormal(); ",
      "use rLMMNormal_reg_known_vcov() or rLMMNormal_reg().",
      call. = FALSE
    )
  }

  .rLMMNormal_reg_run_with_pilot(
    inp            = inp,
    block          = block,
    P              = P,
    dispersion     = dispersion,
    pfamily_list   = pfamily_list,
    pf_summary     = pf_summary,
    start          = start,
    icm_tol        = icm_tol,
    icm_maxit      = icm_maxit,
    progbar        = progbar,
    verbose        = verbose,
    stage_verbose  = stage_verbose,
    gap_tol        = gap_tol,
    mode_gap_max   = mode_gap_max,
    diag_sweeps    = diag_sweeps,
    engine_label   = fn_name,
    result_class   = "rLMMNormal_reg_estimated_vcov",
    cl             = cl
  )
}

#' Replicate-chain Gibbs sampling for Bayesian LMMs with dGamma observation dispersion
#'
#' Matrix-level LMM engine scaffold for a \code{dGamma()} prior on observation-level
#' dispersion (precision), parallel to \code{\link{rindepNormalGamma_reg}} at the
#' GLM level.  \code{prior_list} must be the \code{prior_list} field extracted
#' from a \code{dGamma(Inv_Dispersion = TRUE)} \code{pfamily}.
#'
#' @param P Random-effect prior precision matrix (\code{p_re x p_re}).
#' @param prior_list \code{dGamma()} prior specification (\code{shape}, \code{rate},
#'   \code{beta}, \code{Inv_Dispersion}, optional truncation bounds).
#' @param dispersion_fix Plug-in \eqn{\sigma^2} for ICM and TV calibration.
#' @inheritParams rLMMNormal_reg
#' @return Object of class \code{c("rLMMindepNormalGamma_reg", "list")}.
#' @family simfuncs
#' @seealso \code{\link{rLMMNormal_reg}}, \code{\link{rindepNormalGamma_reg}}
#' @export
rLMMindepNormalGamma_reg <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    P,
    prior_list,
    pfamily_list,
    dispersion_fix,
    start           = NULL,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    m_convergence   = NULL,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()
  inp <- .rLMM_validate_matrix_inputs(
    n, y, x, x_hyper, tv_tol, m_convergence,
    re_coef_names, group_levels, group_name, block
  )
  P <- .rLMM_validate_P(P, length(inp$re_names))
  prior_list <- .rLMM_validate_dGamma_dispersion_prior_list(prior_list)

  if (is.null(dispersion_fix) || !is.numeric(dispersion_fix) ||
      length(dispersion_fix) != 1L || !is.finite(dispersion_fix) ||
      dispersion_fix <= 0) {
    stop(
      "'dispersion_fix' must be a single positive number (plug-in sigma^2).",
      call. = FALSE
    )
  }
  dispersion_fix <- as.numeric(dispersion_fix)

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)

  prior_list_block1_cal <- list(
    P          = P,
    dispersion = dispersion_fix,
    ddef       = FALSE
  )
  .two_block_validate_block1_prior(prior_list_block1_cal, family = gaussian())

  fixef_start <- start
  icm_info    <- NULL
  ranef_mode  <- NULL
  if (is.null(fixef_start)) {
    icm <- .rLMM_icm_at_start(
      y                 = inp$y,
      x                 = inp$x,
      block             = block,
      x_hyper           = inp$x_hyper,
      prior_list_block1 = prior_list_block1_cal,
      pfamily_list      = pfamily_list,
      re_names          = inp$re_names,
      group_levels      = inp$group_levels,
      group_name        = inp$group_name,
      icm_tol           = icm_tol,
      icm_maxit         = icm_maxit,
      verbose           = verbose,
      engine_label      = "rLMMindepNormalGamma_reg"
    )
    fixef_start <- icm$start
    ranef_mode  <- icm$b_start
    icm_info    <- icm$icm
  } else {
    if (!is.list(fixef_start) || is.null(names(fixef_start))) {
      stop("'start' must be a named list or NULL.", call. = FALSE)
    }
    if (!setequal(names(fixef_start), inp$re_names)) {
      stop("names(start) must match re_coef_names.", call. = FALSE)
    }
    fixef_start <- fixef_start[inp$re_names]
  }

  calib <- .rLMM_calibrate_m_convergence(
    x                 = inp$x,
    block             = block,
    x_hyper           = inp$x_hyper,
    prior_list_block1 = prior_list_block1_cal,
    pfamily_list      = pfamily_list,
    group_levels      = inp$group_levels,
    m_convergence     = inp$m_convergence,
    tv_tol            = inp$tv_tol,
    any_non_normal    = pf_summary$any_non_normal,
    engine_label      = "rLMMindepNormalGamma_reg",
    verbose           = verbose
  )
  m_convergence    <- calib$m_convergence
  convergence_info <- calib$convergence_info
  convergence_info$draw_engine <- "rLMMindepNormalGamma_reg_outer"

  if (is.null(ranef_mode)) {
    ranef_mode <- matrix(
      0,
      nrow = length(inp$group_levels),
      ncol = length(inp$re_names),
      dimnames = list(inp$group_levels, inp$re_names)
    )
  }
  b_mat <- ranef_mode

  n_obs <- length(inp$y)
  x_disp <- matrix(1, nrow = n_obs, ncol = 1L)
  wt     <- rep(1, n_obs)

  fixef_cur <- fixef_start
  re_names  <- inp$re_names
  n         <- inp$n

  fixef_draws <- stats::setNames(
    lapply(re_names, function(k) {
      q_k <- length(fixef_cur[[k]])
      matrix(NA_real_, nrow = n, ncol = q_k,
             dimnames = list(NULL, names(fixef_cur[[k]])))
    }),
    re_names
  )
  dispersion_ranef_draws <- numeric(n)
  coef_parts             <- vector("list", n)
  disp_fixef_draws       <- NULL
  iters_fixef_draws      <- NULL

  if (isTRUE(progbar) && n > 1L) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(n)) {
    mu_hat <- .rLMM_observation_mu(
      inp$x, block, b_mat, inp$group_levels
    )
    gamma_out <- rGamma_reg(
      n           = 1L,
      y           = inp$y,
      x           = x_disp,
      prior_list  = prior_list,
      offset      = mu_hat,
      weights     = wt,
      family      = gaussian(),
      progbar     = FALSE,
      verbose     = FALSE
    )
    sigma2_i <- as.numeric(gamma_out$dispersion[1L])
    dispersion_ranef_draws[i] <- sigma2_i

    lmm_i <- rLMMNormal_reg(
      n               = 1L,
      y               = inp$y,
      x               = inp$x,
      block           = block,
      x_hyper         = inp$x_hyper,
      P               = P,
      prior_list      = list(dispersion = sigma2_i),
      pfamily_list    = pfamily_list,
      start           = fixef_cur,
      icm_tol         = icm_tol,
      icm_maxit       = icm_maxit,
      m_convergence   = m_convergence,
      tv_tol          = inp$tv_tol,
      re_coef_names   = re_names,
      group_levels    = inp$group_levels,
      group_name      = inp$group_name,
      progbar         = FALSE,
      verbose         = FALSE
    )

    for (k in re_names) {
      fixef_draws[[k]][i, ] <- lmm_i$fixef[[k]][1L, ]
    }
    coef_parts[[i]] <- lmm_i$coefficients
    fixef_cur <- lapply(lmm_i$fixef, function(m) {
      stats::setNames(m[1L, ], colnames(m))
    })
    b_mat <- .rLMM_b_matrix_from_coefficients(
      lmm_i$coefficients,
      re_names,
      inp$group_levels,
      inp$group_name
    )
    if (i == 1L) {
      disp_fixef_draws  <- lmm_i$fixef.dispersion
      iters_fixef_draws <- lmm_i$fixef.iters
    } else {
      disp_fixef_draws  <- rbind(disp_fixef_draws, lmm_i$fixef.dispersion)
      iters_fixef_draws <- rbind(iters_fixef_draws, lmm_i$fixef.iters)
    }

    if (isTRUE(progbar) && n > 1L) {
      utils::setTxtProgressBar(pb, i / n)
    }
  }

  coefficients <- do.call(rbind, coef_parts)
  rownames(coefficients) <- NULL

  v2_out <- list(
    fixef_draws            = fixef_draws,
    coefficients           = coefficients,
    dispersion_fixef_draws = disp_fixef_draws,
    iters_fixef_draws      = iters_fixef_draws,
    mu_all_last            = build_mu_all(
      list(
        X_hyper       = inp$x_hyper,
        re_coef_names = re_names,
        groups        = factor(block, levels = inp$group_levels)
      ),
      fixef_cur,
      group_levels = inp$group_levels
    )$mu_all
  )

  staged <- .rLMM_format_v2_out(
    v2_out       = v2_out,
    n            = n,
    re_names     = re_names,
    group_levels = inp$group_levels,
    fixef_mode   = fixef_start,
    fixef_init   = fixef_start
  )

  staged$call                  <- cl
  staged$m_convergence         <- m_convergence
  staged$convergence_info      <- convergence_info
  staged$draw_engine           <- "rLMMindepNormalGamma_reg_outer"
  staged$pfamily_list          <- pfamily_list
  staged$prior_list            <- prior_list
  staged$dispersion_prior_list <- prior_list
  staged$dispersion_ranef      <- dispersion_ranef_draws
  staged$dispersion_ranef.mean <- mean(dispersion_ranef_draws)
  staged$family                <- gaussian()
  staged$ranef.mode            <- ranef_mode
  staged$icm_info              <- icm_info
  staged$P                     <- P
  staged$dispersion_fix          <- dispersion_fix
  staged$ptypes         <- pf_summary$ptypes
  staged$any_non_normal <- pf_summary$any_non_normal

  class(staged) <- c("rLMMindepNormalGamma_reg", "list")
  staged
}

#' Format v2 batch output for staged \code{fixef.*} naming
#' @noRd
.rLMM_format_v2_out <- function(
    v2_out,
    n,
    re_names,
    group_levels,
    fixef_mode,
    fixef_init
) {
  x <- list(
    fixef_draws            = v2_out$fixef_draws,
    coefficients           = v2_out$coefficients,
    dispersion_fixef_draws = v2_out$dispersion_fixef_draws,
    iters_fixef_draws      = v2_out$iters_fixef_draws,
    mu_all_last            = v2_out$mu_all_last,
    re_coef_names          = re_names,
    group_levels           = group_levels,
    n                      = n
  )
  .two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}
