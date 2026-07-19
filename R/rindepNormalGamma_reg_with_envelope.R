#' @noRd
.rindepNormalGamma_reg_impl <- function(
    n,
    y,
    x,
    prior_list,
    offset = NULL,
    weights = 1,
    family = gaussian(),
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    progbar = TRUE,
    return_envelope = FALSE,
    call = NULL
) {
  if (is.null(call)) {
    call <- sys.call()
  }
  
  offset2=offset
  wt=weights
  
  if(length(wt)==1) wt=rep(wt,length(y))
  
  ### Initial implementation of Likelihood subgradient Sampling 
  ### Currently uses as single point for conditional tangencis
  ### (at conditional posterior modes)
  ### Verify this yields correct results and then try to implement grid approach
  
  ## Use the prior list to set the prior elements if it is not missing
  ## Error checking to verify that the correct elements are present
  ## Shold be implemented
  
  
  ## Step 1: Validate Prior Specification
  
  if(missing(prior_list)) stop("Prior Specification Missing")
  if(!missing(prior_list)){
    if(!is.null(prior_list$mu)) mu=prior_list$mu
    if(!is.null(prior_list$Sigma)) Sigma=prior_list$Sigma
    if(!is.null(prior_list$dispersion)) dispersion=prior_list$dispersion
    else dispersion=NULL
    if(!is.null(prior_list$shape)) shape=prior_list$shape
    else shape=NULL
    if(!is.null(prior_list$rate)) rate=prior_list$rate
    else rate=NULL
    if (!is.null(prior_list$max_disp_perc)) {
      max_disp_perc <- prior_list$max_disp_perc
    } else {
      max_disp_perc <- 0.99
    }
    
    ## New: extract optional low/upp from prior_list
    if (!is.null(prior_list$disp_lower))  disp_lower <- prior_list$disp_lower  else disp_lower <- NULL
    if (!is.null(prior_list$disp_upper))  disp_upper <- prior_list$disp_upper  else disp_upper <- NULL
    
    ## Validation if both are provided
    if (!is.null(disp_lower) && !is.null(disp_upper)) {
      if (!is.numeric(disp_lower) || !is.numeric(disp_upper)) {
        stop("prior_list$disp_lower and prior_list$disp_upper must be numeric.")
      }
      if (disp_lower <= 0 || disp_upper <= 0) {
        stop("prior_list$disp_lower and prior_list$disp_upper must be positive.")
      }
      if (disp_upper <= disp_lower) {
        stop("prior_list$disp_upper must be strictly greater than prior_list$disp_lower.")
      }
    }
    
  }
  

  
  # Reconstruct P from Sigma
  R <- chol(Sigma)
  P <- chol2inv(R)
  P <- 0.5 * (P + t(P))
  
  
  ##########################  BEGIN *.CPP  MIGRATION   #########################################################

  ## --- NEW: Normalize and validate inputs before calling Rcore ---
  
  # Coerce basic types
  y  <- as.numeric(y)
  x  <- as.matrix(x)
  mu <- as.numeric(mu)
  wt <- as.numeric(wt)
  
  # Normalize offset
  if (is.null(offset2)) offset2 <- rep(0, length(y))
  offset2 <- as.numeric(offset2)
  
  # Normalize weights
  if (length(wt) == 1L) wt <- rep(wt, length(y))
  stopifnot(length(wt) == length(y))
  
  # Dimension checks
  stopifnot(nrow(x) == length(y))
  stopifnot(length(mu) == ncol(x))
  
  ## Prior-vs-data balance guard. The dispersion envelope caps the log-tilt
  ## at n_w/2 (the data contribution to shape2; see Remark 4.1.3 of the ING
  ## vignette), which presumes a likelihood-dominated regime. Under the
  ## Prior_Setup ING calibration shape = (n_prior + 1 + p)/2, so the implied
  ## effective prior sample size must not exceed the weighted observation
  ## count n_w = sum(wt) (equivalently pwt <= 0.5).
  if (!is.null(shape)) {
    .ing_stop_if_prior_exceeds_data(
      shape       = shape,
      p           = ncol(x),
      n_w         = sum(wt),
      detail      = paste0(
        "the data supply only n_w = sum(weights) = ",
        signif(sum(wt), 4)
      ),
      limit_label = "n_w"
    )
  }
  
  # Reconstruct P from Sigma and enforce SPD
  R    <- chol(Sigma)
  Pinv <- chol2inv(R)
  P    <- 0.5 * (Pinv + t(Pinv))
  
  stopifnot(isSymmetric(P))
  
  tol <- 1e-6
  ev  <- eigen(P, symmetric = TRUE)$values
  stopifnot(all(ev >= -tol * abs(ev[1L])))
  
  ## NOTE: this used to require a (numerically) Zellner g-prior for the
  ## coefficient covariance, working around a gap in Chapter A07's Claim 7 /
  ## Remark 5.5.7 (endpoint-only minimization of UB2_j(d) is only exact when
  ## K = Q^{-1/2} P Q^{-1/2} is isotropic). That guard has been removed: the
  ## gap is now fixed directly via exact root-finding for UB2_Min_j in
  ## src/EnvelopeDispersionBuild.cpp::bound_ub2_over_dispersion(), which is
  ## correct for anisotropic priors too. See
  ## data-raw/README_ub2_rootfinding_fix.md.
  
  # dispersion must be numeric scalar or NULL
  if (!is.null(dispersion)) {
    dispersion <- as.numeric(dispersion)
    stopifnot(length(dispersion) == 1L, is.finite(dispersion))
  }
  
  if (is.null(n_envopt)) n_envopt <- n
  n_envopt <- as.integer(n_envopt)

  .lmebayes_check_disp_bounds_or_stop(
    disp_lower, disp_upper, "rindepNormalGamma_reg (pre-.Call)"
  )

  core_out <- if (isTRUE(return_envelope)) {
    .rIndepNormalGammaReg_with_envelope_cpp(
      n,
      y,
      x,
      mu,
      P,
      offset2,
      wt,
      shape,
      rate,
      max_disp_perc,
      disp_lower,
      disp_upper,
      Gridtype,
      n_envopt,
      use_parallel,
      use_opencl,
      verbose,
      progbar
    )
  } else {
    .rIndepNormalGammaReg_cpp(
      n,
      y,
      x,
      mu,
      P,
      offset2,
      wt,
      shape,
      rate,
      max_disp_perc,
      disp_lower,
      disp_upper,
      Gridtype,
      n_envopt,
      use_parallel,
      use_opencl,
      verbose,
      progbar
    )
  }
  


  
  out        <- core_out$out
  betastar   <- core_out$betastar
  disp_out   <- core_out$disp_out
  iters_out  <- core_out$iters_out
  weight_out <- core_out$weight_out
  low        <- core_out$low
  upp        <- core_out$upp
  
  famfunc=glmbayesCore::glmbfamfunc(gaussian())  
  f1=famfunc$f1
  
  R <- chol(Sigma)
  Prec <- chol2inv(R)
  Prec <- 0.5 * (Prec + t(Prec))   # enforce symmetry

  pfamily_obj <- list(
    pfamily = "dIndependent_Normal_Gamma",
    prior_list = list(
      mu = mu,
      Sigma = Sigma,
      dispersion = dispersion,
      shape = shape,
      rate = rate,
      max_disp_perc = max_disp_perc,
      disp_lower = low,
      disp_upper = upp
    )
  )
  attr(pfamily_obj, "Prior Type") <- "dIndependent_Normal_Gamma"
  class(pfamily_obj) <- "pfamily"

  if (isTRUE(return_envelope)) {
    outlist <- list(
      coefficients = t(out),
      coef.mode = betastar,
      dispersion = disp_out,
      Prior = list(
        mean = mu, Sigma = Sigma, shape = shape, rate = rate, Precision = Prec
      ),
      family = gaussian(),
      prior.weights = wt,
      y = y,
      x = x,
      call = call,
      famfunc = famfunc,
      iters = iters_out,
      Envelope = core_out$Env,
      loglike = NULL,
      weight_out = weight_out,
      sim_bounds = list(low = low, upp = upp),
      gamma_list = core_out$gamma_list,
      UB_list = core_out$UB_list,
      diagnostics = core_out$diagnostics
    )
  } else {
    ## Match glmbayesCore::rindepNormalGamma_reg return structure exactly.
    outlist <- list(
      coefficients = t(out),
      coef.mode = betastar,
      dispersion = disp_out,
      Prior = list(
        mean = mu, Sigma = Sigma, shape = shape, rate = rate, Precision = Prec
      ),
      family = gaussian(),
      prior.weights = wt,
      y = y,
      x = x,
      call = call,
      famfunc = famfunc,
      iters = iters_out,
      Envelope = NULL,
      loglike = NULL,
      weight_out = weight_out,
      sim_bounds = list(low = low, upp = upp)
    )
  }

  outlist$pfamily <- pfamily_obj
  
  colnames(outlist$coefficients)<-colnames(x)
  outlist$offset2<-offset2
  class(outlist)<-c(outlist$class,"rglmb")
  
  return(outlist)
}


#' Independent Normal--Gamma regression with envelope artifacts returned
#'
#' Diagnostic wrapper with the same sampling pipeline as
#' \code{glmbayesCore::rindepNormalGamma_reg()}, but also returns the
#' standardized envelope (\code{Envelope}), \code{gamma_list}, \code{UB_list},
#' and \code{diagnostics} used by the joint accept--reject sampler
#' (standardized subproblem; see \code{glmbayesCore::EnvelopeOrchestrator()}).
#' Intended for parity checks and development; production callers should use
#' \code{glmbayesCore::rindepNormalGamma_reg()}.
#'
#' @param n Number of draws to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param y A vector of observations of length \code{m}.
#' @param x A design matrix of dimension \code{m * p}.
#' @param prior_list A list with prior parameters (e.g., shape, rate, beta) used in the simulation.
#' @param offset Optional numeric vector of length \code{m} specifying known components of the linear predictor.
#' @param weights Optional numeric vector of prior weights.
#' @param family A description of the error distribution and link function (see \code{\link{family}}).
#' @param Gridtype Optional integer specifying the method used to construct the envelope function.
#' @param n_envopt Effective sample size passed to \code{EnvelopeOpt} for grid
#'   construction. Defaults to match \code{n}. Larger values encourage tighter
#'   envelopes.
#' @param use_parallel Logical. Whether to use parallel processing.
#' @param use_opencl Logical. Whether to use OpenCL acceleration.
#' @param verbose Logical. Whether to print progress messages.
#' @param progbar Logical. Whether to display a progress bar during simulation.
#'
#' @return An \code{rglmb} object like \code{glmbayesCore::rindepNormalGamma_reg()},
#'   plus \code{Envelope}, \code{gamma_list}, \code{UB_list}, and \code{diagnostics}.
#' @export
#' @rdname rindepNormalGamma_reg_with_envelope
rindepNormalGamma_reg_with_envelope <- function(
    n,
    y,
    x,
    prior_list,
    offset = NULL,
    weights = 1,
    family = gaussian(),
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    progbar = TRUE
) {
  .rindepNormalGamma_reg_impl(
    n = n,
    y = y,
    x = x,
    prior_list = prior_list,
    offset = offset,
    weights = weights,
    family = family,
    Gridtype = Gridtype,
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose,
    progbar = progbar,
    return_envelope = TRUE,
    call = match.call()
  )
}
