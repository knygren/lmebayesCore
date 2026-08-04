## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## Post-hoc, end-of-simulation "revised" total-variation bound for the ING
## known_vcov engine, per inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 17:
##
##   A   = { draws where every group's Lambda + H_j(beta_j) is PD AND the
##           marginal lambda_star at that draw is <= lambda_star_used }
##   A^c = its complement (some group's block not PD, or lambda_star exceeds
##         the value that actually calibrated m_convergence for this run)
##
##   p_hat        = (# main-stage draws in A^c) / n_draws
##   tv_revised       = tv_tol + 2 * p_hat   (Section 17.4, triangle-inequality
##                                            decomposition + same-p_hat
##                                            heuristic for both A^c terms)
##   tv_revised_loose = tv_tol + p_hat       (one-sided sensitivity check)
##
## Deliberately built to REUSE the per-draw scan already computed by
## data-raw/_scratch_lambda_star_marginal_over_draws.R's
## .tmp_lambda_star_marginal_over_draws() (demo("Ex_13b_...")/("Ex_13c_...")
## Section 7d's `res_marg`) rather than recomputing it -- this function only
## adds the lambda_star_used comparison and the tv_tol arithmetic on top of
## that existing result.

## ---------------------------------------------------------------------------
## res_marg        -- output of .tmp_lambda_star_marginal_over_draws() (or the
##                     package's own .two_block_lambda_star_marginal_over_draws()):
##                     must carry $lambda_star_vec, $skipped, $n_draws.
## lambda_star_used -- the lambda_star that actually calibrated this run's
##                     m_convergence: fit$convergence_info$lambda_star_combined
##                     (the rank-matched-max combination of the marginal
##                     safeguard and the disp_upper plug-in envelope; see
##                     .two_block_combine_rate_envelopes()) if
##                     fit$convergence_info$marginal_rate_valid is TRUE,
##                     else fit$convergence_info$lambda_star_upper.
## tv_tol           -- the TV tolerance two_block_l_for_tv() was calibrated
##                      against (package default 0.01; see fit$convergence_info
##                      if the run recorded it, otherwise pass explicitly).
## ---------------------------------------------------------------------------
.tmp_tv_bound_revised <- function(res_marg, lambda_star_used, tv_tol = 0.01) {
  stopifnot(is.finite(lambda_star_used), is.finite(tv_tol))

  n_draws <- res_marg$n_draws
  lv      <- res_marg$lambda_star_vec
  skipped <- res_marg$skipped

  in_Ac <- skipped | (is.finite(lv) & lv > lambda_star_used)
  ## NA (never happens here: skipped draws have lv = NA and are already
  ## caught by `skipped`) would otherwise propagate through `>`.
  in_Ac[is.na(in_Ac)] <- TRUE

  n_exceed <- sum(in_Ac)
  p_hat    <- n_exceed / n_draws

  data.frame(
    lambda_star_used  = lambda_star_used,
    tv_tol            = tv_tol,
    n_draws           = n_draws,
    n_exceed          = n_exceed,
    p_hat             = p_hat,
    tv_revised        = tv_tol + 2 * p_hat,
    tv_revised_loose  = tv_tol + p_hat
  )
}

## ---------------------------------------------------------------------------
## Print a one-line summary in the same spirit as the ellipsoid test /
## marginal-scan cat() reports.
## ---------------------------------------------------------------------------
.tmp_print_tv_bound_revised <- function(res) {
  cat(sprintf(
    paste0(
      "\n=== Split-support revised TV bound (Section 17) ===\n\n",
      "  lambda_star_used (calibrated m_convergence) = %.6f\n",
      "  tv_tol (calibrated)                          = %.4f\n",
      "  draws in A^c (not PD, or lambda_star > used): %d / %d (p_hat = %.4f)\n",
      "  tv_revised        = tv_tol + 2*p_hat = %.4f\n",
      "  tv_revised_loose  = tv_tol +   p_hat = %.4f\n"
    ),
    res$lambda_star_used, res$tv_tol,
    res$n_exceed, res$n_draws, res$p_hat,
    res$tv_revised, res$tv_revised_loose
  ))
  invisible(res)
}

## ---------------------------------------------------------------------------
## Usage: right after Section 7d's res_marg <- .tmp_lambda_star_marginal_over_draws(...)
## in demo("Ex_13b_...")/demo("Ex_13c_..."):
##
##   source("data-raw/_scratch_tv_bound_revised.R")
##   lambda_star_used <- if (isTRUE(fit$convergence_info$marginal_rate_valid)) {
##     fit$convergence_info$lambda_star_combined  # rank-matched max w/ upper bound
##   } else {
##     fit$convergence_info$lambda_star_upper
##   }
##   tv_revised <- .tmp_tv_bound_revised(
##     res_marg, lambda_star_used, tv_tol = fit$convergence_info$tv_tol
##   )
##   .tmp_print_tv_bound_revised(tv_revised)
##
## No demo re-run needed -- this file only *defines* the functions above;
## sourcing it does not execute anything else.
## ---------------------------------------------------------------------------
