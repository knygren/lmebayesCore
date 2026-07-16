# -------------------------------------------------------------------------
#  Rcpp Interface Wrappers for glmbayes
#
#  These functions provide the minimal, strictly positional R → C++ bridges
#  required by the package.  Each wrapper mirrors the exact argument order
#  expected by the corresponding C++ routine and performs no preprocessing,
#  validation, or postprocessing.  Their sole purpose is to ensure that
#  high‑level R code calls the correct compiled symbol with the correct
#  signature.
#
#  All wrappers are internal:
#    - They are not part of the public API.
#    - They exist only to guarantee stable, explicit R–C++ boundaries.
#    - They prevent accidental reliance on .Call() with named arguments,
#      which R ignores, and which can silently break when signatures change.
#
#  Any future C++ interface changes must be reflected here to maintain
#  positional consistency and avoid NULL → double coercion errors.
#
#  Wrappers are organized by tier:
#    Tier 1: Core Simulation   - Main sampling entry points (rNormal_reg, etc.)
#    Tier 2: Envelope          - Envelope build/eval, EnvelopeCentering,
#                                rNormalGLM_std, rIndepNormalGammaReg_std
#    Tier 3: Model Utilities   - Standardization
#    Tier 4: OpenCL/GPU        - Kernel loading, GPU diagnostics
# -------------------------------------------------------------------------


# =============================================================================
#  Tier 1: Core Simulation
#  Callers: rNormal_reg, rNormalGamma_reg, rindepNormalGamma_reg, rGamma_reg,
#           block_rNormalGLM, block_rNormalReg
#  User:    All users – primary paths via rglmb, rlmb, glmb, pfamily
# =============================================================================

#' @noRd
#' @keywords internal
.rNormalGLM_cpp <- function(n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family = "binomial", link = "logit", Gridtype = 2L, n_envopt = -1L, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE) {
  .Call(`_lmebayesCore_rNormalGLM_cpp_export`, n, y, x, mu, P, offset, wt, dispersion, f2, f3, start, family, link, Gridtype, n_envopt, use_parallel, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.rNormalGLMBlocks_cpp <- function(n, y, x, offset, wt, dispersion, mu, P_blocks, prior_by_block, row_blocks, f2, f3, family = "binomial", link = "logit", Gridtype = 2L, n_envopt = -1L, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE) {
  .Call(`_lmebayesCore_rNormalGLMBlocks_cpp_export`, n, y, x, offset, wt, dispersion, mu, P_blocks, prior_by_block, row_blocks, f2, f3, family, link, Gridtype, n_envopt, use_parallel, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.rNormalRegBlocks_cpp <- function(
    n, y, x, offset, wt, dispersion, mu, P_blocks, prior_by_block, row_blocks,
    f2, f3,
    Gridtype = 2L
) {
  .Call(`_lmebayesCore_rNormalRegBlocks_cpp_export`,
    n, y, x, offset, wt, dispersion, mu, P_blocks, prior_by_block, row_blocks,
    f2, f3,
    Gridtype
  )
}

#' @noRd
#' @keywords internal
.block_rNormalReg_cpp <- function(
    n, y, x, block, prior_list, prior_lists, offset, wt, f2, f3, Gridtype = 2L
) {
  .Call(`_lmebayesCore_block_rNormalReg_cpp_export`,
    n, y, x, block, prior_list, prior_lists, offset, wt, f2, f3, Gridtype
  )
}

#' @noRd
#' @keywords internal
.block_rNormalGLM_cpp <- function(
    n, y, x, block, prior_list, prior_lists, offset, wt, f2, f3,
    family = "binomial", link = "logit",
    Gridtype = 2L, n_envopt = -1L,
    use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE
) {
  .Call(`_lmebayesCore_block_rNormalGLM_cpp_export`,
    n, y, x, block, prior_list, prior_lists, offset, wt, f2, f3,
    family, link, Gridtype, n_envopt, use_parallel, use_opencl, verbose
  )
}

#' @noRd
#' @keywords internal
.two_block_rNormal_reg_cpp <- function(
    n, m_convergence, y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    pfamily_list, fixef_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt,
    Gridtype = 2L, n_envopt = 1L,
    use_parallel = TRUE, use_opencl = FALSE,
    verbose = FALSE, progbar = TRUE
) {
  .Call(`_lmebayesCore_two_block_rNormal_reg_v2_cpp_export`,
    n, m_convergence, y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    pfamily_list, fixef_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt, Gridtype, n_envopt,
    use_parallel, use_opencl, verbose, progbar
  )
}

#' @noRd
#' @keywords internal
.two_block_rNormal_reg_staged_cpp <- function(
    n_main, m_convergence_main,
    n_pilot, m_convergence_pilot,
    y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    pfamily_list, fixef_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt,
    Gridtype = 2L, n_envopt = 1L,
    use_parallel = TRUE, use_opencl = FALSE,
    verbose = FALSE,
    progbar_main = TRUE, progbar_pilot = FALSE
) {
  .Call(`_lmebayesCore_two_block_rNormal_reg_staged_cpp_export`,
    n_main, m_convergence_main,
    n_pilot, m_convergence_pilot,
    y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    pfamily_list, fixef_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt, Gridtype, n_envopt,
    use_parallel, use_opencl, verbose, progbar_main, progbar_pilot
  )
}

# =============================================================================
#  Tier 1b: Two-block batch (Block 1 / Block 2 piecewise; rGLMM_sweep)
#  Callers: two_block_batch_gibbs.R, build_mu_all.R
# =============================================================================

#' @noRd
#' @keywords internal
.two_block_build_mu_all_cpp <- function(x_hyper, fixef, re_names, group_levels) {
  .Call(`_lmebayesCore_two_block_build_mu_all_cpp_export`,
    x_hyper, fixef, re_names, group_levels
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_prior_with_tau2_cpp <- function(
    base_prior, tau2_vec, ptypes, re_names, mu_all
) {
  .Call(`_lmebayesCore_two_block_block1_prior_with_tau2_cpp_export`,
    base_prior, tau2_vec, ptypes, re_names, mu_all
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_iters_mean_cpp <- function(block_out) {
  .Call(`_lmebayesCore_two_block_block1_iters_mean_cpp_export`, block_out)
}

#' @noRd
#' @keywords internal
.two_block_batch_fixef_chain_cpp <- function(batch_fixef, chain_i, re_names) {
  .Call(`_lmebayesCore_two_block_batch_fixef_chain_cpp_export`,
    batch_fixef, chain_i, re_names
  )
}

#' @noRd
#' @keywords internal
.two_block_batch_tau2_chain_row_cpp <- function(batch_tau2, chain_i) {
  .Call(`_lmebayesCore_two_block_batch_tau2_chain_row_cpp_export`,
    batch_tau2, chain_i
  )
}

#' @noRd
#' @keywords internal
.two_block_batch_b_assign_slice_cpp <- function(b_store, chain_i, b_draw) {
  .Call(`_lmebayesCore_two_block_batch_b_assign_slice_cpp_export`,
    b_store, chain_i, b_draw
  )
}

#' @noRd
#' @keywords internal
.two_block_batch_iters_ranef_add_cpp <- function(iters_ranef, chain_i, iters_mean) {
  .Call(`_lmebayesCore_two_block_batch_iters_ranef_add_cpp_export`,
    iters_ranef, chain_i, iters_mean
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_one_chain_draw_cpp <- function(
    chain_i, batch_fixef, tau2_i, y, Z, groups, offset, wt, x_hyper,
    re_names, group_levels, ptypes, block1_prior, is_gaussian,
    f2, f3, f2_gauss, f3_gauss, family, link, Gridtype, n_envopt
) {
  .Call(`_lmebayesCore_two_block_block1_one_chain_draw_cpp_export`,
    chain_i, batch_fixef, tau2_i, y, Z, groups, offset, wt, x_hyper,
    re_names, group_levels, ptypes, block1_prior, is_gaussian,
    f2, f3, f2_gauss, f3_gauss, family, link, Gridtype, n_envopt
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_one_chain_cpp <- function(
    chain_i, b_store, iters_ranef, batch_fixef, batch_tau2, design,
    block1_prior, family, ptypes, re_names, group_levels,
    f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add
) {
  .Call(`_lmebayesCore_two_block_block1_one_chain_cpp_export`,
    chain_i, b_store, iters_ranef, batch_fixef, batch_tau2, design,
    block1_prior, family, ptypes, re_names, group_levels,
    f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_one_chain_from_mu_P_cpp <- function(
    mu_all,
    P,
    dispersion,
    ddef,
    design,
    family,
    re_names,
    group_levels,
    f2,
    f3,
    f2_gauss,
    f3_gauss
) {
  .Call(`_lmebayesCore_two_block_block1_one_chain_from_mu_P_cpp_export`,
    mu_all, P, dispersion, ddef, design, family, re_names, group_levels,
    f2, f3, f2_gauss, f3_gauss
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_one_chain_v2_cpp <- function(
    fixef_i, tau2_i, design, block1_prior, family, ptypes,
    re_names, group_levels, f2, f3, f2_gauss, f3_gauss
) {
  .Call(`_lmebayesCore_two_block_block1_one_chain_v2_cpp_export`,
    fixef_i, tau2_i, design, block1_prior, family, ptypes,
    re_names, group_levels, f2, f3, f2_gauss, f3_gauss
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_all_chains_v2_internal_cpp <- function(
    fixef, chain_i, tau2, b, iters_ranef, design, block1_prior, family,
    ptypes, re_names, group_levels, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE
) {
  .Call(`_lmebayesCore_two_block_block1_all_chains_v2_internal_cpp_export`,
    fixef, chain_i, tau2, b, iters_ranef, design, block1_prior, family,
    ptypes, re_names, group_levels, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_all_chains_v2_internal_loop_cpp <- function(
    n, fixef, tau2, b_in_master, iters_ranef_in, design, block1_prior, family,
    ptypes, re_names, group_levels, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE,
    show_bar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE
) {
  .Call(`_lmebayesCore_two_block_block1_all_chains_v2_internal_loop_cpp_export`,
    n, fixef, tau2, b_in_master, iters_ranef_in, design, block1_prior, family,
    ptypes, re_names, group_levels, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add,
    show_bar, progbar_prefix, progbar_finish_newline
  )
}

#' @noRd
#' @keywords internal
.two_block_block1_all_chains_cpp <- function(
    n, fixef, tau2, b, iters_ranef, re_names, group_levels, design,
    block1_prior, family, ptypes, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add,
    show_bar, progbar_prefix, progbar_finish_newline
) {
  .Call(`_lmebayesCore_two_block_block1_all_chains_cpp_export`,
    n, fixef, tau2, b, iters_ranef, re_names, group_levels, design,
    block1_prior, family, ptypes, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add,
    show_bar, progbar_prefix, progbar_finish_newline
  )
}

#' @noRd
#' @keywords internal
.two_block_reorder_b_to_group_levels_cpp <- function(b_draw, block_ids, group_levels) {
  .Call(`_lmebayesCore_two_block_reorder_b_to_group_levels_cpp_export`,
    b_draw, block_ids, group_levels
  )
}

#' @noRd
#' @keywords internal
.two_block_align_b_to_xhyper_cpp <- function(b_vec, X_k, group_levels) {
  .Call(`_lmebayesCore_two_block_align_b_to_xhyper_cpp_export`,
    b_vec, X_k, group_levels
  )
}

#' @noRd
#' @keywords internal
.two_block_block2_one_chain_cpp <- function(
    b_i, fixef_rows, tau2_i, iters_i, x_hyper, group_levels,
    pfamily_list, ptypes, re_names
) {
  .Call(`_lmebayesCore_two_block_block2_one_chain_cpp_export`,
    b_i, fixef_rows, tau2_i, iters_i, x_hyper, group_levels,
    pfamily_list, ptypes, re_names
  )
}

#' @noRd
#' @keywords internal
.rNormalReg_cpp <- function(
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start,
    family = "gaussian",
    link = "identity",
    Gridtype = 2
) {
  .Call(`_lmebayesCore_rNormalReg_cpp_export`,
    n, y, x, mu, P, offset, wt, dispersion,
    f2, f3, start,
    family, link, Gridtype
  )
}

#' @noRd
#' @keywords internal
.rIndepNormalGammaReg_cpp <- function(n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower, disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose, progbar) {
  .Call(`_lmebayesCore_rIndepNormalGammaReg_cpp_export`, n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower, disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose, progbar)
}

#' @noRd
#' @keywords internal
.rIndepNormalGammaReg_with_envelope_cpp <- function(n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower, disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose, progbar) {
  .Call(`_lmebayesCore_rIndepNormalGammaReg_with_envelope_cpp_export`, n, y, x, mu, P, offset, wt, shape, rate, max_disp_perc, disp_lower, disp_upper, Gridtype, n_envopt, use_parallel, use_opencl, verbose, progbar)
}

#' @noRd
#' @keywords internal
.rNormalGammaReg_cpp <- function(n, y, x, mu, P, offset, wt, shape, rate,
                                 max_disp_perc, disp_lower, disp_upper,
                                 verbose = FALSE) {
  .Call(`_lmebayesCore_rNormalGammaReg_cpp_export`,
        n, y, x, mu, P, offset, wt, shape, rate,
        max_disp_perc, disp_lower, disp_upper, verbose)
}

#' @noRd
#' @keywords internal
.rGammaGaussian_cpp <- function(n, y, x, beta, wt, alpha, shape, rate,
                                disp_lower = NULL, disp_upper = NULL,
                                verbose = FALSE) {
  .Call(`_lmebayesCore_rGammaGaussian_cpp_export`,
        n, y, x, beta, wt, alpha, shape, rate,
        disp_lower, disp_upper, verbose)
}

#' @noRd
#' @keywords internal
.rGammaGamma_cpp <- function(n, y, x, beta, wt, alpha, shape, rate,
                             max_disp_perc, disp_lower = NULL,
                             disp_upper = NULL, verbose = FALSE) {
  .Call(`_lmebayesCore_rGammaGamma_cpp_export`,
        n, y, x, beta, wt, alpha, shape, rate,
        max_disp_perc, disp_lower, disp_upper, verbose)
}


# =============================================================================
#  Tier 2: Envelope & Standardization
#  Callers: EnvelopeSize, EnvelopeBuild, EnvelopeEval, EnvelopeDispersionBuild,
#           EnvelopeOrchestrator, EnvelopeCentering, rNormalGLM_std,
#           rIndepNormalGammaReg_std; EnvelopeSet_* are internal
#  User:    Advanced users – understanding algorithm, custom envelope workflows
# =============================================================================

#' @noRd
#' @keywords internal
.rNormalGLM_std_cpp <- function(n, y, x, mu, P, alpha, wt,
                                f2, Envelope,
                                family, link,
                                progbar = 1L,
                                verbose = FALSE) {
  .Call(`_lmebayesCore_rNormalGLM_std_cpp_export`,
        n, y, x, mu, P, alpha, wt,
        f2, Envelope,
        family, link,
        progbar, verbose)
}

#' @noRd
#' @keywords internal
.rIndepNormalGammaReg_std_cpp <- function(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose) {
  .Call(`_lmebayesCore_rIndepNormalGammaReg_std_cpp_export`, n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose)
}

#' @noRd
#' @keywords internal
.rIndepNormalGammaReg_std_parallel_cpp <- function(n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose) {
  .Call(`_lmebayesCore_rIndepNormalGammaReg_std_parallel_cpp_export`, n, y, x, mu, P, alpha, wt, f2, Envelope, gamma_list, UB_list, family, link, progbar, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeCentering_cpp <- function(y, x, mu, P, offset, wt, shape, rate, Gridtype = 2L, verbose = FALSE) {
  .Call(`_lmebayesCore_EnvelopeCentering_cpp_export`, y, x, mu, P, offset, wt, shape, rate, Gridtype, verbose)
}

#' @noRd
#' @keywords internal
.BlockEnvelopeCentering_cpp <- function(
    y, x, block, prior_list, prior_lists,
    offset, wt, shape, rate, max_disp_perc,
    disp_lower = NULL, disp_upper = NULL,
    p_re = -1L, n_rss_iter = 10L, verbose = FALSE
) {
  .Call(`_lmebayesCore_BlockEnvelopeCentering_cpp_export`,
    y, x, block, prior_list, prior_lists,
    offset, wt, shape, rate, max_disp_perc,
    disp_lower, disp_upper, p_re, n_rss_iter, verbose)
}

#' @noRd
#' @keywords internal
.BlockEnvelopeBuild_cpp <- function(
    centering_out, y, x, block, prior_list, prior_lists,
    offset, wt, max_disp_perc,
    disp_lower = NULL, disp_upper = NULL,
    n = 1L, Gridtype = 3L, n_envopt = -1L,
    RSS_ML = NA_real_,
    use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE
) {
  .Call(`_lmebayesCore_BlockEnvelopeBuild_cpp_export`,
    centering_out, y, x, block, prior_list, prior_lists,
    offset, wt, max_disp_perc, disp_lower, disp_upper,
    n, Gridtype, n_envopt, RSS_ML,
    use_parallel, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.BlockEnvelopeDispersionBuild_cpp <- function(
    build_out,
    centering_out,
    y,
    x,
    block,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower = NULL,
    disp_upper = NULL,
    RSS_ML = NA_real_,
    use_parallel = TRUE,
    verbose = FALSE
) {
  .Call(`_lmebayesCore_BlockEnvelopeDispersionBuild_cpp_export`,
    build_out,
    centering_out,
    y,
    x,
    block,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    RSS_ML,
    use_parallel,
    verbose
  )
}

#' @noRd
#' @keywords internal
.BlockEnvelopeSim_cpp <- function(
    build_out,
    n = 1L,
    progbar = FALSE,
    verbose = FALSE
) {
  .Call(`_lmebayesCore_BlockEnvelopeSim_cpp_export`,
    build_out, n, progbar, verbose)
}

#' Block ING envelope pipeline: Centering → Build → DispersionBuild → Sim
#' @noRd
#' @keywords internal
.rIndepNormalGammaRegBlock_cpp <- function(
    n,
    y,
    x,
    block,
    prior_list,
    prior_lists = NULL,
    offset,
    wt,
    p_re = -1L,
    n_rss_iter = 10L,
    Gridtype = 3L,
    n_envopt = -1L,
    RSS_ML = NA_real_,
    use_parallel = TRUE,
    use_opencl = FALSE,
    progbar = FALSE,
    verbose = FALSE,
    group_levels = character(0),
    re_names = character(0)
) {
  .Call(
    `_lmebayesCore_rIndepNormalGammaRegBlock_cpp_export`,
    n, y, x, block, prior_list, prior_lists,
    offset, wt, p_re, n_rss_iter, Gridtype, n_envopt, RSS_ML,
    use_parallel, use_opencl, progbar, verbose,
    group_levels, re_names
  )
}

#' @noRd
#' @keywords internal
.EnvelopeSize_cpp <- function(a, G1, Gridtype, n, n_envopt, use_opencl, verbose) {
  .Call(`_lmebayesCore_EnvelopeSize_cpp_export`, a, G1, Gridtype, n, n_envopt, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeBuild_cpp <- function(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose) {
  .Call(`_lmebayesCore_EnvelopeBuild_cpp_export`, bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeBuild_Ind_Normal_Gamma_cpp <- function(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose) {
  .Call(`_lmebayesCore_EnvelopeBuild_Ind_Normal_Gamma_cpp_export`, bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeEval_cpp <- function(G4, y, x, mu, P, alpha, wt,
                          family, link,
                          use_opencl = FALSE,
                          verbose = FALSE) {
  .Call(`_lmebayesCore_EnvelopeEval_cpp_export`,
        G4, y, x, mu, P, alpha, wt,
        family, link,
        use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeDispersionBuild_cpp <- function(
    Env,
    Shape,
    Rate,
    P,
    y,
    x,
    alpha,
    n_obs,
    RSS_post,
    RSS_ML,
    mu,
    wt,
    max_disp_perc,
    disp_lower = NULL,
    disp_upper = NULL,
    verbose = FALSE,
    use_parallel = TRUE
) {
  .Call(`_lmebayesCore_EnvelopeDispersionBuild_cpp_export`,
    Env,
    Shape,
    Rate,
    P,
    y,
    x,
    alpha,
    n_obs,
    RSS_post,
    RSS_ML,
    mu,
    wt,
    max_disp_perc,
    disp_lower,
    disp_upper,
    verbose,
    use_parallel
  )
}

#' @noRd
#' @keywords internal
.EnvelopeOrchestrator_cpp <- function(bstar2, A, y, x2, mu2, P2, alpha, wt, n, Gridtype, n_envopt, shape, rate, RSS_Post2, RSS_ML, max_disp_perc, disp_lower, disp_upper, use_parallel, use_opencl, verbose) {
  .Call(`_lmebayesCore_EnvelopeOrchestrator_cpp_export`, bstar2, A, y, x2, mu2, P2, alpha, wt, n, Gridtype, n_envopt, shape, rate, RSS_Post2, RSS_ML, max_disp_perc, disp_lower, disp_upper, use_parallel, use_opencl, verbose)
}

#' @noRd
#' @keywords internal
.EnvelopeSet_Grid_cpp <- function(GIndex, cbars, Lint) {
  .Call(`_lmebayesCore_EnvelopeSet_Grid_cpp_export`, GIndex, cbars, Lint)
}

#' @noRd
#' @keywords internal
.EnvelopeSet_LogP_cpp <- function(logP, NegLL, cbars, G3) {
  .Call(`_lmebayesCore_EnvelopeSet_LogP_cpp_export`, logP, NegLL, cbars, G3)
}


# =============================================================================
#  Tier 3: Model Utilities
#  Callers: glmb_Standardize_Model
#  User:    Advanced users – model preparation, standardization
# =============================================================================

#' @noRd
#' @keywords internal
.glmb_Standardize_Model_cpp <- function(y, x, P, bstar, A1) {
  .Call(`_lmebayesCore_glmb_Standardize_Model_cpp_export`, y, x, P, bstar, A1)
}


# =============================================================================
#  Tier 4: OpenCL / GPU
#  Callers: glmbayesCore_has_opencl, gpu_names
#  Kernel loading / core count: opencltools (see ?opencltools::load_kernel_source)
# =============================================================================

#' @noRd
#' @keywords internal
.glmbayesCore_has_opencl_cpp <- function() {
  .Call(`_lmebayesCore_glmbayesCore_has_opencl_cpp_export`)
}

#' @noRd
#' @keywords internal
.gpu_names_cpp <- function() {
  .Call(`_lmebayesCore_gpu_names_cpp_export`)
}


# =============================================================================
#  Independent-block separable-overbound variants (Ind)
#  See Appendix A of inst/BLOCK_ING_RINDEPNORMALGAMMA_REG.md for the theory.
# =============================================================================

#' @noRd
#' @keywords internal
.BlockEnvelopeDispersionBuildInd_cpp <- function(
    build_out,
    centering_out,
    y,
    x,
    block,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower = NULL,
    disp_upper = NULL,
    RSS_ML = NA_real_,
    use_parallel = TRUE,
    verbose = FALSE
) {
  .Call(`_lmebayesCore_BlockEnvelopeDispersionBuildInd_cpp_export`,
    build_out,
    centering_out,
    y,
    x,
    block,
    offset,
    wt,
    shape,
    rate,
    max_disp_perc,
    disp_lower,
    disp_upper,
    RSS_ML,
    use_parallel,
    verbose
  )
}

#' @noRd
#' @keywords internal
.BlockEnvelopeSimInd_cpp <- function(
    build_out,
    n = 1L,
    progbar = FALSE,
    verbose = FALSE
) {
  .Call(`_lmebayesCore_BlockEnvelopeSimInd_cpp_export`,
    build_out,
    n,
    progbar,
    verbose
  )
}

#' Block ING independent-block pipeline: Centering → Build → DispersionBuildInd → SimInd
#' @noRd
#' @keywords internal
.rIndepNormalGammaRegBlockInd_cpp <- function(
    n,
    y,
    x,
    block,
    prior_list,
    prior_lists = NULL,
    offset,
    wt,
    p_re = -1L,
    n_rss_iter = 10L,
    Gridtype = 3L,
    n_envopt = -1L,
    RSS_ML = NA_real_,
    use_parallel = TRUE,
    use_opencl = FALSE,
    progbar = FALSE,
    verbose = FALSE,
    group_levels = character(0),
    re_names = character(0)
) {
  .Call(
    `_lmebayesCore_rIndepNormalGammaRegBlockInd_cpp_export`,
    n, y, x, block, prior_list, prior_lists,
    offset, wt, p_re, n_rss_iter, Gridtype, n_envopt, RSS_ML,
    use_parallel, use_opencl, progbar, verbose,
    group_levels, re_names
  )
}

# =============================================================================
#  Phased Out (no R wrappers; C++ exports may still exist for compatibility)
#  - .rss_face_at_disp_cpp, .UB2_cpp
#  - Former RSS/UB2 minimization callbacks; active path uses closed-form C++ bounds
# =============================================================================
