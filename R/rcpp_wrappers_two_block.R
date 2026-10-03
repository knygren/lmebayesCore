# -------------------------------------------------------------------------
#  Rcpp interface wrappers for the two-block Gibbs engines and the
#  block ING (rIndepNormalGammaRegBlock*) envelope pipeline.
#  Same conventions as rcpp_wrappers.R: strictly positional, no validation.
# -------------------------------------------------------------------------

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
