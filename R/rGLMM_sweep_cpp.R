#' Portable two-block sweep, C++ loop
#'
#' Same arguments and output as \code{.rGLMM_sweep_indexed}. Setup and the
#' final naming are the shared R steps; the sweep loop is
#' \code{two_block_port_sweep_cpp_export}, the C++ port of
#' \code{.two_block_port_sweep_r}. \code{rwt} is accepted and is not used
#' inside the draw.
#' @param index_map Result of \code{.two_block_index_map()}.
#' @noRd
.rGLMM_sweep_cpp <- function(
    n_chains,
    start_fixef,
    inner_sweeps,
    design,
    block1_prior,
    pfamily_list,
    family,
    re_names,
    group_levels,
    collect_block1 = TRUE,
    progbar        = FALSE,
    stage_label    = "",
    diag_sweeps    = FALSE,
    fixef_mode     = NULL,
    b_mode         = NULL,
    b_start        = NULL,
    ptypes         = NULL,
    tau2_start     = NULL,
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE,
    use_cpp_block2 = TRUE,
    rwt = 0.01,
    index_map
) {
  force(rwt)
  force(diag_sweeps)
  force(fixef_mode)
  force(ptypes)
  force(use_cpp_tau2_row)
  force(use_cpp_b_slice)
  force(use_cpp_iters_ranef_add)
  force(use_cpp_block2)

  prep <- .two_block_port_prepare(
    n_chains, start_fixef, design, block1_prior, pfamily_list, family,
    re_names, group_levels, b_start, b_mode, tau2_start, index_map
  )
  state <- two_block_port_sweep_cpp_export(
    spec         = prep$spec,
    b            = prep$state$b,
    fixef        = prep$state$fixef,
    iters        = prep$state$iters,
    iters_ranef  = prep$state$iters_ranef,
    inner_sweeps = as.integer(inner_sweeps),
    use_mode     = FALSE,
    progbar      = isTRUE(progbar),
    stage_label  = if (is.null(stage_label)) "" else as.character(stage_label[1L])
  )
  .two_block_port_finish(
    state, prep$tau2, index_map, design, re_names, collect_block1, prep$spec
  )
}
