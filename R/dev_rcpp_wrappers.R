# -------------------------------------------------------------------------
#  DEV ONLY -- does not migrate to glmbayesCore.
#  Rcpp wrappers for src/dev_exports.cpp (two-block replay / record /
#  prior tracking). Strictly positional, no validation.
# -------------------------------------------------------------------------

#' @noRd
#' @keywords internal
.dev_two_block_log_new_cpp <- function(spec, n, M) {
  .Call(`_lmebayesCore_dev_two_block_log_new_cpp_export`, spec, n, M)
}

#' @noRd
#' @keywords internal
.dev_two_block_sweep_cpp <- function(
    spec, b, fixef, iters, iters_ranef,
    inner_sweeps, use_mode = FALSE, progbar = FALSE, stage_label = "",
    tape = NULL, record = FALSE
) {
  .Call(`_lmebayesCore_dev_two_block_sweep_cpp_export`,
    spec, b, fixef, iters, iters_ranef,
    inner_sweeps, use_mode, progbar, stage_label,
    tape, record
  )
}
