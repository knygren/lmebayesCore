## Isolated repro: replay chain 2's captured Block-1 ING call directly,
## bypassing the full pilot loop, with verbose=TRUE for fine-grained tracing.
suppressPackageStartupMessages({
  root <- Sys.getenv("LMEBAYES_ROOT", unset = normalizePath("../lmebayes"))
  core <- Sys.getenv("GLMBAYESCORE_ROOT", unset = normalizePath("."))
  pkgload::load_all(core, quiet = TRUE)
})

cap <- readRDS(file.path(root, "chain_capture_02.rds"))
cat("chain_i =", cap$chain_i, "\n")
str(cap$prior_list)

cat("\n=== Direct .rIndepNormalGammaRegBlockInd_cpp (verbose=TRUE, use_parallel=FALSE) ===\n\n")
out <- tryCatch({
  glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp(
    n             = 1L,
    y             = cap$y,
    x             = cap$Z,
    block         = cap$groups,
    prior_list    = cap$prior_list,
    prior_lists   = NULL,
    offset        = rep(0, length(cap$y)),
    wt            = rep(1, length(cap$y)),
    p_re          = cap$p_re,
    n_rss_iter    = 10L,
    Gridtype      = 3L,
    n_envopt      = -1L,
    RSS_ML        = NA_real_,
    use_parallel  = FALSE,
    use_opencl    = FALSE,
    progbar       = FALSE,
    verbose       = TRUE,
    group_levels  = cap$group_levels,
    re_names      = cap$re_names
  )
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(out)) {
  cat(sprintf("\nOK: dispersion_ranef=%.4f iters_mean=%.0f\n", out$dispersion_ranef, out$iters_mean))
}
