## Diagnose where Block-1 ING hangs on the 5-school fixture.
## Run from lmebayes root:
##   Rscript ../glmbayesCore/data-raw/diag_small5_block1_one_chain.R

suppressPackageStartupMessages({
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("need pkgload")
  root <- Sys.getenv("LMEBAYES_ROOT", unset = normalizePath("../lmebayes"))
  core <- Sys.getenv("GLMBAYESCORE_ROOT", unset = normalizePath("."))
  pkgload::load_all(core, quiet = TRUE)
  pkgload::load_all(root, quiet = TRUE)
})

source(file.path(root, "tests/manual/_small5_lmerb_fixture.R"))

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
ps <- Prior_Setup_lmebayes(fx$form, data = fx$dat, pwt = 0.01)
Z <- as.matrix(fx$design$Z)
y <- fx$design$y
groups <- fx$design$groups
re_names <- fx$design$re_coef_names
group_levels <- levels(groups)
p_re <- length(re_names)

ing <- ps$ing_prior_measurement
prior_list <- list(
  mu = matrix(0, nrow = p_re, ncol = length(group_levels),
              dimnames = list(re_names, group_levels)),
  Sigma = diag(2),
  shape = ing$shape,
  rate = ing$rate,
  max_disp_perc = 0.99,
  disp_lower = ing$disp_lower,
  disp_upper = ing$disp_upper
)
rownames(prior_list$Sigma) <- re_names
colnames(prior_list$Sigma) <- re_names

cat("\n=== Direct .rIndepNormalGammaRegBlockInd_cpp (n=1, verbose=TRUE) ===\n")
cat("use_parallel=FALSE\n\n")
flush.console()

t0 <- proc.time()
out <- glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp(
  n             = 1L,
  y             = y,
  x             = Z,
  block         = groups,
  prior_list    = prior_list,
  prior_lists   = NULL,
  offset        = rep(0, length(y)),
  wt            = rep(1, length(y)),
  p_re          = p_re,
  n_rss_iter    = 10L,
  Gridtype      = 3L,
  n_envopt      = -1L,
  RSS_ML        = NA_real_,
  use_parallel  = FALSE,
  use_opencl    = FALSE,
  progbar       = FALSE,
  verbose       = TRUE,
  group_levels  = group_levels,
  re_names      = re_names
)
t1 <- proc.time()

cat(sprintf("\nDone: elapsed = %.2f s, iters_mean = %.0f, dispersion = %.6g\n",
            (t1 - t0)["elapsed"], out$iters_mean, out$dispersion_ranef))
