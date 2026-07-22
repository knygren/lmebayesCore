## Stress test: repeatedly call the Block-1 Ind envelope draw at the mode
## start (pwt_measurement = 0.2, worst-case AR acceptance) to try to trigger
## the non-finite-proposal diagnostic many times quickly.
suppressPackageStartupMessages({
  root <- Sys.getenv("LMEBAYES_ROOT", unset = normalizePath("../lmebayes"))
  core <- Sys.getenv("GLMBAYESCORE_ROOT", unset = normalizePath("."))
  pkgload::load_all(core, quiet = TRUE)
  pkgload::load_all(root, quiet = TRUE)
})
source(file.path(root, "tests/manual/_small5_lmerb_fixture.R"))

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
ps <- Prior_Setup_lmebayes(fx$form, data = fx$dat, pwt = 0.01, pwt_measurement = 0.2)
pf <- pfamily_list(ps)
Z <- as.matrix(fx$design$Z)
y <- fx$design$y
groups <- fx$design$groups
re_names <- fx$design$re_coef_names
group_levels <- levels(groups)
p_re <- length(re_names)

m_disp <- ps$ing_prior_measurement
disp_pf <- dGamma(
  shape = m_disp$shape, rate = m_disp$rate,
  beta = matrix(0, 1, 1), Inv_Dispersion = TRUE,
  disp_lower = m_disp$disp_lower, disp_upper = m_disp$disp_upper
)
prior_container <- lmebayesCore::priors_from_pfamily_list(
  pfamily_list = pf, dispersion_ranef = disp_pf,
  design = fx$design, family = gaussian(), fn_name = "diag"
)
block1 <- glmbayesCore:::.lmebayes_block1_prior_list(
  prior_container, dispersion_ranef = m_disp$rate / (m_disp$shape - 1)
)
icm <- glmbayesCore:::.two_block_icm_at_start(
  design = fx$design, prior_list = block1,
  pfamily_list = pf, re_names = re_names, family = gaussian()
)
mu_all <- as.matrix(glmbayesCore::build_mu_all(fx$design, icm$start)$mu_all)

prior_list <- list(
  mu = mu_all, Sigma = prior_container$Sigma_ranef,
  shape = m_disp$shape, rate = m_disp$rate, max_disp_perc = 0.99,
  disp_lower = m_disp$disp_lower, disp_upper = m_disp$disp_upper
)

n_reps <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(n_reps)) n_reps <- 15L
cat("Running", n_reps, "reps...\n")

n_ok <- 0L
n_err <- 0L
for (rep in seq_len(n_reps)) {
  t0 <- proc.time()
  res <- tryCatch({
    out <- glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp(
      n = 1L, y = y, x = Z, block = groups, prior_list = prior_list,
      prior_lists = NULL, offset = rep(0, length(y)), wt = rep(1, length(y)),
      p_re = p_re, n_rss_iter = 10L, Gridtype = 3L, n_envopt = -1L,
      RSS_ML = NA_real_, use_parallel = TRUE, use_opencl = FALSE,
      progbar = FALSE, verbose = FALSE,
      group_levels = group_levels, re_names = re_names
    )
    list(ok = TRUE, iters = out$iters_mean)
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  el <- (proc.time() - t0)["elapsed"]
  if (isTRUE(res$ok)) {
    n_ok <- n_ok + 1L
    cat(sprintf("[%2d] OK   iters=%.0f elapsed=%.1fs\n", rep, res$iters, el))
  } else {
    n_err <- n_err + 1L
    cat(sprintf("[%2d] FAIL elapsed=%.1fs msg=%s\n", rep, el, res$msg))
  }
}
cat(sprintf("\nDone: %d ok, %d fail (of %d)\n", n_ok, n_err, n_reps))
