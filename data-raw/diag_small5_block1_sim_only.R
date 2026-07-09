## Time Block-1 draw with realistic mu_all; 15s timeout.
suppressPackageStartupMessages({
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("need pkgload")
  root <- normalizePath("../lmebayes")
  core <- normalizePath(".")
  pkgload::load_all(core, quiet = TRUE)
  pkgload::load_all(root, quiet = TRUE)
})
source(file.path(root, "tests/manual/_small5_lmerb_fixture.R"))

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
ps <- Prior_Setup_lmebayes(fx$form, data = fx$dat, pwt = 0.01)
pf <- pfamily_list(ps)
Z <- as.matrix(fx$design$Z)
y <- fx$design$y
groups <- fx$design$groups
re_names <- fx$design$re_coef_names
group_levels <- levels(groups)
p_re <- length(re_names)

m_disp <- ps$ing_prior_measurement
disp_lower <- m_disp$disp_lower + 0.25 * (m_disp$disp_upper - m_disp$disp_lower)
disp_upper <- m_disp$disp_upper - 0.25 * (m_disp$disp_upper - m_disp$disp_lower)

prior_container <- glmbayesCore:::.lmebayes_priors_from_pfamily_list(
  pfamily_list = pf,
  dispersion_ranef = dGamma(
    shape = m_disp$shape, rate = m_disp$rate,
    beta = matrix(0, 1, 1), Inv_Dispersion = TRUE,
    disp_lower = disp_lower, disp_upper = disp_upper
  ),
  design = fx$design, family = gaussian(), fn_name = "diag"
)
P <- solve(prior_container$Sigma_ranef)
block1 <- glmbayesCore:::.lmebayes_block1_prior_list(
  prior_container, dispersion_ranef = m_disp$rate / (m_disp$shape - 1)
)
icm <- glmbayesCore:::.two_block_icm_at_start(
  design = fx$design, prior_list = block1,
  pfamily_list = pf, re_names = re_names, family = gaussian()
)
mu_all <- as.matrix(glmbayesCore::build_mu_all(fx$design, icm$start)$mu_all)

prior_list <- list(
  mu = mu_all,
  Sigma = prior_container$Sigma_ranef,
  shape = m_disp$shape, rate = m_disp$rate,
  max_disp_perc = 0.99,
  disp_lower = disp_lower, disp_upper = disp_upper
)

run_one <- function(label, fn) {
  cat("\n", label, " (15s timeout)...\n", sep = "")
  flush.console()
  setTimeLimit(cpu = 15, elapsed = 15, transient = TRUE)
  res <- tryCatch({
    t0 <- proc.time()
    out <- fn()
    list(ok = TRUE, out = out, elapsed = (proc.time() - t0)["elapsed"])
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE)
  if (isTRUE(res$ok)) {
    cat(sprintf("  OK: iters_mean=%.0f elapsed=%.3fs\n",
                res$out$iters_mean, res$elapsed))
  } else {
    cat("  FAIL:", res$msg, "\n")
  }
  invisible(res)
}

args <- list(
  n = 1L, y = y, x = Z, block = groups, prior_list = prior_list,
  prior_lists = NULL, offset = rep(0, length(y)), wt = rep(1, length(y)),
  p_re = p_re, n_rss_iter = 10L, Gridtype = 3L, n_envopt = -1L,
  RSS_ML = NA_real_, use_parallel = FALSE, use_opencl = FALSE,
  progbar = FALSE, verbose = FALSE,
  group_levels = group_levels, re_names = re_names
)

run_one("IND", function() do.call(glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp, args))
run_one("JOINT", function() do.call(glmbayesCore:::.rIndepNormalGammaRegBlock_cpp, args))
