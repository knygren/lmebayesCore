# Two-block stacked Dobson: legacy vs block, n = 10000 draws
pkgload::load_all(quiet = TRUE)

n_draws <- 10000L
n_envopt <- 10000L
seed <- 360L

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group1 <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight1 <- c(ctl, trt)
n1 <- length(weight1)
weight <- c(weight1, weight1)
group_stacked <- factor(c(group1, group1), levels = levels(group1))
p_setup <- Prior_Setup(
  weight ~ group,
  data = data.frame(weight = weight, group = group_stacked),
  family = gaussian()
)
y <- p_setup$y
x_block <- p_setup$x
l1 <- ncol(x_block)
block <- factor(c(rep("B1", n1), rep("B2", n1)), levels = c("B1", "B2"))
x1 <- x_block[seq_len(n1), , drop = FALSE]
zeros <- matrix(0, n1, l1)
x_old <- rbind(cbind(x1, zeros), cbind(zeros, x1))
mu1 <- p_setup$mu
Sigma1 <- p_setup$Sigma
prior_list_old <- list(
  mu = c(mu1, mu1),
  Sigma = rbind(
    cbind(Sigma1, matrix(0, l1, l1)),
    cbind(matrix(0, l1, l1), Sigma1)
  ),
  dispersion = p_setup$dispersion,
  shape = p_setup$shape,
  rate = p_setup$rate,
  Precision = solve(rbind(
    cbind(Sigma1, matrix(0, l1, l1)),
    cbind(matrix(0, l1, l1), Sigma1)
  )),
  max_disp_perc = 0.99
)
prior_list_block <- list(
  mu = mu1,
  Sigma = Sigma1,
  shape = p_setup$shape,
  rate = p_setup$rate,
  max_disp_perc = 0.99
)

cat("=== Envelope / gamma / UB constants (n = 1 build) ===\n")
sim_env <- rindepNormalGamma_reg_with_envelope(
  n = 1L,
  y = y,
  x = x_old,
  prior_list = prior_list_old,
  n_envopt = n_envopt,
  Gridtype = 3L,
  use_parallel = FALSE,
  progbar = FALSE
)
sim_block_build <- glmbayesCore:::.rIndepNormalGammaRegBlock_cpp(
  n = 1L,
  y = y,
  x = x_block,
  block = block,
  prior_list = prior_list_block,
  n_envopt = n_envopt,
  Gridtype = 3L,
  use_parallel = FALSE,
  progbar = FALSE,
  verbose = FALSE,
  offset = rep(0, length(y)),
  wt = rep(1, length(y)),
  p_re = -1L,
  n_rss_iter = 10L,
  RSS_ML = NA_real_,
  use_opencl = FALSE,
  group_levels = character(0),
  re_names = character(0)
)
de <- sim_block_build$build_out$dispersion_envelope
ubL <- sim_env$UB_list
gamL <- sim_env$gamma_list
ubB <- de$UB_list
gamB <- de$gamma_list
keys <- c(
  "lmc1", "lmc2", "lm_log1", "lm_log2", "max_New_LL_UB",
  "max_LL_log_disp", "RSS_Min", "shape3", "rate2", "disp_lower", "disp_upper"
)
for (k in keys) {
  vL <- if (k %in% names(ubL)) ubL[[k]] else gamL[[k]]
  vB <- if (k %in% names(ubB)) {
    ubB[[k]]
  } else if (k == "RSS_Min") {
    de$RSS_Min
  } else {
    gamB[[k]]
  }
  cat(sprintf(
    "%-16s legacy %12g  block %12g  rel_diff %8.3e\n",
    k, vL, vB, (vB - vL) / max(1, abs(vL))
  ))
}
cat("aggregation:", de$cross_face_meta$aggregation, "\n\n")

cat("=== Sampling comparison (n =", n_draws, ", seed =", seed, ") ===\n")
t0 <- Sys.time()
set.seed(seed)
sim_legacy <- rindepNormalGamma_reg(
  n = n_draws,
  y = y,
  x = x_old,
  prior_list = prior_list_old,
  n_envopt = n_envopt,
  Gridtype = 3L,
  use_parallel = FALSE,
  progbar = FALSE
)
t_legacy <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t0 <- Sys.time()
set.seed(seed)
sim_block <- glmbayesCore:::.rIndepNormalGammaRegBlock_cpp(
  n = n_draws,
  y = y,
  x = x_block,
  block = block,
  prior_list = prior_list_block,
  n_envopt = n_envopt,
  Gridtype = 3L,
  use_parallel = FALSE,
  progbar = FALSE,
  verbose = FALSE,
  offset = rep(0, length(y)),
  wt = rep(1, length(y)),
  p_re = -1L,
  n_rss_iter = 10L,
  RSS_ML = NA_real_,
  use_opencl = FALSE,
  group_levels = character(0),
  re_names = character(0)
)
t_block <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

legacy_iters <- sim_legacy$iters
block_iters <- sim_block$iters_out
legacy_dpa <- mean(legacy_iters)
block_dpa <- mean(block_iters)
legacy_rate <- 1 / legacy_dpa
block_rate <- 1 / block_dpa

cat(sprintf(
  "Legacy: mean(iters) = %.4f  accept rate = %.4f  wall_sec = %.1f\n",
  legacy_dpa, legacy_rate, t_legacy
))
cat(sprintf(
  "Block:  mean(iters_out) = %.4f  accept rate = %.4f  wall_sec = %.1f\n",
  block_dpa, block_rate, t_block
))
cat(sprintf("Ratio block/legacy mean(iters): %.3f\n", block_dpa / legacy_dpa))
cat(sprintf(
  "Block iters_out: min=%s max=%s  prop(1)=%.4f\n",
  min(block_iters), max(block_iters), mean(block_iters == 1)
))
cat(sprintf(
  "Legacy iters: min=%s max=%s  median=%s\n",
  min(legacy_iters), max(legacy_iters), median(legacy_iters)
))

cat("\nDispersion summaries:\n")
cat(sprintf(
  "Legacy mean(dispersion) = %.6f  sd = %.6f\n",
  mean(sim_legacy$dispersion), sd(sim_legacy$dispersion)
))
cat(sprintf(
  "Block  mean(dispersion) = %.6f  sd = %.6f\n",
  mean(sim_block$dispersion), sd(sim_block$dispersion)
))
