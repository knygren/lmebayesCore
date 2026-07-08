# Acceptance rate comparison only, n = 10000
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

set.seed(seed)
sim_legacy <- rindepNormalGamma_reg(
  n = n_draws, y = y, x = x_old, prior_list = prior_list_old,
  n_envopt = n_envopt, Gridtype = 3L, use_parallel = FALSE, progbar = FALSE
)
set.seed(seed)
sim_block <- glmbayesCore:::.rIndepNormalGammaRegBlock_cpp(
  n = n_draws, y = y, x = x_block, block = block,
  prior_list = prior_list_block, n_envopt = n_envopt, Gridtype = 3L,
  use_parallel = FALSE, progbar = FALSE, verbose = FALSE,
  offset = rep(0, length(y)), wt = rep(1, length(y)),
  p_re = -1L, n_rss_iter = 10L, RSS_ML = NA_real_, use_opencl = FALSE,
  group_levels = character(0), re_names = character(0)
)

legacy_iters <- sim_legacy$iters
block_iters <- sim_block$iters_out
legacy_dpa <- mean(legacy_iters)
block_dpa <- mean(block_iters)
legacy_rate <- 1 / legacy_dpa
block_rate <- 1 / block_dpa
abs_diff <- block_dpa - legacy_dpa
rel_diff <- abs_diff / legacy_dpa

cat(sprintf("n = %d, seed = %d\n\n", n_draws, seed))
cat(sprintf("Legacy  accept rate     = %.5f  (%.2f%%)\n", legacy_rate, 100 * legacy_rate))
cat(sprintf("Legacy  mean(iters)      = %.3f  (median iters = %g)\n", legacy_dpa, median(legacy_iters)))
cat(sprintf("Block   accept rate     = %.5f  (%.2f%%)\n", block_rate, 100 * block_rate))
cat(sprintf("Block   mean(iters)     = %.3f  (median iters = %g)\n", block_dpa, median(block_iters)))
cat(sprintf("\nAbsolute gap mean(iters) (block - legacy) = %+.3f\n", abs_diff))
cat(sprintf("Block / legacy mean(iters) ratio              = %.3f\n", block_dpa / legacy_dpa))
cat(sprintf("Accept-rate ratio (legacy/block)              = %.3f\n", legacy_rate / block_rate))
