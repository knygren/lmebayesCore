# Mean / SD coefficients and dispersion: two-block stacked Dobson, n = 10000
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
colnames(x_old) <- c(colnames(x_block), paste0(colnames(x_block), ".B2"))
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
  n = n_draws,
  y = y,
  x = x_old,
  prior_list = prior_list_old,
  n_envopt = n_envopt,
  Gridtype = 3L,
  use_parallel = FALSE,
  progbar = FALSE
)
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

summ_vec <- function(x) {
  c(mean = mean(x), sd = sd(x))
}

## Legacy coefficients: n x p matrix
leg_coef <- sim_legacy$coefficients
leg_names <- colnames(leg_coef)
leg_b1 <- leg_coef[, seq_len(l1), drop = FALSE]
leg_b2 <- leg_coef[, seq_len(l1) + l1, drop = FALSE]
cn <- colnames(leg_b1)
if (is.null(cn) || !length(cn)) {
  cn <- paste0("beta", seq_len(l1))
}
block_b1 <- t(sim_block$sim$block_results[[1]]$beta)
block_b2 <- t(sim_block$sim$block_results[[2]]$beta)
nm_vec <- c(paste0("B1.", cn), paste0("B2.", cn))
leg_stack <- cbind(leg_b1, leg_b2)
colnames(leg_stack) <- nm_vec
block_stack <- cbind(block_b1, block_b2)
colnames(block_stack) <- nm_vec

cat(sprintf("=== Coefficients (n = %d, seed = %d) ===\n\n", n_draws, seed))

fmt_row <- function(nm, leg, blk) {
  sprintf(
    "  %-12s  legacy mean %8.4f  sd %7.4f  |  block mean %8.4f  sd %7.4f  |  diff mean %+.4f",
    nm, mean(leg), sd(leg), mean(blk), sd(blk), mean(blk) - mean(leg)
  )
}

cat("Per-parameter mean and SD (B1/B2 aligned):\n")
for (nm in colnames(leg_stack)) {
  cat(fmt_row(nm, leg_stack[, nm], block_stack[, nm]), "\n")
}

cat("\n=== Dispersion ===\n")
leg_disp <- sim_legacy$dispersion
blk_disp <- sim_block$disp_out
cat(sprintf(
  "  legacy  mean = %.6f   sd = %.6f\n  block   mean = %.6f   sd = %.6f\n  diff mean = %+.6f   diff sd = %+.6f\n",
  mean(leg_disp), sd(leg_disp), mean(blk_disp), sd(blk_disp),
  mean(blk_disp) - mean(leg_disp), sd(blk_disp) - sd(leg_disp)
))

cat("\n=== Iteration counts ===\n")
cat(sprintf(
  "  legacy mean(iters) = %.3f   sd = %.3f\n  block  mean(iters) = %.3f   sd = %.3f\n",
  mean(sim_legacy$iters), sd(sim_legacy$iters),
  mean(sim_block$iters_out), sd(sim_block$iters_out)
))

cat("\nCompact table (mean / sd):\n")
tab <- data.frame(
  legacy_mean = colMeans(leg_stack),
  legacy_sd = apply(leg_stack, 2, sd),
  block_mean = colMeans(block_stack),
  block_sd = apply(block_stack, 2, sd),
  row.names = colnames(leg_stack)
)
print(round(tab, 4))

disp_tab <- data.frame(
  legacy_mean = mean(leg_disp),
  legacy_sd = sd(leg_disp),
  block_mean = mean(blk_disp),
  block_sd = sd(blk_disp),
  row.names = "dispersion"
)
print(round(disp_tab, 6))
