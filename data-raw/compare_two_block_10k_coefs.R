# Mean coefficients: two-block stacked Dobson, n = 10000
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

legacy_mean <- colMeans(sim_legacy$coefficients)
cat("Legacy mean coefficients (all accepted draws, n =", n_draws, "):\n")
print(round(legacy_mean, 4))
cat("\nLegacy names:", paste(colnames(sim_legacy$coefficients), collapse = ", "), "\n\n")

b1 <- sim_block$sim$block_results[[1]]$beta
b2 <- sim_block$sim$block_results[[2]]$beta
cn <- rownames(b1)
if (is.null(cn)) {
  cn <- paste0("beta", seq_len(nrow(b1)))
}
block_mean_all <- c(rowMeans(b1), rowMeans(b2))
names(block_mean_all) <- c(paste0("B1.", cn), paste0("B2.", cn))
cat("Block mean coefficients (all proposals, n =", n_draws, "):\n")
print(round(block_mean_all, 4))

acc <- sim_block$iters_out == 1
cat(
  "\nBlock mean coefficients (would-accept only, n_acc =",
  sum(acc), "):\n"
)
block_mean_acc <- c(
  rowMeans(b1[, acc, drop = FALSE]),
  rowMeans(b2[, acc, drop = FALSE])
)
names(block_mean_acc) <- names(block_mean_all)
print(round(block_mean_acc, 4))

leg_names <- colnames(sim_legacy$coefficients)
leg_stack <- legacy_mean
names(leg_stack) <- c(
  paste0("B1.", leg_names[seq_len(l1)]),
  paste0("B2.", leg_names[seq_len(l1) + l1])
)
cmp <- cbind(
  legacy = leg_stack,
  block_all = block_mean_all,
  block_acc = block_mean_acc
)
cat("\nSide-by-side (legacy block-diagonal vs block B1/B2):\n")
print(round(cmp, 4))
