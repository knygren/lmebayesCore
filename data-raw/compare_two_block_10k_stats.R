# Statistical comparison: legacy vs block (n = 10000, seed 360)
# Two independent samples (same seed start, divergent RNG paths).
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

leg <- sim_legacy$coefficients
cn <- colnames(leg)[seq_len(l1)]
nm <- c(paste0("B1.", cn), paste0("B2.", cn))
leg_stack <- cbind(leg[, seq_len(l1)], leg[, seq_len(l1) + l1])
colnames(leg_stack) <- nm
blk_stack <- cbind(
  t(sim_block$sim$block_results[[1]]$beta),
  t(sim_block$sim$block_results[[2]]$beta)
)
colnames(blk_stack) <- nm

## Normal-iid exact CI: (n-1)s^2 / sigma^2 ~ ChiSq_{n-1}
sd_ci_chisq <- function(x, conf = 0.95) {
  n <- length(x)
  v <- var(x)
  alpha <- 1 - conf
  lo <- sqrt(v * (n - 1) / qchisq(1 - alpha / 2, n - 1))
  hi <- sqrt(v * (n - 1) / qchisq(alpha / 2, n - 1))
  c(sd = sqrt(v), sd_lo = lo, sd_hi = hi, var = v,
    var_lo = v * (n - 1) / qchisq(1 - alpha / 2, n - 1),
    var_hi = v * (n - 1) / qchisq(alpha / 2, n - 1))
}

## Ratio of variances CI under independent normals: F = s1^2/s2^2
var_ratio_test <- function(x, y, conf = 0.95) {
  n1 <- length(x)
  n2 <- length(y)
  v1 <- var(x)
  v2 <- var(y)
  ratio <- v1 / v2
  alpha <- 1 - conf
  c(
    ratio = ratio,
    ratio_lo = ratio / qf(1 - alpha / 2, n1 - 1, n2 - 1),
    ratio_hi = ratio / qf(alpha / 2, n1 - 1, n2 - 1),
    p_f = var.test(x, y)$p.value
  )
}

compare_pair <- function(x, y, label) {
  tt <- t.test(x, y, var.equal = FALSE)
  vr <- var_ratio_test(x, y)
  cx <- sd_ci_chisq(x)
  cy <- sd_ci_chisq(y)
  data.frame(
    param = label,
    legacy_mean = mean(x),
    block_mean = mean(y),
    mean_diff = mean(y) - mean(x),
    mean_p_welch = tt$p.value,
    mean_ci_diff_lo = tt$conf.int[1],
    mean_ci_diff_hi = tt$conf.int[2],
    legacy_sd = cx["sd"],
    legacy_sd_lo = cx["sd_lo"],
    legacy_sd_hi = cx["sd_hi"],
    block_sd = cy["sd"],
    block_sd_lo = cy["sd_lo"],
    block_sd_hi = cy["sd_hi"],
    var_ratio = vr["ratio"],
    var_ratio_lo = vr["ratio_lo"],
    var_ratio_hi = vr["ratio_hi"],
    var_p_f = vr["p_f"],
    stringsAsFactors = FALSE
  )
}

cat(sprintf(
  "=== Legacy vs block: statistical comparison (n = %d, seed = %d) ===\n\n",
  n_draws, seed
))
cat(
  "Design: two independent MCMC/envelope samples (not paired draw-by-draw).\n",
  "Mean test: Welch two-sample t-test.\n",
  "Variance test: F-test (var.test).\n",
  "SD CI: exact 95%% chi-square interval under i.i.d. Normal (diagnostic).\n",
  "Variance-ratio CI: 95%% F interval for legacy/block variance ratio.\n\n"
)

rows <- lapply(nm, function(j) compare_pair(leg_stack[, j], blk_stack[, j], j))
tab <- do.call(rbind, rows)
tab_num <- tab
tab_num$param <- NULL
print(cbind(param = tab$param, round(tab_num, 5)), row.names = FALSE)

xd <- sim_legacy$dispersion
yd <- sim_block$disp_out
disp_row <- compare_pair(xd, yd, "dispersion")
cat("\n=== Dispersion row ===\n")
disp_num <- disp_row
disp_num$param <- NULL
print(round(t(disp_num), 3))

cat("\n=== Interpretation at alpha = 0.05 ===\n")
all_p_mean <- c(tab$mean_p_welch, disp_row$mean_p_welch)
all_p_var <- c(tab$var_p_f, disp_row$var_p_f)
cat(sprintf(
  "Mean differences significant: %d / %d\n",
  sum(all_p_mean < 0.05), length(all_p_mean)
))
cat(sprintf(
  "Variance differences significant: %d / %d\n",
  sum(all_p_var < 0.05), length(all_p_var)
))
cat(sprintf(
  "Bonferroni (5 params): mean sig if p < 0.01 -> %d; var sig if p < 0.01 -> %d\n",
  sum(all_p_mean < 0.01), sum(all_p_var < 0.01)
))

cat("\n=== Do SD CIs overlap between samplers? ===\n")
for (i in seq_len(nrow(tab))) {
  ov <- !(tab$legacy_sd_hi[i] < tab$block_sd_lo[i] ||
            tab$block_sd_hi[i] < tab$legacy_sd_lo[i])
  cat(sprintf(
    "  %s: SD CIs overlap = %s (legacy [%.4f,%.4f], block [%.4f,%.4f])\n",
    tab$param[i], ov,
    tab$legacy_sd_lo[i], tab$legacy_sd_hi[i],
    tab$block_sd_lo[i], tab$block_sd_hi[i]
  ))
}
ov_d <- !(disp_row$legacy_sd_hi < disp_row$block_sd_lo ||
            disp_row$block_sd_hi < disp_row$legacy_sd_lo)
cat(sprintf(
  "  dispersion: SD CIs overlap = %s\n",
  ov_d
))
