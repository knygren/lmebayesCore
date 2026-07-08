# Assess sources of remaining legacy vs block acceptance gap (two-block fixture)
pkgload::load_all(quiet = TRUE)

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

sim_env <- rindepNormalGamma_reg_with_envelope(
  n = 1L, y = y, x = x_old, prior_list = prior_list_old,
  n_envopt = n_envopt, Gridtype = 3L, use_parallel = FALSE, progbar = FALSE
)
sim_block <- glmbayesCore:::.rIndepNormalGammaRegBlock_cpp(
  n = 1L, y = y, x = x_block, block = block,
  prior_list = prior_list_block, n_envopt = n_envopt, Gridtype = 3L,
  use_parallel = FALSE, progbar = FALSE, verbose = FALSE,
  offset = rep(0, length(y)), wt = rep(1, length(y)),
  p_re = -1L, n_rss_iter = 10L, RSS_ML = NA_real_, use_opencl = FALSE,
  group_levels = character(0), re_names = character(0)
)

lgL <- sim_env$UB_list$lg_prob_factor
ub2L <- sim_env$UB_list$UB2min
de <- sim_block$build_out$dispersion_envelope
bd <- de$block_dispersion
max_upp <- de$prob_max_upp
max_low <- de$prob_max_low

bd1 <- bd[[1]]
bd2 <- bd[[2]]
gs <- nrow(sim_block$build_out$block_envelopes[[1]]$cbars)

# All 81 product-face constants from block shortcuts
ub2_joint <- lg_joint <- numeric(gs * gs)
idx <- 0L
for (j1 in seq_len(gs)) {
  for (j2 in seq_len(gs)) {
    idx <- idx + 1L
    ub2_low <- bd1$ub2_at_low[j1] + bd2$ub2_at_low[j2]
    ub2_upp <- bd1$ub2_at_upp[j1] + bd2$ub2_at_upp[j2]
    ub2_joint[idx] <- min(ub2_low, ub2_upp)
    upp_sum <- bd1$upp_apprx[j1] + bd2$upp_apprx[j2]
    low_sum <- bd1$low_apprx[j1] + bd2$low_apprx[j2]
    lg_joint[idx] <- max(upp_sum - max_upp, low_sum - max_low)
  }
}

summ_diff <- function(x, y, label) {
  d <- x - y
  cat(sprintf(
    "\n=== %s (n=%d) ===\n  mean diff = %+.6g  max |diff| = %.6g  cor(sorted) = %.5f\n",
    label, length(d), mean(d), max(abs(d)), cor(sort(x), sort(y))
  ))
  print(summary(d))
  invisible(d)
}

cat("=== Methodology note ===\n")
cat(
  "Legacy iters = attempts until accept (resample loop).\n",
  "Block iters_out uses the same convention (starts at 1, increments on reject).\n",
  "Same seed still does not sync proposal streams (joint vs product PLSD).\n\n"
)

cat("=== Global constants (already matched in tests) ===\n")
gL <- sim_env$gamma_list
gB <- de$gamma_list
uL <- sim_env$UB_list
uB <- de$UB_list
keys <- c(
  "lmc1", "lmc2", "lm_log1", "lm_log2", "max_New_LL_UB",
  "max_LL_log_disp", "RSS_Min", "shape3", "rate2"
)
for (k in keys) {
  vL <- if (k %in% names(gL)) gL[[k]] else uL[[k]]
  vB <- if (k %in% names(gB)) gB[[k]] else uB[[k]]
  cat(sprintf("  %-16s legacy %.6g  block %.6g  diff %+.3g\n", k, vL, vB, vB - vL))
}

summ_diff(lg_joint, lgL, "UB3A lg shortcut vs legacy lg_prob_factor (81 product faces)")
summ_diff(ub2_joint, ub2L, "UB2 joint min-of-sums vs legacy UB2min")

# Per-block UB2min sum (old bug) for reference
ub2_sum_bug <- outer(bd1$UB2min, bd2$UB2min, "+")
summ_diff(as.vector(ub2_sum_bug), ub2L, "UB2min SUM of per-block minima (pre-fix bug)")

cat("\n=== Acceptance at n=10000 (seed 360) ===\n")
n_draws <- 10000L
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
legacy_dpa <- mean(sim_legacy$iters)
block_dpa <- mean(sim_block$iters_out)
legacy_rate <- 1 / legacy_dpa
block_rate <- 1 / block_dpa
cat(sprintf(
  "  legacy 1/mean(iters) = %.5f (%.2f%%)\n  block 1/mean(iters) = %.5f (%.2f%%)\n  mean(iters) ratio block/legacy = %.4f  gap = %+.2f\n",
  legacy_rate, 100 * legacy_rate, block_rate, 100 * block_rate,
  block_dpa / legacy_dpa, block_dpa - legacy_dpa
))
legacy_first <- mean(sim_legacy$iters == 1)
block_first <- mean(sim_block$iters_out == 1)
cat(sprintf(
  "  legacy P(iters==1) = %.5f  block P(iters==1) = %.5f\n",
  legacy_first, block_first
))

cat("\n=== PLSD / proposal distribution ===\n")
plsdL <- sim_env$Envelope$PLSD
plsd1 <- sim_block$build_out$block_envelopes[[1]]$PLSD
plsd2 <- sim_block$build_out$block_envelopes[[2]]$PLSD
cat(sprintf(
  "  legacy joint faces = %d  block gs per block = %d, %d\n",
  length(plsdL), length(plsd1), length(plsd2)
))
cat(sprintf(
  "  legacy max PLSD = %.4f  block max PLSD = %.4f, %.4f\n",
  max(plsdL), max(plsd1), max(plsd2)
))
# Product PLSD if independent
plsd_prod <- as.vector(outer(plsd1, plsd2, "*"))
cat(sprintf(
  "  product PLSD entropy ratio H(block prod)/H(legacy) = %.4f\n",
  {
    h <- function(p) -sum(p[p > 0] * log(p[p > 0]))
    h(plsd_prod) / h(plsdL)
  }
))

cat("\n=== Likely residual gap sources (ranked) ===\n")
cat("1. Different proposal RNG paths (joint vs product face draws) — not same proposals.\n")
cat("2. UB2 joint uses rss_face_at_disp endpoints; legacy UB2min uses M_min/M_max (method 2).\n")
cat("3. lg_prob_factor shortcut vs legacy joint prob_factor — small sorted-vector drift.\n")
cat("4. Per-block vs joint envelope build (cbars / face geometry not byte-identical).\n")
