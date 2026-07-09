## Three-way 10k comparison: legacy vs joint-block vs independent-block
## Run via: Rscript data-raw/three_way_10k.R
suppressPackageStartupMessages({
  library(glmbayesCore)
  if (requireNamespace("glmbayes", quietly = TRUE)) library(glmbayes)
})

fix <- local({
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group1  <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight1 <- c(ctl, trt); n1 <- length(weight1)
  weight  <- c(weight1, weight1)
  group_stacked <- factor(c(group1, group1), levels = levels(group1))
  df_stacked <- data.frame(weight = weight, group = group_stacked)
  p_setup <- Prior_Setup(weight ~ group, data = df_stacked, family = gaussian())
  y <- p_setup$y; x_block <- p_setup$x; l1 <- ncol(x_block)
  block <- factor(c(rep("B1", n1), rep("B2", n1)), levels = c("B1", "B2"))
  x1    <- x_block[seq_len(n1), , drop = FALSE]
  zeros <- matrix(0, n1, l1)
  x_old <- rbind(cbind(x1, zeros), cbind(zeros, x1))
  prior_list_old <- list(
    mu   = c(p_setup$mu, p_setup$mu),
    Sigma = rbind(cbind(p_setup$Sigma, matrix(0, l1, l1)),
                  cbind(matrix(0, l1, l1), p_setup$Sigma)),
    dispersion = p_setup$dispersion,
    shape = p_setup$shape, rate = p_setup$rate,
    Precision = solve(rbind(cbind(p_setup$Sigma, matrix(0, l1, l1)),
                            cbind(matrix(0, l1, l1), p_setup$Sigma))),
    max_disp_perc = 0.99
  )
  prior_list_block <- list(
    mu = p_setup$mu, Sigma = p_setup$Sigma,
    shape = p_setup$shape, rate = p_setup$rate, max_disp_perc = 0.99
  )
  list(y = y, x_block = x_block, x_old = x_old, block = block,
       prior_list_old = prior_list_old, prior_list_block = prior_list_block,
       l1 = l1, n_obs = length(y))
})

n_draws  <- 10000L
n_envopt <- 10000L

block_args <- list(
  y = fix$y, x = fix$x_block, block = fix$block,
  prior_list = fix$prior_list_block, prior_lists = NULL,
  offset = rep(0, fix$n_obs), wt = rep(1, fix$n_obs),
  p_re = -1L, n_rss_iter = 10L, Gridtype = 3L, n_envopt = n_envopt,
  RSS_ML = NA_real_, use_parallel = FALSE, use_opencl = FALSE,
  progbar = FALSE, verbose = FALSE,
  group_levels = character(0), re_names = character(0)
)

cat("Running legacy (n =", n_draws, ")...\n")
set.seed(2026)
t_leg <- system.time(sim_leg <- rindepNormalGamma_reg(
  n = n_draws, y = fix$y, x = fix$x_old,
  prior_list = fix$prior_list_old,
  n_envopt = n_envopt, Gridtype = 3L, use_parallel = FALSE, progbar = FALSE
))

cat("Running joint block...\n")
set.seed(2026)
t_jnt <- system.time(sim_jnt <- do.call(
  glmbayesCore:::.rIndepNormalGammaRegBlock_cpp,
  c(list(n = n_draws), block_args)
))

cat("Running independent block...\n")
set.seed(2026)
t_ind <- system.time(sim_ind <- do.call(
  glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp,
  c(list(n = n_draws), block_args)
))

leg_b <- sim_leg$coefficients          # n x p
jnt_b <- rbind(sim_jnt$sim$block_results[[1]]$beta,
               sim_jnt$sim$block_results[[2]]$beta)  # p x n
ind_b <- rbind(sim_ind$sim$block_results[[1]]$beta,
               sim_ind$sim$block_results[[2]]$beta)
leg_d <- sim_leg$dispersion
jnt_d <- sim_jnt$disp_out
ind_d <- sim_ind$disp_out
p <- ncol(leg_b)
bnames <- c("(Intercept)_B1", "TrtTrt_B1", "(Intercept)_B2", "TrtTrt_B2")

cat("\n=== Sampling times (elapsed seconds) ===\n")
cat(sprintf("  Legacy : %.2f s\n", t_leg["elapsed"]))
cat(sprintf("  Joint  : %.2f s\n", t_jnt["elapsed"]))
cat(sprintf("  Ind    : %.2f s\n", t_ind["elapsed"]))

cat("\n=== Acceptance (mean iters per draw) ===\n")
cat(sprintf("  Legacy : %.3f   (accept rate = 1/%.3f = %.4f)\n",
  mean(sim_leg$iters), mean(sim_leg$iters), 1/mean(sim_leg$iters)))
cat(sprintf("  Joint  : %.3f   (accept rate = %.4f)\n",
  mean(sim_jnt$iters_out), 1/mean(sim_jnt$iters_out)))
cat(sprintf("  Ind    : %.3f   (accept rate = %.4f)\n",
  mean(sim_ind$iters_out), 1/mean(sim_ind$iters_out)))
cat(sprintf("  Ind/Joint ratio: %.3f\n",
  mean(sim_ind$iters_out) / mean(sim_jnt$iters_out)))

cat("\n=== Coefficient means and SDs ===\n")
cat(sprintf("%-20s  %8s %8s %8s  |  %8s %8s %8s\n",
  "Param", "Leg_mu", "Jnt_mu", "Ind_mu", "Leg_sd", "Jnt_sd", "Ind_sd"))
for (i in seq_len(p)) {
  cat(sprintf("%-20s  %8.4f %8.4f %8.4f  |  %8.4f %8.4f %8.4f\n",
    bnames[i],
    mean(leg_b[, i]), mean(jnt_b[i, ]), mean(ind_b[i, ]),
    sd(leg_b[, i]),   sd(jnt_b[i, ]),   sd(ind_b[i, ])))
}
cat(sprintf("%-20s  %8.4f %8.4f %8.4f  |  %8.4f %8.4f %8.4f\n",
  "sigma2",
  mean(leg_d), mean(jnt_d), mean(ind_d),
  sd(leg_d),   sd(jnt_d),   sd(ind_d)))

cat("\n=== t-tests: means (joint vs ind) ===\n")
for (i in seq_len(p)) {
  tt <- t.test(jnt_b[i, ], ind_b[i, ])
  cat(sprintf("  %-20s  t=%7.3f  p=%6.4f  diff=%+.5f  95CI=[%+.5f, %+.5f]\n",
    bnames[i], tt$statistic, tt$p.value,
    diff(tt$estimate), tt$conf.int[1], tt$conf.int[2]))
}
tt_d <- t.test(jnt_d, ind_d)
cat(sprintf("  %-20s  t=%7.3f  p=%6.4f  diff=%+.5f  95CI=[%+.5f, %+.5f]\n",
  "sigma2", tt_d$statistic, tt_d$p.value,
  diff(tt_d$estimate), tt_d$conf.int[1], tt_d$conf.int[2]))

cat("\n=== t-tests: means (legacy vs joint) ===\n")
for (i in seq_len(p)) {
  tt <- t.test(leg_b[, i], jnt_b[i, ])
  cat(sprintf("  %-20s  t=%7.3f  p=%6.4f  diff=%+.5f\n",
    bnames[i], tt$statistic, tt$p.value, diff(tt$estimate)))
}
tt_d2 <- t.test(leg_d, jnt_d)
cat(sprintf("  %-20s  t=%7.3f  p=%6.4f  diff=%+.5f\n",
  "sigma2", tt_d2$statistic, tt_d2$p.value, diff(tt_d2$estimate)))

cat("\n=== F-tests: variance ratios (joint vs ind) ===\n")
for (i in seq_len(p)) {
  ft <- var.test(jnt_b[i, ], ind_b[i, ])
  cat(sprintf("  %-20s  F=%7.4f  p=%6.4f  var_ratio=%7.4f  95CI=[%.4f, %.4f]\n",
    bnames[i], ft$statistic, ft$p.value,
    ft$estimate, ft$conf.int[1], ft$conf.int[2]))
}
ft_d <- var.test(jnt_d, ind_d)
cat(sprintf("  %-20s  F=%7.4f  p=%6.4f  var_ratio=%7.4f  95CI=[%.4f, %.4f]\n",
  "sigma2", ft_d$statistic, ft_d$p.value,
  ft_d$estimate, ft_d$conf.int[1], ft_d$conf.int[2]))
