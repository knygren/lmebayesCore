## Timing breakdown: build vs sample
suppressPackageStartupMessages(library(glmbayesCore))

fix <- local({
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group1  <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight1 <- c(ctl, trt); n1 <- length(weight1)
  weight  <- c(weight1, weight1)
  df_stacked <- data.frame(weight = weight,
                           group  = factor(c(group1, group1), levels = levels(group1)))
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
    max_disp_perc = 0.99)
  prior_list_block <- list(mu = p_setup$mu, Sigma = p_setup$Sigma,
    shape = p_setup$shape, rate = p_setup$rate, max_disp_perc = 0.99)
  list(y = y, x_block = x_block, x_old = x_old, block = block,
       prior_list_old = prior_list_old, prior_list_block = prior_list_block,
       l1 = l1, n_obs = length(y))
})

n_envopt <- 10000L
block_base <- list(
  y = fix$y, x = fix$x_block, block = fix$block,
  prior_list = fix$prior_list_block, prior_lists = NULL,
  offset = rep(0, fix$n_obs), wt = rep(1, fix$n_obs),
  p_re = -1L, n_rss_iter = 10L, Gridtype = 3L, n_envopt = n_envopt,
  RSS_ML = NA_real_, use_parallel = FALSE, use_opencl = FALSE,
  progbar = FALSE, verbose = FALSE,
  group_levels = character(0), re_names = character(0))

## ---- Isolate build phase (n=1 draw) ----
set.seed(1)
t_leg_build  <- system.time(do.call(rindepNormalGamma_reg,
  c(list(n=1L, y=fix$y, x=fix$x_old, prior_list=fix$prior_list_old,
         n_envopt=n_envopt, Gridtype=3L, use_parallel=FALSE, progbar=FALSE))))

set.seed(1)
t_jnt_build  <- system.time(do.call(glmbayesCore:::.rIndepNormalGammaRegBlock_cpp,
  c(list(n=1L), block_base)))

set.seed(1)
t_ind_build  <- system.time(do.call(glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp,
  c(list(n=1L), block_base)))

## ---- Full 10k run ----
set.seed(2026)
t_leg_full <- system.time(do.call(rindepNormalGamma_reg,
  c(list(n=10000L, y=fix$y, x=fix$x_old, prior_list=fix$prior_list_old,
         n_envopt=n_envopt, Gridtype=3L, use_parallel=FALSE, progbar=FALSE))))

set.seed(2026)
t_jnt_full <- system.time(do.call(glmbayesCore:::.rIndepNormalGammaRegBlock_cpp,
  c(list(n=10000L), block_base)))

set.seed(2026)
t_ind_full <- system.time(do.call(glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp,
  c(list(n=10000L), block_base)))

cat("\n=== Timing breakdown (elapsed seconds) ===\n")
cat(sprintf("%-12s  %8s  %8s  %8s\n", "Phase", "Legacy", "Joint", "Ind"))
cat(sprintf("%-12s  %8.3f  %8.3f  %8.3f\n",
  "Build (n=1)",  t_leg_build["elapsed"], t_jnt_build["elapsed"], t_ind_build["elapsed"]))
cat(sprintf("%-12s  %8.3f  %8.3f  %8.3f\n",
  "Full (n=10k)", t_leg_full["elapsed"],  t_jnt_full["elapsed"],  t_ind_full["elapsed"]))

sample_leg <- t_leg_full["elapsed"] - t_leg_build["elapsed"]
sample_jnt <- t_jnt_full["elapsed"] - t_jnt_build["elapsed"]
sample_ind <- t_ind_full["elapsed"] - t_ind_build["elapsed"]
cat(sprintf("%-12s  %8.3f  %8.3f  %8.3f\n",
  "~Sample only", sample_leg, sample_jnt, sample_ind))
cat(sprintf("\nPer-draw sampling cost (ms): Legacy=%.3f  Joint=%.3f  Ind=%.3f\n",
  sample_leg*1000/10000, sample_jnt*1000/10000, sample_ind*1000/10000))
