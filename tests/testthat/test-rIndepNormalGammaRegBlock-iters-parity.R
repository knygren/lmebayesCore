## k = 1 parity: BlockEnvelope iters vs rindepNormalGamma_reg accept/reject counts.
## Fixture from inst/examples/Ex_rindepNormalGamma_reg.R.

test_that("BlockEnvelope iters_out matches rindepNormalGamma_reg draws per acceptance (k = 1)", {
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  p_setup <- Prior_Setup(weight ~ group, family = gaussian())

  y <- p_setup$y
  x <- p_setup$x
  n_draws <- 10000L

  prior_list_old <- list(
    mu = p_setup$mu,
    Sigma = p_setup$Sigma,
    dispersion = p_setup$dispersion,
    shape = p_setup$shape,
    rate = p_setup$rate,
    Precision = solve(p_setup$Sigma),
    max_disp_perc = 0.99
  )

  prior_list_block <- prior_list_old
  prior_list_block$dispersion <- NULL

  set.seed(360)
  sim_old <- rindepNormalGamma_reg(
    n = n_draws,
    y = y,
    x = x,
    prior_list = prior_list_old,
    progbar = FALSE
  )

  set.seed(360)
  sim_block <- glmbayesCore:::.rIndepNormalGammaRegBlock_cpp(
    n = n_draws,
    y = y,
    x = x,
    block = factor(rep(1L, length(y))),
    prior_list = prior_list_block,
    prior_lists = NULL,
    offset = rep(0, length(y)),
    wt = rep(1, length(y)),
    p_re = -1L,
    n_rss_iter = 10L,
    Gridtype = 2L,
    n_envopt = n_draws,
    RSS_ML = NA_real_,
    use_parallel = TRUE,
    use_opencl = FALSE,
    progbar = FALSE,
    verbose = FALSE,
    group_levels = character(0),
    re_names = character(0)
  )

  old_draws_per_accept <- mean(sim_old$iters)
  block_iters <- sim_block$iters_out
  block_draws_per_accept <- mean(block_iters)

  expect_length(block_iters, n_draws)
  expect_true(all(block_iters >= 1))
  expect_equal(sim_block$iters_mean, block_iters[[1L]])
  expect_equal(
    sim_block$sim$meta$accept_mode,
    "resample_until_accept_v1"
  )
  expect_true(is.finite(old_draws_per_accept))
  expect_true(old_draws_per_accept >= 1)
  expect_true(is.finite(block_draws_per_accept))
  expect_true(block_draws_per_accept >= 1)

  ## Both paths: iters_out counts candidates until acceptance (starts at 1).
  expect_equal(
    block_draws_per_accept,
    old_draws_per_accept,
    tolerance = 0.3
  )
})
