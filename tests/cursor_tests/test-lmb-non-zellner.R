## Non-Zellner prior: Independent Normal-Gamma with a strongly shrunk diagonal
## covariance (0.001 * diag(diag(ps$Sigma))) instead of full Prior_Setup Sigma.

test_that("rlmb: Independent Normal-Gamma with scaled diagonal Sigma (non-Zellner)", {
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)

  ps <- Prior_Setup(weight ~ group, gaussian())
  Sigma_non_zellner <- 0.001 * diag(diag(ps$Sigma))

  fit <- rlmb(
    n = 500L,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dIndependent_Normal_Gamma(
      ps$mu,
      Sigma_non_zellner,
      shape = ps$shape_ING,
      rate = ps$rate
    ),
    weights = rep(1, length(ps$y)),
    verbose = FALSE,
    use_parallel = FALSE
  )

  expect_s3_class(fit, "rlmb")
  expect_equal(nrow(fit$coefficients), 500L)
  expect_equal(ncol(fit$coefficients), nrow(ps$mu))
})

test_that("rindepNormalGamma_reg rejects ING priors with n_prior > n_w", {
  err <- tryCatch(
    rindepNormalGamma_reg(
      n = 1L,
      y = stats::rnorm(4L),
      x = matrix(1, 4L, 1L),
      prior_list = list(mu = 0, Sigma = matrix(1), shape = 10, rate = 1),
      progbar = FALSE,
      verbose = FALSE,
      use_parallel = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "n_prior <= n_w")
})

test_that("Prior_Setup Gaussian calibration: shape from n_prior, E[sigma^2|y] = dispersion", {
  ctl <- c(4.17, 5.58)
  trt <- c(4.81, 4.17)
  group <- gl(2, 2, 4)
  weight <- c(ctl, trt)
  ps <- Prior_Setup(
    weight ~ group,
    gaussian(),
    pwt = 0.01
  )
  p <- ncol(ps$x)
  n_w <- ps$PriorSettings$n_effective
  n_prior <- ps$PriorSettings$n_prior
  S_marg <- ps$dispersion * (n_w - p)
  expect_equal(ps$shape, (n_prior + 1) / 2)
  post_mean_sigma2 <- (ps$rate + S_marg / 2) / (ps$shape + n_w / 2 - 1)
  expect_equal(post_mean_sigma2, ps$dispersion)
})
