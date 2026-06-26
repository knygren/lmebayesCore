test_that("Bayesian Gaussian regression with Independent Normal-Gamma prior — OpenCL", {
  skip_if(!glmbayesCore_has_opencl(), "OpenCL not enabled in this build of glmbayesCore")
  skip_on_cran()

  library(MASS)

  data("Boston")

  predictors <- setdiff(names(Boston), "medv")
  Boston_centered <- Boston
  Boston_centered[predictors] <- scale(Boston[predictors], center = TRUE, scale = FALSE)

  form <- medv ~
    crim + zn + indus + chas + nox + age + dis + rad + tax + ptratio + black + lstat + rm

  ps <- Prior_Setup(form, gaussian(), data = Boston_centered)

  fit_normal <- rlmb(
    n = 1000,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = ps$dispersion),
    weights = rep(1, length(ps$y))
  )
  expect_s3_class(fit_normal, "rlmb")

  fit_ng <- rlmb(
    n = 1000,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dNormal_Gamma(
      mu = ps$mu,
      Sigma_0 = ps$Sigma_0,
      shape = ps$shape,
      rate = ps$rate
    ),
    weights = rep(1, length(ps$y))
  )
  expect_s3_class(fit_ng, "rlmb")

  fit_ing <- rglmb(
    n = 1000,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dIndependent_Normal_Gamma(
      ps$mu,
      ps$Sigma,
      shape = ps$shape_ING,
      rate = ps$rate
    ),
    family = gaussian(),
    weights = rep(1, length(ps$y)),
    use_parallel = TRUE,
    use_opencl = TRUE,
    verbose = FALSE
  )

  avg_candidates <- mean(fit_ing$iters)
  expect_true(avg_candidates < 400)
})
