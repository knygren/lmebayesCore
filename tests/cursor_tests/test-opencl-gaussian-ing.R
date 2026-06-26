test_that("OpenCL f2_f3_gaussian ING (Boston centered predictors)", {
  skip_if_no_opencl()
  skip_on_cran()

  data("Boston_centered")

  form <- medv ~
    crim + zn + indus + chas + nox + age + dis + rad + tax + ptratio + black + lstat + rm

  ps <- Prior_Setup(form, gaussian(), data = Boston_centered)

  fit <- rglmb(
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
    use_opencl = TRUE
  )

  expect_s3_class(fit, "rglmb")
})
