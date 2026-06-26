test_that("OpenCL f2_f3_gamma (carinsca)", {
  skip_if_no_opencl()
  skip_on_cran()

  data(carinsca)
  carinsca$Merit <- ordered(carinsca$Merit)
  carinsca$Class <- factor(carinsca$Class)
  oldopt <- options(contrasts = c("contr.treatment", "contr.treatment"))
  on.exit(options(oldopt), add = TRUE)

  out <- glm(
    Cost / Claims ~ Merit + Class,
    family = Gamma(link = "log"),
    data = carinsca,
    weights = carinsca$Claims,
    x = TRUE
  )
  disp <- gamma.dispersion(out)
  ps <- Prior_Setup(
    Cost / Claims ~ Merit + Class,
    family = Gamma(link = "log"),
    data = carinsca,
    weights = carinsca$Claims
  )

  fit <- rglmb(
    n = 1000,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = disp),
    family = Gamma(link = "log"),
    weights = carinsca$Claims,
    use_opencl = TRUE
  )

  expect_s3_class(fit, "rglmb")
})
