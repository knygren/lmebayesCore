# inst/examples/Ex_rglmb.R (Dobson Poisson block); only use_opencl = TRUE added

test_that("OpenCL f2_f3_poisson (Dobson RCT)", {
  skip_if_no_opencl()
  skip_on_cran()

  set.seed(333)
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  d.AD <- data.frame(treatment, outcome, counts)

  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson(), data = d.AD)

  fit <- rglmb(
    n = 1000,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    family = poisson(),
    weights = rep(1, nrow(ps$x)),
    use_opencl = TRUE
  )

  expect_s3_class(fit, "rglmb")
})
