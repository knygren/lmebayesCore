menarche_opencl_ag_fit <- function(link) {
  data("menarche", package = "MASS")
  menarche$Age2 <- menarche$Age - 13
  menarche$prop <- menarche$Menarche / menarche$Total
  fam <- binomial(link = link)
  ps <- Prior_Setup(
    prop ~ Age2,
    family = fam,
    data = menarche,
    weights = menarche$Total
  )
  rglmb(
    n = 200,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    family = fam,
    weights = menarche$Total,
    Gridtype = 2,
    use_parallel = TRUE,
    use_opencl = TRUE,
    verbose = FALSE
  )
}

test_that("OpenCL binomial logit with MASS menarche (proportion + trial weights)", {
  skip_if_no_opencl()
  skip_on_cran()

  fit <- menarche_opencl_ag_fit("logit")
  expect_s3_class(fit, "rglmb")
})

test_that("OpenCL binomial probit with MASS menarche (proportion + trial weights)", {
  skip_if_no_opencl()
  skip_on_cran()

  fit <- menarche_opencl_ag_fit("probit")
  expect_s3_class(fit, "rglmb")
})

test_that("OpenCL binomial cloglog with MASS menarche (proportion + trial weights)", {
  skip_if_no_opencl()
  skip_on_cran()

  fit <- suppressWarnings(menarche_opencl_ag_fit("cloglog"))
  expect_s3_class(fit, "rglmb")
})
