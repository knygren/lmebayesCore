## Tests for dBeta pfamily + rBeta_reg + glmbfamfunc(binomial(identity))
## -----------------------------------------------------------------------
## Covers:
##   1. dBeta() constructor validation
##   2. Conjugate draw mean/SD match analytic posterior via rglmb()
##   3. Prior_Setup() produces non-NULL conj_binomial

test_that("dBeta() rejects bad inputs", {
  b <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  expect_error(dBeta(shape1 = -1, shape2 = 2,  beta = b), "positive")
  expect_error(dBeta(shape1 =  2, shape2 = -1, beta = b), "positive")
  expect_error(dBeta(shape1 =  0, shape2 =  2, beta = b), "positive")
  expect_error(dBeta(shape1 = "a", shape2 = 2, beta = b), "non-numeric")
  expect_error(dBeta(shape1 = c(1, 2), shape2 = 2, beta = b), "single")
})

test_that("dBeta() constructor returns valid pfamily object", {
  b  <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  pf <- dBeta(shape1 = 2, shape2 = 2, beta = b)
  expect_s3_class(pf, "pfamily")
  expect_equal(pf$pfamily, "dBeta")
  expect_true("binomial" %in% pf$okfamilies)
  expect_equal(pf$plinks(binomial(link = "identity")), "identity")
  expect_null(pf$plinks(poisson()))
  pl <- pf$prior_list
  expect_equal(as.numeric(pl$mu), 0.5)
  expect_true(as.numeric(pl$Sigma) > 0)
})

test_that("rglmb() with dBeta draws match analytic Beta posterior", {
  set.seed(101)
  n_obs     <- 40L
  theta_true <- 0.3
  y_dat  <- rbinom(n_obs, size = 1, prob = theta_true)

  alpha0 <- 3
  beta0 <- 7

  s1_post <- alpha0 + sum(y_dat)
  s2_post <- beta0  + (n_obs - sum(y_dat))
  post_mean_analytic <- s1_post / (s1_post + s2_post)
  post_sd_analytic   <- sqrt(s1_post * s2_post /
                               ((s1_post + s2_post)^2 * (s1_post + s2_post + 1)))

  b_init <- matrix(alpha0 / (alpha0 + beta0), nrow = 1L, ncol = 1L,
                   dimnames = list(NULL, "(Intercept)"))
  pf <- dBeta(shape1 = alpha0, shape2 = beta0, beta = b_init)

  ps <- Prior_Setup(
    y ~ 1,
    data = data.frame(y = y_dat),
    weights = rep(1L, n_obs),
    family = binomial(link = "identity")
  )

  set.seed(2026)
  fit <- rglmb(
    n = 30000,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = pf,
    family = binomial(link = "identity"),
    weights = rep(1L, n_obs)
  )

  smp <- fit$coefficients[, 1L]

  expect_equal(mean(smp), post_mean_analytic, tolerance = 0.006)
  expect_equal(sd(smp),   post_sd_analytic,   tolerance = 0.006)
})

test_that("Prior_Setup() produces non-NULL conj_binomial for binomial(identity)", {
  y_dat <- c(rep(1, 7), rep(0, 18))
  df    <- data.frame(y = y_dat)

  ps <- Prior_Setup(
    y ~ 1,
    data = df,
    weights = rep(1L, nrow(df)),
    family = binomial(link = "identity"),
    pwt = 0.05
  )

  expect_false(is.null(ps$conj_binomial))
  cb <- ps$conj_binomial
  expect_true(cb$shape1 > 0)
  expect_true(cb$shape2 > 0)
  expect_equal(cb$shape1 / (cb$shape1 + cb$shape2),
               cb$weighted_mean_prop, tolerance = 1e-10)
})
