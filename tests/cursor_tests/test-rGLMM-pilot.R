library(glmbayesCore)

# Minimal random-intercept fixture for matrix-level rGLMM tests
.mini_rGLMM_inputs <- function(family = gaussian()) {
  set.seed(1L)
  J <- 4L
  n_per <- 5L
  g <- factor(rep(seq_len(J), each = n_per))
  n <- length(g)
  re_nm <- "(Intercept)"
  if (identical(family$family, "gaussian")) {
    y <- stats::rnorm(n)
  } else {
    y <- stats::rpois(n, lambda = 2)
  }
  x <- matrix(1, nrow = n, ncol = 1L, dimnames = list(NULL, re_nm))
  gl <- as.character(seq_len(J))
  x_hyper <- stats::setNames(
    list(matrix(1, J, 1L, dimnames = list(gl, re_nm))),
    re_nm
  )
  pf <- dNormal(mu = c(`(Intercept)` = 0), Sigma = matrix(1), dispersion = 1)
  pfamily_list <- stats::setNames(list(pf), re_nm)
  if (identical(family$family, "gaussian")) {
    prior_list <- list(P = matrix(1), dispersion = 1)
  } else {
    prior_list <- list(P = matrix(1))
  }
  list(
    y = y, x = x, block = g, x_hyper = x_hyper,
    prior_list = prior_list,
    pfamily_list = pfamily_list,
    re_names = re_nm,
    group_levels = gl,
    family = family
  )
}

test_that("rGLMM: Gaussian skips pilot even with default gap_tol", {
  inp <- .mini_rGLMM_inputs(gaussian())
  out <- rGLMM_reg(
    n = 2L,
    y = inp$y,
    x = inp$x,
    block = inp$block,
    x_hyper = inp$x_hyper,
    prior_list = inp$prior_list,
    pfamily_list = inp$pfamily_list,
    family = inp$family,
    n_pilot = NULL,
    gap_tol = 0.0196,
    verbose = FALSE,
    progbar = FALSE
  )
  expect_equal(out$n_pilot, 0L)
  expect_null(out$pilot)
  expect_false(is.null(out$fixef))
})

test_that("rGLMM: Poisson runs pilot when n_pilot is explicit", {
  inp <- .mini_rGLMM_inputs(poisson())
  out <- rGLMM_reg(
    n = 2L,
    y = inp$y,
    x = inp$x,
    block = inp$block,
    x_hyper = inp$x_hyper,
    prior_list = inp$prior_list,
    pfamily_list = inp$pfamily_list,
    family = inp$family,
    n_pilot = 2L,
    gap_tol = 0.0196,
    verbose = FALSE,
    progbar = FALSE
  )
  expect_equal(out$n_pilot, 2L)
  expect_false(is.null(out$pilot))
  expect_false(is.null(out$pilot_chisq))
  expect_false(is.null(out$fixef))
  expect_s3_class(out$sweep_history, "two_block_sweep_history")
  expect_s3_class(out$pilot$sweep_history, "two_block_sweep_history")
})

test_that("rGLMM: Poisson skips pilot when n_pilot = 0L", {
  inp <- .mini_rGLMM_inputs(poisson())
  out <- rGLMM_reg(
    n = 2L,
    y = inp$y,
    x = inp$x,
    block = inp$block,
    x_hyper = inp$x_hyper,
    prior_list = inp$prior_list,
    pfamily_list = inp$pfamily_list,
    family = inp$family,
    n_pilot = 0L,
    gap_tol = 0.0196,
    verbose = FALSE,
    progbar = FALSE
  )
  expect_equal(out$n_pilot, 0L)
  expect_null(out$pilot)
})

test_that("rGLMM: Poisson skips pilot when tv_tol and gap_tol are NULL", {
  inp <- .mini_rGLMM_inputs(poisson())
  out <- rGLMM_reg(
    n = 2L,
    y = inp$y,
    x = inp$x,
    block = inp$block,
    x_hyper = inp$x_hyper,
    prior_list = inp$prior_list,
    pfamily_list = inp$pfamily_list,
    family = inp$family,
    n_pilot = NULL,
    gap_tol = NULL,
    tv_tol = NULL,
    verbose = FALSE,
    progbar = FALSE
  )
  expect_equal(out$n_pilot, 0L)
  expect_null(out$pilot)
})

test_that("rGLMM: Poisson derives n_pilot from cost optimization when tv_tol set", {
  inp <- .mini_rGLMM_inputs(poisson())
  out <- rGLMM_reg(
    n = 100L,
    y = inp$y,
    x = inp$x,
    block = inp$block,
    x_hyper = inp$x_hyper,
    prior_list = inp$prior_list,
    pfamily_list = inp$pfamily_list,
    family = inp$family,
    n_pilot = NULL,
    gap_tol = 0.0196,
    verbose = FALSE,
    progbar = FALSE
  )
  expect_true(out$n_pilot > 0L)
  expect_true(out$n_pilot < 10000L)
  expect_equal(out$convergence$n_pilot_source, "cost")
  expect_equal(out$n_pilot, out$convergence$pilot_cost_opt$n_pilot_opt)
  expect_equal(out$m_convergence, out$convergence$pilot_cost_opt$m_convergence_opt)
})
