.mock_two_block_rate <- function(lambda_star, eigenvalues) {
  structure(
    list(lambda_star = lambda_star, eigenvalues = eigenvalues),
    class = "two_block_rate"
  )
}

test_that("two_block_d0_pilot_start decreases in n_pilot", {
  p <- 5L
  d1 <- two_block_d0_pilot_start(10L, p)
  d2 <- two_block_d0_pilot_start(100L, p)
  expect_true(d2 < d1)
  expect_equal(d1, stats::qchisq(0.95, df = p) / 10)
})

test_that("two_block_m_convergence_for_pilot_start is monotone in n_pilot", {
  rate <- .mock_two_block_rate(0.85, c(0.5, 0.85))
  tv_tol <- 0.01
  p <- 2L
  m10 <- two_block_m_convergence_for_pilot_start(rate, 10L, tv_tol, p)
  m100 <- two_block_m_convergence_for_pilot_start(rate, 100L, tv_tol, p)
  m_min <- two_block_m_convergence_for_pilot_start(rate, 100000L, tv_tol, p)
  expect_true(m100$m_convergence <= m10$m_convergence)
  expect_true(m100$m_convergence >= m_min$m_min)
  expect_equal(m10$m_min, m100$m_min)
})

test_that("two_block_optimize_pilot_cost returns sensible optimum", {
  rate <- .mock_two_block_rate(0.85, c(0.5, 0.85))
  opt <- two_block_optimize_pilot_cost(
    n                   = 100L,
    rate                = rate,
    tv_tol              = 0.01,
    m_convergence_pilot = 5L,
    p                   = 2L,
    n_pilot_max         = 500L
  )
  expect_true(opt$n_pilot_opt >= 1L)
  expect_true(opt$m_convergence_opt >= opt$cost_at_opt$m_min)
  expect_equal(opt$total_cost_opt, opt$cost_at_opt$total_cost)
  expect_true(all(opt$cost_curve$total_cost >= opt$total_cost_opt - 1e-12))
})

test_that("rGLMM stores pilot_cost_opt when verbose pilot path runs", {
  set.seed(1L)
  J <- 4L
  g <- factor(rep(seq_len(J), each = 5L))
  y <- stats::rpois(length(g), 2)
  re_nm <- "(Intercept)"
  x <- matrix(1, nrow = length(g), ncol = 1L, dimnames = list(NULL, re_nm))
  gl <- as.character(seq_len(J))
  x_hyper <- stats::setNames(
    list(matrix(1, J, 1L, dimnames = list(gl, re_nm))),
    re_nm
  )
  pf <- dNormal(mu = c(`(Intercept)` = 0), Sigma = matrix(1), dispersion = 1)
  out <- rGLMM_reg(
    n = 2L,
    y = y,
    x = x,
    block = g,
    x_hyper = x_hyper,
    prior_list = list(P = matrix(1)),
    pfamily_list = stats::setNames(list(pf), re_nm),
    family = poisson(),
    n_pilot = 4L,
    verbose = TRUE,
    progbar = FALSE
  )
  expect_equal(out$n_pilot, 4L)
  expect_equal(out$convergence$n_pilot_source, "explicit")
  expect_false(is.null(out$convergence$m_certificate))
})
