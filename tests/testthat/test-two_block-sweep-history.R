library(glmbayesCore)

.mini_v6_inputs <- function(family = gaussian(), inner_sweeps = 3L, n_chains = 2L) {
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
  fixef_mode <- stats::setNames(list(c(`(Intercept)` = 0)), re_nm)
  list(
    n_chains = n_chains,
    start_fixef = fixef_mode,
    inner_sweeps = inner_sweeps,
    design = list(
      y = y,
      Z = x,
      groups = g,
      X_hyper = x_hyper,
      re_coef_names = re_nm,
      group_name = "g"
    ),
    block1_prior = prior_list,
    pfamily_list = pfamily_list,
    re_names = re_nm,
    group_levels = gl,
    family = family,
    fixef_mode = fixef_mode
  )
}

test_that("run_sweep_outer_chains_v6 always returns sweep_history", {
  inp <- .mini_v6_inputs(inner_sweeps = 3L)
  out <- run_sweep_outer_chains_v6(
    n_chains = inp$n_chains,
    start_fixef = inp$start_fixef,
    inner_sweeps = inp$inner_sweeps,
    design = inp$design,
    block1_prior = inp$block1_prior,
    pfamily_list = inp$pfamily_list,
    family = inp$family,
    re_names = inp$re_names,
    group_levels = inp$group_levels,
    fixef_mode = inp$fixef_mode,
    b_mode = matrix(0, length(inp$group_levels), 1L),
    progbar = FALSE,
    stage_label = "main"
  )
  expect_s3_class(out$sweep_history, "two_block_sweep_history")
  expect_equal(out$sweep_history$n_sweeps, 3L)
  expect_equal(out$sweep_history$stage, "main")
  n_params <- length(inp$fixef_mode[[1L]])
  expect_equal(
    nrow(out$sweep_history$table),
    inp$inner_sweeps * n_params + n_params
  )
})

test_that("rGLMM verbose no longer auto-prints sweep history table", {
  inp <- .mini_v6_inputs(family = poisson(), inner_sweeps = 2L)
  out <- capture.output(
    fit <- rGLMM(
      n = 2L,
      y = inp$design$y,
      x = inp$design$Z,
      block = inp$design$groups,
      x_hyper = inp$design$X_hyper,
      prior_list = inp$block1_prior,
      pfamily_list = inp$pfamily_list,
      family = inp$family,
      n_pilot = 0L,
      verbose = TRUE,
      progbar = FALSE
    ),
    type = "output"
  )
  expect_s3_class(fit$sweep_history, "two_block_sweep_history")
  expect_false(any(grepl("fixef by sweep", out, fixed = TRUE)))
})

test_that("print.two_block_sweep_history respects max_sweeps", {
  inp <- .mini_v6_inputs(inner_sweeps = 5L)
  out <- run_sweep_outer_chains_v6(
    n_chains = inp$n_chains,
    start_fixef = inp$start_fixef,
    inner_sweeps = inp$inner_sweeps,
    design = inp$design,
    block1_prior = inp$block1_prior,
    pfamily_list = inp$pfamily_list,
    family = inp$family,
    re_names = inp$re_names,
    group_levels = inp$group_levels,
    fixef_mode = inp$fixef_mode,
    b_mode = matrix(0, length(inp$group_levels), 1L),
    progbar = FALSE,
    stage_label = "main"
  )
  full <- capture.output(print(out$sweep_history))
  short <- capture.output(print(out$sweep_history, max_sweeps = 2L))
  expect_lt(length(short), length(full))
  expect_false(any(grepl("sweep 1", short, fixed = TRUE)))
  expect_true(any(grepl("sweep 5", short, fixed = TRUE)))
})
