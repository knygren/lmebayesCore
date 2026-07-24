## sim_method dispatch for the fixed-dispersion / known-variance-components
## ("lmm_fixed_known") route: rLMMNormal_reg_known_vcov() dispatches between
## the exact-iid engine (rLMMNormal_reg_known_vcov_iid() / rLMMNormal_joint_iid(),
## sim_method = "DEFAULT") and the two-block Gibbs engine
## (rLMMNormal_reg_known_vcov_two_bg(), sim_method = "TWO_BLOCK_GIBBS"). See
## inst/README_KNOWN_VCOV_GAUSSIAN.md \S4.4/\S4.4a.

## Single random-intercept fixture: J groups, fixed per-group dispersion,
## dNormal() Block~2 prior (known tau^2) -- exactly the lmm_fixed_known cell.
.sim_method_fixture <- function() {
  set.seed(20260720)
  J   <- 5L
  n_j <- 8L
  n_obs <- J * n_j
  group_levels <- paste0("g", seq_len(J))
  group <- factor(rep(group_levels, each = n_j), levels = group_levels)
  attr(group, "group_name") <- "group"

  tau2_true   <- 0.7
  sigma2_true <- stats::setNames(c(0.3, 0.6, 0.9, 1.2, 1.5), group_levels)
  b_true <- stats::setNames(stats::rnorm(J, sd = sqrt(tau2_true)), group_levels)
  y <- 2 + b_true[as.character(group)] +
    stats::rnorm(n_obs, sd = sqrt(sigma2_true[as.character(group)]))

  D <- matrix(1, n_obs, 1, dimnames = list(NULL, "(Intercept)"))
  W <- list(
    "(Intercept)" = matrix(1, J, 1, dimnames = list(NULL, "(Intercept)"))
  )
  prior_list <- list(dispersion = unname(sigma2_true))
  pfamily_list <- list(
    "(Intercept)" = glmbayesCore::dNormal(
      mu = 0, Sigma = matrix(100), dispersion = tau2_true
    )
  )

  list(
    n_obs = n_obs, J = J, group = group, y = y,
    D = D, W = W,
    prior_list = prior_list, pfamily_list = pfamily_list,
    sigma2_true = sigma2_true, tau2_true = tau2_true
  )
}

test_that("rLMMNormal_reg_known_vcov() dispatches sim_method and tags sim_method_used", {
  fx <- .sim_method_fixture()

  set.seed(1)
  fit_default <- rLMMNormal_reg_known_vcov(
    n = 20L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE
  )
  expect_identical(fit_default$sim_method_used, "DEFAULT")
  expect_identical(fit_default$m_convergence, 1L)
  expect_identical(fit_default$convergence_info$method, "exact_iid")
  expect_identical(fit_default$draw_engine, "rLMMNormal_joint_iid")

  set.seed(1)
  fit_iid <- rLMMNormal_reg_known_vcov(
    n = 20L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE, sim_method = "DEFAULT"
  )
  expect_identical(fit_iid$sim_method_used, "DEFAULT")
  expect_equal(fit_iid$fixef, fit_default$fixef)

  set.seed(1)
  fit_gibbs <- rLMMNormal_reg_known_vcov(
    n = 20L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE, sim_method = "TWO_BLOCK_GIBBS"
  )
  expect_identical(fit_gibbs$sim_method_used, "TWO_BLOCK_GIBBS")
  expect_gt(fit_gibbs$m_convergence, 1L)
  expect_false(identical(fit_gibbs$draw_engine, "rLMMNormal_joint_iid"))

  ## Direct calls to the two named engines behind the dispatcher.
  set.seed(1)
  fit_iid_direct <- rLMMNormal_reg_known_vcov_iid(
    n = 20L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE
  )
  expect_identical(fit_iid_direct$sim_method_used, "DEFAULT")
  expect_equal(fit_iid_direct$fixef, fit_default$fixef)

  set.seed(1)
  fit_bg_direct <- rLMMNormal_reg_known_vcov_two_bg(
    n = 20L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE
  )
  expect_identical(fit_bg_direct$sim_method_used, "TWO_BLOCK_GIBBS")
  expect_equal(fit_bg_direct$fixef, fit_gibbs$fixef)

  for (cls in c("rLMMNormal_reg_known_vcov", "rLMMNormal_reg", "list")) {
    expect_true(cls %in% class(fit_default))
    expect_true(cls %in% class(fit_gibbs))
  }
})

test_that("rLMMNormal_joint_iid() fixef_mean matches lmerb_posterior_mean() exactly", {
  fx <- .sim_method_fixture()

  fit_iid <- rLMMNormal_joint_iid(
    n = 50L, y = fx$y, x = fx$D, group = fx$group, x_hyper = fx$W,
    prior_list_block1 = list(
      dispersion = unname(fx$sigma2_true), ddef = FALSE
    ),
    pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE
  )

  design <- list(
    y = fx$y, Z = fx$D, groups = fx$group, X_hyper = fx$W,
    re_coef_names = "(Intercept)", group_name = "group"
  )
  ## Sigma_ranef must match what rLMMNormal_joint_iid() now derives
  ## internally from pfamily_list (diag(tau2_k), here tau2_true), not an
  ## independently-supplied P/Sigma.
  mpl <- list(
    Sigma_ranef      = matrix(fx$tau2_true, 1, 1),
    dispersion_ranef = unname(fx$sigma2_true),
    prior_list       = list(
      "(Intercept)" = list(
        mu_fixef = 0, Sigma_fixef = matrix(100),
        dispersion_fixef = fx$tau2_true
      )
    )
  )
  pm <- lmerb_posterior_mean(design, mpl)

  expect_equal(
    as.numeric(fit_iid$fixef_mean[["(Intercept)"]]),
    as.numeric(pm$fixef[["(Intercept)"]])
  )
  expect_identical(fit_iid$m_convergence, 1L)
  expect_identical(class(fit_iid)[1L], "rLMMNormal_joint_iid")
})

test_that("iid and two-block Gibbs engines agree on the posterior mean (Monte Carlo)", {
  fx <- .sim_method_fixture()

  set.seed(2026)
  fit_iid <- rLMMNormal_reg_known_vcov(
    n = 2000L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE, sim_method = "DEFAULT"
  )
  set.seed(2026)
  fit_gibbs <- rLMMNormal_reg_known_vcov(
    n = 2000L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE, sim_method = "TWO_BLOCK_GIBBS"
  )

  mean_iid   <- mean(fit_iid$fixef[["(Intercept)"]][, "(Intercept)"])
  mean_gibbs <- mean(fit_gibbs$fixef[["(Intercept)"]][, "(Intercept)"])
  expect_equal(mean_iid, mean_gibbs, tolerance = 0.1)

  ## Every draw should differ from the previous one for both engines (no
  ## degenerate/constant chain).
  expect_gt(
    stats::sd(fit_iid$fixef[["(Intercept)"]][, "(Intercept)"]), 0
  )
  expect_gt(
    stats::sd(fit_gibbs$fixef[["(Intercept)"]][, "(Intercept)"]), 0
  )
})

test_that("sim_method validation rejects unknown values", {
  fx <- .sim_method_fixture()

  expect_error(
    rLMMNormal_reg_known_vcov(
      n = 5L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
      prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
      progbar = FALSE, sim_method = "bogus"
    ),
    "sim_method"
  )

  expect_error(
    rlmerb(
      n = 5L, design = list(), prior = list(), dispersion_ranef = 1,
      sim_method = "bogus"
    ),
    "sim_method"
  )

  expect_error(
    rglmerb(
      n = 5L, design = list(), prior = list(), family = gaussian(),
      dispersion_ranef = 1, sim_method = "bogus"
    ),
    "sim_method"
  )
})

test_that("rLMMNormal_joint_iid()'s b_mean is the exact posterior mean of beta_j (not b_last)", {
  fx <- .sim_method_fixture()

  out <- rLMMNormal_joint_iid(
    n = 2000L, y = fx$y, x = fx$D, group = fx$group, x_hyper = fx$W,
    prior_list_block1 = list(dispersion = unname(fx$sigma2_true), ddef = FALSE),
    pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE
  )

  ## Manual closed-form E[beta_j | data] using fx's known sigma2_j/tau2 and
  ## rLMMNormal_joint_iid()'s own exact gamma_mean (out$fixef_mean) -- the
  ## same "plug fixef_mean into the conditional-mean map" identity b_mean
  ## itself relies on, computed independently here from first principles.
  group_levels <- levels(fx$group)
  n_j    <- as.integer(table(fx$group)[group_levels])
  ybar_j <- tapply(fx$y, fx$group, mean)[group_levels]
  gamma_hat <- out$fixef_mean[["(Intercept)"]][["(Intercept)"]]
  prec_prior_bj <- 1 / fx$tau2_true
  prec_lik_j    <- n_j / fx$sigma2_true[group_levels]
  exact_beta_mean <- (prec_prior_bj * gamma_hat + prec_lik_j * ybar_j) /
    (prec_prior_bj + prec_lik_j)

  expect_equal(
    unname(out$b_mean[group_levels, "(Intercept)"]),
    as.numeric(exact_beta_mean),
    tolerance = 1e-8
  )

  ## b_mean must be a stable point estimate: it should be far closer to the
  ## across-draw average of fit$coefficients than the single last draw
  ## (b_last) is, on average across groups.
  mcmc_mean_j <- tapply(
    out$coefficients[["(Intercept)"]], out$coefficients$group, mean
  )[group_levels]
  err_b_mean <- mean(abs(out$b_mean[group_levels, "(Intercept)"] - mcmc_mean_j))
  err_b_last <- mean(abs(out$b_last[group_levels, "(Intercept)"] - mcmc_mean_j))
  expect_lt(err_b_mean, err_b_last)
  expect_lt(err_b_mean, 0.05)
})

test_that("rLMMNormal_reg_known_vcov()'s ranef.mode is a stable point estimate for both engines", {
  fx <- .sim_method_fixture()

  set.seed(2026)
  fit_iid <- rLMMNormal_reg_known_vcov(
    n = 3000L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE, sim_method = "DEFAULT"
  )
  set.seed(2026)
  fit_gibbs <- rLMMNormal_reg_known_vcov(
    n = 3000L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
    prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
    progbar = FALSE, verbose = FALSE, sim_method = "TWO_BLOCK_GIBBS"
  )

  ## The two engines target the same exact posterior, so their ranef.mode
  ## (iid: exact b_mean; two-block Gibbs: ICM mode, which coincides with the
  ## exact posterior mean for this fully Gaussian-conjugate model) should
  ## agree closely -- unlike a raw last-draw comparison, which would not.
  expect_equal(
    fit_iid$ranef.mode[, "(Intercept)"],
    fit_gibbs$ranef.mode[, "(Intercept)"],
    tolerance = 1e-6
  )

  ## And each engine's own ranef.mode should track its own MCMC mean.
  for (fit in list(fit_iid, fit_gibbs)) {
    mcmc_mean_j <- tapply(
      fit$coefficients[["(Intercept)"]], fit$coefficients$group, mean
    )[rownames(fit$ranef.mode)]
    expect_equal(
      unname(fit$ranef.mode[, "(Intercept)"]),
      as.numeric(mcmc_mean_j),
      tolerance = 0.05
    )
  }
})

test_that("sim_method is accepted-but-inert on routes with only a two-block Gibbs engine", {
  ## rLMMNormal_reg_estimated_vcov() (estimated variance components -- always
  ## requires at least one non-dNormal Block~2 component, so it can never hit
  ## the lmm_fixed_known route) accepts sim_method for interface parity with
  ## rLMMNormal_reg_known_vcov(), still validates it, but never dispatches on
  ## it (there is no iid engine for estimated variance components).
  fn_args <- formals(rLMMNormal_reg_estimated_vcov)
  expect_true("sim_method" %in% names(fn_args))
  expect_identical(eval(fn_args$sim_method), "DEFAULT")

  fx <- .sim_method_fixture()
  expect_error(
    rLMMNormal_reg_estimated_vcov(
      n = 5L, y = fx$y, D = fx$D, group = fx$group, W = fx$W,
      prior_list = fx$prior_list, pfamily_list = fx$pfamily_list,
      progbar = FALSE, sim_method = "bogus"
    ),
    "sim_method"
  )
})
