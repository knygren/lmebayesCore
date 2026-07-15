## Quick check: lambda* for 5-school dGamma fixture at several pwt_measurement values.
## Run: Rscript data-raw/_diag_lambda_star.R

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    pkgload::load_all(".", quiet = TRUE)
  }
})

if (!requireNamespace("lmebayes", quietly = TRUE)) {
  stop("Need lmebayes sibling package.", call. = FALSE)
}

source("C:/Rpackages/lmebayes/tests/manual/_small5_lmerb_fixture.R")
fx <- .prepare_small5_lmerb_manual(5L)
dat <- fx$dat
form <- fx$form

ms <- lmebayes::model_setup(form, data = dat)
block <- ms$groups
group_levels <- levels(block)

for (pwt_meas in c(0.05, 0.1, 0.49)) {
  cat("\n--- pwt_measurement =", pwt_meas, "---\n")
  ps <- lmebayes::Prior_Setup_lmebayes(
    form, data = dat, pwt = 0.05, pwt_measurement = pwt_meas,
    max_disp_perc = 0.8, dispformula = ~school_id
  )
  pf <- lmebayes::pfamily_list(ps)
  disp_pf <- lmebayes::dGamma_list(ps, max_disp_perc = 0.8, warn_asymmetric = FALSE)
  ing <- ps$ing_prior_measurement_group
  P <- solve(ps$Sigma_ranef)

  ing_prior_list <- list(
    shape_group = vapply(ing, function(g) g$shape_ING, 0),
    rate_group  = vapply(ing, function(g) g$rate_gamma, 0),
    disp_lower_group = vapply(disp_pf, function(p) p$prior_list$disp_lower, 0),
    disp_upper_group = vapply(disp_pf, function(p) p$prior_list$disp_upper, 0)
  )
  names(ing_prior_list$shape_group) <- group_levels
  names(ing_prior_list$rate_group) <- group_levels
  names(ing_prior_list$disp_lower_group) <- group_levels
  names(ing_prior_list$disp_upper_group) <- group_levels

  rate_inputs <- glmbayesCore:::.rLMM_measurement_rate_inputs(
    ing_prior_list, block, group_levels, "diag"
  )
  pl1_rate <- glmbayesCore:::.rLMM_block1_prior_gaussian(P, rate_inputs$dispersion_scalar)

  rate <- tryCatch(
    glmbayesCore::two_block_rate_from_pfamily_list(
      x = ms$Z, block = block, x_hyper = ms$X_hyper,
      prior_list_block1 = pl1_rate, pfamily_list = pf,
      weights = rate_inputs$weights, family = gaussian(),
      group_levels = group_levels
    ),
    error = function(e) e
  )
  if (inherits(rate, "error")) {
    cat("STOPPED:", conditionMessage(rate), "\n")
  } else {
    cat("lambda_star =", format(rate$lambda_star, digits = 17), "\n")
  }
}
