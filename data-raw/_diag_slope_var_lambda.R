## Slope-variance vs lambda* — scratch diagnostic (data-raw, not shipped).
suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    pkgload::load_all(".", quiet = TRUE)
  }
  library(lmebayes)
})

source("C:/Rpackages/lmebayes/tests/manual/_small5_lmerb_fixture.R")
fx <- .prepare_small5_lmerb_manual(5L)
dat <- fx$dat
form <- fx$form
ms <- model_setup(form, data = dat)
block <- ms$groups
group_levels <- levels(block)

raw_lambda <- function(ps) {
  pf <- pfamily_list(ps)
  disp_pf <- dGamma_list(ps, max_disp_perc = 0.8, warn_asymmetric = FALSE)
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
  re_names <- names(ms$X_hyper)
  pfamily_list <- glmbayesCore:::.two_block_validate_pfamily_list(pf, re_names)
  prior_list_block2 <- lapply(pfamily_list, function(pf_i) {
    pl <- pf_i$prior_list
    list(
      mu = pl$mu,
      Sigma = pl$Sigma,
      dispersion = if (identical(pf_i$pfamily, "dNormal")) pl$dispersion else pl$disp_lower
    )
  })
  inp <- glmbayesCore:::.two_block_rate_inputs(
    ms$Z, block, ms$X_hyper, pl1_rate, prior_list_block2,
    rate_inputs$weights, gaussian(), group_levels
  )
  sp <- glmbayesCore:::.two_block_S_P11(inp)
  R <- chol(sp$P11)
  q <- nrow(R)
  Rinv <- backsolve(R, diag(q))
  M <- t(Rinv) %*% sp$S %*% Rinv
  M <- 0.5 * (M + t(M))
  max(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
}

ps_tmb <- Prior_Setup_lmebayes(
  form, data = dat, pwt = 0.05, pwt_measurement = 0.1,
  max_disp_perc = 0.8, dispformula = ~school_id
)
vc_lmer <- 4.7937^2  # lme4 REML on this fixture (see manual reference print)

cat("glmmTMB slope var:", ps_tmb$Sigma_ranef[2, 2],
    " lambda* =", raw_lambda(ps_tmb), "\n")

ps_cf <- ps_tmb
ps_cf$Sigma_ranef[2, 2] <- vc_lmer
ps_cf$sd_tau[2] <- sqrt(vc_lmer)
cat("lme4 slope var:", vc_lmer,
    " lambda* =", raw_lambda(ps_cf), "\n")

ps_pool <- Prior_Setup_lmebayes(
  form, data = dat, pwt = 0.05, pwt_measurement = 0.1,
  max_disp_perc = 0.8, dispformula = ~1
)
cat("pooled lme4 slope var:", ps_pool$Sigma_ranef[2, 2],
    " lambda* =", raw_lambda(ps_pool), "\n")
