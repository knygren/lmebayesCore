## Report Block 1 / Block 2 candidates per stored draw by pwt_measurement.
lmebayes_root <- Sys.getenv("LMEBAYES_ROOT", unset = normalizePath("../lmebayes"))
setwd(lmebayes_root)
source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

N_DGAMMA <- 50L
fx <- .prepare_small5_lmerb_manual(n_schools = 5L)

for (pwt in c(0.01, 0.2, 0.49)) {
  ps <- Prior_Setup_lmebayes(
    fx$form,
    data            = fx$dat,
    pwt             = 0.01,
    pwt_measurement = pwt
  )
  pf <- pfamily_list(ps)
  m  <- ps$ing_prior_measurement
  disp_pf <- dGamma(
    shape          = m$shape,
    rate           = m$rate,
    beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
    Inv_Dispersion = TRUE,
    disp_lower     = m$disp_lower,
    disp_upper     = m$disp_upper
  )
  fit <- lmerb(
    fx$form,
    data             = fx$dat,
    pfamily_list     = pf,
    dispersion_ranef = disp_pf,
    n                = N_DGAMMA,
    progbar          = FALSE,
    verbose          = FALSE
  )
  cat(sprintf(
    "\npwt_measurement = %.2f  window [%.1f, %.1f]  m_convergence = %d\n",
    pwt, m$disp_lower, m$disp_upper, fit$m_convergence
  ))
  cat(sprintf(
    "  Block 1 ranef.iters.mean (candidates per stored draw): %.2f\n",
    fit$ranef.iters.mean
  ))
  if (!is.null(fit$fixef.iters.mean)) {
    cat("  Block 2 fixef.iters.mean (per RE component):\n")
    print(round(fit$fixef.iters.mean, 2))
  }
  cat(sprintf(
    "  mean(ranef.iters) raw total per chain: %.0f\n",
    mean(fit$ranef.iters)
  ))
}
