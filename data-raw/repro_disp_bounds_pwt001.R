# Smoke test: dGamma Block~1, pwt_measurement = 0.01 (weakest default-style prior).
lmebayes_root <- Sys.getenv("LMEBAYES_ROOT", unset = normalizePath("../lmebayes"))
setwd(lmebayes_root)
source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

N_DGAMMA <- 50L

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
dat  <- fx$dat
form <- fx$form

ps <- Prior_Setup_GLMM(
  form,
  data            = dat,
  pwt             = 0.01,
  pwt_measurement = 0.01
)
pf <- pfamily_list(ps)

m_disp <- ps$ing_prior_measurement
cat(sprintf(
  "\npwt_measurement = %.2f  n_prior = %.4g  n_combined = %.4g\n",
  ps$pwt_measurement, m_disp$n_prior, m_disp$n_combined
))
cat(sprintf(
  "ING sigma^2 window [%.1f, %.1f]  ratio = %.2f\n\n",
  m_disp$disp_lower, m_disp$disp_upper,
  m_disp$disp_upper / m_disp$disp_lower
))

disp_pf <- dGamma(
  shape          = m_disp$shape,
  rate           = m_disp$rate,
  beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
  Inv_Dispersion = TRUE,
  disp_lower     = m_disp$disp_lower,
  disp_upper     = m_disp$disp_upper
)

cat("=== lmerb n =", N_DGAMMA, "===\n\n")

t_fit <- system.time({
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pf,
    dispersion_ranef = disp_pf,
    n                = N_DGAMMA,
    progbar          = TRUE,
    verbose          = FALSE
  )
})

cat(sprintf("\nTiming: %.2f s\n", t_fit["elapsed"]))
cat(sprintf("draw_engine: %s\n", fit$convergence_info$draw_engine))
cat(sprintf(
  "sigma2: mean = %.4f, sd = %.4f, range = [%.1f, %.1f]\n",
  fit$group.dispersion.mean, stats::sd(fit$group.dispersion),
  min(fit$group.dispersion), max(fit$group.dispersion)
))

stopifnot(inherits(fit, "lmerb"))
stopifnot(identical(fit$convergence_info$draw_engine, "rGLMM_sweep_ing_block1_ind"))
stopifnot(is.finite(fit$group.dispersion.mean))
stopifnot(all(fit$group.dispersion >= m_disp$disp_lower))
stopifnot(all(fit$group.dispersion <= m_disp$disp_upper))

cat("\ntest_lmerb_dgamma_pwt001: OK\n")
