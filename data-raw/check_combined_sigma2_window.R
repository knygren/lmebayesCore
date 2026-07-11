## Compare sigma^2 truncation windows: prior-only vs combined (n_prior + n).
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
root <- Sys.getenv("LMEBAYES_ROOT", unset = normalizePath("../lmebayes"))
source(file.path(root, "tests/manual/_small5_lmerb_fixture.R"))

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
n <- length(fx$design$y)
p_re <- length(fx$design$re_coef_names)
sigma2_hat <- fx$design$residual_var

prior_only_window <- function(n_prior, d_hat, n_obs, p_re) {
  shape <- (n_prior + 1) / 2 + p_re / 2
  rate  <- d_hat * (n_prior + p_re - 1) / 2
  glmbayesCore:::.lmebayes_ing_prior_quantile_window(shape, rate)
}

cat("sigma2_hat =", sigma2_hat, " n =", n, " p_re =", p_re, "\n\n")
for (pwt in c(0.2, 0.49)) {
  ps <- Prior_Setup_lmebayes(
    fx$form, data = fx$dat, pwt = 0.01, pwt_measurement = pwt
  )
  m <- ps$ing_prior_measurement
  old <- prior_only_window(m$n_prior, sigma2_hat, n, p_re)
  cat(sprintf(
    "pwt=%.2f  n_prior=%.1f  n_combined=%.1f\n",
    pwt, m$n_prior, m$n_combined
  ))
  cat(sprintf(
    "  OLD prior-only : [%.1f, %.1f]  ratio=%.2f\n",
    old$disp_lower, old$disp_upper, old$disp_upper / old$disp_lower
  ))
  cat(sprintf(
    "  NEW combined   : [%.1f, %.1f]  ratio=%.2f\n\n",
    m$disp_lower, m$disp_upper, m$disp_upper / m$disp_lower
  ))
}
