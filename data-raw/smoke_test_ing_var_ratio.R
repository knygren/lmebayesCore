## Fast, fully-synthetic smoke test of plot_sweep_history_var_ratio() wired
## against the REAL rLMMindepNormalGamma_reg_known_vcov()/_estimated_vcov()
## sweeps-outer/chains-inner engines (same code path as Ex_13/Ex_14), bypassing
## Prior_Setup_lmebayes()/dGamma_list() calibration (and bayesrules::big_word_club)
## entirely so this runs in seconds rather than the demos' 10-20+ minutes.
devtools::load_all(".", quiet = TRUE)
grDevices::pdf(tempfile(fileext = ".pdf"))

set.seed(42)
J   <- 6L
n_j <- 10L
n_obs <- J * n_j
group_levels <- paste0("g", seq_len(J))
group <- factor(rep(group_levels, each = n_j), levels = group_levels)
attr(group, "group_name") <- "group"

tau2_true   <- 0.5
sigma2_true <- stats::setNames(rep(1, J), group_levels)
b_true <- stats::setNames(stats::rnorm(J, sd = sqrt(tau2_true)), group_levels)
y <- 2 + b_true[as.character(group)] +
  stats::rnorm(n_obs, sd = sqrt(sigma2_true[as.character(group)]))

D <- matrix(1, n_obs, 1, dimnames = list(NULL, "(Intercept)"))
W <- list("(Intercept)" = matrix(1, J, 1, dimnames = list(NULL, "(Intercept)")))
re_names <- "(Intercept)"
p_re <- 1L

## Per-group dGamma() shape/rate: well-behaved, diffuse prior centered near
## the true sigma^2_j = 1 (shape = 3 => finite variance, rate = shape - 1 => mean 1).
shape_group      <- stats::setNames(rep(3, J), group_levels)
rate_group       <- stats::setNames(rep(2, J), group_levels)
disp_lower_group <- stats::setNames(rep(0.05, J), group_levels)
disp_upper_group <- stats::setNames(rep(20, J), group_levels)

prior_list <- list(
  mu               = matrix(0, nrow = p_re, ncol = 1L, dimnames = list(re_names, NULL)),
  Sigma            = matrix(tau2_true, p_re, p_re, dimnames = list(re_names, re_names)),
  shape_group      = shape_group,
  rate_group       = rate_group,
  disp_lower_group = disp_lower_group,
  disp_upper_group = disp_upper_group
)

## --- Case A: known vcov (all-dNormal Block~2) --------------------------
pf_known <- list(
  "(Intercept)" = glmbayesCore::dNormal(
    mu = 0, Sigma = matrix(100), dispersion = tau2_true
  )
)

cat("\n### rLMMindepNormalGamma_reg_known_vcov() ###\n")
fit_known <- rLMMindepNormalGamma_reg_known_vcov(
  n = 200L, y = y, D = D, group = group, W = W,
  prior_list = prior_list, pfamily_list = pf_known,
  progbar = FALSE, verbose = FALSE
)
cat("fit_known class:", paste(class(fit_known), collapse = ", "), "\n")
cat("has sweep_history:", !is.null(fit_known$sweep_history), "\n")
cat("has pilot sweep_history:", !is.null(fit_known$pilot$sweep_history), "\n")
cat("has cov_by_sweep (main):", !is.null(fit_known$sweep_history$cov_by_sweep), "\n")

cat("\nplot_sweep_history_var_ratio(main, whitened = FALSE)...\n")
plot_sweep_history_var_ratio(fit_known$sweep_history, whitened = FALSE)
cat("plot_sweep_history_var_ratio(main, whitened = TRUE)...\n")
plot_sweep_history_var_ratio(fit_known$sweep_history, whitened = TRUE)
if (!is.null(fit_known$pilot$sweep_history)) {
  cat("plot_sweep_history_var_ratio(pilot, whitened = TRUE)...\n")
  plot_sweep_history_var_ratio(fit_known$pilot$sweep_history, whitened = TRUE)
}
cat("Case A OK\n")

## --- Case B: estimated vcov (at least one non-dNormal Block~2) --------
pf_est <- list(
  "(Intercept)" = glmbayesCore::dIndependent_Normal_Gamma(
    mu = 0, Sigma = matrix(100), shape = 3, rate = 2 * tau2_true
  )
)

cat("\n### rLMMindepNormalGamma_reg_estimated_vcov() ###\n")
fit_est <- rLMMindepNormalGamma_reg_estimated_vcov(
  n = 200L, y = y, D = D, group = group, W = W,
  prior_list = prior_list, pfamily_list = pf_est,
  progbar = FALSE, verbose = FALSE
)
cat("fit_est class:", paste(class(fit_est), collapse = ", "), "\n")
cat("has sweep_history:", !is.null(fit_est$sweep_history), "\n")
cat("has pilot sweep_history:", !is.null(fit_est$pilot$sweep_history), "\n")
cat("has cov_by_sweep (main):", !is.null(fit_est$sweep_history$cov_by_sweep), "\n")

cat("\nplot_sweep_history_var_ratio(main, whitened = FALSE)...\n")
plot_sweep_history_var_ratio(fit_est$sweep_history, whitened = FALSE)
cat("plot_sweep_history_var_ratio(main, whitened = TRUE)...\n")
plot_sweep_history_var_ratio(fit_est$sweep_history, whitened = TRUE)
if (!is.null(fit_est$pilot$sweep_history)) {
  cat("plot_sweep_history_var_ratio(pilot, whitened = TRUE)...\n")
  plot_sweep_history_var_ratio(fit_est$pilot$sweep_history, whitened = TRUE)
}
cat("Case B OK\n")

grDevices::dev.off()
cat("\nALL OK\n")
