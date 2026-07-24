## Smoke test: plot_var_convergence()/plot_mean_convergence()/
## plot_sweep_history_diag() fit-object S3 methods (rLMMNormal_reg,
## rLMMindepNormalGamma_reg, rGLMM_reg, rlmerb, rglmerb).
devtools::load_all(".", quiet = TRUE)

set.seed(42)
J   <- 20L
n_j <- 20L
dat <- data.frame(
  group = factor(rep(paste0("g", seq_len(J)), each = n_j)),
  x1    = stats::rnorm(J * n_j)
)
b0_true <- stats::setNames(stats::rnorm(J, sd = 1.5), levels(dat$group))
b1_true <- stats::setNames(stats::rnorm(J, sd = 0.8), levels(dat$group))
dat$y <- 2 + b0_true[as.character(dat$group)] +
  (1.5 + b1_true[as.character(dat$group)]) * dat$x1 +
  stats::rnorm(nrow(dat), sd = 1)

form_lmer <- y ~ x1 + (1 + x1 || group)

design <- model_setup(form_lmer, data = dat)
grp <- design$groups
attr(grp, "group_name") <- design$group_name
re_names <- design$re_coef_names

ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.05)
pf_known_vcov <- pfamily_list(ps)
pf_est_vcov   <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
prior_list_disp <- list(dispersion = ps$dispersion_ranef)

out_dir <- tempdir()
png_dev <- function(name) grDevices::png(file.path(out_dir, name), width = 900, height = 700)

## --- rLMMNormal_reg_known_vcov(TWO_BLOCK_GIBBS): has sweep_history AND
## qualifies for the *exact* reference (dispersion fixed, vcov known). ------
fit1 <- rLMMNormal_reg_known_vcov(
  n = 30L, y = design$y, D = design$Z, group = grp, W = design$X_hyper,
  prior_list = prior_list_disp, pfamily_list = pf_known_vcov,
  progbar = FALSE, verbose = FALSE, sim_method = "TWO_BLOCK_GIBBS"
)
stopifnot(lmebayesCore:::.lmebayes_convergence_exact_ref_ok(fit1))
png_dev("s3_var_fit1.png");  plot_var_convergence(fit1);  grDevices::dev.off()
png_dev("s3_mean_fit1.png"); plot_mean_convergence(fit1); grDevices::dev.off()
png_dev("s3_diag_fit1.png")
plot_sweep_history_diag(fit1, coef_focus = list(c("(Intercept)", "(Intercept)"), c("x1", "(Intercept)")))
grDevices::dev.off()
cat("fit1 (rLMMNormal_reg_known_vcov, exact ref) OK\n\n")

## --- rLMMNormal_reg_known_vcov(DEFAULT): no sweep_history -> clear error --
fit1b <- rLMMNormal_reg_known_vcov(
  n = 30L, y = design$y, D = design$Z, group = grp, W = design$X_hyper,
  prior_list = prior_list_disp, pfamily_list = pf_known_vcov,
  progbar = FALSE, verbose = FALSE, sim_method = "DEFAULT"
)
err <- tryCatch(plot_var_convergence(fit1b), error = function(e) e)
stopifnot(inherits(err, "error"), grepl("sweep_history", conditionMessage(err)))
cat("fit1b (DEFAULT/iid, no sweep_history) -> clear error OK:\n  ", conditionMessage(err), "\n\n")

## --- rLMMNormal_reg_estimated_vcov(): sweep_history present but vcov
## estimated -> empirical fallback (no exact reference). --------------------
fit2 <- rLMMNormal_reg_estimated_vcov(
  n = 30L, y = design$y, D = design$Z, group = grp, W = design$X_hyper,
  prior_list = prior_list_disp, pfamily_list = pf_est_vcov,
  progbar = FALSE, verbose = FALSE
)
stopifnot(!lmebayesCore:::.lmebayes_convergence_exact_ref_ok(fit2))
png_dev("s3_var_fit2.png"); plot_var_convergence(fit2, n_chains = 5L); grDevices::dev.off()
cat("fit2 (rLMMNormal_reg_estimated_vcov, empirical fallback) OK\n\n")

## --- rlmerb(): same underlying engine as fit1, via the model_setup()/
## Prior_Setup_lmebayes() pipeline; result_class is lost (class == "rlmerb")
## so exact-ref detection must come from fit$Prior$dispersion_mode. ---------
prior_rlmerb <- list(
  Sigma_ranef  = as.matrix(ps$Sigma_ranef),
  prior_list   = ps$prior_list,
  pfamily_list = pf_known_vcov,
  any_non_normal = FALSE
)
fit3 <- rlmerb(
  n = 30L, design = design, prior = prior_rlmerb,
  dispersion_ranef = as.numeric(ps$dispersion_ranef),
  progbar = FALSE, verbose = FALSE, print_icm_table = FALSE,
  sim_method = "TWO_BLOCK_GIBBS"
)
stopifnot(identical(class(fit3), c("rlmerb", "list")))
stopifnot(lmebayesCore:::.lmebayes_convergence_exact_ref_ok(fit3))
png_dev("s3_var_fit3.png"); plot_var_convergence(fit3); grDevices::dev.off()
cat("fit3 (rlmerb(), exact ref via Prior$dispersion_mode) OK\n\n")

## --- rLMMindepNormalGamma_reg_known_vcov(): dispersion always estimated
## per group here -> never an exact reference, even though vcov is known. --
ps_g <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.05, dispformula = ~group)
disp_pf_list <- dGamma_list(ps_g)
group_levels <- levels(grp)
p_re <- length(re_names)
shape_group      <- stats::setNames(numeric(J), group_levels)
rate_group       <- stats::setNames(numeric(J), group_levels)
disp_lower_group <- stats::setNames(numeric(J), group_levels)
disp_upper_group <- stats::setNames(numeric(J), group_levels)
for (lev in group_levels) {
  pl <- disp_pf_list[[lev]]$prior_list
  shape_group[[lev]]      <- pl$shape[1L]
  rate_group[[lev]]       <- pl$rate[1L]
  disp_lower_group[[lev]] <- pl$disp_lower
  disp_upper_group[[lev]] <- pl$disp_upper
}
prior_list_disp_group <- list(
  mu               = matrix(0, nrow = p_re, ncol = 1L, dimnames = list(re_names, NULL)),
  Sigma            = as.matrix(ps_g$Sigma_ranef),
  shape_group      = shape_group,
  rate_group       = rate_group,
  disp_lower_group = disp_lower_group,
  disp_upper_group = disp_upper_group
)
fit4 <- rLMMindepNormalGamma_reg_known_vcov(
  n = 15L, y = design$y, D = design$Z, group = grp, W = design$X_hyper,
  prior_list = prior_list_disp_group, pfamily_list = pf_known_vcov,
  progbar = FALSE, verbose = FALSE
)
stopifnot(!lmebayesCore:::.lmebayes_convergence_exact_ref_ok(fit4))
png_dev("s3_var_fit4.png"); plot_var_convergence(fit4, n_chains = 5L); grDevices::dev.off()
cat("fit4 (rLMMindepNormalGamma_reg_known_vcov, never exact) OK\n\n")

## --- rGLMM_reg_known_vcov(poisson): non-Gaussian -> never an exact
## reference regardless of any_non_normal. ----------------------------------
dat_pois <- data.frame(
  group = factor(rep(paste0("g", seq_len(J)), each = n_j)),
  x1    = stats::rnorm(J * n_j, sd = 0.3)
)
b0_pois <- stats::setNames(stats::rnorm(J, sd = 0.3), levels(dat_pois$group))
b1_pois <- stats::setNames(stats::rnorm(J, sd = 0.15), levels(dat_pois$group))
eta_pois <- 1.5 + b0_pois[as.character(dat_pois$group)] +
  (0.3 + b1_pois[as.character(dat_pois$group)]) * dat_pois$x1
dat_pois$y <- stats::rpois(nrow(dat_pois), lambda = exp(eta_pois))

design_pois <- model_setup(form_lmer, data = dat_pois)
grp_pois <- design_pois$groups
attr(grp_pois, "group_name") <- design_pois$group_name
ps_pois <- Prior_Setup_lmebayes(form_lmer, data = dat_pois, family = poisson(), pwt = 0.05)
pf_pois <- pfamily_list(ps_pois)
fit5 <- rGLMM_reg_known_vcov(
  n = 15L, y = design_pois$y, D = design_pois$Z, group = grp_pois, W = design_pois$X_hyper,
  prior_list = list(), pfamily_list = pf_pois, family = poisson(),
  progbar = FALSE, verbose = FALSE
)
stopifnot(!lmebayesCore:::.lmebayes_convergence_exact_ref_ok(fit5))
png_dev("s3_var_fit5.png"); plot_var_convergence(fit5, n_chains = 5L); grDevices::dev.off()
cat("fit5 (rGLMM_reg_known_vcov(poisson), never exact) OK\n\n")

## --- rglmerb(family = poisson()): same non-Gaussian never-exact rule,
## through the rlmerb/rglmerb-style Prior$dispersion_mode path ("none"). ----
prior_rglmerb <- list(
  Sigma_ranef    = as.matrix(ps_pois$Sigma_ranef),
  prior_list     = ps_pois$prior_list,
  pfamily_list   = pf_pois,
  any_non_normal = FALSE
)
fit6 <- rglmerb(
  n = 15L, design = design_pois, prior = prior_rglmerb, family = poisson(),
  progbar = FALSE, verbose = FALSE
)
stopifnot(identical(class(fit6), c("rglmerb", "list")))
stopifnot(!lmebayesCore:::.lmebayes_convergence_exact_ref_ok(fit6))
png_dev("s3_var_fit6.png"); plot_var_convergence(fit6, n_chains = 5L); grDevices::dev.off()
cat("fit6 (rglmerb(poisson), never exact) OK\n\n")

cat("ALL OK\n")
