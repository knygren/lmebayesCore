## Smoke test: plot_var_convergence()/plot_mean_convergence()/
## plot_sweep_history_diag() fit-object S3 methods (rLMMNormal_reg,
## rLMMindepNormalGamma_reg, rGLMM_reg, rlmerb, rglmerb), plus the
## split/max_whitened per-group-plot behaviour (each split group is its own
## chart/page, never stacked via layout()/mfrow across groups) and its
## '...' forwarding from the fit-object methods down to the .default
## methods via .lmebayes_convergence_from_fit().
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
grp <- design$group
attr(grp, "group_name") <- design$group_name
re_names <- design$re_coef_names

ps <- Prior_Setup_GLMM(form_lmer, data = dat, pwt = 0.05)
pf_known_vcov <- pfamily_list(ps)
pf_est_vcov   <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
prior_list_disp <- list(dispersion = ps$dispersion_ranef)

out_dir <- tempdir()
png_dev <- function(name) grDevices::png(file.path(out_dir, name), width = 900, height = 700)

## --- rLMMNormal_reg_known_vcov(TWO_BLOCK_GIBBS): has sweep_history AND
## qualifies for the *exact* reference (dispersion fixed, vcov known). ------
fit1 <- rLMMNormal_reg_known_vcov(
  n = 30L, y = design$y, D = design$D, group = grp, W = design$W,
  pfamily_list = pf_known_vcov, dispprior_list = prior_list_disp,
  progbar = FALSE, verbose = FALSE, sim_method = "TWO_BLOCK_GIBBS"
)
stopifnot(lmebayesCore:::.lmebayes_convergence_exact_ref_ok(fit1))
png_dev("s3_var_fit1.png");  plot_var_convergence(fit1);  grDevices::dev.off()
png_dev("s3_mean_fit1.png"); plot_mean_convergence(fit1); grDevices::dev.off()
png_dev("s3_diag_fit1.png")
plot_sweep_history_diag(fit1, coef_focus = list(c("(Intercept)", "(Intercept)"), c("x1", "(Intercept)")))
grDevices::dev.off()
cat("fit1 (rLMMNormal_reg_known_vcov, exact ref) OK\n\n")

## --- split = "auto" (default) vs split = "none": fit1 has one
## Intercept-side and one slope-side Block-2 fixed effect, so the default
## split produces two separate chart pages (via .lmebayes_convergence_
## from_fit()'s '...' forwarding down to plot_var_convergence.default()),
## while split = "none" restores the single-combined-chart page. Uses a
## multi-page pdf() device (unlike png()) so both pages are actually kept.
pdf(file.path(out_dir, "s3_var_fit1_split_auto.pdf"))
plot_var_convergence(fit1)                  # split = "auto" (default) -> 2 pages
grDevices::dev.off()
pdf(file.path(out_dir, "s3_var_fit1_split_none.pdf"))
plot_var_convergence(fit1, split = "none")  # single combined chart -> 1 page
grDevices::dev.off()
pdf(file.path(out_dir, "s3_mean_fit1_split_auto.pdf"))
plot_mean_convergence(fit1)
grDevices::dev.off()
cat("fit1 split = 'auto' vs 'none' (via '...' forwarding on the S3 method) OK\n\n")

## --- .convergence_split_coef_focus(): grouping logic on synthetic keys ---
split_coef_focus <- lmebayesCore:::.convergence_split_coef_focus
keys_mixed <- data.frame(
  re_component = c("(Intercept)", "(Intercept)", "x1", "x2"),
  covariate    = c("(Intercept)", "z1", "(Intercept)", "(Intercept)"),
  stringsAsFactors = FALSE
)
g_mixed <- split_coef_focus(keys_mixed)
stopifnot(
  identical(names(g_mixed), c("Intercept predictors", "Slope predictors")),
  length(g_mixed[["Intercept predictors"]]) == 2L,
  length(g_mixed[["Slope predictors"]]) == 2L
)
keys_icpt_only <- keys_mixed[keys_mixed$re_component == "(Intercept)", , drop = FALSE]
g_icpt_only <- split_coef_focus(keys_icpt_only)
stopifnot(identical(names(g_icpt_only), ""), length(g_icpt_only[[1L]]) == 2L)
cat(".convergence_split_coef_focus() grouping OK\n\n")

## --- .convergence_split_whitened(): batches of at most max_per_plot -----
split_whitened <- lmebayesCore:::.convergence_split_whitened
series7 <- stats::setNames(as.list(seq_len(7L)), paste0("var", seq_len(7L)))
g7 <- split_whitened(series7, 4L)
stopifnot(
  identical(names(g7), c("var1-var4", "var5-var7")),
  length(g7[["var1-var4"]]) == 4L,
  length(g7[["var5-var7"]]) == 3L
)
series3 <- stats::setNames(as.list(seq_len(3L)), paste0("var", seq_len(3L)))
g3 <- split_whitened(series3, 4L)
stopifnot(identical(names(g3), ""), length(g3[[1L]]) == 3L)
cat(".convergence_split_whitened() batching OK\n\n")

## --- rLMMNormal_reg_known_vcov(DEFAULT): no sweep_history -> clear error --
fit1b <- rLMMNormal_reg_known_vcov(
  n = 30L, y = design$y, D = design$D, group = grp, W = design$W,
  pfamily_list = pf_known_vcov, dispprior_list = prior_list_disp,
  progbar = FALSE, verbose = FALSE, sim_method = "DEFAULT"
)
err <- tryCatch(plot_var_convergence(fit1b), error = function(e) e)
stopifnot(inherits(err, "error"), grepl("sweep_history", conditionMessage(err)))
cat("fit1b (DEFAULT/iid, no sweep_history) -> clear error OK:\n  ", conditionMessage(err), "\n\n")

## --- rLMMNormal_reg_estimated_vcov(): sweep_history present but vcov
## estimated -> empirical fallback (no exact reference). --------------------
fit2 <- rLMMNormal_reg_estimated_vcov(
  n = 30L, y = design$y, D = design$D, group = grp, W = design$W,
  pfamily_list = pf_est_vcov, dispprior_list = prior_list_disp,
  progbar = FALSE, verbose = FALSE
)
stopifnot(!lmebayesCore:::.lmebayes_convergence_exact_ref_ok(fit2))
png_dev("s3_var_fit2.png"); plot_var_convergence(fit2, n_chains = 5L); grDevices::dev.off()
cat("fit2 (rLMMNormal_reg_estimated_vcov, empirical fallback) OK\n\n")

## --- whitened = TRUE end-to-end, forcing a split via max_whitened (P = 2
## Block-2 fixed effects here, so max_whitened = 1 forces exactly 2 pages,
## each its own separate chart -- never one page with two stacked panels). -
pdf(file.path(out_dir, "s3_var_fit2_whitened_split.pdf"))
plot_var_convergence(fit2, whitened = TRUE, n_chains = 5L, max_whitened = 1L)
grDevices::dev.off()
cat("fit2 whitened = TRUE, max_whitened = 1 (forced split) OK\n\n")

## --- component = "precision": fit2 (estimated vcov, 2 ING components) has
## a populated disp_table; fit1 (known vcov, TWO_BLOCK_GIBBS) has none. -----
disp2 <- fit2$sweep_history$disp_table
stopifnot(
  is.data.frame(disp2), nrow(disp2) > 0L,
  setequal(unique(disp2$re_component), re_names),
  identical(unique(disp2$covariate), "precision"),
  all(disp2$mean > 0), all(disp2$sd[disp2$sweep > 0L] >= 0)
)
cat("fit2 sweep_history$disp_table populated (", nrow(disp2), "rows) OK\n\n")

disp1 <- fit1$sweep_history$disp_table
stopifnot(is.data.frame(disp1), nrow(disp1) == 0L)
cat("fit1 sweep_history$disp_table empty (known vcov) OK\n\n")

png_dev("s3_var_fit2_precision.png")
plot_var_convergence(fit2, component = "precision", n_chains = 5L)
grDevices::dev.off()
png_dev("s3_mean_fit2_precision.png")
plot_mean_convergence(fit2, component = "precision", n_chains = 5L)
grDevices::dev.off()
cat("fit2 component = 'precision' (named, non-whitened) plots OK\n\n")

err_empty <- tryCatch(plot_var_convergence(fit1, component = "precision"), error = function(e) e)
stopifnot(inherits(err_empty, "error"), grepl("disp_table", conditionMessage(err_empty)))
cat("fit1 component = 'precision' (empty disp_table) -> clear error OK:\n  ",
    conditionMessage(err_empty), "\n\n")

err_whitened <- tryCatch(
  plot_var_convergence(fit2, component = "precision", whitened = TRUE, n_chains = 5L),
  error = function(e) e
)
stopifnot(inherits(err_whitened, "error"), grepl("whitened", conditionMessage(err_whitened)))
cat("fit2 component = 'precision', whitened = TRUE -> clear error OK:\n  ",
    conditionMessage(err_whitened), "\n\n")

err_mean_whitened <- tryCatch(
  plot_mean_convergence(fit2, component = "precision", whitened = TRUE, n_chains = 5L),
  error = function(e) e
)
stopifnot(inherits(err_mean_whitened, "error"), grepl("whitened", conditionMessage(err_mean_whitened)))
cat("fit2 plot_mean_convergence component = 'precision', whitened = TRUE -> clear error OK\n\n")

print(fit2$sweep_history)
cat("print(fit2$sweep_history) shows RE precision table OK\n\n")

## --- rlmerb(): same underlying engine as fit1, via the model_setup()/
## Prior_Setup_GLMM() pipeline; result_class is lost (class == "rlmerb")
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
ps_g <- Prior_Setup_GLMM(form_lmer, data = dat, pwt = 0.05, dispformula = ~group)
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
  n = 15L, y = design$y, D = design$D, group = grp, W = design$W,
  pfamily_list = pf_known_vcov, dispprior_list = prior_list_disp_group,
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
grp_pois <- design_pois$group
attr(grp_pois, "group_name") <- design_pois$group_name
ps_pois <- Prior_Setup_GLMM(form_lmer, data = dat_pois, family = poisson(), pwt = 0.05)
pf_pois <- pfamily_list(ps_pois)
fit5 <- rGLMM_reg_known_vcov(
  n = 15L, y = design_pois$y, D = design_pois$D, group = grp_pois, W = design_pois$W,
  pfamily_list = pf_pois, dispprior_list = list(), family = poisson(),
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
