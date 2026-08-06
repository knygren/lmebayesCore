## Smoke test: $design is now populated on the 6 fit-object exports touched
## by the "convergence plot methods" plan step 4, and rGLMM_reg_known_vcov()/
## rGLMM_reg_estimated_vcov() share a new "rGLMM_reg" parent class.
##
## Uses the real Prior_Setup_GLMM()/pfamily_list()/dGamma_list() pipeline
## (same pattern as demo/Ex_10...Ex_14), on a small synthetic random-intercept
## dataset, so priors are well-scaled to the data (no hand-rolled prior_list/
## pfamily_list numbers) and the whole thing still runs in seconds.
devtools::load_all(".", quiet = TRUE)

## y ~ x1 + (1 + x1 || group): x1 is a fixed effect *and* a random slope, so
## it is an allowed "population mean slope" (a plain level-1 fixed effect
## with no matching random slope is rejected by model_setup(); see
## extract_re_hyper_matrices()'s "Fixed effects must be constant within
## <group>..." check) -- mirrors how every Ex_10...Ex_14 demo's fixed effects
## are either group-level covariates or population mean slopes matching a
## random-slope name (e.g. distracted_ppvt/distracted_a1 in the BigWordClub
## model).
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

## --- Known observation dispersion + known/estimated vcov priors ------
ps <- Prior_Setup_GLMM(form_lmer, data = dat, pwt = 0.05)
pf_known_vcov <- pfamily_list(ps)
pf_est_vcov   <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
prior_list_disp <- list(dispersion = ps$dispersion_ranef)

## --- Estimated (per-group) observation dispersion prior ---------------
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

check_design <- function(fit, label, expect_classes) {
  cat(sprintf("--- %s ---\n", label))
  cat("class:", paste(class(fit), collapse = ", "), "\n")
  has_design <- !is.null(fit$design)
  cat("has $design:", has_design, "\n")
  stopifnot(has_design)
  stopifnot(all(c("y", "D", "group", "W", "re_coef_names") %in% names(fit$design)))
  for (cls in expect_classes) {
    stopifnot(cls %in% class(fit))
  }
  cat("OK\n\n")
}

## --- rLMMNormal_reg_known_vcov(sim_method = "TWO_BLOCK_GIBBS") --------
fit1 <- rLMMNormal_reg_known_vcov(
  n = 20L, y = design$y, D = design$D, group = grp, W = design$W,
  prior_list = prior_list_disp, pfamily_list = pf_known_vcov,
  progbar = FALSE, verbose = FALSE, sim_method = "TWO_BLOCK_GIBBS"
)
check_design(fit1, "rLMMNormal_reg_known_vcov(TWO_BLOCK_GIBBS)", c("rLMMNormal_reg_known_vcov", "rLMMNormal_reg"))

## --- rLMMNormal_reg_known_vcov(sim_method = "DEFAULT", the actual default) ---
## .rLMMNormal_reg_run_iid() (the exact-iid engine behind "DEFAULT") was
## missing 'design' -- only the "TWO_BLOCK_GIBBS" route (.rLMMNormal_reg_run())
## had it. Exercise this route explicitly so a regression here is caught.
fit1b <- rLMMNormal_reg_known_vcov(
  n = 20L, y = design$y, D = design$D, group = grp, W = design$W,
  prior_list = prior_list_disp, pfamily_list = pf_known_vcov,
  progbar = FALSE, verbose = FALSE, sim_method = "DEFAULT"
)
check_design(fit1b, "rLMMNormal_reg_known_vcov(DEFAULT)", c("rLMMNormal_reg_known_vcov", "rLMMNormal_reg"))

## --- rLMMNormal_reg_estimated_vcov() ----------------------------------
fit2 <- rLMMNormal_reg_estimated_vcov(
  n = 20L, y = design$y, D = design$D, group = grp, W = design$W,
  prior_list = prior_list_disp, pfamily_list = pf_est_vcov,
  progbar = FALSE, verbose = FALSE
)
check_design(fit2, "rLMMNormal_reg_estimated_vcov()", c("rLMMNormal_reg_estimated_vcov", "rLMMNormal_reg"))

## --- rLMMindepNormalGamma_reg_known_vcov()/_estimated_vcov() ----------
fit3 <- rLMMindepNormalGamma_reg_known_vcov(
  n = 20L, y = design$y, D = design$D, group = grp, W = design$W,
  prior_list = prior_list_disp_group, pfamily_list = pf_known_vcov,
  progbar = FALSE, verbose = FALSE
)
check_design(fit3, "rLMMindepNormalGamma_reg_known_vcov()", c("rLMMindepNormalGamma_reg_known_vcov", "rLMMindepNormalGamma_reg"))

fit4 <- rLMMindepNormalGamma_reg_estimated_vcov(
  n = 20L, y = design$y, D = design$D, group = grp, W = design$W,
  prior_list = prior_list_disp_group, pfamily_list = pf_est_vcov,
  progbar = FALSE, verbose = FALSE
)
check_design(fit4, "rLMMindepNormalGamma_reg_estimated_vcov()", c("rLMMindepNormalGamma_reg_estimated_vcov", "rLMMindepNormalGamma_reg"))

## --- rGLMM_reg_known_vcov()/_estimated_vcov() -------------------------
## glmerb() (lmebayes) only ever routes non-Gaussian families to
## rGLMM_reg_known_vcov()/_estimated_vcov() -- Gaussian goes through
## rLMMNormal_reg_*()/rLMMindepNormalGamma_reg_*() instead (fit1...fit4
## above). Use a Poisson GLMM here so this matches how glmerb() actually
## calls these two exports (see inst/ARCHITECTURE_glmerb.md).
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

ps_pois        <- Prior_Setup_GLMM(form_lmer, data = dat_pois, family = poisson(), pwt = 0.05)
pf_known_pois  <- pfamily_list(ps_pois)
pf_est_pois    <- pfamily_list(ps_pois, ptypes = "dIndependent_Normal_Gamma")

fit5 <- rGLMM_reg_known_vcov(
  n = 20L, y = design_pois$y, D = design_pois$D, group = grp_pois, W = design_pois$W,
  prior_list = list(), pfamily_list = pf_known_pois,
  family = poisson(), progbar = FALSE, verbose = FALSE
)
check_design(fit5, "rGLMM_reg_known_vcov(poisson)", c("rGLMM_reg_known_vcov", "rGLMM_reg"))

fit6 <- rGLMM_reg_estimated_vcov(
  n = 20L, y = design_pois$y, D = design_pois$D, group = grp_pois, W = design_pois$W,
  prior_list = list(), pfamily_list = pf_est_pois,
  family = poisson(), progbar = FALSE, verbose = FALSE
)
check_design(fit6, "rGLMM_reg_estimated_vcov(poisson)", c("rGLMM_reg_estimated_vcov", "rGLMM_reg"))

cat("ALL OK\n")
