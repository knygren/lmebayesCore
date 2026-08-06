## TEMPORARY test runner (not part of the deliverable scratch diagnostic) --
## reproduces Ex_13c's setup/fit quietly (verbose = FALSE, no plots/prints)
## to smoke-test _scratch_lambda_star_marginal_over_draws.R quickly.

devtools::load_all(".", quiet = TRUE)

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design_all <- model_setup(form_lmer, data = dat)
full_rank_schools <- names(design_all$re_rank)[design_all$re_rank]
dat <- subset(dat, school_id %in% full_rank_schools)
dat$school_id <- droplevels(dat$school_id)

temp_drop_schools <- c("18", "2", "6")
drop <- intersect(temp_drop_schools, levels(dat$school_id))
dat <- subset(dat, !as.character(school_id) %in% drop)
dat$school_id <- droplevels(dat$school_id)

design <- model_setup(form_lmer, data = dat)
ps <- Prior_Setup_GLMM(
  form_lmer, data = dat, pwt = 0.01, dispformula = ~school_id,
  max_disp_perc_measurement = 0.8, pwt_measurement = 0.1
)
pf <- pfamily_list(ps)

max_disp_perc <- 0.8
block_formula <- ps$block_formula
sd_tau        <- sqrt(diag(ps$Sigma_ranef))
re_names_all  <- design$re_coef_names
group_levels0 <- levels(dat$school_id)

part_vi_group <- stats::setNames(
  lapply(group_levels0, function(lev) {
    dat_j <- dat[dat$school_id == lev, , drop = FALSE]
    inp <- lmebayesCore:::.lmebayes_ing_prior_measurement_group_glm_inputs(
      lev = lev, dat_j = dat_j, block_formula = block_formula, sd_tau = sd_tau,
      family = gaussian(), intercept_source = "null_model", effects_source = "null_effects"
    )
    n_j          <- inp$n_j
    n_prior_j    <- ps$ing_prior_measurement_group[[lev]]$n_prior
    n_combined_j <- n_prior_j + n_j
    p_re         <- length(sd_tau)

    pwt_j <- diag(inp$V0)
    pwt_j <- pwt_j / (pwt_j + inp$sd_vec^2)
    if (length(pwt_j) == 1L) {
      Sigma_j <- ((1 - pwt_j) / pwt_j) * inp$V0
    } else {
      scale_vec <- sqrt((1 - pwt_j) / pwt_j)
      Sigma_j <- inp$V0 * outer(scale_vec, scale_vec)
    }

    Omega_j <- matrix(0, nrow = length(inp$var_names), ncol = length(inp$var_names),
                       dimnames = list(inp$var_names, inp$var_names))
    for (k in re_names_all) {
      Wk_row <- design$W[[k]][lev, , drop = FALSE]
      Sigma_fixef_k <- ps$prior_list[[k]]$Sigma_fixef
      Omega_j[k, k] <- as.numeric(Wk_row %*% Sigma_fixef_k %*% t(Wk_row))
    }

    cal <- lmebayesCore:::.lmebayes_compute_ing_prior_cal_from_sigma(
      inp, Sigma_j + Omega_j, n_prior_j
    )

    shape_w <- (n_combined_j + 1) / 2 + p_re / 2
    rate_w  <- cal$dispersion * (n_combined_j + p_re - 1) / 2
    disp_lower <- 1 / qgamma(max_disp_perc,     shape = shape_w, rate = rate_w)
    disp_upper <- 1 / qgamma(1 - max_disp_perc, shape = shape_w, rate = rate_w)

    list(
      sigma2_hat = cal$dispersion, shape = cal$shape_ING, rate = cal$rate,
      disp_lower = disp_lower, disp_upper = disp_upper, omega_j = Omega_j
    )
  }),
  group_levels0
)

grp <- design$group
attr(grp, "group_name") <- design$group_name
group_levels <- levels(grp)
re_names     <- design$re_coef_names
p_re         <- length(re_names)

shape_group      <- stats::setNames(numeric(length(group_levels)), group_levels)
rate_group       <- stats::setNames(numeric(length(group_levels)), group_levels)
disp_lower_group <- stats::setNames(numeric(length(group_levels)), group_levels)
disp_upper_group <- stats::setNames(numeric(length(group_levels)), group_levels)
for (lev in group_levels) {
  g <- part_vi_group[[lev]]
  shape_group[[lev]]      <- g$shape
  rate_group[[lev]]       <- g$rate
  disp_lower_group[[lev]] <- g$disp_lower
  disp_upper_group[[lev]] <- g$disp_upper
}

prior_list <- list(
  mu               = matrix(0, nrow = p_re, ncol = 1L, dimnames = list(re_names, NULL)),
  Sigma            = as.matrix(ps$Sigma_ranef),
  shape_group      = shape_group,
  rate_group       = rate_group,
  disp_lower_group = disp_lower_group,
  disp_upper_group = disp_upper_group
)

cat("Fitting (verbose = FALSE, quiet)...\n")
fit <- rLMMindepNormalGamma_reg_known_vcov(
  n            = 1000L,
  y            = design$y,
  D            = design$D,
  group        = grp,
  W            = design$W,
  prior_list   = prior_list,
  pfamily_list = pf,
  progbar      = FALSE,
  verbose      = FALSE
)
cat("Fit done.\n")

n_draws <- nrow(fit$fixef[[re_names[1L]]])
n_group <- stats::setNames(as.numeric(table(grp)), group_levels)

prior_list_block1_rate <- list(
  Sigma      = as.matrix(ps$Sigma_ranef),
  dispersion = disp_upper_group[group_levels]
)
prior_list_block2_rate <- lapply(pf, function(pfk) {
  pl <- pfk$prior_list
  list(mu = pl$mu, Sigma = pl$Sigma, dispersion = pl$dispersion)
})

source("data-raw/_scratch_rss_ellipsoid_test.R")
tab <- .tmp_rss_ellipsoid_test(
  fit = fit, D = design$D, y = design$y, group = grp,
  group_name = design$group_name, re_coef_names = re_names,
  rate_group = rate_group
)
cat("\n--- RSS ellipsoid test (sanity check) ---\n")
print(tab, row.names = FALSE, digits = 4)

source("data-raw/_scratch_lambda_star_marginal_over_draws.R")
inp <- lmebayesCore:::.two_block_rate_inputs(
  x = design$D, block = grp, x_hyper = design$W,
  prior_list_block1 = prior_list_block1_rate,
  prior_list_block2 = prior_list_block2_rate
)
blocks <- lmebayesCore:::.two_block_S_P11(inp)
group_setup <- .tmp_marginal_group_setup(
  D = design$D, y = design$y, group = grp, group_levels = group_levels,
  re_coef_names = re_names, shape_group = shape_group, rate_group = rate_group
)
res <- .tmp_lambda_star_marginal_over_draws(
  fit = fit, n_draws = n_draws, y = design$y,
  group_name = design$group_name, group_setup = group_setup,
  inp = inp, blocks = blocks
)
.tmp_print_marginal_over_draws_summary(res)

cat(sprintf("\nBase lambda_star (fixed-Omega corner) for comparison = %.6f\n",
            lmebayesCore:::.two_block_gen_eigen(blocks$S, blocks$P11, strict = FALSE)[1L]))
