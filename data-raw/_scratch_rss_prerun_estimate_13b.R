## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## Same closed-form Gaussian pre-run estimate as
## data-raw/_scratch_rss_prerun_estimate.R, but built from Ex_13b's Part VI
## per-group calibration (model-derived Omega_j folded into sigma2_hat/rate)
## instead of Ex_13's plain dGamma_list(ps, disp_upper_anchor = "symmetric").
## Everything else (design, ps, fit_ref, Sigma_ranef) is identical to Ex_13 --
## only shape_group/rate_group and the sigma2_hat plug-in differ.

library(lmebayesCore)

## ---------------------------------------------------------------------------
## 1. Ex_13b setup only (through Section 2b's Part VI calibration) -- NO
##    sampler call. Mirrors Ex_13b lines 1-252.
## ---------------------------------------------------------------------------
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

temp_drop_schools <- c("18", "2")
drop <- intersect(temp_drop_schools, levels(dat$school_id))
dat <- subset(dat, !as.character(school_id) %in% drop)
dat$school_id <- droplevels(dat$school_id)

design <- model_setup(form_lmer, data = dat)

ps <- Prior_Setup_lmebayes(
  form_lmer,
  data            = dat,
  pwt             = 0.01,
  dispformula     = ~school_id,
  max_disp_perc_measurement = 0.8,
  pwt_measurement = 0.1
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
      sigma2_hat = cal$dispersion,
      shape      = cal$shape_ING,
      rate       = cal$rate,
      disp_lower = disp_lower,
      disp_upper = disp_upper,
      omega_j    = Omega_j
    )
  }),
  group_levels0
)

grp <- design$group
attr(grp, "group_name") <- design$group_name
group_levels <- levels(grp)
re_names <- design$re_coef_names

shape_group <- stats::setNames(numeric(length(group_levels)), group_levels)
rate_group  <- stats::setNames(numeric(length(group_levels)), group_levels)
for (lev in group_levels) {
  g <- part_vi_group[[lev]]
  shape_group[[lev]] <- g$shape
  rate_group[[lev]]  <- g$rate
}

cat(sprintf("\n=== Ex_13b (Part VI) setup complete: %d groups ===\n", length(group_levels)))

## ---------------------------------------------------------------------------
## 2. Same pre-run Monte Carlo estimator as Ex_13's scratch script (gamma_hat
##    mapped per RE component, plugged-in Omega_hat_j -- now from part_vi_group
##    instead of ps$ing_prior_measurement_group).
## ---------------------------------------------------------------------------
gamma_hat_raw <- lmebayesCore:::.lmebayes_reference_fixef(ps$fit_ref)
Sigma_ranef <- as.matrix(ps$Sigma_ranef)
Psi_inv <- diag(1 / diag(Sigma_ranef))
dimnames(Psi_inv) <- dimnames(Sigma_ranef)

.tmp_fe_name <- function(k, nm, fe_names) {
  if (identical(k, "(Intercept)") && identical(nm, "(Intercept)")) return("(Intercept)")
  if (identical(nm, "(Intercept)")) return(k)
  if (identical(k, "(Intercept)")) return(nm)
  cand <- c(paste0(nm, ":", k), paste0(k, ":", nm))
  hit <- cand[cand %in% fe_names]
  if (length(hit)) hit[1L] else NA_character_
}

gamma_hat_by_component <- stats::setNames(
  lapply(re_names, function(k) {
    Wk <- design$W[[k]]
    nm <- colnames(Wk)
    fe_nm <- vapply(nm, .tmp_fe_name, character(1), k = k, fe_names = names(gamma_hat_raw))
    if (anyNA(fe_nm)) stop("Could not map hyper-predictor(s) for '", k, "'.", call. = FALSE)
    stats::setNames(unname(gamma_hat_raw[fe_nm]), nm)
  }),
  re_names
)

.tmp_rss_prerun_estimate <- function(D, y, group, group_levels, re_coef_names,
                                      rate_group, ing_grp, x_hyper, gamma_hat_by_component,
                                      Psi_inv, n_sim = 20000L, seed = 1L) {
  set.seed(seed)
  group_chr <- as.character(group)

  rows <- lapply(group_levels, function(lev) {
    idx <- which(group_chr == lev)
    D_j <- D[idx, re_coef_names, drop = FALSE]
    y_j <- y[idx]
    DtD <- crossprod(D_j)
    beta_ols <- as.vector(solve(DtD, crossprod(D_j, y_j)))
    RSS_ols  <- sum((y_j - D_j %*% beta_ols)^2)
    thresh   <- 2 * rate_group[[lev]] + RSS_ols

    Omega_hat_j <- 1 / ing_grp[[lev]]$sigma2_hat
    mu_j <- vapply(re_coef_names, function(k) {
      Wk <- x_hyper[[k]]
      as.numeric(Wk[lev, , drop = FALSE] %*% gamma_hat_by_component[[k]])
    }, numeric(1))

    Sigma_j_inv <- Omega_hat_j * DtD + Psi_inv
    Sigma_j <- solve(Sigma_j_inv)
    beta_bar_j <- as.vector(Sigma_j %*% (Omega_hat_j * crossprod(D_j, y_j) + Psi_inv %*% mu_j))

    L <- chol(Sigma_j)
    Z <- matrix(rnorm(n_sim * length(re_coef_names)), nrow = n_sim)
    draws <- Z %*% L
    draws <- sweep(draws, 2, beta_bar_j, "+")

    diffs <- sweep(draws, 2, beta_ols, "-")
    q_draws <- rowSums((diffs %*% DtD) * diffs)

    d_mean <- beta_bar_j - beta_ols
    q_mean <- as.numeric(t(d_mean) %*% DtD %*% d_mean)

    data.frame(
      group = lev, n_j = length(idx), RSS_ols = RSS_ols, threshold = thresh,
      q_mean_prerun = q_mean, inside_at_mean_prerun = q_mean <= thresh,
      pct_outside_prerun = 100 * mean(q_draws > thresh)
    )
  })
  do.call(rbind, rows)
}

tab_prerun_13b <- .tmp_rss_prerun_estimate(
  D = design$D, y = design$y, group = grp, group_levels = group_levels,
  re_coef_names = re_names, rate_group = rate_group,
  ing_grp = part_vi_group, x_hyper = design$W,
  gamma_hat_by_component = gamma_hat_by_component, Psi_inv = Psi_inv
)

cat("\n=== Pre-run estimate: Ex_13b (Part VI) ===\n\n")
print(tab_prerun_13b, row.names = FALSE, digits = 4)

## ---------------------------------------------------------------------------
## 3. Compare to Ex_13's own pre-run estimate (Ex_13 vs Ex_13b, same groups).
##    Values below are Ex_13's pre-run pct_outside_prerun from
##    _scratch_rss_prerun_estimate.R's own run (recorded so this script does
##    not need to re-run Ex_13's setup from scratch).
## ---------------------------------------------------------------------------
ex13_pct_prerun <- c(
  `3` = 0.055, `5` = 0.145, `6` = 11.075, `9` = 0.800, `10` = 0.120, `11` = 1.045,
  `12` = 3.615, `13` = 1.050, `14` = 1.130, `15` = 0.375, `17` = 0.650, `19` = 0.705,
  `20` = 0.475, `22` = 0.370, `24` = 0.375, `25` = 0.080, `26` = 1.710, `29` = 0.005,
  `30` = 2.670, `31` = 0.140, `33` = 4.610, `35` = 0.475, `36` = 0.110, `37` = 0.150,
  `38` = 1.055, `39` = 0.035, `40` = 1.210, `42` = 0.385, `43` = 0.000, `44` = 0.250,
  `45` = 0.425, `47` = 0.790
)
cmp <- data.frame(
  group = tab_prerun_13b$group,
  pct_prerun_Ex13  = ex13_pct_prerun[tab_prerun_13b$group],
  pct_prerun_Ex13b = tab_prerun_13b$pct_outside_prerun
)
cmp$improved <- cmp$pct_prerun_Ex13b < cmp$pct_prerun_Ex13
cat("\n=== Ex_13 vs Ex_13b (Part VI) pre-run estimate, side by side ===\n\n")
print(cmp, row.names = FALSE, digits = 4)
cat(sprintf(
  "\nGroups improved (lower pct_outside_prerun) under Part VI: %d / %d\n",
  sum(cmp$improved), nrow(cmp)
))
cat(sprintf("Mean pct_outside_prerun -- Ex_13: %.3f%%, Ex_13b (Part VI): %.3f%%\n",
            mean(cmp$pct_prerun_Ex13), mean(cmp$pct_prerun_Ex13b)))
cat(sprintf("Max  pct_outside_prerun -- Ex_13: %.3f%% (group %s), Ex_13b (Part VI): %.3f%% (group %s)\n",
            max(cmp$pct_prerun_Ex13), cmp$group[which.max(cmp$pct_prerun_Ex13)],
            max(cmp$pct_prerun_Ex13b), cmp$group[which.max(cmp$pct_prerun_Ex13b)]))
