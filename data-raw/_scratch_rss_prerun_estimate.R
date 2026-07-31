## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## Can we predict, BEFORE running the (expensive) two-block Gibbs sampler,
## something close to the empirical "pct_draws_outside" column from the RSS
## ellipsoid test?
##
## Idea: for a plug-in point estimate of the group's own dispersion precision
## Omega_hat_j (already available pre-run from Prior_Setup_lmebayes()'s own
## calibration, ps$ing_prior_measurement_group[[j]]$sigma2_hat) and a plug-in
## gamma_hat (the glmmTMB/lmer reference fit's fixef, also already available
## pre-run), the model's OWN Block~1 full conditional for beta_j is EXACTLY
## Gaussian (Prior_Setup_lmebayes.R's documented formula):
##
##   beta_j | Omega_hat_j, gamma_hat  ~  N(beta_bar_j, Sigma_j)
##   Sigma_j^{-1} = Omega_hat_j * D_j'D_j + Psi^{-1}      (Psi = Sigma_ranef, known)
##   beta_bar_j   = Sigma_j (Omega_hat_j * D_j'y_j + Psi^{-1} mu_j),  mu_j = W_j gamma_hat
##
## So instead of running MCMC, just draw a large number of samples directly
## from this closed-form Gaussian (near-instant) and compute the fraction
## with q_j(beta) > threshold_j -- a fast PRE-RUN estimate of the sampler's
## own empirical pct_draws_outside. This necessarily ignores the sampler's
## additional Omega_j/gamma uncertainty (both held fixed at their plug-in
## values here), so it is an approximation, not a certified bound.

## ---------------------------------------------------------------------------
## 1. Ex_13 setup only (model_setup/Prior_Setup_lmebayes/dGamma_list + the
##    prior_list/shape_group/rate_group construction) -- NO sampler call.
## ---------------------------------------------------------------------------
library(lmebayesCore)

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

disp_pf_list <- dGamma_list(
  ps,
  max_disp_perc_measurement = 0.8,
  disp_upper_anchor = "symmetric"
)

grp <- design$group
attr(grp, "group_name") <- design$group_name
group_levels <- levels(grp)
re_names <- design$re_coef_names
p_re <- length(re_names)

shape_group <- stats::setNames(numeric(length(group_levels)), group_levels)
rate_group  <- stats::setNames(numeric(length(group_levels)), group_levels)
for (lev in group_levels) {
  pl <- disp_pf_list[[lev]]$prior_list
  shape_group[[lev]] <- pl$shape[1L]
  rate_group[[lev]]  <- pl$rate[1L]
}

cat(sprintf("\n=== Ex_13 setup complete: %d groups, p_re = %d ===\n", length(group_levels), p_re))

## ---------------------------------------------------------------------------
## 2. Pre-run Monte Carlo estimate of pct_draws_outside
## ---------------------------------------------------------------------------
gamma_hat_raw <- lmebayesCore:::.lmebayes_reference_fixef(ps$fit_ref)
Sigma_ranef <- as.matrix(ps$Sigma_ranef)
Psi_inv <- diag(1 / diag(Sigma_ranef))
dimnames(Psi_inv) <- dimnames(Sigma_ranef)

## Map each RE component k's hyper-predictor name 'nm' to the correct entry
## of fixef(fit_ref): a main effect for intercept-associated predictors, but
## the OBSERVATION-LEVEL INTERACTION TERM for a non-intercept RE component's
## own hyper-predictors -- same mapping Ex_13's own Section 7 table uses.
.tmp_fe_name <- function(k, nm, fe_names) {
  if (identical(k, "(Intercept)") && identical(nm, "(Intercept)")) {
    return("(Intercept)")
  }
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
    if (anyNA(fe_nm)) {
      stop("Could not map hyper-predictor(s) for '", k, "' to fixef(fit_ref).", call. = FALSE)
    }
    stats::setNames(unname(gamma_hat_raw[fe_nm]), nm)
  }),
  re_names
)
cat("\n=== gamma_hat, correctly mapped per RE component ===\n\n")
print(gamma_hat_by_component)

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

    ## Plug-ins, all available BEFORE running the sampler:
    Omega_hat_j <- 1 / ing_grp[[lev]]$sigma2_hat
    mu_j <- vapply(re_coef_names, function(k) {
      Wk <- x_hyper[[k]]
      as.numeric(Wk[lev, , drop = FALSE] %*% gamma_hat_by_component[[k]])
    }, numeric(1))

    Sigma_j_inv <- Omega_hat_j * DtD + Psi_inv
    Sigma_j <- solve(Sigma_j_inv)
    beta_bar_j <- as.vector(Sigma_j %*% (Omega_hat_j * crossprod(D_j, y_j) + Psi_inv %*% mu_j))

    ## Closed-form Gaussian full conditional -> instant simulation (no MCMC).
    L <- chol(Sigma_j)  ## upper-triangular; t(L) %*% L = Sigma_j
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

tab_prerun <- .tmp_rss_prerun_estimate(
  D = design$D, y = design$y, group = grp, group_levels = group_levels,
  re_coef_names = re_names, rate_group = rate_group,
  ing_grp = ps$ing_prior_measurement_group, x_hyper = design$W,
  gamma_hat_by_component = gamma_hat_by_component, Psi_inv = Psi_inv
)

cat("\n=== Pre-run (closed-form Gaussian plug-in) estimate of pct_draws_outside ===\n\n")
print(tab_prerun, row.names = FALSE, digits = 4)

## ---------------------------------------------------------------------------
## 3. Compare to the ACTUAL post-sampler empirical table already reported.
## ---------------------------------------------------------------------------
actual_pct <- c(
  `3` = 0.0, `5` = 0.1, `6` = 8.3, `9` = 0.6, `10` = 0.3, `11` = 0.6, `12` = 3.2,
  `13` = 1.2, `14` = 1.5, `15` = 0.5, `17` = 0.8, `19` = 0.2, `20` = 0.2, `22` = 0.1,
  `24` = 0.1, `25` = 0.2, `26` = 1.0, `29` = 0.0, `30` = 1.9, `31` = 0.3, `33` = 4.2,
  `35` = 0.6, `36` = 0.1, `37` = 0.0, `38` = 0.6, `39` = 0.0, `40` = 1.4, `42` = 0.6,
  `43` = 0.0, `44` = 0.5, `45` = 0.8, `47` = 0.7
)
cmp <- data.frame(
  group = tab_prerun$group,
  pct_outside_prerun = tab_prerun$pct_outside_prerun,
  pct_draws_outside_actual = actual_pct[tab_prerun$group]
)
cat("\n=== Pre-run estimate vs actual post-sampler empirical pct ===\n\n")
print(cmp, row.names = FALSE, digits = 4)
cat(sprintf("\nCorrelation (pre-run estimate, actual): %.3f\n",
            cor(cmp$pct_outside_prerun, cmp$pct_draws_outside_actual)))
cat(sprintf("Rank correlation (Spearman): %.3f\n",
            cor(cmp$pct_outside_prerun, cmp$pct_draws_outside_actual, method = "spearman")))
