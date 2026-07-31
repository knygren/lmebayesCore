## Scratch: end-to-end sanity check that rLMMindepNormalGamma_reg_known_vcov()
## still runs cleanly (envelope construction, accept/reject, pilot + main
## stage) against the migrated Prior_Setup_lmebayes()/dGamma_list() bounds.
## Small n, verbose = FALSE -- speed/log-flooding is not the point here,
## just "does it run without error, on the exact Ex_13 fixture/priors".
devtools::load_all(".", quiet = TRUE)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

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
dat <- subset(dat, !as.character(school_id) %in% c("18", "2"))
dat$school_id <- droplevels(dat$school_id)

design <- model_setup(form_lmer, data = dat)
stopifnot(all(design$re_rank))

ps <- Prior_Setup_lmebayes(
  form_lmer, data = dat, pwt = 0.01, dispformula = ~school_id,
  max_disp_perc_measurement = 0.8, pwt_measurement = 0.1
)
pf <- pfamily_list(ps)
disp_pf_list <- dGamma_list(ps, max_disp_perc_measurement = 0.8)

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
  pl <- disp_pf_list[[lev]]$prior_list
  shape_group[[lev]]      <- pl$shape[1L]
  rate_group[[lev]]       <- pl$rate[1L]
  disp_lower_group[[lev]] <- pl$disp_lower
  disp_upper_group[[lev]] <- pl$disp_upper
}

prior_list <- list(
  mu               = matrix(0, nrow = p_re, ncol = 1L, dimnames = list(re_names, NULL)),
  Sigma            = as.matrix(ps$Sigma_ranef),
  shape_group      = shape_group,
  rate_group       = rate_group,
  disp_lower_group = disp_lower_group,
  disp_upper_group = disp_upper_group
)

cat(sprintf(
  "=== Running rLMMindepNormalGamma_reg_known_vcov(n = 50) on %d groups ===\n",
  length(group_levels)
))

t0 <- Sys.time()
fit <- rLMMindepNormalGamma_reg_known_vcov(
  n            = 50L,
  y            = design$y,
  D            = design$D,
  group        = grp,
  W            = design$W,
  prior_list   = prior_list,
  pfamily_list = pf,
  progbar      = FALSE,
  verbose      = FALSE
)
cat(sprintf("Completed in %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))

stopifnot(is.matrix(fit$dispersion_ranef))
stopifnot(all(is.finite(fit$dispersion_ranef)), all(fit$dispersion_ranef > 0))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(!is.null(fit$sweep_history))
n_draws <- nrow(fit$fixef[[re_names[1L]]])
stopifnot(n_draws == 50L)

cat(sprintf(
  "OK: %d draws, dispersion_ranef range [%.3f, %.3f], m_convergence = %s\n",
  n_draws, min(fit$dispersion_ranef), max(fit$dispersion_ranef),
  if (!is.null(fit$convergence_info$m_convergence)) fit$convergence_info$m_convergence else NA
))
if (!is.null(fit$convergence_info$lambda_star_marginal)) {
  cat(sprintf(
    "lambda_star_marginal (pilot safeguard) = %.4f, valid = %s\n",
    fit$convergence_info$lambda_star_marginal,
    fit$convergence_info$marginal_rate_valid
  ))
}
cat("\nSampler runs cleanly against the migrated bounds: OK\n")
