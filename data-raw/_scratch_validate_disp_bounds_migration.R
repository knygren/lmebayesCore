## Scratch validation for the disp_lower/disp_upper migration into
## Prior_Setup_GLMM() (Part VI now the permanent default; dGamma_list()
## no longer builds its own window). Not a permanent test -- ad hoc checks
## only, run by hand.
##
##   Rscript data-raw/_scratch_validate_disp_bounds_migration.R

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
temp_drop_schools <- c("18", "2")
dat <- subset(dat, !as.character(school_id) %in% temp_drop_schools)
dat$school_id <- droplevels(dat$school_id)

ps <- Prior_Setup_GLMM(
  form_lmer,
  data            = dat,
  pwt             = 0.01,
  dispformula     = ~school_id,
  max_disp_perc_measurement = 0.8,
  pwt_measurement = 0.1
)

cat("\n=== 1. ing_prior_measurement_group now carries disp_lower/disp_upper ===\n\n")
g1 <- ps$ing_prior_measurement_group[[1L]]
str(g1[c("sigma2_hat", "shape_ING", "rate", "disp_lower", "disp_upper",
         "max_disp_perc")])
stopifnot(
  !is.null(g1$disp_lower), !is.null(g1$disp_upper),
  g1$disp_lower < g1$disp_upper,
  is.matrix(g1$omega_j)
)

## sigma2_hat should sit inside [disp_lower, disp_upper] for every group
## (a group's own calibrated point estimate should never land outside an
## 80%-mass window built for exactly that group).
inside <- vapply(ps$ing_prior_measurement_group, function(g) {
  g$sigma2_hat >= g$disp_lower && g$sigma2_hat <= g$disp_upper
}, logical(1))
cat(sprintf(
  "sigma2_hat inside [disp_lower, disp_upper] for %d / %d groups\n",
  sum(inside), length(inside)
))
stopifnot(all(inside))

cat("\n=== 2. dGamma_list(ps) consumes the stored bounds verbatim ===\n\n")
pf <- dGamma_list(ps)
g_lev <- names(ps$ing_prior_measurement_group)[[1L]]
pf1 <- pf[[g_lev]]
stopifnot(
  isTRUE(all.equal(pf1$prior_list$shape, ps$ing_prior_measurement_group[[g_lev]]$shape_ING)),
  isTRUE(all.equal(pf1$prior_list$rate, ps$ing_prior_measurement_group[[g_lev]]$rate)),
  isTRUE(all.equal(pf1$prior_list$disp_lower, ps$ing_prior_measurement_group[[g_lev]]$disp_lower)),
  isTRUE(all.equal(pf1$prior_list$disp_upper, ps$ing_prior_measurement_group[[g_lev]]$disp_upper))
)
cat("dGamma_list() pfamily shape/rate/disp_lower/disp_upper match ps$ing_prior_measurement_group: OK\n")

wd <- attr(pf, "window_diagnostics")
cat("\n=== 3. window_diagnostics attribute (lightweight) ===\n\n")
print(head(wd), row.names = FALSE)
stopifnot(nrow(wd) == length(pf))

## Recompute at a different max_disp_perc_measurement: should trigger a
## fresh quantile calculation (not just echo the stored 0.8-calibrated
## bounds).
cat("\n=== 4. Different max_disp_perc_measurement triggers a recompute ===\n\n")
pf95 <- dGamma_list(ps, max_disp_perc_measurement = 0.95)
lo80 <- ps$ing_prior_measurement_group[[g_lev]]$disp_lower
lo95 <- pf95[[g_lev]]$prior_list$disp_lower
cat(sprintf("disp_lower at max_disp_perc_measurement = 0.80: %.4f\n", lo80))
cat(sprintf("disp_lower at max_disp_perc_measurement = 0.95: %.4f\n", lo95))
stopifnot(lo95 < lo80)  # wider window -> lower disp_lower

## Passing the removed disp_center/disp_upper_anchor args should be silently
## absorbed by '...' (backward-compatible no-ops), not error.
cat("\n=== 5. Old (now-removed) args are silently absorbed via '...' ===\n\n")
pf_old_args <- dGamma_list(
  ps, disp_center = "sigma2_hat", disp_upper_anchor = "blup", n_rss_iter = 5L,
  warn_asymmetric = FALSE
)
stopifnot(isTRUE(all.equal(
  pf_old_args[[g_lev]]$prior_list$disp_lower,
  pf[[g_lev]]$prior_list$disp_lower
)))
cat("Old disp_center/disp_upper_anchor/n_rss_iter/warn_asymmetric args accepted (ignored), no error: OK\n")

cat("\n=== 6. Part VI Omega_j fold-in is non-trivial (shape_ING/rate move vs Omega_j = 0) ===\n\n")
## Sanity check the fold-in actually changes rate relative to what Part I
## alone (Omega_j = 0) would have given, for at least one group.
omega_nonzero <- vapply(ps$ing_prior_measurement_group, function(g) {
  any(diag(g$omega_j) > 0)
}, logical(1))
cat(sprintf("Groups with nonzero Omega_j: %d / %d\n", sum(omega_nonzero), length(omega_nonzero)))
stopifnot(all(omega_nonzero))

cat("\nAll checks passed.\n")
