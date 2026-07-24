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

design <- model_setup(form_lmer, data = dat)
ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)
pf <- pfamily_list(ps)
grp <- design$groups
attr(grp, "group_name") <- design$group_name
re_names <- design$re_coef_names
prior_list_block1 <- list(
  P          = lmebayesCore:::.rLMM_P_from_pfamily_list(pf, re_names),
  dispersion = ps$dispersion_ranef,
  ddef       = FALSE
)

full_rank_levs <- names(design$re_rank)[design$re_rank]
cat(sprintf(
  "\nFull rank: %d of %d groups (%d rank-deficient)\n",
  length(full_rank_levs), length(design$re_rank),
  sum(!design$re_rank)
))

## (a) FULL model: all 47 groups (exactly what Ex_10 actually samples).
rate_full <- two_block_rate_from_pfamily_list(
  x                 = design$Z,
  block             = grp,
  x_hyper           = design$X_hyper,
  prior_list_block1 = prior_list_block1,
  pfamily_list      = pf,
  family            = gaussian(),
  group_levels      = levels(grp)
)
cat("\n=== (a) ALL 47 groups (incl. 13 rank-deficient) ===\n")
print(rate_full)

## (b) Restricted model: drop the 13 rank-deficient groups entirely (both
## from the data AND the group_levels used to build W_j's rows) -- tests
## whether their presence in (a) meaningfully raises lambda_star relative to
## a model built from only the well-identified groups.
dat_fr <- subset(dat, school_id %in% full_rank_levs)
dat_fr$school_id <- droplevels(dat_fr$school_id)
design_fr <- model_setup(form_lmer, data = dat_fr)
stopifnot(all(design_fr$re_rank))
grp_fr <- design_fr$groups
attr(grp_fr, "group_name") <- design_fr$group_name

rate_fr <- two_block_rate_from_pfamily_list(
  x                 = design_fr$Z,
  block             = grp_fr,
  x_hyper           = design_fr$X_hyper,
  prior_list_block1 = list(
    P          = prior_list_block1$P,
    dispersion = prior_list_block1$dispersion,
    ddef       = FALSE
  ),
  pfamily_list      = pf,
  family            = gaussian(),
  group_levels      = levels(grp_fr)
)
cat("\n=== (b) Only the 34 full-rank groups ===\n")
print(rate_fr)

cat("\n=== Comparison ===\n")
cat(sprintf("lambda_star (a) all 47 groups      : %.6f\n", rate_full$lambda_star))
cat(sprintf("lambda_star (b) 34 full-rank groups: %.6f\n", rate_fr$lambda_star))
cat(sprintf("difference (a - b)                 : %+.6f\n", rate_full$lambda_star - rate_fr$lambda_star))

cat("\nFull eigenvalue spectra:\n")
print(round(rbind(
  all_47_groups     = sort(rate_full$eigenvalues, decreasing = TRUE),
  full_rank_34_only = sort(rate_fr$eigenvalues, decreasing = TRUE)
), 6))
