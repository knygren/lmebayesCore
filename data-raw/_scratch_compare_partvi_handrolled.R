## Scratch: compare Ex_13b's hand-rolled Part VI derivation against the new
## Prior_Setup_GLMM() default (should now match on shape_ING/rate/
## sigma2_hat -- only the disp_lower/disp_upper window construction changed).
devtools::load_all(".", quiet = TRUE)

data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) & !is.na(invalid_ppvt) & invalid_ppvt == 0L &
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

ps <- Prior_Setup_GLMM(
  form_lmer, data = dat, pwt = 0.01, dispformula = ~school_id,
  max_disp_perc_measurement = 0.8, pwt_measurement = 0.1
)

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
    n_prior_j <- ps$ing_prior_measurement_group[[lev]]$n_prior
    pwt_j <- diag(inp$V0)
    pwt_j <- pwt_j / (pwt_j + inp$sd_vec^2)
    if (length(pwt_j) == 1L) {
      Sigma_j <- ((1 - pwt_j) / pwt_j) * inp$V0
    } else {
      scale_vec <- sqrt((1 - pwt_j) / pwt_j)
      Sigma_j <- inp$V0 * outer(scale_vec, scale_vec)
    }
    Omega_j <- matrix(
      0, nrow = length(inp$var_names), ncol = length(inp$var_names),
      dimnames = list(inp$var_names, inp$var_names)
    )
    for (k in re_names_all) {
      Wk_row <- design$W[[k]][lev, , drop = FALSE]
      Sigma_fixef_k <- ps$prior_list[[k]]$Sigma_fixef
      Omega_j[k, k] <- as.numeric(Wk_row %*% Sigma_fixef_k %*% t(Wk_row))
    }
    cal <- lmebayesCore:::.lmebayes_compute_ing_prior_cal_from_sigma(
      inp, Sigma_j + Omega_j, n_prior_j
    )
    list(sigma2_hat = cal$dispersion, shape = cal$shape_ING, rate = cal$rate)
  }),
  group_levels0
)

diffs <- sapply(group_levels0, function(lev) {
  a <- part_vi_group[[lev]]
  b <- ps$ing_prior_measurement_group[[lev]]
  c(
    shape_diff   = a$shape - b$shape_ING,
    rate_diff    = a$rate - b$rate,
    sigma2_diff  = a$sigma2_hat - unname(b$sigma2_hat)
  )
})
cat("Max abs diffs (hand-rolled Ex_13b Part VI vs new package default):\n")
print(apply(abs(diffs), 1, max))
stopifnot(all(apply(abs(diffs), 1, max) < 1e-8))
cat("\nMatch confirmed: new Prior_Setup_GLMM() default reproduces Ex_13b's\n")
cat("hand-rolled Part VI shape_ING/rate/sigma2_hat exactly.\n")
