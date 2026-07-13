## Runnable companion to inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md
##
## Reproduces, on the `bayesrules::big_word_club` all-full-rank fixture
## (39 schools, 2 uncorrelated RE components):
##   1. The current dGamma_list() per-group truncation window (mean-matched at
##      sigma2_hat, asymmetric BLUP-inflated upper tail).
##   2. A proposed refinement that mean-matches at a dispersion estimate which
##      integrates over random-effect uncertainty (EnvelopeCentering-style
##      trace correction), with a worked illustration of the double-counting
##      pitfall the naive version of this refinement falls into.
##   3. An independent classical validation via nlme::lme() with a per-group
##      residual-variance (varIdent) structure.
##
## Requires: bayesrules, nlme (both Suggests-only; script is not run by
## R CMD check). See inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md for the derivation
## this script checks numerically.

devtools::load_all(".", quiet = TRUE)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install 'bayesrules' to run this script.", call. = FALSE)
}

## ---------------------------------------------------------------------------
## 0. Fixture: all full-rank schools, 2 uncorrelated RE components
## ---------------------------------------------------------------------------
data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c("score_ppvt", "distracted_ppvt", "school_id")])
)

form <- score_ppvt ~ 1 + distracted_ppvt + (1 + distracted_ppvt || school_id)

design_scr <- model_setup(form, data = dat)
full_rank_schools <- names(design_scr$re_rank)[design_scr$re_rank]
dropped <- names(design_scr$re_rank)[!design_scr$re_rank]
if (length(dropped)) {
  message(
    "Dropping ", length(dropped), " rank-deficient school_id level(s): ",
    paste(dropped, collapse = ", ")
  )
  dat <- subset(dat, school_id %in% full_rank_schools)
  dat$school_id <- droplevels(dat$school_id)
}

design <- model_setup(form, data = dat)
stopifnot(all(design$re_rank))
group_levels <- levels(design$groups)
p_re <- length(design$re_coef_names)
max_disp_perc <- 0.8

message(sprintf(
  "%d schools, %d obs (mean %.1f obs/school)",
  nlevels(design$groups), nrow(dat), mean(table(dat$school_id))
))

## ---------------------------------------------------------------------------
## 1. Current dGamma_list() bounds (mean-matched at sigma2_hat, BLUP upper tail)
## ---------------------------------------------------------------------------
ps <- Prior_Setup_lmebayes(
  form, data = dat, pwt = 0.01, pwt_measurement = 0.1,
  max_disp_perc = max_disp_perc
)
disp_pf_current <- dGamma_list(ps, max_disp_perc = max_disp_perc, warn_asymmetric = FALSE)
window_diag <- attr(disp_pf_current, "window_diagnostics")

ing_grp       <- ps$ing_prior_measurement_group
block_formula <- ps$block_formula
Sigma_ranef   <- ps$Sigma_ranef
P             <- solve(Sigma_ranef)
RA            <- chol(P)
fit_ref       <- ps$fit_ref
group_name    <- ps$design$group_name
beta_blup     <- coef(fit_ref)[[group_name]]
sigma2_pooled_lmer <- stats::sigma(fit_ref)^2

## ---------------------------------------------------------------------------
## 2. Proposed refinement: mean-match at an EnvelopeCentering-style dispersion
##    estimate that integrates over b_j uncertainty via a trace correction.
##
##    IMPORTANT (the double-counting pitfall, see README Part III): the
##    dispersion2 fixed-point below uses ONLY the group's own n_j/RSS_precomputed
##    at each iteration -- it must NOT be seeded by feeding the group's own
##    shape_ING/rate_gamma (already calibrated from the same n_j observations)
##    back in as a second "prior" contribution. Doing so double-counts the
##    group's data and produces artificially tight windows (verified below to
##    push at least one group's BLUP point estimate outside its own window).
## ---------------------------------------------------------------------------
group_calc <- function(lev) {
  idx   <- design$groups == lev
  dat_j <- dat[idx, , drop = FALSE]
  X <- stats::model.matrix(block_formula, data = dat_j)
  Y <- stats::model.response(stats::model.frame(block_formula, data = dat_j))
  n_j <- nrow(X)

  mu_j <- glmbayesCore:::.lmebayes_block_formula_prior_mu(
    block_formula = block_formula, dat_j = dat_j,
    intercept_source = "null_model", effects_source = "null_effects"
  )
  mu_vec <- as.numeric(mu_j)

  fit0 <- lm.fit(X, Y)
  rss0 <- sum(fit0$residuals^2)
  sigma2_ols <- rss0 / (n_j - ncol(X))

  z_bot <- RA %*% mu_vec
  XtX   <- t(X) %*% X

  dispersion2 <- sigma2_ols
  for (iter in 1:10) {
    s <- 1 / sqrt(dispersion2)
    W <- rbind(s * X, RA)
    z <- c(s * Y, z_bot)
    Sigma_post  <- solve(t(W) %*% W)
    b2          <- Sigma_post %*% (t(W) %*% z)
    r           <- Y - X %*% b2
    rss_at_mean <- sum(r^2)
    trace_term  <- sum(diag(XtX %*% Sigma_post))
    RSS_precomputed <- rss_at_mean + trace_term
    dispersion2 <- RSS_precomputed / (n_j - ncol(X))
  }

  g <- ing_grp[[lev]]
  n_combined   <- g$n_combined
  shape_w_prop <- (n_combined + 1) / 2 + p_re / 2
  rate_w_prop  <- dispersion2 * (n_combined + p_re - 1) / 2  # mean-matched, symmetric

  disp_lower_prop <- 1 / qgamma(max_disp_perc,     shape = shape_w_prop, rate = rate_w_prop)
  disp_upper_prop <- 1 / qgamma(1 - max_disp_perc, shape = shape_w_prop, rate = rate_w_prop)

  beta_lev   <- as.numeric(beta_blup[lev, colnames(X), drop = FALSE])
  RSS_blup   <- sum((Y - X %*% beta_lev)^2)
  df_resid   <- n_j - ncol(X)

  list(
    n_j = n_j, dispersion2 = dispersion2,
    sigma2_ols  = sigma2_ols,
    sigma2_blup = RSS_blup / df_resid,
    disp_lower_prop = disp_lower_prop,
    disp_upper_prop = disp_upper_prop
  )
}

rows <- lapply(group_levels, function(lev) {
  rp <- group_calc(lev)
  pf <- disp_pf_current[[lev]]
  wdiag <- window_diag[window_diag$group == lev, ]
  data.frame(
    group = lev, n_j = rp$n_j,
    sigma2_ols = rp$sigma2_ols, sigma2_blup = rp$sigma2_blup,
    blup_infl = wdiag$blup_infl, dispersion2_prop = rp$dispersion2,
    disp_lower_cur = pf$prior_list$disp_lower, disp_upper_cur = pf$prior_list$disp_upper,
    disp_lower_prop = rp$disp_lower_prop, disp_upper_prop = rp$disp_upper_prop,
    blup_inside_prop = rp$sigma2_blup >= rp$disp_lower_prop & rp$sigma2_blup <= rp$disp_upper_prop,
    blup_inside_cur  = rp$sigma2_blup >= pf$prior_list$disp_lower & rp$sigma2_blup <= pf$prior_list$disp_upper,
    stringsAsFactors = FALSE
  )
})
tab <- do.call(rbind, rows)
tab$pct_delta_lower <- 100 * (tab$disp_lower_prop - tab$disp_lower_cur) / tab$disp_lower_cur
tab$pct_delta_upper <- 100 * (tab$disp_upper_prop - tab$disp_upper_cur) / tab$disp_upper_cur

cat("\n=== Current vs. proposed per-group dGamma() bounds (max_disp_perc = 0.8) ===\n\n")
print(round(tab[order(-tab$blup_infl), setdiff(names(tab), "group")], 2), row.names = FALSE)

cat(sprintf(
  "\nGroups where the BLUP point estimate falls OUTSIDE the proposed window: %s\n",
  paste(tab$group[!tab$blup_inside_prop], collapse = ", ")
))
cat(sprintf(
  "Median %%delta lower = %.1f%%, median %%delta upper = %.1f%%\n",
  median(tab$pct_delta_lower), median(tab$pct_delta_upper)
))

## ---------------------------------------------------------------------------
## 3. Independent classical validation: nlme::lme() with per-group residual
##    variance (varIdent), vs. lme4's single pooled sigma^2.
## ---------------------------------------------------------------------------
if (requireNamespace("nlme", quietly = TRUE)) {
  fit_hetero <- nlme::lme(
    score_ppvt ~ 1 + distracted_ppvt,
    random  = list(school_id = nlme::pdDiag(~ 1 + distracted_ppvt)),
    weights = nlme::varIdent(form = ~ 1 | school_id),
    data    = dat,
    control = nlme::lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim")
  )
  base_sigma <- fit_hetero$sigma
  vw <- nlme::varWeights(fit_hetero$modelStruct$varStruct)
  dat$.sigma_j <- base_sigma / vw
  sigma2_nlme <- tapply(dat$.sigma_j^2, dat$school_id, function(x) x[1])
  tab$sigma2_nlme <- unname(sigma2_nlme[tab$group])

  lrt_stat <- 2 * (as.numeric(logLik(fit_hetero)) - as.numeric(logLik(fit_ref)))
  df_added <- length(group_levels) - 1
  cat(sprintf(
    "\n=== nlme::lme() heteroscedastic validation ===\nLRT (pooled lme4 vs. per-group nlme) = %.2f on %d df, p = %.4g\n",
    lrt_stat, df_added, pchisq(lrt_stat, df_added, lower.tail = FALSE)
  ))
  cat(sprintf(
    "cor(sigma2_ols, sigma2_nlme) = %.3f | cor(sigma2_hat_calibrated, sigma2_nlme) = %.3f\n",
    cor(tab$sigma2_ols, tab$sigma2_nlme), cor(window_diag$sigma2_hat, tab$sigma2_nlme)
  ))
} else {
  message("Install 'nlme' to run the classical validation step.")
}
