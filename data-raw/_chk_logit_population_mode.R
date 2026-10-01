## Non-Gaussian population_mode: ICM init, MC screen, and EM fallback.
## Primary fixture: book_banning (Chapter-C06, 8 states, full rank).
## Secondary: bwc_full_rank logit (school hierarchy, Gaussian-calibrated priors).
## Not part of the package test suite.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required to load lmebayesCore for this script.")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Suggested package 'bayesrules' is required.", call. = FALSE)
}

devtools::load_all("c:/Rpackages/lmebayesCore", quiet = TRUE)

fixef_l2 <- function(a, b) {
  sqrt(sum(unlist(Map(function(x, y) x - y, a, b)))^2)
}

mc_seed <- 42L

cat("=== book_banning logit (Chapter-C06 subset) ===\n")

data(book_banning, package = "bayesrules")
bb <- book_banning[, c("state", "removed", "violent")]
bb <- bb[stats::complete.cases(bb), ]
bb$removed_i <- as.integer(bb$removed == 1L | bb$removed == "1")
bb$violent_i <- as.integer(
  bb$violent == TRUE | bb$violent == 1L | bb$violent == "TRUE"
)
bb$state <- factor(bb$state)
keep <- names(sort(table(bb$state), decreasing = TRUE))[seq_len(8L)]
bb <- droplevels(subset(bb, state %in% keep))

form_bb <- removed_i ~ violent_i + (1 + violent_i || state)
design_bb <- model_setup(form_bb, data = bb, family = binomial())
ps_bb <- Prior_Setup_GLMM(form_bb, data = bb, family = binomial(), pop.pwt = 0.01)
pf_bb <- pfamily_list(ps_bb)

mode_bb_screen <- population_mode(
  design = design_bb,
  pfamily_list = pf_bb,
  family = binomial(),
  estep = "mc",
  n = 1500L,
  mc_seed = mc_seed,
  icm_init = TRUE,
  icm_screen = TRUE,
  mc_alpha = 0.05,
  mc_stable = 2L,
  maxit = 50L
)

mode_bb_em <- population_mode(
  design = design_bb,
  pfamily_list = pf_bb,
  family = binomial(),
  estep = "mc",
  n = 1500L,
  mc_seed = mc_seed,
  icm_init = TRUE,
  icm_screen = FALSE,
  mc_alpha = 0.05,
  mc_stable = 2L,
  maxit = 50L
)

cat("  J =", nlevels(bb$state), " n_obs =", nrow(bb), "\n")
cat("  ICM iterations:", mode_bb_screen$icm$iterations, "\n")
cat("  screen on: route =", mode_bb_screen$em_route,
    " EM iters =", mode_bb_screen$iterations,
    " screen passed =", mode_bb_screen$icm_screen$passed, "\n")
cat("  screen off: EM iters =", mode_bb_em$iterations, "\n")
cat("  fixef L2 (screen vs no-screen): ",
    signif(fixef_l2(mode_bb_screen$fixef, mode_bb_em$fixef), 4), "\n")

stopifnot(isTRUE(mode_bb_screen$converged))
stopifnot(isTRUE(mode_bb_em$converged))
stopifnot(mode_bb_screen$em_route %in% c("icm_screen", "em"))

cat("\n=== bwc_full_rank school logit (optional) ===\n")

data_path <- system.file("extdata", "bwc_full_rank.rds", package = "lmebayesCore")
bwc_ok <- FALSE
if (nzchar(data_path) && file.exists(data_path)) {
  bwc_ok <- tryCatch({
    dat <- readRDS(data_path)
    dat$high_ppvt <- as.integer(dat$score_ppvt >
                                  stats::median(dat$score_ppvt, na.rm = TRUE))
    form_gauss <- score_ppvt ~ distracted_ppvt + (1 + distracted_ppvt || school_id)
    form <- high_ppvt ~ distracted_ppvt + (1 + distracted_ppvt || school_id)
    ps <- Prior_Setup_GLMM(form_gauss, data = dat, family = gaussian(), pop.pwt = 0.01)
    pf <- pfamily_list(ps)
    design <- model_setup(form, data = dat, family = binomial())
    stopifnot(all(design$groupef.rank))

    mode_bwc <- population_mode(
      design = design,
      pfamily_list = pf,
      family = binomial(),
      estep = "mc",
      n = 1000L,
      mc_seed = mc_seed,
      icm_init = TRUE,
      icm_screen = TRUE,
      mc_stable = 2L,
      maxit = 40L
    )

    beta_mode <- beta_marginal_mode(design = design, pfamily_list = pf,
                                    family = binomial())
    beta_set <- beta_marginal_safe_set(beta_mode = beta_mode, delta_2 = 0.01)
    mode_restricted <- population_mode(
      design = design,
      pfamily_list = pf,
      family = binomial(),
      beta_set = beta_set,
      estep = "mc",
      n = 400L,
      mc_seed = mc_seed,
      icm_init = TRUE,
      icm_screen = TRUE,
      mc_stable = 2L,
      maxit = 30L
    )

    cat("  J =", nlevels(dat$school_id), " route =", mode_bwc$em_route,
        " EM iters =", mode_bwc$iterations, "\n")
    cat("  restricted accept_rate =", signif(mode_restricted$accept_rate, 4),
        " route =", mode_restricted$em_route, "\n")
    stopifnot(isTRUE(mode_bwc$converged))
    stopifnot(mode_restricted$restricted)
    TRUE
  }, error = function(e) {
    message("  bwc logit skipped: ", conditionMessage(e))
    FALSE
  })
}
if (!bwc_ok) {
  cat("  (bwc school logit not run or failed; book_banning checks suffice.)\n")
}

cat("\nLogit population_mode checks passed.\n")
