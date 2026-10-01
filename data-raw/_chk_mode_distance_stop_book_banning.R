## EM trace via population_mode(verbose=TRUE) on book_banning logit.
## Not part of the package test suite.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Suggested package 'bayesrules' is required.", call. = FALSE)
}

pkg_root <- if (dir.exists("R")) "." else if (dir.exists("../R")) ".." else
  "c:/Rpackages/lmebayesCore"
devtools::load_all(pkg_root, quiet = TRUE)

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

beta_set_bb <- beta_marginal_safe_set(
  delta_2 = 0.05,
  beta_mode = beta_marginal_mode(
    design = design_bb, pfamily_list = pf_bb, family = binomial()
  )
)

cat("=== book_banning logit (8 states) ===\n")
cat("J =", nlevels(bb$state), " n_obs =", nrow(bb), "\n\n")

mode <- population_mode(
  design = design_bb,
  pfamily_list = pf_bb,
  family = binomial(),
  beta_set = beta_set_bb,
  n = 800L,
  mc_seed = 42L,
  icm_init = TRUE,
  icm_screen = FALSE,
  mc_stable = 1L,
  maxit = 20L,
  verbose = TRUE
)

cat(
  "\nDone. route =", mode$em_route,
  " EM iters =", mode$iterations,
  " converged =", mode$converged, "\n"
)
