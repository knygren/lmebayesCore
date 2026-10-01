## Quick check: population_mode() em_trace + verbose on book_banning.
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("bayesrules required.", call. = FALSE)
}
devtools::load_all("c:/Rpackages/lmebayesCore", quiet = TRUE)

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
  "em_route=", mode$em_route,
  " iters=", mode$iterations,
  " converged=", mode$converged, "\n"
)
stopifnot(!is.null(mode$em_trace))
stopifnot(!is.null(mode$em_trace$em_start))
stopifnot(is.finite(mode$em_trace$em_start$log_post))
stopifnot(length(mode$em_trace$history) >= 1L)
cat("population_mode em_trace check passed.\n")
