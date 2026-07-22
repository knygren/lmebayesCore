pkgload::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
pkgload::load_all("c:/Rpackages/lmebayes", quiet = TRUE)
library(bayesrules)

dat <- book_banning[, c("state", "removed", "violent")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)
form_glmerb <- removed_i ~ violent_i + (1 + violent_i || state)
design <- model_setup(form_glmerb, dat, binomial(), fit_mer = FALSE)
ps <- Prior_Setup_lmebayes(form_glmerb, dat, binomial(), pwt = 0.01)
prior <- lmebayesCore::priors_from_pfamily_list(
  pfamily_list(ps), ps$dispersion_ranef, design, binomial(), "glmerb"
)
block1_prior <- lmebayes:::.lmebayes_block1_prior_list(prior)
re_names <- design$re_coef_names
group_levels <- levels(design$groups)
pm <- glmerb_posterior_mode(design, binomial(), prior)
fixef_start <- pm$fixef
b_mode <- pm$b_mean
m_pilot <- 19L
n_pilot <- 20L

run_pilot <- function(engine) {
  if (engine == "R") {
    set.seed(42L)
    sample.int(.Machine$integer.max - 1L, 1L)
    out <- glmbayesCore::rGLMM_sweep(
      n_chains = n_pilot, start_fixef = fixef_start, inner_sweeps = m_pilot,
      design = design, block1_prior = block1_prior,
      pfamily_list = prior$pfamily_list, family = binomial(),
      re_names = re_names, group_levels = group_levels,
      b_start = b_mode, progbar = FALSE, stage_label = "pilot"
    )
  } else {
    set.seed(42L)
    sample.int(.Machine$integer.max - 1L, 1L)
    out <- glmbayesCore::two_block_rNormal_reg_v5(
      n = n_pilot, y = design$y, x = design$Z, block = design$groups,
      x_hyper = design$X_hyper, prior_list_block1 = block1_prior,
      pfamily_list = prior$pfamily_list, fixef_start = fixef_start,
      re_coef_names = re_names, group_levels = group_levels,
      group_name = design$group_name, family = binomial(),
      m_convergence = m_pilot, use_parallel = FALSE, seed = NULL,
      progbar = FALSE, stage_label = "pilot"
    )
  }
  glmbayesCore:::.two_block_fixef_colmeans(out$fixef_draws, re_names, fixef_start)
}

cm_r <- run_pilot("R")
cm_cpp <- run_pilot("cpp")

cat("ICM mode:", paste(round(unlist(fixef_start), 4), collapse = ", "), "\n")
cat("R pilot colMeans:", paste(round(unlist(cm_r), 4), collapse = ", "), "\n")
cat("Cpp pilot colMeans:", paste(round(unlist(cm_cpp), 4), collapse = ", "), "\n")
cat("Diff:", paste(round(unlist(cm_cpp) - unlist(cm_r), 4), collapse = ", "), "\n")
