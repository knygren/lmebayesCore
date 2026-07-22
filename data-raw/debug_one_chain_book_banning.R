pkgload::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
pkgload::load_all("c:/Rpackages/lmebayes", quiet = TRUE)
library(bayesrules)

dat <- book_banning[, c("state", "removed", "violent")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)
form <- removed_i ~ violent_i + (1 + violent_i || state)
design <- model_setup(form, dat, binomial(), fit_mer = FALSE)
ps <- Prior_Setup_lmebayes(form, dat, binomial(), pwt = 0.01)
prior <- lmebayesCore::priors_from_pfamily_list(
  pfamily_list(ps), ps$dispersion_ranef, design, binomial(), "glmerb"
)
block1 <- lmebayes:::.lmebayes_block1_prior_list(prior)
re_names <- design$re_coef_names
group_levels <- levels(design$groups)
pm <- glmerb_posterior_mode(design, binomial(), prior)

one_chain <- function() {
  set.seed(42L)
  sample.int(.Machine$integer.max - 1L, 1L)
  glmbayesCore::rGLMM_sweep(
    n_chains = 1L, start_fixef = pm$fixef, inner_sweeps = 19L,
    design = design, block1_prior = block1, pfamily_list = prior$pfamily_list,
    family = binomial(), re_names = re_names, group_levels = group_levels,
    b_start = pm$b_mean, progbar = FALSE
  )$fixef_draws
}

one_cpp <- function() {
  set.seed(42L)
  sample.int(.Machine$integer.max - 1L, 1L)
  glmbayesCore::two_block_rNormal_reg_v5(
    n = 1L, y = design$y, x = design$Z, block = design$groups,
    x_hyper = design$X_hyper, prior_list_block1 = block1,
    pfamily_list = prior$pfamily_list, fixef_start = pm$fixef,
    re_coef_names = re_names, group_levels = group_levels,
    group_name = design$group_name, family = binomial(), m_convergence = 19L,
    use_parallel = FALSE, seed = NULL, progbar = FALSE
  )$fixef_draws
}

r <- unlist(one_chain())
c <- unlist(one_cpp())
cat("R:  ", r, "\n")
cat("Cpp:", c, "\n")
cat("Diff:", c - r, "\n")
