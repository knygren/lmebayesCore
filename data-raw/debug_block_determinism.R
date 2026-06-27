pkgload::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
pkgload::load_all("c:/Rpackages/lmebayes", quiet = TRUE)
library(bayesrules)

dat <- book_banning[, c("state", "removed", "violent")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)
design <- model_setup(
  removed_i ~ violent_i + (1 + violent_i || state), dat, binomial(), fit_mer = FALSE
)
ps <- Prior_Setup_lmebayes(
  removed_i ~ violent_i + (1 + violent_i || state), dat, binomial(), pwt = 0.01
)
prior <- lmebayes:::.lmebayes_priors_from_pfamily_list(
  pfamily_list(ps), ps$dispersion_ranef, design, binomial(), "glmerb"
)
block1 <- lmebayes:::.lmebayes_block1_prior_list(prior)
pm <- glmerb_posterior_mode(design, binomial(), prior)
ptypes <- vapply(prior$pfamily_list, function(pf) pf$pfamily, character(1))
tau2 <- glmbayesCore:::.two_block_tau2_start_from_pfamily(
  prior$pfamily_list, design$re_coef_names
)
mu_all <- as.matrix(build_mu_all(design, pm$fixef, levels(design$groups))$mu_all)
pl <- glmbayesCore:::.two_block_block1_prior_with_tau2(
  block1, tau2, ptypes, design$re_coef_names, mu_all
)
f <- glmbfamfunc(binomial())
args <- list(
  n = 1L, y = design$y, x = design$Z, block = design$groups,
  prior_list = pl, prior_lists = NULL,
  offset = rep(0, length(design$y)), wt = rep(1, length(design$y)),
  f2 = f$f2, f3 = f$f3, family = "binomial", link = "logit",
  Gridtype = 2L, n_envopt = 1L, use_parallel = FALSE,
  use_opencl = FALSE, verbose = FALSE
)
set.seed(1L)
a <- do.call(glmbayesCore:::.block_rNormalGLM_cpp, args)$coefficients
set.seed(1L)
b <- do.call(glmbayesCore:::.block_rNormalGLM_cpp, args)$coefficients
cat("repeat cpp max diff:", max(abs(a - b)), "\n")

# Compare wrapper vs cpp same seed
set.seed(1L)
w <- block_rNormalGLM(
  n = 1L, y = design$y, x = design$Z, block = design$groups,
  prior_list = pl, family = binomial(), use_parallel = FALSE,
  verbose = FALSE, progbar = FALSE
)$coefficients
set.seed(1L)
c <- do.call(glmbayesCore:::.block_rNormalGLM_cpp, args)$coefficients
cat("wrapper vs cpp max diff:", max(abs(w - c)), "\n")
