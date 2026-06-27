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
prior <- lmebayes:::.lmebayes_priors_from_pfamily_list(
  pfamily_list(ps), ps$dispersion_ranef, design, binomial(), "glmerb"
)
block1 <- lmebayes:::.lmebayes_block1_prior_list(prior)
group_levels <- levels(design$groups)
pm <- glmerb_posterior_mode(design, binomial(), prior)
fixef <- pm$fixef

ptypes <- vapply(prior$pfamily_list, function(pf) pf$pfamily, character(1))
tau2 <- glmbayesCore:::.two_block_tau2_start_from_pfamily(
  prior$pfamily_list, design$re_coef_names
)
mu_all <- as.matrix(build_mu_all(design, fixef, group_levels)$mu_all)
prior_list <- glmbayesCore:::.two_block_block1_prior_with_tau2(
  block1, tau2, ptypes, design$re_coef_names, mu_all
)

set.seed(42L)
ss <- sample.int(.Machine$integer.max - 1L, 1L)
set.seed(ss + 1L)

b1_r <- block_rNormalGLM(
  n = 1L, y = design$y, x = design$Z, block = design$groups,
  prior_list = prior_list, family = binomial(),
  use_parallel = FALSE, verbose = FALSE, progbar = FALSE
)
b_r <- b1_r$coefficients
rn <- rownames(b_r)
ord <- match(group_levels, rn)
b_r <- b_r[ord, , drop = FALSE]

set.seed(ss + 1L)
b1_cpp <- glmbayesCore:::.block_rNormalGLM_cpp(
  n = 1L, y = design$y, x = design$Z, block = design$groups,
  prior_list = prior_list, prior_lists = NULL,
  offset = rep(0, length(design$y)), wt = rep(1, length(design$y)),
  f2 = glmbfamfunc(binomial())$f2, f3 = glmbfamfunc(binomial())$f3,
  family = "binomial", link = "logit", Gridtype = 2L, n_envopt = 1L,
  use_parallel = FALSE, use_opencl = FALSE, verbose = FALSE
)
ids <- b1_cpp$block_info$ids
b_cpp <- b1_cpp$coefficients
ord2 <- match(group_levels, ids)
b_cpp <- b_cpp[ord2, , drop = FALSE]

cat("Block1 draw diff max:", max(abs(b_r - b_cpp)), "\n")
cat("ids match group_levels:", identical(as.character(ids[ord2]), group_levels), "\n")

# Block 2 component 1
k <- design$re_coef_names[1]
X_k <- as.matrix(design$X_hyper[[k]])
y_k <- glmbayesCore:::.two_block_align_b_to_xhyper(b_r[, k], X_k, group_levels)
pf <- prior$pfamily_list[[k]]
set.seed(999L)
fit_r <- rglmb(1L, y_k, X_k, gaussian(), pf, use_parallel = FALSE, verbose = FALSE)
set.seed(999L)
pl <- pf$prior_list
R <- chol(pl$Sigma)
P <- 0.5 * (chol2inv(R) + t(chol2inv(R)))
out_cpp <- glmbayesCore:::.rNormalReg_cpp(
  1L, y_k, X_k, pl$mu, P, rep(0, nrow(X_k)), rep(1, nrow(X_k)),
  pl$dispersion, glmbfamfunc(gaussian())$f2, glmbfamfunc(gaussian())$f3,
  pl$mu, "gaussian", "identity", 2L
)
cat("Block2 coef.mode diff:", max(abs(fit_r$coef.mode - out_cpp$coef.mode)), "\n")
