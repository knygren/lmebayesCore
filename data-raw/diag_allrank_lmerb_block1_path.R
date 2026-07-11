## Reproduce lmerb Block1 ING path on all-rank (mu_all 2 x 39).
devtools::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
setwd("c:/Rpackages/lmebayes")
source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

fx <- .prepare_small5_all_full_rank_manual()
ps <- Prior_Setup_lmebayes(
  fx$form, data = fx$dat, pwt = 0.01, pwt_measurement = 0.49
)
pf <- pfamily_list(ps)
design <- ps$design
ing <- pf$ing_prior_list

# Minimal batch like pilot chain 1
re_names <- design$re_coef_names
group_levels <- levels(design$groups)
p_re <- length(re_names)
n_chains <- 4L

fixef <- lapply(re_names, function(r) {
  setNames(rep(0, ncol(design$X_hyper[[r]])), colnames(design$X_hyper[[r]]))
})
names(fixef) <- re_names
tau2 <- matrix(1, n_chains, p_re, dimnames = list(NULL, re_names))
batch <- list(
  n = n_chains,
  fixef = fixef,
  tau2 = tau2,
  re_names = re_names,
  group_levels = group_levels,
  b = matrix(0, n_chains, p_re * length(group_levels)),
  iters_ranef = numeric(n_chains)
)

for (i in seq_len(n_chains)) {
  cat("\n=== chain", i, "===\n")
  prep <- glmbayesCore:::.two_block_block1_prep_one_chain(
    batch = batch, i = i, design = design,
    block1_prior = pf$block1_prior, ptypes = pf$ptypes,
    use_cpp_mu_all = FALSE, use_cpp_prior_tau2 = FALSE
  )
  prior_list <- glmbayesCore:::.two_block_block1_ing_prior_list_one_chain(
    prep$mu_all, ing
  )
  cat("mu_all dim:", paste(dim(prior_list$mu), collapse = "x"), "\n")
  tryCatch({
    out <- glmbayesCore:::.two_block_block1_envelope_draw_one_chain(
      y = design$y,
      Z = as.matrix(design$Z),
      groups = design$groups,
      prior_list = prior_list,
      p_re = p_re,
      re_names = re_names,
      group_levels = group_levels
    )
    cat("dispersion_ranef length:", length(out$dispersion_ranef), "\n")
    cat("iters_mean:", out$iters_mean, "\n")
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
}
