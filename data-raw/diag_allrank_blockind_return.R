devtools::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
setwd("c:/Rpackages/lmebayes")
source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

fx <- .prepare_small5_all_full_rank_manual()
ps <- Prior_Setup_lmebayes(
  fx$form, data = fx$dat, pwt = 0.01, pwt_measurement = 0.49
)
m <- ps$ing_prior_measurement
design <- ps$design
Z <- as.matrix(design$Z)
y <- design$y
groups <- design$groups
re_names <- design$re_coef_names
group_levels <- levels(groups)
p_re <- length(re_names)

prior_list <- list(
  shape = m$shape, rate = m$rate, max_disp_perc = 0.99,
  disp_lower = m$disp_lower, disp_upper = m$disp_upper,
  mu = matrix(0, p_re, 1, dimnames = list(re_names, NULL)),
  P = diag(p_re), Inv_Dispersion = TRUE
)

raw <- glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp(
  n = 1L, y = y, x = Z, block = groups, prior_list = prior_list,
  offset = rep(0, length(y)), wt = rep(1, length(y)), p_re = p_re,
  use_parallel = TRUE, verbose = FALSE,
  group_levels = group_levels, re_names = re_names
)

cat("names(raw):\n")
print(names(raw))
cat("\ndispersion_ranef:\n")
print(raw$dispersion_ranef)
cat("length:", length(raw$dispersion_ranef), "\n")
cat("\ndisp_out:\n")
print(raw$disp_out)
cat("\niters_mean:\n")
print(raw$iters_mean)
cat("\nsubset check:\n")
sub <- raw[c("b", "dispersion_ranef", "iters_mean")]
str(sub, max.level = 1)
