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
re_names <- design$re_coef_names
group_levels <- levels(design$groups)
pm <- glmerb_posterior_mode(design, binomial(), prior)

cat("re_names:", paste(re_names, collapse = ", "), "\n")
cat("Z colnames:", paste(colnames(design$Z), collapse = ", "), "\n")
cat("fixef names:", paste(names(pm$fixef), collapse = ", "), "\n")
cat("ICM:\n")
print(unlist(pm$fixef))

# Route Block 2 through R rglmb from C++ by running one R-batch sweep manually
# vs checking mu_all orientation: rownames should match re_names
mu <- build_mu_all(design, pm$fixef, group_levels)$mu_all
cat("\nmu_all rownames:", paste(rownames(mu), collapse = ", "), "\n")
cat("mu_all col 1 (first state):", paste(round(mu[, 1], 4), collapse = ", "), "\n")
