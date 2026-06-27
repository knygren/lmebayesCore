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
ps <- Prior_Setup_lmebayes(form, dat, binomial(), pwt = 0.01)

set.seed(42L)
fit <- glmerb(
  form, data = dat, family = binomial(), pfamily_list = pfamily_list(ps),
  n = 500L, mode_gap_max = 1.0, progbar = FALSE, verbose = FALSE
)

init <- unlist(fit$fixef.init)
mode <- unlist(fit$coef.mode)
cat("draw_engine:", fit$convergence$draw_engine, "\n")
cat("m_pilot:", fit$convergence$m_convergence_pilot,
    "m_main:", fit$convergence$m_convergence,
    "n_pilot:", fit$convergence$n_pilot, "\n")
cat("ICM mode:  ", paste(round(mode, 4), collapse = ", "), "\n")
cat("fixef.init:", paste(round(init, 4), collapse = ", "), "\n")
cat("target ~:  -1.11, 0.39\n")
