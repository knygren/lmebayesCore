pkgload::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
pkgload::load_all("c:/Rpackages/lmebayes", quiet = TRUE)
library(bayesrules)

dat <- book_banning[, c("state", "removed", "violent")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)
ps <- Prior_Setup_lmebayes(
  removed_i ~ violent_i + (1 + violent_i || state), dat, binomial(), pwt = 0.01
)
pf <- pfamily_list(ps)
for (k in names(pf)) {
  S <- pf[[k]]$prior_list$Sigma
  R <- chol(S)
  P_r <- 0.5 * (chol2inv(R) + t(chol2inv(R)))
  P_cpp <- 0.5 * (solve(S) + t(solve(S)))
  cat(k, "max|P diff|", max(abs(P_r - P_cpp)), "\n")
}
