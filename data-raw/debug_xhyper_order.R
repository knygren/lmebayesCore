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
group_levels <- levels(design$groups)
re_names <- design$re_coef_names

for (k in re_names) {
  X_k <- as.matrix(design$X_hyper[[k]])
  rn <- rownames(X_k)
  cat("\nComponent:", k, "\n")
  cat("  nrow:", nrow(X_k), " group_levels:", length(group_levels), "\n")
  if (!is.null(rn)) {
    cat("  X row order matches group_levels:", identical(rn, group_levels), "\n")
    if (!identical(rn, group_levels)) {
      cat("  first 5 X rn:", paste(head(rn, 5), collapse = ", "), "\n")
      cat("  first 5 grp:", paste(head(group_levels, 5), collapse = ", "), "\n")
    }
  } else {
    cat("  X has NO rownames\n")
  }
}
