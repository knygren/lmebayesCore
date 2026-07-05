pkgload::load_all(".", quiet = TRUE)
pkgload::load_all("../lmebayes", quiet = TRUE)
data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)
form <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)
design <- model_setup(form, data = dat)
ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
re <- design$re_coef_names
for (k in re) {
  cat(k, "X_hyper dim", paste(dim(design$X_hyper[[k]]), collapse = "x"), "\n")
}
cat("J groups", nlevels(design$groups), "\n")

fit1 <- lmerb(form, dat, pfamily_list(ps), ps$dispersion_ranef, simulate = FALSE)
fit2 <- lmerb(form, dat, pf_ing, ps$dispersion_ranef, simulate = FALSE)
b <- fit2$ranef.mode
g <- fit2$fixef.mode
b1 <- fit1$ranef.mode
g1 <- fit1$fixef.mode
gl <- levels(design$groups)

for (k in re) {
  pl <- pf_ing[[k]]$prior_list
  Xk <- design$X_hyper[[k]]
  bk <- b[, k]
  gk <- g[[k]]
  bal <- glmbayesCore:::.two_block_align_b_to_xhyper(bk, Xk, gl)
  eta <- as.numeric(Xk %*% gk)
  rss <- sum((bal - eta)^2)
  m <- length(bal)
  shape_post <- pl$shape + m / 2
  rate_post <- pl$rate + rss / 2
  tau2_mode <- rate_post / (shape_post - 1)
  cat("\n", k, ":\n", sep = "")
  cat("  m=", m, " rss=", rss, " sd(b)=", sd(bal), "\n", sep = "")
  cat("  shape=", pl$shape, " rate=", pl$rate, "\n", sep = "")
  cat("  plug-in=", ps$prior_list[[k]]$dispersion_fixef,
      " formula=", tau2_mode, " fit=", fit2$tau2.mode[[k]], "\n", sep = "")
  cat("  case1 sd(b)=", sd(b1[, k]), " case2 sd(b)=", sd(bal), "\n", sep = "")
}

cat("\n--- rglmb Block2 mode at case1 b (Gibbs path) ---\n")
for (k in re) {
  Xk <- design$X_hyper[[k]]
  yk <- glmbayesCore:::.two_block_align_b_to_xhyper(b1[, k], Xk, gl)
  fit_k <- rglmb(
    n = 1L, y = yk, x = Xk, family = stats::gaussian(),
    pfamily = pf_ing[[k]], verbose = FALSE, use_parallel = FALSE
  )
  cat(k, ": rglmb tau2=", fit_k$dispersion[1L],
      " gamma=", paste(round(fit_k$coef.mode, 4), collapse = ","), "\n", sep = "")
}

cat("\n--- rglmb Block2 mode at case2 converged b ---\n")
for (k in re) {
  Xk <- design$X_hyper[[k]]
  yk <- glmbayesCore:::.two_block_align_b_to_xhyper(b[, k], Xk, gl)
  fit_k <- rglmb(
    n = 1L, y = yk, x = Xk, family = stats::gaussian(),
    pfamily = pf_ing[[k]], verbose = FALSE, use_parallel = FALSE
  )
  cat(k, ": rglmb tau2=", fit_k$dispersion[1L], "\n", sep = "")
}

