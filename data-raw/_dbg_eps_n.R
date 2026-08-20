devtools::load_all(".", quiet = TRUE)

compare_eps <- function(label, mode, ns = c(200L, 10000L), mc_seed = 42L) {
  cat("\n=== ", label, " ===\n", sep = "")
  cat(
    "mode: n =", mode$n,
    " converged:", mode$converged,
    " iters:", mode$iterations,
    " delta:", signif(mode$delta, 4),
    if (is.finite(mode$tol_eff)) paste0(" tol_eff:", signif(mode$tol_eff, 4)) else "",
    "\n"
  )
  eps_cl <- epsilon_star(mode, method = "closure")$eps_star
  cat("closure eps* (from tilde_J at mode n):", signif(eps_cl, 6), "\n")
  out <- lapply(ns, function(n) {
    opt <- epsilon_optimize(mode, n = n, mc_seed = mc_seed)
    c(eps = opt$eps_star, certified = opt$certified, g = opt$g_opt)
  })
  for (i in seq_along(ns)) {
    cat(
      "optimize eps* n =", ns[[i]], ":",
      signif(out[[i]]["eps"], 6),
      " certified:", out[[i]]["certified"],
      " g_opt:", signif(out[[i]]["g"], 5), "\n"
    )
  }
  if (length(ns) == 2L) {
    e200 <- out[[1]]["eps"]
    e10k <- out[[2]]["eps"]
    cat(
      "|optimize_200 - optimize_10k|:", signif(abs(e200 - e10k), 4),
      " rel:", signif(abs(e200 - e10k) / e10k, 4), "\n"
    )
  }
  invisible(list(closure = eps_cl, optimize = out))
}

# --- Gaussian (exact mode; optimize n should not matter) -----------------
data(big_word_club, package = "bayesrules")
dat_g <- big_word_club
dat_g$school_id <- factor(dat_g$school_id)
dat_g <- subset(
  dat_g,
  !is.na(score_ppvt) & !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat_g[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)
form_g <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)
design0 <- model_setup(form_g, data = dat_g)
dat_g <- subset(dat_g, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
dat_g$school_id <- droplevels(dat_g$school_id)
design_g <- model_setup(form_g, data = dat_g)
ps_g <- Prior_Setup_GLMM(form_g, data = dat_g, pop.pwt = 0.01)
pf_g <- pfamily_list(ps_g)

mode_g <- population_mode(
  design_g, pf_g, gaussian(),
  dispprior_list = list(dispersion = ps_g$group.dispersion),
  estep = "exact", icm_init = TRUE
)
compare_eps("Gaussian big_word_club (exact mode)", mode_g)

# --- Binomial (MC mode + MC optimize integral) ---------------------------
data(book_banning, package = "bayesrules")
dat_b <- book_banning[, c("state", "removed", "violent")]
dat_b <- dat_b[stats::complete.cases(dat_b), ]
dat_b$removed_i <- as.integer(dat_b$removed == 1L | dat_b$removed == "1")
dat_b$violent_i <- as.integer(
  dat_b$violent == TRUE | dat_b$violent == 1L | dat_b$violent == "TRUE"
)
dat_b$state <- factor(dat_b$state)
keep <- names(sort(table(dat_b$state), decreasing = TRUE))[seq_len(8L)]
dat_b <- droplevels(subset(dat_b, state %in% keep))
form_b <- removed_i ~ violent_i + (1 + violent_i || state)
design_b <- model_setup(form_b, data = dat_b, family = binomial())
ps_b <- Prior_Setup_GLMM(form_b, data = dat_b, family = binomial(), pop.pwt = 0.01)
pf_b <- pfamily_list(ps_b)

run_binomial <- function(n, label) {
  mode <- population_mode(
    design_b, pf_b, binomial(),
    estep = "mc", n = n, mc_seed = 42L,
    icm_init = TRUE, maxit = 200L
  )
  res <- compare_eps(label, mode)
  res
}

cat("\n\n######## Binomial book_banning: mode built at each n ########\n")
r200 <- run_binomial(200L, "Binomial mode n=200")
r10k <- run_binomial(10000L, "Binomial mode n=10000")

cat(
  "\nCross-mode closure eps* gap:",
  signif(abs(r200$closure - r10k$closure), 4),
  " rel:", signif(abs(r200$closure - r10k$closure) / r10k$closure, 4), "\n"
)

cat("\nSame mode (n=10000), optimize integral at n=200 vs n=10000:\n")
mode10k <- population_mode(
  design_b, pf_b, binomial(),
  estep = "mc", n = 10000L, mc_seed = 42L,
  icm_init = TRUE, maxit = 200L
)
compare_eps("Binomial fixed mode n=10000", mode10k)
