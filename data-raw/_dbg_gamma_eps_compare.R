devtools::load_all(".", quiet = TRUE)
data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) & !is.na(invalid_ppvt) & invalid_ppvt == 0L &
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
design0 <- model_setup(form, data = dat)
dat <- subset(dat, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
dat$school_id <- droplevels(dat$school_id)
design <- model_setup(form, data = dat)
ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
pf <- pfamily_list(ps)
dispprior_list <- list(dispersion = ps$group.dispersion)

.priors_from_pfamily_list <- getFromNamespace(".priors_from_pfamily_list", "lmebayesCore")
prior_pack <- .priors_from_pfamily_list(pf, ps$group.dispersion, design, gaussian(), "cmp")
mpl <- list(
  group.Sigma = prior_pack$group.Sigma,
  group.dispersion = prior_pack$group.dispersion,
  pop.prior_list = prior_pack$pop.prior_list
)

mode_c05 <- population_mode(
  design, pf, gaussian(), dispprior_list,
  estep = "exact", icm_init = TRUE
)
ref <- lmerb_posterior_mean(design, mpl)
pt <- rlmerb(
  design = design,
  pfamily_list = pf,
  dispprior_list = dispprior_list,
  simulate = FALSE,
  verbose = FALSE,
  print_icm_table = FALSE
)

flatten_fixef <- function(fixef) {
  unlist(lapply(names(fixef), function(k) {
    setNames(fixef[[k]], paste(k, names(fixef[[k]]), sep = " | "))
  }))
}

v_mode <- flatten_fixef(mode_c05$fixef)
v_ref <- flatten_fixef(ref$fixef)
v_rlmerb <- flatten_fixef(pt$popef.mode)

cat("=== Stacked gamma (fixef) vectors ===\n")
cmp <- data.frame(
  mode = v_mode,
  lmerb = v_ref,
  rlmerb = v_rlmerb,
  gap_mode_lmerb = v_mode - v_ref,
  gap_mode_rlmerb = v_mode - v_rlmerb
)
print(cmp, digits = 8)

cat("\nmax |mode - lmerb|:", max(abs(cmp$gap_mode_lmerb)), "\n")
cat("max |mode - rlmerb|:", max(abs(cmp$gap_mode_rlmerb)), "\n")
cat("max |lmerb - rlmerb|:", max(abs(v_ref - v_rlmerb)), "\n")

cat("\ngamma_star length:", length(mode_c05$gamma_star), " q:", mode_c05$q, "\n")
cat("max |gamma_star - stack(mode$fixef)|:",
    max(abs(mode_c05$gamma_star - as.numeric(v_mode[names(v_mode)]))), "\n")

cat("\n=== epsilon* diagnostics at mode ===\n")
cat("kappa (EM rate):", signif(mode_c05$kappa, 5), "\n")
cat("eps_star_closure:", signif(mode_c05$eps_star_closure, 6), "\n")
cat("stationarity:", signif(mode_c05$stationarity, 6), "\n")
eps_opt <- epsilon_optimize(mode_c05)
cat("eps_star optimize:", signif(eps_opt$eps_star, 6),
    " gamma_prime at opt (first 3):",
    paste(signif(eps_opt$gamma_prime[1:3], 5), collapse = ", "), "\n")
cat("gamma_star (first 3):",
    paste(signif(mode_c05$gamma_star[1:3], 5), collapse = ", "), "\n")
cat("|gamma_prime* - gamma_star| at optimizer arg:",
    signif(sqrt(sum((eps_opt$gamma_prime - mode_c05$gamma_star)^2)), 6), "\n")
cat("g_opt (log eps*):", signif(eps_opt$g_opt, 5), "  exp(g_opt):", signif(exp(eps_opt$g_opt), 6), "\n")

# two_block_rate comparison if available
if (exists("two_block_rate", mode = "function")) {
  rate <- tryCatch(two_block_rate(design, mpl), error = function(e) NULL)
  if (!is.null(rate)) {
    cat("\ntwo_block_rate lambda*:", signif(rate$lambda_star, 5), "\n")
  }
}
