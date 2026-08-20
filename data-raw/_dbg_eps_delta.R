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

mode <- population_mode(
  design, pf, gaussian(), dispprior_list,
  estep = "exact", icm_init = TRUE
)
eps_star <- mode$eps_star_closure

for (delta in c(0.01, 0.1)) {
  c <- epsilon(eps_star, delta, mode)
  cat(
    "delta =", delta,
    "  d =", signif(c$d, 5),
    "  exp(-d) =", signif(exp(-c$d), 4),
    "  escape =", signif(c$cal$escape_mass, 4),
    "  eps_d =", signif(c$eps_d, 4),
    "  Q_lb =", signif(c$Q_mass_lb, 4),
    "  eps =", signif(c$eps, 4),
    "\n"
  )
}
