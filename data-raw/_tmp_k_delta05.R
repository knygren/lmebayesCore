devtools::load_all("c:/Rpackages/lmebayesCore", quiet = TRUE)

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
form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 + free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)
design0 <- model_setup(form_lmer, data = dat_g)
dat_g <- subset(dat_g, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
dat_g$school_id <- droplevels(dat_g$school_id)
design <- model_setup(form_lmer, data = dat_g)
ps <- Prior_Setup_GLMM(form_lmer, data = dat_g, pop.pwt = 0.01)
pf <- pfamily_list(ps)
dispprior_list <- list(dispersion = ps$group.dispersion)

run_case <- function(delta_2, inner_tol, label) {
  cert <- gamma_beta_tv_certificate(
    design = design,
    pfamily_list = pf,
    family = gaussian(),
    delta_2 = delta_2,
    dispprior_list = dispprior_list,
    tol = inner_tol,
    kappa_method = "none",
    estep = "exact"
  )
  ros <- cert$rosenthal
  full <- cert$full_pi_gamma$full_bound
  cat(sprintf("\n=== %s ===\n", label))
  cat("delta_2:", delta_2, "  r_Gauss:", signif(cert$beta_set$level$r_gauss, 5), "\n")
  cat("kappa_max^LB:", signif(cert$floor_spectrum$kappa_max_lb, 5), "\n")
  cat("k sweeps:", cert$sweeps$k, "\n")
  cat("inner bound:", signif(cert$sweeps$bound, 6), "\n")
  cat("full bound (inner + delta_2):", signif(full, 6), "\n")
  cat("alpha:", signif(ros$alpha, 6),
      "  minorization:", signif(ros$minorization, 4),
      "  drift:", signif(ros$drift, 6), "\n")
  invisible(cert)
}

# Literal reading: delta_2 = 0.05, inner = 0.05 => combined = 0.10
run_case(0.05, 0.05, "delta_2=0.05, inner target 0.05 (combined = 0.10)")

# Combined = 0.01 requires inner = 0.01 - delta_2
cat("\n--- Arithmetic check: full = inner + delta_2 ---\n")
cat("If delta_2 = 0.05 and full = 0.01, inner would need to be", 0.01 - 0.05, "\n")
cat("If delta_2 = 0.005 and full = 0.01, inner target =", 0.01 - 0.005, "\n\n")

# Likely intent: combined = 0.01 with delta_2 = 0.005, inner = 0.005
run_case(0.005, 0.005, "delta_2=0.005, inner target 0.005 (combined = 0.01)")

# Alternative: delta_2 = 0.05 but only ask inner small enough that combined <= 0.01?
# impossible unless inner negative

# Maybe they meant delta_2=0.05 and want combined 0.05? inner=0?
# Or delta_2=0.05, inner=0, full=0.05 - trivial

# Scenario: delta_2=0.05, inner target chosen so full=0.01 -> impossible
# Scenario: delta_2=0.05, inner as small as possible at some k - just report at k for inner=0.05
