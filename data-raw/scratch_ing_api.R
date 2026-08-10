## Scratch: current API for ING Block 2 and per-group dispersion.
suppressMessages(library(lmebayes))

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01,
                       pop.dispersion.pwt = 0.2)
cat("=== names(ps) ===\n"); print(names(ps))
cat("\n=== ps$pop.ing_prior ===\n"); print(ps$pop.ing_prior)
cat("\n=== ps$group.ing_prior ===\n"); print(ps$group.ing_prior)
cat("\n=== ps$group.dispersion ===\n"); print(ps$group.dispersion)

pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
cat("\n=== pfamily_list ING ===\n"); print(pf_ing)

t1 <- system.time({
  fit_ing <- lmerb(form, data = dat, pfamily_list = pf_ing,
                   dispersion_ranef = ps$group.dispersion,
                   n = 100L, progbar = FALSE)
})
cat("\nlmerb ING n=100:", t1[["elapsed"]], "s\n")
cat("names(fit_ing):\n"); print(names(fit_ing))
cat("\npilot:\n"); print(fit_ing$pilot)
cat("\nconvergence:\n"); print(fit_ing$convergence)
cat("\nsweep_history names:\n"); print(names(fit_ing$sweep_history))
cat("\npopef.dispersion head:\n"); print(utils::head(fit_ing$popef.dispersion))
cat("\nsummary:\n"); print(summary(fit_ing))

## Per-group observation dispersion
cat("\n=== dGamma_list ===\n")
print(args(dGamma_list))
