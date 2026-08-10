## Scratch part 3: shape of the stored draw containers on an lmerb fit.
suppressMessages(library(lmebayes))

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

ps  <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
set.seed(1)
fit <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
             dispersion_ranef = ps$group.dispersion,
             n = 50L, progbar = FALSE)

cat("=== str(groupef) ===\n"); print(utils::str(fit$groupef))
cat("\n=== head(groupef) ===\n"); print(utils::head(fit$groupef))
cat("\n=== groupef.mode ===\n"); print(utils::head(fit$groupef.mode))
cat("\n=== popef.mode ===\n"); print(fit$popef.mode)
cat("\n=== popef.means ===\n"); print(fit$popef.means)
cat("\n=== str(popef) ===\n"); print(utils::str(fit$popef))
cat("\n=== model_setup names ===\n"); print(names(fit$model_setup))
cat("\n=== groupef.names / group_name ===\n")
print(fit$model_setup$groupef.names); print(fit$model_setup$group_name)
cat("\n=== summary_sigma2 ===\n"); print(summary_sigma2(fit))
cat("\n=== convergence ===\n"); print(fit$convergence)
cat("\n=== m_convergence ===\n"); print(fit$m_convergence)
cat("\n=== sim_method_used / draw_engine ===\n")
print(fit$sim_method_used); print(fit$draw_engine)
