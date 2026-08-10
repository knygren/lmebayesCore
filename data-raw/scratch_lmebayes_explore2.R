## Scratch part 2: accessor methods on an lmerb fit.
suppressMessages(library(lmebayes))
suppressMessages(library(lme4))

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

ps  <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
set.seed(1)
fit <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
             dispersion_ranef = ps$group.dispersion,
             n = 200L, progbar = FALSE)

cat("=== fixef ===\n"); print(fixef(fit))
cat("\n=== ranef ===\n"); print(ranef(fit))
cat("\n=== class(ranef) ===\n"); print(class(ranef(fit)))
cat("\n=== summary(ranef) ===\n"); print(summary(ranef(fit)))
cat("\n=== coef ===\n"); print(utils::head(coef(fit)))
cat("\n=== summary_sigma2 ===\n"); print(summary_sigma2(fit))
cat("\n=== m_convergence ===\n"); print(fit$m_convergence)
cat("\n=== convergence ===\n"); print(fit$convergence)
cat("\n=== sim_method_used / draw_engine ===\n")
print(fit$sim_method_used); print(fit$draw_engine)
cat("\n=== head(popef) ===\n"); print(utils::head(fit$popef))
cat("\n=== head(groupef) ===\n"); print(utils::head(fit$groupef))
cat("\n=== sweep_history null? ===\n"); print(is.null(fit$sweep_history))
