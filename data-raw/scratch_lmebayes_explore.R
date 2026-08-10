## Scratch: inspect lmerb() output structures for the lmebayes vignettes.
suppressMessages(library(lmebayes))

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

design <- model_setup(form, data = dat)
cat("=== model_setup ===\n"); print(design)

ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
cat("\n=== Prior_Setup_GLMM ===\n"); print(ps)
cat("\n=== names(ps) ===\n"); print(names(ps))

pf <- pfamily_list(ps)
cat("\n=== pfamily_list ===\n"); print(pf)

set.seed(1)
fit <- lmerb(form, data = dat, pfamily_list = pf,
             dispersion_ranef = ps$group.dispersion,
             n = 200L, progbar = FALSE)

cat("\n=== class ===\n"); print(class(fit))
cat("\n=== names(fit) ===\n"); print(names(fit))
cat("\n=== print(fit) ===\n"); print(fit)
cat("\n=== summary(fit) ===\n"); print(summary(fit))
cat("\n=== fixef ===\n"); print(fixef(fit))
cat("\n=== ranef ===\n"); print(ranef(fit))
cat("\n=== coef ===\n"); print(head(coef(fit)))
cat("\n=== summary_sigma2 ===\n"); print(summary_sigma2(fit))
cat("\n=== m_convergence / convergence ===\n")
print(fit$m_convergence)
print(fit$convergence)
print(fit$convergence_info)
cat("\n=== sweep_history ===\n"); print(is.null(fit$sweep_history))
