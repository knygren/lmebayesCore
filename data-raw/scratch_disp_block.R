## Scratch: per-group dispersion (dGamma_list) and the row-block interface.
suppressMessages(library(lmebayes))

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01,
                       dispformula = ~Subject)

cat("=== dGamma_list(ps) ===\n")
dl <- dGamma_list(ps)
print(class(dl))
print(utils::head(names(dl)))
print(dl[[1L]])

t1 <- system.time({
  fit <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
               dispersion_ranef = dl, dispformula = ~Subject,
               n = 100L, progbar = FALSE)
})
cat("\nlmerb dGamma_list n=100:", t1[["elapsed"]], "s\n")
cat("sim_method:", fit$sim_method_used, " m_conv:", fit$m_convergence, "\n")
cat("\n=== summary_sigma2 ===\n")
print(summary_sigma2(fit))

## --- Row-block interface --------------------------------------------------
cat("\n=== lmbBlock ===\n")
psb <- Prior_SetupGroup(Reaction ~ Days_c, block = dat$Subject, data = dat)
print(names(psb))
t2 <- system.time({
  fb <- lmbBlock(Reaction ~ Days_c, block = dat$Subject, data = dat,
                 pfamily_list = pfamily_list(psb), n = 50L)
})
cat("lmbBlock n=50:", t2[["elapsed"]], "s\n")
print(class(fb))
print(fb)
