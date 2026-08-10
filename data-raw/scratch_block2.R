## Scratch: the row-block interface (Prior_SetupGroup + lmbBlock).
suppressMessages(library(lmebayes))

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)

psb <- Prior_SetupGroup(Reaction ~ Days_c, group = dat$Subject, data = dat,
                        pwt = 0.01)
cat("=== class / names ===\n")
print(class(psb)); print(names(psb))
cat("\n=== print ===\n"); print(psb)

pfb <- pfamily_list(psb)
cat("\n=== pfamily_list length ===\n"); print(length(pfb)); print(names(pfb)[1:3])

t2 <- system.time({
  fb <- lmbBlock(Reaction ~ Days_c, block = dat$Subject, data = dat,
                 pfamily_list = pfb, n = 50L)
})
cat("\nlmbBlock n=50:", t2[["elapsed"]], "s\n")
cat("class:", class(fb), "\n")
cat("names:\n"); print(names(fb))
cat("\n=== print(fb) ===\n"); print(fb)
cat("\n=== summary(fb) ===\n"); print(summary(fb))
