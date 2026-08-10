## Scratch: verify the Chapter-B08 worked-failure example.
library(lmebayesCore)

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)

bad <- dat[!(dat$Subject == "308" & dat$Days > 0), ]
bad$Subject <- droplevels(bad$Subject)
cat("rows for 308:", sum(bad$Subject == "308"), "\n")

d <- suppressWarnings(
  model_setup(Reaction ~ Days_c + (1 + Days_c || Subject), data = bad)
)
cat("---- groupef.rank[1:4] ----\n"); print(d$groupef.rank[1:4])
cat("---- table(groupef.rank) ----\n"); print(table(d$groupef.rank))
cat("---- groupef.estimable[1:4] ----\n"); print(d$groupef.estimable[1:4])
cat("---- popef.rank ----\n"); print(d$popef.rank)
cat("---- popef.rank_ok ----\n"); print(d$popef.rank_ok)
cat("---- print(d) ----\n"); print(d)
