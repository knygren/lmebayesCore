suppressMessages(library(lmebayesCore))
data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)
design <- model_setup(form, data = dat)
ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)

set.seed(1)
fit <- rlmerb(n = 50L, design = design,
              pfamily_list = pfamily_list(ps),
              dispprior_list = list(dispersion = ps$group.dispersion),
              sim_method = "TWO_BLOCK_GIBBS",
              progbar = FALSE, verbose = FALSE, print_icm_table = FALSE)

cat("names(fit):\n"); print(names(fit))
cat("\nnames(sweep_history):\n"); print(names(fit$sweep_history))
cat("\nclass(sweep_history):\n"); print(class(fit$sweep_history))
cat("\nclass(sweep_history$main):\n"); print(class(fit$sweep_history$main))
cat("\nstr(sweep_history, 2):\n"); str(fit$sweep_history, max.level = 2)
