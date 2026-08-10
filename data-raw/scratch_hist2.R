suppressMessages(library(lmebayesCore))
data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)
design <- model_setup(form, data = dat)
ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01,
                       pop.dispersion.pwt = 0.2)
pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

set.seed(1)
fit <- rlmerb(n = 100L, design = design,
              pfamily_list = pf_ing,
              dispprior_list = list(dispersion = ps$group.dispersion),
              progbar = FALSE, verbose = FALSE, print_icm_table = FALSE)

cat("class(sweep_history):\n"); print(class(fit$sweep_history))
cat("names(sweep_history):\n"); print(names(fit$sweep_history))
cat("stage:", fit$sweep_history$stage, "\n")
cat("\nany pilot fields in fit?\n")
print(grep("pilot", names(fit), value = TRUE, ignore.case = TRUE))
cat("\nconvergence_info:\n"); str(fit$convergence_info, max.level = 1)
cat("\nsweep_history print:\n"); print(fit$sweep_history)
