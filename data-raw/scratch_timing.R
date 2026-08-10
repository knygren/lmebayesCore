## Scratch: how long do the candidate vignette fits take?
suppressMessages(library(lmebayes))

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

t0 <- system.time({
  ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
})
cat("Prior_Setup_GLMM:", t0[["elapsed"]], "s\n")

t1 <- system.time({
  fit <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
               dispersion_ranef = ps$group.dispersion,
               n = 200L, progbar = FALSE)
})
cat("lmerb n=200 DEFAULT:", t1[["elapsed"]], "s\n")

t2 <- system.time({
  fit2 <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
                dispersion_ranef = ps$group.dispersion,
                n = 200L, progbar = FALSE, sim_method = "TWO_BLOCK_GIBBS")
})
cat("lmerb n=200 TWO_BLOCK_GIBBS:", t2[["elapsed"]], "s\n")
cat("m_convergence:", fit2$m_convergence, "\n")
cat("sweep_history null?", is.null(fit2$sweep_history), "\n")

## Does lmerb forward print_icm_table through ...?
t3 <- try(lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
                dispersion_ranef = ps$group.dispersion,
                n = 10L, progbar = FALSE, print_icm_table = FALSE),
          silent = TRUE)
cat("print_icm_table forwarded:", !inherits(t3, "try-error"), "\n")
