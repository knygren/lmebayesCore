## Scratch: convergence objects, plots, and simulate = FALSE.
suppressMessages(library(lmebayes))

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)
ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)

invisible(utils::capture.output(
  fit <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
               dispersion_ranef = ps$group.dispersion, n = 200L,
               progbar = FALSE, sim_method = "TWO_BLOCK_GIBBS")
))

cat("=== convergence ===\n"); print(fit$convergence)
cat("\n=== m_convergence ===\n"); print(fit$m_convergence)
cat("\n=== sweep_history names ===\n"); print(names(fit$sweep_history))
cat("\n=== sweep_history$main ===\n"); print(fit$sweep_history$main)
cat("\n=== class(sweep_history$main) ===\n"); print(class(fit$sweep_history$main))

cat("\n=== plot_mean_convergence ===\n")
p1 <- try(plot_mean_convergence(fit$sweep_history$main), silent = TRUE)
print(class(p1))
cat("\n=== plot_var_convergence ===\n")
p2 <- try(plot_var_convergence(fit$sweep_history$main), silent = TRUE)
print(class(p2))

cat("\n=== tv_tol sensitivity ===\n")
for (tol in c(0.1, 0.01, 1e-4)) {
  invisible(utils::capture.output(
    f <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
               dispersion_ranef = ps$group.dispersion, n = 20L,
               progbar = FALSE, sim_method = "TWO_BLOCK_GIBBS", tv_tol = tol)
  ))
  cat(sprintf("  tv_tol = %-8g  m_convergence = %d\n", tol, f$m_convergence))
}

cat("\n=== simulate = FALSE ===\n")
invisible(utils::capture.output(
  pt <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
              dispersion_ranef = ps$group.dispersion, simulate = FALSE)
))
print(names(pt))
print(pt)
