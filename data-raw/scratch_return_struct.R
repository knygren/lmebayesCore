## Scratch: dump the returned structure of rlmerb()/rglmerb() so the @return
## sections can be rewritten against what the functions actually return.
library(lmebayesCore)

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

design <- model_setup(form, data = dat)
ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)

fit <- rlmerb(
  n = 50L, design = design,
  pfamily_list   = pfamily_list(ps),
  dispprior_list = list(dispersion = ps$group.dispersion),
  progbar = FALSE, verbose = FALSE, print_icm_table = FALSE
)

cat("==== class ====\n"); print(class(fit))
cat("\n==== names(fit) ====\n"); print(names(fit))

cat("\n==== shapes ====\n")
for (nm in names(fit)) {
  v <- fit[[nm]]
  d <- if (is.null(v)) "NULL"
       else if (is.function(v)) "function"
       else if (is.list(v)) paste0("list[", length(v), "]: ",
                                   paste(utils::head(names(v), 8), collapse = ", "))
       else if (!is.null(dim(v))) paste0(class(v)[1], " ", paste(dim(v), collapse = "x"))
       else paste0(class(v)[1], " len=", length(v))
  cat(sprintf("  %-28s %s\n", nm, d))
}

cat("\n==== popef ====\n"); str(fit$popef, max.level = 2)
cat("\n==== groupef ====\n"); str(fit$groupef, max.level = 1)
cat("\n==== convergence_info ====\n"); str(fit$convergence_info, max.level = 2)
cat("\n==== m_convergence ====\n"); print(fit$m_convergence)

cat("\n\n######## GLMM ########\n")
ps_p <- Prior_Setup_GLMM(form, data = dat, family = poisson(), pop.pwt = 0.01)
