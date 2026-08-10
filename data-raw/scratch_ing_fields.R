suppressMessages(library(lmebayes))
data(sleepstudy, package = "lme4")
d <- sleepstudy
d$Days_c <- d$Days - mean(d$Days)
ps <- Prior_Setup_GLMM(Reaction ~ Days_c + (1 + Days_c || Subject),
                       data = d, pop.pwt = 0.01, pop.dispersion.pwt = 0.2)
cat("class(pop.ing_prior):", class(ps$pop.ing_prior), "\n")
str(ps$pop.ing_prior)
cat("\nnames(pop.prior_list[[1]]):\n")
print(names(ps$pop.prior_list[[1]]))
cat("\ndispersion per component:\n")
print(vapply(ps$pop.prior_list, function(p) p$dispersion, numeric(1)))
