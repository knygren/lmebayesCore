## Scratch: B06 -- prior-mean plug-in vs disp_lower corner, and the u = 0 floor.
library(lmebayesCore)

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

design <- model_setup(form, data = dat)
ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01,
                       pop.dispersion.pwt = 0.2)
pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
pl1 <- list(Sigma = ps$group.Sigma, dispersion = ps$group.dispersion)

u_ref <- lmebayesCore:::.lmebayes_reference_u(
  design$lmer, design$group_name, design$groupef.names
)

make_lambda_ing <- function(lambda_fun, u_scale = 1) {
  stats::setNames(lapply(names(ps$pop.ing_prior), function(k) {
    g <- ps$pop.ing_prior[[k]]
    list(lambda = lambda_fun(g),
         shape  = g$shape,
         u      = stats::setNames(u_scale * u_ref[[k]], rownames(u_ref)))
  }), names(ps$pop.ing_prior))
}

run <- function(li) {
  r <- lmebayesCore:::two_block_rate_ing(
    x = design$D, group = design$group, x_hyper = design$W,
    prior_list_block1 = pl1, prior_list_block2 = ps$pop.prior_list,
    lambda_ing = li, warn_slow = FALSE
  )
  c(base = r$lambda_star_base, ext = r$lambda_star)
}

cat("---- prior-mean plug-in lambda = (shape-1)/rate ----\n")
print(run(make_lambda_ing(function(g) (g$shape - 1) / g$rate)))

cat("---- disp_lower corner lambda = 1/disp_lower ----\n")
print(run(make_lambda_ing(function(g) 1 / g$disp_lower)))

cat("---- prior-mean lambda, u = 0 (should equal base) ----\n")
print(run(make_lambda_ing(function(g) (g$shape - 1) / g$rate, u_scale = 0)))

cat("---- prior-mean lambda, u scaled 0, .25, .5, 1, 2 ----\n")
for (s in c(0, 0.25, 0.5, 1, 2)) {
  r <- run(make_lambda_ing(function(g) (g$shape - 1) / g$rate, u_scale = s))
  cat(sprintf("  u_scale = %-5s base = %.6f  ext = %.6f\n", s, r[["base"]], r[["ext"]]))
}

cat("---- tau2 plug-ins for reference ----\n")
print(vapply(ps$pop.ing_prior, function(g) g$rate / (g$shape - 1), numeric(1)))
print(vapply(ps$pop.ing_prior, function(g) g$disp_lower, numeric(1)))
