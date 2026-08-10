## Scratch: verify the Chapter-B06 extended-rate example runs as written.
library(lmebayesCore)

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)
form <- Reaction ~ Days_c + (1 + Days_c || Subject)

design <- model_setup(form, data = dat)
ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01,
                       pop.dispersion.pwt = 0.2)

cat("---- names(ps) ----\n"); print(names(ps))
cat("---- ps$pop.ing_prior ----\n"); print(ps$pop.ing_prior)
cat("---- ps$pop.prior_list ----\n"); print(ps$pop.prior_list)
cat("---- design$groupef.names / group_name ----\n")
print(design$groupef.names); print(design$group_name)

pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
cat("---- pf_ing ----\n"); print(pf_ing)

cat("---- group.Sigma / group.dispersion ----\n")
print(ps$group.Sigma); print(ps$group.dispersion)

cat("---- base certified rate ----\n")
rate_base <- two_block_rate_from_pfamily_list(
  x = design$D, group = design$group, x_hyper = design$W,
  prior_list_block1 = list(Sigma = ps$group.Sigma,
                           dispersion = ps$group.dispersion),
  pfamily_list = pf_ing
)
print(rate_base$lambda_star)

cat("---- reference u ----\n")
u_ref <- lmebayesCore:::.lmebayes_reference_u(
  design$lmer, design$group_name, design$groupef.names
)
print(utils::head(u_ref))

cat("---- build lambda_ing ----\n")
lambda_ing <- stats::setNames(lapply(names(ps$pop.ing_prior), function(k) {
  g <- ps$pop.ing_prior[[k]]
  list(lambda = 1 / g$disp_lower,
       shape  = g$shape,
       u      = stats::setNames(u_ref[[k]], rownames(u_ref)))
}), names(ps$pop.ing_prior))
str(lambda_ing, max.level = 2)

cat("---- extended rate ----\n")
rate_ext <- lmebayesCore:::two_block_rate_ing(
  x = design$D, group = design$group, x_hyper = design$W,
  prior_list_block1 = list(Sigma = ps$group.Sigma,
                           dispersion = ps$group.dispersion),
  prior_list_block2 = ps$pop.prior_list,
  lambda_ing = lambda_ing
)
print(rate_ext)
cat("---- fields ----\n")
print(c(base = rate_ext$lambda_star_base, extended = rate_ext$lambda_star))
print(rate_ext$ext_names)
