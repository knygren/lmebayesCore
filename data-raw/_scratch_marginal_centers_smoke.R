devtools::load_all(quiet = TRUE)
library(lme4)
data(sleepstudy, package = "lme4")

ps <- Prior_Setup_GLMM(
  Reaction ~ Days + (Days || Subject),
  data = sleepstudy
)
cat("source=", attr(ps$group.dispersion, "source"), "\n")
cat(
  "sigma2_marg=", ps$group.dispersion,
  " classical=", attr(ps$group.dispersion, "classical"), "\n"
)
cat(
  "ing rate/(shape-1)=",
  ps$group.ing_prior$rate / (ps$group.ing_prior$shape - 1), "\n"
)
for (k in names(ps$pop.prior_list)) {
  d <- ps$pop.prior_list[[k]]$dispersion
  ing <- ps$pop.ing_prior[[k]]
  cat(
    k, " tau_marg=", d,
    " ref=", ps$pop.dispersion.ref[[k]],
    " ing_mean=", ing$rate / (ing$shape - 1), "\n"
  )
}
pf <- pfamily_list(ps)
cat(
  "pfamily ok:",
  paste(vapply(pf, function(x) x$pfamily, character(1)), collapse = ", "),
  "\n"
)

pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
cat(
  "ING pfamily ok:",
  paste(vapply(pf_ing, function(x) x$pfamily, character(1)), collapse = ", "),
  "\n"
)
pl0 <- pf_ing[["(Intercept)"]]$prior_list
cat(
  "ING intercept rate/(shape-1)=",
  pl0$rate / (pl0$shape - 1),
  " vs prior_list=",
  ps$pop.prior_list[["(Intercept)"]]$dispersion, "\n"
)
