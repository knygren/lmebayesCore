## pfamily_list() — population / group-effect pfamilies (Ex_13b data)

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design0 <- model_setup(form_lmer, data = dat)
dat <- subset(dat, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
dat$school_id <- droplevels(dat$school_id)

ps <- Prior_Setup_GLMM(form_lmer, data = dat, pop.pwt = 0.01)

## Default: dNormal() per component (known tau^2_k)
pf_n <- pfamily_list(ps)
names(pf_n)   # population models: (Intercept), slopes, ...

## Print a subset of population models
print(pf_n, components = c("(Intercept)", "distracted_a1"))

## Same via subsetting (class preserved)
print(pf_n[c("(Intercept)", "distracted_ppvt")])

## Variation: Independent Normal-Gamma for every component
pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
print(pf_ing, components = "distracted_ppvt")

## Variation: mix types by name
pf_mix <- pfamily_list(
  ps,
  ptypes = c(
    "(Intercept)" = "dNormal",
    distracted_ppvt = "dIndependent_Normal_Gamma",
    distracted_a1 = "dNormal"
  )
)
print(pf_mix, components = c("distracted_ppvt", "distracted_a1"))
