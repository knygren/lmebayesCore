## Temporary smoke test for the pop.*/group.* rename of Prior_Setup_lmebayes().
## Not package code -- investigation/verification only.

devtools::load_all(".")

dat <- lme4::sleepstudy

## --- pooled mode ------------------------------------------------------------
ps_pooled <- Prior_Setup_lmebayes(
  Reaction ~ Days + (Days || Subject),
  data = dat,
  pop.pwt = 0.01
)
stopifnot(!is.null(ps_pooled$group.ing_prior))
stopifnot(!lmebayesCore:::.lmebayes_ing_prior_is_grouped(ps_pooled$group.ing_prior))
stopifnot(!is.null(ps_pooled$pop.prior_list))
stopifnot(!is.null(ps_pooled$pop.ing_prior))
stopifnot(is.numeric(ps_pooled$group.pwt) && length(ps_pooled$group.pwt) == 1L)
stopifnot(is.numeric(ps_pooled$group.n_prior) && length(ps_pooled$group.n_prior) == 1L)
print(ps_pooled)

pf_pooled <- pfamily_list(ps_pooled)
stopifnot(is.list(pf_pooled))
tryCatch(
  dGamma_list(ps_pooled),
  error = function(e) cat("Expected error on pooled dGamma_list():", conditionMessage(e), "\n")
)

## --- group mode --------------------------------------------------------------
ps_group <- Prior_Setup_lmebayes(
  Reaction ~ Days + (Days || Subject),
  data                = dat,
  pop.pwt             = 0.01,
  group.pwt           = 0.1,
  dispformula         = ~Subject,
  group.max_disp_perc = 0.9,
  group.alpha_target  = 0.02
)
stopifnot(!is.null(ps_group$group.ing_prior))
stopifnot(lmebayesCore:::.lmebayes_ing_prior_is_grouped(ps_group$group.ing_prior))
stopifnot(length(ps_group$group.ing_prior) == nlevels(dat$Subject))
stopifnot(!is.null(ps_group$group.pwt_calibration))
stopifnot(identical(ps_group$group.alpha_target, 0.02))
print(ps_group)

disp_pf <- dGamma_list(ps_group)
stopifnot(length(disp_pf) == nlevels(dat$Subject))
stopifnot(all(vapply(disp_pf, inherits, logical(1L), "pfamily")))

pf_group <- pfamily_list(ps_group, ptypes = "dIndependent_Normal_Gamma")
stopifnot(is.list(pf_group))

cat("\nALL SMOKE CHECKS PASSED\n")
