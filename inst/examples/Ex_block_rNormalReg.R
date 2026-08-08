## block_rNormalReg() — independent Gaussian regressions by school
##
## Same big_word_club data and first 3 school_id levels as printed in
## Ex_rLMM_reg, but each school is its own model via Prior_SetupBlock()
## (not Prior_Setup_GLMM / hierarchical Psi). Within-school formula is
## intercept-only so every small school stays full rank. Requires bayesrules.

if (requireNamespace("bayesrules", quietly = TRUE)) {
  data(big_word_club, package = "bayesrules")
  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  dat <- subset(
    dat,
    !is.na(score_ppvt) &
      !is.na(invalid_ppvt) & invalid_ppvt == 0L &
      complete.cases(dat[, c("score_ppvt", "school_id")])
  )
  keep <- levels(dat$school_id)[1:3]
  dat <- droplevels(subset(dat, school_id %in% keep))

  form_school <- score_ppvt ~ 1
  ps_block <- Prior_SetupBlock(
    form_school,
    block = "school_id",
    data = dat,
    family = gaussian(),
    pwt = 0.01
  )
  prior_lists <- lapply(ps_block, function(ps) {
    list(
      mu = as.numeric(ps$mu),
      Sigma = ps$Sigma,
      dispersion = ps$dispersion,
      ddef = FALSE
    )
  })

  mf <- model.frame(form_school, data = dat)
  y <- model.response(mf)
  x <- model.matrix(form_school, data = mf)
  block <- dat$school_id

  set.seed(1)
  out <- block_rNormalReg(
    n = 1L,
    y = y,
    x = x,
    block = block,
    prior_lists = prior_lists
  )

  out$coefficients
  out$coef.mode
}
