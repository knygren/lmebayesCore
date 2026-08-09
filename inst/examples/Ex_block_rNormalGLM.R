## block_rNormalGLM() — independent binomial regressions by state
##
## Same book_banning data / formula path as Ex_rGLMM_reg. Prior_Setup_GLMM()
## flags rank and estimability; only estimable states are then fit as
## independent GLMs via Prior_SetupBlock() + block_rNormalGLM().
## Binomial responses use cbind(success, failure). Requires bayesrules.

if (requireNamespace("bayesrules", quietly = TRUE)) {
  data(book_banning, package = "bayesrules")
  dat <- book_banning[, c("state", "removed", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$success <- as.integer(dat$removed == 1L | dat$removed == "1")
  dat$failure <- 1L - dat$success
  dat$violent_i <- as.integer(
    dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
  )
  dat$state <- factor(dat$state)

  ## Keep the largest states (same subset as Ex_rGLMM_reg).
  keep <- names(sort(table(dat$state), decreasing = TRUE))[seq_len(8L)]
  dat <- droplevels(subset(dat, state %in% keep))

  ## Hierarchical setup for rank / estimability (same RE structure as Ex_rGLMM_reg).
  form <- success ~ violent_i + (1 + violent_i || state)
  ps <- Prior_Setup_GLMM(
    form, data = dat, family = binomial(), pop.pwt = 0.01
  )
  est <- names(ps$design$groupef.estimable)[ps$design$groupef.estimable]
  dat <- droplevels(subset(dat, state %in% est))

  form_state <- cbind(success, failure) ~ violent_i
  ps_block <- Prior_SetupBlock(
    form_state,
    block = "state",
    data = dat,
    family = binomial(),
    pwt = 0.01
  )
  ## Default dNormal() per block; block_rNormalGLM takes the prior_list payload.
  pf <- pfamily_list(ps_block)
  prior_lists <- lapply(pf, `[[`, "prior_list")

  ## Prior_SetupBlock takes cbind(success, failure); block_rNormalGLM /
  ## rglmb expect proportions + trial weights.
  mf <- model.frame(form_state, data = dat)
  Y <- model.response(mf)
  weights <- rowSums(Y)
  y <- Y[, 1L] / weights
  x <- model.matrix(form_state, data = mf)

  set.seed(1)
  out <- block_rNormalGLM(
    n = 1L,
    y = y,
    x = x,
    block = dat$state,
    prior_lists = prior_lists,
    weights = weights,
    family = binomial(),
    use_parallel = FALSE
  )

  out$coefficients
  out$coef.mode
}
