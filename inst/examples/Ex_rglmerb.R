## rglmerb() — matrix-level GLMM sampler (binomial) on book_banning
##
## Probability that a challenged book is removed, with a random intercept and
## a `violent` slope by state.  Restricted to the eight largest states so the
## example stays quick under R CMD check.  Requires suggested package
## bayesrules.  For the formula interface see lmebayes::glmerb().

if (requireNamespace("bayesrules", quietly = TRUE)) {
  data(book_banning, package = "bayesrules")
  dat <- book_banning[, c("state", "removed", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
  dat$violent_i <- as.integer(
    dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
  )
  dat$state <- factor(dat$state)

  keep <- names(sort(table(dat$state), decreasing = TRUE))[seq_len(8L)]
  dat <- droplevels(subset(dat, state %in% keep))

  form <- removed_i ~ violent_i + (1 + violent_i || state)

  design <- model_setup(form, data = dat, family = binomial())
  ps <- Prior_Setup_GLMM(form, data = dat, family = binomial(), pop.pwt = 0.01)

  ## Point estimates only: joint posterior mode by ICM, no sampling.
  pt <- rglmerb(
    design       = design,
    pfamily_list = pfamily_list(ps),
    family       = binomial(),
    simulate     = FALSE,
    verbose      = FALSE
  )
  pt$popef.mode

  ## Draws.  Non-Gaussian families always run a pilot stage before the main
  ## chains, and dispprior_list must stay NULL for binomial().
  set.seed(1)
  fit <- rglmerb(
    n            = 10L,
    design       = design,
    pfamily_list = pfamily_list(ps),
    family       = binomial(),
    verbose      = FALSE,
    progbar      = FALSE
  )

  print_groupef(fit, draws = 1:3, groups = 1:3)
}
