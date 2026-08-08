## rGLMM_reg_known_vcov() — binomial, dNormal population priors on book_banning
##
## Same data / formula path as lmebayes demo
## demo("Ex_16_glmerb_book_banning", package = "lmebayes"), restricted to the
## eight largest states so the example stays quick under R CMD check.
## Requires suggested package bayesrules.

if (requireNamespace("bayesrules", quietly = TRUE)) {
  data(book_banning, package = "bayesrules")
  dat <- book_banning[, c("state", "removed", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
  dat$violent_i <- as.integer(
    dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
  )
  dat$state <- factor(dat$state)

  ## Keep the largest states (avoids singular glmer fits on tiny subsets).
  keep <- names(sort(table(dat$state), decreasing = TRUE))[seq_len(8L)]
  dat <- droplevels(subset(dat, state %in% keep))

  form <- removed_i ~ violent_i + (1 + violent_i || state)

  design <- model_setup(form, data = dat, family = binomial())
  ps <- Prior_Setup_GLMM(form, data = dat, family = binomial(), pop.pwt = 0.01)
  pf <- pfamily_list(ps)

  set.seed(1)
  fit <- rGLMM_reg_known_vcov(
    n = 10L,
    y = design$y,
    D = design$D,
    group = design$group,
    W = design$W,
    pfamily_list = pf,
    dispprior_list = list(),
    family = binomial(),
    progbar = FALSE,
    verbose = FALSE
  )

  print(fit, draws = 1:5)
  print_groupef(fit, draws = 1:5, groups = 1:3)
}
