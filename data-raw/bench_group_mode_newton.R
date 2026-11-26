# Newton vs vmmin group mode finder (glmbayesCore group_mode_newton.cpp).
#
# Run by hand after installing glmbayesCore, then lmebayesCore, then
# lmebayes:
#   Rscript data-raw/bench_group_mode_newton.R
# No set.seed(): agreement in part 3 is judged against Monte Carlo error.
#
# Model: the book_banning logit glmerb of lmebayes tests/manual
# (_glmerb_book_banning_ridge.R), with few draws and rwt = 0.75 so each fit
# needs few sweeps.
#   1. Accuracy: options(lmebayesCore.mode_method) "newton_compare" and
#      "newton_f4_compare" keep vmmin downstream and report the
#      "[... Newton compare]" message: expect flagged 0, max d well under
#      0.01, max rel Hessian diff around 1e-8 or less, 1-3 iterations per
#      call; the f2_f3_f4 kernel line should be near machine precision.
#   2. Speed: "vmmin", "newton" and "newton_f4" with
#      options(lmebayesCore.block_timing = TRUE); compare the "mode parts"
#      lines and the elapsed times.
#   3. Downstream: population-effect means and SDs of the Newton runs
#      against the vmmin run.

library(lmebayes)
data("book_banning", package = "bayesrules")

n_draws <- 200L
rwt <- 0.75

dat <- book_banning
st_vars <- c("political_value_index", "median_income", "hs_grad_rate")
dat <- dat[stats::complete.cases(dat[, c(
  "removed", "explicit", "language", "violent", "state", st_vars
)]), ]
dat$removed_i  <- as.integer(as.character(dat$removed) == "1")
dat$explicit_i <- as.integer(as.character(dat$explicit) == "1")
dat$language_i <- as.integer(as.character(dat$language) == "1")
dat$violent_i  <- as.integer(as.character(dat$violent) == "1")
st <- unique(dat[, c("state", st_vars)])
z <- function(v) {
  s <- as.numeric(scale(st[[v]]))
  s[match(dat$state, st$state)]
}
dat$pvi_c    <- z("political_value_index")
dat$income_c <- z("median_income")
dat$hs_c     <- z("hs_grad_rate")
dat$state <- droplevels(factor(dat$state))

form <- removed_i ~ explicit_i + language_i + violent_i +
  pvi_c + income_c + hs_c +
  (1 + explicit_i + language_i + violent_i || state)
ps <- Prior_Setup_GLMM(form, data = dat, family = binomial(), pop.pwt = 0.01)

run_fit <- function(mode_method, block_timing = FALSE) {
  old <- options(lmebayesCore.mode_method = mode_method,
                 lmebayesCore.block_timing = block_timing)
  on.exit(options(old))
  cat("\n=== mode_method = \"", mode_method, "\" ===\n", sep = "")
  t0 <- proc.time()[["elapsed"]]
  fit <- glmerb(form, data = dat, family = binomial(),
                pfamily_list = pfamily_list(ps), n = n_draws, rwt = rwt)
  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("elapsed %.1f s\n", elapsed))
  list(fit = fit, elapsed = elapsed)
}

pop_summary <- function(fit) {
  pe <- fit$popef
  mats <- if (is.list(pe) && !is.data.frame(pe)) pe else list(pe)
  m <- unlist(lapply(mats, function(x) colMeans(as.matrix(x))))
  s <- unlist(lapply(mats, function(x) apply(as.matrix(x), 2, stats::sd)))
  list(mean = m, sd = s)
}

## 1. Accuracy ---------------------------------------------------------------

invisible(run_fit("newton_compare"))
invisible(run_fit("newton_f4_compare"))

## 2. Speed ------------------------------------------------------------------

runs <- list(
  vmmin     = run_fit("vmmin", block_timing = TRUE),
  newton    = run_fit("newton", block_timing = TRUE),
  newton_f4 = run_fit("newton_f4", block_timing = TRUE)
)
cat("\nelapsed (s):\n")
print(vapply(runs, function(r) r$elapsed, numeric(1)))

## 3. Downstream agreement ---------------------------------------------------

ref <- pop_summary(runs$vmmin$fit)
for (nm in c("newton", "newton_f4")) {
  s <- pop_summary(runs[[nm]]$fit)
  zdiff <- (s$mean - ref$mean) / (ref$sd * sqrt(2 / n_draws))
  cat(sprintf("\n%s vs vmmin: population-effect means, difference in MC SEs\n",
              nm))
  print(round(cbind(vmmin = ref$mean, other = s$mean, z = zdiff), 3))
  cat(sprintf("max |z| %.2f (up to about 3 is unremarkable over %d terms)\n",
              max(abs(zdiff)), length(zdiff)))
}
