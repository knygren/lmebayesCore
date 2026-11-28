# Block-2 (population) draws: glmbayesCore rNormalReg per draw vs a factor
# of the posterior precision built once per sweep call
# (src/two_block_population_chol.cpp).
#
# Run by hand after installing lmebayesCore, then lmebayes, one part at a
# time from an R session (part 2 runs two fits).
#
# Model: the book_banning logit glmerb of lmebayes tests/manual
# (_glmerb_book_banning_ridge.R), group modes by "newton_f4".
#   1. Accuracy: options(lmebayesCore.pop_draw = "compare") keeps the
#      rNormalReg draws and makes the fixed-factor draw from the same R RNG
#      state; the "[... pop compare]" message should show flagged 0 and mean
#      and draw differences near machine precision. Matching draws (not just
#      means) also confirm the same posterior covariance.
#   2. Speed: "reg" and "chol" with options(lmebayesCore.block_timing = TRUE);
#      compare the block2 times and the elapsed times.
#   3. Downstream: population-effect means and SDs of the "chol" run against
#      the "reg" run (Monte Carlo agreement; block 1 uses the thread-local
#      RNG, so the runs are not draw-for-draw equal).

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

run_fit <- function(pop_draw, block_timing = FALSE) {
  old <- options(lmebayesCore.pop_draw = pop_draw,
                 lmebayesCore.mode_method = "newton_f4",
                 lmebayesCore.block_timing = block_timing)
  on.exit(options(old))
  cat("\n=== pop_draw = \"", pop_draw, "\" ===\n", sep = "")
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

invisible(run_fit("compare"))

## 2. Speed ------------------------------------------------------------------

runs <- list(
  reg  = run_fit("reg", block_timing = TRUE),
  chol = run_fit("chol", block_timing = TRUE)
)
cat("\nelapsed (s):\n")
print(vapply(runs, function(r) r$elapsed, numeric(1)))

## 3. Downstream agreement ---------------------------------------------------

ref <- pop_summary(runs$reg$fit)
s <- pop_summary(runs$chol$fit)
zm <- (s$mean - ref$mean) / (ref$sd * sqrt(2 / n_draws))
zs <- (s$sd / ref$sd - 1) / sqrt(1 / n_draws)
cat("\nchol vs reg: population effects, differences in MC SEs\n")
print(round(cbind(reg_mean = ref$mean, chol_mean = s$mean, z_mean = zm,
                  sd_ratio = s$sd / ref$sd, z_sd = zs), 3))
cat(sprintf("max |z| mean %.2f, sd %.2f (up to about 3 is unremarkable)\n",
            max(abs(zm)), max(abs(zs))))
