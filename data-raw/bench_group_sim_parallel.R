# Serial vs parallel group simulation (glmbayesCore
# group_simulate_parallel.cpp).
#
# Run by hand after installing glmbayesCore, then lmebayesCore, then
# lmebayes, one part at a time from an R session (part 2 runs three fits).
# The parallel draws use the thread-local RNG, so no seed fixes them;
# agreement in part 3 is judged against Monte Carlo error.
#
# Model: the book_banning logit glmerb of lmebayes tests/manual
# (_glmerb_book_banning_ridge.R), group modes by "newton_f4".
#   1. Accuracy: options(glmbayesCore.group_sim = "compare") keeps the serial
#      draws and reports the "[... sim compare]" message: expect flagged 0,
#      f2 and back-transform differences near machine precision, and similar
#      candidates per draw for both paths.
#   2. Speed: "serial", "parallel", and "parallel" on one thread, with
#      options(lmebayesCore.block_timing = TRUE); compare the "simulate"
#      share and the elapsed times. The one-thread run separates removing
#      the per-draw R overhead from the threading gain.
#   3. Downstream: population-effect means and SDs of the parallel runs
#      against the serial run.

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

run_fit <- function(group_sim, block_timing = FALSE, threads = NULL) {
  old <- options(glmbayesCore.group_sim = group_sim,
                 lmebayesCore.mode_method = "newton_f4",
                 lmebayesCore.block_timing = block_timing)
  on.exit(options(old))
  if (!is.null(threads)) {
    RcppParallel::setThreadOptions(numThreads = threads)
    on.exit(RcppParallel::setThreadOptions(numThreads = "auto"), add = TRUE)
  }
  label <- if (is.null(threads)) group_sim else
    sprintf("%s, %d thread(s)", group_sim, threads)
  cat("\n=== group_sim = \"", label, "\" ===\n", sep = "")
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
  serial     = run_fit("serial", block_timing = TRUE),
  parallel   = run_fit("parallel", block_timing = TRUE),
  parallel_1 = run_fit("parallel", block_timing = TRUE, threads = 1L)
)
cat("\nelapsed (s):\n")
print(vapply(runs, function(r) r$elapsed, numeric(1)))

## 3. Downstream agreement ---------------------------------------------------

ref <- pop_summary(runs$serial$fit)
for (nm in c("parallel", "parallel_1")) {
  s <- pop_summary(runs[[nm]]$fit)
  zm <- (s$mean - ref$mean) / (ref$sd * sqrt(2 / n_draws))
  zs <- (s$sd / ref$sd - 1) / sqrt(1 / n_draws)
  cat(sprintf("\n%s vs serial: population effects, differences in MC SEs\n",
              nm))
  print(round(cbind(serial_mean = ref$mean, other_mean = s$mean,
                    z_mean = zm, sd_ratio = s$sd / ref$sd, z_sd = zs), 3))
  cat(sprintf("max |z| mean %.2f, sd %.2f (up to about 3 is unremarkable)\n",
              max(abs(zm)), max(abs(zs))))
}
