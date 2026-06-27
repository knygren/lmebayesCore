pkgload::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
pkgload::load_all("c:/Rpackages/lmebayes", quiet = TRUE)
library(bayesrules)

dat <- book_banning[, c("state", "removed", "violent")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)
form <- removed_i ~ violent_i + (1 + violent_i || state)
ps <- Prior_Setup_lmebayes(form, dat, binomial(), pwt = 0.01)

run_fit <- function(engine) {
  assign(".rglmerb_engine", engine, envir = asNamespace("lmebayes"))
  on.exit(
    assign(".rglmerb_engine", "cpp_engine", envir = asNamespace("lmebayes")),
    add = TRUE
  )
  set.seed(42L)
  fit <- glmerb(
    form, data = dat, family = binomial(), pfamily_list = pfamily_list(ps),
    n = 500L, mode_gap_max = 1.0, progbar = FALSE, verbose = FALSE
  )
  init <- unlist(fit$fixef.init)
  list(
    engine = fit$convergence$draw_engine,
    m_pilot = fit$convergence$m_convergence_pilot,
    m_main = fit$convergence$m_convergence,
    n_pilot = fit$convergence$n_pilot,
    init = init,
    mode = unlist(fit$coef.mode)
  )
}

r_out <- run_fit("R_engine")
cpp_out <- run_fit("cpp_engine")

cat("R_engine:", r_out$engine, "m_pilot", r_out$m_pilot, "n_pilot", r_out$n_pilot, "\n")
cat("  init:", paste(round(r_out$init, 4), collapse = ", "), "\n")
cat("cpp_engine:", cpp_out$engine, "m_pilot", cpp_out$m_pilot, "n_pilot", cpp_out$n_pilot, "\n")
cat("  init:", paste(round(cpp_out$init, 4), collapse = ", "), "\n")
cat("mode:", paste(round(r_out$mode, 4), collapse = ", "), "\n")
