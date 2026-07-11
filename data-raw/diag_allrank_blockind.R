## Isolate BlockInd crash on all-rank fixture.
devtools::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
setwd("c:/Rpackages/lmebayes")
source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

fx <- .prepare_small5_all_full_rank_manual()
form <- fx$form
dat  <- fx$dat

ps <- Prior_Setup_lmebayes(
  form,
  data            = dat,
  pwt             = 0.01,
  pwt_measurement = 0.49
)
m <- ps$ing_prior_measurement
disp_pf <- dGamma(
  shape          = m$shape,
  rate           = m$rate,
  beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
  Inv_Dispersion = TRUE,
  disp_lower     = m$disp_lower,
  disp_upper     = m$disp_upper
)

design <- ps$design
Z <- as.matrix(design$Z)
y <- design$y
groups <- design$groups
re_names <- design$re_coef_names
group_levels <- levels(groups)
p_re <- length(re_names)

cat("k =", length(group_levels), " n =", length(y), " p_re =", p_re, "\n")
cat("re_rank all TRUE:", all(design$re_rank), "\n")
cat("obs per school: min =", min(table(groups)), "max =", max(table(groups)), "\n\n")

prior_list <- list(
  shape          = m$shape,
  rate           = m$rate,
  max_disp_perc  = 0.99,
  disp_lower     = m$disp_lower,
  disp_upper     = m$disp_upper,
  mu             = matrix(0, p_re, 1, dimnames = list(re_names, NULL)),
  P              = diag(p_re),
  Inv_Dispersion = TRUE
)

cat("=== Direct BlockInd n=1 verbose ===\n")
tryCatch(
  glmbayesCore:::.rIndepNormalGammaRegBlockInd_cpp(
    n             = 1L,
    y             = y,
    x             = Z,
    block         = groups,
    prior_list    = prior_list,
    offset        = rep(0, length(y)),
    wt            = rep(1, length(y)),
    p_re          = p_re,
    verbose       = TRUE,
    progbar       = FALSE,
    use_parallel  = FALSE,
    group_levels  = group_levels,
    re_names      = re_names
  ),
  error = function(e) {
    cat("\nERROR:", conditionMessage(e), "\n")
    traceback()
  }
)
