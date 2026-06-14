## Smoke test for two_block_rNormal_reg_staged_cpp_export (Phase 2a).
## Run: Rscript data-raw/test_two_block_staged_cpp.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

set.seed(11)
J <- 6L
n_per <- 15L
grp <- factor(rep(sprintf("g%02d", seq_len(J)), each = n_per))
group_levels <- levels(grp)
w_j <- round(rnorm(J), 2)

z1 <- rnorm(J * n_per)
x_re <- cbind(`(Intercept)` = 1, slope = z1)
re_names <- colnames(x_re)

X_int <- cbind(1, w_j)
rownames(X_int) <- group_levels
colnames(X_int) <- c("(Intercept)", "w")
X_slp <- matrix(1, J, 1L, dimnames = list(NULL, "(Intercept)"))
x_hyper <- list(`(Intercept)` = X_int, slope = X_slp)

gamma_int <- c(1.0, 0.5)
gamma_slp <- 0.8
b_int <- as.numeric(X_int %*% gamma_int) + rnorm(J, sd = 0.4)
b_slp <- as.numeric(X_slp %*% gamma_slp) + rnorm(J, sd = 0.4)
eta <- b_int[as.integer(grp)] + b_slp[as.integer(grp)] * z1

pl1_gauss <- list(Sigma = diag(0.25, 2L), dispersion = 0.25)
pfam_list <- list(
  `(Intercept)` = dNormal(mu = c(0, 0), Sigma = diag(4, 2L), dispersion = 0.16),
  slope         = dNormal(mu = 0, Sigma = diag(4, 1L), dispersion = 0.16)
)
fixef_start <- list(
  `(Intercept)` = stats::setNames(c(0, 0), colnames(X_int)),
  slope         = stats::setNames(0, colnames(X_slp))
)

n_main <- 20L
m_main <- 4L
fam <- gaussian()
famfunc <- glmbfamfunc(fam)
famfunc_g <- glmbfamfunc(gaussian())
disp_block1 <- pl1_gauss$dispersion
ddef_block1 <- NULL

## No pilot: staged main should match a direct v2 call (same seed).
y_gauss <- eta + rnorm(length(eta), sd = 0.5)
set.seed(515)
staged0 <- glmbayesCore:::.two_block_rNormal_reg_staged_cpp(
  n_main = n_main,
  m_convergence_main = m_main,
  n_pilot = 0L,
  m_convergence_pilot = 1L,
  y = y_gauss, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  dispersion_block1 = disp_block1,
  ddef_block1 = ddef_block1,
  pfamily_list = pfam_list,
  fixef_start = fixef_start,
  group_levels = group_levels,
  family = fam$family, link = fam$link,
  f2 = famfunc$f2, f3 = famfunc$f3,
  f2_gauss = famfunc_g$f2, f3_gauss = famfunc_g$f3,
  offset = rep(0, length(y_gauss)),
  wt = rep(1, length(y_gauss)),
  progbar_main = FALSE,
  progbar_pilot = FALSE
)

set.seed(515)
v2_only <- glmbayesCore:::.two_block_rNormal_reg_v2_cpp(
  n = n_main, m_convergence = m_main,
  y = y_gauss, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  dispersion_block1 = disp_block1,
  ddef_block1 = ddef_block1,
  pfamily_list = pfam_list,
  fixef_start = fixef_start,
  group_levels = group_levels,
  family = fam$family, link = fam$link,
  f2 = famfunc$f2, f3 = famfunc$f3,
  f2_gauss = famfunc_g$f2, f3_gauss = famfunc_g$f3,
  offset = rep(0, length(y_gauss)),
  wt = rep(1, length(y_gauss)),
  progbar = FALSE
)

for (i in seq_along(re_names)) {
  d <- max(abs(staged0$fixef_draws[[i]] - v2_only$fixef_draws[[i]]))
  if (!is.finite(d) || d > 1e-12) {
    stop("no-pilot staged vs v2 mismatch for ", re_names[i], ": max = ", d)
  }
}
stopifnot(is.null(staged0$pilot))
stopifnot(staged0$m_convergence_used == m_main)

## With pilot: chi-squared structure present (UB deferred to Phase 2b R wrapper).
set.seed(919)
staged_p <- glmbayesCore:::.two_block_rNormal_reg_staged_cpp(
  n_main = n_main,
  m_convergence_main = m_main,
  n_pilot = 8L,
  m_convergence_pilot = 3L,
  y = y_gauss, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  dispersion_block1 = disp_block1,
  ddef_block1 = ddef_block1,
  pfamily_list = pfam_list,
  fixef_start = fixef_start,
  group_levels = group_levels,
  family = fam$family, link = fam$link,
  f2 = famfunc$f2, f3 = famfunc$f3,
  f2_gauss = famfunc_g$f2, f3_gauss = famfunc_g$f3,
  offset = rep(0, length(y_gauss)),
  wt = rep(1, length(y_gauss)),
  progbar_main = FALSE,
  progbar_pilot = FALSE
)

stopifnot(!is.null(staged_p$pilot))
stopifnot(!is.null(staged_p$pilot_chisq))
stopifnot(is.finite(staged_p$pilot_chisq$Q))
stopifnot(staged_p$pilot_chisq$df == sum(vapply(fixef_start, length, integer(1L))))
stopifnot(staged_p$m_convergence_used == m_main)
stopifnot(is.list(staged_p$fixef_main_start))
stopifnot(!is.null(staged_p$pilot$b_draws))

cat("test_two_block_staged_cpp.R: OK\n")
