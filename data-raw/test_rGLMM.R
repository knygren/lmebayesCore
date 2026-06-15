## Integration test for rGLMM.
## Run: Rscript data-raw/test_rGLMM.R

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
start <- list(
  `(Intercept)` = stats::setNames(c(0, 0), colnames(X_int)),
  slope         = stats::setNames(0, colnames(X_slp))
)

n_main <- 12L
m_main <- 4L
y_gauss <- eta + rnorm(length(eta), sd = 0.5)
fam <- gaussian()

common <- list(
  y = y_gauss, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list = pl1_gauss,
  pfamily_list = pfam_list,
  start = start,
  group_levels = group_levels,
  family = fam,
  m_convergence = m_main,
  progbar = FALSE
)

## 1. No pilot: rGLMM matches v2 on same seed.
set.seed(515)
fit0 <- do.call(glmbayesCore::rGLMM, c(list(n = n_main, n_pilot = 0L), common))
set.seed(515)
v2_common <- common
v2_common$fixef_start <- v2_common$start
v2_common$start <- NULL
v2_common$prior_list_block1 <- v2_common$prior_list
v2_common$prior_list <- NULL
v2_only <- do.call(glmbayesCore::two_block_rNormal_reg_v2,
                   c(list(n = n_main), v2_common))

for (k in re_names) {
  d <- max(abs(fit0$fixef[[k]] - v2_only$fixef_draws[[k]]))
  if (!is.finite(d) || d > 1e-12) {
    stop("rGLMM (no pilot) vs v2 mismatch for ", k, ": max = ", d)
  }
}
stopifnot(is.null(fit0$pilot))
stopifnot(inherits(fit0, "rGLMM"))

## 2. Pilot without UB: C++ staged path (tv_tol NULL).
set.seed(919)
fit_p <- do.call(
  glmbayesCore::rGLMM,
  c(list(n = n_main, n_pilot = 6L, m_convergence_pilot = 3L), common)
)
stopifnot(!is.null(fit_p$pilot))
stopifnot(!is.null(fit_p$pilot_chisq))
stopifnot(is.finite(fit_p$pilot_chisq$p_value))
stopifnot(is.list(fit_p$fixef.init))
stopifnot(!is.null(fit_p$pilot$coefficients))

## 3. Pilot with UB: tv_tol triggers R calibration path.
set.seed(707)
fit_ub <- do.call(
  glmbayesCore::rGLMM,
  c(
    list(
      n = n_main,
      n_pilot = 5L,
      m_convergence_pilot = 3L,
      tv_tol = 1e-2
    ),
    common
  )
)
stopifnot(!is.null(fit_ub$pilot_ub))
stopifnot(is.finite(fit_ub$pilot_ub$m_min_upper))
stopifnot(fit_ub$m_convergence >= m_main)
stopifnot(!is.null(fit_ub$pilot_chisq))

cat("test_rGLMM.R: OK\n")
