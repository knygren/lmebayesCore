devtools::load_all(quiet = TRUE)
set.seed(42L)
J <- 4L
n_per <- 5L
g <- factor(rep(seq_len(J), each = n_per))
n_obs <- length(g)
re_nm <- "(Intercept)"
y <- stats::rnorm(n_obs)
x <- matrix(1, nrow = n_obs, ncol = 1L, dimnames = list(NULL, re_nm))
gl <- as.character(seq_len(J))
x_hyper <- stats::setNames(
  list(matrix(1, J, 1L, dimnames = list(gl, re_nm))),
  re_nm
)
P <- matrix(1)
ing <- list(
  mu = matrix(0, 1, 1, dimnames = list(re_nm, NULL)),
  Sigma = matrix(1),
  shape = 2,
  rate = 1,
  max_disp_perc = 0.99,
  disp_lower = stats::qgamma(0.01, shape = 2, rate = 1),
  disp_upper = stats::qgamma(0.99, shape = 2, rate = 1)
)
pf <- dNormal(mu = c(`(Intercept)` = 0), Sigma = matrix(1), dispersion = 1)
pfamily_list <- stats::setNames(list(pf), re_nm)
disp_fix <- 1

message("known vcov...")
k <- rLMMindepNormalGamma_reg_known_vcov(
  n = 2L,
  y = y,
  x = x,
  block = g,
  x_hyper = x_hyper,
  P = P,
  prior_list = ing,
  pfamily_list = pfamily_list,
  dispersion_fix = disp_fix,
  m_convergence = 2L,
  progbar = FALSE,
  verbose = FALSE
)
stopifnot(
  length(k$dispersion_ranef) == 2L,
  length(k$dispersion_ranef.mean) == 1L,
  is.null(k$dispersion_ranef_draws)
)
message("known OK: dispersion_ranef.mean = ", k$dispersion_ranef.mean)

pf_ing <- dIndependent_Normal_Gamma(
  mu = c(`(Intercept)` = 0),
  Sigma = matrix(1),
  shape = 2,
  rate = 1,
  disp_lower = 0.01,
  disp_upper = 100
)
pfamily_ing <- stats::setNames(list(pf_ing), re_nm)
fixef_start <- stats::setNames(list(c(`(Intercept)` = 0)), re_nm)
message("estimated vcov...")
e <- rLMMindepNormalGamma_reg_estimated_vcov(
  n = 1L,
  y = y,
  x = x,
  block = g,
  x_hyper = x_hyper,
  P = P,
  prior_list = ing,
  pfamily_list = pfamily_ing,
  dispersion_fix = disp_fix,
  start = fixef_start,
  m_convergence = 1L,
  tv_tol = 0.05,
  progbar = FALSE,
  verbose = FALSE,
  gap_tol = NULL
)
stopifnot(
  length(e$dispersion_ranef) == 1L,
  length(e$dispersion_ranef.mean) == 1L,
  is.null(e$dispersion_ranef_draws)
)
message("estimated OK")
