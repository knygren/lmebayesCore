devtools::load_all(quiet = TRUE)
set.seed(42L)
J <- 4L
n_per <- 5L
g <- factor(rep(seq_len(J), each = n_per))
re_nm <- "(Intercept)"
y <- stats::rnorm(length(g))
x <- matrix(1, nrow = length(g), ncol = 1L, dimnames = list(NULL, re_nm))
gl <- as.character(seq_len(J))
x_hyper <- stats::setNames(
  list(matrix(1, J, 1L, dimnames = list(gl, re_nm))),
  re_nm
)
ing <- list(
  mu = matrix(0, 1, 1, dimnames = list(re_nm, NULL)),
  Sigma = matrix(1),
  shape = 2,
  rate = 1,
  max_disp_perc = 0.99
)
pf_ing <- dIndependent_Normal_Gamma(
  mu = c(`(Intercept)` = 0),
  Sigma = matrix(1),
  shape = 2,
  rate = 1,
  disp_lower = 0.01,
  disp_upper = 100
)
pfamily_ing <- stats::setNames(list(pf_ing), re_nm)
fixef_mode <- stats::setNames(list(c(`(Intercept)` = 0)), re_nm)
b0 <- matrix(0, J, 1, dimnames = list(gl, re_nm))
design <- list(
  y = y, Z = x, groups = g, X_hyper = x_hyper,
  re_coef_names = re_nm, group_name = "g"
)
block1_prior <- list(P = matrix(1), dispersion = 1, ddef = FALSE)

out <- glmbayesCore:::.rGLMM_sweep_ing_block1(
  n_chains = 1L,
  start_fixef = fixef_mode,
  inner_sweeps = 1L,
  design = design,
  block1_prior = block1_prior,
  ing_prior_list = ing,
  pfamily_list = pfamily_ing,
  family = gaussian(),
  re_names = re_nm,
  group_levels = gl,
  b_start = b0,
  ptypes = c(`(Intercept)` = "dIndependent_Normal_Gamma"),
  use_cpp_block2 = FALSE,
  progbar = FALSE
)
stopifnot(
  length(out$dispersion_ranef) == 1L,
  is.finite(out$dispersion_ranef)
)
message("estimated sweep (no pilot) OK")
