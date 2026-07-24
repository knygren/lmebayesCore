devtools::load_all(".", quiet = TRUE)

set.seed(20260720)
J   <- 5L
n_j <- 8L
n_obs <- J * n_j
group_levels <- paste0("g", seq_len(J))
group <- factor(rep(group_levels, each = n_j), levels = group_levels)
attr(group, "group_name") <- "group"

tau2_true   <- 0.7
sigma2_true <- stats::setNames(c(0.3, 0.6, 0.9, 1.2, 1.5), group_levels)
b_true <- stats::setNames(stats::rnorm(J, sd = sqrt(tau2_true)), group_levels)
y <- 2 + b_true[as.character(group)] +
  stats::rnorm(n_obs, sd = sqrt(sigma2_true[as.character(group)]))

D <- matrix(1, n_obs, 1, dimnames = list(NULL, "(Intercept)"))
W <- list("(Intercept)" = matrix(1, J, 1, dimnames = list(NULL, "(Intercept)")))

design <- list(
  y = y, Z = D, groups = group, X_hyper = W,
  re_coef_names = "(Intercept)", group_name = "group"
)
mpl <- list(
  Sigma_ranef      = matrix(tau2_true, 1, 1),
  dispersion_ranef = unname(sigma2_true),
  prior_list       = list(
    "(Intercept)" = list(
      mu_fixef = 0, Sigma_fixef = matrix(100),
      dispersion_fixef = tau2_true
    )
  )
)

pm  <- lmerb_posterior_mean(design, mpl)
cat("posterior mean fixef:\n")
print(pm$fixef)

Sigma <- lmerb_posterior_covariance(design, mpl)
cat("\nposterior covariance:\n")
print(Sigma)
cat("\ndim:", paste(dim(Sigma), collapse = " x "), "\n")
cat("dimnames:", paste(rownames(Sigma), collapse = ", "), "\n")
cat("symmetric:", isSymmetric(Sigma), "\n")
cat("positive-definite (all eigenvalues > 0):", all(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values > 0), "\n")

cat("\nALL OK\n")
