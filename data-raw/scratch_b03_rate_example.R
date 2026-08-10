## Scratch: verify the small runnable two_block_rate() example used in
## vignettes/Chapter-B02.Rmd and Chapter-B03.Rmd.
library(lmebayesCore)

set.seed(1)
J     <- 12L                      # groups
n_j   <- 8L                       # observations per group
sigma2 <- 1.0                     # observation variance
tau2   <- 0.5                     # random-intercept variance
V_gamma <- 100                    # Block 2 prior variance on gamma

group <- factor(rep(seq_len(J), each = n_j))
Z     <- matrix(1, nrow = J * n_j, ncol = 1L,
                dimnames = list(NULL, "(Intercept)"))

x_hyper <- list("(Intercept)" = matrix(1, nrow = J, ncol = 1L,
                                       dimnames = list(levels(group),
                                                       "(Intercept)")))

rate <- two_block_rate(
  x                 = Z,
  group             = group,
  x_hyper           = x_hyper,
  prior_list_block1 = list(Sigma = matrix(tau2), dispersion = sigma2),
  prior_list_block2 = list("(Intercept)" = list(mu = 0,
                                                Sigma = matrix(V_gamma)))
)
print(rate)

cat("\nlambda_star  =", rate$lambda_star, "\n")
cat("closed form  =", (1 / tau2) / (1 / tau2 + n_j / sigma2), "\n")

cat("\nTV bound at l = 1..8 (theorem3):\n")
print(signif(lmebayesCore:::two_block_tv_bound(rate, 1:8), 4))

cat("\nsweeps for tv_tol = 0.01:",
    lmebayesCore:::two_block_l_for_tv(rate, 0.01), "\n")

cat("\npilot cost optimum:\n")
opt <- lmebayesCore:::two_block_optimize_pilot_cost(
  n = 1000L, rate = rate, tv_tol = 0.01,
  m_convergence_pilot = 5L, p = 1L
)
str(opt[c("n_pilot_opt", "m_convergence_opt", "total_cost_opt")])
