## Validate group_precision_floor() internals against _chk_group_marginal_bound.R
## scratch fixture (logit, J = 10, p = 3, weak data).

devtools::load_all("c:/Rpackages/lmebayesCore", quiet = TRUE)

set.seed(2026)
J <- 10L; p <- 3L; nobs <- 10L
idx <- function(j) ((j - 1L) * p + 1L):(j * p)

Zl <- lapply(1:J, function(j)
  cbind(1, matrix(rnorm(nobs * (p - 1L)), nobs, p - 1L)))

Rb     <- matrix(c(1, .30, -.20, .30, 1, .25, -.20, .25, 1), p, p)
sb     <- c(.75, .55, .45)
Sig_b  <- diag(sb) %*% Rb %*% diag(sb)
Lam_b1 <- solve(Sig_b)
mu0    <- c(-0.4, 0.3, -0.2)
gam_tr <- mu0 + 0.25 * rnorm(p)
b_tr   <- as.vector(vapply(1:J, function(j)
  gam_tr + drop(t(chol(Sig_b)) %*% rnorm(p)), numeric(p)))

n_tr0 <- 3L; lg0 <- 0.05; eps <- 0.01
set.seed(7)
yl <- lapply(1:J, function(j)
  rbinom(nobs, n_tr0, plogis(drop(Zl[[j]] %*% b_tr[idx(j)]))))

Lam_g <- diag(rep(lg0, p))
Mcpl  <- Lam_b1 %*% solve(Lam_g + J * Lam_b1) %*% Lam_b1
Lam_b <- kronecker(diag(J), Lam_b1) - kronecker(matrix(1, J, J), Mcpl)
mu_b  <- rep(mu0, J)

family_hook <- .group_floor_family(stats::binomial())

group_data <- list(
  group_levels = as.character(seq_len(J)),
  Z_list = Zl,
  y_list = yl,
  eta_offset = rep(list(rep(0, nobs)), J),
  size_list = rep(list(rep(n_tr0, nobs)), J),
  wt_list = rep(list(rep(1, nobs)), J),
  n_obs = rep(nobs, J),
  dispersion = 1,
  J = J,
  p_re = p
)

prior_int <- list(
  Lambda_beta = Lam_b,
  mu_beta = mu_b,
  P_b = Lam_b1,
  p_re = p,
  J = J,
  n_dim = J * p
)

engine <- .group_floor_engine(
  prior_int = prior_int,
  group_data = group_data,
  family_hook = family_hook
)

Gbar <- .group_floor_build_gbar(group_data, family_hook)
kappa <- vapply(seq_len(J), function(j) {
  .group_floor_kappa_crude(j, Lam_b, Gbar, p, J)
}, numeric(1L))

r_grp <- .gamma_level_for_budget(p, eps / J)
R_crude <- r_grp + 2 * max(kappa)

kappa_lap <- vapply(seq_len(J), function(j) {
  .group_floor_kappa_laplace(j, engine, engine$emb(j, Zl[[j]][1, ]), r_grp)
}, numeric(1L))
R_lap <- r_grp + 2 * max(kappa_lap)

fi_crude <- lapply(seq_len(J), function(j) {
  .group_floor_one_group(j, engine, group_data, family_hook, R_crude)
})
fi_lap <- lapply(seq_len(J), function(j) {
  .group_floor_one_group(j, engine, group_data, family_hook, R_lap)
})

ePb <- eigen(Lam_b1, symmetric = TRUE)
Pmh <- ePb$vectors %*% diag(1 / sqrt(ePb$values)) %*% t(ePb$vectors)
omega_of <- function(Gm) 1 / (1 + min(eigen(Pmh %*% Gm %*% Pmh,
                                            symmetric = TRUE)$values))

cat("=== crude kappa level (matches scratch grp_crude) ===\n")
cat(sprintf("R = %.4f  eta_star_max = %.4f  wbar_min = %.3e  omega_max = %.6f\n",
            R_crude,
            max(vapply(fi_crude, function(z) max(z$es), 0)),
            min(vapply(fi_crude, function(z) min(z$wbar), 0)),
            max(vapply(fi_crude, function(z) omega_of(z$Gamma), 0))))

cat("\n=== laplace kappa level (matches scratch grp_laplace) ===\n")
cat(sprintf("R = %.4f  eta_star_max = %.4f  wbar_min = %.3e  omega_max = %.6f\n",
            R_lap,
            max(vapply(fi_lap, function(z) max(z$es), 0)),
            min(vapply(fi_lap, function(z) min(z$wbar), 0)),
            max(vapply(fi_lap, function(z) omega_of(z$Gamma), 0))))
cat("\nScratch reference grp_laplace: eta~11.63, wbar~2.67e-5, omega~0.9978.\n")
