############################### Start of EnvelopeOrchestrator example ####################

# This example demonstrates calling EnvelopeOrchestrator directly for Gaussian
# regression with an independent Normal-Gamma prior. It mirrors the algorithm
# path used inside rIndepNormalGammaReg:
#   - Step A: Initial dispersion via weighted lm.wfit residual variance
#   - Step B: EnvelopeCentering loop (sample beta, compute RSS_Post2, update
#             dispersion2 via Gamma posterior) to anchor the envelope
#   - Step C: Coefficient posterior mode optimization (optim + f2/f3)
#   - Step D: Standardize the model (glmb_Standardize_Model)
#   - Step E: EnvelopeOrchestrator (EnvelopeBuild + EnvelopeDispersionBuild
#             + EnvelopeSort in one call)
# It stops after envelope construction (no sampling).

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

ps <- Prior_Setup(weight ~ group, gaussian())

x <- as.matrix(ps$x)
y <- as.vector(ps$y)
mu <- ps$mu
Sigma <- ps$Sigma
shape <- ps$shape
rate <- ps$rate

n_obs <- length(y)
wt <- rep(1, n_obs)
offset2 <- rep(0, n_obs)

# Reconstruct coefficient precision P (matches rindepNormalGamma_reg)
Rchol <- chol(Sigma)
Pinv <- chol2inv(Rchol)
P <- 0.5 * (Pinv + t(Pinv))

famfunc <- glmbfamfunc(gaussian())
f2 <- famfunc$f2
f3 <- famfunc$f3

Gridtype_core <- as.integer(2)

###############################################################################
# Step A: Initial dispersion2 via weighted lm.wfit residual variance
###############################################################################
y_star <- y - offset2
fit0 <- lm.wfit(x = x, y = y_star, w = wt)
res0 <- fit0$residuals
RSS0 <- sum(res0^2)
p_rank <- as.integer(fit0$rank)

dispersion2 <- RSS0 / (n_obs - p_rank)

###############################################################################
# Step B: EnvelopeCentering loop (RSS_Post2 and dispersion2 anchoring)
###############################################################################
n_beta_draws <- as.integer(10000)
n_rss_iter <- as.integer(10)

n_w <- sum(wt)

RSS_Post2 <- NA_real_
cpp_out <- NULL

for (j in seq_len(n_rss_iter)) {
  prior_list_loop <- list(
    mu = mu,
    P = P,
    dispersion = dispersion2
  )

  cpp_out <- rNormal_reg(
    n = n_beta_draws,
    y = y,
    x = x,
    prior_list = prior_list_loop,
    offset = offset2,
    weights = wt,
    family = gaussian(),
    Gridtype = Gridtype_core,
    use_parallel = FALSE,
    use_opencl = FALSE,
    verbose = FALSE,
    progbar = FALSE
  )

  beta_draws <- cpp_out$coefficients

  lp_mat <- beta_draws %*% t(x)
  eta_mat <- lp_mat + matrix(offset2, nrow = n_beta_draws, ncol = n_obs, byrow = TRUE)
  diff_mat <- eta_mat - matrix(y, nrow = n_beta_draws, ncol = n_obs, byrow = TRUE)
  res_sq <- diff_mat * diff_mat
  res_sq_weighted <- res_sq * matrix(wt, nrow = n_beta_draws, ncol = n_obs, byrow = TRUE)
  RSS_temp <- rowSums(res_sq_weighted)
  RSS_Post2 <- mean(RSS_temp)

  shape2 <- shape + n_w / 2.0
  rate2 <- rate + RSS_Post2 / 2.0
  dispersion2 <- rate2 / (shape2 - 1.0)
}

###############################################################################
# Step C: Coefficient posterior mode optimization (optim + f2/f3)
###############################################################################
dispstar <- dispersion2
wt2_opt <- wt / dispstar
alpha <- as.vector(x %*% as.vector(mu) + offset2)

mu2_opt <- rep(0, length(as.vector(mu)))
parin <- rep(0, length(as.vector(mu)))

opt_out <- optim(
  par = parin,
  fn = f2,
  gr = f3,
  y = as.vector(y),
  x = as.matrix(x),
  mu = as.vector(mu2_opt),
  P = as.matrix(P),
  alpha = as.vector(alpha),
  wt = as.vector(wt2_opt),
  method = "BFGS",
  hessian = TRUE
)

bstar <- opt_out$par
A1 <- opt_out$hessian

###############################################################################
# Step D: Standardize model (glmb_Standardize_Model)
###############################################################################
Standard_Mod <- glmb_Standardize_Model(
  y = as.vector(y),
  x = as.matrix(x),
  P = as.matrix(P),
  bstar = as.matrix(bstar, ncol = 1),
  A1 = as.matrix(A1)
)

bstar2 <- Standard_Mod$bstar2
A <- Standard_Mod$A
x2_std <- Standard_Mod$x2
mu2_std <- Standard_Mod$mu2
P2_std <- Standard_Mod$P2

###############################################################################
# Step E: EnvelopeOrchestrator (EnvelopeBuild + EnvelopeDispersionBuild
#         + EnvelopeSort in one call)
###############################################################################
max_disp_perc <- 0.99
n_env <- as.integer(200)
Gridtype_env <- as.integer(3)  # EnvelopeOrchestrator overrides to 3 for unknown dispersion

env_out <- EnvelopeOrchestrator(
  bstar2 = as.vector(bstar2),
  A = as.matrix(A),
  y = as.vector(y),
  x2 = as.matrix(x2_std),
  mu2 = as.matrix(mu2_std, ncol = 1),
  P2 = as.matrix(P2_std),
  alpha = as.vector(alpha),
  wt = as.vector(wt),
  n = n_env,
  Gridtype = Gridtype_env,
  n_envopt = as.integer(1),
  shape = shape,
  rate = rate,
  RSS_Post2 = RSS_Post2,
  RSS_ML = NA_real_,
  max_disp_perc = max_disp_perc,
  disp_lower = NULL,
  disp_upper = NULL,
  use_parallel = TRUE,
  use_opencl = FALSE,
  verbose = FALSE
)

# Output structure matches that from the step-by-step Ex_EnvelopeDispersionBuild
print(env_out$low)
print(env_out$upp)
print(env_out$gamma_list[c("shape3", "rate2")])

env_out

###############################################################################
# End: envelope construction only
###############################################################################
