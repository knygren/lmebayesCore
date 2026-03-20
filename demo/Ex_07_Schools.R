## -----------------------------
## Data
## -----------------------------
school   <- c("A","B","C","D","E","F","G","H")
estimate <- c(28.39, 7.94, -2.75, 6.82, -0.64, 0.63, 18.01, 12.16)
sd_obs   <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)

J           <- length(school)
theta       <- estimate                     # initialize at MLEs
sigma_y_sq  <- sd_obs^2                     # known sampling variances

## -----------------------------
## Independent Normal-Gamma prior
## -----------------------------
mu_mu    <- mean(estimate)
sigma_mu <- var(estimate)

n_prior  <- 0.5
disp_ML  <- var(estimate)

shape <- n_prior / 2
rate  <- disp_ML * shape

## -----------------------------
## MCMC settings
## -----------------------------
n_burn <- 1000
n_sim  <- 1000

theta_out        <- matrix(0, nrow = n_sim, ncol = J)
mu_out           <- numeric(n_sim)
sigma_theta_out  <- numeric(n_sim)
iters_out1       <- numeric(n_burn)
iters_out2       <- numeric(n_sim)

set.seed(102)
x_one <- matrix(1, nrow = J, ncol = 1)

## -----------------------------
## Burn-in
## -----------------------------
burn_time <- system.time({
  for (k in seq_len(n_burn)) {
    
    out_pop <- rlmb(
      1,
      y       = theta,
      x       = x_one,
      pfamily = dIndependent_Normal_Gamma(mu_mu, sigma_mu, shape, rate)
    )
    
    mu_theta       <- out_pop$coefficients[1, 1]
    sigma_theta_sq <- out_pop$dispersion
    
    for (j in seq_len(J)) {
      theta[j] <- rlmb(
        1,
        y       = estimate[j],
        x       = matrix(1),
        pfamily = dNormal(mu_theta, Sigma = sigma_theta_sq, dispersion = sigma_y_sq[j])
      )$coefficients[1, 1]
    }
    
    iters_out1[k] <- out_pop$iters
  }
})

## -----------------------------
## Main simulation
## -----------------------------
sim_time <- system.time({
  for (k in seq_len(n_sim)) {
    
    out_pop <- rlmb(
      1,
      y       = theta,
      x       = x_one,
      pfamily = dIndependent_Normal_Gamma(mu_mu, sigma_mu, shape, rate)
    )
    
    mu_theta       <- out_pop$coefficients[1, 1]
    sigma_theta_sq <- out_pop$dispersion
    
    for (j in seq_len(J)) {
      theta[j] <- rlmb(
        1,
        y       = estimate[j],
        x       = matrix(1),
        pfamily = dNormal(mu_theta, Sigma = sigma_theta_sq, dispersion = sigma_y_sq[j])
      )$coefficients[1, 1]
    }
    
    theta_out[k, ]       <- theta
    mu_out[k]            <- mu_theta
    sigma_theta_out[k]   <- sqrt(sigma_theta_sq)
    iters_out2[k]        <- out_pop$iters
  }
})


cat("Burn-in time:\n")
print(burn_time)

cat("\nSimulation time:\n")
print(sim_time)

## -----------------------------
## Posterior summaries
## -----------------------------
colMeans(theta_out)
mean(mu_out)
mean(sigma_theta_out)
sqrt(diag(var(theta_out)))

## -----------------------------
## Acceptance diagnostics
## -----------------------------
mean_iters_out1 <- mean(iters_out1)
cat("\nMean iters_out1 (burn-in; 1 / acceptance rate):\n")
print(mean_iters_out1)
cat("\nEstimated burn-in acceptance rate:\n")
print(1 / mean_iters_out1)

mean_iters_out2 <- mean(iters_out2)
cat("\nMean iters_out2 (1 / acceptance rate):\n")
print(mean_iters_out2)

cat("\nEstimated acceptance rate:\n")
print(1 / mean_iters_out2)

cat("\nCODA diagnostics (Independent NG):\n")

mcmc_ind_ng <- coda::mcmc(cbind(theta_out,
                               mu = mu_out,
                               sigma_theta = sigma_theta_out))
print(summary(mcmc_ind_ng))

cat("\nEffective sample size:\n")
print(coda::effectiveSize(mcmc_ind_ng))

cat("\nGeweke z-scores (first 10 parameters shown):\n")
gwe <- coda::geweke.diag(mcmc_ind_ng)$z
print(utils::head(gwe, 10))

mcmc_iters <- coda::mcmc(iters_out2)
cat("\nLag-1 autocorrelation (mu, sigma_theta, iters_out2):\n")
print(c(mu = coda::autocorr(coda::mcmc(mu_out), lag = 1),
        sigma_theta = coda::autocorr(coda::mcmc(sigma_theta_out), lag = 1),
        iters_out2 = coda::autocorr(mcmc_iters, lag = 1)))


readline("\n--- End of Independent NG sampler block. Press <Enter> to continue... ---")

######################################## Normal_Gamma prior ################################

## Prior values (conjugate Normal-Gamma)

mu_mu    <- mean(estimate)
sigma_mu <- var(estimate)

n_prior  <- 0.5
disp_ML  <- var(estimate)

shape <- n_prior / 2
rate  <- disp_ML * shape

## MCMC settings (conjugate NG)

n_burn_ng <- 1000
n_sim_ng  <- 1000

theta_ng            <- estimate
theta_out_ng        <- matrix(0, nrow = n_sim_ng, ncol = J)
mu_out_ng           <- numeric(n_sim_ng)
sigma_theta_out_ng  <- numeric(n_sim_ng)
iters_out1_ng       <- numeric(n_burn_ng)
iters_out2_ng       <- numeric(n_sim_ng)

set.seed(101)
x_one <- matrix(1, nrow = J, ncol = 1)

## -----------------------------
## Burn-in (timed)
## -----------------------------
burn_time_ng <- system.time({
  for (k in seq_len(n_burn_ng)) {
    
    out_pop <- rlmb(
      1,
      y       = theta_ng,
      x       = x_one,
      pfamily = dNormal_Gamma(mu_mu, sigma_mu / disp_ML, shape = shape, rate = rate)
    )
    
    mu_theta       <- out_pop$coefficients[1, 1]
    sigma_theta_sq <- out_pop$dispersion
    
    for (j in seq_len(J)) {
      theta_ng[j] <- rlmb(
        1,
        y       = estimate[j],
        x       = matrix(1),
        pfamily = dNormal(mu_theta, Sigma = sigma_theta_sq, dispersion = sigma_y_sq[j])
      )$coefficients[1, 1]
    }
    
    iters_out1_ng[k] <- out_pop$iters
  }
})

## -----------------------------
## Main simulation (timed)
## -----------------------------
sim_time_ng <- system.time({
  for (k in seq_len(n_sim_ng)) {
    
    out_pop <- rlmb(
      1,
      y       = theta_ng,
      x       = x_one,
      pfamily = dNormal_Gamma(mu_mu, sigma_mu / disp_ML, shape = shape, rate = rate)
    )
    
    mu_theta       <- out_pop$coefficients[1, 1]
    sigma_theta_sq <- out_pop$dispersion
    
    for (j in seq_len(J)) {
      theta_ng[j] <- rlmb(
        1,
        y       = estimate[j],
        x       = matrix(1),
        pfamily = dNormal(mu_theta, Sigma = sigma_theta_sq, dispersion = sigma_y_sq[j])
      )$coefficients[1, 1]
    }
    
    theta_out_ng[k, ]       <- theta_ng
    mu_out_ng[k]            <- mu_theta
    sigma_theta_out_ng[k]   <- sqrt(sigma_theta_sq)
    iters_out2_ng[k]        <- out_pop$iters
  }
})

## Print timing
cat("Conjugate NG burn-in time:\n")
print(burn_time_ng)

cat("\nConjugate NG simulation time:\n")
print(sim_time_ng)
## Posterior summaries (conjugate NG)

colMeans(theta_out_ng)
mean(mu_out_ng)
mean(sigma_theta_out_ng)
sqrt(diag(var(theta_out_ng)))

## -----------------------------
## Acceptance diagnostics (conjugate NG)
## -----------------------------
mean_iters_out1_ng <- mean(iters_out1_ng)
cat("\nMean iters_out1_ng (burn-in; 1 / acceptance rate):\n")
print(mean_iters_out1_ng)
cat("\nEstimated burn-in acceptance rate:\n")
print(1 / mean_iters_out1_ng)

mean_iters_out2_ng <- mean(iters_out2_ng)
cat("\nMean iters_out2_ng (1 / acceptance rate):\n")
print(mean_iters_out2_ng)
cat("\nEstimated acceptance rate:\n")
print(1 / mean_iters_out2_ng)

cat("\nCODA diagnostics (Conjugate NG):\n")
mcmc_conj_ng <- coda::mcmc(cbind(theta_out_ng,
                                mu = mu_out_ng,
                                sigma_theta = sigma_theta_out_ng))
print(summary(mcmc_conj_ng))

cat("\nEffective sample size:\n")
print(coda::effectiveSize(mcmc_conj_ng))

cat("\nGeweke z-scores (first 10 parameters shown):\n")
gwe <- coda::geweke.diag(mcmc_conj_ng)$z
print(utils::head(gwe, 10))

cat("\nLag-1 autocorrelation (mu, sigma_theta):\n")
print(c(mu = coda::autocorr(coda::mcmc(mu_out_ng), lag = 1),
        sigma_theta = coda::autocorr(coda::mcmc(sigma_theta_out_ng), lag = 1)))


readline("\n--- End of Conjugate NG sampler block. Press <Enter> to continue... ---")


######################################## Normal prior (fixed sigma_theta^2) ###############

## Data
school   <- c("A","B","C","D","E","F","G","H")
estimate <- c(28.39, 7.94, -2.75, 6.82, -0.64, 0.63, 18.01, 12.16)
sd_obs   <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)

J          <- length(school)
y          <- estimate
sigma_y_sq <- sd_obs^2

## Prior: mu ~ N(mu_mu, sigma_mu), theta_j | mu ~ N(mu, sigma_theta^2)
mu_mu    <- mean(estimate)
sigma_mu <- var(estimate)

## Fixed group-level variance from independent NG sampler
sigma_theta_fixed <- mean(sigma_theta_out^2)

## MCMC settings
n_burn_n <- 1000
n_sim_n  <- 1000

theta_n      <- y
mu_n         <- mu_mu
theta_out_n  <- matrix(0, nrow = n_sim_n, ncol = J)
mu_out_n     <- numeric(n_sim_n)
iters_out1_n <- numeric(n_burn_n)
iters_out2_n <- numeric(n_sim_n)

set.seed(104)
x_one <- matrix(1, nrow = J, ncol = 1)

## -----------------------------
## Burn-in (timed)
## -----------------------------
burn_time_n <- system.time({
  for (k in seq_len(n_burn_n)) {
    
    ## Population block: update mu | theta
    out_pop <- rlmb(
      1,
      y       = theta_n,
      x       = x_one,
      pfamily = dNormal(mu_mu, Sigma = sigma_mu, dispersion = sigma_theta_fixed)
    )
    mu_n <- out_pop$coefficients[1, 1]
    
    ## School-level block: update theta_j | mu, y_j
    for (j in seq_len(J)) {
      theta_n[j] <- rlmb(
        1,
        y       = y[j],
        x       = matrix(1),
        pfamily = dNormal(mu_n, Sigma = sigma_theta_fixed, dispersion = sigma_y_sq[j])
      )$coefficients[1, 1]
    }
    
    iters_out1_n[k] <- out_pop$iters
  }
})

## -----------------------------
## Main simulation (timed)
## -----------------------------
sim_time_n <- system.time({
  for (k in seq_len(n_sim_n)) {
    
    ## Population block
    out_pop <- rlmb(
      1,
      y       = theta_n,
      x       = x_one,
      pfamily = dNormal(mu_mu, Sigma = sigma_mu, dispersion = sigma_theta_fixed)
    )
    mu_n <- out_pop$coefficients[1, 1]
    
    ## School-level block
    for (j in seq_len(J)) {
      theta_n[j] <- rlmb(
        1,
        y       = y[j],
        x       = matrix(1),
        pfamily = dNormal(mu_n, Sigma = sigma_theta_fixed, dispersion = sigma_y_sq[j])
      )$coefficients[1, 1]
    }
    
    theta_out_n[k, ] <- theta_n
    mu_out_n[k]      <- mu_n
    iters_out2_n[k]  <- out_pop$iters
  }
})

## Print timing
cat("Fixed-sigma Normal burn-in time:\n")
print(burn_time_n)

cat("\nFixed-sigma Normal simulation time:\n")
print(sim_time_n)

## Summaries
colMeans(theta_out_n)
mean(mu_out_n)
rep(sqrt(sigma_theta_fixed), J)
sqrt(diag(var(theta_out_n)))

## -----------------------------
## Acceptance diagnostics (Fixed-sigma Normal)
## -----------------------------
mean_iters_out1_n <- mean(iters_out1_n)
cat("\nMean iters_out1_n (burn-in; 1 / acceptance rate):\n")
print(mean_iters_out1_n)
cat("\nEstimated burn-in acceptance rate:\n")
print(1 / mean_iters_out1_n)

mean_iters_out2_n <- mean(iters_out2_n)
cat("\nMean iters_out2_n (1 / acceptance rate):\n")
print(mean_iters_out2_n)
cat("\nEstimated acceptance rate:\n")
print(1 / mean_iters_out2_n)

cat("\nCODA diagnostics (Fixed-sigma Normal):\n")
mcmc_fix_n <- coda::mcmc(cbind(theta_out_n, mu = mu_out_n))
print(summary(mcmc_fix_n))

cat("\nEffective sample size:\n")
print(coda::effectiveSize(mcmc_fix_n))

cat("\nGeweke z-scores (first 10 parameters shown):\n")
gwe <- coda::geweke.diag(mcmc_fix_n)$z
print(utils::head(gwe, 10))

cat("\nLag-1 autocorrelation (mu):\n")
print(coda::autocorr(coda::mcmc(mu_out_n), lag = 1))


