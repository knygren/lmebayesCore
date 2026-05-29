# Benchmark the BikeSharing two-block Gibbs sampler (Chapter 18 / Ex_09).
# Block 1: rglmb on population (theta ~ X_train).
# Block 2: rNormalGLM_reg_block() (new) or per-observation rglmb loop (legacy).
#
#   Rscript data-raw/benchmark_BikeSharing_rNormalGLM_reg_block.R
#   Rscript data-raw/benchmark_BikeSharing_rNormalGLM_reg_block.R quick
#   Rscript data-raw/benchmark_BikeSharing_rNormalGLM_reg_block.R legacy
#   Rscript data-raw/benchmark_BikeSharing_rNormalGLM_reg_block.R quick legacy
#
# Default: n_burn = 200, n_sim = 1000 (same as demo("Ex_09_BikeSharingPoisson")).
# Append `quick` for n_burn = 5, n_sim = 10 (smoke test only).
#
# Does not modify inst/extdata/BikeSharing_ch14_gibbs.rds.

args <- commandArgs(trailingOnly = TRUE)
run_legacy <- any(tolower(args) %in% c("legacy", "--legacy", "-l"))
run_quick <- any(tolower(args) %in% c("quick", "--quick", "-q"))
path_args <- args[!tolower(args) %in% c(
  "legacy", "--legacy", "-l",
  "quick", "--quick", "-q"
)]

root <- if (length(path_args) >= 1L) {
  normalizePath(path_args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload, e.g. install.packages('pkgload')")
}

pkgload::load_all(export_all = FALSE)

summ_time <- function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
}

fmt_hms <- function(secs) {
  secs <- as.numeric(secs)
  if (!is.finite(secs) || secs < 0) {
    secs <- 0
  }
  h <- floor(secs / 3600)
  rem <- secs - h * 3600
  m <- floor(rem / 60)
  s <- rem - m * 60
  sprintf("%d h %d min %.2f s", h, m, s)
}

make_mcmc_main <- function(beta_out, sigma_out, beta_names) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Install coda for coefficient summaries, e.g. install.packages('coda')")
  }
  mcmc <- coda::mcmc(cbind(beta_out, sigma_theta = sigma_out))
  colnames(mcmc) <- c(beta_names, "sigma_theta")
  mcmc
}

summarize_gibbs_coefficients <- function(beta_out, sigma_out, beta_names, label = "") {
  mcmc_main <- make_mcmc_main(beta_out, sigma_out, beta_names)
  if (nzchar(label)) {
    message("\n--- CODA summary (", label, ") ---")
  } else {
    message("\n--- CODA summary (beta + sigma_theta) ---")
  }
  print(summary(mcmc_main))
  message("\n--- posterior means ---")
  print(colMeans(as.matrix(mcmc_main)))
  if (nrow(beta_out) >= 50L) {
    es <- coda::effectiveSize(mcmc_main)
    message("\n--- effective sample size ---")
    print(es)
  } else {
    message("(skipped effectiveSize: n_sim < 50)")
  }
  invisible(mcmc_main)
}

compare_coefficients_to_vignette <- function(beta_out, sigma_out, beta_names) {
  ch14_path <- file.path(root, "inst", "extdata", "BikeSharing_ch14_gibbs.rds")
  if (!file.exists(ch14_path)) {
    ch14_path <- system.file("extdata", "BikeSharing_ch14_gibbs.rds", package = "glmbayes")
  }
  if (!nzchar(ch14_path) || !file.exists(ch14_path)) {
    message("Vignette reference BikeSharing_ch14_gibbs.rds not found; skip comparison.")
    return(invisible(NULL))
  }
  ch14 <- readRDS(ch14_path)
  mcmc_new <- make_mcmc_main(beta_out, sigma_out, beta_names)
  mcmc_vig <- ch14$mcmc_main
  stopifnot(identical(colnames(mcmc_new), colnames(mcmc_vig)))

  mean_new <- colMeans(as.matrix(mcmc_new))
  mean_vig <- colMeans(as.matrix(mcmc_vig))
  diff_mean <- mean_new - mean_vig

  message("\n--- compare posterior means vs vignette (legacy rglmb Block 2) ---")
  cmp <- data.frame(
    parameter = names(diff_mean),
    mean_new = mean_new,
    mean_vignette = mean_vig,
    diff = diff_mean,
    row.names = NULL
  )
  print(cmp, digits = 4, row.names = FALSE)
  message("max |diff| in posterior means: ", signif(max(abs(diff_mean)), 4))

  q_new <- apply(as.matrix(mcmc_new), 2, quantile, probs = c(0.025, 0.5, 0.975))
  q_vig <- apply(as.matrix(mcmc_vig), 2, quantile, probs = c(0.025, 0.5, 0.975))
  message("max |diff| in 2.5% / 50% / 97.5% quantiles: ",
          signif(max(abs(q_new - q_vig)), 4))

  invisible(list(compare_means = cmp, mcmc_new = mcmc_new, mcmc_vignette = mcmc_vig))
}

n_burn <- if (run_quick) 5L else 200L
n_sim  <- if (run_quick) 10L else 1000L
n_gibbs <- n_burn + n_sim

## --- Data / priors (same as demo Ex_09 and Chapter-18) -----------------------
data("BikeSharing", package = "glmbayes")

cont_vars <- c(
  "temp", "atemp", "hum", "windspeed",
  "hr_sin", "hr_cos", "mon_sin", "mon_cos"
)
BikeSharing_c <- BikeSharing
BikeSharing_c[cont_vars] <- scale(BikeSharing[cont_vars], center = TRUE, scale = FALSE)

form2 <- cnt ~ part_of_day + quarter + holiday + workingday + weathersit +
  hr_sin + hr_cos + mon_sin + mon_cos

pct_train <- 0.01
set.seed(42)
n <- nrow(BikeSharing_c)
idx_train <- sample(n, size = round(pct_train * n))

Bike_train <- BikeSharing_c[idx_train, ]
X_train <- model.matrix(form2, data = Bike_train)
y_train <- Bike_train$cnt
n_train <- length(y_train)
p <- ncol(X_train)

theta <- log(y_train + 0.5)
data_pop <- data.frame(theta = theta, Bike_train)
form_pop <- theta ~ part_of_day + quarter + holiday + workingday + weathersit +
  hr_sin + hr_cos + mon_sin + mon_cos
ps_pop <- Prior_Setup(form_pop, family = gaussian(), data = data_pop)

x_one <- matrix(1, n_train, 1)
colnames(x_one) <- "(Intercept)"

pfamily_pop <- dNormal_Gamma(
  ps_pop$mu, Sigma_0 = ps_pop$Sigma_0,
  ps_pop$shape, ps_pop$rate
)

message("BikeSharing train n = ", n_train, ", p = ", p)
message("Block Gibbs: n_burn = ", n_burn, ", n_sim = ", n_sim,
        " (", n_gibbs, " full iterations)")
if (run_quick) {
  message("Quick mode: not comparable to Chapter-18 / Ex_09 timings.")
}
if (!run_legacy) {
  message("Legacy full-Gibbs cross-check skipped (append 'legacy').")
}

block2_prior_lists <- function(mu_all, sigma_theta_sq) {
  lapply(mu_all, function(m) {
    list(
      mu = m,
      Sigma = matrix(sigma_theta_sq, 1, 1),
      dispersion = 1,
      ddef = FALSE
    )
  })
}

block1_update <- function(theta) {
  out_pop <- rglmb(
    1L, theta, X_train, family = gaussian(),
    pfamily = pfamily_pop,
    use_parallel = FALSE,
    verbose = FALSE
  )
  beta <- as.vector(out_pop$coefficients[1, ])
  sigma_theta_sq <- out_pop$dispersion[1]
  mu_all <- as.vector(X_train %*% beta)
  list(
    beta = beta,
    sigma_theta_sq = sigma_theta_sq,
    sigma_theta = sqrt(sigma_theta_sq),
    mu_all = mu_all
  )
}

block2_theta_reg_block <- function(mu_all, sigma_theta_sq, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- rNormalGLM_reg_block(
    n = 1L,
    y = y_train,
    x = x_one,
    block = seq_len(n_train),
    prior_lists = block2_prior_lists(mu_all, sigma_theta_sq),
    family = poisson(),
    Gridtype = 2L,
    n_envopt = 1L,
    use_parallel = FALSE,
    use_opencl = FALSE,
    verbose = FALSE,
    progbar = FALSE
  )
  as.vector(out$coefficients[, 1])
}

block2_theta_rglmb_loop <- function(mu_all, sigma_theta_sq, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  theta <- numeric(n_train)
  for (i in seq_len(n_train)) {
    theta[i] <- rglmb(
      1L, y = y_train[i], x = matrix(1, 1, 1),
      family = poisson(),
      pfamily = dNormal(mu = mu_all[i], Sigma = sigma_theta_sq),
      Gridtype = 2L,
      use_parallel = FALSE,
      use_opencl = FALSE,
      verbose = FALSE
    )$coefficients[1, 1]
  }
  theta
}

run_block_gibbs <- function(block2_fun, store = FALSE) {
  beta_out <- if (store) matrix(0, nrow = n_sim, ncol = p) else NULL
  sigma_out <- if (store) numeric(n_sim) else NULL
  theta_out <- if (store) matrix(0, nrow = n_sim, ncol = n_train) else NULL

  burn_time <- system.time({
    for (k in seq_len(n_burn)) {
      b1 <- block1_update(theta)
      theta <- block2_fun(b1$mu_all, b1$sigma_theta_sq)
    }
  })

  sim_time <- system.time({
    for (k in seq_len(n_sim)) {
      b1 <- block1_update(theta)
      theta <- block2_fun(b1$mu_all, b1$sigma_theta_sq)
      if (store) {
        beta_out[k, ] <- b1$beta
        sigma_out[k] <- b1$sigma_theta
        theta_out[k, ] <- theta
      }
    }
  })

  list(
    burn_time = burn_time,
    sim_time = sim_time,
    beta_out = beta_out,
    sigma_out = sigma_out,
    theta_out = theta_out,
    theta_final = theta
  )
}

# =============================================================================
# (1) Full two-block Gibbs — Block 2 via rNormalGLM_reg_block
# =============================================================================

message("\n========== two-block Gibbs (Block 2: rNormalGLM_reg_block) ==========")
message("Started: ", format(Sys.time(), usetz = TRUE))

set.seed(123)
gibbs_blk <- run_block_gibbs(block2_theta_reg_block, store = TRUE)

t_burn <- as.numeric(gibbs_blk$burn_time["elapsed"])
t_sim  <- as.numeric(gibbs_blk$sim_time["elapsed"])
t_total <- t_burn + t_sim
t_per_iter <- t_total / n_gibbs

message("\n--- timing: full Block Gibbs sampler ---")
message("burn-in (", n_burn, " iterations): ", fmt_hms(t_burn),
        " (", signif(t_burn, 4), " s)")
message("main (", n_sim, " iterations):     ", fmt_hms(t_sim),
        " (", signif(t_sim, 4), " s)")
message("TOTAL (", n_gibbs, " iterations):   ", fmt_hms(t_total),
        " (", signif(t_total, 4), " s)")
message("mean seconds per full Gibbs iteration: ", signif(t_per_iter, 4))

if (!run_quick) {
  message("\n--- final theta (after chain) ---")
  message("mean(theta) = ", signif(mean(gibbs_blk$theta_final), 4),
          ", sd(theta) = ", signif(sd(gibbs_blk$theta_final), 4))
}

beta_names <- colnames(X_train)
mcmc_blk <- summarize_gibbs_coefficients(
  gibbs_blk$beta_out,
  gibbs_blk$sigma_out,
  beta_names,
  label = "Block 2: rNormalGLM_reg_block"
)
cmp_vig <- compare_coefficients_to_vignette(
  gibbs_blk$beta_out,
  gibbs_blk$sigma_out,
  beta_names
)

benchmark <- list(
  n_train = n_train,
  p = p,
  n_burn = n_burn,
  n_sim = n_sim,
  gibbs_iterations = n_gibbs,
  quick_mode = run_quick,
  timing_seconds = list(
    burn_in = t_burn,
    main = t_sim,
    total = t_total,
    per_iteration_mean = t_per_iter
  ),
  block2_method = "rNormalGLM_reg_block",
  beta_names = beta_names,
  posterior_mean = colMeans(as.matrix(mcmc_blk)),
  compare_vignette_max_abs_mean_diff = if (!is.null(cmp_vig)) {
    max(abs(cmp_vig$compare_means$diff))
  } else {
    NA_real_
  },
  timestamp = Sys.time()
)

# =============================================================================
# (2) LEGACY: full Gibbs with vignette rglmb loop (+ optional checks)
# =============================================================================

if (run_legacy) {
  message("\n========== two-block Gibbs (Block 2: rglmb loop, legacy) ==========")
  message("Started: ", format(Sys.time(), usetz = TRUE))

  set.seed(123)
  gibbs_rglmb <- run_block_gibbs(block2_theta_rglmb_loop, store = FALSE)

  t_burn_l <- as.numeric(gibbs_rglmb$burn_time["elapsed"])
  t_sim_l  <- as.numeric(gibbs_rglmb$sim_time["elapsed"])
  t_total_l <- t_burn_l + t_sim_l

  message("\n--- timing: legacy full Block Gibbs ---")
  message("burn-in: ", fmt_hms(t_burn_l), " (", signif(t_burn_l, 4), " s)")
  message("main:    ", fmt_hms(t_sim_l), " (", signif(t_sim_l, 4), " s)")
  message("TOTAL:   ", fmt_hms(t_total_l), " (", signif(t_total_l, 4), " s)")
  message("speedup block vs legacy (total): ",
          signif(t_total_l / t_total, 3), "x")

  message("\n--- one Block-2 update (not full Gibbs), for reference ---")
  b1_ref <- block1_update(theta)
  n_time <- 3L
  time_blk <- numeric(n_time)
  time_rglmb <- numeric(n_time)
  for (t in seq_len(n_time)) {
    time_blk[t] <- system.time({
      block2_theta_reg_block(b1_ref$mu_all, b1_ref$sigma_theta_sq)
    })["elapsed"]
    time_rglmb[t] <- system.time({
      block2_theta_rglmb_loop(b1_ref$mu_all, b1_ref$sigma_theta_sq)
    })["elapsed"]
  }
  s_blk2 <- summ_time(time_blk)
  s_rglmb2 <- summ_time(time_rglmb)
  message("rNormalGLM_reg_block (one Block-2 step):")
  print(round(s_blk2, 3))
  message("rglmb loop (one Block-2 step):")
  print(round(s_rglmb2, 3))

  benchmark$legacy <- list(
    block2_method = "rglmb_loop",
    timing_seconds = list(
      burn_in = t_burn_l,
      main = t_sim_l,
      total = t_total_l,
      per_iteration_mean = t_total_l / n_gibbs,
      block2_only_rNormalGLM_reg_block = s_blk2,
      block2_only_rglmb_loop = s_rglmb2
    ),
    speedup_total_block_vs_legacy = t_total_l / t_total
  )
}

out_path <- file.path(root, "data-raw", "BikeSharing_block_reg_benchmark.rds")
saveRDS(benchmark, out_path)
message("\nWrote summary: ", out_path)
message("Done.")
