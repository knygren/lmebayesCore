# Helper: find best alignment key for product vs legacy faces
pkgload::load_all(quiet = TRUE)
n_envopt <- 10000L

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group1 <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight1 <- c(ctl, trt)
n1 <- length(weight1)
weight <- c(weight1, weight1)
group_stacked <- factor(c(group1, group1), levels = levels(group1))
p_setup <- Prior_Setup(
  weight ~ group,
  data = data.frame(weight = weight, group = group_stacked),
  family = gaussian()
)
y <- p_setup$y
x_block <- p_setup$x
l1 <- ncol(x_block)
block <- factor(c(rep("B1", n1), rep("B2", n1)), levels = c("B1", "B2"))
x1 <- x_block[seq_len(n1), , drop = FALSE]
zeros <- matrix(0, n1, l1)
x_old <- rbind(cbind(x1, zeros), cbind(zeros, x1))
mu1 <- p_setup$mu
Sigma1 <- p_setup$Sigma
prior_list_old <- list(
  mu = c(mu1, mu1),
  Sigma = rbind(
    cbind(Sigma1, matrix(0, l1, l1)),
    cbind(matrix(0, l1, l1), Sigma1)
  ),
  dispersion = p_setup$dispersion,
  shape = p_setup$shape,
  rate = p_setup$rate,
  Precision = solve(rbind(
    cbind(Sigma1, matrix(0, l1, l1)),
    cbind(matrix(0, l1, l1), Sigma1)
  )),
  max_disp_perc = 0.99
)
prior_list_block <- list(
  mu = mu1,
  Sigma = Sigma1,
  shape = p_setup$shape,
  rate = p_setup$rate,
  max_disp_perc = 0.99
)

sim_env <- rindepNormalGamma_reg_with_envelope(
  n = 1L, y = y, x = x_old, prior_list = prior_list_old,
  n_envopt = n_envopt, Gridtype = 3L, use_parallel = FALSE, progbar = FALSE
)
sim_block <- glmbayesCore:::.rIndepNormalGammaRegBlock_cpp(
  n = 1L, y = y, x = x_block, block = block,
  prior_list = prior_list_block, n_envopt = n_envopt, Gridtype = 3L,
  use_parallel = FALSE, progbar = FALSE, verbose = FALSE,
  offset = rep(0, length(y)), wt = rep(1, length(y)),
  p_re = -1L, n_rss_iter = 10L, RSS_ML = NA_real_, use_opencl = FALSE,
  group_levels = character(0), re_names = character(0)
)

envL <- sim_env$Envelope
env1 <- sim_block$build_out$block_envelopes[[1]]
env2 <- sim_block$build_out$block_envelopes[[2]]

match_row <- function(target, pool, tol = 1e-6) {
  which(apply(pool, 1, function(r) max(abs(r - target)) < tol))
}

gs1 <- nrow(env1$cbars)
gs2 <- nrow(env2$cbars)
flat <- function(j1, j2) (j1 - 1L) * gs2 + j2

keys <- list(
  cbars = matrix(NA, gs1 * gs2, 2L * l1),
  loglt = matrix(NA, gs1 * gs2, 2L * l1),
  logrt = matrix(NA, gs1 * gs2, 2L * l1),
  thetabars = matrix(NA, gs1 * gs2, 2L * l1)
)
for (j1 in seq_len(gs1)) for (j2 in seq_len(gs2)) {
  k <- flat(j1, j2)
  keys$cbars[k, ] <- c(env1$cbars[j1, ], env2$cbars[j2, ])
  keys$loglt[k, ] <- c(env1$loglt[j1, ], env2$loglt[j2, ])
  keys$logrt[k, ] <- c(env1$logrt[j1, ], env2$logrt[j2, ])
  keys$thetabars[k, ] <- c(env1$thetabars[j1, ], env2$thetabars[j2, ])
}

legacy_pools <- list(
  cbars = envL$cbars,
  loglt = envL$loglt,
  logrt = envL$logrt,
  thetabars = envL$thetabars
)

for (nm in names(keys)) {
  idx <- integer(gs1 * gs2)
  for (k in seq_len(gs1 * gs2)) {
    hit <- match_row(keys[[nm]][k, ], legacy_pools[[nm]])
    idx[k] <- if (length(hit) == 1L) hit else NA_integer_
  }
  cat(sprintf("%s: matched %d / %d\n", nm, sum(!is.na(idx)), length(idx)))
}

# If loglt/logrt match, compare PLSD
idx <- integer(gs1 * gs2)
for (k in seq_len(gs1 * gs2)) {
  hit <- match_row(keys$loglt[k, ], legacy_pools$loglt)
  if (length(hit) != 1L) hit <- match_row(keys$logrt[k, ], legacy_pools$logrt)
  idx[k] <- if (length(hit) == 1L) hit else NA_integer_
}
cat(sprintf("loglt primary match: %d / %d\n", sum(!is.na(idx)), length(idx)))

if (sum(!is.na(idx)) == gs1 * gs2) {
  plsdL <- envL$PLSD[idx]
  plsdP <- outer(env1$PLSD, env2$PLSD)[cbind(
    ((seq_len(81) - 1L) %/% gs2) + 1L,
    ((seq_len(81) - 1L) %% gs2) + 1L
  )]
  diff <- plsdP - plsdL
  cat(sprintf("PLSD TV (loglt-matched) = %.6f\n", 0.5 * sum(abs(diff))))
  cat(sprintf("PLSD max |diff| = %.6f  cor = %.5f\n", max(abs(diff)), cor(plsdP, plsdL)))
}
