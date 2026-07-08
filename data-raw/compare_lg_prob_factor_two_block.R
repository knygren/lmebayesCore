# Compare lg_prob_factor: legacy joint (81 faces) vs block per-block (9+9)
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

lgL <- sim_env$UB_list$lg_prob_factor
cbL <- sim_env$Envelope$cbars
ub2L <- sim_env$UB_list$UB2min

de <- sim_block$build_out$dispersion_envelope
bd <- de$block_dispersion
lgB1 <- bd[[1]]$lg_prob_factor
lgB2 <- bd[[2]]$lg_prob_factor
cbB1 <- sim_block$build_out$block_envelopes[[1]]$cbars
cbB2 <- sim_block$build_out$block_envelopes[[2]]$cbars
ub2B1 <- bd[[1]]$UB2min
ub2B2 <- bd[[2]]$UB2min

cat("=== Storage layout ===\n")
cat(sprintf(
  "Legacy: lg_prob_factor length %d (= gs joint faces); also UB2min length %d\n",
  length(lgL), length(ub2L)
))
cat(sprintf(
  "        in UB_list (post-EnvelopeSort, aligned with Env$cbars / PLSD)\n"
))
cat(sprintf(
  "Block B1: lg_prob_factor length %d; B2 length %d (per-block EDB, own max_upp)\n",
  length(lgB1), length(lgB2)
))
cat("\nPer-block lg_prob_factor summaries:\n")
print(summary(lgB1))
cat("B2 identical to B1 on this duplicate-data fixture:",
    isTRUE(all.equal(lgB1, lgB2)), "\n")
cat("\nLegacy joint lg_prob_factor summary:\n")
print(summary(lgL))

# Match joint faces by concatenated cbars (B1 cols || B2 cols)
gs1 <- nrow(cbB1)
gs2 <- nrow(cbB2)
prod_sum <- numeric(gs1 * gs2)
prod_cb <- matrix(NA_real_, gs1 * gs2, 2L * l1)
idx <- 0L
for (j1 in seq_len(gs1)) {
  for (j2 in seq_len(gs2)) {
    idx <- idx + 1L
    prod_sum[idx] <- lgB1[j1] + lgB2[j2]
    prod_cb[idx, ] <- c(cbB1[j1, ], cbB2[j2, ])
  }
}

match_legacy <- function(cbars_target, cbars_pool, tol = 1e-8) {
  vapply(seq_len(nrow(cbars_pool)), function(i) {
    max(abs(cbars_pool[i, ] - cbars_target)) < tol
  }, logical(1))
}

legacy_idx <- integer(gs1 * gs2)
for (k in seq_along(prod_sum)) {
  hit <- which(match_legacy(prod_cb[k, ], cbL))
  if (length(hit) != 1L) {
    legacy_idx[k] <- NA_integer_
  } else {
    legacy_idx[k] <- hit
  }
}

n_match <- sum(!is.na(legacy_idx))
cat(sprintf(
  "\n=== Product-face matching via cbars (gs=%d x gs=%d -> %d) ===\n",
  gs1, gs2, gs1 * gs2
))
cat(sprintf("Matched %d / %d product faces to a unique legacy face\n", n_match, gs1 * gs2))

if (n_match == gs1 * gs2) {
  lg_legacy_at_prod <- lgL[legacy_idx]
  diff <- prod_sum - lg_legacy_at_prod
  cat("\nBlock sum lg_prob_factor(j1)+lg_prob_factor(j2) minus legacy joint:\n")
  print(summary(diff))
  cat(sprintf(
    "max abs diff = %.6g  RMSE = %.6g\n",
    max(abs(diff)), sqrt(mean(diff^2))
  ))
  worst <- order(-abs(diff))[seq_len(min(5L, length(diff)))]
  cat("\nLargest |diff| product faces (j1,j2):\n")
  for (w in worst) {
    j1 <- ((w - 1L) %/% gs2) + 1L
    j2 <- ((w - 1L) %% gs2) + 1L
    cat(sprintf(
      "  (j1=%d,j2=%d): block_sum=%.6f legacy=%.6f diff=%+.6f\n",
      j1, j2, prod_sum[w], lg_legacy_at_prod[w], diff[w]
    ))
  }

  ub2_sum <- outer(ub2B1, ub2B2, "+")[cbind(
    ((seq_along(prod_sum) - 1L) %/% gs2) + 1L,
    ((seq_along(prod_sum) - 1L) %% gs2) + 1L
  )]
  # outer gives matrix; flatten same odometer order
  ub2_prod <- numeric(gs1 * gs2)
  idx <- 0L
  for (j1 in seq_len(gs1)) {
    for (j2 in seq_len(gs2)) {
      idx <- idx + 1L
      ub2_prod[idx] <- ub2B1[j1] + ub2B2[j2]
    }
  }
  ub2_legacy_at_prod <- ub2L[legacy_idx]
  diff2 <- ub2_prod - ub2_legacy_at_prod
  cat("\nUB2min sum (block) minus legacy joint:\n")
  print(summary(diff2))
  cat(sprintf("max abs diff = %.6g\n", max(abs(diff2))))
}

cat("\nLegacy lg_prob_factor range vs block single-block range:\n")
cat(sprintf("  legacy [%.4f, %.4f]  block [%.4f, %.4f]\n",
            min(lgL), max(lgL), min(lgB1), max(lgB1)))

prod_sum <- as.vector(outer(lgB1, lgB2, "+"))
cat("\n=== If block UB3A used lg_B1(j1)+lg_B2(j2) on product faces ===\n")
cat("Product-sum lg_prob_factor (81 values):\n")
print(summary(prod_sum))
cat(sprintf(
  "Mean product-sum = %.4f vs legacy joint mean = %.4f (diff %+.4f)\n",
  mean(prod_sum), mean(lgL), mean(prod_sum) - mean(lgL)
))
cat(sprintf(
  "Sorted-vector correlation(product-sum, legacy) = %.4f\n",
  cor(sort(prod_sum), sort(lgL))
))
cat(
  "\nNote: cbars live in different standardizations (joint 4-D vs per-block 2-D),\n",
  "so face indices are not directly alignable by cbars concat on this fixture.\n",
  "Legacy prob_factor uses joint max_upp; block uses per-block max_upp in EDB.\n"
)
