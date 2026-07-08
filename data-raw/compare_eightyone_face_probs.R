# Compare all 81 product-face probabilities: legacy joint vs block (independence product)
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

logp_vec <- function(logP) {
  if (is.matrix(logP) || is.data.frame(logP)) {
    as.numeric(logP[, 1L])
  } else {
    as.numeric(logP)
  }
}

match_legacy_by_cbars <- function(cbars_target, cbars_pool, tol = 1e-8) {
  vapply(seq_len(nrow(cbars_pool)), function(i) {
    max(abs(cbars_pool[i, ] - cbars_target)) < tol
  }, logical(1))
}

## Legacy joint (81 faces, 4-D envelope)
plsdL <- sim_env$Envelope$PLSD
logpL <- logp_vec(sim_env$Envelope$logP)
cbL <- sim_env$Envelope$cbars
lgL <- sim_env$UB_list$lg_prob_factor
ub2L <- sim_env$UB_list$UB2min

## Block: per-block EDB envelopes (9 faces each); sim draws J1, J2 independently
env1 <- sim_block$build_out$block_envelopes[[1]]
env2 <- sim_block$build_out$block_envelopes[[2]]
plsd1 <- env1$PLSD
plsd2 <- env2$PLSD
logp1 <- logp_vec(env1$logP)
logp2 <- logp_vec(env2$logP)
cb1 <- env1$cbars
cb2 <- env2$cbars
bd1 <- sim_block$build_out$dispersion_envelope$block_dispersion[[1]]
bd2 <- sim_block$build_out$dispersion_envelope$block_dispersion[[2]]
lg1 <- bd1$lg_prob_factor
lg2 <- bd2$lg_prob_factor
ub2_1 <- bd1$UB2min
ub2_2 <- bd2$UB2min
joint_lg <- sim_block$build_out$dispersion_envelope$joint_lg_prob_factor
joint_ub2 <- sim_block$build_out$dispersion_envelope$joint_ub2min_product

gs1 <- length(plsd1)
gs2 <- length(plsd2)
stopifnot(gs1 * gs2 == length(plsdL))

# Odometer: last block (B2) varies fastest -> flat = j1 * gs2 + j2 (0-based: j1*gs2+j2)
product_flat <- function(j1, j2, gs2) (j1 - 1L) * gs2 + j2

plsd_prod <- numeric(gs1 * gs2)
logp_prod <- numeric(gs1 * gs2)
lg_prod_indep <- numeric(gs1 * gs2)
ub2_prod_indep <- numeric(gs1 * gs2)
logw_prod_indep <- numeric(gs1 * gs2)
logw_joint_product <- numeric(gs1 * gs2)
norm2_prod <- numeric(gs1 * gs2)
face_j1 <- integer(gs1 * gs2)
face_j2 <- integer(gs1 * gs2)

softmax <- function(x) {
  z <- exp(x - max(x))
  z / sum(z)
}

norm2_row <- function(cbars, row) sum(cbars[row, ]^2)

legacy_idx <- integer(gs1 * gs2)
prod_cb <- matrix(NA_real_, gs1 * gs2, ncol(cbL))

for (j1 in seq_len(gs1)) {
  for (j2 in seq_len(gs2)) {
    k <- product_flat(j1, j2, gs2)
    plsd_prod[k] <- plsd1[j1] * plsd2[j2]
    logp_prod[k] <- logp1[j1] + logp2[j2]
    lg_prod_indep[k] <- lg1[j1] + lg2[j2]
    ub2_prod_indep[k] <- ub2_1[j1] + ub2_2[j2]  # old wrong sum-of-minima
    face_j1[k] <- j1
    face_j2[k] <- j2
    prod_cb[k, ] <- c(cb1[j1, ], cb2[j2, ])

    hit <- which(match_legacy_by_cbars(prod_cb[k, ], cbL))
    legacy_idx[k] <- if (length(hit) == 1L) hit else NA_integer_

    # EDB log-weight: logP + 0.5||c||^2 + (lg - UB2min); PLSD = softmax(logw)
    n2_1 <- norm2_row(cb1, j1)
    n2_2 <- norm2_row(cb2, j2)
    logw1 <- logp1[j1] + 0.5 * n2_1 + (lg1[j1] - ub2_1[j1])
    logw2 <- logp2[j2] + 0.5 * n2_2 + (lg2[j2] - ub2_2[j2])
    logw_prod_indep[k] <- logw1 + logw2
    norm2_prod[k] <- n2_1 + n2_2
    logw_joint_product[k] <- logp_prod[k] + 0.5 * norm2_prod[k] +
      (joint_lg[k] - joint_ub2[k])
  }
}

plsd_joint_on_product <- softmax(logw_joint_product)

# Legacy logw at each legacy face index
logw_legacy <- numeric(length(plsdL))
for (k in seq_along(plsdL)) {
  n2 <- norm2_row(cbL, k)
  logw_legacy[k] <- logpL[k] + 0.5 * n2 + (lgL[k] - ub2L[k])
}

n_match <- sum(!is.na(legacy_idx))
cat(sprintf(
  "Matched %d / %d product faces to legacy via cbars concat\n",
  n_match, gs1 * gs2
))

# cbars-aligned vectors (product face k -> legacy face legacy_idx[k])
plsdL_at_prod <- plsdL[legacy_idx]
logw_legacy_at_prod <- logw_legacy[legacy_idx]
lgL_at_prod <- lgL[legacy_idx]
ub2L_at_prod <- ub2L[legacy_idx]

cat("=== Face-count / normalization ===\n")
cat(sprintf("Legacy PLSD sum = %.10f  (n = %d)\n", sum(plsdL), length(plsdL)))
cat(sprintf("Block PLSD B1 sum = %.10f  B2 sum = %.10f\n", sum(plsd1), sum(plsd2)))
cat(sprintf("Product PLSD sum  = %.10f  (independence: should be 1)\n", sum(plsd_prod)))

cat("\n=== Three proposal laws on the same 81 (j1,j2) product faces ===\n")
cat("  (1) indep: PLSD_B1[j1] * PLSD_B2[j2]  — block sim face draw today\n")
cat("  (2) joint-on-product: softmax(logP_sum + 0.5||c||^2 + joint_lg - joint_ub2)\n")
cat("      same per-block logP/cbars; joint UB3A/UB2 slack (accept tables)\n")
cat("  (3) legacy: PLSD from single 4-D envelope (different face geometry)\n")
diff_indep_joint <- plsd_prod - plsd_joint_on_product
cat(sprintf(
  "  TV(indep vs joint-on-product) = %.6f\n",
  0.5 * sum(abs(diff_indep_joint))
))
cat(sprintf(
  "  max |PLSD_indep - PLSD_joint_on_product| = %.6f  cor = %.5f\n",
  max(abs(diff_indep_joint)), cor(plsd_prod, plsd_joint_on_product)
))
cat(sprintf(
  "  TV(legacy vs joint-on-product) sorted-index = %.6f  cor(sort) = %.5f\n",
  0.5 * sum(abs(sort(plsd_joint_on_product) - sort(plsdL))),
  cor(sort(plsd_joint_on_product), sort(plsdL))
))

# Optimal face matching: assign each product face to a unique legacy face
# by minimum |logw_joint_product - logw_legacy| (same weight formula, different geometry)
if (requireNamespace("clue", quietly = TRUE)) {
  cost <- outer(logw_joint_product, logw_legacy, FUN = function(a, b) abs(a - b))
  assign <- clue::solve_LSAP(cost)
  plsdL_matched <- plsdL[assign]
  logwL_matched <- logw_legacy[assign]
  diff_plsd_matched <- plsd_joint_on_product - plsdL_matched
  cat(sprintf(
    "\n  Optimal logw-matching (clue::solve_LSAP): TV(joint-on-product vs legacy) = %.6f\n",
    0.5 * sum(abs(diff_plsd_matched))
  ))
  cat(sprintf(
    "    cor(PLSD) = %.5f  max |diff| = %.6f\n",
    cor(plsd_joint_on_product, plsdL_matched), max(abs(diff_plsd_matched))
  ))
  cat(sprintf(
    "    TV(indep vs legacy, same matching) = %.6f\n",
    0.5 * sum(abs(plsd_prod - plsdL_matched))
  ))
  use_matched_legacy <- TRUE
} else {
  assign <- order(order(logw_joint_product))  # fallback: sort-rank alignment
  plsdL_matched <- plsdL[assign]
  diff_plsd_matched <- plsd_joint_on_product - plsdL_matched
  use_matched_legacy <- FALSE
  cat("\n  (Install 'clue' for optimal logw matching; using sort-rank fallback)\n")
}

cat("\n=== PLSD: legacy joint vs block independence product ===\n")
if (n_match == gs1 * gs2) {
  diff_plsd <- plsd_prod - plsdL_at_prod
  cat("  (cbars-matched: same geometric face, different proposal weights)\n")
  cat(sprintf(
    "  mean diff (prod - legacy) = %+.6g\n  max |diff| = %.6g\n  cor = %.5f\n",
    mean(diff_plsd), max(abs(diff_plsd)), cor(plsd_prod, plsdL_at_prod)
  ))
  cat(sprintf(
    "  TV distance = %.6f (0.5 * sum |p_prod - p_legacy|)\n",
    0.5 * sum(abs(diff_plsd))
  ))
  cat(sprintf(
    "  legacy max PLSD = %.6f  product max = %.6f\n",
    max(plsdL_at_prod), max(plsd_prod)
  ))
} else {
  diff_plsd <- plsd_prod - plsdL[seq_along(plsd_prod)]
  cat("  (WARNING: cbars match incomplete; comparing raw index 1:81)\n")
  cat(sprintf(
    "  mean diff = %+.6g  max |diff| = %.6g  cor(sort) = %.5f\n",
    mean(diff_plsd), max(abs(diff_plsd)), cor(sort(plsd_prod), sort(plsdL))
  ))
  cat(sprintf(
    "  TV distance (raw index) = %.6f\n",
    0.5 * sum(abs(diff_plsd))
  ))
}

cat("\n=== logP (log envelope mass): sum of block logP vs legacy ===\n")
if (n_match == gs1 * gs2) {
  logpL_at_prod <- logpL[legacy_idx]
  cat(sprintf(
    "  logp_prod vs logpL (matched): mean diff = %+.4f  max |diff| = %.4f\n",
    mean(logp_prod - logpL_at_prod), max(abs(logp_prod - logpL_at_prod))
  ))
} else {
  cat(sprintf(
    "  logp_prod vs logpL: cor(sort) = %.5f  max |diff| on sorted = %.4f\n",
    cor(sort(logp_prod), sort(logpL)),
    max(abs(sort(logp_prod) - sort(logpL)))
  ))
}

cat("\n=== lg_prob_factor (UB3A slack): three block constructions ===\n")
if (n_match == gs1 * gs2) {
  cat("  (a) sum per-block lg (independence anchors)\n")
  cat(sprintf(
    "    mean diff vs legacy = %+.4f  max |diff| = %.4f\n",
    mean(lg_prod_indep - lgL_at_prod), max(abs(lg_prod_indep - lgL_at_prod))
  ))
  cat("  (b) block joint product table (joint anchors, used at sim accept)\n")
  cat(sprintf(
    "    mean diff vs legacy = %+.4f  max |diff| = %.4f\n",
    mean(joint_lg - lgL_at_prod), max(abs(joint_lg - lgL_at_prod))
  ))
} else {
  cat(sprintf(
    "  sum lg: cor(sort) = %.5f; joint lg: cor(sort) = %.5f\n",
    cor(sort(lg_prod_indep), sort(lgL)), cor(sort(joint_lg), sort(lgL))
  ))
}

cat("\n=== UB2min: sum-of-block-min vs joint-min vs legacy ===\n")
ub2_sum_bug <- numeric(gs1 * gs2)
for (j1 in seq_len(gs1)) for (j2 in seq_len(gs2)) {
  k <- product_flat(j1, j2, gs2)
  ub2_sum_bug[k] <- ub2_1[j1] + ub2_2[j2]
}
cat(sprintf(
  "  sum UB2min (bug): mean diff = %+.4f  max |diff| = %.4f\n",
  mean(ub2_sum_bug - ub2L), max(abs(ub2_sum_bug - ub2L))
))
if (n_match == gs1 * gs2) {
  cat(sprintf(
    "  joint_ub2 table:  mean diff = %+.4f  max |diff| = %.4f\n",
    mean(joint_ub2 - ub2L_at_prod), max(abs(joint_ub2 - ub2L_at_prod))
  ))
} else {
  cat(sprintf(
    "  joint_ub2 table:  cor(sort) = %.5f\n",
    cor(sort(joint_ub2), sort(ub2L))
  ))
}

cat("\n=== Unnormalized log-weights (logP + 0.5||c||^2 + lg - UB2min) ===\n")
if (n_match == gs1 * gs2) {
  cat(sprintf(
    "  indep product logw vs legacy logw (cbars-matched):\n"
  ))
  cat(sprintf(
    "    mean diff = %+.4f  max |diff| = %.4f  cor = %.5f\n",
    mean(logw_prod_indep - logw_legacy_at_prod),
    max(abs(logw_prod_indep - logw_legacy_at_prod)),
    cor(logw_prod_indep, logw_legacy_at_prod)
  ))
  lw <- logw_prod_indep - max(logw_prod_indep)
  w_indep <- exp(lw)
  plsd_from_logw <- w_indep / sum(w_indep)
  cat(sprintf(
    "  PLSD from indep logw renormalize vs legacy: TV = %.6f  cor = %.5f\n",
    0.5 * sum(abs(plsd_from_logw - plsdL_at_prod)), cor(plsd_from_logw, plsdL_at_prod)
  ))
} else {
  cat(sprintf(
    "    mean diff (raw index) = %+.4f  cor(sort) = %.5f\n",
    mean(logw_prod_indep - logw_legacy[seq_along(logw_prod_indep)]),
    cor(sort(logw_prod_indep), sort(logw_legacy))
  ))
  lw <- logw_prod_indep - max(logw_prod_indep)
  plsd_from_logw <- w_indep <- exp(lw) / sum(exp(lw))
}
cat(sprintf(
  "  PLSD from indep logw vs product PLSD1*PLSD2: max |diff| = %.2e (should ~0)\n",
  max(abs(plsd_from_logw - plsd_prod))
))

cat("\n=== All 81 product faces (j1,j2): proposal probabilities ===\n")
face_tbl <- data.frame(
  j1 = face_j1,
  j2 = face_j2,
  PLSD_indep = plsd_prod,
  PLSD_joint_on_product = plsd_joint_on_product,
  diff_indep_joint = diff_indep_joint,
  PLSD_legacy_matched = plsdL_matched,
  diff_joint_legacy = plsd_joint_on_product - plsdL_matched
)
face_tbl <- face_tbl[order(-abs(face_tbl$diff_indep_joint)), ]
print(face_tbl, row.names = FALSE, digits = 5)

cat("\n=== Marginal over block-1 face (sum over j2) ===\n")
marg_indep <- rowsum(plsd_prod, face_j1, reorder = FALSE)
marg_joint <- rowsum(plsd_joint_on_product, face_j1, reorder = FALSE)
marg_legacy <- rowsum(plsdL_matched, face_j1, reorder = FALSE)
marg_tbl <- data.frame(
  j1 = seq_len(gs1),
  PLSD_B1 = plsd1,
  marginal_indep = as.numeric(marg_indep),
  marginal_joint = as.numeric(marg_joint),
  marginal_legacy_matched = as.numeric(marg_legacy)
)
print(marg_tbl, row.names = FALSE, digits = 5)
cat(sprintf(
  "  Note: marginal_indep equals PLSD_B1 exactly (max |diff| = %.2e)\n",
  max(abs(marg_tbl$marginal_indep - marg_tbl$PLSD_B1))
))

cat("\n=== All 81 faces: PLSD legacy vs product (cbars-matched) ===\n")
if (n_match == gs1 * gs2) {
  face_tbl <- data.frame(
    j1 = face_j1,
    j2 = face_j2,
    PLSD_legacy = plsdL_at_prod,
    PLSD_product = plsd_prod,
    diff = diff_plsd,
    ratio = plsd_prod / plsdL_at_prod
  )
  face_tbl <- face_tbl[order(-abs(face_tbl$diff)), ]
  print(head(face_tbl, 12L), row.names = FALSE)
  cat(sprintf(
    "\n  Summary ratio (product/legacy): mean=%.4f  min=%.4f  max=%.4f\n",
    mean(face_tbl$ratio), min(face_tbl$ratio), max(face_tbl$ratio)
  ))
} else {
  cat("  (skipped — cbars matching incomplete)\n")
}

cat("\n=== Worst PLSD mismatches (product - legacy) ===\n")
worst <- order(-abs(diff_plsd))[seq_len(min(8L, length(diff_plsd)))]
for (w in worst) {
  leg <- if (n_match == gs1 * gs2) plsdL_at_prod[w] else plsdL[w]
  cat(sprintf(
    "  face (j1=%d,j2=%d): legacy=%.6f  product=%.6f  diff=%+.6f\n",
    face_j1[w], face_j2[w], leg, plsd_prod[w], diff_plsd[w]
  ))
}

cat("\n=== Interpretation ===\n")
cat(
  "Block sim (k>1) draws one flat index from build_out$dispersion_envelope$joint_PLSD,\n",
  "then decodes (j1,j2). PLSD_indep = PLSD_B1*PLSD_B2 is no longer used at sim.\n",
  "joint_PLSD matches joint-on-product softmax (legacy weight multiset, product indexing).\n"
)
