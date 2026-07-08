# Product-face slack vs legacy after v2 fix
pkgload::load_all(quiet = TRUE)
source("data-raw/compare_lg_prob_factor_two_block.R", local = TRUE)

de <- sim_block$build_out$dispersion_envelope
lgB <- de$joint_lg_prob_factor
ub2B <- de$joint_ub2min_product

cat("\n=== Product-face slack vs legacy (n=", length(lgB), " faces) ===\n", sep = "")
diff_lg <- lgB - lgL
diff_ub2 <- ub2B - ub2L
cat("lg:  mean diff", mean(diff_lg), " max |diff|", max(abs(diff_lg)),
    " cor(sorted)", cor(sort(lgB), sort(lgL)), "\n")
cat("ub2: mean diff", mean(diff_ub2), " max |diff|", max(abs(diff_ub2)),
    " cor(sorted)", cor(sort(ub2B), sort(ub2L)), "\n")
cat("legacy ub2 range", range(ub2L), " block product ub2 range", range(ub2B), "\n")

# Old shortcut from per-block stored vectors (recompute)
bd1 <- bd[[1]]; bd2 <- bd[[2]]
max_upp <- de$prob_max_upp
max_low <- de$prob_max_low
gs <- 9
lg_old <- ub2_old <- numeric(81)
idx <- 0
for (j1 in 1:gs) for (j2 in 1:gs) {
  idx <- idx + 1
  upp_sum <- bd1$upp_apprx[j1] + bd2$upp_apprx[j2]
  low_sum <- bd1$low_apprx[j1] + bd2$low_apprx[j2]
  lg_old[idx] <- max(upp_sum - max_upp, low_sum - max_low)
  ub2_l <- bd1$ub2_at_low[j1] + bd2$ub2_at_low[j2]
  ub2_u <- bd1$ub2_at_upp[j1] + bd2$ub2_at_upp[j2]
  ub2_old[idx] <- min(ub2_l, ub2_u)
}
cat("\nTable vs per-draw shortcut recompute:\n")
cat("lg  max |diff|", max(abs(lgB - lg_old)), "\n")
cat("ub2 max |diff|", max(abs(ub2B - ub2_old)), "\n")
