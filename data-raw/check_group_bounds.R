devtools::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)

n_j <- c("4" = 14, "10" = 16, "20" = 15, "30" = 15, "43" = 15)
sigma2_hat_j <- c("4" = 157.67, "10" = 428.47, "20" = 316.08, "30" = 579.97, "43" = 352.89)
p_re <- 2
pwt <- 0.2

for (lev in names(n_j)) {
  design_j <- list(re_coef_names = c("a", "b"), y = seq_len(n_j[[lev]]))
  n_prior_j <- pwt / (1 - pwt) * n_j[[lev]]
  gp <- glmbayesCore:::.lmebayes_calibrate_ing_prior_measurement(
    design_j, sigma2_hat_j[[lev]], n_prior_j, 0.99
  )
  cat(sprintf(
    "school %s: shape=%.3f rate=%.1f disp_lower=%.3f disp_upper=%.1f ratio=%.1f\n",
    lev, gp$shape, gp$rate, gp$disp_lower, gp$disp_upper, gp$disp_upper / gp$disp_lower
  ))
}
