## k = 2 stacked Dobson: compare legacy joint envelope artifacts
## (rindepNormalGamma_reg_with_envelope) vs block pipeline build_out.

.ing_legacy_ub_summary <- function(sim) {
  ub <- sim$UB_list
  gam <- sim$gamma_list
  c(
    lmc1 = ub$lmc1,
    lmc2 = ub$lmc2,
    lm_log1 = ub$lm_log1,
    lm_log2 = ub$lm_log2,
    max_New_LL_UB = ub$max_New_LL_UB,
    max_LL_log_disp = ub$max_LL_log_disp,
    RSS_Min = ub$RSS_Min,
    RSS_ML = ub$RSS_ML,
    shape3 = gam$shape3,
    rate2 = gam$rate2,
    disp_lower = gam$disp_lower,
    disp_upper = gam$disp_upper,
    n_faces = length(sim$Envelope$PLSD)
  )
}

.block_pooled_ub_summary <- function(sim_block) {
  de <- sim_block$build_out$dispersion_envelope
  ub <- de$UB_list
  gam <- de$gamma_list
  meta <- de$cross_face_meta
  list(
    numeric = c(
      lmc1 = ub$lmc1,
      lmc2 = ub$lmc2,
      lm_log1 = ub$lm_log1,
      lm_log2 = ub$lm_log2,
      max_New_LL_UB = ub$max_New_LL_UB,
      max_LL_log_disp = ub$max_LL_log_disp,
      RSS_Min = de$RSS_Min,
      RSS_ML = de$RSS_ML,
      shape3 = gam$shape3,
      rate2 = gam$rate2,
      disp_lower = gam$disp_lower,
      disp_upper = gam$disp_upper,
      n_identifiable = meta$n_identifiable
    ),
    aggregation = meta$aggregation
  )
}

.block_per_block_ub_rows <- function(sim_block) {
  de <- sim_block$build_out$dispersion_envelope
  blocks <- de$block_dispersion
  gs <- de$cross_face_meta$gs_per_block
  do.call(
    rbind,
    lapply(seq_along(blocks), function(j) {
      ub <- blocks[[j]]$UB_list
      gam <- blocks[[j]]$gamma_list
      data.frame(
        block_id = blocks[[j]]$block_id,
        gs = gs[[j]],
        lmc1 = ub$lmc1,
        lmc2 = ub$lmc2,
        max_New_LL_UB = ub$max_New_LL_UB,
        max_LL_log_disp = ub$max_LL_log_disp,
        RSS_Min = ub$RSS_Min,
        shape3 = gam$shape3,
        rate2 = gam$rate2,
        stringsAsFactors = FALSE
      )
    })
  )
}

.format_num_vec <- function(x) {
  paste(names(x), signif(as.numeric(x), 5), sep = " = ", collapse = "\n")
}

.legacy_vs_pooled_message <- function(legacy, pooled_num) {
  keys <- c(
    "disp_lower", "disp_upper", "RSS_Min", "rate2",
    "lmc1", "lmc2", "lm_log1", "lm_log2",
    "max_New_LL_UB", "max_LL_log_disp", "shape3"
  )
  lines <- vapply(keys, function(k) {
    sprintf(
      "%s: legacy %s | block pooled %s",
      k,
      signif(as.numeric(legacy[k]), 5),
      signif(as.numeric(pooled_num[k]), 5)
    )
  }, character(1))
  paste(lines, collapse = "\n")
}

.dobson_plant_two_block_fixture <- function() {
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group1 <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight1 <- c(ctl, trt)
  n1 <- length(weight1)

  weight <- c(weight1, weight1)
  group_stacked <- factor(c(group1, group1), levels = levels(group1))
  df_stacked <- data.frame(weight = weight, group = group_stacked)

  p_setup <- glmbayesCore::Prior_Setup(
    weight ~ group,
    data = df_stacked,
    family = gaussian()
  )

  y <- p_setup$y
  x_block <- p_setup$x
  l1 <- ncol(x_block)

  block <- factor(
    c(rep("B1", n1), rep("B2", n1)),
    levels = c("B1", "B2")
  )

  x1 <- x_block[seq_len(n1), , drop = FALSE]
  zeros <- matrix(0, n1, l1)
  x_old <- rbind(
    cbind(x1, zeros),
    cbind(zeros, x1)
  )

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

  list(
    y = y,
    x_block = x_block,
    x_old = x_old,
    block = block,
    prior_list_old = prior_list_old,
    prior_list_block = prior_list_block,
    l1 = l1,
    k = 2L,
    n_obs = length(y)
  )
}

.ing_rindepNormalGamma_reg_glmbayes_names <- function() {
  c(
    "coefficients", "coef.mode", "dispersion", "Prior", "family",
    "prior.weights", "y", "x", "call", "famfunc", "iters", "Envelope",
    "loglike", "weight_out", "sim_bounds", "pfamily", "offset2"
  )
}

test_that("two-block stacked Dobson: legacy vs block envelope / gamma / UB constants", {
  fix <- .dobson_plant_two_block_fixture()
  n_envopt <- 10000L
  n_draws <- 1L

  sim_std <- glmbayesCore::rindepNormalGamma_reg(
    n = n_draws,
    y = fix$y,
    x = fix$x_old,
    prior_list = fix$prior_list_old,
    n_envopt = n_envopt,
    Gridtype = 3L,
    use_parallel = FALSE,
    progbar = FALSE
  )
  expect_identical(names(sim_std), .ing_rindepNormalGamma_reg_glmbayes_names())
  expect_identical(class(sim_std), "rglmb")
  expect_null(sim_std$Envelope)
  expect_null(sim_std$loglike)
  expect_false("gamma_list" %in% names(sim_std))
  expect_false("UB_list" %in% names(sim_std))
  expect_false("diagnostics" %in% names(sim_std))

  sim_env <- rindepNormalGamma_reg_with_envelope(
    n = n_draws,
    y = fix$y,
    x = fix$x_old,
    prior_list = fix$prior_list_old,
    n_envopt = n_envopt,
    Gridtype = 3L,
    use_parallel = FALSE,
    progbar = FALSE
  )

  sim_block <- lmebayesCore:::.rIndepNormalGammaRegBlock_cpp(
    n = n_draws,
    y = fix$y,
    x = fix$x_block,
    block = fix$block,
    prior_list = fix$prior_list_block,
    prior_lists = NULL,
    offset = rep(0, fix$n_obs),
    wt = rep(1, fix$n_obs),
    p_re = -1L,
    n_rss_iter = 10L,
    Gridtype = 3L,
    n_envopt = n_envopt,
    RSS_ML = NA_real_,
    use_parallel = FALSE,
    use_opencl = FALSE,
    progbar = FALSE,
    verbose = FALSE,
    group_levels = character(0),
    re_names = character(0)
  )

  legacy <- .ing_legacy_ub_summary(sim_env)
  pooled <- .block_pooled_ub_summary(sim_block)
  pooled_num <- pooled$numeric
  per_block <- .block_per_block_ub_rows(sim_block)

  message("legacy joint envelope (p = ", 2L * fix$l1, ", one EnvelopeOrchestrator):\n",
          .format_num_vec(legacy))
  message("block pooled globals (k = ", fix$k, "):\n",
          .format_num_vec(pooled_num),
          "\naggregation = ", pooled$aggregation)
  message("legacy vs block pooled (matched quantities first):\n",
          .legacy_vs_pooled_message(legacy, pooled_num))
  message("legacy joint grid n_faces = ", legacy["n_faces"],
          "; per-block gs = ", paste(per_block$gs, collapse = ", "))

  expect_true(is.list(sim_env$Envelope))
  expect_true(is.list(sim_env$gamma_list))
  expect_true(is.list(sim_env$UB_list))
  expect_equal(sim_block$k, fix$k)

  ## Shared dispersion-interval bounds and data-side RSS floor (should align).
  expect_equal(unname(pooled_num["disp_lower"]), unname(legacy["disp_lower"]), tolerance = 1e-4)
  expect_equal(unname(pooled_num["disp_upper"]), unname(legacy["disp_upper"]), tolerance = 1e-4)
  expect_equal(unname(pooled_num["RSS_Min"]), unname(legacy["RSS_Min"]), tolerance = 1e-3)
  expect_equal(unname(pooled_num["rate2"]), unname(legacy["rate2"]), tolerance = 1e-3)

  expect_equal(as.integer(unname(pooled_num["n_identifiable"])), fix$k)

  ## Joint face-product globals should match legacy EnvelopeDispersionBuild constants.
  expect_equal(pooled$aggregation, "joint_face_product_edb")
  keys <- c(
    "lmc1", "lmc2", "lm_log1", "lm_log2",
    "max_New_LL_UB", "max_LL_log_disp", "shape3"
  )
  for (k in keys) {
    expect_equal(
      unname(pooled_num[k]),
      unname(legacy[k]),
      tolerance = 0.02,
      info = paste("mismatch on", k)
    )
  }

  ## Per-block EDB rows still use block-local geometry (not expected to match legacy joint).
  expect_equal(as.integer(per_block$gs), rep(3^fix$l1, fix$k))
  expect_equal(unname(legacy["n_faces"]), 3^(2L * fix$l1))

  ## With matched globals, block resample-until-accept iters should track legacy.
  set.seed(360)
  sim_legacy_ar <- glmbayesCore::rindepNormalGamma_reg(
    n = 500L,
    y = fix$y,
    x = fix$x_old,
    prior_list = fix$prior_list_old,
    n_envopt = n_envopt,
    Gridtype = 3L,
    use_parallel = FALSE,
    progbar = FALSE
  )
  set.seed(360)
  sim_block_ar <- lmebayesCore:::.rIndepNormalGammaRegBlock_cpp(
    n = 500L,
    y = fix$y,
    x = fix$x_block,
    block = fix$block,
    prior_list = fix$prior_list_block,
    n_envopt = n_envopt,
    Gridtype = 3L,
    use_parallel = FALSE,
    progbar = FALSE,
    verbose = FALSE,
    offset = rep(0, fix$n_obs),
    wt = rep(1, fix$n_obs),
    p_re = -1L,
    n_rss_iter = 10L,
    RSS_ML = NA_real_,
    use_opencl = FALSE,
    group_levels = character(0),
    re_names = character(0)
  )
  legacy_dpa <- mean(sim_legacy_ar$iters)
  block_dpa <- mean(sim_block_ar$iters_out)
  message(sprintf(
    "acceptance (n = 500): legacy mean(iters) %.3f vs block mean(iters) %.3f (ratio %.2f)",
    legacy_dpa, block_dpa, block_dpa / legacy_dpa
  ))
  expect_equal(
    sim_block_ar$sim$meta$accept_mode,
    "resample_until_accept_joint_product_slack_v2"
  )
  expect_gt(1 / block_dpa, 0.005)
  expect_equal(block_dpa, legacy_dpa, tolerance = 0.35)
})

# ---------------------------------------------------------------------------
# Three-way comparison: legacy vs joint-block vs independent-block
# For the stacked Dobson fixture (two identical blocks), Appendix A of
# BLOCK_ING_RINDEPNORMALGAMMA_REG.md shows that the lg-overbound slack is
# zero for every product face, so the independent sampler should give
# acceptance rates essentially identical to the joint sampler.
# ---------------------------------------------------------------------------
test_that("two-block Dobson: three-way means/SDs/acceptance — legacy, joint, ind", {
  fix <- .dobson_plant_two_block_fixture()
  n_envopt <- 10000L
  n_draws  <- 500L

  block_args <- list(
    y            = fix$y,
    x            = fix$x_block,
    block        = fix$block,
    prior_list   = fix$prior_list_block,
    prior_lists  = NULL,
    offset       = rep(0, fix$n_obs),
    wt           = rep(1, fix$n_obs),
    p_re         = -1L,
    n_rss_iter   = 10L,
    Gridtype     = 3L,
    n_envopt     = n_envopt,
    RSS_ML       = NA_real_,
    use_parallel = FALSE,
    use_opencl   = FALSE,
    progbar      = FALSE,
    verbose      = FALSE,
    group_levels = character(0),
    re_names     = character(0)
  )

  set.seed(2026)
  sim_legacy <- glmbayesCore::rindepNormalGamma_reg(
    n            = n_draws,
    y            = fix$y,
    x            = fix$x_old,
    prior_list   = fix$prior_list_old,
    n_envopt     = n_envopt,
    Gridtype     = 3L,
    use_parallel = FALSE,
    progbar      = FALSE
  )

  set.seed(2026)
  sim_joint <- do.call(
    lmebayesCore:::.rIndepNormalGammaRegBlock_cpp,
    c(list(n = n_draws), block_args)
  )

  set.seed(2026)
  sim_ind <- do.call(
    lmebayesCore:::.rIndepNormalGammaRegBlockInd_cpp,
    c(list(n = n_draws), block_args)
  )

  # --- accept_mode tags ---
  expect_equal(
    sim_joint$sim$meta$accept_mode,
    "resample_until_accept_joint_product_slack_v2"
  )
  expect_equal(
    sim_ind$sim$meta$accept_mode,
    "resample_until_accept_ind_v1"
  )
  expect_equal(
    sim_ind$sim$meta$face_draw_mode,
    "per_block_plsd_ind_v1"
  )

  # --- Coefficient matrices ---
  # legacy$coefficients is n x p  (each row = one draw)
  # block beta is l1 x n          (each row = one parameter, each col = one draw)
  legacy_coefs  <- sim_legacy$coefficients           # n x (2*l1)
  joint_coefs   <- rbind(                            # (2*l1) x n
    sim_joint$sim$block_results[[1]]$beta,
    sim_joint$sim$block_results[[2]]$beta
  )
  ind_coefs     <- rbind(
    sim_ind$sim$block_results[[1]]$beta,
    sim_ind$sim$block_results[[2]]$beta
  )

  legacy_means  <- colMeans(legacy_coefs)            # p-element vector
  joint_means   <- rowMeans(joint_coefs)             # p-element vector
  ind_means     <- rowMeans(ind_coefs)               # p-element vector
  legacy_sds    <- apply(legacy_coefs, 2L, sd)
  joint_sds     <- apply(joint_coefs,  1L, sd)
  ind_sds       <- apply(ind_coefs,    1L, sd)

  # legacy stores dispersion draws in $dispersion; block samplers use $disp_out
  legacy_disp_mean <- mean(sim_legacy$dispersion)
  joint_disp_mean  <- mean(sim_joint$disp_out)
  ind_disp_mean    <- mean(sim_ind$disp_out)
  legacy_disp_sd   <- sd(sim_legacy$dispersion)
  joint_disp_sd    <- sd(sim_joint$disp_out)
  ind_disp_sd      <- sd(sim_ind$disp_out)

  legacy_dpa <- mean(sim_legacy$iters)
  joint_dpa  <- mean(sim_joint$iters_out)
  ind_dpa    <- mean(sim_ind$iters_out)

  coef_labels <- paste0("b", seq_len(ncol(legacy_coefs)))
  message(sprintf(
    "Three-way comparison (n = %d):", n_draws
  ))
  message(sprintf(
    "  acceptance (mean iters): legacy %.3f | joint %.3f | ind %.3f",
    legacy_dpa, joint_dpa, ind_dpa
  ))
  for (i in seq_along(coef_labels)) {
    message(sprintf(
      "  %s mean: legacy %.4f | joint %.4f | ind %.4f",
      coef_labels[i], legacy_means[i], joint_means[i], ind_means[i]
    ))
    message(sprintf(
      "  %s sd:   legacy %.4f | joint %.4f | ind %.4f",
      coef_labels[i], legacy_sds[i], joint_sds[i], ind_sds[i]
    ))
  }
  message(sprintf(
    "  disp mean: legacy %.4f | joint %.4f | ind %.4f",
    legacy_disp_mean, joint_disp_mean, ind_disp_mean
  ))
  message(sprintf(
    "  disp sd:   legacy %.4f | joint %.4f | ind %.4f",
    legacy_disp_sd, joint_disp_sd, ind_disp_sd
  ))

  # --- Acceptance rate sanity ---
  # The Ind sampler is valid but less efficient: drawing block faces
  # independently can pair an "upp"-bound face from block 1 with a
  # "low"-bound face from block 2, producing non-zero slack even for
  # identical blocks (Appendix A §A.6 of BLOCK_ING_RINDEPNORMALGAMMA_REG.md).
  # For this two-block Dobson fixture the observed ratio is roughly 1.3–1.5x.
  expect_gt(1 / ind_dpa,   0.005)   # ind acceptance rate must be finite
  expect_gt(1 / joint_dpa, 0.005)  # joint acceptance rate must be finite
  expect_lt(ind_dpa / joint_dpa, 2.5)  # ind not more than 2.5x worse than joint

  # --- Coefficient means: each sampler draws from the same posterior,
  #     so E[b] should agree within Monte Carlo error (~3 * sd/sqrt(n)). ---
  expect_equal(ind_means, joint_means, tolerance = 0.20,
    info = "ind vs joint coefficient means"
  )
  expect_equal(ind_means, unname(legacy_means), tolerance = 0.20,
    info = "ind vs legacy coefficient means"
  )

  # --- Coefficient SDs: should agree within 30% relative. ---
  expect_equal(ind_sds, joint_sds, tolerance = 0.30,
    info = "ind vs joint coefficient SDs"
  )
  expect_equal(ind_sds, unname(legacy_sds), tolerance = 0.30,
    info = "ind vs legacy coefficient SDs"
  )

  # --- Dispersion mean and SD. ---
  expect_equal(ind_disp_mean, joint_disp_mean, tolerance = 0.20,
    info = "ind vs joint dispersion mean"
  )
  expect_equal(ind_disp_mean, legacy_disp_mean, tolerance = 0.20,
    info = "ind vs legacy dispersion mean"
  )
  expect_equal(ind_disp_sd, joint_disp_sd, tolerance = 0.30,
    info = "ind vs joint dispersion SD"
  )
})
