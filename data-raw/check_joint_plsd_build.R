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
block <- factor(c(rep("B1", n1), rep("B2", n1)), levels = c("B1", "B2"))
prior_list_block <- list(
  mu = p_setup$mu,
  Sigma = p_setup$Sigma,
  shape = p_setup$shape,
  rate = p_setup$rate,
  max_disp_perc = 0.99
)
sim_block <- glmbayesCore:::.rIndepNormalGammaRegBlock_cpp(
  n = 1L, y = y, x = x_block, block = block,
  prior_list = prior_list_block, n_envopt = n_envopt, Gridtype = 3L,
  use_parallel = FALSE, progbar = FALSE, verbose = FALSE,
  offset = rep(0, length(y)), wt = rep(1, length(y)),
  p_re = -1L, n_rss_iter = 10L, RSS_ML = NA_real_, use_opencl = FALSE,
  group_levels = character(0), re_names = character(0)
)
de <- sim_block$build_out$dispersion_envelope
jp <- de$joint_PLSD
cat(sprintf("joint_PLSD: n=%d sum=%.10f max=%.6f\n", length(jp), sum(jp), max(jp)))
cat(sprintf(
  "face_draw_mode (build) = %s\nsim meta = %s\n",
  de$cross_face_meta$face_draw_mode,
  sim_block$sim$meta$face_draw_mode
))
