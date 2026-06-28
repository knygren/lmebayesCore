## Parity: C++ two_block_block1_one_chain vs v6 R prep + draw (piecewise C++ TRUE).

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists("two_block_block1_one_chain_cpp_export", mode = "function"))

compare_prep <- function(prep_r, cpp_out) {
  stopifnot(all.equal(prep_r$mu_all, cpp_out$mu_all, tolerance = 0))
  stopifnot(all.equal(prep_r$prior_list$P, cpp_out$prior_list$P, tolerance = 0))
  stopifnot(identical(prep_r$prior_list$ddef, cpp_out$prior_list$ddef))
  invisible(TRUE)
}

one_chain_cpp_args <- function(i) {
  c(
    list(
      chain_i      = i,
      batch_fixef  = batch$fixef,
      tau2_i       = batch$tau2[i, ],
      y            = as.numeric(design$y),
      Z            = as.matrix(design$Z),
      groups       = design$groups,
      offset       = offset,
      wt           = wt,
      x_hyper      = x_hyper,
      re_names     = batch$re_names,
      group_levels = batch$group_levels,
      ptypes       = ptypes,
      block1_prior = block1_prior,
      is_gaussian  = FALSE,
      f2           = fam$f2,
      f3           = fam$f3,
      f2_gauss     = fam_g$f2,
      f3_gauss     = fam_g$f3,
      family       = "binomial",
      link         = "logit",
      Gridtype     = 2L,
      n_envopt     = 1L
    )
  )
}

J <- 4L
n_chains <- 3L
group_levels <- c("d", "a", "c", "b")
re_names <- c("(Intercept)", "x")
ptypes <- c("dNormal", "dIndependent_Normal_Gamma")
X_int <- matrix(1, nrow = J, ncol = 1)
X_slope <- cbind(x = c(0.1, 0.2, 0.3, 0.4), z = c(1, 0, 1, 0))
rownames(X_int) <- sort(unique(group_levels))
rownames(X_slope) <- sort(unique(group_levels))
design_Z <- matrix(rnorm(J * 2L * 2L), ncol = 2L)
colnames(design_Z) <- re_names
design <- list(
  y = rbinom(J * 2L, 1, 0.3),
  Z = design_Z,
  groups = factor(rep(group_levels, each = 2L), levels = sort(unique(group_levels))),
  X_hyper = list("(Intercept)" = X_int, "x" = X_slope),
  re_coef_names = re_names
)
block1_prior <- list(P = diag(c(10, 20)), ddef = TRUE)
batch <- glmbayesCore:::.two_block_batch_init(
  n_chains     = n_chains,
  start_fixef  = list(
    "(Intercept)" = c("(Intercept)" = 0.5),
    "x" = c(x = 0.1, z = -0.2)
  ),
  b_start      = matrix(0, nrow = J, ncol = length(re_names),
                        dimnames = list(group_levels, re_names)),
  tau2_start   = c("(Intercept)" = 1, "x" = 2),
  re_names     = re_names,
  group_levels = group_levels
)

l2 <- length(design$y)
offset <- rep(0, l2)
wt <- rep(1, l2)
x_hyper <- lapply(design$X_hyper, as.matrix)
fam <- glmbayesCore::glmbfamfunc(stats::binomial())
fam_g <- glmbayesCore::glmbfamfunc(stats::gaussian())
family <- stats::binomial()

for (i in seq_len(n_chains)) {
  prep_r <- glmbayesCore:::.two_block_block1_prep_one_chain(
    batch, i, design, block1_prior, ptypes,
    use_cpp_mu_all = TRUE,
    use_cpp_prior_tau2 = TRUE
  )
  cpp_out <- do.call(
    two_block_block1_one_chain_cpp_export,
    one_chain_cpp_args(i)
  )
  compare_prep(prep_r, cpp_out)
  stopifnot(identical(dim(prep_r$mu_all), dim(cpp_out$mu_all)))
  stopifnot(identical(dim(cpp_out$b), c(J, length(re_names))))
  stopifnot(identical(rownames(cpp_out$b), group_levels))
  stopifnot(identical(colnames(cpp_out$b), re_names))
  stopifnot(is.finite(cpp_out$iters_mean))

  set.seed(1000L + i)
  prep_draw <- glmbayesCore:::.two_block_block1_prep_one_chain(
    batch, i, design, block1_prior, ptypes,
    use_cpp_mu_all = TRUE,
    use_cpp_prior_tau2 = TRUE
  )
  draw_r <- glmbayesCore:::.two_block_block1_draw_one_chain(
    prior_list      = prep_draw$prior_list,
    design          = design,
    family          = family,
    is_gaussian     = FALSE,
    group_levels    = batch$group_levels,
    use_cpp_reorder = TRUE,
    use_cpp_iters   = TRUE
  )

  set.seed(1000L + i)
  cpp_draw <- do.call(
    two_block_block1_one_chain_cpp_export,
    one_chain_cpp_args(i)
  )

  stopifnot(all.equal(draw_r$b, cpp_draw$b, tolerance = 0))
  stopifnot(all.equal(draw_r$iters_mean, cpp_draw$iters_mean, tolerance = 0))
  stopifnot(identical(rownames(draw_r$b), rownames(cpp_draw$b)))
  stopifnot(identical(colnames(draw_r$b), colnames(cpp_draw$b)))
}

message("test_block1_one_chain_cpp.R: OK")
