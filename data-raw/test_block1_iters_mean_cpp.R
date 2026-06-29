## Parity: C++ two_block_block1_iters_mean vs R reference (.two_block_block1_iters_mean_r).
## Does not use v5 drivers — compares implementations on the same block_out object.

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists(".two_block_block1_iters_mean_r", mode = "function", where = asNamespace("glmbayesCore")))
stopifnot(exists(".two_block_block1_iters_mean_cpp", mode = "function",
                 where = asNamespace("glmbayesCore")))

## Synthetic block_out shapes (no sampler required)
empty_br <- list(block_results = list())
stopifnot(identical(
  glmbayesCore:::.two_block_block1_iters_mean_r(empty_br),
  as.numeric(glmbayesCore:::.two_block_block1_iters_mean_cpp(empty_br))
))

br_one <- list(block_results = list(list(iters = matrix(3, 1, 1))))
stopifnot(identical(
  glmbayesCore:::.two_block_block1_iters_mean_r(br_one),
  as.numeric(glmbayesCore:::.two_block_block1_iters_mean_cpp(br_one))
))

br_mix <- list(block_results = list(
  list(iters = matrix(2, 1, 1)),
  list(iters = NULL),
  list(iters = 4)
))
r_val <- glmbayesCore:::.two_block_block1_iters_mean_r(br_mix)
cpp_val <- as.numeric(glmbayesCore:::.two_block_block1_iters_mean_cpp(br_mix))
stopifnot(identical(r_val, cpp_val))

## Real Block 1 draw (binomial GLMM block path) if lmebayes available
if (requireNamespace("lmebayes", quietly = TRUE)) {
  data(book_banning, package = "lmebayes")
  ps <- lmebayes::Prior_Setup_lmebayes(
    formula = cbind(banned, unbanned) ~ 1 + logpop,
    data = book_banning,
    family = binomial(),
    re_formula = ~ 1,
    group = county
  )
  design <- ps$design
  prior <- ps$prior
  block1_prior <- lmebayes:::.lmebayes_block1_prior_list(prior, dispersion_ranef = NULL)
  fixef <- lapply(ps$fixef_init, function(x) x[1, , drop = TRUE])
  mu_all <- as.matrix(glmbayesCore::build_mu_all(design, fixef)$mu_all)
  prior_list <- glmbayesCore:::.two_block_block1_prior_with_tau2(
    block1_prior, rep(1, length(ps$re_names)), ps$ptypes, ps$re_names, mu_all
  )
  block_out <- glmbayesCore::block_rNormalGLM(
    n = 1L,
    y = design$y,
    x = design$Z,
    block = design$groups,
    prior_list = prior_list,
    family = binomial(),
    use_parallel = FALSE,
    verbose = FALSE,
    progbar = FALSE
  )
  r_live <- glmbayesCore:::.two_block_block1_iters_mean_r(block_out)
  cpp_live <- as.numeric(glmbayesCore:::.two_block_block1_iters_mean_cpp(block_out))
  if (!identical(r_live, cpp_live)) {
    stop(
      "Block 1 iters_mean mismatch on live draw: R=", r_live, " C++=", cpp_live,
      call. = FALSE
    )
  }
  message("Live block_rNormalGLM draw: R and C++ iters_mean = ", r_live)
}

message("test_block1_iters_mean_cpp.R: OK")
