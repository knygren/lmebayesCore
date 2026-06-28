## Parity: C++ two_block_block1_prior_with_tau2 vs R reference.

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists(".two_block_block1_prior_with_tau2_r", mode = "function",
                 where = asNamespace("glmbayesCore")))
stopifnot(exists("two_block_block1_prior_with_tau2_cpp_export", mode = "function"))

compare_prior <- function(base_prior, tau2_vec, ptypes, re_names, mu_all) {
  r_out <- glmbayesCore:::.two_block_block1_prior_with_tau2_r(
    base_prior, tau2_vec, ptypes, re_names, mu_all
  )
  cpp_out <- glmbayesCore:::.two_block_block1_prior_with_tau2(
    base_prior, tau2_vec, ptypes, re_names, mu_all, use_cpp = TRUE
  )
  stopifnot(identical(r_out$ddef, cpp_out$ddef))
  stopifnot(identical(r_out$dispersion, cpp_out$dispersion))
  stopifnot(all.equal(r_out$mu, cpp_out$mu, tolerance = 0))
  stopifnot(all.equal(r_out$P, cpp_out$P, tolerance = 0))
  invisible(list(r = r_out, cpp = cpp_out))
}

re_names <- c("(Intercept)", "violent_i")
ptypes_all_ing <- c("dIndependent_Normal_Gamma", "dIndependent_Normal_Gamma")
ptypes_mixed <- c("dNormal", "dIndependent_Normal_Gamma")
mu_all <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 2, byrow = TRUE)
rownames(mu_all) <- re_names
colnames(mu_all) <- c("a", "b", "c")

base_P <- diag(c(10, 20))
tau2 <- c(0.25, 4)

## glmerb / lmebayes shape: list(P, ddef = TRUE) with no "dispersion" name at all.
base_prior_glmerb <- list(P = base_P, ddef = TRUE)
compare_prior(base_prior_glmerb, tau2, ptypes_all_ing, re_names, mu_all)

base_prior <- list(
  P          = base_P,
  dispersion = NULL,
  ddef       = TRUE
)

## Key case: ddef = TRUE (Poisson/binomial glmerb) must still refresh ING rows.
compare_prior(base_prior, tau2, ptypes_all_ing, re_names, mu_all)

base_prior$ddef <- FALSE
compare_prior(base_prior, tau2, ptypes_all_ing, re_names, mu_all)
compare_prior(base_prior, tau2, ptypes_mixed, re_names, mu_all)

ptypes_dnormal <- c("dNormal", "dNormal")
compare_prior(base_prior, tau2, ptypes_dnormal, re_names, mu_all)

## Direct export vs R reference
r_out <- glmbayesCore:::.two_block_block1_prior_with_tau2_r(
  base_prior, tau2, ptypes_all_ing, re_names, mu_all
)
cpp_out <- two_block_block1_prior_with_tau2_cpp_export(
  base_prior = base_prior,
  tau2_vec   = tau2,
  ptypes     = ptypes_all_ing,
  re_names   = re_names,
  mu_all     = mu_all
)
stopifnot(all.equal(r_out$P, cpp_out$P, tolerance = 0))

message("test_block1_prior_with_tau2_cpp.R: OK")
