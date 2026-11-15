# Legacy / optional code (not part of the installed package)

Files here are excluded from `R CMD build` via `.Rbuildignore`.
Copy or source manually when experimenting; they are not exported from **glmbayes**.

## `rIndepNormalGammaReg_with_envelope/`

lmebayesCore's former copy of `rIndepNormalGammaReg()` (C++, with an extra
`return_envelope` flag) and the internal `rindepNormalGamma_reg_with_envelope()`
diagnostic sampler built on it, with its `.Call` wrapper and Rcpp export.
Moved out of `R/` and `src/` so lmebayesCore uses glmbayesCore's
`rIndepNormalGammaReg`. The first block of
`tests/testthat/test-rIndepNormalGammaRegBlock-iters-parity-two-block.R` is
skipped while this lives here.

## `pp_check.glmb.R`

Former S3 method `pp_check.glmb()` (wrapper around **bayesplot** for `glmb` fits).
Removed from `R/` so **glmbayes** does not depend on **bayesplot** (~41 packages on R-universe). **bayesplot** is not in **Imports** or **Suggests**.

To use locally:

1. `install.packages("bayesplot")`
2. `source("legacy_code/pp_check.glmb.R")` from the package root (after `library(glmbayes)`)
3. Call `bayesplot::pp_check(fit, ...)` or register the method by sourcing the file

Vignette chunks that previously called this method are commented out; see
`vignettes/Chapter-05.Rmd`, `Chapter-06.Rmd`, and `Chapter-12.Rmd` (in **glmbayes**).
