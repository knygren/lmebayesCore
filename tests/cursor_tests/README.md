# cursor_tests

Tests moved here were added during Cursor-assisted development and are **not** part of the official `R CMD check` / `testthat` suite.

The package runs only `tests/testthat/test-prior-setup-poisson-conj.R` via `test_check("glmbayesCore")`.

To run these files locally (optional):

```r
testthat::test_dir("tests/cursor_tests", package = "glmbayesCore")
```

Do not add new files here without explicit approval; prefer extending `tests/testthat/` after review.
