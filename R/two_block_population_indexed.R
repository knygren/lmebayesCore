## Population block of the indexed two-block sweep. Each step has a C++ twin
## in src/two_block_population_indexed.cpp
## (glmbayes::sim::twoblock::population_*) with the same inputs, output
## shape, and loop order. `spec` is from .two_block_spec().

#' Population response for component \code{k}, in the row order of \code{W_k}
#' @return Numeric vector of length \code{nrow(W_k)}.
#' @noRd
.two_block_population_response <- function(spec, b_i, k) {
  b_i[spec$b_for_wrow[[k]], k]
}

#' Population draw for component \code{k}
#' @return List with \code{coef} and \code{coef_mode}, numeric vectors in
#'   \code{W_k} column order.
#' @noRd
.two_block_population_draw <- function(spec, y_k, k) {
  n_w <- length(y_k)
  out <- .rNormalReg_cpp(
    n          = 1L,
    y          = y_k,
    x          = spec$W[[k]],
    mu         = spec$mu2[[k]],
    P          = spec$P2[[k]],
    offset     = numeric(n_w),
    wt         = rep(1, n_w),
    dispersion = spec$disp2[k],
    f2         = spec$f2_gauss,
    f3         = spec$f3_gauss,
    start      = spec$mu2[[k]],
    family     = "gaussian",
    link       = "identity",
    Gridtype   = spec$gridtype
  )
  list(
    coef      = as.numeric(out$coefficients[1L, ]),
    coef_mode = as.numeric(out$coef.mode)
  )
}

#' Write a population draw into row \code{i} of the stored fixef matrix
#' @return \code{fixef_k} with \code{fixef_k[i, w_col_k] <- coef}.
#' @noRd
.two_block_population_write <- function(spec, fixef_k, i, k, coef) {
  fixef_k[i, spec$w_col[[k]]] <- coef
  fixef_k
}
