## Groups block of the indexed two-block sweep. Each step has a C++ twin in
## src/two_block_groups_indexed.cpp (glmbayes::sim::twoblock::groups_*) with
## the same inputs, output shape, and loop order. `spec` is from
## .two_block_spec(); index vectors in it are 1-based.

#' Groups prior mean for chain \code{i}
#'
#' Row \code{d_col[k]}, column \code{j} is
#' \code{sum_c W_k[w_for_group_k[j], c] * fixef_k[i, w_col_k[c]]}, summed over
#' \code{c} in increasing order.
#' @return Numeric matrix \code{ncol(D) x J}.
#' @noRd
.two_block_groups_prior_mean <- function(spec, fixef, i) {
  l1 <- ncol(spec$D)
  J <- length(spec$obs_rows)
  mu_D <- matrix(0, l1, J)
  for (k in seq_along(spec$d_col)) {
    Wk <- spec$W[[k]]
    wg <- spec$w_for_group[[k]]
    wc <- spec$w_col[[k]]
    s <- numeric(J)
    for (c in seq_len(ncol(Wk))) {
      s <- s + Wk[wg, c] * fixef[[k]][i, wc[c]]
    }
    mu_D[spec$d_col[k], ] <- s
  }
  mu_D
}

#' Groups draw for one chain
#' @return List with \code{coef} and \code{coef_mode} (both \code{J x ncol(D)})
#'   and \code{iters} (mean envelope candidates per group).
#' @noRd
.two_block_groups_draw <- function(spec, mu_D) {
  out <- .rNormalGLMGroups_cpp(
    n              = 1L,
    y              = spec$y,
    x              = spec$D,
    offset         = spec$offset1,
    wt             = spec$wt1,
    dispersion     = spec$disp1,
    mu             = mu_D,
    P_blocks       = list(spec$P_D),
    prior_by_block = TRUE,
    row_blocks     = spec$obs_rows,
    f2             = spec$f2,
    f3             = spec$f3,
    family         = spec$family,
    link           = spec$link,
    Gridtype       = spec$gridtype,
    n_envopt       = spec$n_envopt,
    use_parallel   = FALSE,
    use_opencl     = FALSE,
    verbose        = FALSE
  )
  it <- vapply(out$group_results, function(g) as.numeric(g$iters)[1L], numeric(1))
  list(
    coef      = out$coefficients,
    coef_mode = out$coef.mode,
    iters     = sum(it) / length(it)
  )
}

#' Columns of a groups draw in stored \code{b} order
#' @return Numeric matrix \code{J x p_re}.
#' @noRd
.two_block_groups_place <- function(spec, coef) {
  out <- coef[, spec$d_col, drop = FALSE]
  dimnames(out) <- NULL
  out
}
