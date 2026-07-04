#' Convergence rate from a \code{pfamily_list} Block~2 spec
#'
#' Thin wrapper around \code{\link{two_block_rate}} that accepts the Block~2
#' priors as \code{pfamily} objects.  For \code{dNormal} components the fixed
#' \code{dispersion} is used; for \code{dIndependent_Normal_Gamma} components
#' the conservative \code{disp_lower} plug-in is used (the rate is then an
#' upper bound over the truncated tau^2 range).
#'
#' @param x,block,x_hyper,prior_list_block1,weights,family,group_levels As in
#'   \code{\link{two_block_rate}}.
#' @param pfamily_list Named list of \code{pfamily} objects (one per RE
#'   component), as in \code{\link{two_block_rNormal_reg}}.
#' @return Object of class \code{"two_block_rate"}.
#' @family simfuncs
#' @seealso \code{\link{two_block_rate}},
#'   \code{\link{two_block_rNormal_reg}}
#' @export
two_block_rate_from_pfamily_list <- function(x,
                                              block,
                                              x_hyper,
                                              prior_list_block1,
                                              pfamily_list,
                                              weights = NULL,
                                              family = gaussian(),
                                              group_levels = levels(block)) {
  re_names <- names(x_hyper)
  pfamily_list <- .two_block_validate_pfamily_list(pfamily_list, re_names)
  prior_list_block2 <- lapply(pfamily_list, function(pf) {
    pl <- pf$prior_list
    list(
      mu = pl$mu,
      Sigma = pl$Sigma,
      dispersion = if (identical(pf$pfamily, "dNormal")) {
        pl$dispersion
      } else {
        pl$disp_lower
      }
    )
  })
  two_block_rate(
    x = x, block = block, x_hyper = x_hyper,
    prior_list_block1 = prior_list_block1,
    prior_list_block2 = prior_list_block2,
    weights = weights, family = family, group_levels = group_levels
  )
}
