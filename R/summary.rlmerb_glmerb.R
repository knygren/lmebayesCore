## summary() for rglmerb. The rlmerb method lives in glmbayesCore.

#' Summary for an rglmerb fit
#'
#' Uses \code{summary.rlmerb()} from \pkg{glmbayesCore}, which accepts an
#' \code{rglmerb} fit as well as an \code{rlmerb} fit.
#'
#' @param object An \code{rglmerb} fit.
#' @param groups Optional groups to include in the group-effect detail.
#'   Passed through to \code{summary.rlmerb()}.
#' @param ... Passed to \code{summary.rlmerb()}.
#' @return An object of class \code{"summary.rlmerb"}, as returned by
#'   \code{summary.rlmerb()} in \pkg{glmbayesCore}.
#' @seealso \code{\link[glmbayesCore]{summary.rlmerb}}
#' @export
#' @method summary rglmerb
summary.rglmerb <- function(object, groups = NULL, ...) {
  glmbayesCore:::summary.rlmerb(object, groups = groups, ...)
}
