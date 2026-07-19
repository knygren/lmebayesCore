#' Build a named list of dGamma measurement-dispersion priors
#'
#' Generic for constructing a named list of \code{\link[glmbayesCore]{dGamma}} prior
#' objects from a prior-specification object.  For mixed models, one
#' \code{dGamma()} per group level supplies Block~1 observation-level
#' \eqn{\sigma^2_j} priors for \code{lmerb(..., dispersion_ranef = ...)}.
#'
#' @param object A prior-specification object.
#' @param ... Additional arguments passed to methods (e.g.
#'   \code{max_disp_perc}, \code{disp_upper_anchor}).
#' @return A named list whose elements are objects of class \code{"pfamily"}.
#' @seealso \code{\link[glmbayesCore]{dGamma}}, \code{\link{Prior_Setup_lmebayes}},
#'   \code{\link{pfamily_list}}
#' @export
dGamma_list <- function(object, ...) UseMethod("dGamma_list")
