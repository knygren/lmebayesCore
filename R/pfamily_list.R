#' Build a named list of pfamily objects
#'
#' Generic for constructing a named list of \code{\link[glmbayesCore]{pfamily}} prior
#' objects from a prior-specification object.  Packages that define prior
#' setup containers (e.g. \code{\link{Prior_Setup_lmebayes}}) provide methods
#' that map each component of the container to a \code{pfamily} constructor such as \code{\link[glmbayesCore]{dNormal}} or
#' \code{\link[glmbayesCore]{dIndependent_Normal_Gamma}}.
#'
#' @param object A prior-specification object.
#' @param ... Additional arguments passed to methods (e.g. \code{ptypes}).
#' @return A named list whose elements are objects of class
#'   \code{"pfamily"}.
#' @seealso \code{\link[glmbayesCore]{pfamily}}, \code{\link[glmbayesCore]{dNormal}},
#'   \code{\link[glmbayesCore]{dIndependent_Normal_Gamma}}, \code{\link{Prior_Setup_lmebayes}}
#' @export
pfamily_list <- function(object, ...) UseMethod("pfamily_list")
