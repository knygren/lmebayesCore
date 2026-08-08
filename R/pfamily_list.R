#' Build a named list of pfamily objects
#'
#' Generic for constructing a named list of \code{\link[glmbayesCore]{pfamily}}
#' prior objects from a prior-specification object.  For
#' \code{\link{Prior_Setup_GLMM}} objects, each population / group-effect
#' coefficient (e.g. \code{"(Intercept)"}, slope names) is mapped to a
#' \code{\link[glmbayesCore]{dNormal}} or
#' \code{\link[glmbayesCore]{dIndependent_Normal_Gamma}} prior according to
#' \code{ptypes}.
#'
#' @param object A prior-specification object (typically from
#'   \code{\link{Prior_Setup_GLMM}}).
#' @param ptypes Character: either a single string applied to every
#'   group-effect coefficient, or a character vector / list with one
#'   string per coefficient.  For \code{\link{Prior_Setup_GLMM}}, allowed
#'   values are \code{"dNormal"} (default; known \eqn{\tau^2_k}) and
#'   \code{"dIndependent_Normal_Gamma"} (Gamma prior on precision
#'   \eqn{1/\tau^2_k}).  A vector may be named with the group-effect
#'   coefficient names (any order); unnamed vectors are matched
#'   positionally against \code{names(object$pop.prior_list)}.
#' @param ... Additional arguments passed to methods. For
#'   \code{print.pfamily_list()}, passed to \code{print.pfamily}.
#'
#' @return For \code{pfamily_list()}, an object of class
#'   \code{"pfamily_list"} (also inherits from \code{"list"}): a named
#'   list whose elements are objects of class \code{"pfamily"}.
#'   Attribute \code{"ptypes"} records the resolved prior-family name
#'   for each component.  For \code{print.pfamily_list()}, \code{x}
#'   invisibly.
#'
#' @seealso \code{\link{Prior_Setup_GLMM}},
#'   \code{\link[glmbayesCore]{pfamily}}, \code{\link[glmbayesCore]{dNormal}},
#'   \code{\link[glmbayesCore]{dIndependent_Normal_Gamma}},
#'   \code{\link{dGamma_list}}
#' @example inst/examples/Ex_pfamily_list.R
#' @order 1
#' @export
pfamily_list <- function(object, ptypes = "dNormal", ...) {
  UseMethod("pfamily_list")
}

#' @rdname pfamily_list
#' @method print pfamily_list
#' @param x An object of class \code{"pfamily_list"}.
#' @param components For \code{print.pfamily_list()}: \code{NULL}
#'   (default; print all components), a character vector of component
#'   names (\code{names(x)}, e.g. \code{"(Intercept)"}, slopes), or
#'   integer indices into \code{x}. Each selected component is printed
#'   with \code{\link[glmbayesCore]{print.pfamily}}.
#' @order 3
#' @export
print.pfamily_list <- function(x, components = NULL, ...) {
  sel <- .lmebayes_select_named_list_keys(
    x, components, arg = "components", what = "component"
  )
  for (nm in sel) {
    cat(sprintf("\n[[%s]]\n", nm))
    print(x[[nm]], ...)
  }
  invisible(x)
}
