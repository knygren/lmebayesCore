#' Build a named list of pfamily objects
#'
#' Generic for constructing a named list of \code{\link[glmbayesCore]{pfamily}}
#' prior objects from a prior-specification object.  Methods exist for
#' \code{\link{Prior_Setup_GLMM}} (one pfamily per population / group-effect
#' coefficient) and \code{\link{Prior_SetupGroup}} (one pfamily per row
#' block).
#'
#' @param object A prior-specification object (typically from
#'   \code{\link{Prior_Setup_GLMM}} or \code{\link{Prior_SetupGroup}}).
#' @param ptypes Character: either a single string applied to every
#'   component, or a character vector / list with one string per
#'   component.  See the method-specific Details for allowed values and
#'   naming (RE coefficient names for \code{Prior_Setup_GLMM}; block IDs
#'   for \code{Prior_SetupGroup}).  For both methods,
#'   \code{ptypes = NULL} (default) resolves to \code{"dNormal"}.
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
#' @seealso \code{\link{Prior_Setup_GLMM}}, \code{\link{Prior_SetupGroup}},
#'   \code{\link[glmbayesCore]{pfamily}}, \code{\link[glmbayesCore]{dNormal}},
#'   \code{\link[glmbayesCore]{dNormal_Gamma}},
#'   \code{\link[glmbayesCore]{dIndependent_Normal_Gamma}},
#'   \code{\link{dGamma_list}}
#' @example inst/examples/Ex_pfamily_list.R
#' @order 1
#' @export
pfamily_list <- function(object, ptypes = NULL, ...) {
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
#' @order 4
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
