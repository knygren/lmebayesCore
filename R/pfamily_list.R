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
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class \code{"pfamily_list"} (also inherits from
#'   \code{"list"}): a named list whose elements are objects of class
#'   \code{"pfamily"}.  Attribute \code{"ptypes"} records the resolved
#'   prior-family name for each component.  Printing uses
#'   \code{\link{print.pfamily_list}}.
#'
#' @seealso \code{\link{Prior_Setup_GLMM}}, \code{\link{print.pfamily_list}},
#'   \code{\link[glmbayesCore]{pfamily}}, \code{\link[glmbayesCore]{dNormal}},
#'   \code{\link[glmbayesCore]{dIndependent_Normal_Gamma}},
#'   \code{\link{dGamma_list}}
#' @example inst/examples/Ex_pfamily_list.R
#' @order 1
#' @export
pfamily_list <- function(object, ptypes = "dNormal", ...) {
  UseMethod("pfamily_list")
}

#' Subset a pfamily_list (keeps class and \code{ptypes})
#'
#' @param x An object of class \code{"pfamily_list"}.
#' @param i Character component names (e.g. \code{"(Intercept)"}, slope
#'   names), integer indices, or a logical vector.
#' @param j Ignored (list subsetting).
#' @param ... Ignored.
#' @param drop Ignored (always returns a list).
#' @return A \code{"pfamily_list"} with the selected population / group-effect
#'   components.
#' @export
#' @method [ pfamily_list
`[.pfamily_list` <- function(x, i, j, ..., drop = TRUE) {
  if (!missing(j)) {
    stop("pfamily_list objects are one-dimensional; do not use [i, j].",
         call. = FALSE)
  }
  all_names <- names(x)
  if (is.null(all_names)) {
    all_names <- as.character(seq_along(x))
    names(x) <- all_names
  }
  if (missing(i)) {
    return(x)
  }
  if (is.logical(i)) {
    if (length(i) != length(x)) {
      stop("Logical 'i' must have length ", length(x), ".", call. = FALSE)
    }
    keep <- all_names[i]
  } else {
    keep <- .lmebayes_select_named_list_keys(
      x, i, arg = "i", what = "component"
    )
  }

  out <- unclass(x)[keep]
  ptypes <- attr(x, "ptypes")
  if (!is.null(ptypes)) {
    attr(out, "ptypes") <- ptypes[keep]
  }
  class(out) <- c("pfamily_list", "list")
  out
}

#' Print a pfamily_list
#'
#' Prints each selected population / group-effect component with
#' \code{\link[glmbayesCore]{print.pfamily}}, one after another.
#' Component names are \code{names(x)} (e.g. \code{"(Intercept)"},
#' random-slope names).
#'
#' @param x An object of class \code{"pfamily_list"}.
#' @param components \code{NULL} (default; print all components), a
#'   character vector of component names (\code{names(x)}), or integer
#'   indices into \code{x}.
#' @param ... Passed to \code{print.pfamily}.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{pfamily_list}}, \code{\link[glmbayesCore]{pfamily}},
#'   \code{\link{print.dGamma_list}}
#' @export
#' @method print pfamily_list
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
