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
#' @details
#' A \code{pfamily_list} is the prior argument every mixed-model sampler
#' takes. It exists because the model requires a \emph{diagonal}
#' random-effect covariance \eqn{\Psi}: with \eqn{\Psi} diagonal, the
#' population (Block~2) update factorizes across random-effect
#' coefficients, so each coefficient carries its own prior object rather
#' than sharing one joint prior over all of \eqn{\gamma}.
#'
#' @section What the choice of family controls:
#' The family assigned to component \eqn{k} decides whether that
#' component's variance component \eqn{\tau^2_k} is a constant or a
#' parameter, and that single decision propagates through the whole engine:
#' \tabular{lll}{
#'   \strong{Family} \tab \strong{\eqn{\tau^2_k}} \tab \strong{Consequence} \cr
#'   \code{"dNormal"} \tab treated as known \tab Posterior stays Gaussian; draws are exact \cr
#'   \code{"dIndependent_Normal_Gamma"} \tab estimated, Gamma prior on \eqn{1/\tau^2_k} \tab Posterior is non-Gaussian; draws come from the final sweep of a calibrated chain \cr
#' }
#' Declaring \eqn{\tau^2_k} known asserts that you already know how widely
#' groups scatter on coefficient \eqn{k}; estimating it lets the data
#' revise that. Both are calibrated to the same center --- by construction
#' \eqn{\mathrm{rate}_k/(\mathrm{shape}_k - 1)} equals the \code{dNormal}
#' dispersion --- so switching \code{ptypes} changes whether a component can
#' move without changing where its prior sits.
#'
#' A single \code{dIndependent_Normal_Gamma} component is enough to make the
#' whole posterior non-Gaussian; there is no partial-exactness discount.
#' Mixing families is nonetheless coherent, and is specified by passing a
#' \emph{named} \code{ptypes} vector keyed by random-effect coefficient
#' name.
#'
#' @section Where the precision comes from:
#' \eqn{\Psi^{-1}} -- the group-level prior precision the samplers call
#' \code{P} -- is never supplied by the caller. Every routed export derives
#' it from the validated \code{pfamily_list} via one \eqn{\tau^2_k} plug-in
#' per component (the fixed \code{dispersion} for \code{dNormal}, the prior
#' mean \eqn{\mathrm{rate}/(\mathrm{shape}-1)} for
#' \code{dIndependent_Normal_Gamma}), and rejects a caller-supplied
#' \code{P} or \code{Sigma}. The prior list is therefore the single source
#' of truth for the population covariance.
#'
#' @section Overriding a component by hand:
#' The result is a plain named list, so an element can be replaced with any
#' \code{pfamily} object of matching dimension. Three quantities are easy to
#' confuse: \code{mu} is the prior mean of the \emph{population}
#' coefficient \eqn{\gamma_k}, \code{Sigma} is the uncertainty about that
#' population coefficient, and \code{dispersion} is \eqn{\tau^2_k}, the
#' scatter of \emph{individual group} coefficients around
#' \eqn{\mathcal{W}_j\gamma}. They are independent statements.
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
