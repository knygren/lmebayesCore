#' Model Formulae for \code{summary.rglmb} Objects
#'
#' Extract a formula for a \code{summary.rglmb} object by refitting a reference
#' \code{\link[stats]{glm}} with the stored response and design matrix.
#'
#' @param x an object of class \code{summary.rglmb}, typically from
#'   \code{\link{summary.rglmb}}.
#' @param ... further arguments passed to or from other methods.
#' @return A model formula.
#' @seealso \code{\link{rglmb}}, \code{\link{summary.rglmb}}, \code{\link{rlmb}},
#'   \code{\link[stats]{formula}}.
#' @export
#' @method formula summary.rglmb
formula.summary.rglmb <- function(x, ...) {
  stats::formula(stats::glm(x$y ~ x$x - 1, family = stats::family(x)))
}
