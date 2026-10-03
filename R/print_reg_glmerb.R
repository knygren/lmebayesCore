## Print methods for rGLMM_reg / rglmerb (glmerb only).
## Must collate after print_reg.R (aliases print_groupef.default).

#' @rdname rGLMM_reg
#' @order 2
#' @method print rGLMM_reg
#' @param x For \code{print} methods, a simfunc return object.
#' @param digits Number of significant digits for printed numeric values.
#' @param components For \code{print} methods: \code{NULL} (default; print all
#'   population components in \code{names(x$popef)}), a character vector of
#'   component names (e.g. \code{"(Intercept)"}, slopes), or integer indices.
#' @param draws For \code{print} methods and \code{print_groupef}: \code{NULL}
#'   (default; all draws) or integer indices / logical of length \code{n}, as in
#'   \code{rnorm(n)[1:10]}. Only the printed rows are subset; the object is
#'   unchanged. Example: \code{print(x, draws = 1:10)}.
#' @param ... further arguments passed to or from other methods.
#' @export
print.rGLMM_reg <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    components = NULL,
    draws = NULL,
    ...
) {
  .lmebayes_print_reg_popef(
    x, digits = digits, components = components, draws = draws, ...
  )
  invisible(x)
}

#' @rdname print_groupef
#' @export
print_groupef.rGLMM_reg <- print_groupef.default

#' @rdname print_groupef
#' @export
print_groupef.rglmerb <- print_groupef.default
