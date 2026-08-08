#' Build a named list of dGamma measurement-dispersion priors
#'
#' Generic for constructing a named list of \code{\link[glmbayesCore]{dGamma}} prior
#' objects from a prior-specification object.  For mixed models, one
#' \code{dGamma()} per group level supplies observation-level
#' \eqn{\sigma^2_j} priors for \code{lmerb(..., group.dispersion = ...)}.
#'
#' @param object A prior-specification object.
#' @param ... Additional arguments passed to methods (e.g.
#'   \code{max_disp_perc_measurement}). For \code{print.dGamma_list()},
#'   passed to \code{print.pfamily}.
#' @return For \code{dGamma_list()}, an object of class
#'   \code{"dGamma_list"} (also inherits from \code{"list"}): a named
#'   list of \code{"pfamily"} objects (one per group).  Attributes such
#'   as \code{window_diagnostics} and \code{measurement_prior_group}
#'   are retained for downstream use but are not shown by
#'   \code{print.dGamma_list()}.  For \code{print.dGamma_list()},
#'   \code{x} invisibly.
#' @seealso \code{\link[glmbayesCore]{dGamma}}, \code{\link{Prior_Setup_GLMM}},
#'   \code{\link{pfamily_list}}
#' @example inst/examples/Ex_dGamma_list.R
#' @order 1
#' @export
dGamma_list <- function(object, ...) UseMethod("dGamma_list")

#' Resolve names/indices for a named list (dGamma_list / pfamily_list)
#' @noRd
.lmebayes_select_named_list_keys <- function(x, keys, arg = "groups",
                                             what = "name") {
  all_names <- names(x)
  if (is.null(all_names)) {
    all_names <- as.character(seq_along(x))
  }
  if (is.null(keys)) {
    return(all_names)
  }
  if (is.numeric(keys)) {
    if (any(is.na(keys)) || any(keys < 1) || any(keys > length(x)) ||
        any(keys != as.integer(keys))) {
      stop(
        "'", arg, "' numeric indices must be integers in 1:", length(x), ".",
        call. = FALSE
      )
    }
    return(all_names[as.integer(keys)])
  }
  if (!is.character(keys) || length(keys) < 1L || anyNA(keys)) {
    stop(
      "'", arg, "' must be NULL, a character vector of ", what, "s, ",
      "or integer indices.",
      call. = FALSE
    )
  }
  keys <- as.character(keys)
  missing <- setdiff(keys, all_names)
  if (length(missing) > 0L) {
    stop(
      "'", arg, "' contains unknown ", what, "(s): ",
      paste(missing, collapse = ", "),
      ".\nAvailable: ", paste(all_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  keys
}

#' @rdname dGamma_list
#' @method print dGamma_list
#' @param x An object of class \code{"dGamma_list"}.
#' @param groups For \code{print.dGamma_list()}: \code{NULL} (default;
#'   print all groups), a character vector of group names
#'   (\code{names(x)}), or integer indices into \code{x}. Each selected
#'   group's \code{dGamma()} is printed with
#'   \code{\link[glmbayesCore]{print.pfamily}}. List attributes
#'   (\code{window_diagnostics}, \code{group.dispersion.fit}, etc.) are
#'   omitted; use \code{attr(x, "...")} to inspect them.
#' @order 3
#' @export
print.dGamma_list <- function(x, groups = NULL, ...) {
  sel <- .lmebayes_select_named_list_keys(
    x, groups, arg = "groups", what = "group name"
  )
  for (nm in sel) {
    cat(sprintf("\n[[%s]]\n", nm))
    print(x[[nm]], ...)
  }
  invisible(x)
}
