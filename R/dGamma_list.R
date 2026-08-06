#' Build a named list of dGamma measurement-dispersion priors
#'
#' Generic for constructing a named list of \code{\link[glmbayesCore]{dGamma}} prior
#' objects from a prior-specification object.  For mixed models, one
#' \code{dGamma()} per group level supplies observation-level
#' \eqn{\sigma^2_j} priors for \code{lmerb(..., group.dispersion = ...)}.
#'
#' @param object A prior-specification object.
#' @param ... Additional arguments passed to methods (e.g.
#'   \code{max_disp_perc_measurement}).
#' @return An object of class \code{"dGamma_list"} (also inherits from
#'   \code{"list"}): a named list of \code{"pfamily"} objects (one per
#'   group).  Attributes such as \code{window_diagnostics} and
#'   \code{measurement_prior_group} are retained for downstream use but
#'   are not shown by \code{\link{print.dGamma_list}}.
#' @seealso \code{\link[glmbayesCore]{dGamma}}, \code{\link{Prior_Setup_GLMM}},
#'   \code{\link{pfamily_list}}, \code{\link{print.dGamma_list}}
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

#' Subset a dGamma_list (keeps class and trims per-group attributes)
#'
#' @param x An object of class \code{"dGamma_list"}.
#' @param i Character group names, integer indices, or a logical vector.
#' @param j Ignored (list subsetting).
#' @param ... Ignored.
#' @param drop Ignored (always returns a list).
#' @return A \code{"dGamma_list"} with the selected groups.
#' @export
#' @method [ dGamma_list
`[.dGamma_list` <- function(x, i, j, ..., drop = TRUE) {
  if (!missing(j)) {
    stop("dGamma_list objects are one-dimensional; do not use [i, j].",
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
      x, i, arg = "i", what = "group name"
    )
  }

  out <- unclass(x)[keep]
  ## Preserve / trim attributes
  a <- attributes(x)
  a$names <- keep
  a$class <- c("dGamma_list", "list")

  wd <- a$window_diagnostics
  if (is.data.frame(wd) && "group" %in% names(wd)) {
    a$window_diagnostics <- wd[as.character(wd$group) %in% keep, , drop = FALSE]
    rownames(a$window_diagnostics) <- NULL
  }

  mpg <- a$measurement_prior_group
  if (is.list(mpg)) {
    a$measurement_prior_group <- lapply(mpg, function(v) {
      if (!is.null(names(v))) v[keep] else v
    })
  }

  attributes(out) <- a
  out
}

#' Print a dGamma_list
#'
#' Prints each selected group's \code{dGamma()} prior with
#' \code{\link[glmbayesCore]{print.pfamily}}, one after another.
#' List attributes (\code{window_diagnostics}, \code{group.dispersion.fit},
#' etc.) are omitted from the printout; use \code{attr(x, "...")} to
#' inspect them. Group names are \code{names(x)}.
#'
#' @param x An object of class \code{"dGamma_list"}.
#' @param groups \code{NULL} (default; print all groups), a character
#'   vector of group names (\code{names(x)}), or integer indices into
#'   \code{x}.
#' @param ... Passed to \code{print.pfamily}.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{dGamma_list}}, \code{\link[glmbayesCore]{pfamily}}
#' @export
#' @method print dGamma_list
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
