#' Resolve names/indices for a named list (dGamma_list / pfamily_list printers)
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
    idx <- as.integer(keys)
    if (anyNA(idx) || any(idx < 1L) || any(idx > length(x))) {
      stop(
        sprintf("Invalid %s index in '%s'.", what, arg),
        call. = FALSE
      )
    }
    return(all_names[idx])
  }
  if (!is.character(keys)) {
    stop(sprintf("'%s' must be NULL, character, or numeric.", arg), call. = FALSE)
  }
  bad <- setdiff(keys, all_names)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "Unknown %s(s) in '%s': %s.",
        what, arg, paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  keys
}
