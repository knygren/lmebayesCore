## summary() method for rglmerb (glmerb only); shares summary.rlmerb().
## Must collate after summary.rlmerb.R so the shared Rd keeps \name{summary.rlmerb}.

#' @rdname summary.rlmerb
#' @export
#' @method summary rglmerb
summary.rglmerb <- function(object, groups = NULL, ...) {
  summary.rlmerb(object, groups = groups, ...)
}
