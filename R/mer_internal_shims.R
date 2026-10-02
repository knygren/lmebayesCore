## Internal helpers still resolved from the lmebayesCore namespace (lmebayes :::, tests).

#' @noRd
.lmebayes_mer_optional_args <- function(...) {
  glmbayesCore:::.lmebayes_mer_optional_args(...)
}

#' @noRd
.lmebayes_normalize_weights <- function(...) {
  glmbayesCore:::.lmebayes_normalize_weights(...)
}

#' @noRd
.lmebayes_normalize_offset <- function(...) {
  glmbayesCore:::.lmebayes_normalize_offset(...)
}

#' @noRd
.lmebayes_stop_if_nondefault_weights_offset <- function(...) {
  glmbayesCore:::.lmebayes_stop_if_nondefault_weights_offset(...)
}

#' @noRd
extract_mer_variance_components <- function(...) {
  glmbayesCore:::extract_mer_variance_components(...)
}

#' @noRd
.lmebayes_ing_prior_is_grouped <- function(x) {
  glmbayesCore:::.lmebayes_ing_prior_is_grouped(x)
}
