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

#' @noRd
.ing_n_prior_from_shape <- function(...) {
  glmbayesCore:::.ing_n_prior_from_shape(...)
}

#' @noRd
.ing_stop_if_prior_exceeds_data <- function(...) {
  glmbayesCore:::.ing_stop_if_prior_exceeds_data(...)
}

## glmmTMB / merMod reference-fit extractors (Tier 2; live in glmbayesCore).

#' @noRd
.lmebayes_dispformula_kind <- function(...) {
  glmbayesCore:::.lmebayes_dispformula_kind(...)
}

#' @noRd
.lmebayes_fit_glmmtmb_reference <- function(...) {
  glmbayesCore:::.lmebayes_fit_glmmtmb_reference(...)
}

#' @noRd
.lmebayes_glmmtmb_convergence_issues <- function(...) {
  glmbayesCore:::.lmebayes_glmmtmb_convergence_issues(...)
}

#' @noRd
.lmebayes_reference_fixef <- function(...) {
  glmbayesCore:::.lmebayes_reference_fixef(...)
}

#' @noRd
.lmebayes_reference_vcov <- function(...) {
  glmbayesCore:::.lmebayes_reference_vcov(...)
}

#' @noRd
.lmebayes_reference_coef <- function(...) {
  glmbayesCore:::.lmebayes_reference_coef(...)
}

#' @noRd
.lmebayes_extract_reference_variance_components <- function(...) {
  glmbayesCore:::.lmebayes_extract_reference_variance_components(...)
}

#' @noRd
extract_glmmtmb_variance_components <- function(...) {
  glmbayesCore:::extract_glmmtmb_variance_components(...)
}

#' @noRd
.lmebayes_glmmtmb_group_sigma2 <- function(...) {
  glmbayesCore:::.lmebayes_glmmtmb_group_sigma2(...)
}
