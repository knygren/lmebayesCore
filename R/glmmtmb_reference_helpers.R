## Internal glmmTMB reference-fit helpers for Prior_Setup_lmebayes() when
## dispformula requests per-group measurement dispersion. These fit and read
## from a glmmTMB::glmmTMB() reference model as an alternative to the lme4
## merMod reference used for dispformula = ~1; see
## inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md for the calibration this feeds.

#' Fit a \code{glmmTMB::glmmTMB()} reference model with per-group residual
#' dispersion.
#'
#' Used by \code{\link{Prior_Setup_lmebayes}} as the calibration reference
#' (fixed effects, RE variances, per-group coefficients and dispersion) when
#' \code{dispformula} requests per-group measurement dispersion. The
#' lme4-embedded \code{lmer}/\code{glmer} reference fit
#' (\code{design$lmer_fit} / \code{design$glmer_fit}) is fit separately and
#' unconditionally by \code{\link{model_setup}} for backward compatibility.
#'
#' @param formula Mixed-model formula (same as \code{Prior_Setup_lmebayes}).
#' @param data Data frame.
#' @param family A \code{\link[stats]{family}} object.
#' @param dispformula One-sided formula for \code{glmmTMB}'s dispersion model.
#' @param REML Logical, passed to \code{glmmTMB}.
#' @param control Optional \code{glmmTMB::glmmTMBControl()}; omitted (glmmTMB
#'   default) when \code{NULL}.
#' @param ... Passed to \code{glmmTMB::glmmTMB()}.
#' @return A fitted \code{"glmmTMB"} object.
#' @keywords internal
#' @noRd
.lmebayes_fit_glmmtmb_reference <- function(
    formula, data, family, dispformula, REML = TRUE, control = NULL, ...
) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required to calibrate Prior_Setup_lmebayes() ",
      "priors when 'dispformula' requests per-group measurement dispersion. ",
      "Install it with install.packages(\"glmmTMB\").",
      call. = FALSE
    )
  }
  tmb_args <- c(
    list(
      formula     = formula,
      data        = data,
      family      = family,
      dispformula = dispformula,
      REML        = isTRUE(REML)
    ),
    if (!is.null(control)) list(control = control),
    list(...)
  )
  do.call(glmmTMB::glmmTMB, tmb_args)
}

## Return character issue messages when a glmmTMB fit failed to converge
## (nlminb convergence code != 0) or has a non-positive-definite Hessian.
## Empty character() = OK. Mirrors .lmebayes_mer_convergence_issues().
#' @noRd
.lmebayes_glmmtmb_convergence_issues <- function(fit, label = "reference fit") {
  if (is.null(fit)) {
    return(sprintf("%s: fit is NULL", label))
  }
  if (!inherits(fit, "glmmTMB")) {
    return(sprintf("%s: not a glmmTMB object", label))
  }
  issues <- character(0)
  conv_code <- fit$fit$convergence
  if (!is.null(conv_code) && !identical(as.numeric(conv_code), 0)) {
    issues <- c(
      issues,
      sprintf("%s: optimizer did not converge (convergence = %s)", label, conv_code)
    )
  }
  pd_hess <- fit$sdr$pdHess
  if (!is.null(pd_hess) && !isTRUE(pd_hess)) {
    issues <- c(
      issues,
      sprintf("%s: Hessian is not positive-definite (pdHess = FALSE)", label)
    )
  }
  issues
}

## Dispatch helpers so downstream calibration code can treat a merMod or a
## glmmTMB reference fit uniformly.
#' @noRd
.lmebayes_reference_fixef <- function(fit) {
  if (inherits(fit, "glmmTMB")) {
    glmmTMB::fixef(fit)$cond
  } else {
    lme4::fixef(fit)
  }
}

#' @noRd
.lmebayes_reference_vcov <- function(fit) {
  ## glmmTMB does not export its own vcov(); vcov.glmmTMB is a registered
  ## S3 method on the stats::vcov generic (returns a $cond/$disp/... list).
  if (inherits(fit, "glmmTMB")) {
    as.matrix(stats::vcov(fit)$cond)
  } else {
    as.matrix(stats::vcov(fit))
  }
}

## coef(merMod) and coef.glmmTMB(fit)$cond both return a named list, one
## data.frame per grouping factor, in the same layout (rows = levels, columns
## = RE coefficient names); this returns that list uniformly. coef.glmmTMB is
## a registered S3 method on the stats::coef generic, not a glmmTMB export.
#' @noRd
.lmebayes_reference_coef <- function(fit) {
  if (inherits(fit, "glmmTMB")) {
    stats::coef(fit)$cond
  } else {
    stats::coef(fit)
  }
}

#' Extract random-effect variance components from a glmmTMB reference fit
#'
#' Mirrors \code{\link{extract_mer_variance_components}} for a
#' \code{glmmTMB} fit: reads \code{RE} standard deviations from
#' \code{glmmTMB::VarCorr(fit)$cond[[group_name]]}'s \code{"stddev"}
#' attribute. The residual (observation-level) variance is not a single
#' scalar under a per-group \code{dispformula} and is reported as \code{NA}
#' here; see \code{\link{.lmebayes_glmmtmb_group_sigma2}} for the per-group
#' values.
#'
#' @param fit A fitted \code{"glmmTMB"} object.
#' @param re_coef_names Random coefficient names (as in
#'   \code{design$re_coef_names}).
#' @param group_name Name of the random-effects grouping factor.
#' @return List with \code{varcorr} (raw \code{VarCorr.glmmTMB} object),
#'   \code{vcov_re} (named numeric vector, one entry per \code{re_coef_names}),
#'   and \code{residual_var} (\code{NA_real_}).
#' @keywords internal
#' @noRd
extract_glmmtmb_variance_components <- function(fit, re_coef_names, group_name) {
  if (!inherits(fit, "glmmTMB")) {
    stop("fit must be a glmmTMB object.", call. = FALSE)
  }

  vc    <- glmmTMB::VarCorr(fit)
  block <- vc$cond[[group_name]]
  if (is.null(block)) {
    stop(
      "Could not find grouping factor '", group_name,
      "' in glmmTMB VarCorr(fit)$cond.",
      call. = FALSE
    )
  }

  sd_vec <- attr(block, "stddev")
  missing_coefs <- setdiff(re_coef_names, names(sd_vec))
  if (length(missing_coefs) > 0L) {
    stop(
      "Could not find variance components for: ",
      paste(missing_coefs, collapse = ", "),
      call. = FALSE
    )
  }

  vcov_re <- (sd_vec[re_coef_names])^2
  names(vcov_re) <- re_coef_names

  list(varcorr = vc, vcov_re = vcov_re, residual_var = NA_real_)
}

#' Per-group observation-level dispersion from a glmmTMB reference fit
#'
#' Reads \code{predict(fit, type = "disp")} (constant within each level of
#' \code{group_name} when \code{dispformula = ~group_name}) and aggregates by
#' group.
#'
#' @param fit A fitted \code{"glmmTMB"} object with a per-group
#'   \code{dispformula}.
#' @param group_name Name of the random-effects grouping factor.
#' @param group_levels Character vector of group levels, in the desired
#'   output order.
#' @return Named numeric vector (length \code{length(group_levels)}) of
#'   per-group residual variances \eqn{\sigma^2_j}.
#' @keywords internal
#' @noRd
.lmebayes_glmmtmb_group_sigma2 <- function(fit, group_name, group_levels) {
  if (!inherits(fit, "glmmTMB")) {
    stop("fit must be a glmmTMB object.", call. = FALSE)
  }

  disp_hat <- stats::predict(fit, type = "disp")
  grp      <- as.character(stats::model.frame(fit)[[group_name]])

  sigma2_group <- vapply(group_levels, function(lev) {
    vals <- disp_hat[grp == lev]
    if (length(vals) == 0L) {
      return(NA_real_)
    }
    if (diff(range(vals)) > sqrt(.Machine$double.eps) * max(abs(vals), 1)) {
      warning(
        "glmmTMB per-observation dispersion is not constant within group '",
        lev, "'; using the mean.",
        call. = FALSE
      )
    }
    mean(vals)
  }, numeric(1L))

  stats::setNames(sigma2_group, group_levels)
}
