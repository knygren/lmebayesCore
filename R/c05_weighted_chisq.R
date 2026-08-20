## Weighted chi-square tails for C05 deficiency calibration.

#' Validate weights for weighted chi-square routines.
#' @noRd
.c05_validate_weights <- function(weights, fn_name = "weighted chi-square") {
  weights <- as.numeric(weights)
  if (length(weights) < 1L) {
    stop("'", fn_name, "' requires at least one positive weight.", call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("'", fn_name, "' weights must be finite and non-negative.", call. = FALSE)
  }
  weights <- weights[weights > 0]
  if (length(weights) < 1L) {
    stop("'", fn_name, "' requires at least one strictly positive weight.",
         call. = FALSE)
  }
  weights
}

#' Upper tail of a weighted sum of central chi-square_1 variates.
#' @noRd
.c05_quadratic_form_survival <- function(q, weights) {
  if (requireNamespace("CompQuadForm", quietly = TRUE)) {
    return(
      CompQuadForm::imhof(
        q,
        weights,
        epsabs = 1e-8,
        epsrel = 1e-8,
        limit = 10000L
      )$Qq
    )
  }
  if (requireNamespace("mgcv", quietly = TRUE)) {
    return(
      mgcv::psum.chisq(
        q,
        weights,
        lower.tail = FALSE,
        tol = 1e-8,
        nlim = 100000L
      )
    )
  }
  stop(
    "Weighted chi-square tails require suggested package 'CompQuadForm' ",
    "(or 'mgcv').",
    call. = FALSE
  )
}

#' Lower tail Pr(sum weights_i Z_i^2 <= q).
#' @noRd
.c05_weighted_chisq_cdf <- function(q, weights) {
  if (!is.finite(q)) {
    stop("'q' must be finite.", call. = FALSE)
  }
  weights <- .c05_validate_weights(weights)
  if (q <= 0) {
    return(0)
  }

  if (length(weights) == 1L) {
    return(stats::pchisq(q / weights, df = 1, lower.tail = TRUE))
  }

  wu <- unique(weights)
  if (length(wu) == 1L) {
    return(stats::pchisq(q / wu, df = length(weights), lower.tail = TRUE))
  }

  1 - .c05_quadratic_form_survival(q, weights)
}

#' Upper tail Pr(sum weights_i Z_i^2 > q).
#' @noRd
.c05_weighted_chisq_tail <- function(q, weights) {
  if (!is.finite(q)) {
    stop("'q' must be finite.", call. = FALSE)
  }
  weights <- .c05_validate_weights(weights)
  if (q <= 0) {
    return(1)
  }

  if (length(weights) == 1L) {
    return(stats::pchisq(q / weights, df = 1, lower.tail = FALSE))
  }

  wu <- unique(weights)
  if (length(wu) == 1L) {
    return(stats::pchisq(q / wu, df = length(weights), lower.tail = FALSE))
  }

  .c05_quadratic_form_survival(q, weights)
}

#' Validate delta in (0, 1).
#' @noRd
.c05_validate_delta <- function(delta) {
  if (!(delta > 0 && delta < 1)) {
    stop("'delta' must lie in (0, 1).", call. = FALSE)
  }
}

#' Invert Pr(sum w_i Z_i^2 > r^2) = delta for reference radius r and d = r^2/2.
#'
#' Weights \code{w_i} are the \eqn{B}-spectrum (\eqn{\kappa_i/(1-\kappa_i)}).
#' Univariate closed form: \code{r^2 = w * qchisq(1-delta, 1)}, \code{d = r^2/2}.
#' @noRd
.c05_r_from_delta <- function(delta, weights) {
  .c05_validate_delta(delta)
  weights <- .c05_validate_weights(weights, fn_name = "r_from_delta")

  if (length(weights) == 1L) {
    r2 <- weights[[1L]] * stats::qchisq(1 - delta, df = 1)
    return(list(r2 = r2, r = sqrt(r2), d = 0.5 * r2))
  }

  wu <- unique(weights)
  if (length(wu) == 1L) {
    r2 <- wu[[1L]] * stats::qchisq(1 - delta, df = length(weights))
    return(list(r2 = r2, r = sqrt(r2), d = 0.5 * r2))
  }

  tail_at <- function(r2) {
    .c05_weighted_chisq_tail(r2, weights)
  }

  r2_hi <- max(weights) * stats::qchisq(1 - delta, df = 1)
  if (!is.finite(r2_hi) || r2_hi <= 0) {
    r2_hi <- max(weights)
  }
  while (tail_at(r2_hi) > delta && r2_hi < 1e300) {
    r2_hi <- r2_hi * 2
  }
  if (tail_at(r2_hi) > delta) {
    stop(
      "Could not bracket reference radius r for delta = ", delta, ".",
      call. = FALSE
    )
  }

  r2 <- stats::uniroot(
    function(r2) tail_at(r2) - delta,
    lower = 0,
    upper = r2_hi,
    tol = .Machine$double.eps^0.5
  )$root

  list(r2 = r2, r = sqrt(r2), d = 0.5 * r2)
}
