## GLM family hooks for group_precision_floor(): NLL, Fisher weights, interval minima.

#' @noRd
.group_floor_l1pe <- function(e) {
  ifelse(e > 30, e, log1p(exp(e)))
}

#' @noRd
.group_floor_family <- function(family) {
  if (!inherits(family, "family")) {
    family <- stats::family(family)
  }
  fam <- family$family
  link <- family$link

  if (identical(fam, "gaussian") && identical(link, "identity")) {
    return(list(
      family = family,
      name = "gaussian",
      endpoint_rule = "constant",
      gaussian = TRUE,
      weight_at = function(eta, size, wt) wt,
      weight_min_on_interval = function(lo, hi, size, wt) wt,
      w_max = function(size, wt) max(wt),
      negll_j = function(y, eta, size, wt, dispersion) {
        0.5 * sum(wt * (y - eta)^2 / dispersion)
      },
      Gmat_j = function(Z, eta, size, wt, dispersion) {
        crossprod(Z, wt * Z) / dispersion
      }
    ))
  }

  if (identical(fam, "binomial") && identical(link, "logit")) {
    wof <- function(e, n) n * exp(-abs(e) - 2 * .group_floor_l1pe(-abs(e)))
    return(list(
      family = family,
      name = "binomial_logit",
      endpoint_rule = "log_concave",
      gaussian = FALSE,
      weight_at = function(eta, size, wt) wof(eta, size * wt),
      weight_min_on_interval = function(lo, hi, size, wt) {
        w1 <- wof(lo, size * wt)
        w2 <- wof(hi, size * wt)
        pmin(w1, w2)
      },
      w_max = function(size, wt) max(size * wt) / 4,
      negll_j = function(y, eta, size, wt, dispersion) {
        n <- size * wt
        -sum(y * eta - n * .group_floor_l1pe(eta))
      },
      Gmat_j = function(Z, eta, size, wt, dispersion) {
        n <- size * wt
        p <- plogis(eta)
        crossprod(Z, (n * p * (1 - p)) * Z)
      }
    ))
  }

  if (identical(fam, "poisson") && identical(link, "log")) {
    return(list(
      family = family,
      name = "poisson_log",
      endpoint_rule = "monotone_inc",
      gaussian = FALSE,
      weight_at = function(eta, size, wt) size * wt * exp(eta),
      weight_min_on_interval = function(lo, hi, size, wt) {
        size * wt * exp(lo)
      },
      w_max = function(size, wt) max(size * wt),
      negll_j = function(y, eta, size, wt, dispersion) {
        n <- size * wt
        -sum(y * eta - n * exp(eta))
      },
      Gmat_j = function(Z, eta, size, wt, dispersion) {
        n <- size * wt
        crossprod(Z, (n * exp(eta)) * Z)
      }
    ))
  }

  stop(
    "group_precision_floor() supports gaussian(identity), ",
    "binomial(logit), and poisson(log) only; got family = ", fam,
    "(), link = ", link, ".",
    call. = FALSE
  )
}
