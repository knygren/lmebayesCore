## Rosenthal (1995) TV bound for the beta-restricted gamma chain.

#' @noRd
.c05_rosenthal_drift_block <- function(spectrum,
                                       gamma_0 = NULL,
                                       mode = NULL,
                                       display_mode = c("sharp", "general"),
                                       d = NULL,
                                       V_sup = NULL) {
  display_mode <- match.arg(display_mode)
  kappa_lb <- spectrum$kappa_lb
  kappa_max <- spectrum$kappa_max_lb
  lambda <- spectrum$lambda_lb
  b <- spectrum$b_drift
  q <- spectrum$q

  if (identical(display_mode, "sharp")) {
    V0 <- 1
    V_sup_use <- 1
  } else {
    if (is.null(gamma_0) && !is.null(mode)) {
      gamma_0 <- unlist(mode$fixef, use.names = FALSE)
    }
    if (is.null(gamma_0) || is.null(mode) || is.null(mode$p11)) {
      stop(
        "display_mode = \"general\" requires 'gamma_0' and 'mode' with 'p11'.",
        call. = FALSE
      )
    }
    gamma_star <- unlist(mode$fixef, use.names = FALSE)
    diff <- gamma_0 - gamma_star
    V0 <- 1 + 0.5 * as.numeric(
      crossprod(diff, mode$p11$P11 %*% diff)
    )
    if (is.null(V_sup)) {
      if (is.null(d)) {
        stop("'d' or 'V_sup' required for display_mode = \"general\".", call. = FALSE)
      }
      V_sup_use <- 1 + d
    } else {
      V_sup_use <- V_sup
    }
  }

  drift_num <- 1 + 2 * b + lambda * V0
  drift_den <- 1 + 2 * b / (1 - lambda)
  U <- 1 + 2 * b + lambda * V_sup_use

  list(
    kappa_lb = kappa_lb,
    kappa_max_lb = kappa_max,
    lambda_lb = lambda,
    b_drift = b,
    C_beta_plus = spectrum$C_beta_plus,
    q = q,
    V_gamma_0 = V0,
    V_sup = V_sup_use,
    U = U,
    drift_numerator = drift_num,
    drift_denominator = drift_den,
    display_mode = display_mode
  )
}

#' @noRd
.c05_rosenthal_bound_at_alpha <- function(k, alpha, eps, drift) {
  if (!(alpha > drift$lambda_lb && alpha < 1)) {
    return(Inf)
  }
  minor <- (1 - eps)^floor(alpha * k)
  drift_term <- (drift$U / alpha) *
    (drift$drift_numerator / drift$drift_denominator) *
    alpha^k
  minor + drift_term
}

#' @noRd
.c05_rosenthal_optimize_alpha <- function(k, eps, drift, n_grid = 40L) {
  lambda <- drift$lambda_lb
  if (!(lambda < 1)) {
    stop("lambda_lb must be strictly less than 1.", call. = FALSE)
  }
  lo <- lambda + max(sqrt(.Machine$double.eps), (1 - lambda) * 1e-6)
  hi <- 1 - sqrt(.Machine$double.eps)
  grid <- unique(c(seq(lo, hi, length.out = n_grid), (lo + hi) / 2))
  vals <- vapply(grid, function(a) {
    .c05_rosenthal_bound_at_alpha(k, a, eps, drift)
  }, numeric(1L))
  i <- which.min(vals)
  list(alpha = grid[i], bound = vals[i], grid = grid, values = vals)
}

#' Rosenthal total-variation bound for the beta-restricted gamma chain.
#'
#' Evaluates the Rosenthal (1995) bound from
#' \code{inst/GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md} section 3.1 at
#' sweep count \code{k}, minorization \code{eps}, and floor spectrum
#' \code{kappa_lb}. With \code{display_mode = "sharp"}, uses
#' \eqn{\varepsilon(\gamma^\star)}, \eqn{V_{\sup}=V(\gamma^\star)=1}, and
#' \eqn{\gamma_0=\gamma^\star}.
#'
#' @param k Number of Gibbs sweeps (\eqn{k \ge 1}).
#' @param eps Minorization constant (typically \code{epsilon_star()$eps_star}
#'   for sharp display).
#' @param spectrum Output of \code{\link{floor_coupling_spectrum}}.
#' @param alpha Tuning parameter in \eqn{(\lambda^{\mathrm{LB}}, 1)}; if
#'   \code{NULL}, optimized at \code{k}.
#' @param gamma_0 Start state for the drift Lyapunov factor; \code{NULL} uses
#'   \eqn{\gamma^\star} from \code{mode} when supplied.
#' @param mode Optional \code{\link{population_mode}} result for general display.
#' @param display_mode \code{"sharp"} (default) or \code{"general"}.
#' @param d Deficiency level for \eqn{V_{\sup}(d)} in general mode when
#'   \code{V_sup} is not supplied.
#' @param V_sup Override for \eqn{V_{\sup}} in general mode.
#' @return A list with \code{bound}, \code{minorization}, \code{drift},
#'   \code{alpha}, \code{drift_block}, and \code{display_mode}.
#' @seealso \code{\link{gamma_beta_tv_certificate}}, \code{\link{floor_coupling_spectrum}}
#' @export
rosenthal_tv_bound <- function(k,
                               eps,
                               spectrum,
                               alpha = NULL,
                               gamma_0 = NULL,
                               mode = NULL,
                               display_mode = c("sharp", "general"),
                               d = NULL,
                               V_sup = NULL) {
  display_mode <- match.arg(display_mode)
  k <- as.integer(k)
  if (k < 1L) {
    stop("'k' must be at least 1.", call. = FALSE)
  }
  if (!(eps > 0 && eps <= 1)) {
    stop("'eps' must lie in (0, 1].", call. = FALSE)
  }
  if (is.null(spectrum$kappa_lb)) {
    stop("'spectrum' must be from floor_coupling_spectrum().", call. = FALSE)
  }

  drift <- .c05_rosenthal_drift_block(
    spectrum = spectrum,
    gamma_0 = gamma_0,
    mode = mode,
    display_mode = display_mode,
    d = d,
    V_sup = V_sup
  )

  if (is.null(alpha)) {
    opt <- .c05_rosenthal_optimize_alpha(k, eps, drift)
    alpha <- opt$alpha
    bound <- opt$bound
  } else {
    if (!(alpha > drift$lambda_lb && alpha < 1)) {
      stop(
        "'alpha' must lie in (lambda_lb, 1); lambda_lb = ",
        drift$lambda_lb,
        ".",
        call. = FALSE
      )
    }
    bound <- .c05_rosenthal_bound_at_alpha(k, alpha, eps, drift)
  }

  minor <- (1 - eps)^floor(alpha * k)
  drift_val <- bound - minor

  list(
    k = k,
    eps = eps,
    alpha = alpha,
    bound = bound,
    minorization = minor,
    drift = drift_val,
    drift_block = drift,
    display_mode = display_mode,
    lambda_lb = drift$lambda_lb
  )
}

#' Invert Rosenthal bound for sweep count at tolerance (inner bound only).
#' @noRd
.c05_rosenthal_sweeps_for_tol <- function(tol,
                                          eps,
                                          spectrum,
                                          display_mode = "sharp",
                                          mode = NULL,
                                          d = NULL,
                                          alpha = NULL) {
  if (!(tol > 0 && tol < 1)) {
    stop("'tol' must lie in (0, 1).", call. = FALSE)
  }
  k <- 1L
  repeat {
    ros <- rosenthal_tv_bound(
      k = k,
      eps = eps,
      spectrum = spectrum,
      alpha = alpha,
      mode = mode,
      display_mode = display_mode,
      d = d
    )
    if (ros$bound <= tol) {
      return(list(k = k, bound = ros$bound, rosenthal = ros))
    }
    if (k > 1e7) {
      stop("Could not reach 'tol' within 1e7 sweeps.", call. = FALSE)
    }
    k <- k + 1L
  }
}
