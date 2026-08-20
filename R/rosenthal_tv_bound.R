## Rosenthal drift constants and TV bound (categories 7).

#' Foster / Rosenthal drift constants from floor coupling eigenvalues.
#'
#' @param eigenvalues Output of \code{\link{floor_coupling_eigenvalues}} (or any
#'   list with \code{kappa_lb} and \code{q}).
#' @param display_mode \code{"sharp"} (default) or \code{"general"}.
#' @param gamma_0 Start state for \eqn{V(\gamma_0)} in general mode.
#' @param mode Optional \code{\link{population_mode}} result for general mode.
#' @param d Deficiency level for \eqn{V_{\sup}(d)} when \code{V_sup} is omitted.
#' @param V_sup Override for \eqn{V_{\sup}} in general mode.
#' @return A list with \code{kappa_lb}, \code{lambda_lb}, \code{b_drift},
#'   \code{C_beta_plus}, \code{U}, drift numerator/denominator, and Lyapunov
#'   values.
#' @seealso \code{\link{rosenthal_tv_bound}}, \code{\link{floor_coupling_eigenvalues}}
#' @export
rosenthal_drift_constants <- function(eigenvalues,
                                      gamma_0 = NULL,
                                      mode = NULL,
                                      display_mode = c("sharp", "general"),
                                      d = NULL,
                                      V_sup = NULL) {
  display_mode <- match.arg(display_mode)
  if (is.null(eigenvalues$kappa_lb)) {
    stop("'eigenvalues' must contain 'kappa_lb'.", call. = FALSE)
  }

  kappa_lb <- eigenvalues$kappa_lb
  kappa_max <- if (is.null(eigenvalues$kappa_max_lb)) {
    max(kappa_lb)
  } else {
    eigenvalues$kappa_max_lb
  }
  q <- if (is.null(eigenvalues$q)) length(kappa_lb) else eigenvalues$q
  lambda <- kappa_max^2
  C_beta_plus <- 0.5 * sum(kappa_lb)
  b <- 1 - lambda + q / 2 + C_beta_plus

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
    C_beta_plus = C_beta_plus,
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

#' Rosenthal total-variation bound at fixed tuning parameter alpha.
#'
#' Evaluates the Rosenthal (1995) bound at sweep count \code{k}, minorization
#' \code{eps}, floor eigenvalues, and fixed \code{alpha}. To optimize
#' \code{alpha} or \code{k}, use \code{\link{optimal_rosenthal_tv_bound}}.
#'
#' @param k Number of Gibbs sweeps (\eqn{k \ge 1}).
#' @param alpha Tuning parameter in \eqn{(\lambda^{\mathrm{LB}}, 1)} (required).
#' @param eps Minorization constant (typically \code{epsilon_star()$eps_star}).
#' @param eigenvalues Output of \code{\link{floor_coupling_eigenvalues}}.
#' @param drift Optional output of \code{\link{rosenthal_drift_constants}}; if
#'   \code{NULL}, computed from \code{eigenvalues}.
#' @inheritParams rosenthal_drift_constants
#' @return A list with \code{bound}, \code{minorization}, \code{drift},
#'   \code{k}, \code{alpha}, \code{eps}, and \code{drift_constants}.
#' @seealso \code{\link{optimal_rosenthal_tv_bound}}
#' @export
rosenthal_tv_bound <- function(k,
                                 alpha,
                                 eps,
                                 eigenvalues,
                                 drift = NULL,
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
  if (missing(alpha) || is.null(alpha)) {
    stop("'alpha' is required; use optimal_rosenthal_tv_bound() to optimize it.",
         call. = FALSE)
  }
  if (!(eps > 0 && eps <= 1)) {
    stop("'eps' must lie in (0, 1].", call. = FALSE)
  }

  if (is.null(drift)) {
    drift <- rosenthal_drift_constants(
      eigenvalues = eigenvalues,
      gamma_0 = gamma_0,
      mode = mode,
      display_mode = display_mode,
      d = d,
      V_sup = V_sup
    )
  }

  if (!(alpha > drift$lambda_lb && alpha < 1)) {
    stop(
      "'alpha' must lie in (lambda_lb, 1); lambda_lb = ",
      drift$lambda_lb,
      ".",
      call. = FALSE
    )
  }

  bound <- .c05_rosenthal_bound_at_alpha(k, alpha, eps, drift)
  minor <- (1 - eps)^floor(alpha * k)

  list(
    k = k,
    eps = eps,
    alpha = alpha,
    bound = bound,
    minorization = minor,
    drift = bound - minor,
    drift_constants = drift,
    display_mode = display_mode,
    lambda_lb = drift$lambda_lb
  )
}
