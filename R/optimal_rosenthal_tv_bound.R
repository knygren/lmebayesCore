## Optimal Rosenthal TV tuning (category 8).

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

#' @noRd
.c05_rosenthal_sweeps_for_tol <- function(tol,
                                          eps,
                                          eigenvalues,
                                          drift,
                                          display_mode = "sharp",
                                          mode = NULL,
                                          d = NULL,
                                          alpha = NULL,
                                          optimize_alpha = TRUE) {
  if (!(tol > 0 && tol < 1)) {
    stop("'tol' must lie in (0, 1).", call. = FALSE)
  }
  k <- 1L
  repeat {
    a_use <- alpha
    if (is.null(a_use) && isTRUE(optimize_alpha)) {
      opt <- .c05_rosenthal_optimize_alpha(k, eps, drift)
      a_use <- opt$alpha
    } else if (is.null(a_use)) {
      stop("'alpha' is required when optimize_alpha = FALSE.", call. = FALSE)
    }
    ros <- rosenthal_tv_bound(
      k = k,
      alpha = a_use,
      eps = eps,
      eigenvalues = eigenvalues,
      drift = drift,
      mode = mode,
      display_mode = display_mode,
      d = d
    )
    if (ros$bound <= tol) {
      return(list(k = k, bound = ros$bound, rosenthal = ros, alpha = a_use))
    }
    if (k > 1e7) {
      stop("Could not reach 'tol' within 1e7 sweeps.", call. = FALSE)
    }
    k <- k + 1L
  }
}

#' Optimize Rosenthal TV bound over alpha and/or sweep count k.
#'
#' @param mode A \code{\link{population_mode}} result.
#' @param eigenvalues Output of \code{\link{floor_coupling_eigenvalues}}.
#' @param eps Minorization constant for the bound.
#' @param k If set, evaluate at this sweep count only.
#' @param alpha If set, skip alpha optimization.
#' @param inner_tol Target for the inner Rosenthal bound.
#' @param total_tol Target for full \eqn{\pi_\gamma} bound when
#'   \code{include_full_pi_gamma = TRUE} (uses \code{total_tol - delta_2}).
#' @param delta_2 Asymptotic tail allowance for the full marginal.
#' @param include_full_pi_gamma If \code{TRUE}, set \code{full_bound} to inner
#'   bound plus \code{delta_2}.
#' @param optimize_alpha If \code{TRUE}, minimize over \code{alpha}.
#' @inheritParams rosenthal_drift_constants
#' @return A list of class \code{"optimal_rosenthal_tv_bound"}.
#' @seealso \code{\link{rosenthal_tv_bound}}, \code{\link{gamma_beta_tv_certificate}}
#' @export
optimal_rosenthal_tv_bound <- function(mode,
                                       eigenvalues,
                                       eps,
                                       k = NULL,
                                       alpha = NULL,
                                       inner_tol = NULL,
                                       total_tol = NULL,
                                       delta_2 = 0,
                                       include_full_pi_gamma = TRUE,
                                       optimize_alpha = TRUE,
                                       display_mode = c("sharp", "general"),
                                       d = NULL,
                                       V_sup = NULL) {
  display_mode <- match.arg(display_mode)
  if (!(eps > 0 && eps <= 1)) {
    stop("'eps' must lie in (0, 1].", call. = FALSE)
  }

  drift <- rosenthal_drift_constants(
    eigenvalues = eigenvalues,
    mode = mode,
    display_mode = display_mode,
    d = d,
    V_sup = V_sup
  )

  target <- inner_tol
  if (is.null(target) && !is.null(total_tol)) {
    if (!(total_tol > delta_2)) {
      stop("'total_tol' must exceed 'delta_2'.", call. = FALSE)
    }
    target <- total_tol - delta_2
  }

  rosenthal <- NULL
  sweeps <- NULL

  if (!is.null(k)) {
    if (is.null(alpha) && isTRUE(optimize_alpha)) {
      opt <- .c05_rosenthal_optimize_alpha(k, eps, drift)
      alpha <- opt$alpha
    } else if (is.null(alpha)) {
      stop("'alpha' is required when optimize_alpha = FALSE.", call. = FALSE)
    }
    rosenthal <- rosenthal_tv_bound(
      k = k,
      alpha = alpha,
      eps = eps,
      eigenvalues = eigenvalues,
      drift = drift,
      mode = mode,
      display_mode = display_mode,
      d = d,
      V_sup = V_sup
    )
  } else if (!is.null(target)) {
    sweeps <- .c05_rosenthal_sweeps_for_tol(
      tol = target,
      eps = eps,
      eigenvalues = eigenvalues,
      drift = drift,
      display_mode = display_mode,
      mode = mode,
      d = d,
      alpha = alpha,
      optimize_alpha = optimize_alpha
    )
    rosenthal <- sweeps$rosenthal
    k <- sweeps$k
    alpha <- sweeps$alpha
  } else {
    stop("Set 'k' or 'inner_tol' / 'total_tol'.", call. = FALSE)
  }

  inner_bound <- rosenthal$bound
  full_bound <- if (isTRUE(include_full_pi_gamma)) {
    inner_bound + delta_2
  } else {
    inner_bound
  }

  structure(
    list(
      k = k,
      alpha = alpha,
      eps = eps,
      inner_tol = target,
      total_tol = total_tol,
      delta_2 = delta_2,
      rosenthal = rosenthal,
      drift = drift,
      eigenvalues = eigenvalues,
      inner_bound = inner_bound,
      full_bound = full_bound,
      include_full_pi_gamma = isTRUE(include_full_pi_gamma),
      display_mode = display_mode,
      sweeps = sweeps,
      call = match.call()
    ),
    class = "optimal_rosenthal_tv_bound"
  )
}

#' @export
print.optimal_rosenthal_tv_bound <- function(x, digits = 4, ...) {
  cat("Optimal Rosenthal TV bound\n")
  cat("  k:", x$k, "  alpha:", signif(x$alpha, digits), "\n")
  cat("  inner bound:", signif(x$inner_bound, digits), "\n")
  if (isTRUE(x$include_full_pi_gamma)) {
    cat("  full pi_gamma bound (inner + delta_2):",
        signif(x$full_bound, digits), "\n")
  }
  invisible(x)
}
