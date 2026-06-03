#' Directional Tail Diagnostic
#'
#' Computes the directional tail probability based on posterior draws and prior mean,
#' using whitening transformation and projection onto the direction of disagreement.
#' This diagnostic identifies directional disagreement between posterior and prior,
#' and is especially useful for visualizing rejection regions in whitened space.
#' The whitening uses Mahalanobis distance \insertCite{Mahalanobis1936}{glmbayesCore} in
#' posterior-precision-scaled coordinates.
#'
#' @name directional_tail
#' @param fit A fitted model object of class \code{rglmb} or \code{rlmb}
#' @param mu0 An optional argument containing a reference vector relative to which the directional tail is computed. Defaults to the prior mean.
#' @param x An object of class \code{directional_tail}
#' @param ... Additional arguments passed to or from other methods.
#' @return An object of class 'directional_tail' containing:
#'   \item{mahalanobis_shift}{Measures the standardized Mahalanobis distance between the posterior and prior means,
#' using posterior precision for scaling. In the Gaussian case, this directly determines the directional tail probability via Phi(-||w||).}
#'   \item{p_directional}{Directional tail probability (proportion of draws in the direction of disagreement)}
#'   \item{delta}{Mean shift in whitened space}
#'   \item{draws}{List containing whitened draws, raw draws, and tail flags}
#' @details Whitening is performed using the posterior precision matrix.
#' The direction vector is computed as the mean shift in whitened space.
#' Tail probability is the proportion of draws with negative projection onto this direction.
#' For theory, interpretation, and relation to t/F statistics, see
#' \insertCite{glmbayesChapterA04}{glmbayesCore}.
#' @seealso \code{\link{rglmb}}, \code{\link{rlmb}}, \code{\link{directional_tail}}
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_directional_tail.R
#' @export
#' @keywords diagnostic geometry Bayesian

directional_tail <- function(fit, mu0 = NULL) {
  main_class <- class(fit)[1]

  B <- as.matrix(fit$coefficients)
  if (is.null(mu0)) {
    mu0 <- as.numeric(fit$Prior$mean)
  }

  if (!is.null(fit$Prior$Precision)) {
    Prec_prior <- fit$Prior$Precision
    Prec_prior <- 0.5 * (Prec_prior + t(Prec_prior))
  } else if (!is.null(fit$Prior$Variance)) {
    V0 <- fit$Prior$Variance
    Prec_prior <- chol2inv(chol(V0))
    Prec_prior <- 0.5 * (Prec_prior + t(Prec_prior))
  } else if (!is.null(fit$Prior$Sigma)) {
    V0 <- fit$Prior$Sigma
    Prec_prior <- chol2inv(chol(V0))
    Prec_prior <- 0.5 * (Prec_prior + t(Prec_prior))
  } else {
    stop("Could not recover prior precision from fit$Prior.", call. = FALSE)
  }

  if (main_class == "rglmb") {
    f7 <- fit$famfunc$f7
    bstar <- fit$coef.mode
    y <- fit$y
    x <- fit$x
    P0 <- Prec_prior
    alpha <- tryCatch(fit$offset, error = function(e) NULL)
    if (is.null(alpha)) alpha <- 0
    wt <- fit$prior.weights
    Prec_lik <- f7(bstar, y, x, mu0, P0, alpha, wt)
  } else if (main_class == "rlmb") {
    x <- fit$x
    wt <- fit$prior.weights
    disp <- fit$dispersion
    if (length(disp) > 1L) {
      disp <- disp[1L]
    }
    Xw <- sqrt(wt) * x
    Prec_lik <- crossprod(Xw) / disp
  } else {
    stop("Unsupported model class: expected rglmb or rlmb.", call. = FALSE)
  }

  Prec_post <- Prec_lik + Prec_prior

  P <- Prec_post
  eig <- eigen(P, symmetric = TRUE)
  D_eig <- diag(sqrt(eig$values), nrow = length(eig$values))
  W <- eig$vectors %*% D_eig %*% t(eig$vectors)
  Z <- t(W %*% t(sweep(B, 2, mu0)))

  delta <- colMeans(Z)
  w <- delta

  proj <- as.vector(Z %*% w)
  flag <- proj <= 0
  p_tail <- mean(flag)

  out <- list(
    mahalanobis_shift = sqrt(sum(delta^2)),
    p_directional = p_tail,
    delta = delta,
    draws = list(
      is_tail = flag,
      Z = Z,
      B = B
    )
  )

  class(out) <- "directional_tail"
  out
}

#' @rdname directional_tail
#' @method print directional_tail
#' @export

print.directional_tail <- function(x, ...) {
  n_draws <- nrow(x$draws$Z)

  cat("Bayesian Estimates Based on", n_draws, "iid draws\n")
  cat("--------------------------------------------------\n")
  cat("Standardized Prior-Posterior Mahalanobis distance, and\n")
  cat("associated tail probability (P[delta^T * Z <= 0]):\n\n")
  cat("  Mahalanobis Prior-Posterior Distance    :", formatC(x$mahalanobis_shift, digits = 4, format = "f"), "\n")
  cat("  Associated Directional Tail Probability :", formatC(x$p_directional, digits = 4, format = "f"), "\n\n")

  Z <- x$draws$Z
  mcse_delta <- apply(Z, 2, function(zj) sd(zj) / sqrt(length(zj)))

  cat("Standardized Prior-posterior shifts mcses: \n")

  delta_table <- data.frame(
    Delta = round(x$delta, 4),
    MCSE = round(mcse_delta, 4),
    row.names = colnames(Z)
  )
  print(delta_table)

  cat("\nPosterior draws available:\n")
  cat("  - `x$draws$Z` for whitened draws\n")
  cat("  - `x$draws$B` for raw coefficient draws\n")
  cat("  - `x$draws$is_tail` to flag tail draws\n")

  invisible(x)
}
