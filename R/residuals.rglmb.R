#' Deviance Residuals for \code{rglmb} and \code{summary.rglmb} Objects
#'
#' Returns a matrix of deviance residuals across posterior draws, using the
#' fitted object's \code{family$dev.resids} function as in
#' \code{\link[stats]{residuals.glm}}.
#'
#' @param object an object of class \code{rglmb}, \code{rlmb}, or
#'   \code{summary.rglmb}.
#' @param ysim optional simulated responses (matrix with one row per draw).
#' @param ... further arguments (currently unused).
#' @return A numeric matrix of deviance residuals with one row per draw.
#' @seealso \code{\link{rglmb}}, \code{\link{summary.rglmb}},
#'   \code{\link[stats]{residuals.glm}}.
#' @export
#' @method residuals rglmb
residuals.rglmb <- function(object, ysim = NULL, ...) {
  .residuals_rglmb_draws(object, ysim = ysim)
}

#' @rdname residuals.rglmb
#' @export
#' @method residuals summary.rglmb
residuals.summary.rglmb <- function(object, ysim = NULL, ...) {
  .residuals_rglmb_draws(object, ysim = ysim)
}

#' @rdname residuals.rglmb
#' @export
#' @method residuals rlmb
residuals.rlmb <- function(object, ysim = NULL, ...) {
  .residuals_rglmb_draws(object, ysim = ysim)
}

.residuals_rglmb_draws <- function(object, ysim = NULL) {
  y <- object$y
  n <- nrow(object$coefficients)
  wts <- object$prior.weights

  if (!is.null(object$fitted.values)) {
    fv_mat <- object$fitted.values
  } else {
    lp_mat <- t(object$x %*% t(object$coefficients))
    fv_mat <- object$family$linkinv(lp_mat)
  }

  devfun <- object$family$dev.resids
  DevRes <- matrix(0, nrow = n, ncol = length(y))

  for (i in seq_len(n)) {
    mu_vec <- if (is.null(ysim)) fv_mat[i, ] else ysim[i, ]
    DevRes[i, ] <- sign(y - mu_vec) * sqrt(devfun(y, mu_vec, wts))
  }

  colnames(DevRes) <- names(y)
  DevRes
}
