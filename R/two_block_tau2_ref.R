#' Plug-in tau^2_k from a Block~2 pfamily prior spec (not simulation state)
#'
#' For \code{dNormal}, returns fixed \code{dispersion}.  For ING, returns
#' \code{rate / shape} = \eqn{1/E[1/\tau^2]} (harmonic-mean / precision-mean
#' plug-in), since \eqn{\tau^2} enters Block~1 and Block~2 conditionals through
#' precision \eqn{1/\tau^2}.
#' @noRd
.two_block_tau2_ref_from_pfamily <- function(pf) {
  pl <- pf$prior_list
  if (identical(pf$pfamily, "dNormal")) {
    return(as.numeric(pl$dispersion))
  }
  if (identical(pf$pfamily, "dIndependent_Normal_Gamma")) {
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (is.finite(shape) && shape > 0 && is.finite(rate) && rate > 0) {
      return(rate / shape)
    }
    stop(
      "pfamily dIndependent_Normal_Gamma requires positive shape and rate ",
      "so the precision-mean plug-in tau^2 = rate/shape is defined.",
      call. = FALSE
    )
  }
  stop("Unsupported pfamily: ", pf$pfamily, call. = FALSE)
}

#' Named plug-in tau^2 vector from a validated pfamily_list
#' @noRd
.two_block_tau2_ref_vector <- function(pfamily_list, re_names) {
  stats::setNames(
    vapply(re_names, function(k) {
      .two_block_tau2_ref_from_pfamily(pfamily_list[[k]])
    }, numeric(1)),
    re_names
  )
}

#' Plug-in tau^2 vector from pilot dispersion draws (precision-mean plug-in)
#'
#' Returns \code{1 / colMeans(1 / dispersion_fixef_draws)} per RE component —
#' the sample analogue of \eqn{1/E[1/\tau^2]} for main-stage chain starts.
#' @noRd
.two_block_tau2_start_from_dispersion_draws <- function(
    dispersion_fixef_draws,
    re_names = colnames(dispersion_fixef_draws)
) {
  if (is.null(dispersion_fixef_draws) || !is.matrix(dispersion_fixef_draws)) {
    stop(
      "'dispersion_fixef_draws' must be a matrix of positive tau^2 values.",
      call. = FALSE
    )
  }
  if (is.null(re_names) || length(re_names) != ncol(dispersion_fixef_draws)) {
    stop(
      "'re_names' must match ncol(dispersion_fixef_draws).",
      call. = FALSE
    )
  }
  tau2 <- vapply(seq_along(re_names), function(j) {
    x <- dispersion_fixef_draws[, j]
    if (!all(is.finite(x)) || any(x <= 0)) {
      stop(
        "dispersion_fixef_draws[*, ", re_names[j],
        "] must be finite and strictly positive.",
        call. = FALSE
      )
    }
    1 / mean(1 / x)
  }, numeric(1))
  stats::setNames(tau2, re_names)
}
