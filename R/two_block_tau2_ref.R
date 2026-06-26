#' Plug-in tau^2_k from a Block~2 pfamily prior spec (not simulation state)
#'
#' For \code{dNormal}, returns fixed \code{dispersion}.  For ING, returns the
#' inverse-Gamma prior mean \code{rate / (shape - 1)} from the prior fields
#' \code{shape} and \code{rate} only.
#' @noRd
.two_block_tau2_ref_from_pfamily <- function(pf) {
  pl <- pf$prior_list
  if (identical(pf$pfamily, "dNormal")) {
    return(as.numeric(pl$dispersion))
  }
  if (identical(pf$pfamily, "dIndependent_Normal_Gamma")) {
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (is.finite(shape) && shape > 1 && is.finite(rate) && rate > 0) {
      return(rate / (shape - 1))
    }
    stop(
      "pfamily dIndependent_Normal_Gamma requires shape > 1 with positive rate ",
      "so the prior-mean plug-in tau^2 is defined.",
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
