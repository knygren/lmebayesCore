#' Reference tau^2_k for ICM and replicate-chain starts (not disp_lower)
#'
#' For \code{dNormal}, returns fixed \code{dispersion}.  For ING, uses
#' internal \code{prior_list$tau2_ref} when set by the sampler, else the
#' inverse-Gamma prior mean \code{rate / (shape - 1)}.
#' @noRd
.two_block_tau2_ref_from_pfamily <- function(pf) {
  pl <- pf$prior_list
  if (identical(pf$pfamily, "dNormal")) {
    return(as.numeric(pl$dispersion))
  }
  if (identical(pf$pfamily, "dIndependent_Normal_Gamma")) {
    if (!is.null(pl$tau2_ref)) {
      return(as.numeric(pl$tau2_ref))
    }
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (is.finite(shape) && shape > 1 && is.finite(rate) && rate > 0) {
      return(rate / (shape - 1))
    }
    stop(
      "pfamily dIndependent_Normal_Gamma requires shape > 1 or an internal ",
      "reference tau^2 (set by the two-block engine).",
      call. = FALSE
    )
  }
  stop("Unsupported pfamily: ", pf$pfamily, call. = FALSE)
}

#' Named vector of reference tau^2 values from a pfamily_list
#' @noRd
.two_block_tau2_ref_vector <- function(pfamily_list, re_names) {
  stats::setNames(
    vapply(re_names, function(k) {
      .two_block_tau2_ref_from_pfamily(pfamily_list[[k]])
    }, numeric(1)),
    re_names
  )
}

#' Copy pfamily_list with updated ING \code{tau2_ref} values
#' @noRd
.two_block_pfamily_with_tau2_ref <- function(pfamily_list, re_names, tau2_ref) {
  out <- pfamily_list
  for (k in re_names) {
    pf <- out[[k]]
    if (!identical(pf$pfamily, "dIndependent_Normal_Gamma")) {
      next
    }
    pl <- pf$prior_list
    pl$tau2_ref <- as.numeric(tau2_ref[[k]])
    pf$prior_list <- pl
    out[[k]] <- pf
  }
  out
}
