#' Build measurement-prior list from matrix-level two-block inputs
#'
#' Converts Block~1 \code{prior_list} and Block~2 \code{pfamily_list} into the
#' \code{measurement_prior_list} structure required by
#' \code{\link{lmerb_posterior_mean}} and \code{\link{glmerb_posterior_mode}}.
#'
#' @param prior_list_block1 Block~1 prior with \code{P} or \code{Sigma} and,
#'   for Gaussian models, \code{dispersion}.
#' @param pfamily_list Named list of \code{pfamily} objects (validated).
#' @param re_names Random-effect component names.
#' @param x_hyper Named list of hyper design matrices (column names for
#'   Block~2 parameters).
#' @param family Response \code{family} object.
#' @return List with \code{dispersion_ranef}, \code{Sigma_ranef}, and
#'   \code{prior_list}.
#' @noRd
.two_block_measurement_prior_list <- function(
    prior_list_block1,
    pfamily_list,
    re_names,
    x_hyper,
    family
) {
  family <- .two_block_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")

  if (!is.null(prior_list_block1$Sigma)) {
    Sigma_ranef <- as.matrix(prior_list_block1$Sigma)
  } else if (!is.null(prior_list_block1$P)) {
    Sigma_ranef <- solve(as.matrix(prior_list_block1$P))
  } else {
    stop(
      "prior_list must contain 'P' or 'Sigma' for ICM.",
      call. = FALSE
    )
  }

  dispersion_ranef <- if (is_gaussian) {
    d <- prior_list_block1$dispersion
    if (is.null(d)) {
      stop(
        "prior_list must contain 'dispersion' for gaussian() ICM.",
        call. = FALSE
      )
    }
    as.numeric(d)
  } else {
    NULL
  }

  prior_list <- stats::setNames(
    vector("list", length(re_names)),
    re_names
  )

  for (k in re_names) {
    pf <- pfamily_list[[k]]
    pl <- pf$prior_list
    par_names <- colnames(x_hyper[[k]])
    if (is.null(par_names) || !length(par_names)) {
      stop(
        "x_hyper[[\"", k, "\"]] must have colnames for ICM.",
        call. = FALSE
      )
    }

    mu_k <- as.numeric(pl$mu)
    if (length(mu_k) != length(par_names)) {
      stop(
        "pfamily_list[[\"", k, "\"]]$prior_list$mu length (", length(mu_k),
        ") must match ncol(x_hyper[[\"", k, "\"]]) (", length(par_names), ").",
        call. = FALSE
      )
    }
    mu_nms <- rownames(pl$mu)
    if (!is.null(mu_nms) && all(nzchar(mu_nms))) {
      if (!setequal(mu_nms, par_names)) {
        stop(
          "Parameter names of pfamily_list[[\"", k,
          "\"]] do not match x_hyper[[\"", k, "\"]] columns.",
          call. = FALSE
        )
      }
      ord <- match(par_names, mu_nms)
      mu_k <- mu_k[ord]
      Sigma_k <- as.matrix(pl$Sigma)[ord, ord, drop = FALSE]
    } else {
      Sigma_k <- as.matrix(pl$Sigma)
    }
    names(mu_k) <- par_names
    dimnames(Sigma_k) <- list(par_names, par_names)

    d_k <- if (identical(pf$pfamily, "dNormal")) {
      pl$dispersion
    } else {
      pl$disp_lower
    }
    if (is.null(d_k) || !is.numeric(d_k) || length(d_k) != 1L ||
        !is.finite(d_k) || d_k <= 0) {
      stop(
        "pfamily_list[[\"", k, "\"]] must supply a positive dispersion ",
        "(dNormal) or disp_lower (dIndependent_Normal_Gamma) for ICM.",
        call. = FALSE
      )
    }

    prior_list[[k]] <- list(
      mu_fixef         = mu_k,
      Sigma_fixef      = Sigma_k,
      dispersion_fixef = as.numeric(d_k)
    )
  }

  list(
    dispersion_ranef = dispersion_ranef,
    Sigma_ranef      = Sigma_ranef,
    prior_list       = prior_list
  )
}

#' Iterated conditional means/modes at matrix-level two-block inputs
#'
#' @param design Design list for ICM (\code{y}, \code{Z}, \code{groups},
#'   \code{X_hyper}, \code{re_coef_names}).
#' @param prior_list Block~1 prior list.
#' @param pfamily_list Validated Block~2 \code{pfamily_list}.
#' @param re_names Random-effect component names.
#' @param family Response family (default \code{gaussian()}).
#' @param tol,maxit ICM convergence controls.
#' @return List with \code{start}, \code{b_start}, and \code{icm}.
#' @noRd
.two_block_icm_at_start <- function(
    design,
    prior_list,
    pfamily_list,
    re_names,
    family = gaussian(),
    tol   = 1e-10,
    maxit = 200L
) {
  family <- .two_block_normalize_family(family)
  mpl <- .two_block_measurement_prior_list(
    prior_list_block1 = prior_list,
    pfamily_list      = pfamily_list,
    re_names          = re_names,
    x_hyper           = design$X_hyper,
    family            = family
  )

  pm <- if (identical(family$family, "gaussian")) {
    lmerb_posterior_mean(
      design,
      mpl,
      tol   = tol,
      maxit = maxit
    )
  } else {
    glmerb_posterior_mode(
      design,
      family,
      mpl,
      tol   = tol,
      maxit = maxit
    )
  }

  list(
    start   = pm$fixef,
    b_start = pm$b_mean,
    icm     = list(
      converged  = pm$converged,
      iterations = pm$iterations,
      delta      = pm$delta
    )
  )
}

#' Validate gap tolerance for pilot chain count derivation
#' @noRd
.two_block_validate_gap_tol <- function(gap_tol) {
  if (is.null(gap_tol)) {
    return(invisible(NULL))
  }
  if (!is.numeric(gap_tol) || length(gap_tol) != 1L ||
      !is.finite(gap_tol) || gap_tol <= 0 || gap_tol >= 1) {
    stop("'gap_tol' must be NULL or a single value in (0, 1).", call. = FALSE)
  }
  invisible(as.numeric(gap_tol))
}

#' Resolve pilot chain count from family, explicit \code{n_pilot}, and \code{gap_tol}
#'
#' Gaussian models never run a pilot stage.  For non-Gaussian families,
#' an explicit non-\code{NULL} \code{n_pilot} (including \code{0L}) is used
#' as-is; when \code{n_pilot} is \code{NULL}, \code{gap_tol} derives the count
#' (or \code{0L} when \code{gap_tol} is \code{NULL}).
#' @noRd
.two_block_resolve_n_pilot <- function(family, n_pilot, gap_tol) {
  family <- .two_block_normalize_family(family)
  if (identical(family$family, "gaussian")) {
    return(0L)
  }
  if (!is.null(n_pilot)) {
    n_pilot <- as.integer(n_pilot[1L])
    if (n_pilot < 0L) {
      stop("'n_pilot' must be non-negative.", call. = FALSE)
    }
    return(n_pilot)
  }
  gap_tol <- .two_block_validate_gap_tol(gap_tol)
  if (is.null(gap_tol)) {
    return(0L)
  }
  as.integer(ceiling((stats::qnorm(0.975) / gap_tol)^2))
}
