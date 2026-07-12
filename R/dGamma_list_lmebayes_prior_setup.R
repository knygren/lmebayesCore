#' Build per-group dGamma priors from a Prior_Setup_lmebayes object
#'
#' Converts the per-group Block~1 measurement-dispersion calibration stored
#' in a \code{\link{Prior_Setup_lmebayes}} object into a named list of
#' \code{\link{dGamma}} \code{pfamily} objects, one per group level.
#'
#' Prior density (\code{shape_ING}, \code{rate_gamma}) comes from
#' \code{object$ing_prior_measurement_group} (calibrated in
#' \code{Prior_Setup_lmebayes()} via \code{\link{compute_gaussian_prior}} with
#' shared population \code{sd_tau}).  Truncation bounds use an approximate
#' posterior at \eqn{n_{\mathrm{combined},j} = n_{\mathrm{prior},j} + n_j}:
#' \code{disp_lower} is OLS-anchored at \code{sigma2_hat_j}; \code{disp_upper}
#' widens by the BLUP/OLS residual RSS inflation ratio when
#' \code{disp_upper_anchor = "blup"} (default).
#'
#' @param object An object of class \code{"lmebayes_prior_setup"} as returned
#'   by \code{\link{Prior_Setup_lmebayes}} (Gaussian models only).
#' @param max_disp_perc Scalar in \eqn{(0.5, 1)}; defaults to
#'   \code{object$max_disp_perc}.
#' @param disp_upper_anchor Character: \code{"blup"} (default) scales the
#'   upper-bound rate by \eqn{\mathrm{RSS}_{\mathrm{blup}}/\mathrm{RSS}_{\mathrm{ols}}}
#'   per group; \code{"symmetric"} uses the same rate as the lower bound.
#' @param ... Currently ignored.
#'
#' @return A named list of \code{"pfamily"} objects keyed by group levels,
#'   suitable for \code{lmerb(..., dispersion_ranef = dGamma_list(ps))}.
#'
#' @seealso \code{\link{Prior_Setup_lmebayes}}, \code{\link{dGamma_list}},
#'   \code{\link{dGamma}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("bayesrules", quietly = TRUE)) {
#'   data(big_word_club, package = "bayesrules")
#'   dat <- big_word_club
#'   dat$school_id <- factor(dat$school_id)
#'   dat <- subset(dat, !is.na(score_ppvt))
#'
#'   ps <- Prior_Setup_lmebayes(
#'     score_ppvt ~ private_school + (1 | school_id),
#'     data = dat,
#'     pwt_measurement = 0.01
#'   )
#'   disp_pf <- dGamma_list(ps)
#'   print(disp_pf[[1L]])
#' }
#' }
#'
#' @export
#' @method dGamma_list lmebayes_prior_setup
dGamma_list.lmebayes_prior_setup <- function(
    object,
    max_disp_perc = NULL,
    disp_upper_anchor = c("blup", "symmetric"),
    ...
) {
  disp_upper_anchor <- match.arg(disp_upper_anchor)

  if (!identical(object$family$family, "gaussian")) {
    stop(
      "dGamma_list() for lmebayes_prior_setup requires family = gaussian().",
      call. = FALSE
    )
  }

  ing_grp <- object$ing_prior_measurement_group
  if (is.null(ing_grp)) {
    stop(
      "object has no ing_prior_measurement_group; ",
      "call Prior_Setup_lmebayes() on a Gaussian model first.",
      call. = FALSE
    )
  }

  if (is.null(max_disp_perc)) {
    max_disp_perc <- object$max_disp_perc
  }
  if (is.null(max_disp_perc)) {
    max_disp_perc <- 0.99
  }
  if (!is.numeric(max_disp_perc) || length(max_disp_perc) != 1L ||
      is.na(max_disp_perc) || max_disp_perc <= 0.5 || max_disp_perc >= 1) {
    stop("'max_disp_perc' must be a scalar in (0.5, 1).", call. = FALSE)
  }

  group_levels <- names(ing_grp)
  if (is.null(group_levels)) {
    group_levels <- levels(object$design$groups)
    names(ing_grp) <- group_levels
  }

  blup_inflation <- if (disp_upper_anchor == "blup") {
    .lmebayes_group_blup_rss_inflation(
      data          = object$data,
      block_formula = object$block_formula,
      fit_ref       = object$fit_ref,
      groups        = object$design$groups,
      group_levels  = group_levels,
      group_name    = object$design$group_name
    )
  } else {
    stats::setNames(rep(1, length(group_levels)), group_levels)
  }

  stats::setNames(
    lapply(group_levels, function(lev) {
      g <- ing_grp[[lev]]
      p_re <- g$p_re
      n_combined <- g$n_combined
      sigma2_hat <- g$sigma2_hat

      shape_w <- (n_combined + 1) / 2 + p_re / 2
      rate_w  <- sigma2_hat * (n_combined + p_re - 1) / 2
      rate_u  <- if (disp_upper_anchor == "blup") {
        rate_w * blup_inflation[[lev]]
      } else {
        rate_w
      }

      if (rate_u < rate_w - sqrt(.Machine$double.eps) * max(rate_w, 1)) {
        stop(
          "Group '", lev, "': upper-bound rate (", rate_u,
          ") is below lower-bound rate (", rate_w, ").",
          call. = FALSE
        )
      }

      win <- .lmebayes_ing_prior_quantile_window_asymmetric(
        shape      = shape_w,
        rate_lower = rate_w,
        rate_upper = rate_u,
        max_disp_perc = max_disp_perc
      )

      dGamma(
        shape          = g$shape_ING,
        rate           = g$rate_gamma,
        beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
        Inv_Dispersion = TRUE,
        max_disp_perc  = max_disp_perc,
        disp_lower     = win$disp_lower,
        disp_upper     = win$disp_upper
      )
    }),
    group_levels
  )
}
