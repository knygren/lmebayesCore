#' Build per-group dGamma priors from a Prior_Setup_lmebayes object
#'
#' Converts the per-group Block~1 measurement-dispersion calibration stored
#' in a \code{\link{Prior_Setup_lmebayes}} object into a named list of
#' \code{\link[glmbayesCore]{dGamma}} \code{pfamily} objects, one per group level.
#'
#' Prior density (\code{shape_ING}, \code{rate}) and truncation bounds
#' (\code{disp_lower}, \code{disp_upper}) both come from
#' \code{object$ing_prior_measurement_group}, calibrated once in
#' \code{\link{Prior_Setup_lmebayes}()} via
#' \code{\link[glmbayesCore]{compute_gaussian_prior}()} with the Part VI
#' extension of \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} (a
#' model-derived \code{Omega_j} folded into \code{Sigma_j}, so \code{rate}/
#' \code{sigma2_hat} integrate out both the random effects \code{b_j} and
#' the fixed effects \code{gamma}). \code{disp_lower}/\code{disp_upper} are
#' literal quantiles of that same \code{Gamma(shape_ING, rate)} marginal
#' (the \eqn{(1-\mathrm{max\_disp\_perc})}/\code{max_disp_perc} quantiles),
#' i.e. a truncated version of the actual sampling prior rather than a
#' separately-constructed window. See
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} for the full derivation.
#'
#' @param object An object of class \code{"lmebayes_prior_setup"} as returned
#'   by \code{\link{Prior_Setup_lmebayes}} (Gaussian models only).
#' @param max_disp_perc Scalar in \eqn{(0.5, 1)}; defaults to
#'   \code{object$max_disp_perc} (the value the stored \code{disp_lower}/
#'   \code{disp_upper} bounds were calibrated with). Supplying a different
#'   value recomputes the bounds as fresh quantiles of the same
#'   \code{Gamma(shape_ING, rate)} for every group.
#' @param ... Currently ignored.
#'
#' @return A named list of \code{"pfamily"} objects keyed by group levels,
#'   suitable for \code{lmerb(..., dispersion_ranef = dGamma_list(ps))}, with
#'   attributes \code{"window_diagnostics"} (data frame, one row per group,
#'   with \code{sigma2_hat}, \code{disp_lower}, \code{disp_upper}),
#'   \code{"dispersion_fit"} (\code{object$dispersion_fit}, the \code{glmmTMB}
#'   reference fit used for calibration -- \code{lmerb()}/\code{glmerb()}
#'   reuse it as their own \code{dispersion_fit} instead of re-fitting
#'   \code{glmmTMB} when \code{dispersion_ranef} carries this attribute), and
#'   \code{"calibration_source"} (\code{object$calibration_source}).
#'
#' @seealso \code{\link{Prior_Setup_lmebayes}}, \code{\link{dGamma_list}},
#'   \code{\link[glmbayesCore]{dGamma}}
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
#'     pwt_measurement = 0.01,
#'     dispformula = ~school_id
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
    ...
) {
  if (!identical(object$family$family, "gaussian")) {
    stop(
      "dGamma_list() for lmebayes_prior_setup requires family = gaussian().",
      call. = FALSE
    )
  }

  ing_grp <- object$ing_prior_measurement_group
  if (is.null(ing_grp)) {
    grp_nm <- object$design$group_name
    stop(
      "object has no ing_prior_measurement_group; call Prior_Setup_lmebayes(",
      "..., dispformula = ~", if (!is.null(grp_nm)) grp_nm else "<group_name>",
      ") on a Gaussian model to calibrate per-group measurement-dispersion ",
      "priors (dispformula = ~1, the default, skips this calibration).",
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
    group_levels <- levels(object$design$group)
    names(ing_grp) <- group_levels
  }

  ## The bounds stored on ing_grp were computed at object$max_disp_perc; if
  ## the caller asks for a different max_disp_perc here, recompute them as
  ## fresh quantiles of the same Gamma(shape_ING, rate) rather than reusing
  ## the stored ones.
  recompute <- !isTRUE(all.equal(max_disp_perc, object$max_disp_perc))

  diag_rows <- vector("list", length(group_levels))
  out <- stats::setNames(
    lapply(seq_along(group_levels), function(i) {
      lev <- group_levels[[i]]
      g <- ing_grp[[lev]]

      if (recompute || is.null(g$disp_lower) || is.null(g$disp_upper)) {
        win <- .lmebayes_ing_prior_quantile_window(
          g$shape_ING, g$rate, max_disp_perc
        )
        disp_lower <- win$disp_lower
        disp_upper <- win$disp_upper
      } else {
        disp_lower <- g$disp_lower
        disp_upper <- g$disp_upper
      }

      diag_rows[[i]] <<- data.frame(
        group      = lev,
        n_j        = g$n_j,
        sigma2_hat = unname(g$sigma2_hat),
        shape_ING  = unname(g$shape_ING),
        rate       = unname(g$rate),
        disp_lower = disp_lower,
        disp_upper = disp_upper,
        stringsAsFactors = FALSE
      )

      glmbayesCore::dGamma(
        shape          = g$shape_ING,
        rate           = g$rate,
        beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
        Inv_Dispersion = TRUE,
        max_disp_perc  = max_disp_perc,
        disp_lower     = disp_lower,
        disp_upper     = disp_upper
      )
    }),
    group_levels
  )

  window_diagnostics <- do.call(rbind, diag_rows)
  rownames(window_diagnostics) <- NULL
  attr(out, "window_diagnostics")  <- window_diagnostics
  ## Carried forward so lmerb()/glmerb() can reuse this glmmTMB reference fit
  ## as their diagnostic-only dispersion_fit instead of re-fitting it; see
  ## lmebayes:::.lmebayes_fit_glmmtmb_dispersion().
  attr(out, "dispersion_fit")      <- object$dispersion_fit
  attr(out, "calibration_source")  <- object$calibration_source

  out
}
