#' Build per-group dGamma priors from a Prior_Setup_lmebayes object
#'
#' Converts the per-group Block~1 measurement-dispersion calibration stored
#' in a \code{\link{Prior_Setup_lmebayes}} object into a named list of
#' \code{\link[glmbayesCore]{dGamma}} \code{pfamily} objects, one per group level.
#'
#' Prior density (\code{shape_ING}, \code{rate}) both come from
#' \code{object$group.ing_prior} (the per-group shape, when
#' \code{dispformula} requested per-group dispersion), calibrated once in
#' \code{\link{Prior_Setup_lmebayes}()} via
#' \code{\link[glmbayesCore]{compute_gaussian_prior}()} with the Part VI
#' extension of \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} (a
#' model-derived \code{Omega_j} folded into \code{Sigma_j}, so \code{rate}/
#' \code{sigma2_hat} integrate out both the random effects \code{b_j} and
#' the fixed effects \code{gamma}). Truncation bounds (\code{disp_lower},
#' \code{disp_upper}) are the \eqn{(1-\mathrm{max\_disp\_perc})}/
#' \code{max_disp_perc} quantiles of \code{Gamma(shape_ING + n_j/2,
#' rate_post)}, \code{rate_post} being \code{sigma2_hat} mean-matched at
#' that inflated shape -- i.e. the window tracks the \emph{posterior}
#' spread the sampler's own envelope machinery actually draws
#' \eqn{\sigma^2_j} from each sweep (\code{EnvelopeDispersionBuild.cpp}'s
#' own \code{shape2 = Shape + n_w/2} fallback), not the prior
#' \code{Gamma(shape_ING, rate)} alone. \code{shape_ING}/\code{rate}
#' themselves -- what is actually fed to the sampler as the Gamma prior --
#' are unaffected by this widening. See
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} for the full derivation.
#'
#' @param object An object of class \code{"lmebayes_prior_setup"} as returned
#'   by \code{\link{Prior_Setup_lmebayes}} (Gaussian models only).
#' @param max_disp_perc_measurement \code{NULL} (default), a scalar in
#'   \eqn{(0.5, 1)}, or a named/positional length-\eqn{J} vector (one value
#'   per group level). \code{NULL} reuses each group's own stored
#'   \code{disp_lower}/\code{disp_upper} (the value
#'   \code{object$group.ing_prior[[lev]]$max_disp_perc} was
#'   calibrated with in \code{\link{Prior_Setup_lmebayes}}). Supplying a
#'   scalar or vector recomputes the bounds as fresh quantiles of the same
#'   posterior-shape \code{Gamma(shape_ING + n_j/2, rate_post)} (see above),
#'   per group, for every group whose resolved value differs from its own
#'   stored one.
#' @param ... Currently ignored.
#'
#' @return A named list of \code{"pfamily"} objects keyed by group levels,
#'   suitable for \code{lmerb(..., dispersion_ranef = dGamma_list(ps))}
#'   \emph{or} passed directly as \code{prior_list} to
#'   \code{\link{rLMMindepNormalGamma_reg_known_vcov}}/
#'   \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}} (both routes
#'   accept a named list of \code{dGamma()} pfamilies directly; \code{mu}/
#'   \code{Sigma} are derived internally from \code{pfamily_list} and are
#'   never read from this object), with attributes
#'   \code{"window_diagnostics"} (data frame, one row per group, with
#'   \code{sigma2_hat}, \code{disp_lower}, \code{disp_upper}),
#'   \code{"measurement_prior_group"} (list of four named-by-group-level
#'   numeric vectors -- \code{shape_group}, \code{rate_group},
#'   \code{disp_lower_group}, \code{disp_upper_group} -- the same
#'   per-group quantities as \code{"window_diagnostics"}'s \code{shape_ING}/
#'   \code{rate}/\code{disp_lower}/\code{disp_upper} columns, reshaped as
#'   vectors; a convenience view, not needed for the primary use above),
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
#'     group.dispersion.pwt = 0.01,
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
    max_disp_perc_measurement = NULL,
    ...
) {
  if (!identical(object$family$family, "gaussian")) {
    stop(
      "dGamma_list() for lmebayes_prior_setup requires family = gaussian().",
      call. = FALSE
    )
  }

  ing_grp <- object$group.ing_prior
  if (is.null(ing_grp) || !.lmebayes_ing_prior_is_grouped(ing_grp)) {
    grp_nm <- object$design$group_name
    stop(
      "object's group.ing_prior is not per-group; call Prior_Setup_lmebayes(",
      "..., dispformula = ~", if (!is.null(grp_nm)) grp_nm else "<group_name>",
      ") on a Gaussian model to calibrate per-group measurement-dispersion ",
      "priors (dispformula = ~1, the default, skips this calibration).",
      call. = FALSE
    )
  }

  group_levels <- names(ing_grp)
  if (is.null(group_levels)) {
    group_levels <- levels(object$design$group)
    names(ing_grp) <- group_levels
  }

  ## NULL: reuse each group's own stored max_disp_perc (no override, no
  ## recompute). Supplied: resolve to a per-group vector (scalar recycled,
  ## or a length-J vector) and recompute per group below.
  mdp_override <- if (!is.null(max_disp_perc_measurement)) {
    .lmebayes_expand_scalar_or_vector(
      max_disp_perc_measurement, group_levels, "max_disp_perc_measurement"
    )
  } else {
    NULL
  }

  diag_rows <- vector("list", length(group_levels))
  out <- stats::setNames(
    lapply(seq_along(group_levels), function(i) {
      lev <- group_levels[[i]]
      g <- ing_grp[[lev]]

      ## The bounds stored on g were computed at g$max_disp_perc; if the
      ## caller asks for a different value for this group, recompute fresh
      ## quantiles of Gamma(shape_ING + n_j/2, rate_post) -- the posterior-
      ## shape window (see Prior_Setup_lmebayes()'s group.ing_prior
      ## docs), NOT g$shape_ING/g$rate directly (those are the prior fed to
      ## the sampler, unchanged) -- rather than reusing the stored bounds.
      mdp_lev <- if (!is.null(mdp_override)) {
        unname(mdp_override[[lev]])
      } else {
        g$max_disp_perc
      }
      recompute <- !isTRUE(all.equal(mdp_lev, g$max_disp_perc))

      if (recompute || is.null(g$disp_lower) || is.null(g$disp_upper)) {
        shape_post_lev <- g$shape_ING + g$n_j / 2
        rate_post_lev  <- g$sigma2_hat * (shape_post_lev - 1)
        win <- .lmebayes_ing_prior_quantile_window(
          shape_post_lev, rate_post_lev, mdp_lev
        )
        disp_lower <- win$disp_lower
        disp_upper <- win$disp_upper
      } else {
        disp_lower <- g$disp_lower
        disp_upper <- g$disp_upper
      }

      diag_rows[[i]] <<- data.frame(
        group         = lev,
        n_j           = g$n_j,
        sigma2_hat    = unname(g$sigma2_hat),
        shape_ING     = unname(g$shape_ING),
        rate          = unname(g$rate),
        max_disp_perc = mdp_lev,
        disp_lower    = disp_lower,
        disp_upper    = disp_upper,
        stringsAsFactors = FALSE
      )

      glmbayesCore::dGamma(
        shape          = g$shape_ING,
        rate           = g$rate,
        beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
        Inv_Dispersion = TRUE,
        max_disp_perc  = mdp_lev,
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

  ## Flattened view of the SAME four per-group quantities already collected
  ## into window_diagnostics above (shape_ING/rate/disp_lower/disp_upper),
  ## reshaped as named-by-group-level vectors instead of data-frame columns
  ## -- ready to pass straight through as 'prior_list' to
  ## rLMMindepNormalGamma_reg_known_vcov()/_estimated_vcov() (which now
  ## accept this list argument-name-for-argument-name; note 'mu'/'Sigma' are
  ## deliberately NOT included here: 'mu' is never read by those routes --
  ## each group's mean is always W_j %*% gamma, recomputed every sweep -- and
  ## 'Sigma' is the Block~2-derived random-effects covariance, which those
  ## routes now always derive internally from 'pfamily_list' as
  ## solve(P) rather than accepting it here, so it can never drift out of
  ## sync with pfamily_list). All four vectors share identical names/order
  ## (= group_levels) by construction: they come from the same per-group
  ## loop above, never independently subset afterward.
  attr(out, "measurement_prior_group") <- list(
    shape_group      = stats::setNames(window_diagnostics$shape_ING,  window_diagnostics$group),
    rate_group       = stats::setNames(window_diagnostics$rate,       window_diagnostics$group),
    disp_lower_group = stats::setNames(window_diagnostics$disp_lower, window_diagnostics$group),
    disp_upper_group = stats::setNames(window_diagnostics$disp_upper, window_diagnostics$group)
  )

  out
}
