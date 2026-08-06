## Internal helpers for matrix-level mixed-model samplers (rlmerb / rglmerb).
## lmebayes calls these via lmebayesCore::: from lmerb() / glmerb() only.

#' @noRd
.lmebayes_resolve_group_dispersion <- function(
    group.dispersion,
    family,
    design = NULL,
    fn_name = "lmerb"
) {
  has_dispersion <- family$family %in%
    c("gaussian", "Gamma", "quasipoisson", "quasibinomial")

  if (!has_dispersion) {
    if (!is.null(group.dispersion)) {
      stop(
        "'group.dispersion' must be NULL for family = ", family$family,
        "() (no observation-level dispersion).",
        call. = FALSE
      )
    }
    return(list(
      mode                   = "none",
      dispersion_fix         = NULL,
      dispersion_prior_list  = NULL,
      dispersion_pfamily     = NULL
    ))
  }

  if (is.list(group.dispersion) && !inherits(group.dispersion, "pfamily")) {
    return(.lmebayes_resolve_group_dispersion_group_list(
      group.dispersion = group.dispersion,
      design           = design,
      fn_name          = fn_name
    ))
  }

  if (inherits(group.dispersion, "pfamily")) {
    if (!identical(group.dispersion$pfamily, "dGamma")) {
      stop(
        fn_name, "(): 'group.dispersion' pfamily must be dGamma(); got ",
        group.dispersion$pfamily, ". RE priors belong in 'pfamily_list'.",
        call. = FALSE
      )
    }
    pl <- group.dispersion$prior_list
    if (!isTRUE(pl$Inv_Dispersion)) {
      stop(
        fn_name, "(): dGamma() observation-dispersion prior requires ",
        "Inv_Dispersion = TRUE.",
        call. = FALSE
      )
    }
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
      stop(
        fn_name, "(): dGamma() group.dispersion prior_list requires positive ",
        "'shape' and 'rate'.",
        call. = FALSE
      )
    }
    return(list(
      mode                  = "gamma",
      dispersion_fix        = shape / rate,
      dispersion_prior_list = pl,
      dispersion_pfamily    = group.dispersion
    ))
  }

  if (is.numeric(group.dispersion) && length(group.dispersion) > 1L) {
    return(.lmebayes_resolve_group_dispersion_fixed_vector(
      group.dispersion = group.dispersion,
      design           = design,
      fn_name          = fn_name
    ))
  }

  if (is.null(group.dispersion) || !is.numeric(group.dispersion) ||
      length(group.dispersion) != 1L || !is.finite(group.dispersion) ||
      group.dispersion <= 0) {
    stop(
      "'group.dispersion' must be a single positive number, a named numeric ",
      "vector of positive per-group values, a dGamma() pfamily, or a named ",
      "list of dGamma() pfamilies (one per group) for family = ",
      family$family, "().",
      call. = FALSE
    )
  }
  list(
    mode                  = "fixed",
    dispersion_fix        = as.numeric(group.dispersion),
    dispersion_prior_list = NULL,
    dispersion_pfamily    = NULL
  )
}

#' Resolve a fixed per-group numeric vector for \code{group.dispersion}
#'
#' Fourth \code{group.dispersion} option alongside a fixed scalar, a single
#' (pooled) \code{dGamma()}, and a named list of \code{dGamma()} pfamilies: a
#' plain named numeric vector with one known, fixed dispersion value per
#' group. Unlike \code{dGamma_list(...)} (\code{"gamma_list"}), this is not a
#' prior to be sampled -- each group's \eqn{\sigma^2_j} is treated as known
#' for the duration of sampling, exactly like the pooled \code{"fixed"} mode
#' but allowed to vary by group. No \code{glmmTMB} reference fit or full-rank
#' requirement applies (there is no ING accept/reject envelope to build).
#' @noRd
.lmebayes_resolve_group_dispersion_fixed_vector <- function(
    group.dispersion,
    design,
    fn_name = "lmerb"
) {
  if (is.null(design) || is.null(design$group)) {
    stop(
      fn_name, "(): a named numeric vector for 'group.dispersion' requires ",
      "'design' with grouping information.",
      call. = FALSE
    )
  }
  group_levels <- levels(design$group)
  J <- length(group_levels)

  if (length(group.dispersion) != J) {
    stop(
      fn_name, "(): 'group.dispersion' is a vector of length ",
      length(group.dispersion), " but there are ", J, " group level(s) (",
      paste(group_levels, collapse = ", "), "). Supply exactly one fixed ",
      "dispersion value per group.",
      call. = FALSE
    )
  }
  nms <- names(group.dispersion)
  if (is.null(nms) || any(!nzchar(nms)) || !setequal(nms, group_levels)) {
    stop(
      fn_name, "(): names(group.dispersion) must match the group levels (",
      paste(group_levels, collapse = ", "), ") exactly.",
      call. = FALSE
    )
  }
  group.dispersion <- stats::setNames(
    as.numeric(group.dispersion[group_levels]),
    group_levels
  )

  if (any(!is.finite(group.dispersion)) || any(group.dispersion <= 0)) {
    stop(
      fn_name, "(): 'group.dispersion' values must all be positive and ",
      "finite.",
      call. = FALSE
    )
  }

  list(
    mode                  = "fixed_vector",
    dispersion_fix        = group.dispersion,
    dispersion_prior_list = NULL,
    dispersion_pfamily    = NULL
  )
}

#' Resolve a per-group list of \code{dGamma()} pfamilies for \code{group.dispersion}
#'
#' Third \code{group.dispersion} option alongside a fixed scalar and a single
#' (pooled) \code{dGamma()}: a named list with one \code{dGamma()} pfamily per
#' group level. Each entry keeps its own \code{shape}/\code{rate}/\code{disp_lower}/
#' \code{disp_upper} -- there is no requirement that groups share the same
#' \code{shape}/\code{rate} (\code{Prior_Setup_lmebayes()} may choose to build
#' them from a shared hyperprior, but the engine itself is agnostic).
#' Requires every group to be full column rank (rank-deficient groups are not
#' yet supported for this option).
#' @noRd
.lmebayes_resolve_group_dispersion_group_list <- function(
    group.dispersion,
    design,
    fn_name = "lmerb"
) {
  if (is.null(design) || is.null(design$group)) {
    stop(
      fn_name, "(): a list of dGamma() priors for 'group.dispersion' requires ",
      "'design' with grouping information.",
      call. = FALSE
    )
  }
  group_levels <- levels(design$group)
  J <- length(group_levels)

  if (length(group.dispersion) != J) {
    stop(
      fn_name, "(): 'group.dispersion' is a list of length ",
      length(group.dispersion), " but there are ", J, " group level(s) (",
      paste(group_levels, collapse = ", "), "). Supply exactly one dGamma() ",
      "pfamily per group.",
      call. = FALSE
    )
  }
  nms <- names(group.dispersion)
  if (is.null(nms) || any(!nzchar(nms)) || !setequal(nms, group_levels)) {
    stop(
      fn_name, "(): names(group.dispersion) must match the group levels (",
      paste(group_levels, collapse = ", "), ") exactly.",
      call. = FALSE
    )
  }
  window_diagnostics <- attr(group.dispersion, "window_diagnostics")
  group.dispersion <- group.dispersion[group_levels]

  if (is.null(design$groupef.rank) || !all(design$groupef.rank[group_levels])) {
    deficient <- if (!is.null(design$groupef.rank)) {
      group_levels[!design$groupef.rank[group_levels]]
    } else {
      group_levels
    }
    stop(
      fn_name, "(): a list of per-group dGamma() dispersion priors currently ",
      "requires every group to be full column rank; rank-deficient group(s): ",
      paste(deficient, collapse = ", "), ". Use a single dGamma() or a fixed ",
      "scalar 'group.dispersion' for models with rank-deficient groups.",
      call. = FALSE
    )
  }

  shape_group      <- stats::setNames(numeric(J), group_levels)
  rate_group       <- stats::setNames(numeric(J), group_levels)
  disp_lower_group <- stats::setNames(numeric(J), group_levels)
  disp_upper_group <- stats::setNames(numeric(J), group_levels)

  for (lev in group_levels) {
    pf <- group.dispersion[[lev]]
    if (!inherits(pf, "pfamily") || !identical(pf$pfamily, "dGamma")) {
      stop(
        fn_name, "(): group.dispersion[[\"", lev, "\"]] must be a dGamma() ",
        "pfamily.",
        call. = FALSE
      )
    }
    pl <- pf$prior_list
    if (!isTRUE(pl$Inv_Dispersion)) {
      stop(
        fn_name, "(): group.dispersion[[\"", lev, "\"]] dGamma() prior ",
        "requires Inv_Dispersion = TRUE.",
        call. = FALSE
      )
    }
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
      stop(
        fn_name, "(): group.dispersion[[\"", lev, "\"]] must have positive, ",
        "finite 'shape' and 'rate'.",
        call. = FALSE
      )
    }
    lo <- pl$disp_lower
    hi <- pl$disp_upper
    if (is.null(lo) || is.null(hi) ||
        !is.numeric(lo) || !is.numeric(hi) ||
        length(lo) != 1L || length(hi) != 1L ||
        !is.finite(lo) || !is.finite(hi) ||
        lo <= 0 || hi <= lo) {
      stop(
        fn_name, "(): group.dispersion[[\"", lev, "\"]] must supply finite ",
        "'disp_lower' and 'disp_upper' with 0 < disp_lower < disp_upper -- a ",
        "list of dGamma() priors requires explicit per-group truncation bounds.",
        call. = FALSE
      )
    }
    shape_group[[lev]]      <- shape
    rate_group[[lev]]       <- rate
    disp_lower_group[[lev]] <- as.numeric(lo)
    disp_upper_group[[lev]] <- as.numeric(hi)
  }

  list(
    mode                  = "gamma_list",
    dispersion_fix        = mean(shape_group / rate_group),
    dispersion_prior_list = list(
      shape_group      = shape_group,
      rate_group       = rate_group,
      disp_lower_group = disp_lower_group,
      disp_upper_group = disp_upper_group
    ),
    dispersion_pfamily    = group.dispersion,
    window_diagnostics    = window_diagnostics
  )
}

#' @noRd
.lmebayes_validate_group_dispersion_arg <- function(
    group.dispersion,
    family,
    fn_name = "lmerb"
) {
  resolved <- .lmebayes_resolve_group_dispersion(
    group.dispersion = group.dispersion,
    family           = family,
    design           = NULL,
    fn_name          = fn_name
  )
  resolved$dispersion_fix
}

#' Normalize a Block~2 \code{pfamily_list} (+ \code{group.dispersion}) into a sampler \code{prior}
#'
#' @description
#' Shared front-door normalizer used by \code{lmerb()}/\code{glmerb()} (via
#' \code{lmebayesCore::priors_from_pfamily_list()} in \strong{lmebayes}) to
#' turn the two independent user-supplied prior specs -- the Block~2
#' \code{pfamily_list} (one prior per random-effect coefficient) and the
#' Block~1 \code{group.dispersion} (observation-level dispersion) -- into the
#' single flat \code{prior} object consumed by \code{\link{matrix_args_lmm}},
#' by \code{rlmerb()}/\code{rglmerb()}'s ICM/reporting code, and by the
#' routed \code{rLMM_reg}/\code{rGLMM_reg} exports.
#'
#' @details
#' Concretely, this function:
#' \enumerate{
#'   \item Resolves \code{group.dispersion} via
#'     \code{.lmebayes_resolve_group_dispersion()} (dispatches on its type/
#'     length to one of the \code{"none"}/\code{"fixed"}/\code{"fixed_vector"}/
#'     \code{"gamma"}/\code{"gamma_list"} Block~1 dispersion modes described
#'     under \code{dispersion_mode} below).
#'   \item Validates that \code{pfamily_list} has exactly one entry per
#'     random-effect coefficient in \code{design$groupef.names} (by name, not
#'     position) and reorders it to that canonical order.
#'   \item For each random-effect component, checks that its \code{pfamily}
#'     is \code{dNormal} or \code{dIndependent_Normal_Gamma} (any other
#'     \code{pfamily} is rejected), checks that \code{prior_list$mu}/
#'     \code{prior_list$Sigma} conform to the corresponding per-group
#'     hyper-design \code{design$W[[k]]}, and reorders/relabels them to
#'     that hyper-design's column order (so the \code{pfamily} objects
#'     returned in \code{pfamily_list} are safe to pass straight to the
#'     matrix-level samplers).
#'   \item Computes each component's Block~2 variance-component plug-in
#'     \eqn{\tau^2_k} via \code{.two_block_tau2_plug_in_from_pfamily()} (the
#'     \code{dNormal()} dispersion for \code{dNormal} components, or an
#'     ICM-style plug-in for \code{dIndependent_Normal_Gamma} components) and
#'     assembles \code{group.Sigma} and \code{pop.prior_list} from those.
#' }
#' It currently combines all of the above (dispersion resolution,
#' \code{pfamily_list} validation/reordering, and deriving
#' \code{group.Sigma}/\code{pop.prior_list}/\code{ptypes}) into one function and
#' is a refactor candidate; its argument list and return shape are not yet
#' considered stable.
#'
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects, one per
#'   random-effect coefficient in \code{design$groupef.names}.
#' @param group.dispersion Block~1 dispersion spec, as accepted by
#'   \code{.lmebayes_resolve_group_dispersion()} (\code{NULL}, a single
#'   scalar, a named/unnamed numeric vector, a \code{dGamma()} pfamily, or a
#'   named list of \code{dGamma()} pfamily objects).
#' @param design A \code{model_setup} object.
#' @param family A \code{\link[stats]{family}} object.
#' @param fn_name Character scalar used to prefix error messages
#'   (e.g. \code{"lmerb"} or \code{"glmerb"}).
#' @return A list (the \code{prior} object) with elements:
#'   \describe{
#'     \item{\code{pfamily_list}}{The input \code{pfamily_list}, reordered to
#'       \code{design$groupef.names} and with each component's
#'       \code{prior_list$mu}/\code{prior_list$Sigma} realigned to the
#'       column order of the corresponding \code{design$W[[k]]}. Safe
#'       to pass straight through to the matrix-level samplers.}
#'     \item{\code{group.dispersion}}{The \emph{resolved} Block~1 dispersion
#'       value (i.e. \code{disp_res$dispersion_fix}, not the raw input):
#'       \code{NULL} when \code{family} has no dispersion parameter
#'       (\code{dispersion_mode == "none"}); a single positive scalar for
#'       \code{"fixed"} (the value supplied) or for \code{"gamma"}/
#'       \code{"gamma_list"} (a plug-in point estimate -- \code{shape/rate},
#'       or the across-group mean of \code{shape_group/rate_group} --
#'       \emph{not} the prior itself, see \code{dispersion_prior_list}); or a
#'       named numeric vector, one entry per group level, for
#'       \code{"fixed_vector"}.}
#'     \item{\code{dispersion_mode}}{Character scalar: one of \code{"none"}
#'       (no observation-level dispersion for this \code{family}),
#'       \code{"fixed"} (single known scalar \eqn{\sigma^2}),
#'       \code{"fixed_vector"} (one known, fixed \eqn{\sigma^2_j} per group,
#'       from a named numeric vector), \code{"gamma"} (a single pooled
#'       \code{dGamma()} prior on the observation precision/dispersion, to be
#'       estimated -- ING), or \code{"gamma_list"} (one \code{dGamma()} prior
#'       per group, to be estimated -- ING).}
#'     \item{\code{dispersion_pfamily}}{\code{NULL} for
#'       \code{"none"}/\code{"fixed"}/\code{"fixed_vector"}. For
#'       \code{"gamma"}, the original \code{dGamma()} pfamily object passed
#'       as \code{group.dispersion}. For \code{"gamma_list"}, the original
#'       named list of per-group \code{dGamma()} pfamily objects.}
#'     \item{\code{dispersion_prior_list}}{\code{NULL} for
#'       \code{"none"}/\code{"fixed"}/\code{"fixed_vector"} (there is no
#'       Block~1 prior to carry -- the dispersion is a known constant). For
#'       \code{"gamma"}, the pooled \code{dGamma()} prior's
#'       \code{prior_list} (\code{shape}, \code{rate}, \code{disp_lower},
#'       \code{disp_upper}, \code{Inv_Dispersion}, \ldots). For
#'       \code{"gamma_list"}, a list with named numeric vectors
#'       \code{shape_group}, \code{rate_group}, \code{disp_lower_group}, and
#'       \code{disp_upper_group} (one value per group level).}
#'     \item{\code{window_diagnostics}}{Usually \code{NULL}. For
#'       \code{"gamma_list"}, the \code{"window_diagnostics"} attribute
#'       carried on the \code{group.dispersion} list (if any), describing how
#'       each group's \code{disp_lower}/\code{disp_upper} truncation window
#'       was calibrated (e.g. by \code{Prior_Setup_lmebayes()}).}
#'     \item{\code{group.Sigma}}{A \code{p_re x p_re} diagonal matrix (row/
#'       column names \code{design$groupef.names}) holding each component's
#'       Block~2 variance-component plug-in \eqn{\tau^2_k} on the diagonal --
#'       the random-effect prior covariance implied by \code{pfamily_list}.
#'       Its inverse is the Block~2 random-effect prior precision the
#'       matrix-level samplers derive internally from \code{pfamily_list}
#'       (there is no separate \code{P} argument).}
#'     \item{\code{pop.prior_list}}{A named list (one entry per
#'       \code{design$groupef.names}), each a list with \code{mu}
#'       (numeric vector, the Block~2 hyperparameter prior mean),
#'       \code{Sigma} (matrix, the Block~2 hyperparameter prior
#'       covariance), and \code{dispersion} (scalar, the same
#'       \eqn{\tau^2_k} plug-in stored on the diagonal of
#'       \code{group.Sigma}). A restructured, per-component echo of
#'       \code{pfamily_list}'s contents keyed for direct use elsewhere (e.g.
#'       ICM reporting).}
#'     \item{\code{ptypes}}{A named character vector (names
#'       \code{design$groupef.names}), each entry either \code{"dNormal"} or
#'       \code{"dIndependent_Normal_Gamma"} -- the \code{pfamily} type of the
#'       corresponding Block~2 component.}
#'     \item{\code{any_non_normal}}{Logical scalar: \code{TRUE} if any
#'       \code{ptypes} entry is not \code{"dNormal"} (i.e. at least one
#'       Block~2 component is \code{dIndependent_Normal_Gamma} and must be
#'       estimated rather than known). Drives the \code{"known"} vs.
#'       \code{"estimated"} route choice in \code{\link{matrix_args_lmm}} /
#'       \code{.lmebayes_reg_route_key()}.}
#'   }
#' @keywords internal
#' @export
priors_from_pfamily_list <- function(pfamily_list,
                                      group.dispersion,
                                      design,
                                      family,
                                      fn_name = "lmerb") {

  re_names <- design$groupef.names
  p_re     <- length(re_names)

  ## --- group.dispersion (Block 1 measurement dispersion) -------------------
  disp_res <- .lmebayes_resolve_group_dispersion(
    group.dispersion = group.dispersion,
    family           = family,
    design           = design,
    fn_name          = fn_name
  )
  group.dispersion <- disp_res$dispersion_fix

  ## --- pfamily_list ---------------------------------------------------------
  if (!is.list(pfamily_list) || length(pfamily_list) != p_re) {
    stop(
      "'pfamily_list' must be a list with one pfamily per random-effect ",
      "component (", p_re, " expected: ", paste(re_names, collapse = ", "),
      "). Build it with pfamily_list(Prior_Setup_lmebayes(...)).",
      call. = FALSE
    )
  }
  if (is.null(names(pfamily_list)) || !setequal(names(pfamily_list), re_names)) {
    stop(
      "Names of 'pfamily_list' must match the random-effect coefficient ",
      "names: ", paste(re_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  pfamily_list <- pfamily_list[re_names]

  prior_list <- stats::setNames(vector("list", p_re), re_names)
  tau2   <- stats::setNames(numeric(p_re), re_names)
  ptypes <- stats::setNames(character(p_re), re_names)

  for (k in re_names) {
    pf <- pfamily_list[[k]]
    if (!inherits(pf, "pfamily")) {
      stop("pfamily_list[[\"", k, "\"]] must be a pfamily object.",
           call. = FALSE)
    }
    if (!pf$pfamily %in% c("dNormal", "dIndependent_Normal_Gamma")) {
      stop(
        fn_name, "() supports only dNormal and dIndependent_Normal_Gamma ",
        "pfamilies in 'pfamily_list'; component \"", k, "\" is ",
        pf$pfamily, ".",
        call. = FALSE
      )
    }
    ptypes[[k]] <- pf$pfamily

    par_names <- colnames(design$W[[k]])
    q_k <- length(par_names)

    mu_k <- as.numeric(pf$prior_list$mu)
    if (length(mu_k) != q_k) {
      stop(
        "pfamily_list[[\"", k, "\"]]$prior_list$mu has length ",
        length(mu_k), " but the hyper design has ", q_k, " column(s): ",
        paste(par_names, collapse = ", "), ".",
        call. = FALSE
      )
    }
    mu_nms <- rownames(pf$prior_list$mu)
    if (!is.null(mu_nms) && all(nzchar(mu_nms))) {
      if (!setequal(mu_nms, par_names)) {
        stop(
          "Parameter names of pfamily_list[[\"", k, "\"]] (",
          paste(mu_nms, collapse = ", "), ") do not match the hyper design ",
          "columns (", paste(par_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      ord <- match(par_names, mu_nms)
      mu_k <- mu_k[ord]
      Sigma_k <- as.matrix(pf$prior_list$Sigma)[ord, ord, drop = FALSE]
    } else {
      Sigma_k <- as.matrix(pf$prior_list$Sigma)
    }
    names(mu_k) <- par_names
    dimnames(Sigma_k) <- list(par_names, par_names)

    ## Keep the pfamily object itself aligned with the hyper-design column
    ## order: it is passed straight to the v2 sampler as the Block 2 source
    ## of truth, so its mu/Sigma must match x_hyper[[k]].
    pfamily_list[[k]]$prior_list$mu <-
      matrix(mu_k, ncol = 1L, dimnames = list(par_names, NULL))
    pfamily_list[[k]]$prior_list$Sigma <- Sigma_k

    if (identical(pf$pfamily, "dNormal")) {
      d_k <- pf$prior_list$dispersion
      if (isTRUE(pf$prior_list$ddef)) {
        warning(
          fn_name, ": pfamily_list[[\"", k, "\"]] uses the default ",
          "dispersion = 1 (none was supplied to dNormal()); the Block 1 ",
          "random-effect variance tau^2 for \"", k, "\" is therefore 1.",
          call. = FALSE
        )
      }
    } else {
      ## ING: disp_lower/disp_upper fix the truncation window and lambda*
      ## calibration only; ICM plug-in tau^2 comes from the pfamily spec.
      d_k <- pf$prior_list$disp_lower
      if (is.null(d_k) || !is.numeric(d_k) || length(d_k) != 1L ||
          !is.finite(d_k) || d_k <= 0) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a positive scalar ",
          "'disp_lower' (lower dispersion truncation) for lambda* calibration.",
          call. = FALSE
        )
      }
      u_k <- pf$prior_list$disp_upper
      if (is.null(u_k) || !is.numeric(u_k) || length(u_k) != 1L ||
          !is.finite(u_k) || u_k <= as.numeric(d_k)) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a finite scalar ",
          "'disp_upper' > 'disp_lower' (upper dispersion truncation), so ",
          "the tau^2 truncation window is fixed across Gibbs sweeps. ",
          "pfamily_list(Prior_Setup_lmebayes(...)) sets both bounds to the ",
          "0.01/0.99 prior dispersion quantiles by default.",
          call. = FALSE
        )
      }
    }

    tau2_k <- .two_block_tau2_plug_in_from_pfamily(pf)
    tau2[[k]] <- tau2_k
    prior_list[[k]] <- list(
      mu         = mu_k,
      Sigma      = Sigma_k,
      dispersion = tau2_k
    )
  }

  group.Sigma <- diag(unname(tau2), nrow = p_re, ncol = p_re)
  dimnames(group.Sigma) <- list(re_names, re_names)

  list(
    pfamily_list          = pfamily_list,
    group.dispersion      = group.dispersion,
    dispersion_mode       = disp_res$mode,
    dispersion_pfamily    = disp_res$dispersion_pfamily,
    dispersion_prior_list = disp_res$dispersion_prior_list,
    window_diagnostics    = disp_res$window_diagnostics,
    group.Sigma           = group.Sigma,
    pop.prior_list        = prior_list,
    ptypes         = ptypes,
    any_non_normal = any(ptypes != "dNormal")
  )
}

#' Expand a scalar-or-length-\eqn{N} argument into a named length-\eqn{N} vector
#'
#' Shared "scalar or one value per named unit" resolver used for arguments
#' that can be supplied either as a single number (recycled to every unit in
#' \code{names_ref}) or as a length-\eqn{N} numeric vector (named, matching
#' \code{names_ref} in any order, or positional). Generalizes the \code{expand()}
#' closure in \code{.lmebayes_resolve_disp_prior()} and the \code{check_w()} /
#' recycling logic in \code{.lmebayes_resolve_measurement_disp_prior_group()}
#' so both Block~1 (per-group) and Block~2 (per-RE-component) "scalar or
#' vector" arguments -- and \code{dGamma_list()}'s override -- share one
#' validation path.
#' @noRd
.lmebayes_expand_scalar_or_vector <- function(x, names_ref, what,
                                               range = c(0.5, 1)) {
  n <- length(names_ref)
  check_range <- function(v) {
    if (!is.numeric(v) || anyNA(v) || any(v <= range[1]) || any(v >= range[2])) {
      stop(
        sprintf(
          "'%s' must be numeric with all values in (%s, %s).",
          what, range[1], range[2]
        ),
        call. = FALSE
      )
    }
  }

  if (length(x) == 1L) {
    check_range(x)
    v <- rep(as.numeric(x), n)
  } else if (length(x) == n) {
    check_range(x)
    v <- as.numeric(x)
    nms <- names(x)
    if (!is.null(nms) && any(nzchar(nms))) {
      if (!setequal(nms, names_ref)) {
        stop(
          sprintf(
            "Names of '%s' must match: %s.",
            what, paste(names_ref, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      names(v) <- nms
      v <- v[names_ref]
    }
  } else {
    stop(
      sprintf("'%s' must have length 1 or %d.", what, n),
      call. = FALSE
    )
  }

  stats::setNames(v, names_ref)
}

#' Resolve Block~1 \eqn{\sigma^2} prior weight into observation-scale \code{n_prior}.
#'
#' \eqn{n_{\mathrm{prior}} = w/(1-w)\times n} with \eqn{w =} \code{pwt_measurement}.
#' Independent of Block~2 fixef \code{pwt} and Block~2 \eqn{\tau^2}
#' \code{pwt_dispersion}.
#' @noRd
.lmebayes_resolve_measurement_disp_prior <- function(
    pwt_measurement,
    n_prior_measurement,
    n_obs
) {
  if (!is.null(pwt_measurement) && !is.null(n_prior_measurement)) {
    stop(
      "Supply at most one of 'pwt_measurement' and 'n_prior_measurement'.",
      call. = FALSE
    )
  }
  if (!is.numeric(n_obs) || length(n_obs) != 1L || !is.finite(n_obs) ||
      n_obs <= 0) {
    stop("'n_obs' must be a positive finite scalar.", call. = FALSE)
  }

  if (!is.null(pwt_measurement)) {
    if (!is.numeric(pwt_measurement) || length(pwt_measurement) != 1L) {
      stop(
        "'pwt_measurement' for the pooled Block~1 path must be a scalar; ",
        "supply a length-J vector for per-group calibration via ",
        "dGamma_list() only.",
        call. = FALSE
      )
    }
    if (is.na(pwt_measurement) || pwt_measurement <= 0 || pwt_measurement >= 1) {
      stop(
        "'pwt_measurement' must be a scalar in (0, 1).",
        call. = FALSE
      )
    }
    w <- as.numeric(pwt_measurement)
    n_prior <- w / (1 - w) * n_obs
    src <- "user-supplied (group.dispersion.pwt)"
  } else if (!is.null(n_prior_measurement)) {
    if (!is.numeric(n_prior_measurement) || length(n_prior_measurement) != 1L ||
        is.na(n_prior_measurement) || n_prior_measurement <= 0 ||
        !is.finite(n_prior_measurement)) {
      stop(
        "'n_prior_measurement' must be a positive finite scalar.",
        call. = FALSE
      )
    }
    n_prior <- as.numeric(n_prior_measurement)
    w <- n_prior / (n_prior + n_obs)
    src <- "user-supplied (group.dispersion.nprior)"
  } else {
    w <- 0.01
    n_prior <- w / (1 - w) * n_obs
    src <- "default (group.dispersion.pwt = 0.01)"
  }

  if (n_prior > n_obs) {
    stop(
      "Measurement dispersion prior requires n_prior <= n (equivalently ",
      "pwt_measurement <= 0.5); got n_prior = ", signif(n_prior, 4),
      ", n = ", n_obs, ".",
      call. = FALSE
    )
  }

  list(
    pwt_measurement     = w,
    n_prior_measurement = n_prior,
    source              = src
  )
}

#' Resolve per-group Block~1 \eqn{\sigma^2} prior weights into \code{n_prior_j}
#'
#' \eqn{n_{\mathrm{prior},j} = w_j/(1-w_j)\times n_j} for each group level.
#' When \code{pwt_measurement} is a scalar, the same weight applies to every
#' group.  When it is a length-\eqn{J} vector, names must match
#' \code{group_levels} if supplied.
#' @noRd
.lmebayes_resolve_measurement_disp_prior_group <- function(
    pwt_measurement,
    n_prior_measurement,
    n_j,
    group_levels
) {
  if (!is.null(pwt_measurement) && !is.null(n_prior_measurement)) {
    stop(
      "Supply at most one of 'pwt_measurement' and 'n_prior_measurement'.",
      call. = FALSE
    )
  }

  J <- length(group_levels)
  n_j <- as.integer(n_j)
  if (length(n_j) != J || anyNA(n_j) || any(n_j <= 0L)) {
    stop(
      "'n_j' must be a positive integer vector of length J (one per group level).",
      call. = FALSE
    )
  }
  names(n_j) <- group_levels

  check_w <- function(v, what) {
    if (!is.numeric(v) || anyNA(v) || any(v <= 0) || any(v >= 1)) {
      stop(sprintf("%s must be numeric with all values in (0, 1).", what),
           call. = FALSE)
    }
  }

  if (!is.null(pwt_measurement)) {
    if (length(pwt_measurement) == 1L) {
      check_w(pwt_measurement, "'pwt_measurement'")
      w <- rep(as.numeric(pwt_measurement), J)
    } else {
      if (length(pwt_measurement) != J) {
        stop(
          sprintf(
            "'pwt_measurement' vector must have length J = %d (number of group levels).",
            J
          ),
          call. = FALSE
        )
      }
      check_w(pwt_measurement, "'pwt_measurement'")
      w <- as.numeric(pwt_measurement)
      nms <- names(pwt_measurement)
      if (!is.null(nms) && any(nzchar(nms))) {
        if (!setequal(nms, group_levels)) {
          stop(
            "Names of 'pwt_measurement' must match group levels: ",
            paste(group_levels, collapse = ", "),
            call. = FALSE
          )
        }
        names(w) <- nms
        w <- w[group_levels]
      } else {
        names(w) <- group_levels
      }
    }
    n_prior <- w / (1 - w) * n_j
    src <- if (length(pwt_measurement) == 1L) {
      "user-supplied scalar (group.dispersion.pwt)"
    } else {
      "user-supplied vector (group.dispersion.pwt)"
    }
  } else {
    w <- rep(0.01, J)
    n_prior <- w / (1 - w) * n_j
    src <- if (!is.null(n_prior_measurement)) {
      "default per group (group.dispersion.pwt = 0.01; scalar group.dispersion.nprior applies to pooled path only)"
    } else {
      "default (group.dispersion.pwt = 0.01 per group)"
    }
  }

  names(w) <- group_levels
  names(n_prior) <- group_levels

  if (any(n_prior > n_j)) {
    bad <- names(n_prior)[n_prior > n_j]
    stop(
      "Per-group measurement dispersion prior requires n_prior_j <= n_j for every group; ",
      "failed for: ", paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    pwt_measurement     = w,
    n_prior_measurement = n_prior,
    source              = src
  )
}

#' Within-group Block~1 formula from random-coefficient names
#'
#' Per-group \eqn{\sigma^2} calibration (\code{\link[glmbayesCore]{Prior_Setup}} parity) fits
#' only predictors that enter the within-group likelihood---the population-mean
#' structure aligned with \code{design$groupef.names}.  Level-2 hyper covariates
#' and cross-level moderation terms in the full mixed-model formula are excluded.
#' @noRd
.lmebayes_block_formula_from_re <- function(formula, groupef.names) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (length(groupef.names) < 1L || anyNA(groupef.names)) {
    stop(
      "'groupef.names' must be a non-empty character vector.",
      call. = FALSE
    )
  }

  resp <- all.vars(formula)[1L]
  slope_terms <- setdiff(groupef.names, "(Intercept)")
  rhs <- if (length(slope_terms) == 0L) {
    "1"
  } else {
    paste(c("1", slope_terms), collapse = " + ")
  }

  stats::as.formula(paste(resp, "~", rhs))
}

#' Prior mean vector for block-formula Gaussian calibration (Prior_Setup parity)
#'
#' Matches \code{\link[glmbayesCore]{Prior_Setup}} defaults on a group subset:
#' intercept from an intercept-only \code{lm()} when
#' \code{intercept_source = "null_model"}, slopes zero when
#' \code{effects_source = "null_effects"}.
#' @noRd
.lmebayes_block_formula_prior_mu <- function(
    block_formula,
    dat_j,
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model")
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)

  X         <- stats::model.matrix(block_formula, data = dat_j)
  var_names <- colnames(X)
  mu        <- rep(0, length(var_names))
  names(mu) <- var_names

  if ("(Intercept)" %in% var_names) {
    if (intercept_source == "null_model") {
      resp <- all.vars(block_formula)[1L]
      null_fit <- stats::lm(
        stats::as.formula(paste(resp, "~ 1")),
        data = dat_j
      )
      mu["(Intercept)"] <- unname(stats::coef(null_fit)["(Intercept)"])
    } else {
      full_fit <- stats::lm(block_formula, data = dat_j)
      mu["(Intercept)"] <- unname(stats::coef(full_fit)["(Intercept)"])
    }
  }

  if (effects_source == "full_model") {
    full_fit <- stats::lm(block_formula, data = dat_j)
    for (nm in setdiff(var_names, "(Intercept)")) {
      mu[nm] <- unname(stats::coef(full_fit)[nm])
    }
  }

  matrix(mu, ncol = 1L, dimnames = list(var_names, "mu"))
}

#' Per-group Gaussian measurement-dispersion calibration (Block~1 dGamma density)
#'
#' Within-group glm inputs for per-group measurement-dispersion ING calibration.
#' @noRd
.lmebayes_ing_prior_measurement_group_glm_inputs <- function(
    lev,
    dat_j,
    block_formula,
    sd_tau,
    family = gaussian(),
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model")
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)

  mf <- stats::model.frame(block_formula, data = dat_j)
  X  <- stats::model.matrix(block_formula, data = dat_j)
  Y  <- stats::model.response(mf)
  var_names <- colnames(X)
  nvar <- ncol(X)
  n_j  <- nrow(X)
  weights <- rep(1, n_j)
  offset  <- rep(0, n_j)

  glm_full <- stats::glm.fit(
    x = X,
    y = Y,
    weights = weights,
    family = family
  )
  glm_full$weights <- weights
  class(glm_full) <- c("glm", "lm")

  V0 <- stats::vcov(glm_full)
  if (anyNA(V0)) {
    XtW <- sweep(X, 1, weights, `*`)
    Gm  <- crossprod(XtW, X)
    Ginv <- tryCatch(
      solve(Gm),
      error = function(e) {
        stop(
          "Group '", lev, "': vcov(glm) is NA and (X'WX) is singular.",
          call. = FALSE
        )
      }
    )
    res <- Y - X %*% coef(glm_full)
    rss <- sum(weights * res^2)
    if (n_j <= nvar || !is.finite(rss) || rss <= 0) {
      stop(
        "Group '", lev, "': cannot recover vcov for rank-deficient glm fit.",
        call. = FALSE
      )
    }
    d_v0 <- rss / (n_j - nvar)
    V0 <- d_v0 * Ginv
    dimnames(V0) <- list(var_names, var_names)
  }

  V0_diag <- diag(V0)
  if (any(V0_diag <= 0)) {
    stop(
      "Group '", lev, "': diagonal entries of V0 must be positive.",
      call. = FALSE
    )
  }

  sd_vec <- sd_tau[var_names]
  if (anyNA(sd_vec)) {
    stop(
      "Group '", lev, "': block_formula coefficients must align with sd_tau names.",
      call. = FALSE
    )
  }

  bhat <- coef(glm_full)
  res  <- residuals(glm_full, type = "response")
  rss  <- sum(weights * res^2)
  if (n_j <= nvar || !is.finite(rss) || rss <= 0) {
    stop(
      "Group '", lev, "': Gaussian dispersion requires n_j > p.",
      call. = FALSE
    )
  }
  dispersion_classical <- rss / (n_j - nvar)
  mu <- .lmebayes_block_formula_prior_mu(
    block_formula    = block_formula,
    dat_j            = dat_j,
    intercept_source = intercept_source,
    effects_source   = effects_source
  )

  list(
    X                    = X,
    Y                    = Y,
    weights              = weights,
    offset               = offset,
    V0                   = V0,
    bhat                 = bhat,
    dispersion_classical = dispersion_classical,
    mu_vec               = as.numeric(mu),
    var_names            = var_names,
    nvar                 = nvar,
    n_j                  = n_j,
    sd_vec               = sd_vec
  )
}

#' Pack \code{compute_gaussian_prior()} output for measurement-dispersion lists.
#' @noRd
.lmebayes_ing_prior_list_from_cal <- function(
    cal,
    n_prior_j,
    n_j,
    p_re,
    pwt_record,
    pwt_group_j
) {
  sh <- cal$shape_ING
  rt <- cal$rate_gamma
  list(
    sigma2_hat  = cal$dispersion,
    shape       = cal$shape,
    shape_ING   = sh,
    rate        = cal$rate,
    rate_gamma  = rt,
    E_sigma2    = if (is.finite(sh) && sh > 1 && is.finite(rt) && rt > 0) {
      rt / (sh - 1)
    } else {
      NA_real_
    },
    inv_E       = if (is.finite(sh) && sh > 0 && is.finite(rt) && rt > 0) {
      rt / sh
    } else {
      NA_real_
    },
    n_prior     = n_prior_j,
    n_j         = n_j,
    n_combined  = n_prior_j + n_j,
    p_re        = p_re,
    pwt         = pwt_record,
    pwt_group   = pwt_group_j
  )
}

#' \code{compute_gaussian_prior()} on within-group glm inputs and \code{Sigma}.
#' @noRd
.lmebayes_compute_ing_prior_cal_from_sigma <- function(inp, Sigma, n_prior_j) {
  Sigma_0 <- Sigma / inp$dispersion_classical
  glmbayesCore::compute_gaussian_prior(
    X           = inp$X,
    Y           = inp$Y,
    weights     = inp$weights,
    offset      = inp$offset,
    dispersion  = NULL,
    n_effective = inp$n_j,
    bhat        = inp$bhat,
    mu          = inp$mu_vec,
    Sigma_0     = Sigma_0,
    Sigma       = Sigma,
    n_prior     = n_prior_j,
    k           = 1
  )
}

#' Per-group Block~1 measurement-dispersion calibration for \code{dGamma_list()}.
#'
#' \code{sigma2_hat}, \code{shape_ING}, and \code{rate_gamma} from shared
#' population \code{sd_tau} coefficient shrinkage (\eqn{V_0} scaled by
#' per-coefficient \code{pwt_j}). Also stores \code{rate} (A12 3.3.4
#' \eqn{S_{\mathrm{marg}}}) for dev comparison only.
#'
#' Also folds in the Part VI extension of
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} -- "also integrating out
#' the prior mean \code{mu_j}" -- with a model-derived \code{Omega_j}: each
#' RE component's \code{mu_j[k]} stands in for the model's own conditional
#' prior mean \code{W_j[[k]] \%*\% gamma_k} (the Block~2 hyper-regression),
#' and \code{gamma_k}'s own calibrated uncertainty (\code{prior_list[[k]]$Sigma},
#' already computed above in \code{Prior_Setup_lmebayes()}) propagates
#' through group \code{j}'s own hyper-design row,
#' \code{Omega_j[k, k] = W_j[[k]] \%*\% Sigma_k \%*\% t(W_j[[k]])}
#' (diagonal across RE components -- each \code{gamma_k} is calibrated
#' independently). \code{Sigma_j' = Sigma_j + Omega_j} is then used in place
#' of \code{Sigma_j} for \code{compute_gaussian_prior()}, so the resulting
#' \code{rate}/\code{sigma2_hat} integrate out both \code{b_j} (random
#' effects, as before) and \code{gamma} (fixed effects, via \code{Omega_j}).
#' This is now the permanent default (not opt-in) -- see
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part VI.
#'
#' \code{disp_lower}/\code{disp_upper} are computed here too, as literal
#' quantiles of this same \code{Gamma(shape_ING, rate)} marginal (via
#' \code{\link{.lmebayes_ing_prior_quantile_window}}) -- the same construction
#' already used for the pooled \code{ing_prior_measurement} case -- rather
#' than \code{dGamma_list()}'s former, decoupled \code{n_combined}-based
#' window.
#' @noRd
.lmebayes_calibrate_ing_prior_measurement_group <- function(
    design,
    data,
    block_formula,
    sd_tau,
    pwt_group,
    n_prior_group,
    group_levels,
    prior_list,
    max_disp_perc_group,
    family = gaussian(),
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model")
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)
  p_re <- length(design$groupef.names)
  if (p_re < 1L) {
    stop(
      "Per-group measurement dispersion calibration requires at least one random coefficient.",
      call. = FALSE
    )
  }
  if (length(sd_tau) != p_re || anyNA(sd_tau) || any(sd_tau <= 0)) {
    stop(
      "'sd_tau' must be a named numeric vector of positive RE standard deviations.",
      call. = FALSE
    )
  }
  re_names <- design$groupef.names

  stats::setNames(
    lapply(group_levels, function(lev) {
      idx   <- design$group == lev
      dat_j <- data[idx, , drop = FALSE]
      n_prior_j <- unname(n_prior_group[[lev]])

      inp <- .lmebayes_ing_prior_measurement_group_glm_inputs(
        lev              = lev,
        dat_j            = dat_j,
        block_formula    = block_formula,
        sd_tau           = sd_tau,
        family           = family,
        intercept_source = intercept_source,
        effects_source   = effects_source
      )

      pwt_j <- diag(inp$V0)
      pwt_j <- pwt_j / (pwt_j + inp$sd_vec^2)
      names(pwt_j) <- inp$var_names

      if (length(pwt_j) == 1L) {
        Sigma <- ((1 - pwt_j) / pwt_j) * inp$V0
      } else {
        scale_vec <- sqrt((1 - pwt_j) / pwt_j)
        Sigma <- inp$V0 * outer(scale_vec, scale_vec)
      }

      ## Part VI: model-derived Omega_j (fixed-effect/gamma uncertainty
      ## about b_j's prior mean), diagonal across RE components.
      Omega_j <- matrix(
        0, nrow = length(inp$var_names), ncol = length(inp$var_names),
        dimnames = list(inp$var_names, inp$var_names)
      )
      for (k in re_names) {
        Wk_row        <- design$W[[k]][lev, , drop = FALSE]
        Sigma_k <- prior_list[[k]]$Sigma
        Omega_j[k, k] <- as.numeric(Wk_row %*% Sigma_k %*% t(Wk_row))
      }

      cal <- .lmebayes_compute_ing_prior_cal_from_sigma(
        inp, Sigma + Omega_j, n_prior_j
      )

      .ing_stop_if_prior_exceeds_data(
        shape       = cal$shape_ING,
        p           = inp$nvar,
        n_w         = inp$n_j,
        detail      = paste0("group '", lev, "' has n_j = ", inp$n_j),
        limit_label = "n_j",
        prefix      = "Per-group measurement dispersion: "
      )

      ## The window must bound the POSTERIOR spread the sampler's own
      ## envelope machinery actually draws sigma2_j from (EnvelopeDispersionBuild.cpp:
      ## shape2 = Shape + n_w/2), not the prior alone -- shape_ING,j (cal$shape_ING)
      ## is the PRIOR shape (n_prior,j-only) fed to the sampler as-is (unchanged
      ## below); only the window's own shape/rate are inflated by n_j/2, mean-matched
      ## at the same sigma2_hat,j. See inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md Part II
      ## (algebraically shape_post,j == shape_w,j = (n_combined,j+1)/2 + p_re/2).
      mdp_j <- unname(max_disp_perc_group[[lev]])
      shape_post_j <- cal$shape_ING + inp$n_j / 2
      rate_post_j  <- cal$dispersion * (shape_post_j - 1)
      win <- .lmebayes_ing_prior_quantile_window(
        shape_post_j, rate_post_j, mdp_j
      )

      out <- .lmebayes_ing_prior_list_from_cal(
        cal         = cal,
        n_prior_j   = n_prior_j,
        n_j         = inp$n_j,
        p_re        = p_re,
        pwt_record  = pwt_j,
        pwt_group_j = unname(pwt_group[[lev]])
      )
      out$disp_lower    <- win$disp_lower
      out$disp_upper    <- win$disp_upper
      out$max_disp_perc <- mdp_j
      out$omega_j       <- Omega_j
      out
    }),
    group_levels
  )
}

#' Print \code{rate_gamma} (A12 3.3.5, downstream) vs \code{rate} (A12 3.3.4
#' \eqn{S_{\mathrm{marg}}}) from the same per-group calibration.
#' @noRd
.lmebayes_print_ing_prior_measurement_group_compare <- function(
    existing,
    digits = 4
) {
  if (is.null(existing)) {
    return(invisible(NULL))
  }
  grp <- names(existing)
  if (is.null(grp) || length(grp) < 1L) {
    return(invisible(NULL))
  }

  pct_diff <- function(new, old) {
    if (!is.finite(old) || old == 0) {
      return(NA_real_)
    }
    100 * (new - old) / old
  }

  inv_E_rate <- function(sh, rt) {
    if (is.finite(sh) && sh > 0 && is.finite(rt) && rt > 0) {
      rt / sh
    } else {
      NA_real_
    }
  }

  tab <- do.call(rbind, lapply(grp, function(g) {
    ex <- existing[[g]]
    rt334 <- ex$rate
    data.frame(
      group          = g,
      n_j            = ex$n_j,
      n_prior        = ex$n_prior,
      shape_ING      = ex$shape_ING,
      rate_gamma     = ex$rate_gamma,
      rate           = rt334,
      pct_rate       = pct_diff(rt334, ex$rate_gamma),
      inv_E_gamma    = ex$inv_E,
      inv_E_rate     = inv_E_rate(ex$shape_ING, rt334),
      pct_inv_E      = pct_diff(
        inv_E_rate(ex$shape_ING, rt334),
        ex$inv_E
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(tab) <- NULL

  cat(
    "\n--- Per-group Block~1 gamma: rate_gamma (A12 3.3.5) vs rate (A12 3.3.4 S_marg) ---\n",
    "  dGamma_list downstream uses rate (S_marg); rate_gamma retained for comparison.\n",
    "  sigma2_hat and truncation bounds unchanged.\n\n",
    sep = ""
  )
  num_cols <- vapply(tab, is.numeric, logical(1))
  if (any(num_cols)) {
    tab[num_cols] <- lapply(tab[num_cols], round, digits = digits)
  }
  print(tab, row.names = FALSE)
  invisible(tab)
}

#' Central 98% prior-mass \eqn{\sigma^2}/\eqn{\tau^2} window from calibrated precision prior
#'
#' Precision \eqn{1/\sigma^2 \sim \mathrm{Gamma}(\code{shape}, \code{rate})};
#' bounds are 0.01/0.99 quantiles inverted to the variance scale.
#' @noRd
.lmebayes_ing_prior_quantile_window <- function(shape, rate, max_disp_perc = 0.99) {
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      "ING prior quantile window requires positive finite shape and rate.",
      call. = FALSE
    )
  }
  list(
    disp_lower = 1 / stats::qgamma(max_disp_perc,       shape = shape, rate = rate),
    disp_upper = 1 / stats::qgamma(1 - max_disp_perc,   shape = shape, rate = rate)
  )
}

#' Prospective \code{dGamma()} measurement \eqn{\sigma^2} calibration from setup
#'
#' Mean-matched inverse-Gamma hyperparameters for Block~1 ING (same algebra as
#' \code{ing_prior} for \eqn{\tau^2_k}, with \eqn{\hat\sigma^2} =
#' \code{group.dispersion}, \eqn{p = p_{\mathrm{re}}}, and
#' \eqn{n_{\mathrm{prior}} = \mathrm{pwt\_measurement}/(1-\mathrm{pwt\_measurement})\times n} on the total
#' observation count).  Truncation bounds are the central 98% prior-mass
#' interval from the same \code{shape}/\code{rate}.
#' @noRd
.lmebayes_calibrate_ing_prior_measurement <- function(
    design,
    group.dispersion,
    n_prior,
    max_disp_perc = 0.99
) {
  p_re <- length(design$groupef.names)
  n    <- length(design$y)
  if (p_re < 1L) {
    stop(
      "Measurement dispersion calibration requires at least one random coefficient.",
      call. = FALSE
    )
  }

  if (!is.numeric(n_prior) || length(n_prior) != 1L || !is.finite(n_prior) ||
      n_prior <= 0) {
    stop(
      "'n_prior' must be a positive finite scalar for measurement dispersion calibration.",
      call. = FALSE
    )
  }
  if (n_prior > n) {
    stop(
      "Measurement dispersion prior requires n_prior <= n; got n_prior = ",
      signif(n_prior, 4), ", n = ", n, ".",
      call. = FALSE
    )
  }

  shape <- (n_prior + 1) / 2 + p_re / 2
  rate  <- as.numeric(group.dispersion) * (n_prior + p_re - 1) / 2
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      "Measurement dispersion ING calibration produced non-positive shape/rate.",
      call. = FALSE
    )
  }

  win <- .lmebayes_ing_prior_quantile_window(shape, rate, max_disp_perc)

  list(
    sigma2_hat    = as.numeric(group.dispersion),
    shape         = shape,
    rate          = rate,
    disp_lower    = win$disp_lower,
    disp_upper    = win$disp_upper,
    max_disp_perc = max_disp_perc,
    n_prior       = n_prior,
    n_effective = n,
    p_re        = p_re
  )
}

#' Mean-match one per-group Block~1 \code{group.ing_prior} entry to a fixed \eqn{\sigma^2_j}.
#' @noRd
.lmebayes_pin_ing_prior_measurement_group_entry <- function(ing, sigma2, p_re) {
  n_prior_j <- ing$n_prior
  shape_ING <- ing$shape_ING
  rate_gamma <- sigma2 * (n_prior_j + p_re - 1) / 2
  shape_post <- shape_ING + ing$n_j / 2
  rate_post <- sigma2 * (shape_post - 1)
  win <- .lmebayes_ing_prior_quantile_window(
    shape_post, rate_post, ing$max_disp_perc
  )
  ing$sigma2_hat <- sigma2
  ing$rate_gamma <- rate_gamma
  if ("E_sigma2" %in% names(ing)) {
    ing$E_sigma2 <- if (is.finite(shape_ING) && shape_ING > 1 &&
                        is.finite(rate_gamma) && rate_gamma > 0) {
      rate_gamma / (shape_ING - 1)
    } else {
      NA_real_
    }
  }
  if ("inv_E" %in% names(ing)) {
    ing$inv_E <- if (is.finite(shape_ING) && shape_ING > 0 &&
                     is.finite(rate_gamma) && rate_gamma > 0) {
      rate_gamma / shape_ING
    } else {
      NA_real_
    }
  }
  ing$disp_lower <- win$disp_lower
  ing$disp_upper <- win$disp_upper
  ing
}

#' Mean-match every per-group Block~1 \code{group.ing_prior} entry from \code{group.dispersion}.
#' @noRd
.lmebayes_pin_ing_prior_measurement_group <- function(ing_grp, override_vec, p_re) {
  stats::setNames(
    lapply(names(ing_grp), function(lev) {
      .lmebayes_pin_ing_prior_measurement_group_entry(
        ing_grp[[lev]], unname(override_vec[[lev]]), p_re
      )
    }),
    names(ing_grp)
  )
}

#' Limiting-posterior \eqn{\sigma^2}/\eqn{\tau^2} truncation window (lmebayes default)
#'
#' Central 98% mass of \code{Gamma((J+1)/2, d_hat*(J-1)/2)} inverted to the
#' variance scale; see \code{inst/ING_TRUNCATION_WINDOW.md} in \pkg{lmebayes}.
#' @noRd
.lmebayes_ing_limiting_posterior_window <- function(d_hat, J, max_disp_perc = 0.99) {
  if (!is.numeric(d_hat) || length(d_hat) != 1L || !is.finite(d_hat) ||
      d_hat <= 0) {
    stop(
      "'d_hat' must be a positive finite scalar (classical variance plug-in).",
      call. = FALSE
    )
  }
  J <- as.integer(J[1L])
  if (!is.finite(J) || J < 1L) {
    stop("'J' must be a positive integer (number of groups).", call. = FALSE)
  }
  a_inf <- (J + 1) / 2
  b_inf <- as.numeric(d_hat) * (J - 1) / 2
  if (b_inf <= 0) {
    stop(
      "Limiting-posterior ING window requires J >= 2 (got J = ", J, ").",
      call. = FALSE
    )
  }
  list(
    disp_lower = 1 / stats::qgamma(max_disp_perc,     shape = a_inf, rate = b_inf),
    disp_upper = 1 / stats::qgamma(1 - max_disp_perc, shape = a_inf, rate = b_inf)
  )
}

#' Build shared ING Block~1 measurement \code{prior_list} for lmebayes glue
#' @noRd
.lmebayes_ing_measurement_prior_list <- function(prior, disp_info, design) {
  re_names <- design$groupef.names
  p_re     <- length(re_names)
  pl       <- disp_info$dispersion_prior_list
  if (is.null(pl$shape) || is.null(pl$rate)) {
    stop(
      "dGamma() group.dispersion prior_list must contain 'shape' and 'rate'.",
      call. = FALSE
    )
  }
  mu <- matrix(
    0,
    nrow = p_re,
    ncol = 1L,
    dimnames = list(re_names, NULL)
  )
  Sigma <- as.matrix(prior$group.Sigma)
  if (nrow(Sigma) != p_re || ncol(Sigma) != p_re) {
    stop(
      "prior$group.Sigma must be ", p_re, " x ", p_re, ".",
      call. = FALSE
    )
  }
  out <- list(
    mu            = mu,
    Sigma         = Sigma,
    shape         = pl$shape,
    rate          = pl$rate,
    max_disp_perc = if (!is.null(pl$max_disp_perc)) pl$max_disp_perc else 0.99
  )
  if (!is.null(pl$disp_lower)) out$disp_lower <- pl$disp_lower
  if (!is.null(pl$disp_upper)) out$disp_upper <- pl$disp_upper
  if (is.null(out$disp_lower) || is.null(out$disp_upper)) {
    win <- .lmebayes_ing_prior_quantile_window(
      shape = as.numeric(pl$shape[1L]),
      rate  = as.numeric(pl$rate[1L])
    )
    if (is.null(out$disp_lower)) out$disp_lower <- win$disp_lower
    if (is.null(out$disp_upper)) out$disp_upper <- win$disp_upper
  }
  if (out$disp_upper <= out$disp_lower) {
    stop(
      "dGamma() measurement prior: implied disp_upper must exceed disp_lower.",
      call. = FALSE
    )
  }
  out
}

#' Build per-group ING Block~1 measurement \code{prior_list} for lmebayes glue
#'
#' Third \code{group.dispersion} option: a named list of \code{dGamma()}
#' pfamilies, one per group. Unlike \code{.lmebayes_ing_measurement_prior_list()}
#' (single pooled \code{shape}/\code{rate}), each group keeps its own
#' \code{shape}/\code{rate}/\code{disp_lower}/\code{disp_upper}.
#' @noRd
.lmebayes_ing_measurement_prior_list_group <- function(prior, disp_info, design) {
  re_names <- design$groupef.names
  p_re     <- length(re_names)
  pl       <- disp_info$dispersion_prior_list
  req <- c("shape_group", "rate_group", "disp_lower_group", "disp_upper_group")
  miss <- req[!req %in% names(pl)]
  if (length(miss)) {
    stop(
      "Internal error: per-group dGamma() group.dispersion prior_list must ",
      "contain ", paste(req, collapse = ", "), ".",
      call. = FALSE
    )
  }
  mu <- matrix(
    0,
    nrow = p_re,
    ncol = 1L,
    dimnames = list(re_names, NULL)
  )
  Sigma <- as.matrix(prior$group.Sigma)
  if (nrow(Sigma) != p_re || ncol(Sigma) != p_re) {
    stop(
      "prior$group.Sigma must be ", p_re, " x ", p_re, ".",
      call. = FALSE
    )
  }
  list(
    mu               = mu,
    Sigma            = Sigma,
    shape_group      = pl$shape_group,
    rate_group       = pl$rate_group,
    disp_lower_group = pl$disp_lower_group,
    disp_upper_group = pl$disp_upper_group,
    max_disp_perc    = if (!is.null(pl$max_disp_perc)) pl$max_disp_perc else 0.99
  )
}

#' Assemble the flat argument list for the routed LMM \code{rLMM_reg} export
#'
#' @description
#' Shared front-door helper that turns \code{design} + the \code{prior}
#' returned by \code{\link{priors_from_pfamily_list}} + \code{disp_info}
#' (from \code{.lmebayes_resolve_group_dispersion()}) into the flat,
#' named argument list passed via \code{do.call()} to whichever
#' \code{rLMM_reg} export \code{.lmebayes_reg_route_fn()} resolves to. There
#' are exactly four possible targets, chosen by \code{disp_info$mode} and
#' \code{prior$any_non_normal}: \code{\link{rLMMNormal_reg_known_vcov}},
#' \code{rLMMNormal_reg_estimated_vcov}, \code{rLMMindepNormalGamma_reg_known_vcov},
#' or \code{rLMMindepNormalGamma_reg_estimated_vcov} (the actual dispatch
#' happens one level up, in \code{.lmebayes_run_lmm_engine()} via
#' \code{.lmebayes_reg_route_key()}/\code{.lmebayes_reg_route_fn()} --
#' this function only builds the argument list, it does not pick the route).
#'
#' @details
#' Concretely, this function does three things:
#' \enumerate{
#'   \item Builds a base argument list that is always the same shape:
#'     the response/design matrices and grouping structure taken from
#'     \code{design}, \code{prior$pfamily_list} unchanged (the routed export
#'     derives its own Block~2 random-effect prior precision from
#'     \code{pfamily_list} internally, so it is not built here), and
#'     the naming/control arguments (\code{tv_tol}, \code{progbar},
#'     \code{verbose}), with \code{group_name} attached to \code{group} as
#'     an attribute rather than a separate argument.
#'   \item Builds \code{args$prior_list} (the Block~1/measurement prior),
#'     whose \emph{shape} depends on \code{disp_info$mode}: a plain known
#'     \code{dispersion} value for \code{"none"}/\code{"fixed"}/
#'     \code{"fixed_vector"}, or a shared/per-group ING measurement prior
#'     (built by \code{.lmebayes_ing_measurement_prior_list()}/
#'     \code{.lmebayes_ing_measurement_prior_list_group()}) for
#'     \code{"gamma"}/\code{"gamma_list"}.
#'   \item Conditionally adds route-specific controls: pilot-stage controls
#'     (\code{gap_tol}, \code{mode_gap_max}, \code{diag_sweeps},
#'     \code{stage_verbose}) when \code{prior$any_non_normal} is
#'     \code{TRUE} (an ING route, which needs a pilot chain), or
#'     \code{sim_method} when the route is fixed \emph{and} known-vcov
#'     (the only route with more than one sampling engine).
#' }
#' It currently combines all of the above (base argument assembly,
#' route-specific \code{prior_list} branching, and conditional pilot-stage/
#' \code{sim_method} fields) into a single flat list and is a refactor
#' candidate; its argument list and return shape are not yet considered
#' stable.
#'
#' @param n Number of stored draws.
#' @param design A \code{model_setup} object.
#' @param prior The list returned by \code{\link{priors_from_pfamily_list}}.
#' @param disp_info Resolved dispersion info, as returned by
#'   \code{.lmebayes_resolve_group_dispersion()} (must supply
#'   \code{mode} and \code{dispersion_fix}).
#' @param tv_tol Total-variation tolerance passed through to the routed
#'   export.
#' @param progbar,verbose Passed through to the routed export.
#' @param gap_tol,mode_gap_max,diag_sweeps Pilot-stage controls, added to the
#'   argument list only when \code{prior$any_non_normal} is \code{TRUE}
#'   (i.e. an ING route).
#' @param sim_method Sampling engine (\code{"DEFAULT"} or
#'   \code{"TWO_BLOCK_GIBBS"}), added to the argument list only on the
#'   fixed+known-vcov (\code{rLMMNormal_reg_known_vcov}) route.
#' @return A flat named list (\code{args}) suitable for
#'   \code{do.call(route$export_fn, args)}, with elements:
#'   \describe{
#'     \item{\code{n}}{Number of stored draws, unchanged from the \code{n}
#'       argument.}
#'     \item{\code{y}}{\code{design$y}, the response vector.}
#'     \item{\code{D}}{\code{design$D}, the level-1 (\eqn{l_2 \times p_{re}})
#'       random-effect design matrix.}
#'     \item{\code{group}}{\code{design$group}, the grouping factor, with
#'       \code{attr(group, "group_name")} set to \code{design$group_name}
#'       (the routed export has no \code{group_name} formal and resolves it
#'       from this attribute, since \code{design$group} is never a bare
#'       variable at the routed export's call site).}
#'     \item{\code{W}}{\code{design$W}, the named list of
#'       group-level hyper-design matrices (one per random-effect
#'       coefficient).}
#'     \item{\code{pfamily_list}}{\code{prior$pfamily_list} unchanged. The
#'       routed export derives its own Block~2 random-effect prior
#'       precision from this internally; it is not built or passed here.}
#'     \item{\code{tv_tol}}{The \code{tv_tol} argument, unchanged.}
#'     \item{\code{progbar}, \code{verbose}}{The \code{progbar}/\code{verbose}
#'       arguments, unchanged.}
#'     \item{\code{pop.prior_list}}{Always present, but its shape depends on
#'       \code{disp_info$mode}: for \code{"none"}/\code{"fixed"}/
#'       \code{"fixed_vector"} it is \code{list(dispersion = disp_info$dispersion_fix)}
#'       (a scalar or named per-group numeric vector, passed through
#'       unchanged -- a \emph{known}, fixed observation dispersion); for
#'       \code{"gamma"} it is the pooled ING measurement prior list from
#'       \code{.lmebayes_ing_measurement_prior_list()} (\code{mu},
#'       \code{Sigma}, \code{shape}, \code{rate}, \code{max_disp_perc},
#'       \code{disp_lower}, \code{disp_upper}); for \code{"gamma_list"} it
#'       is the per-group ING measurement prior list from
#'       \code{.lmebayes_ing_measurement_prior_list_group()} (\code{mu},
#'       \code{Sigma}, \code{shape_group}, \code{rate_group},
#'       \code{disp_lower_group}, \code{disp_upper_group},
#'       \code{max_disp_perc}).}
#'     \item{\code{gap_tol}, \code{mode_gap_max}, \code{diag_sweeps},
#'       \code{stage_verbose}}{Present \strong{only} when
#'       \code{prior$any_non_normal} is \code{TRUE} (an ING/estimated-vcov
#'       route that needs a pilot chain): the \code{gap_tol}/
#'       \code{mode_gap_max}/\code{diag_sweeps} arguments unchanged, and
#'       \code{stage_verbose} set to \code{verbose}. Absent for
#'       \code{"none"}/\code{"fixed"}/\code{"fixed_vector"} dispersion
#'       modes.}
#'     \item{\code{sim_method}}{Present \strong{only} when
#'       \code{prior$any_non_normal} is \code{FALSE} \emph{and}
#'       \code{disp_info$mode} is not \code{"gamma"}/\code{"gamma_list"}
#'       (i.e. only on the route that resolves to
#'       \code{\link{rLMMNormal_reg_known_vcov}}, the only export with a
#'       real \code{sim_method} dispatch between exact iid draws and
#'       two-block Gibbs sweeps): the \code{sim_method} argument, unchanged.
#'       Every other route ignores \code{sim_method} (if it even accepts
#'       it), so it is only forwarded here.}
#'   }
#' @keywords internal
#' @export
matrix_args_lmm <- function(
    n,
    design,
    prior,
    disp_info,
    tv_tol        = 0.01,
    progbar       = TRUE,
    verbose       = FALSE,
    gap_tol       = 0.0196,
    mode_gap_max  = 1.0,
    diag_sweeps   = FALSE,
    sim_method    = "DEFAULT"
) {
  ## The routed export has no 'group_name' formal: attach it to 'group'
  ## itself (design$group is never a bare variable here, so the export's
  ## substitute()-based fallback could not resolve it anyway).
  grp <- design$group
  attr(grp, "group_name") <- design$group_name

  args <- list(
    n             = n,
    y             = design$y,
    D             = design$D,
    group         = grp,
    W             = design$W,
    pfamily_list  = prior$pfamily_list,
    tv_tol        = tv_tol,
    progbar       = progbar,
    verbose       = verbose
  )

  if (identical(disp_info$mode, "gamma")) {
    args$prior_list <- .lmebayes_ing_measurement_prior_list(
      prior     = prior,
      disp_info = disp_info,
      design    = design
    )
  } else if (identical(disp_info$mode, "gamma_list")) {
    args$prior_list <- .lmebayes_ing_measurement_prior_list_group(
      prior     = prior,
      disp_info = disp_info,
      design    = design
    )
  } else {
    args$prior_list <- list(dispersion = disp_info$dispersion_fix)
  }

  if (isTRUE(prior$any_non_normal)) {
    args$gap_tol       <- gap_tol
    args$mode_gap_max  <- mode_gap_max
    args$diag_sweeps   <- diag_sweeps
    args$stage_verbose <- verbose
  } else if (!identical(disp_info$mode, "gamma") &&
             !identical(disp_info$mode, "gamma_list")) {
    ## lmm_fixed_known route only: rLMMNormal_reg_known_vcov() is the only
    ## export with a real sim_method dispatch (exact iid vs. two-block
    ## Gibbs). Every other route accepts the argument (if at all) as a
    ## no-op, so it is only forwarded here.
    args$sim_method <- sim_method
  }

  args
}

#' Shared matrix-level arguments for GLMM reg routes
#' @noRd
.lmebayes_matrix_args_glmm <- function(
    n,
    design,
    prior,
    family,
    gap_tol       = 0.0196,
    tv_tol        = 0.01,
    mode_gap_max  = 1.0,
    verbose       = FALSE,
    progbar       = FALSE,
    collect_block1 = TRUE
) {
  block1_prior <- .lmebayes_block1_prior_list(prior, group.dispersion = NULL)

  ## The routed export has no 'group_name' formal: attach it to 'group'
  ## itself (design$group is never a bare variable here, so the export's
  ## substitute()-based fallback could not resolve it anyway).
  grp <- design$group
  attr(grp, "group_name") <- design$group_name

  list(
    n               = n,
    y               = design$y,
    D               = design$D,
    group           = grp,
    W               = design$W,
    prior_list      = block1_prior,
    pfamily_list    = prior$pfamily_list,
    family          = family,
    gap_tol         = gap_tol,
    tv_tol          = tv_tol,
    mode_gap_max    = mode_gap_max,
    verbose         = verbose,
    progbar         = progbar,
    stage_verbose   = verbose,
    collect_block1  = collect_block1
  )
}

#' @noRd
.lmebayes_run_lmm_engine <- function(
    n,
    design,
    prior,
    disp_info,
    tv_tol        = 0.01,
    progbar       = TRUE,
    verbose       = FALSE,
    gap_tol             = 0.0196,
    mode_gap_max        = 1.0,
    diag_sweeps         = FALSE,
    sim_method          = "DEFAULT"
) {
  route_key <- .lmebayes_reg_route_key(
    family         = gaussian(),
    disp_mode      = disp_info$mode,
    any_non_normal = prior$any_non_normal
  )
  route <- .lmebayes_reg_route_fn(route_key)
  args  <- matrix_args_lmm(
    n             = n,
    design        = design,
    prior         = prior,
    disp_info     = disp_info,
    tv_tol        = tv_tol,
    progbar       = progbar,
    verbose       = verbose,
    gap_tol       = gap_tol,
    mode_gap_max  = mode_gap_max,
    diag_sweeps   = diag_sweeps,
    sim_method    = sim_method
  )
  out <- do.call(route$export_fn, args)
  if (is.null(out$sim_method_used)) {
    out$sim_method_used <- "TWO_BLOCK_GIBBS"
  }
  .lmebayes_attach_sigma2(out, disp_info)
}

#' @noRd
.lmebayes_run_glmm_engine <- function(
    n,
    design,
    prior,
    family,
    gap_tol       = 0.0196,
    tv_tol        = 0.01,
    mode_gap_max  = 1.0,
    verbose       = FALSE,
    progbar       = FALSE,
    collect_block1 = TRUE
) {
  route_key <- .lmebayes_reg_route_key(
    family         = family,
    disp_mode      = "none",
    any_non_normal = prior$any_non_normal
  )
  route <- .lmebayes_reg_route_fn(route_key)
  args  <- .lmebayes_matrix_args_glmm(
    n              = n,
    design         = design,
    prior          = prior,
    family         = family,
    gap_tol        = gap_tol,
    tv_tol         = tv_tol,
    mode_gap_max   = mode_gap_max,
    verbose        = verbose,
    progbar        = progbar,
    collect_block1 = collect_block1
  )
  out <- do.call(route$export_fn, args)
  disp_none <- list(mode = "none")
  .lmebayes_attach_sigma2(out, disp_none)
}

#' Build a P/Sigma-free Block~1 prior list (\code{dispersion}/\code{ddef} only)
#'
#' The Block~1 random-effect prior precision is now always derived
#' internally (from \code{pfamily_list}) by the routed \code{rGLMM_reg}/
#' \code{rLMM_reg} exports, which reject a caller-supplied \code{P}/
#' \code{Sigma}; this helper therefore no longer computes or returns one
#' (previously \code{solve(measurement_prior_list$group.Sigma)}).
#' @noRd
.lmebayes_block1_prior_list <- function(
    measurement_prior_list,
    group.dispersion = NULL
) {
  dispersion <- if (!is.null(group.dispersion)) {
    group.dispersion
  } else {
    measurement_prior_list$group.dispersion
  }
  if (is.null(dispersion)) {
    list(ddef = TRUE)
  } else {
    list(dispersion = dispersion, ddef = FALSE)
  }
}

#' Attach \code{sigma2} / \code{sigma2.mean} from dispersion mode and sampler draws.
#'
#' Fixed measurement dispersion returns a scalar; a fixed per-group vector
#' returns that same named length-\code{J} vector (constant, not sampled); a
#' single \code{dGamma()} returns the length-\code{n} vector from the final
#' inner sweep (\code{group.dispersion}); a list of per-group \code{dGamma()}
#' pfamilies returns an \code{n x J} matrix (one column per group); families
#' without observation-level dispersion get \code{NULL}.
#' @noRd
.lmebayes_attach_sigma2 <- function(out, disp_info) {
  mode <- disp_info$mode
  if (identical(mode, "none")) {
    out$sigma2 <- NULL
    out$sigma2.mean <- NULL
    return(out)
  }
  if (identical(mode, "fixed")) {
    val <- as.numeric(disp_info$dispersion_fix)
    out$sigma2 <- val
    out$sigma2.mean <- val
    return(out)
  }
  if (identical(mode, "fixed_vector")) {
    val <- disp_info$dispersion_fix
    out$sigma2 <- val
    out$sigma2.mean <- val
    return(out)
  }
  if (identical(mode, "gamma")) {
    ## Fit draw field remains dispersion_ranef until a later rename pass.
    dr <- out$dispersion_ranef
    if (is.null(dr)) {
      stop(
        "Internal error: dGamma measurement dispersion requires ",
        "'dispersion_ranef' draws on the sampler output.",
        call. = FALSE
      )
    }
    out$sigma2 <- as.numeric(dr)
    out$sigma2.mean <- mean(out$sigma2)
    return(out)
  }
  if (identical(mode, "gamma_list")) {
    dr <- out$dispersion_ranef
    if (is.null(dr)) {
      stop(
        "Internal error: per-group dGamma() measurement dispersion requires ",
        "'dispersion_ranef' draws on the sampler output.",
        call. = FALSE
      )
    }
    out$sigma2 <- as.matrix(dr)
    out$sigma2.mean <- colMeans(out$sigma2)
    return(out)
  }
  stop("Unknown dispersion mode: ", mode, call. = FALSE)
}

#' @noRd
.lmebayes_add_fixef_summaries <- function(x) {
  if (!is.null(x$fixef)) {
    x$fixef.means <- lapply(x$fixef, colMeans)
  }
  if (!is.null(x$fixef.dispersion)) {
    x$fixef.dispersion.mean <- colMeans(x$fixef.dispersion)
  }
  if (!is.null(x$fixef.iters) && !is.null(x$m_convergence)) {
    x$fixef.iters.mean <- colMeans(x$fixef.iters) / x$m_convergence
  }
  if (!is.null(x$ranef.iters) && !is.null(x$m_convergence)) {
    x$ranef.iters.mean <- mean(x$ranef.iters) / x$m_convergence
  }
  if (!is.null(x$sigma2.iters) && !is.null(x$m_convergence)) {
    x$sigma2.iters.mean <- colMeans(x$sigma2.iters) / x$m_convergence
  }
  x
}

#' @noRd
.lmebayes_block2_icm_labels <- function(prior, family = gaussian()) {
  any_ing <- isTRUE(prior$any_non_normal)
  is_gauss <- is.null(family) || identical(family$family, "gaussian")
  ref_label <- "prior mean"
  if (any_ing) {
    icm_label   <- "gamma @ lmer tau2"
    icm_verbose <- "Block 2 start at lmer tau^2 plug-in"
    conv_label  <- "Plug-in fixed point"
  } else if (is_gauss) {
    icm_label   <- "ICM mean"
    icm_verbose <- "ICM posterior mean"
    conv_label  <- "ICM"
  } else {
    icm_label   <- "ICM mode"
    icm_verbose <- "ICM posterior mode"
    conv_label  <- "ICM"
  }
  list(
    ref_label   = ref_label,
    icm_label   = icm_label,
    icm_verbose = icm_verbose,
    conv_label  = conv_label
  )
}

#' Per-group full column-rank flag for Block~1 \code{Z_j} (same rule as
#' \code{model_setup()$groupef.rank}).
#' @noRd
.lmebayes_groupef_rank_from_Z <- function(Z, groups, group_levels = NULL) {
  Z <- as.matrix(Z)
  g_chr <- as.character(groups)
  levs <- if (is.null(group_levels)) {
    unique(g_chr)
  } else {
    as.character(group_levels)
  }
  p_re <- ncol(Z)
  stats::setNames(
    vapply(
      levs,
      function(lev) {
        rows <- which(g_chr == lev)
        Z_j  <- Z[rows, , drop = FALSE]
        nrow(Z_j) >= p_re &&
          Matrix::rankMatrix(Z_j, method = "qr")[1L] == p_re
      },
      logical(1L)
    ),
    levs
  )
}

#' @noRd
.lmebayes_print_icm_fixef_table <- function(
    prior_list,
    re_names,
    fixef_icm,
    icm_info,
    ref_label,
    icm_label,
    conv_label = "ICM",
    header,
    verbose
) {
  if (!isTRUE(verbose) || is.null(fixef_icm)) {
    return(invisible(NULL))
  }
  fixef_ref <- lapply(prior_list, `[[`, "mu")
  names(fixef_ref) <- re_names
  hdr <- sprintf("  %-18s  %-30s  %14s  %18s",
                 "RE component", "parameter", ref_label, icm_label)
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat(header, "\n")
  cat(hdr, "\n")
  cat(sep, "\n")
  for (k in re_names) {
    nms_k  <- names(fixef_ref[[k]])
    ref_v  <- fixef_ref[[k]]
    icm_v  <- fixef_icm[[k]]
    for (nm in nms_k) {
      cat(sprintf("  %-18s  %-30s  %14.4f  %18.4f\n",
                  k, nm, ref_v[[nm]], icm_v[[nm]]))
    }
  }
  if (!is.null(icm_info)) {
    cat(sprintf("  (%s converged: %s, %d iter, delta = %.2e)\n\n",
                conv_label,
                icm_info$converged, icm_info$iterations, icm_info$delta))
  } else {
    cat("\n")
  }
  invisible(NULL)
}

#' @noRd
.lmebayes_print_ranef_mode_reference <- function(
    ranef_mode,
    re_names,
    group_levels,
    verbose
) {
  invisible(NULL)
}

#' @noRd
.lmebayes_print_fixef_init <- function(
    fixef_init,
    re_names,
    verbose,
    header = "--- main-stage fixef.init (pilot colMeans) ---"
) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat(header, "\n")
  for (k in re_names) {
    for (nm in names(fixef_init[[k]])) {
      cat(sprintf("  %-18s  %-30s  %12.4f\n",
                  k, nm, fixef_init[[k]][[nm]]))
    }
  }
  cat("\n")
  invisible(NULL)
}
