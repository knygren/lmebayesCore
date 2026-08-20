## Conditional means E[beta_j | gamma, y] for the C05 population mean map.

#' Conditional means of group random effects given population effects
#'
#' Computes the C05 EM E-step: per-group conditional means
#' \eqn{E[\beta_j \mid \gamma, y]} (and posterior covariances \eqn{V_j}) at a
#' supplied population location \code{fixef}. For \code{gaussian()} with
#' \code{estep = "exact"} and no \code{beta_set}, the posterior mean is closed
#' form. Otherwise means are Monte Carlo averages.
#'
#' @details
#' When \code{beta_set} (\eqn{\widetilde B(\delta_2)}) is supplied, means are
#' Monte Carlo averages of accepted joint draws from
#' \code{\link{rNormal_reg_group_safe}} or
#' \code{\link{rNormalGLM_reg_group_safe}} (by \code{family}).
#'
#' @param design A \code{\link{model_setup}} list (required if \code{mode} is
#'   \code{NULL}).
#' @param pfamily_list Block~2 prior list (required if \code{mode} is
#'   \code{NULL}).
#' @param family A \code{\link[stats]{family}} object (required if \code{mode}
#'   is \code{NULL}).
#' @param fixef Named list of population fixed effects (\eqn{\gamma} in C05
#'   notation). Ignored when \code{mode} is supplied.
#' @param mode Optional result of \code{\link{population_mode}}; supplies
#'   \code{fixef}, \code{design}, \code{p11}, and model context.
#' @param dispprior_list Optional Block~1 dispersion prior for \code{gaussian()}.
#' @param beta_set Optional \code{\link{beta_marginal_safe_set}} object
#'   (\eqn{\widetilde B(\delta_2)}); enables joint-draw rejection (see Details).
#' @param estep \code{"exact"} (Gaussian closed form when unrestricted),
#'   \code{"mc"} (simulated means), or \code{"aghq"} (not yet implemented).
#'   Ignored when \code{beta_set} is supplied (always uses rejection MC).
#' @param n Target number of **accepted joint draws** when \code{beta_set} is
#'   set; otherwise number of \eqn{\beta_j} draws per group when
#'   \code{estep = "mc"} (required for non-Gaussian families).
#' @param mc_seed Optional RNG seed for the MC E-step.
#' @return An object of class \code{"group_effects_conditional_mean"} with
#'   \code{b_mean}, \code{V_list}, \code{b_mc_se} (when MC is used),
#'   \code{fixef}, \code{estep}, \code{n}, \code{n_target}, \code{n_tried},
#'   \code{accept_rate} (when \code{beta_set} is set), \code{beta_set},
#'   \code{restricted}, and model context (\code{design}, \code{family},
#'   \code{p11}).
#' @seealso \code{\link{population_mode}}, \code{\link{beta_marginal_safe_set}},
#'   \code{\link{rNormal_reg_group_safe}}, \code{\link{rNormalGLM_reg_group_safe}}
#' @export
group_effects_conditional_mean <- function(design = NULL,
                                           pfamily_list = NULL,
                                           family = NULL,
                                           fixef = NULL,
                                           mode = NULL,
                                           dispprior_list = NULL,
                                           beta_set = NULL,
                                           estep = c("exact", "aghq", "mc"),
                                           n = 10000L,
                                           mc_seed = NULL) {
  estep <- match.arg(estep)
  n <- as.integer(n)
  if (!(n >= 1L)) {
    stop("'n' must be at least 1.", call. = FALSE)
  }

  if (!is.null(beta_set) && !inherits(beta_set, "beta_marginal_safe_set")) {
    stop("'beta_set' must be from beta_marginal_safe_set().", call. = FALSE)
  }

  if (!is.null(mode)) {
    if (!is.null(fixef)) {
      stop("Supply only one of 'fixef' and 'mode'.", call. = FALSE)
    }
    fixef <- mode$fixef
    design <- mode$design
    family <- mode$family
    mpl <- mode$measurement_prior_list
    p11 <- mode$p11
  } else {
    if (is.null(design) || is.null(pfamily_list) || is.null(family)) {
      stop(
        "'design', 'pfamily_list', and 'family' are required when 'mode' is NULL.",
        call. = FALSE
      )
    }
    if (is.null(fixef)) {
      stop("'fixef' is required when 'mode' is NULL.", call. = FALSE)
    }
    prep <- .c05_validate(design, pfamily_list, family, dispprior_list)
    mpl <- prep$measurement_prior_list
    p11 <- .c05_p11(
      design,
      prep$prior_pack,
      prep$group_levels,
      measurement_prior_list = mpl,
      family = family
    )
  }

  estep_use <- estep
  if (!is.null(beta_set)) {
    estep_use <- "mc"
  } else if (!identical(family$family, "gaussian") && identical(estep, "exact")) {
    estep_use <- "mc"
    message(
      "Non-Gaussian family: using estep = \"mc\" with n = ", n, " draws per group.",
      appendLF = FALSE
    )
  }

  estep_out <- .c05_estep(
    design = design,
    fixef = fixef,
    p11 = p11,
    measurement_prior_list = mpl,
    family = family,
    estep = estep_use,
    n = n,
    mc_seed = mc_seed,
    beta_set = beta_set
  )

  n_report <- if (isTRUE(estep_out$restricted)) {
    estep_out$n
  } else if (identical(estep_out$estep, "mc")) {
    n
  } else {
    NA_integer_
  }

  structure(
    list(
      b_mean = estep_out$b_mean,
      b_mc_se = estep_out$b_mc_se,
      V_list = estep_out$V_list,
      fixef = fixef,
      estep = estep_out$estep,
      n = n_report,
      n_target = estep_out$n_target,
      n_tried = estep_out$n_tried,
      accept_rate = estep_out$accept_rate,
      beta_set = beta_set,
      restricted = isTRUE(estep_out$restricted),
      design = design,
      family = family,
      measurement_prior_list = mpl,
      p11 = p11,
      call = match.call()
    ),
    class = "group_effects_conditional_mean"
  )
}

#' @export
print.group_effects_conditional_mean <- function(x, digits = 4, ...) {
  cat("Group effects conditional mean E[beta_j | gamma, y]\n\n")
  cat("  estep: ", x$estep, sep = "")
  if (isTRUE(x$restricted)) {
    cat("  restricted to widetilde B\n")
    cat("  accepted: ", x$n,
        if (!is.na(x$n_target)) paste0(" / target ", x$n_target) else "",
        if (is.finite(x$accept_rate)) {
          paste0("  rate: ", signif(x$accept_rate, digits))
        } else {
          ""
        },
        "\n", sep = "")
  } else if (identical(x$estep, "mc") && !is.na(x$n)) {
    cat("  n per group: ", x$n, "\n", sep = "")
  } else {
    cat("\n")
  }
  if (!is.null(x$beta_set) && !isTRUE(x$restricted)) {
    cat("  beta_set delta_2: ", x$beta_set$delta_2, "\n", sep = "")
  }
  cat("  ||b_mean||: ",
      signif(sqrt(sum(x$b_mean^2)), digits), "\n", sep = "")
  invisible(x)
}
