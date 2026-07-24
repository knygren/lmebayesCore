#' Plot Block~2 sweep-history diagnostics (cross-chain mean or SD)
#'
#' Plots cross-chain Block~2 hyperparameter summaries stored on an object of
#' class \code{"two_block_sweep_history"} (typically
#' \code{fit$sweep_history$pilot} or \code{fit$sweep_history$main} from
#' \code{\link{rlmerb}}, \code{\link{rglmerb}}, or formula drivers such as
#' \code{lmerb()} / \code{glmerb()} in \pkg{lmebayes}; see
#' \code{\link{print.two_block_sweep_history}}).
#'
#' @param hist Object of class \code{"two_block_sweep_history"} (see
#'   \code{\link{print.two_block_sweep_history}}).
#' @param coef_focus List of length-2 character vectors
#'   \code{c(re_component, covariate)}, matching rows of \code{hist$table}.
#'   Example: \code{list(c("(Intercept)", "(Intercept)"), c("violent_i", "(Intercept)"))}.
#' @param what One or both of \code{"sd"} and \code{"mean"} (cross-chain
#'   summary after each inner sweep). Default both (\code{sd} first).
#' @param engine \code{"base"} for one panel per coefficient (default), or
#'   \code{"ggplot"} for a single faceted figure (requires \pkg{ggplot2}).
#' @param stage_label Character label for titles; defaults to \code{hist$stage}.
#' @return \code{hist} invisibly.
#' @seealso \code{\link{rlmerb}}, \code{\link{rglmerb}},
#'   \code{\link{print.two_block_sweep_history}}
#' @export
plot_sweep_history_diag <- function(
    hist,
    coef_focus,
    what = c("sd", "mean"),
    engine = c("base", "ggplot"),
    stage_label = hist$stage
) {
  if (!inherits(hist, "two_block_sweep_history")) {
    stop(
      "'hist' must be a two_block_sweep_history object ",
      "(e.g. fit$sweep_history$main).",
      call. = FALSE
    )
  }
  if (!is.list(coef_focus) || !length(coef_focus)) {
    stop("'coef_focus' must be a non-empty list of c(re_component, covariate) pairs.",
         call. = FALSE)
  }

  what <- match.arg(what, c("sd", "mean"), several.ok = TRUE)
  what <- what[order(match(what, c("sd", "mean")))]
  engine <- match.arg(engine)
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) {
    stage_label <- if (!is.null(hist$stage)) hist$stage else "stage"
  }

  if (identical(engine, "ggplot") && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("'engine = \"ggplot\"' requires the ggplot2 package.", call. = FALSE)
  }

  sh_tab <- hist$table
  sh_tab$re_component <- as.character(sh_tab$re_component)
  sh_tab$covariate <- as.character(sh_tab$covariate)
  sh_sweeps <- sh_tab[sh_tab$sweep > 0L, , drop = FALSE]
  if (!nrow(sh_sweeps)) {
    warning("No sweep rows in sweep history for stage ", stage_label, call. = FALSE)
    return(invisible(hist))
  }

  sh_plot <- do.call(rbind, lapply(coef_focus, function(cc) {
    if (length(cc) < 2L) {
      stop("Each element of 'coef_focus' must be c(re_component, covariate).",
           call. = FALSE)
    }
    sh_sweeps[
      sh_sweeps$re_component == as.character(cc[1L]) &
        sh_sweeps$covariate == as.character(cc[2L]),
      ,
      drop = FALSE
    ]
  }))
  rownames(sh_plot) <- NULL

  plot_one <- function(re_comp, cov, metric, ylab) {
    sub <- sh_sweeps[
      sh_sweeps$re_component == re_comp & sh_sweeps$covariate == cov,
      ,
      drop = FALSE
    ]
    sub <- sub[order(sub$sweep), , drop = FALSE]
    if (!nrow(sub)) {
      warning("No sweep rows for ", re_comp, " | ", cov, call. = FALSE)
      return(invisible(FALSE))
    }
    y <- if (metric == "sd") sub$sd else sub$mean
    if (!any(is.finite(y))) {
      warning("No finite ", metric, " values for ", re_comp, " | ", cov,
              call. = FALSE)
      return(invisible(FALSE))
    }
    graphics::plot(
      sub$sweep, y,
      type = "b", pch = 16,
      xlab = "Inner sweep", ylab = ylab,
      main = paste(re_comp, cov, sep = " | ")
    )
    if (metric == "mean") {
      mode_val <- sh_tab$mean[
        sh_tab$re_component == re_comp &
          sh_tab$covariate == cov &
          sh_tab$sweep == 0L
      ]
      if (length(mode_val) == 1L && is.finite(mode_val)) {
        graphics::abline(h = mode_val, lty = 2, col = "gray40")
      }
    }
    invisible(TRUE)
  }

  for (metric in what) {
    ylab <- if (metric == "sd") "Cross-chain SD" else "Cross-chain mean"
    cat(sprintf(
      "\n=== %s sweep history (%s; %s) ===\n\n",
      stage_label, ylab, engine
    ))

    if (identical(engine, "base")) {
      op <- graphics::par(
        mfrow = c(length(coef_focus), 1L),
        mar = c(4, 4, 2.5, 1),
        oma = c(0, 0, 2, 0)
      )
      n_plotted <- 0L
      for (cc in coef_focus) {
        if (isTRUE(plot_one(as.character(cc[1L]), as.character(cc[2L]), metric, ylab))) {
          n_plotted <- n_plotted + 1L
        }
      }
      graphics::mtext(
        sprintf("%s Block 2 fixef: cross-chain %s by inner sweep", stage_label, metric),
        outer = TRUE, line = 0.5, cex = 0.95
      )
      if (metric == "mean") {
        graphics::mtext(
          "Dashed line = ICM mode (sweep 0)",
          outer = TRUE, line = -1.5, cex = 0.85
        )
      }
      graphics::par(op)
      if (n_plotted == 0L) {
        warning(
          "No panels drawn for ", metric, ". Check coef_focus against ",
          "unique(hist$table[, c('re_component', 'covariate')]).",
          call. = FALSE
        )
      }
    } else if (nrow(sh_plot)) {
      sh_plot$coef <- interaction(
        sh_plot$re_component, sh_plot$covariate, sep = " | "
      )
      y_var <- if (metric == "sd") "sd" else "mean"
      aes_y <- switch(
        y_var,
        sd = sh_plot$sd,
        mean = sh_plot$mean
      )
      p <- ggplot2::ggplot(
        sh_plot,
        ggplot2::aes(sweep, aes_y, group = coef, colour = coef)
      ) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(~ coef, scales = "free_y") +
        ggplot2::labs(
          x = "Inner sweep",
          y = ylab,
          title = sprintf(
            "%s Block 2 fixef - cross-chain %s by sweep",
            stage_label, metric
          )
        ) +
        ggplot2::theme(legend.position = "none")
      if (metric == "mean") {
        mode_df <- do.call(rbind, lapply(coef_focus, function(cc) {
          sh_tab[
            sh_tab$re_component == as.character(cc[1L]) &
              sh_tab$covariate == as.character(cc[2L]) &
              sh_tab$sweep == 0L,
            ,
            drop = FALSE
          ]
        }))
        if (nrow(mode_df)) {
          mode_df$coef <- interaction(
            mode_df$re_component, mode_df$covariate, sep = " | "
          )
          mode_df$yint <- mode_df$mean
          p <- p + ggplot2::geom_hline(
            ggplot2::aes(yintercept = yint, linetype = "ICM mode"),
            data = mode_df,
            colour = "gray40"
          ) +
            ggplot2::scale_linetype_manual(
              name = NULL, values = c("ICM mode" = "dashed")
            )
        }
      }
      print(p)
    }
  }

  invisible(hist)
}

#' Combined Var/Var_final ratio plot for Block~2 fixed effects (Claim 3)
#'
#' Tracks convergence of the Block~2 fixed-effects (hyperparameter) posterior
#' variance across inner Gibbs sweeps as a single combined chart (all traces
#' overlaid on one axis, following \code{lmebayes::Ex_16}'s style), per the
#' two-block Gibbs ergodicity result (Claim~3 of the package's ergodicity
#' reference, \code{inst/BLOCK_GIBBS_ERGODICITY.md}): after \code{l} sweeps,
#' the cross-chain covariance \eqn{\Sigma^{(l)}_{11}} satisfies
#' \deqn{
#'   \Sigma_{11}^{-1/2} \Sigma^{(l)}_{11} \Sigma_{11}^{-1/2} = I - A^{2l}
#' }
#' for a convergent (symmetric, eigenvalues in \eqn{[0, 1)}) matrix \eqn{A}.
#' \code{Var(l) / Var_ref} ratios are therefore bounded above by 1 and
#' increase toward it as \code{l} grows.
#'
#' @param hist Object of class \code{"two_block_sweep_history"} (see
#'   \code{\link{print.two_block_sweep_history}}). For \code{whitened = TRUE},
#'   \code{hist} must also carry \code{cov_by_sweep}/\code{coef_index}, which
#'   are only populated by sweeps-outer/chains-inner engines (currently the
#'   \verb{_run_with_pilot()} family: \code{rLMMNormal_reg_estimated_vcov()},
#'   \code{rLMMindepNormalGamma_reg_known_vcov()},
#'   \code{rLMMindepNormalGamma_reg_estimated_vcov()}) -- not
#'   \code{rLMMNormal_reg_known_vcov(sim_method = "TWO_BLOCK_GIBBS")}, whose
#'   engine runs entirely inside compiled code and captures no per-sweep
#'   history yet.
#' @param coef_focus Optional list of length-2 character vectors
#'   \code{c(re_component, covariate)} to restrict the plot to; default
#'   (\code{NULL}) uses every coefficient tracked in \code{hist}.
#' @param design,measurement_prior_list When both are supplied, the *exact*
#'   reference covariance is resolved automatically via
#'   \code{\link{lmerb_posterior_covariance}} and used as the denominator --
#'   this is not an opt-in toggle, it is used whenever the ingredients for it
#'   exist. Only valid when dispersion and the RE variance-covariance are
#'   both fixed (not sampled). When omitted (the only option for estimated
#'   dispersion or estimated-vcov models), the denominator falls back to the
#'   empirical last-sweep cross-chain variance/covariance ("Var_final").
#' @param whitened When \code{FALSE} (default), plots \code{Var(l) / Var_ref}
#'   per named coefficient in \code{coef_focus}. When \code{TRUE}, whitens
#'   each sweep's full cross-chain covariance by the reference covariance
#'   (\eqn{\Sigma_{ref}^{-1/2} \mathrm{Cov}(l) \Sigma_{ref}^{-1/2}}) and plots
#'   its eigenvalues -- basis-invariant, so no eigenvectors of the
#'   theoretical \eqn{A} are ever needed -- labeled \code{var1, var2, ...}
#'   instead of named coefficients (per Claim~3, these eigen-directions mix
#'   the original coefficients together, so named labels no longer apply).
#' @param engine \code{"base"} (default) or \code{"ggplot"} (requires
#'   \pkg{ggplot2}); both draw one combined chart with all traces overlaid.
#' @param n_chains Optional integer, the number of independent chains each
#'   plotted \code{Var(l)} (and \code{Var_ref}, when it is the empirical
#'   \code{Var_final} rather than an exact reference) was computed across --
#'   i.e. the same \code{n_chains}/\code{n_pilot}/\code{n} passed to the
#'   sweep-generating call (see Details). When supplied, draws a single pair
#'   of horizontal dotted reference lines: a naive \code{conf_level} band
#'   around 1 for \code{Ratio(l)}, under the null that the true ratio at that
#'   sweep already equals 1 (unrealistic at small \code{l}, but useful as a
#'   ruler for how much of the departure from 1 is just Monte Carlo noise).
#'   The band's distribution depends on whether \code{Var_ref} is exact or
#'   empirical (see \code{design}/\code{measurement_prior_list} above):
#'   \itemize{
#'     \item \code{ref_source = "exact"}: only \code{Var(l)} (the numerator)
#'       has sampling error, so
#'       \eqn{\mathrm{Ratio}(l) \sim \chi^2_{n_{chains}-1} / (n_{chains}-1)}.
#'     \item \code{ref_source = "empirical (Var_final)"}: \code{Var_ref}
#'       itself is a sample variance/covariance across the same
#'       \code{n_chains} chains (just at the final sweep), so both
#'       numerator and denominator carry sampling error and the classic
#'       two-sample variance-ratio result applies instead:
#'       \eqn{\mathrm{Ratio}(l) \sim F(n_{chains}-1,\, n_{chains}-1)}. This
#'       is *wider* than the chi-squared band, and reduces to it exactly as
#'       \code{Var_final}'s own degrees of freedom go to infinity
#'       (\code{stats::qf(p, df1, Inf) == stats::qchisq(p, df1) / df1}), so
#'       both cases are computed via a single \code{stats::qf()} call with
#'       \code{df2 = Inf} for the exact case.
#'   }
#'   Not series-specific (a single band, not one per coefficient/eigenvalue)
#'   since it depends only on \code{n_chains}/\code{ref_source}. For
#'   \code{whitened = TRUE} this is an additional approximation regardless of
#'   \code{ref_source}: whitened eigenvalues are eigenvalues of a
#'   Wishart-distributed matrix, not marginally chi-squared/F, and tend to
#'   spread out more than this band implies (eigenvalue repulsion), so treat
#'   it as an illustrative ruler, not an exact whitened confidence band. Also
#'   note that \code{Var(l)} and \code{Var_final} come from the *same*
#'   \code{n_chains} chains at two different sweep indices, not independent
#'   samples the way a textbook two-sample F-test assumes -- another reason
#'   this is a naive/illustrative band rather than an exact one. Omit (the
#'   default, \code{NULL}) to skip the band entirely.
#' @param conf_level Confidence level for the \code{n_chains} band (default
#'   \code{0.95}); ignored when \code{n_chains} is \code{NULL}.
#' @param stage_label Character label for the title; defaults to \code{hist$stage}.
#' @return \code{hist} invisibly.
#' @details
#' \code{n_chains} is not stored on \code{hist} itself, so it must be
#' supplied from the fitted object that produced it: for the main-stage
#' history (\code{fit$sweep_history}) use \code{nrow(fit$fixef[[k]])} for
#' any RE component \code{k} (all have the same number of rows, the \code{n}
#' passed to the \verb{_run_with_pilot()} call); for the pilot-stage history
#' (\code{fit$pilot$sweep_history}) use \code{fit$pilot_chisq$n_pilot}.
#' @seealso \code{\link{plot_sweep_history_diag}}, \code{\link{lmerb_posterior_covariance}}
#' @export
plot_sweep_history_var_ratio <- function(
    hist,
    coef_focus = NULL,
    design = NULL,
    measurement_prior_list = NULL,
    whitened = FALSE,
    engine = c("base", "ggplot"),
    n_chains = NULL,
    conf_level = 0.95,
    stage_label = hist$stage
) {
  if (!inherits(hist, "two_block_sweep_history")) {
    stop(
      "'hist' must be a two_block_sweep_history object ",
      "(e.g. fit$sweep_history$main).",
      call. = FALSE
    )
  }
  engine <- match.arg(engine)
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) {
    stage_label <- if (!is.null(hist$stage)) hist$stage else "stage"
  }
  if (identical(engine, "ggplot") && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("'engine = \"ggplot\"' requires the ggplot2 package.", call. = FALSE)
  }

  sh_tab <- hist$table
  sh_tab$re_component <- as.character(sh_tab$re_component)
  sh_tab$covariate <- as.character(sh_tab$covariate)
  sh_sweeps <- sh_tab[sh_tab$sweep > 0L, , drop = FALSE]
  if (!nrow(sh_sweeps)) {
    warning("No sweep rows in sweep history for stage ", stage_label, call. = FALSE)
    return(invisible(hist))
  }

  exact_ref <- NULL
  ref_source <- "empirical (Var_final)"
  if (!is.null(design) && !is.null(measurement_prior_list)) {
    exact_ref <- lmerb_posterior_covariance(design, measurement_prior_list)
    ref_source <- "exact"
  }
  ## Exact reference => only Var(l) has sampling error (chi-sq/df); empirical
  ## Var_final reference => both numerator and denominator are sample
  ## variances across the same n_chains chains (two-sample F). The former is
  ## the df2 = Inf limit of the latter, so one qf() call covers both.
  band <- .sweep_var_ratio_naive_band(n_chains, conf_level, identical(ref_source, "exact"))

  if (isTRUE(whitened)) {
    var_ratio <- .sweep_var_ratio_whitened_series(hist, coef_focus, exact_ref)
    ylab <- "Whitened variance ratio (eigenvalue)"
    sub <- sprintf(
      "(whitened cross-chain covariance eigenvalues; reference = %s; Claim 3: I - A^(2l))",
      ref_source
    )
  } else {
    var_ratio <- .sweep_var_ratio_named_series(sh_sweeps, coef_focus, exact_ref)
    ylab <- if (identical(ref_source, "exact")) {
      "Var / Var_exact ratio"
    } else {
      "Var / Var_final ratio"
    }
    sub <- sprintf("(cross-chain variance from sweep history; reference = %s)", ref_source)
  }
  if (!is.null(band)) {
    ref_note <- if (identical(ref_source, "exact")) {
      "chi-sq"
    } else {
      "F, Var_ref also estimated"
    }
    sub <- paste0(
      sub,
      sprintf(
        "\ndotted: naive %.0f%% band around 1 (%s; n_chains = %d%s)",
        conf_level * 100, ref_note, n_chains,
        if (isTRUE(whitened)) "; approx. for eigenvalues" else ""
      )
    )
  }

  cat(sprintf("\n=== %s sweep history (%s; %s) ===\n\n", stage_label, ylab, engine))
  .sweep_var_ratio_render(var_ratio, stage_label, ylab, sub, engine, band)

  invisible(hist)
}

#' Naive band around 1 for a variance-ratio estimated from \code{n_chains}
#' draws, against either an exact or an empirical reference variance
#'
#' Under the (naive) null that the true ratio at a sweep already equals 1:
#' \itemize{
#'   \item \code{exact_ref = TRUE} (\code{Var_ref} is a fixed/known number,
#'     e.g. \code{lmerb_posterior_covariance()}'s \code{M^{-1}}): only the
#'     numerator \code{Var(l)} has sampling error, so
#'     \eqn{(n_{chains}-1) \cdot \mathrm{Ratio}(l) \sim \chi^2_{n_{chains}-1}},
#'     i.e. \code{Ratio(l)} itself is distributed as
#'     \eqn{\chi^2_{n_{chains}-1} / (n_{chains}-1)}.
#'   \item \code{exact_ref = FALSE} (\code{Var_ref} is the empirical
#'     \code{Var_final}, itself a sample variance across the same
#'     \code{n_chains} chains, just at the last sweep): both numerator and
#'     denominator are sample variances of the same true (co)variance, so
#'     the classic two-sample variance-ratio result applies instead:
#'     \eqn{\mathrm{Ratio}(l) \sim F(n_{chains}-1,\, n_{chains}-1)}.
#' }
#' These are computed via a single \code{stats::qf()} call: the exact case
#' is simply \code{df2 = Inf}, since
#' \code{stats::qf(p, df1, Inf) == stats::qchisq(p, df1) / df1} exactly (a
#' known limiting identity -- a fixed/known reference variance behaves like
#' a denominator sample variance with infinite degrees of freedom). Both are
#' exact for the named-coefficient (non-whitened) series, and only an
#' approximation for whitened eigenvalues (see
#' \code{\link{plot_sweep_history_var_ratio}}'s \code{n_chains} docs).
#' @noRd
.sweep_var_ratio_naive_band <- function(n_chains, conf_level, exact_ref) {
  if (is.null(n_chains)) {
    return(NULL)
  }
  n_chains <- as.integer(n_chains)[1L]
  if (is.na(n_chains) || n_chains <= 1L) {
    stop("'n_chains' must be a single integer > 1.", call. = FALSE)
  }
  conf_level <- as.numeric(conf_level)[1L]
  if (is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single value strictly between 0 and 1.", call. = FALSE)
  }
  df1 <- n_chains - 1L
  df2 <- if (isTRUE(exact_ref)) Inf else df1
  alpha <- 1 - conf_level
  stats::qf(c(alpha / 2, 1 - alpha / 2), df1 = df1, df2 = df2)
}

#' Resolve \code{coef_focus} into a (re_component, covariate) key data frame
#' @noRd
.sweep_var_ratio_resolve_keys <- function(coef_focus, available) {
  if (is.null(coef_focus)) {
    keys <- unique(available[, c("re_component", "covariate")])
    rownames(keys) <- NULL
    return(keys)
  }
  if (!is.list(coef_focus) || !length(coef_focus)) {
    stop("'coef_focus' must be a non-empty list of c(re_component, covariate) pairs.",
         call. = FALSE)
  }
  keys <- do.call(rbind, lapply(coef_focus, function(cc) {
    if (length(cc) < 2L) {
      stop("Each element of 'coef_focus' must be c(re_component, covariate).",
           call. = FALSE)
    }
    data.frame(
      re_component = as.character(cc[1L]),
      covariate    = as.character(cc[2L]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(keys) <- NULL
  keys
}

#' \code{"re_component | covariate"} labels for a key data frame
#' @noRd
.sweep_var_ratio_labels <- function(keys) {
  paste(keys$re_component, keys$covariate, sep = " | ")
}

#' Check that every label in \code{lbl} is a row/col name of \code{mat}
#' @noRd
.sweep_var_ratio_check_labels <- function(mat, lbl, what) {
  missing <- setdiff(lbl, rownames(mat))
  if (length(missing)) {
    stop(
      what, " has no entry for: ", paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Named-coefficient \code{Var(l) / Var_ref} series (non-whitened mode)
#' @noRd
.sweep_var_ratio_named_series <- function(sh_sweeps, coef_focus, exact_ref) {
  keys <- .sweep_var_ratio_resolve_keys(coef_focus, sh_sweeps)
  lbl  <- .sweep_var_ratio_labels(keys)
  if (!is.null(exact_ref)) {
    .sweep_var_ratio_check_labels(exact_ref, lbl, "Exact reference covariance")
  }

  var_list <- stats::setNames(
    lapply(seq_len(nrow(keys)), function(i) {
      sub <- sh_sweeps[
        sh_sweeps$re_component == keys$re_component[i] &
          sh_sweeps$covariate == keys$covariate[i],
        ,
        drop = FALSE
      ]
      if (!nrow(sub)) {
        stop("No sweep rows for ", lbl[i], call. = FALSE)
      }
      sub <- sub[order(sub$sweep), , drop = FALSE]
      stats::setNames(sub$sd^2, sub$sweep)
    }),
    lbl
  )

  var_ref <- stats::setNames(vapply(seq_along(var_list), function(i) {
    if (!is.null(exact_ref)) {
      exact_ref[lbl[[i]], lbl[[i]]]
    } else {
      v <- var_list[[i]]
      m_final <- max(as.integer(names(v)))
      as.numeric(v[as.character(m_final)])
    }
  }, numeric(1L)), lbl)

  stats::setNames(
    lapply(seq_along(var_list), function(i) var_list[[i]] / var_ref[[i]]),
    lbl
  )
}

#' Whitened \code{Var(l) / Var_ref} eigenvalue series (whitened mode)
#' @noRd
.sweep_var_ratio_whitened_series <- function(hist, coef_focus, exact_ref) {
  if (is.null(hist$cov_by_sweep) || !length(hist$cov_by_sweep)) {
    stop(
      "'hist' has no per-sweep covariance ('cov_by_sweep') -- whitened = TRUE ",
      "requires a sweep history built by a sweeps-outer/chains-inner engine ",
      "(e.g. rLMMNormal_reg_estimated_vcov(), rLMMindepNormalGamma_reg_known_vcov(), ",
      "rLMMindepNormalGamma_reg_estimated_vcov()).",
      call. = FALSE
    )
  }
  coef_index <- hist$coef_index
  keys <- .sweep_var_ratio_resolve_keys(coef_focus, coef_index)
  lbl  <- .sweep_var_ratio_labels(keys)

  n_sweeps <- length(hist$cov_by_sweep)
  .sweep_var_ratio_check_labels(hist$cov_by_sweep[[n_sweeps]], lbl, "Sweep-history covariance")

  Sigma_ref <- if (!is.null(exact_ref)) {
    .sweep_var_ratio_check_labels(exact_ref, lbl, "Exact reference covariance")
    exact_ref[lbl, lbl, drop = FALSE]
  } else {
    hist$cov_by_sweep[[n_sweeps]][lbl, lbl, drop = FALSE]
  }

  eg <- eigen(Sigma_ref, symmetric = TRUE)
  d_inv_sqrt <- ifelse(eg$values > .Machine$double.eps, 1 / sqrt(eg$values), 0)
  Sigma_ref_inv_sqrt <- eg$vectors %*% (d_inv_sqrt * t(eg$vectors))

  P <- length(lbl)
  eig_by_sweep <- matrix(NA_real_, nrow = n_sweeps, ncol = P)
  for (m in seq_len(n_sweeps)) {
    Cov_l <- hist$cov_by_sweep[[m]][lbl, lbl, drop = FALSE]
    W_l <- Sigma_ref_inv_sqrt %*% Cov_l %*% Sigma_ref_inv_sqrt
    W_l <- (W_l + t(W_l)) / 2
    eig_by_sweep[m, ] <- eigen(W_l, symmetric = TRUE, only.values = TRUE)$values
  }

  var_names <- paste0("var", seq_len(P))
  stats::setNames(
    lapply(seq_len(P), function(j) {
      stats::setNames(eig_by_sweep[, j], seq_len(n_sweeps))
    }),
    var_names
  )
}

#' Render one combined Var-ratio chart (base or ggplot)
#' @noRd
.sweep_var_ratio_render <- function(var_ratio, stage_label, ylab, sub, engine, band = NULL) {
  if (!length(var_ratio)) {
    warning("No variance-ratio series to plot for stage ", stage_label, call. = FALSE)
    return(invisible(NULL))
  }
  if (identical(engine, "base")) {
    .sweep_var_ratio_plot_base(var_ratio, stage_label, ylab, sub, band)
  } else {
    .sweep_var_ratio_plot_ggplot(var_ratio, stage_label, ylab, sub, band)
  }
}

#' Recycled palette of maximally-distinguishable point shapes
#'
#' \code{graphics}' \code{pch} and \code{ggplot2}'s \code{shape} aesthetic
#' both accept the same integer codes, so this one vector drives point
#' shapes in both \code{\link{.sweep_var_ratio_plot_base}} and
#' \code{\link{.sweep_var_ratio_plot_ggplot}}. Colour alone (via
#' \code{grDevices::hcl.colors()}) can be hard to tell apart across many
#' similarly-hued series -- shape gives a second, colour-blind-friendly
#' channel to distinguish traces/legend entries. Recycles (with a
#' \code{warning()}) if there are more series than shapes.
#' @noRd
.sweep_var_ratio_shapes <- function(n) {
  ## filled circle, triangle, square, diamond, plus, cross, asterisk, open
  ## circle, open triangle, open square, open diamond, inverted triangle.
  shapes <- c(16L, 17L, 15L, 18L, 3L, 4L, 8L, 1L, 2L, 0L, 5L, 6L)
  rep(shapes, length.out = n)
}

#' Estimate how many legend columns fit the current device width
#'
#' \code{graphics::legend(..., ncol = k)} does not auto-wrap based on
#' available width -- if the labels don't actually fit \code{k} across, the
#' extra columns silently spill off the left/right of the device/panel
#' instead of wrapping to more rows (this is exactly what happened with a
#' hardcoded \code{ncol = 6}: long \code{"re_component | covariate"} labels,
#' e.g. \code{"distracted_a1 | free_reduced_lunch"}, don't fit six-across on
#' an ordinary device, and the overflow columns ended up positioned outside
#' the plot area). Estimates the widest label's rendered width (in inches,
#' via \code{graphics::strwidth()}, so it reflects the actual font/cex/
#' device), adds a fixed allowance for the point symbol and inter-column
#' gap, and divides that into the current device width to get a column
#' count that actually fits.
#' @noRd
.sweep_var_ratio_legend_ncol <- function(labels, cex = 0.8, max_ncol = 6L) {
  label_w_in <- max(graphics::strwidth(labels, units = "inches", cex = cex))
  ## Allowance for the pch point plus legend()'s own internal symbol/text
  ## gap (a generous fixed estimate; legend() itself does not expose this).
  col_w_in   <- label_w_in + 0.35
  avail_w_in <- graphics::par("din")[1L] * 0.94
  ncol_fit   <- max(1L, floor(avail_w_in / col_w_in))
  max(1L, min(max_ncol, ncol_fit, length(labels)))
}

#' Base-R combined Var-ratio chart (all series on one panel + legend below)
#'
#' The legend lives in its own \code{\link[graphics]{layout}} panel below the
#' main plot, not hand-placed via \code{par("usr")} data-range arithmetic --
#' the previous approach positioned the legend at a fixed *fraction of the
#' plot's own y data range* below \code{usr[3]}, which (for any data range
#' larger than a couple of units) lands well outside the actual figure/device
#' area and is silently clipped, i.e. never drawn at all. A dedicated layout
#' panel is robust to the data range, the device size, and the number of
#' series (wraps to multiple legend rows via \code{ncol}).
#' @param band Optional length-2 numeric vector \code{c(lower, upper)} from
#'   \code{\link{.sweep_var_ratio_naive_band}}, drawn as a pair of horizontal
#'   dotted reference lines (not a per-series/legend entry).
#' @noRd
.sweep_var_ratio_plot_base <- function(var_ratio, stage_label, ylab, sub, band = NULL) {
  sweeps <- sort(unique(as.integer(unlist(lapply(var_ratio, names)))))
  if (!length(sweeps)) {
    warning("No sweep rows to plot for stage ", stage_label, call. = FALSE)
    return(invisible(NULL))
  }
  y_top <- max(1, unlist(var_ratio), band, na.rm = TRUE)
  if (!is.finite(y_top) || y_top <= 0) {
    y_top <- 1
  }
  ylim <- c(0, y_top * 1.05)

  n_series <- length(var_ratio)
  cols   <- grDevices::hcl.colors(n_series, palette = "Dark 3")
  shapes <- .sweep_var_ratio_shapes(n_series)
  y_first <- var_ratio[[1L]][match(sweeps, names(var_ratio[[1L]]))]

  has_legend <- n_series > 1L
  legend_ncol <- if (has_legend) {
    .sweep_var_ratio_legend_ncol(names(var_ratio), cex = 0.8)
  } else {
    1L
  }
  n_legend_rows <- if (has_legend) ceiling(n_series / legend_ncol) else 0L

  ## 'sub' may carry multiple '\n'-separated lines (e.g. the naive-band
  ## note appended below the base caption) -- mtext() does not honour
  ## embedded newlines the way text() does, so split and place each line at
  ## its own 'line=' offset, sizing the bottom margin to fit them all.
  sub_lines <- if (nzchar(sub)) strsplit(sub, "\n", fixed = TRUE)[[1L]] else character(0L)
  bottom_mar <- 3.2 + 1.3 * max(1L, length(sub_lines))

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1L)
    graphics::par(old_par)
  }, add = TRUE)

  if (has_legend) {
    ## Reserve a short bottom panel sized to the number of legend rows
    ## needed (each row ~0.09 of the plot panel's height), so the legend
    ## always has room regardless of how many series are being plotted.
    graphics::layout(matrix(1:2, nrow = 2L), heights = c(1, 0.12 + 0.09 * n_legend_rows))
    graphics::par(mar = c(bottom_mar, 4, 3, 1))
  } else {
    graphics::par(mar = c(bottom_mar, 4, 3, 1))
  }

  graphics::plot(
    sweeps, y_first,
    type = "n",
    xlab = "Inner sweep",
    ylab = ylab,
    ylim = ylim,
    main = paste0(stage_label, ": ", ylab, " vs sweep")
  )
  graphics::grid()
  graphics::abline(h = 1, lty = 2, col = "gray40")
  if (!is.null(band)) {
    graphics::abline(h = band, lty = 3, col = "gray50")
  }
  for (i in seq_len(n_series)) {
    v <- var_ratio[[i]]
    x <- as.integer(names(v))
    graphics::lines(x, v, type = "b", pch = shapes[i], col = cols[i])
  }
  for (i in seq_along(sub_lines)) {
    graphics::mtext(sub_lines[i], side = 1, line = 2.2 + 1.3 * i, cex = 0.8)
  }

  if (has_legend) {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::legend(
      "center", legend = names(var_ratio), col = cols, pch = shapes,
      ncol = legend_ncol, bty = "n", cex = 0.8
    )
  }
  invisible(var_ratio)
}

#' ggplot2 combined Var-ratio chart (all series on one panel)
#'
#' Shape (via \code{\link{.sweep_var_ratio_shapes}}, the same palette
#' \code{\link{.sweep_var_ratio_plot_base}} uses) is mapped alongside colour,
#' so series stay distinguishable even where the \code{hcl.colors()} hues
#' are close. \code{coef} is kept as a factor in \code{names(var_ratio)}'s
#' own order (not ggplot's default alphabetical order) so its levels line up
#' 1:1 with \code{shapes}.
#' @param band Optional length-2 numeric vector \code{c(lower, upper)} from
#'   \code{\link{.sweep_var_ratio_naive_band}}, drawn as a pair of horizontal
#'   dotted reference lines (not a per-series/legend entry).
#' @noRd
.sweep_var_ratio_plot_ggplot <- function(var_ratio, stage_label, ylab, sub, band = NULL) {
  coef_levels <- names(var_ratio)
  shapes <- .sweep_var_ratio_shapes(length(coef_levels))

  df <- do.call(rbind, lapply(coef_levels, function(nm) {
    v <- var_ratio[[nm]]
    data.frame(
      sweep = as.integer(names(v)), value = as.numeric(v), coef = nm,
      stringsAsFactors = FALSE
    )
  }))
  df$coef <- factor(df$coef, levels = coef_levels)

  ## Columns referenced via !!as.name(...) (injected symbols), not bare names,
  ## so aes() sees ordinary symbols resolved against 'df' at plot-build time --
  ## no "no visible binding for global variable" NOTE, and no need for a
  ## bare (unqualified) '.data' pronoun, which only resolves correctly when
  ## written literally inside aes() (see .sweep_var_ratio_plot_ggplot() docs).
  aes_sweep <- as.name("sweep")
  aes_value <- as.name("value")
  aes_coef  <- as.name("coef")
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = !!aes_sweep, y = !!aes_value,
      group = !!aes_coef, colour = !!aes_coef, shape = !!aes_coef
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "gray40") +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(
      x = "Inner sweep", y = ylab,
      title = paste0(stage_label, ": ", ylab, " vs sweep"),
      subtitle = sub, colour = NULL, shape = NULL
    )
  if (!is.null(band)) {
    p <- p + ggplot2::geom_hline(yintercept = band, linetype = "dotted", colour = "gray50")
  }
  print(p)
  invisible(var_ratio)
}
