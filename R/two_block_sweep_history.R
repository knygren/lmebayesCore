#' Build structured two-block sweep history from per-sweep fixef snapshots
#' @param stage_label Character stage label (e.g. \code{"main"}).
#' @param sweep_stats Nested list from \code{.two_block_snapshot_fixef_stats}.
#' @param fixef_mode Named list of ICM mode vectors per RE component.
#' @param re_names Character vector of RE component names.
#' @param sweep_cov Optional list (length \code{length(sweep_stats)}) of
#'   \code{.two_block_snapshot_fixef_cov()} results, one per inner sweep. When
#'   supplied, stored as \code{cov_by_sweep} (list of covariance matrices) and
#'   \code{coef_index} (data frame of stacking order) on the returned object,
#'   enabling \code{\link{plot_var_convergence}}'s whitened mode.
#' @param disp_stats Optional list (length \code{length(sweep_stats)}) of
#'   \code{\link{.two_block_snapshot_disp_stats}} results, one per inner
#'   sweep. When supplied (and non-empty), stored as \code{disp_table} --
#'   same 5-column shape as \code{table} (\code{covariate} is always
#'   \code{"precision"}) -- enabling \code{\link{plot_var_convergence}}'s
#'   \code{component = "precision"} mode. \code{NULL}/empty (known-vcov
#'   engines, or zero ING Block~2 components) yields a 0-row \code{disp_table}.
#' @return Object of class \code{"two_block_sweep_history"}.
#' @noRd
.two_block_build_sweep_history <- function(
    stage_label,
    sweep_stats,
    fixef_mode,
    re_names,
    sweep_cov = NULL,
    disp_stats = NULL
) {
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) {
    stage_label <- "stage"
  }
  n_sweep <- length(sweep_stats)
  rows <- list()

  for (k in re_names) {
    cn <- if (n_sweep >= 1L) {
      names(sweep_stats[[1L]][[k]]$mean)
    } else {
      character(0L)
    }
    if (is.null(cn) || !length(cn)) {
      cn <- names(fixef_mode[[k]])
    }
    if (is.null(cn)) {
      cn <- character(0L)
    }
    for (nm in cn) {
      mode_v <- .two_block_fixef_mode_at(fixef_mode, k, nm, cn)
      rows[[length(rows) + 1L]] <- data.frame(
        re_component = k,
        covariate    = nm,
        sweep        = 0L,
        mean         = mode_v,
        sd           = NA_real_,
        stringsAsFactors = FALSE
      )
      for (m in seq_len(n_sweep)) {
        rows[[length(rows) + 1L]] <- data.frame(
          re_component = k,
          covariate    = nm,
          sweep        = as.integer(m),
          mean         = sweep_stats[[m]][[k]]$mean[[nm]],
          sd           = sweep_stats[[m]][[k]]$sd[[nm]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  table <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      re_component = character(),
      covariate    = character(),
      sweep        = integer(),
      mean         = numeric(),
      sd           = numeric(),
      stringsAsFactors = FALSE
    )
  }

  cov_by_sweep <- NULL
  coef_index   <- NULL
  if (!is.null(sweep_cov) && length(sweep_cov)) {
    cov_by_sweep <- lapply(sweep_cov, function(sc) sc$cov)
    coef_index   <- sweep_cov[[1L]]$coef_index
  }

  disp_ing_names <- if (!is.null(disp_stats) && length(disp_stats)) {
    names(disp_stats[[1L]])
  } else {
    character(0L)
  }
  disp_rows <- list()
  for (k in disp_ing_names) {
    for (m in seq_len(n_sweep)) {
      disp_rows[[length(disp_rows) + 1L]] <- data.frame(
        re_component = k,
        covariate    = "precision",
        sweep        = as.integer(m),
        mean         = disp_stats[[m]][[k]]$mean,
        sd           = disp_stats[[m]][[k]]$sd,
        stringsAsFactors = FALSE
      )
    }
  }
  disp_table <- if (length(disp_rows)) {
    do.call(rbind, disp_rows)
  } else {
    data.frame(
      re_component = character(),
      covariate    = character(),
      sweep        = integer(),
      mean         = numeric(),
      sd           = numeric(),
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      stage        = stage_label,
      n_sweeps     = as.integer(n_sweep),
      re_names     = re_names,
      fixef_mode   = fixef_mode,
      table        = table,
      cov_by_sweep = cov_by_sweep,
      coef_index   = coef_index,
      disp_table   = disp_table
    ),
    class = c("two_block_sweep_history", "list")
  )
}


#' Print stage-end table via structured sweep history (legacy helper)
#' @noRd
.two_block_print_sweep_history_tables <- function(
    stage_label,
    sweep_stats,
    fixef_mode,
    re_names,
    ...
) {
  hist <- .two_block_build_sweep_history(
    stage_label = stage_label,
    sweep_stats = sweep_stats,
    fixef_mode  = fixef_mode,
    re_names    = re_names
  )
  print(hist, ...)
}
