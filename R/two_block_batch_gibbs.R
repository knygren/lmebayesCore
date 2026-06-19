# Pure-R batch driver for two-block Gibbs sampling (sweep-outer layout).
#
# Two-block Gibbs alternates:
#   Block 1 — given hyperparameters (fixef gamma, tau2), sample random effects b
#             per neighborhood via block_rNormalGLM / block_rNormalReg.
#   Block 2 — given b, treat each RE column as a Gaussian pseudo-response and
#             update hyperparameters via rglmb (one call per RE component).
#
# Sweep-outer: ALL chains run Block 1, then ALL chains run Block 2, per inner
# sweep. Block 1 is split into prep (mu_all + prior_list per chain) and draw
# (block_rNormalGLM / block_rNormalReg per chain); both phases are embarrassingly
# parallel over chains (optional n_cores on Unix/macOS).


#' Text progress bar matching glmbayesCore C++ style
#' @param current Completed step (1-based, up to \code{total}).
#' @param total Total number of steps.
#' @noRd
.two_block_progress_bar <- function(current, total) {
  if (total <= 0L) {
    return(invisible())
  }
  totaldotz <- 40L
  fraction  <- current / total
  dotz      <- round(fraction * totaldotz)
  cat("\r", strrep(" ", 80L), "\r", sep = "")
  cat(sprintf("%3.0f%% [", fraction * 100), sep = "")
  cat(paste0(rep("=", dotz), collapse = ""))
  cat(paste0(rep(" ", totaldotz - dotz), collapse = ""))
  cat("]", sep = "")
  utils::flush.console()
}

#' Finish a progress bar started by \code{.two_block_progress_bar}
#' @noRd
.two_block_progress_bar_finish <- function() {
  cat("\n")
}

#' Print Block 1 prep/draw sub-phase boundary with wall-clock timestamp
#' @noRd
.two_block_print_block1_phase <- function(phase, boundary, n_chains) {
  phase <- as.character(phase)[1L]
  boundary <- as.character(boundary)[1L]
  action <- if (identical(boundary, "enter")) "Entering" else "Exiting"
  phase_label <- if (identical(phase, "prep")) {
    "Block1 prep (mu_all + prior_list)"
  } else if (identical(phase, "draw")) {
    "Block1 draw (block_rNormalGLM / block_rNormalReg)"
  } else {
    phase
  }
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf(
    "[Block1] %s %s, n=%d chains (%s)\n",
    action, phase_label, as.integer(n_chains)[1L], ts
  ))
  utils::flush.console()
  invisible(NULL)
}

#' Print sweep/block enter or exit line with wall-clock timestamp
#' @noRd
.two_block_print_sweep_boundary <- function(
    stage_label,
    sweep,
    inner_sweeps,
    phase,
    boundary
) {
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) stage_label <- "stage"
  phase_label <- if (identical(phase, "Block1")) {
    "random effects update"
  } else {
    "fixed effects update"
  }
  action <- if (identical(boundary, "enter")) "Entering" else "Exiting"
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf(
    "[%s] Sweep %d / %d: %s %s (%s)\n",
    stage_label, sweep, inner_sweeps, action, phase_label, ts
  ))
  utils::flush.console()
  invisible(NULL)
}

#' Lookup ICM mode for one fixef (re_name, covariate) pair
#' @noRd
.two_block_fixef_mode_at <- function(fixef_mode, re_name, cov_name, col_names) {
  if (is.null(fixef_mode) || is.null(fixef_mode[[re_name]])) {
    return(NA_real_)
  }
  mode_k <- fixef_mode[[re_name]]
  if (!is.null(names(mode_k)) && cov_name %in% names(mode_k)) {
    return(unname(mode_k[[cov_name]]))
  }
  j <- match(cov_name, col_names)
  if (is.finite(j) && length(mode_k) >= j) {
    return(unname(mode_k[j]))
  }
  NA_real_
}

#' Snapshot chain colMeans and SDs of Block 2 fixef after one sweep
#' @noRd
.two_block_snapshot_fixef_stats <- function(batch, re_names) {
  out <- list()
  for (k in re_names) {
    mat <- batch$fixef[[k]]
    cn  <- colnames(mat)
    if (is.null(cn)) {
      cn <- paste0("V", seq_len(ncol(mat)))
    }
    out[[k]] <- list(
      mean = stats::setNames(colMeans(mat), cn),
      sd   = stats::setNames(apply(mat, 2, stats::sd), cn)
    )
  }
  out
}

#' Print stage-end tables: fixef mode + per-sweep means; then per-sweep SDs
#' @noRd
.two_block_print_sweep_history_tables <- function(
    stage_label,
    sweep_stats,
    fixef_mode,
    re_names
) {
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) stage_label <- "stage"
  n_sweep <- length(sweep_stats)
  if (n_sweep < 1L) {
    return(invisible(NULL))
  }

  row_keys <- list()
  for (k in re_names) {
    cn <- names(sweep_stats[[1L]][[k]]$mean)
    if (is.null(cn)) {
      cn <- names(sweep_stats[[1L]][[k]]$sd)
    }
    for (nm in cn) {
      row_keys[[length(row_keys) + 1L]] <- list(re = k, cov = nm, cn = cn)
    }
  }

  sweep_hdr <- paste(
    vapply(seq_len(n_sweep), function(m) {
      sprintf("%12s", paste0("sweep ", m))
    }, character(1)),
    collapse = "  "
  )

  cat(sprintf(
    "\n--- two-block [%s stage summary: fixef means by sweep (%d sweeps)] ---\n",
    stage_label, n_sweep
  ))
  cat("  Block 2 fixed effects (chain colMeans after each sweep):\n")
  hdr_mean <- sprintf(
    "  %-18s  %-30s  %12s  %s",
    "Random effect", "Covariate", "mode", sweep_hdr
  )
  cat(hdr_mean, "\n")
  cat(paste0("  ", strrep("-", nchar(hdr_mean) - 2L)), "\n")

  for (row in row_keys) {
    k  <- row$re
    nm <- row$cov
    mode_v <- .two_block_fixef_mode_at(fixef_mode, k, nm, row$cn)
    mean_cols <- vapply(seq_len(n_sweep), function(m) {
      sweep_stats[[m]][[k]]$mean[[nm]]
    }, numeric(1))
    mean_fmt <- paste(sprintf("%12.4f", mean_cols), collapse = "  ")
    cat(sprintf(
      "  %-18s  %-30s  %12.4f  %s\n",
      k, nm, mode_v, mean_fmt
    ))
  }

  cat(sprintf(
    "\n--- two-block [%s stage summary: fixef sd by sweep (%d sweeps)] ---\n",
    stage_label, n_sweep
  ))
  cat("  Block 2 fixed effects (chain SD after each sweep):\n")
  hdr_sd <- sprintf(
    "  %-18s  %-30s  %s",
    "Random effect", "Covariate", sweep_hdr
  )
  cat(hdr_sd, "\n")
  cat(paste0("  ", strrep("-", nchar(hdr_sd) - 2L)), "\n")

  for (row in row_keys) {
    k  <- row$re
    nm <- row$cov
    sd_cols <- vapply(seq_len(n_sweep), function(m) {
      sweep_stats[[m]][[k]]$sd[[nm]]
    }, numeric(1))
    sd_fmt <- paste(sprintf("%12.4f", sd_cols), collapse = "  ")
    cat(sprintf("  %-18s  %-30s  %s\n", k, nm, sd_fmt))
  }
  cat("\n")
  invisible(NULL)
}

#' Print per-sweep block diagnostics (fixef table across chains; b vs mode optional)
#' @noRd
.two_block_print_block_diag <- function(
    stage_label,
    sweep,
    inner_sweeps,
    phase,
    batch,
    fixef_mode,
    b_mode,
    re_names,
    group_levels
) {
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) stage_label <- "stage"
  phase_note <- if (identical(phase, "Block1")) {
    " (fixef unchanged this phase)"
  } else {
    ""
  }
  cat(sprintf(
    "--- two-block [%s sweep %d / %d after %s, n=%d]%s ---\n",
    stage_label, sweep, inner_sweeps, phase, batch$n, phase_note
  ))

  hdr <- sprintf(
    "  %-18s  %-30s  %12s  %12s  %12s",
    "Random effect", "Covariate", "mode", "mean", "sd"
  )
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat("  Block 2 fixed effects (chain draws):\n")
  cat(hdr, "\n")
  cat(sep, "\n")

  for (k in re_names) {
    mat <- batch$fixef[[k]]
    cn  <- colnames(mat)
    if (is.null(cn)) {
      cn <- paste0("V", seq_len(ncol(mat)))
    }
    for (nm in cn) {
      fe_mean <- mean(mat[, nm])
      fe_sd   <- stats::sd(mat[, nm])
      fe_mode <- .two_block_fixef_mode_at(fixef_mode, k, nm, cn)
      cat(sprintf(
        "  %-18s  %-30s  %12.4f  %12.4f  %12.4f\n",
        k, nm, fe_mode, fe_mean, fe_sd
      ))
    }
  }

  # Random-effects vs mode (disabled; keep for debugging)
  # if (!is.null(b_mode)) {
  #   b_mean <- apply(batch$b, c(1, 2), mean)
  #   if (is.null(colnames(b_mean))) colnames(b_mean) <- re_names
  #   if (is.null(rownames(b_mean))) rownames(b_mean) <- group_levels
  #   cat(sprintf(
  #     "  b mean vs mode (%d group x RE coefficients):\n",
  #     length(b_mean)
  #   ))
  #   for (g in rownames(b_mean)) {
  #     for (k in colnames(b_mean)) {
  #       bm <- b_mean[g, k]
  #       md <- b_mode[g, k]
  #       cat(sprintf(
  #         "    %s::%-18s  mean %12.4f  mode %12.4f  delta %+.4f\n",
  #         g, k, bm, md, bm - md
  #       ))
  #     }
  #   }
  # }

  cat("\n")
  invisible(NULL)
}

#' Initialize batch state for sweep-outer R Gibbs driver
#' @noRd
.two_block_batch_init <- function(
    n_chains,
    start_fixef,
    b_start,
    tau2_start,
    re_names,
    group_levels
) {
  p_re <- length(re_names)
  J    <- length(group_levels)

  fixef <- lapply(start_fixef, function(beta0) {
    matrix(
      beta0,
      nrow = n_chains,
      ncol = length(beta0),
      byrow = TRUE,
      dimnames = list(NULL, names(beta0))
    )
  })
  names(fixef) <- re_names

  tau2 <- matrix(
    tau2_start,
    nrow = n_chains,
    ncol = p_re,
    byrow = TRUE,
    dimnames = list(NULL, re_names)
  )

  b <- array(b_start, dim = c(J, p_re, n_chains),
             dimnames = list(group_levels, re_names, NULL))

  iters <- matrix(0, nrow = n_chains, ncol = p_re)
  colnames(iters) <- re_names

  list(
    n     = n_chains,
    fixef = fixef,
    tau2  = tau2,
    b     = b,
    iters = iters,
    re_names     = re_names,
    group_levels = group_levels
  )
}

#' Extract chain-i fixef list from batch state
#' @noRd
.two_block_batch_fixef_chain <- function(batch, i) {
  lapply(batch$fixef, function(mat) {
    v <- mat[i, , drop = TRUE]
    stats::setNames(as.numeric(v), colnames(mat))
  })
}

#' Refresh Block 1 prior precision for ING components (mirrors C++ twoBlockGibbs)
#' @noRd
.two_block_block1_prior_with_tau2 <- function(
    base_prior,
    tau2_vec,
    ptypes,
    re_names,
    mu_all
) {
  out <- list(
    mu         = mu_all,
    dispersion = base_prior$dispersion,
    ddef       = base_prior$ddef
  )

  if (isTRUE(base_prior$ddef)) {
    out$P <- base_prior$P
    return(out)
  }

  if (!any(ptypes == "dIndependent_Normal_Gamma")) {
    out$P <- base_prior$P
    return(out)
  }

  P1 <- base_prior$P
  for (k in seq_along(re_names)) {
    if (ptypes[[k]] != "dIndependent_Normal_Gamma") next
    P1[k, ] <- 0
    P1[, k] <- 0
    P1[k, k] <- 1 / tau2_vec[k]
  }
  out$P <- P1
  out
}

#' One-chain Block 1 prep: fixef -> mu_all -> prior_list (no sampling)
#' @noRd
.two_block_block1_prep_one_chain <- function(
    batch,
    i,
    design,
    block1_prior,
    ptypes
) {
  fixef_i <- .two_block_batch_fixef_chain(batch, i)
  mu_all  <- as.matrix(build_mu_all(
    design, fixef_i, batch$group_levels
  )$mu_all)
  tau2_i  <- batch$tau2[i, ]
  prior_list <- .two_block_block1_prior_with_tau2(
    block1_prior, tau2_i, ptypes, batch$re_names, mu_all
  )
  list(mu_all = mu_all, prior_list = prior_list)
}

#' All-chain Block 1 prep: mu_all and prior_list for every chain
#' @noRd
.two_block_block1_prep_all_chains <- function(
    batch,
    design,
    block1_prior,
    ptypes,
    n_cores = NULL,
    progbar = FALSE
) {
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)

  prep_i <- function(i) {
    if (show_bar) .two_block_progress_bar(i, n)
    .two_block_block1_prep_one_chain(
      batch        = batch,
      i            = i,
      design       = design,
      block1_prior = block1_prior,
      ptypes       = ptypes
    )
  }

  prep_list <- .two_block_lapply_chains(n, prep_i, n_cores = n_cores)
  if (show_bar) .two_block_progress_bar_finish()

  structure(
    list(
      mu_all      = lapply(prep_list, `[[`, "mu_all"),
      prior_lists = lapply(prep_list, `[[`, "prior_list")
    ),
    class = "two_block_block1_prep"
  )
}

#' One-chain Block 1 draw given a prepared prior_list
#' @noRd
.two_block_block1_draw_one_chain <- function(
    prior_list,
    design,
    family,
    is_gaussian,
    group_levels
) {
  if (is_gaussian) {
    block_out <- block_rNormalReg(
      n          = 1L,
      y          = design$y,
      x          = design$Z,
      block      = design$groups,
      prior_list = prior_list
    )
  } else {
    block_out <- block_rNormalGLM(
      n            = 1L,
      y            = design$y,
      x            = design$Z,
      block        = design$groups,
      prior_list   = prior_list,
      family       = family,
      use_parallel = FALSE,
      verbose      = FALSE,
      progbar      = FALSE
    )
  }

  b_draw <- block_out$coefficients
  rn <- rownames(b_draw)
  if (!is.null(rn)) {
    ord <- match(group_levels, rn)
    if (any(is.na(ord))) {
      stop("Block 1 group ids do not match group_levels.", call. = FALSE)
    }
    b_draw <- b_draw[ord, , drop = FALSE]
  }
  b_draw
}

#' All-chain Block 1 draw from prepared prior_lists
#' @noRd
.two_block_block1_draw_all_chains <- function(
    batch,
    prep,
    design,
    family,
    n_cores = NULL,
    progbar = FALSE
) {
  is_gaussian <- identical(family$family, "gaussian")
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)
  prior_lists <- prep$prior_lists
  if (length(prior_lists) != n) {
    stop("length(prep$prior_lists) must equal batch$n.", call. = FALSE)
  }

  draw_i <- function(i) {
    if (show_bar) .two_block_progress_bar(i, n)
    .two_block_block1_draw_one_chain(
      prior_list   = prior_lists[[i]],
      design       = design,
      family       = family,
      is_gaussian  = is_gaussian,
      group_levels = batch$group_levels
    )
  }

  b_draws <- .two_block_lapply_chains(n, draw_i, n_cores = n_cores)
  if (show_bar) .two_block_progress_bar_finish()

  for (i in seq_len(n)) {
    batch$b[, , i] <- b_draws[[i]]
  }
  batch
}

#' Apply FUN to each chain index, optionally in parallel (Unix/macOS only)
#' @noRd
.two_block_lapply_chains <- function(n, FUN, n_cores = NULL) {
  idx <- seq_len(n)
  if (is.null(n_cores)) {
    return(lapply(idx, FUN))
  }
  n_cores <- as.integer(n_cores[1L])
  if (!is.finite(n_cores) || n_cores < 2L) {
    return(lapply(idx, FUN))
  }
  n_cores <- min(n_cores, n)
  if (.Platform$OS.type == "windows") {
    warning(
      "Chain-parallel Block 1 (n_cores > 1) is not supported on Windows; ",
      "using sequential lapply.",
      call. = FALSE
    )
    return(lapply(idx, FUN))
  }
  parallel::mclapply(idx, FUN, mc.cores = n_cores)
}

#' Align b vector to X_hyper row order (group levels)
#' @noRd
.two_block_align_b_to_xhyper <- function(b_vec, X_k, group_levels) {
  rn <- rownames(X_k)
  if (is.null(rn)) {
    if (length(b_vec) != nrow(X_k)) {
      stop(
        "length(b) (", length(b_vec), ") must equal nrow(X_hyper) (",
        nrow(X_k), ") when X_hyper has no rownames.",
        call. = FALSE
      )
    }
    return(b_vec)
  }
  if (!is.null(names(b_vec))) {
    miss <- setdiff(rn, names(b_vec))
    if (length(miss) > 0L) {
      stop(
        "Group level(s) missing from b: ", paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    return(unname(b_vec[rn]))
  }
  if (length(b_vec) != length(group_levels) ||
      length(b_vec) != nrow(X_k)) {
    stop(
      "b and X_hyper row counts disagree (b: ", length(b_vec),
      ", X_hyper: ", nrow(X_k), ", group_levels: ", length(group_levels), ").",
      call. = FALSE
    )
  }
  names(b_vec) <- group_levels
  miss <- setdiff(rn, group_levels)
  if (length(miss) > 0L) {
    stop(
      "X_hyper rownames do not match group_levels; missing in groups: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  b_vec[rn]
}

#' One-chain Block 2 update (writes batch$fixef, batch$tau2, batch$iters)
#' @noRd
.two_block_block2_one_chain <- function(
    batch,
    i,
    design,
    pfamily_list,
    ptypes
) {
  b_i <- batch$b[, , i, drop = FALSE]
  b_i <- matrix(b_i, nrow = nrow(b_i), ncol = ncol(b_i),
                dimnames = dimnames(b_i)[1:2])

  for (k in batch$re_names) {
    X_k <- as.matrix(design$X_hyper[[k]])
    y_k <- .two_block_align_b_to_xhyper(
      b_vec        = b_i[, k],
      X_k          = X_k,
      group_levels = batch$group_levels
    )
    pf  <- pfamily_list[[k]]

    fit_k <- rglmb(
      n       = 1L,
      y       = y_k,
      x       = X_k,
      family  = stats::gaussian(),
      pfamily = pf,
      verbose = FALSE
    )

    cn <- colnames(batch$fixef[[k]])
    coef_k <- fit_k$coef.mode
    if (!is.null(names(coef_k))) {
      batch$fixef[[k]][i, names(coef_k)] <- coef_k
    } else {
      batch$fixef[[k]][i, ] <- coef_k
    }

    if (ptypes[[k]] == "dIndependent_Normal_Gamma") {
      batch$tau2[i, k] <- fit_k$dispersion[1L]
      it_k <- if (!is.null(fit_k$iters)) fit_k$iters[1L, 1L] else 1L
      batch$iters[i, k] <- batch$iters[i, k] + it_k
    } else {
      batch$iters[i, k] <- batch$iters[i, k] + 1L
    }
  }
  batch
}

#' Block 1 batch: update random effects for all chains
#' @noRd
.two_block_block1_all_chains <- function(
    batch,
    design,
    block1_prior,
    family,
    ptypes,
    n_cores = NULL,
    progbar = FALSE
) {
  n <- batch$n
  .two_block_print_block1_phase("prep", "enter", n)
  prep <- .two_block_block1_prep_all_chains(
    batch        = batch,
    design       = design,
    block1_prior = block1_prior,
    ptypes       = ptypes,
    n_cores      = n_cores,
    progbar      = FALSE
  )
  .two_block_print_block1_phase("prep", "exit", n)
  .two_block_print_block1_phase("draw", "enter", n)
  batch <- .two_block_block1_draw_all_chains(
    batch   = batch,
    prep    = prep,
    design  = design,
    family  = family,
    n_cores = n_cores,
    progbar = progbar
  )
  .two_block_print_block1_phase("draw", "exit", n)
  batch
}

#' Block 2 batch: update fixed effects / tau2 for all chains
#' @noRd
.two_block_block2_all_chains <- function(
    batch,
    design,
    pfamily_list,
    ptypes,
    progbar = FALSE
) {
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L

  for (i in seq_len(n)) {
    if (show_bar) .two_block_progress_bar(i, n)
    batch <- .two_block_block2_one_chain(
      batch        = batch,
      i            = i,
      design       = design,
      pfamily_list = pfamily_list,
      ptypes       = ptypes
    )
  }
  if (show_bar) .two_block_progress_bar_finish()
  batch
}

#' Starting tau2 vector from pfamily_list (plug-in dispersions)
#' @noRd
.two_block_tau2_start_from_pfamily <- function(pfamily_list, re_names) {
  vapply(re_names, function(k) {
    pl <- pfamily_list[[k]]$prior_list
    pf <- pfamily_list[[k]]$pfamily
    if (pf == "dNormal") {
      pl$dispersion
    } else if (pf == "dIndependent_Normal_Gamma") {
      pl$disp_lower
    } else {
      stop("Unsupported pfamily: ", pf, call. = FALSE)
    }
  }, numeric(1))
}

#' Pack batch state into replicate-chain draw list
#' @noRd
.two_block_pack_batch_draws <- function(
    batch,
    design,
    collect_block1 = TRUE
) {
  re_names     <- batch$re_names
  group_levels <- batch$group_levels
  n            <- batch$n
  J            <- length(group_levels)

  fixef_draws <- lapply(batch$fixef, function(mat) {
    out <- mat
    rownames(out) <- NULL
    out
  })
  names(fixef_draws) <- re_names

  if (isTRUE(collect_block1)) {
    grp_col <- design$group_name
    if (is.null(grp_col) || !nzchar(grp_col)) grp_col <- "group"
    coef_rows <- vector("list", n)
    for (i in seq_len(n)) {
      draw_df <- data.frame(
        draw = rep(i, J),
        stringsAsFactors = FALSE
      )
      draw_df[[grp_col]] <- group_levels
      for (k in re_names) {
        draw_df[[k]] <- batch$b[, k, i]
      }
      coef_rows[[i]] <- draw_df
    }
    coefficients <- do.call(rbind, coef_rows)
    rownames(coefficients) <- NULL
  } else {
    coefficients <- NULL
  }

  fixef_mean <- lapply(batch$fixef, colMeans)
  mu_all_last <- as.matrix(build_mu_all(
    design, fixef_mean, group_levels
  )$mu_all)

  list(
    fixef_draws            = fixef_draws,
    dispersion_fixef_draws = batch$tau2,
    iters_fixef_draws      = batch$iters,
    coefficients           = coefficients,
    mu_all_last            = mu_all_last
  )
}
