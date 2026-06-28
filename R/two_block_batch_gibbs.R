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


#' Build a progress-bar prefix for sweep-outer Gibbs stages
#' @noRd
.two_block_progbar_prefix <- function(stage_label, sweep, inner_sweeps, phase) {
  phase_label <- if (identical(phase, "Block1")) {
    "RE"
  } else {
    "fixef"
  }
  stage_label <- as.character(stage_label)[1L]
  if (nzchar(stage_label)) {
    sprintf("[%s] sweep %d/%d %s: ", stage_label, sweep, inner_sweeps, phase_label)
  } else {
    sprintf("sweep %d/%d %s: ", sweep, inner_sweeps, phase_label)
  }
}

#' Text progress bar matching glmbayesCore C++ style
#' @param current Completed step (1-based, up to \code{total}).
#' @param total Total number of steps.
#' @param prefix Optional label printed before the bar (not cleared by \code{\r}).
#' @noRd
.two_block_progress_bar <- function(current, total, prefix = "") {
  total <- as.integer(total[1L])
  current <- as.integer(current[1L])
  if (!is.finite(total) || total <= 0L) {
    return(invisible())
  }
  if (!is.finite(current)) {
    current <- 0L
  }
  current <- max(0L, min(current, total))
  totaldotz <- 40L
  fraction  <- current / total
  if (!is.finite(fraction)) {
    return(invisible())
  }
  dotz      <- round(fraction * totaldotz)
  cat("\r", strrep(" ", 100L), "\r", sep = "")
  if (nzchar(prefix)) {
    cat(prefix)
  }
  cat(sprintf("%3.0f%% [", fraction * 100), sep = "")
  cat(paste0(rep("=", dotz), collapse = ""))
  cat(paste0(rep(" ", totaldotz - dotz), collapse = ""))
  cat("]", sep = "")
  utils::flush.console()
}

#' Finish a progress bar started by \code{.two_block_progress_bar}
#' @param newline If \code{FALSE}, leave the completed bar on the current
#'   line so the next bar can overwrite it with \code{\r}.
#' @noRd
.two_block_progress_bar_finish <- function(newline = TRUE) {
  if (isTRUE(newline)) {
    cat("\n")
  }
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

#' Envelope candidate count from a single-draw \code{rglmb} fit (\code{n = 1})
#' @noRd
.two_block_rglmb_iter_count <- function(fit_k) {
  if (is.null(fit_k$iters)) {
    return(1L)
  }
  it <- fit_k$iters
  if (is.matrix(it)) {
    return(as.integer(it[1L, 1L]))
  }
  as.integer(as.numeric(it)[1L])
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

  iters_ranef <- numeric(n_chains)

  list(
    n     = n_chains,
    fixef = fixef,
    tau2  = tau2,
    b     = b,
    iters = iters,
    iters_ranef = iters_ranef,
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

#' Refresh Block 1 prior precision for ING components (R reference; v5 any_ing)
#' @noRd
.two_block_block1_prior_with_tau2_r <- function(
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

#' Refresh Block 1 prior precision for ING components (mirrors C++ twoBlockGibbs)
#' @param use_cpp If \code{TRUE} (default), use the C++ implementation.
#' @noRd
.two_block_block1_prior_with_tau2 <- function(
    base_prior,
    tau2_vec,
    ptypes,
    re_names,
    mu_all,
    use_cpp = TRUE
) {
  if (isTRUE(use_cpp)) {
    return(two_block_block1_prior_with_tau2_cpp_export(
      base_prior = base_prior,
      tau2_vec   = tau2_vec,
      ptypes     = ptypes,
      re_names   = re_names,
      mu_all     = mu_all
    ))
  }
  .two_block_block1_prior_with_tau2_r(
    base_prior, tau2_vec, ptypes, re_names, mu_all
  )
}

#' One-chain Block 1 prep: fixef -> mu_all -> prior_list (no sampling)
#' @noRd
.two_block_block1_prep_one_chain <- function(
    batch,
    i,
    design,
    block1_prior,
    ptypes,
    use_cpp_mu_all = TRUE,
    use_cpp_prior_tau2 = TRUE
) {
  fixef_i <- .two_block_batch_fixef_chain(batch, i)
  mu_all  <- as.matrix(build_mu_all(
    design, fixef_i, batch$group_levels, use_cpp = use_cpp_mu_all
  )$mu_all)
  tau2_i  <- batch$tau2[i, ]
  prior_list <- .two_block_block1_prior_with_tau2(
    block1_prior, tau2_i, ptypes, batch$re_names, mu_all,
    use_cpp = use_cpp_prior_tau2
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
    progbar = FALSE,
    progbar_finish_newline = TRUE,
    use_cpp_mu_all = TRUE,
    use_cpp_prior_tau2 = TRUE
) {
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)

  prep_i <- function(i) {
    if (show_bar) .two_block_progress_bar(i, n)
    .two_block_block1_prep_one_chain(
      batch              = batch,
      i                  = i,
      design             = design,
      block1_prior       = block1_prior,
      ptypes             = ptypes,
      use_cpp_mu_all     = use_cpp_mu_all,
      use_cpp_prior_tau2 = use_cpp_prior_tau2
    )
  }

  prep_list <- .two_block_lapply_chains(n, prep_i, n_cores = n_cores)
  if (show_bar) .two_block_progress_bar_finish(newline = progbar_finish_newline)

  structure(
    list(
      mu_all      = lapply(prep_list, `[[`, "mu_all"),
      prior_lists = lapply(prep_list, `[[`, "prior_list")
    ),
    class = "two_block_block1_prep"
  )
}

#' Mean envelope candidates per group from a Block~1 draw (R reference)
#' @noRd
.two_block_block1_iters_mean_r <- function(block_out) {
  br <- block_out$block_results
  if (is.null(br) || !length(br)) {
    return(1)
  }
  vals <- vapply(br, function(b) {
    if (is.null(b$iters)) {
      return(1)
    }
    it <- b$iters
    if (is.matrix(it)) {
      return(as.numeric(it[1, 1]))
    }
    as.numeric(it[1])
  }, numeric(1))
  mean(vals)
}

#' Mean envelope candidates per group from a Block~1 draw
#' @param block_out Output from \code{block_rNormalGLM} or \code{block_rNormalReg}.
#' @param use_cpp If \code{TRUE} (default), use the C++ implementation.
#' @noRd
.two_block_block1_iters_mean <- function(block_out, use_cpp = TRUE) {
  if (isTRUE(use_cpp)) {
    return(as.numeric(two_block_block1_iters_mean_cpp_export(block_out)))
  }
  .two_block_block1_iters_mean_r(block_out)
}

#' Reorder Block~1 \code{b} rows to \code{group_levels} (R reference)
#' @noRd
.two_block_block1_reorder_b_r <- function(b_draw, group_levels) {
  rn <- rownames(b_draw)
  if (is.null(rn)) {
    return(b_draw)
  }
  ord <- match(group_levels, rn)
  if (any(is.na(ord))) {
    stop("Block 1 group ids do not match group_levels.", call. = FALSE)
  }
  b_draw[ord, , drop = FALSE]
}

#' Reorder Block~1 \code{b} rows to \code{group_levels}
#' @param b_draw Coefficient matrix from a Block~1 draw.
#' @param group_levels Target row order.
#' @param block_ids Row ids for \code{b_draw}; defaults to \code{rownames(b_draw)}.
#' @param use_cpp If \code{TRUE} (default), use the C++ implementation.
#' @noRd
.two_block_block1_reorder_b <- function(
    b_draw,
    group_levels,
    block_ids = rownames(b_draw),
    use_cpp = TRUE
) {
  if (isTRUE(use_cpp) && !is.null(block_ids)) {
    return(two_block_reorder_b_to_group_levels_cpp_export(
      b_draw, block_ids, group_levels
    ))
  }
  .two_block_block1_reorder_b_r(b_draw, group_levels)
}

#' One-chain Block 1 draw given a prepared prior_list
#' @noRd
.two_block_block1_draw_one_chain <- function(
    prior_list,
    design,
    family,
    is_gaussian,
    group_levels,
    use_cpp_reorder = TRUE,
    use_cpp_iters = TRUE
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
  b_draw <- .two_block_block1_reorder_b(
    b_draw        = b_draw,
    group_levels  = group_levels,
    use_cpp       = use_cpp_reorder
  )
  list(
    b          = b_draw,
    iters_mean = .two_block_block1_iters_mean(block_out, use_cpp = use_cpp_iters)
  )
}

#' All-chain Block 1 draw from prepared prior_lists
#' @noRd
.two_block_block1_draw_all_chains <- function(
    batch,
    prep,
    design,
    family,
    n_cores = NULL,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE,
    use_cpp_reorder = TRUE,
    use_cpp_iters = TRUE
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
    if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
    .two_block_block1_draw_one_chain(
      prior_list       = prior_lists[[i]],
      design           = design,
      family           = family,
      is_gaussian      = is_gaussian,
      group_levels     = batch$group_levels,
      use_cpp_reorder  = use_cpp_reorder,
      use_cpp_iters    = use_cpp_iters
    )
  }

  b_draws <- .two_block_lapply_chains(n, draw_i, n_cores = n_cores)
  if (show_bar) .two_block_progress_bar_finish(newline = progbar_finish_newline)

  for (i in seq_len(n)) {
    batch$b[, , i] <- b_draws[[i]]$b
    batch$iters_ranef[i] <- batch$iters_ranef[i] + b_draws[[i]]$iters_mean
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

#' Align random-effect vector to \code{X_hyper} row order
#'
#' Maps one column of Block~1 random effects (\code{b}, in \code{group_levels}
#' order) to the row order of \code{X_hyper[[k]]} for Block~2 pseudo-response
#' \code{y}. See \code{inst/ARCHITECTURE_glmerb.md}.
#'
#' @param b_vec Length-\code{J} vector for one RE component (named by
#'   \code{group_levels} or positional).
#' @param X_k Group-level design matrix (\code{J x q_k}).
#' @param group_levels Character vector defining Block~1 row order of \code{b}.
#' @return Numeric vector of length \code{nrow(X_k)} in \code{X_hyper} row order.
#' @export
two_block_align_b_to_xhyper <- function(b_vec, X_k, group_levels) {
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

#' @rdname two_block_align_b_to_xhyper
#' @export
two_block_align_b_to_xhyper_cpp <- function(b_vec, X_k, group_levels) {
  two_block_align_b_to_xhyper_cpp_export(b_vec, X_k, group_levels)
}

#' @keywords internal
#' @export
.two_block_align_b_to_xhyper <- two_block_align_b_to_xhyper

#' One-chain Block 2 update (writes \code{batch$fixef}, \code{batch$tau2},
#' \code{batch$iters})
#'
#' Given current random effects for replicate chain \code{i}, update fixed
#' effects (and ING dispersion) via one \code{rglmb()} call per RE component.
#'
#' @param batch Batch list from \code{.two_block_batch_init()}.
#' @param i Chain index (\code{1..batch$n}).
#' @param design Model design list (\code{X_hyper}, etc.).
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects.
#' @param ptypes Named character vector of \code{pfamily} types.
#' @return Updated \code{batch}.
#' @seealso \code{\link{two_block_block2_one_chain_cpp}},
#'   \code{\link{run_sweep_outer_chains_v6}}
#' @export
two_block_block2_one_chain <- function(
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
      n            = 1L,
      y            = y_k,
      x            = X_k,
      family       = stats::gaussian(),
      pfamily      = pf,
      verbose      = FALSE,
      use_parallel = FALSE
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
      it_k <- .two_block_rglmb_iter_count(fit_k)
      batch$iters[i, k] <- batch$iters[i, k] + it_k
    } else {
      batch$iters[i, k] <- batch$iters[i, k] + 1L
    }
  }
  batch
}

#' @keywords internal
#' @export
.two_block_block2_one_chain <- two_block_block2_one_chain

#' Block 2 one-chain update via C++ (native align + \code{rglmb})
#'
#' Same semantics as \code{\link{two_block_block2_one_chain}}; \code{b} is
#' aligned to \code{X_hyper} rows in C++ before each \code{rglmb()} call.
#'
#' @inheritParams two_block_block2_one_chain
#' @return Updated \code{batch}.
#' @export
two_block_block2_one_chain_cpp <- function(
    batch,
    i,
    design,
    pfamily_list,
    ptypes
) {
  b_i <- batch$b[, , i, drop = FALSE]
  b_i <- matrix(
    b_i, nrow = nrow(b_i), ncol = ncol(b_i),
    dimnames = dimnames(b_i)[1:2]
  )
  fixef_rows <- lapply(batch$re_names, function(k) batch$fixef[[k]][i, ])
  names(fixef_rows) <- batch$re_names
  x_hyper <- lapply(design$X_hyper, as.matrix)
  out <- two_block_block2_one_chain_cpp_export(
    b_i            = b_i,
    fixef_rows     = fixef_rows,
    tau2_i         = batch$tau2[i, ],
    iters_i        = batch$iters[i, ],
    x_hyper        = x_hyper,
    group_levels   = batch$group_levels,
    pfamily_list   = pfamily_list,
    ptypes         = ptypes,
    re_names       = batch$re_names
  )
  for (k in batch$re_names) {
    batch$fixef[[k]][i, ] <- out$fixef[[k]]
  }
  batch$tau2[i, ] <- out$tau2
  batch$iters[i, ] <- out$iters
  batch
}

#' Block 1 one-chain prep + draw (R reference; all piece flags FALSE)
#' @noRd
.two_block_block1_one_chain_r <- function(
    batch,
    i,
    design,
    block1_prior,
    family,
    ptypes
) {
  prep <- .two_block_block1_prep_one_chain(
    batch           = batch,
    i               = i,
    design          = design,
    block1_prior    = block1_prior,
    ptypes          = ptypes,
    use_cpp_mu_all  = FALSE,
    use_cpp_prior_tau2 = FALSE
  )
  draw <- .two_block_block1_draw_one_chain(
    prior_list      = prep$prior_list,
    design          = design,
    family          = family,
    is_gaussian     = identical(family$family, "gaussian"),
    group_levels    = batch$group_levels,
    use_cpp_reorder = FALSE,
    use_cpp_iters   = FALSE
  )
  c(prep, draw)
}

#' Shared \code{.Call} arguments for Block~1 C++ one-chain / all-chains exports
#' @noRd
.two_block_block1_cpp_call_args <- function(
    batch,
    design,
    block1_prior,
    family,
    ptypes
) {
  l2 <- length(design$y)
  offset <- design$offset
  if (is.null(offset)) {
    offset <- rep(0, l2)
  } else {
    offset <- as.numeric(offset)
    if (length(offset) == 1L) offset <- rep(offset, l2)
  }
  wt <- design$weights
  if (is.null(wt)) {
    wt <- rep(1, l2)
  } else {
    wt <- as.numeric(wt)
    if (length(wt) == 1L) wt <- rep(wt, l2)
  }
  is_gaussian <- identical(family$family, "gaussian")
  fam <- glmbfamfunc(family)
  fam_g <- glmbfamfunc(stats::gaussian())
  list(
    batch_fixef  = batch$fixef,
    batch_tau2   = batch$tau2,
    y            = as.numeric(design$y),
    Z            = as.matrix(design$Z),
    groups       = design$groups,
    offset       = offset,
    wt           = wt,
    x_hyper      = lapply(design$X_hyper, as.matrix),
    re_names     = batch$re_names,
    group_levels = batch$group_levels,
    ptypes       = ptypes,
    block1_prior = block1_prior,
    is_gaussian  = is_gaussian,
    f2           = fam$f2,
    f3           = fam$f3,
    f2_gauss     = fam_g$f2,
    f3_gauss     = fam_g$f3,
    family       = family$family,
    link         = family$link,
    Gridtype     = 2L,
    n_envopt     = 1L
  )
}

#' Block 1 one-chain prep + draw via C++ (native piecewise helpers)
#'
#' Same semantics as \code{.two_block_block1_prep_one_chain} followed by
#' \code{.two_block_block1_draw_one_chain} with all piecewise C++ flags
#' \code{TRUE} (\code{use_cpp_mu_all}, \code{use_cpp_prior_tau2},
#' \code{use_cpp_reorder}, \code{use_cpp_iters}).
#' Updates \code{batch$b[,, i]} and \code{batch$iters_ranef[i]}.
#'
#' @param batch Batch state from \code{.two_block_batch_init}.
#' @param i Chain index (\code{1..batch$n}).
#' @param design Model design list.
#' @param block1_prior Block~1 prior list.
#' @param family Response \code{family} object.
#' @param ptypes Named \code{pfamily} type vector.
#' @return Updated \code{batch} (invisibly).
#' @keywords internal
#' @export
two_block_block1_one_chain_cpp <- function(
    batch,
    i,
    design,
    block1_prior,
    family,
    ptypes
) {
  args <- .two_block_block1_cpp_call_args(
    batch, design, block1_prior, family, ptypes
  )
  args$batch_tau2 <- NULL
  out <- do.call(two_block_block1_one_chain_cpp_export, c(
    list(chain_i = as.integer(i), tau2_i = batch$tau2[i, ]),
    args
  ))

  batch$b[, , i] <- out$b
  batch$iters_ranef[i] <- batch$iters_ranef[i] + out$iters_mean
  batch
}

#' Block 1 batch: update random effects for all chains
#' @param use_cpp_block1 When \code{TRUE}, an R loop calls
#'   \code{\link{two_block_block1_one_chain_cpp}} (one C++ \code{.Call} per
#'   chain). When \code{FALSE}, use the R prep/draw loops with optional
#'   piecewise C++ helpers (\code{use_cpp_mu_all}, etc.). Default \code{TRUE}.
#' @noRd
.two_block_block1_all_chains <- function(
    batch,
    design,
    block1_prior,
    family,
    ptypes,
    n_cores = NULL,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE,
    use_cpp_block1 = TRUE,
    use_cpp_reorder = TRUE,
    use_cpp_iters = TRUE,
    use_cpp_mu_all = TRUE,
    use_cpp_prior_tau2 = TRUE
) {
  n <- batch$n

  if (isTRUE(use_cpp_block1)) {
    show_bar <- isTRUE(progbar) && n > 1L &&
      (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)
    for (i in seq_len(n)) {
      if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
      batch <- two_block_block1_one_chain_cpp(
        batch        = batch,
        i            = i,
        design       = design,
        block1_prior = block1_prior,
        family       = family,
        ptypes       = ptypes
      )
    }
    if (show_bar) .two_block_progress_bar_finish(newline = progbar_finish_newline)
    return(batch)
  }

  # .two_block_print_block1_phase("prep", "enter", n)
  prep <- .two_block_block1_prep_all_chains(
    batch                = batch,
    design               = design,
    block1_prior         = block1_prior,
    ptypes               = ptypes,
    n_cores              = n_cores,
    progbar              = FALSE,
    use_cpp_mu_all       = use_cpp_mu_all,
    use_cpp_prior_tau2   = use_cpp_prior_tau2
  )
  # .two_block_print_block1_phase("prep", "exit", n)
  # .two_block_print_block1_phase("draw", "enter", n)
  batch <- .two_block_block1_draw_all_chains(
    batch                   = batch,
    prep                    = prep,
    design                  = design,
    family                  = family,
    n_cores                 = n_cores,
    progbar                 = progbar,
    progbar_prefix          = progbar_prefix,
    progbar_finish_newline  = progbar_finish_newline,
    use_cpp_reorder         = use_cpp_reorder,
    use_cpp_iters           = use_cpp_iters
  )
  # .two_block_print_block1_phase("draw", "exit", n)
  batch
}

#' Block 2 batch: update fixed effects / tau2 for all chains
#' @noRd
.two_block_block2_all_chains <- function(
    batch,
    design,
    pfamily_list,
    ptypes,
    use_cpp_block2 = TRUE,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE
) {
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L
  block2_fn <- if (isTRUE(use_cpp_block2)) {
    two_block_block2_one_chain_cpp
  } else {
    two_block_block2_one_chain
  }

  for (i in seq_len(n)) {
    if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
    batch <- block2_fn(
      batch        = batch,
      i            = i,
      design       = design,
      pfamily_list = pfamily_list,
      ptypes       = ptypes
    )
  }
  if (show_bar) .two_block_progress_bar_finish(newline = progbar_finish_newline)
  batch
}

#' Starting tau2 vector from pfamily prior fields (plug-in dispersions)
#' @noRd
.two_block_tau2_start_from_pfamily <- function(pfamily_list, re_names) {
  vapply(re_names, function(k) {
    .two_block_tau2_ref_from_pfamily(pfamily_list[[k]])
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
    iters_ranef_draws      = batch$iters_ranef,
    coefficients           = coefficients,
    mu_all_last            = mu_all_last
  )
}
