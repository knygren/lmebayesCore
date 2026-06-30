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
#' @param fixef Named list of fixed-effect matrices (\code{n x p_k} per RE block).
#' @param re_names Character vector of random-effect block names.
#' @return Named list of \code{list(mean = ..., sd = ...)} per \code{re_names} entry.
#' @noRd
.two_block_snapshot_fixef_stats <- function(fixef, re_names) {
  out <- list()
  for (k in re_names) {
    mat <- fixef[[k]]
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

#' Initialize sweep state for \code{\link{rGLMM_sweep}}
#' @noRd
.rGLMM_sweep_initialize <- function(
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
    return(.two_block_block1_prior_with_tau2_cpp(
      base_prior, tau2_vec, ptypes, re_names, mu_all
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
    return(as.numeric(.two_block_block1_iters_mean_cpp(block_out)))
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
    return(.two_block_reorder_b_to_group_levels_cpp(
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
  .two_block_align_b_to_xhyper_cpp(b_vec, X_k, group_levels)
}

#' @keywords internal
#' @export
.two_block_align_b_to_xhyper <- two_block_align_b_to_xhyper

#' One-chain Block 2 update (fixed effects, \eqn{\tau^2}, iteration counts)
#'
#' Given current random effects for replicate chain \code{i}, update fixed
#' effects (and ING dispersion) via one \code{rglmb()} call per RE component.
#'
#' @param i Chain index (\code{1..n}).
#' @param b \code{J x p_re x n} array of random effects (read-only).
#' @param fixef Named list of fixed-effect matrices (\code{n x p_k} per RE block).
#' @param tau2 \code{n x p_re} matrix of random-effect variances.
#' @param iters \code{n x p_re} matrix of Block~2 iteration counts.
#' @param re_names Character vector of random-effect block names.
#' @param group_levels Character vector of group level labels.
#' @param design Model design list (\code{X_hyper}, etc.).
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects.
#' @param ptypes Named character vector of \code{pfamily} types.
#' @return List with \code{fixef}, \code{tau2}, and \code{iters} (each fully
#'   updated for chain \code{i}).
#' @seealso \code{\link{two_block_block2_one_chain_cpp}},
#'   \code{\link{rGLMM_sweep}}
#' @export
two_block_block2_one_chain <- function(
    i,
    b,
    fixef,
    tau2,
    iters,
    re_names,
    group_levels,
    design,
    pfamily_list,
    ptypes
) {
  b_i <- b[, , i, drop = FALSE]
  b_i <- matrix(b_i, nrow = nrow(b_i), ncol = ncol(b_i),
                dimnames = dimnames(b_i)[1:2])

  fixef_out <- fixef
  tau2_out <- tau2
  iters_out <- iters

  for (k in re_names) {
    X_k <- as.matrix(design$X_hyper[[k]])
    y_k <- .two_block_align_b_to_xhyper(
      b_vec        = b_i[, k],
      X_k          = X_k,
      group_levels = group_levels
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

    coef_k <- fit_k$coef.mode
    if (!is.null(names(coef_k))) {
      fixef_out[[k]][i, names(coef_k)] <- coef_k
    } else {
      fixef_out[[k]][i, ] <- coef_k
    }

    if (ptypes[[k]] == "dIndependent_Normal_Gamma") {
      tau2_out[i, k] <- fit_k$dispersion[1L]
      it_k <- .two_block_rglmb_iter_count(fit_k)
      iters_out[i, k] <- iters_out[i, k] + it_k
    } else {
      iters_out[i, k] <- iters_out[i, k] + 1L
    }
  }

  list(
    fixef = fixef_out,
    tau2  = tau2_out,
    iters = iters_out
  )
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
#' @return List with \code{fixef}, \code{tau2}, and \code{iters}.
#' @export
two_block_block2_one_chain_cpp <- function(
    i,
    b,
    fixef,
    tau2,
    iters,
    re_names,
    group_levels,
    design,
    pfamily_list,
    ptypes
) {
  b_i <- b[, , i, drop = FALSE]
  b_i <- matrix(
    b_i, nrow = nrow(b_i), ncol = ncol(b_i),
    dimnames = dimnames(b_i)[1:2]
  )
  fixef_rows <- lapply(re_names, function(k) fixef[[k]][i, ])
  names(fixef_rows) <- re_names
  x_hyper <- lapply(design$X_hyper, as.matrix)
  out <- .two_block_block2_one_chain_cpp(
    b_i            = b_i,
    fixef_rows     = fixef_rows,
    tau2_i         = tau2[i, ],
    iters_i        = iters[i, ],
    x_hyper        = x_hyper,
    group_levels   = group_levels,
    pfamily_list   = pfamily_list,
    ptypes         = ptypes,
    re_names       = re_names
  )

  fixef_out <- fixef
  for (k in re_names) {
    fixef_out[[k]][i, ] <- out$fixef[[k]]
  }
  tau2_out <- tau2
  tau2_out[i, ] <- out$tau2
  iters_out <- iters
  iters_out[i, ] <- out$iters

  list(
    fixef = fixef_out,
    tau2  = tau2_out,
    iters = iters_out
  )
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

#' Extract chain-i tau2 row from batch state (R reference)
#' @noRd
.two_block_batch_tau2_chain_row_r <- function(batch_tau2, chain_i) {
  batch_tau2[chain_i, , drop = TRUE]
}

#' Extract chain-i tau2 row from batch state (all-chains step A)
#'
#' Default \code{TRUE} on \code{\link{rGLMM_sweep}}; pass
#' \code{use_cpp_tau2_row = FALSE} to revert to \code{batch$tau2[i, ]} R indexing.
#'
#' @param batch_tau2 \code{n x p_re} matrix (\code{batch$tau2}).
#' @param chain_i Chain index (\code{1..n}).
#' @param use_cpp_tau2_row When \code{TRUE}, use C++ row extract; \code{FALSE}
#'   uses the R reference (\code{batch$tau2[i, ]}).
#' @noRd
.two_block_batch_tau2_chain_row <- function(
    batch_tau2,
    chain_i,
    use_cpp_tau2_row = TRUE
) {
  if (isTRUE(use_cpp_tau2_row)) {
    return(.two_block_batch_tau2_chain_row_cpp(
      batch_tau2,
      as.integer(chain_i)
    ))
  }
  .two_block_batch_tau2_chain_row_r(batch_tau2, chain_i)
}

#' Assign one chain slice of \code{batch$b} (R reference; all-chains step C)
#' @noRd
.two_block_batch_b_assign_slice_r <- function(b_store, chain_i, b_draw) {
  b_store[, , chain_i] <- b_draw
  b_store
}

#' Assign \code{batch$b[, , chain_i] <- b_draw} (all-chains step C)
#'
#' Default \code{TRUE} on \code{\link{rGLMM_sweep}}; pass
#' \code{use_cpp_b_slice = FALSE} to revert to R \code{batch$b[, , i] <- b_draw}.
#'
#' @param b_store 3-D array \code{batch$b} (\code{dim = c(J, p_re, n)}).
#' @param chain_i Chain index (\code{1..n}).
#' @param b_draw \code{J x p_re} matrix (\code{out$b} from one-chain export).
#' @param use_cpp_b_slice When \code{TRUE}, use C++ slice export; \code{FALSE}
#'   uses the R reference subassignment.
#' @noRd
.two_block_batch_b_assign_slice <- function(
    b_store,
    chain_i,
    b_draw,
    use_cpp_b_slice = FALSE
) {
  if (isTRUE(use_cpp_b_slice)) {
    return(.two_block_batch_b_assign_slice_cpp(
      b_store,
      as.integer(chain_i),
      b_draw
    ))
  }
  .two_block_batch_b_assign_slice_r(b_store, chain_i, b_draw)
}

#' Add one-chain Block~1 envelope iters to \code{batch$iters_ranef} (R reference; step D)
#' @noRd
.two_block_batch_iters_ranef_add_r <- function(iters_ranef, chain_i, iters_mean) {
  stopifnot(length(iters_ranef) >= chain_i)
  iters_ranef[chain_i] <- iters_ranef[chain_i] + as.numeric(iters_mean)
  iters_ranef
}

#' \code{batch$iters_ranef[chain_i] <- batch$iters_ranef[chain_i] + iters_mean}
#' (all-chains step D)
#'
#' Default \code{TRUE} on \code{\link{rGLMM_sweep}}; pass
#' \code{use_cpp_iters_ranef_add = FALSE} to revert to R accumulation.
#'
#' @param iters_ranef Length-\code{n} vector (\code{batch$iters_ranef}).
#' @param chain_i Chain index (\code{1..n}).
#' @param iters_mean Scalar from \code{out$iters_mean} (one-chain export).
#' @param use_cpp_iters_ranef_add When \code{TRUE}, use C++ add export; \code{FALSE}
#'   uses the R reference.
#' @noRd
.two_block_batch_iters_ranef_add <- function(
    iters_ranef,
    chain_i,
    iters_mean,
    use_cpp_iters_ranef_add = FALSE
) {
  if (isTRUE(use_cpp_iters_ranef_add)) {
    return(.two_block_batch_iters_ranef_add_cpp(
      iters_ranef,
      as.integer(chain_i),
      as.numeric(iters_mean)
    ))
  }
  .two_block_batch_iters_ranef_add_r(iters_ranef, chain_i, iters_mean)
}

#' Resolve Block~1 \code{f2}/\code{f3} closures once per batch (mirrors
#' \code{two_block_rNormal_reg}).
#' @noRd
.two_block_block1_glmbfamfunc <- function(family) {
  is_gaussian <- identical(family$family, "gaussian")
  famfunc_block1 <- glmbfamfunc(if (is_gaussian) gaussian() else family)
  famfunc_gauss <- glmbfamfunc(gaussian())
  list(
    f2       = famfunc_block1$f2,
    f3       = famfunc_block1$f3,
    f2_gauss = famfunc_gauss$f2,
    f3_gauss = famfunc_gauss$f3
  )
}

#' Block 1 batch: update random effects for all chains (C++ per-chain loop)
#'
#' Reads batch slice inputs, returns only Block~1 outputs (\code{b},
#' \code{iters_ranef}). Caller assigns \code{b} and \code{iters_ranef} into batch state.
#'
#' @param n Number of replicate chains.
#' @param fixef Named list of fixed-effect matrices (read-only for \code{mu_all}).
#' @param tau2 \code{n x p_re} matrix (read per chain).
#' @param b \code{J x p_re x n} array of random effects.
#' @param iters_ranef Length-\code{n} vector of random-effect iteration counts.
#' @param re_names Character vector of random-effect block names.
#' @param group_levels Character vector of group level labels.
#' @param use_cpp_tau2_row Step A (\code{tau2[i, ]}).
#' @param use_cpp_b_slice Step C (\code{b[, , i] <- out$b}).
#' @param use_cpp_iters_ranef_add Step D (\code{iters_ranef[i] += iters_mean}).
#' @return List with \code{b} and \code{iters_ranef}.
#' @noRd
.two_block_block1_all_chains <- function(
    n,
    fixef,
    tau2,
    b,
    iters_ranef,
    re_names,
    group_levels,
    design,
    block1_prior,
    family,
    ptypes,
    n_cores = NULL,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE,
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE
) {
  if (!is.null(n_cores) && as.integer(n_cores[1L]) >= 2L) {
    warning(
      "Chain-parallel Block 1 (n_cores > 1) is not supported; ",
      "running sequential loop.",
      call. = FALSE
    )
  }
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)

  fam_f23 <- .two_block_block1_glmbfamfunc(family)

  b_out <- b
  iters_ranef_out <- iters_ranef

  for (i in seq_len(n)) {
    if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
    chain_out <- .two_block_block1_one_chain_cpp(
      chain_i                 = i,
      b_store                 = b_out,
      iters_ranef             = iters_ranef_out,
      batch_fixef             = fixef,
      batch_tau2              = tau2,
      design                  = design,
      block1_prior            = block1_prior,
      family                  = family,
      ptypes                  = ptypes,
      re_names                = re_names,
      group_levels            = group_levels,
      f2                      = fam_f23$f2,
      f3                      = fam_f23$f3,
      f2_gauss                = fam_f23$f2_gauss,
      f3_gauss                = fam_f23$f3_gauss,
      use_cpp_tau2_row        = use_cpp_tau2_row,
      use_cpp_b_slice         = use_cpp_b_slice,
      use_cpp_iters_ranef_add = use_cpp_iters_ranef_add
    )
    b_out <- chain_out$b
    iters_ranef_out <- chain_out$iters_ranef
  }
  if (show_bar) {
    .two_block_progress_bar_finish(newline = progbar_finish_newline)
  }

  list(
    b           = b_out,
    iters_ranef = iters_ranef_out
  )
}

#' Block 1 batch: all chains via C++ loop (mirrors \code{.two_block_block1_all_chains})
#'
#' Single \code{.Call} entry; C++ loop mirrors \code{.two_block_block1_all_chains}.
#'
#' @inheritParams .two_block_block1_all_chains
#' @return List with \code{b} and \code{iters_ranef}.
#' @noRd
.two_block_block1_all_chains_via_cpp <- function(
    n,
    fixef,
    tau2,
    b,
    iters_ranef,
    re_names,
    group_levels,
    design,
    block1_prior,
    family,
    ptypes,
    n_cores = NULL,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE,
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE
) {
  if (!is.null(n_cores) && as.integer(n_cores[1L]) >= 2L) {
    warning(
      "Chain-parallel Block 1 (n_cores > 1) is not supported; ",
      "running sequential loop.",
      call. = FALSE
    )
  }
  fam_f23 <- .two_block_block1_glmbfamfunc(family)
  .two_block_block1_all_chains_cpp(
    b_store                 = b,
    iters_ranef             = iters_ranef,
    batch_fixef             = fixef,
    batch_tau2              = tau2,
    design                  = design,
    block1_prior            = block1_prior,
    family                  = family,
    ptypes                  = ptypes,
    re_names                = re_names,
    group_levels            = group_levels,
    f2                      = fam_f23$f2,
    f3                      = fam_f23$f3,
    f2_gauss                = fam_f23$f2_gauss,
    f3_gauss                = fam_f23$f3_gauss,
    use_cpp_tau2_row        = use_cpp_tau2_row,
    use_cpp_b_slice         = use_cpp_b_slice,
    use_cpp_iters_ranef_add = use_cpp_iters_ranef_add,
    progbar                 = progbar,
    progbar_prefix          = progbar_prefix,
    progbar_finish_newline  = progbar_finish_newline
  )
}

#' Block 2 batch: update fixed effects / tau2 for all chains
#'
#' Reads batch slice inputs, returns only Block~2 outputs (\code{fixef},
#' \code{tau2}, \code{iters}). Caller assigns \code{fixef}, \code{tau2}, and \code{iters}
#'   into batch state.
#'
#' @param n Number of replicate chains.
#' @param b \code{J x p_re x n} array of random effects (read-only).
#' @param fixef Named list of fixed-effect matrices.
#' @param tau2 \code{n x p_re} matrix of random-effect variances.
#' @param iters \code{n x p_re} matrix of Block~2 iteration counts.
#' @param re_names Character vector of random-effect block names.
#' @param group_levels Character vector of group level labels.
#' @param design Model design list.
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects.
#' @param ptypes Named character vector of \code{pfamily} types.
#' @param use_cpp_block2 When \code{TRUE}, use
#'   \code{\link{two_block_block2_one_chain_cpp}}; otherwise R reference.
#' @return List with \code{fixef}, \code{tau2}, and \code{iters}.
#' @noRd
.two_block_block2_all_chains <- function(
    n,
    b,
    fixef,
    tau2,
    iters,
    re_names,
    group_levels,
    design,
    pfamily_list,
    ptypes,
    use_cpp_block2 = TRUE,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE
) {
  show_bar <- isTRUE(progbar) && n > 1L
  block2_fn <- if (isTRUE(use_cpp_block2)) {
    two_block_block2_one_chain_cpp
  } else {
    two_block_block2_one_chain
  }

  fixef_out <- fixef
  tau2_out <- tau2
  iters_out <- iters

  for (i in seq_len(n)) {
    if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
    chain_out <- block2_fn(
      i            = i,
      b            = b,
      fixef        = fixef_out,
      tau2         = tau2_out,
      iters        = iters_out,
      re_names     = re_names,
      group_levels = group_levels,
      design       = design,
      pfamily_list = pfamily_list,
      ptypes       = ptypes
    )
    fixef_out <- chain_out$fixef
    tau2_out <- chain_out$tau2
    iters_out <- chain_out$iters
  }
  if (show_bar) .two_block_progress_bar_finish(newline = progbar_finish_newline)

  list(
    fixef = fixef_out,
    tau2  = tau2_out,
    iters = iters_out
  )
}

#' Starting tau2 vector from pfamily prior fields (plug-in dispersions)
#' @noRd
.two_block_tau2_start_from_pfamily <- function(pfamily_list, re_names) {
  vapply(re_names, function(k) {
    .two_block_tau2_ref_from_pfamily(pfamily_list[[k]])
  }, numeric(1))
}

#' Pack \code{rGLMM_sweep} outputs (draws, coefficients, \code{mu_all_last})
#' @param n Number of replicate chains.
#' @param fixef Named list of fixed-effect matrices.
#' @param tau2 \code{n x p_re} matrix of random-effect variances.
#' @param b \code{J x p_re x n} array of random effects.
#' @param iters \code{n x p_re} matrix of Block~2 iteration counts.
#' @param iters_ranef Length-\code{n} vector of random-effect iteration counts.
#' @param re_names Character vector of random-effect block names.
#' @param group_levels Character vector of group level labels.
#' @param design Model design list.
#' @param collect_block1 When \code{TRUE}, build long \code{coefficients} table from \code{b}.
#' @return List with \code{fixef_draws}, \code{dispersion_fixef_draws},
#'   \code{iters_fixef_draws}, \code{iters_ranef_draws}, \code{coefficients},
#'   and \code{mu_all_last}.
#' @noRd
.rGLMM_sweep_save <- function(
    n,
    fixef,
    tau2,
    b,
    iters,
    iters_ranef,
    re_names,
    group_levels,
    design,
    collect_block1 = TRUE
) {
  J <- length(group_levels)

  fixef_draws <- lapply(fixef, function(mat) {
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
        draw_df[[k]] <- b[, k, i]
      }
      coef_rows[[i]] <- draw_df
    }
    coefficients <- do.call(rbind, coef_rows)
    rownames(coefficients) <- NULL
  } else {
    coefficients <- NULL
  }

  fixef_mean <- lapply(fixef, colMeans)
  mu_all_last <- as.matrix(build_mu_all(
    design, fixef_mean, group_levels
  )$mu_all)

  list(
    fixef_draws            = fixef_draws,
    dispersion_fixef_draws = tau2,
    iters_fixef_draws      = iters,
    iters_ranef_draws      = iters_ranef,
    coefficients           = coefficients,
    mu_all_last            = mu_all_last
  )
}
