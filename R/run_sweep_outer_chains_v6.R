#' Run independent short Gibbs chains via sweep-outer R driver
#'
#' Executes \code{n_chains} independent two-block Gibbs chains, each
#' initialised at \code{start_fixef} and run for \code{inner_sweeps} sweeps.
#' Loop order is sweep-outer: for each inner sweep, all chains run Block~1,
#' then all chains run Block~2. Each chain stores exactly one draw
#' (\code{n = 1L} per chain), so the result is a collection of
#' \code{n_chains} approximately independent draws from the target posterior
#' (for large enough \code{inner_sweeps}).
#'
#' @param n_chains Integer. Number of independent chains to run.
#' @param start_fixef Named list of starting hyper-parameter vectors, one
#'   named numeric vector per RE component (matches \code{re_names}).
#' @param inner_sweeps Integer. Number of inner Gibbs sweeps per chain.
#' @param design List with components \code{y}, \code{Z}, \code{groups},
#'   \code{X_hyper}, \code{re_coef_names}, and optional \code{group_name}.
#' @param block1_prior Block~1 prior list (\code{P} or \code{Sigma},
#'   \code{dispersion}, optional \code{ddef}).
#' @param pfamily_list Named list of \code{\link{pfamily}} objects (one per
#'   RE component).
#' @param family A \code{\link[stats]{family}} object for the response model.
#' @param re_names Character vector of RE coefficient names.
#' @param group_levels Character vector of group levels.
#' @param collect_block1 Logical. If \code{TRUE}, collect and rbind Block~1
#'   (\code{coefficients}) draws from every chain. Default \code{TRUE}.
#' @param progbar Logical. When \code{TRUE} (or when \code{diag_sweeps} is
#'   \code{TRUE}), show text progress bars over chains during each inner sweep.
#'   Default \code{FALSE}.
#' @param stage_label Character label for sweep diagnostics (e.g.
#'   \code{"pilot"}).
#' @param diag_sweeps If \code{TRUE}, print stage-end sweep history tables
#'   (fixef means and SDs by sweep) after the inner loop completes. Per-sweep
#'   boundary and fixef tables are suppressed by default.
#' @param fixef_mode ICM mode reference for fixef diagnostics.
#' @param b_mode ICM mode reference for random-effect diagnostics.
#' @param b_start Initial random-effect matrix for all chains (\code{J x p_re}).
#'   Defaults to \code{b_mode} when \code{NULL}.
#' @param ptypes Per-component pfamily names (optional; derived from
#'   \code{pfamily_list} when \code{NULL}).
#' @return A list with components \code{fixef_draws}, \code{dispersion_fixef_draws},
#'   \code{iters_fixef_draws}, \code{coefficients}, and \code{mu_all_last}.
#' @family simfuncs
#' @seealso \code{\link{two_block_rNormal_reg_v2}}, \code{\link{rGLMM}}
#' @export
run_sweep_outer_chains_v6 <- function(
    n_chains,
    start_fixef,
    inner_sweeps,
    design,
    block1_prior,
    pfamily_list,
    family,
    re_names,
    group_levels,
    collect_block1 = TRUE,
    progbar        = FALSE,
    stage_label    = "",
    diag_sweeps    = FALSE,
    fixef_mode     = NULL,
    b_mode         = NULL,
    b_start        = NULL,
    ptypes         = NULL
) {
  if (is.null(ptypes)) {
    ptypes <- vapply(pfamily_list, function(pf) pf$pfamily, character(1))
    names(ptypes) <- re_names
  }

  tau2_start <- .two_block_tau2_start_from_pfamily(pfamily_list, re_names)
  if (is.null(b_start)) {
    if (is.null(b_mode)) {
      stop("'b_start' or 'b_mode' required for batch init.", call. = FALSE)
    }
    b_start <- b_mode
  }

  batch <- .two_block_batch_init(
    n_chains     = n_chains,
    start_fixef  = start_fixef,
    b_start      = b_start,
    tau2_start   = tau2_start,
    re_names     = re_names,
    group_levels = group_levels
  )

  verbose_block_diag <- isTRUE(diag_sweeps)
  progbar_use <- isTRUE(progbar) || verbose_block_diag
  sweep_stats <- if (verbose_block_diag) {
    vector("list", inner_sweeps)
  } else {
    NULL
  }

  for (m in seq_len(inner_sweeps)) {
    prefix_b1 <- if (progbar_use) {
      .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block1")
    } else {
      ""
    }
    prefix_b2 <- if (progbar_use) {
      .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block2")
    } else {
      ""
    }
    # .two_block_print_sweep_boundary(
    #   stage_label  = stage_label,
    #   sweep        = m,
    #   inner_sweeps = inner_sweeps,
    #   phase        = "Block1",
    #   boundary     = "enter"
    # )
    batch <- .two_block_block1_all_chains(
      batch                  = batch,
      design                 = design,
      block1_prior           = block1_prior,
      family                 = family,
      ptypes                 = ptypes,
      progbar                = progbar_use,
      progbar_prefix         = prefix_b1,
      progbar_finish_newline = FALSE
    )
    # .two_block_print_sweep_boundary(
    #   stage_label  = stage_label,
    #   sweep        = m,
    #   inner_sweeps = inner_sweeps,
    #   phase        = "Block1",
    #   boundary     = "exit"
    # )
    # Fixef table only after Block 2 (gamma updated); skip after Block 1.
    # if (verbose_block_diag) {
    #   .two_block_print_block_diag(
    #     stage_label  = stage_label,
    #     sweep        = m,
    #     inner_sweeps = inner_sweeps,
    #     phase        = "Block1",
    #     batch        = batch,
    #     fixef_mode   = fixef_mode,
    #     b_mode       = b_mode,
    #     re_names     = re_names,
    #     group_levels = group_levels
    #   )
    # }

    # .two_block_print_sweep_boundary(
    #   stage_label  = stage_label,
    #   sweep        = m,
    #   inner_sweeps = inner_sweeps,
    #   phase        = "Block2",
    #   boundary     = "enter"
    # )
    batch <- .two_block_block2_all_chains(
      batch                  = batch,
      design                 = design,
      pfamily_list           = pfamily_list,
      ptypes                 = ptypes,
      progbar                = progbar_use,
      progbar_prefix         = prefix_b2,
      progbar_finish_newline = (m == inner_sweeps)
    )
    # .two_block_print_sweep_boundary(
    #   stage_label  = stage_label,
    #   sweep        = m,
    #   inner_sweeps = inner_sweeps,
    #   phase        = "Block2",
    #   boundary     = "exit"
    # )
    if (verbose_block_diag) {
      sweep_stats[[m]] <- .two_block_snapshot_fixef_stats(batch, re_names)
      # .two_block_print_block_diag(
      #   stage_label  = stage_label,
      #   sweep        = m,
      #   inner_sweeps = inner_sweeps,
      #   phase        = "Block2",
      #   batch        = batch,
      #   fixef_mode   = fixef_mode,
      #   b_mode       = b_mode,
      #   re_names     = re_names,
      #   group_levels = group_levels
      # )
    }
    if (progbar_use && n_chains <= 1L) {
      prefix_sweep <- if (nzchar(stage_label)) {
        sprintf("[%s] sweep %d/%d: ", stage_label, m, inner_sweeps)
      } else {
        sprintf("sweep %d/%d: ", m, inner_sweeps)
      }
      .two_block_progress_bar(m, inner_sweeps, prefix = prefix_sweep)
      .two_block_progress_bar_finish(newline = (m == inner_sweeps))
    }
  }

  if (verbose_block_diag) {
    .two_block_print_sweep_history_tables(
      stage_label = stage_label,
      sweep_stats = sweep_stats,
      fixef_mode  = fixef_mode,
      re_names    = re_names
    )
  }

  .two_block_pack_batch_draws(
    batch          = batch,
    design         = design,
    collect_block1 = collect_block1
  )
}
