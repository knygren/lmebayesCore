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
#' @param progbar Logical. When \code{TRUE}, show text progress bars over chains
#'   during each inner sweep. Default \code{FALSE}.
#' @param stage_label Character label stored on \code{$sweep_history} (e.g.
#'   \code{"pilot"} or \code{"main"}).
#' @param diag_sweeps When \code{TRUE}, print one combined Block~2 fixef
#'   chain-mean table for the stage when the inner-sweep loop finishes (same
#'   layout as \code{print()} on \code{$sweep_history} with all sweeps).
#'   Sweep history is always collected; use \code{print()} on
#'   \code{$sweep_history} to display tables later.
#' @param fixef_mode ICM mode reference stored on \code{$sweep_history}.
#' @param b_mode ICM mode reference for random-effect diagnostics.
#' @param b_start Initial random-effect matrix for all chains (\code{J x p_re}).
#'   Defaults to \code{b_mode} when \code{NULL}.
#' @param ptypes Per-component pfamily names (optional; derived from
#'   \code{pfamily_list} when \code{NULL}).
#' @param tau2_start Optional named numeric vector of plug-in
#'   \eqn{\tau^2_k} values for chain initialisation (one per \code{re_names}).
#'   When \code{NULL}, derived from \code{pfamily_list} prior fields
#'   (\code{dNormal} dispersion or ING \code{rate/shape} = \eqn{1/E[1/\tau^2]}).
#' @param use_cpp_block2 When \code{TRUE}, Block~2 uses
#'   \code{\link{two_block_block2_one_chain_cpp}} (native C++ align + \code{rglmb})
#'   instead of the pure-R reference.
#' @param use_cpp_block1 When \code{TRUE} (default), Block~1 uses an R loop over
#'   \code{two_block_block1_one_chain_cpp} (one \code{.Call} per chain). Set
#'   \code{use_cpp_block1_all_chains = TRUE} for a single all-chains \code{.Call}.
#'   When \code{use_cpp_block1 = FALSE}, use the R prep/draw loops (reference oracle).
#' @param use_cpp_block1_all_chains When \code{TRUE}, Block~1 uses
#'   \code{two_block_block1_all_chains_cpp_export} instead of the R chain loop.
#'   Default \code{FALSE}.
#' @param use_cpp_tau2_row Step A (\code{batch$tau2[i, ]}): \code{NULL} uses
#'   \code{getOption("glmbayesCore.use_cpp_tau2_row", TRUE)} (C++ export).
#'   Pass \code{FALSE} or set the option to \code{FALSE} for pure R row extract.
#' @param use_cpp_b_slice Step C (\code{batch$b[, , i] <- out$b}): \code{NULL} uses
#'   \code{getOption("glmbayesCore.use_cpp_b_slice", FALSE)} (pure R default).
#' @return A list with components \code{fixef_draws}, \code{dispersion_fixef_draws},
#'   \code{iters_fixef_draws}, \code{iters_ranef_draws}, \code{coefficients},
#'   \code{mu_all_last}, and \code{sweep_history} (class
#'   \code{"two_block_sweep_history"}).
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
    ptypes         = NULL,
    tau2_start     = NULL,
    use_cpp_block1 = TRUE,
    use_cpp_block1_all_chains = FALSE,
    use_cpp_tau2_row = NULL,
    use_cpp_b_slice = NULL,
    use_cpp_block2 = TRUE
) {
  if (is.null(ptypes)) {
    ptypes <- vapply(pfamily_list, function(pf) pf$pfamily, character(1))
    names(ptypes) <- re_names
  }

  if (is.null(tau2_start)) {
    tau2_start <- .two_block_tau2_start_from_pfamily(pfamily_list, re_names)
  } else {
    if (is.null(names(tau2_start)) || !setequal(names(tau2_start), re_names)) {
      stop("'tau2_start' must be a named vector with names(re_names).",
           call. = FALSE)
    }
    tau2_start <- as.numeric(tau2_start[re_names])
    names(tau2_start) <- re_names
  }
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

  progbar_use <- isTRUE(progbar)
  sweep_stats <- vector("list", inner_sweeps)

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
      progbar_finish_newline = FALSE,
      use_cpp_block1         = use_cpp_block1,
      use_cpp_block1_all_chains = use_cpp_block1_all_chains,
      use_cpp_tau2_row       = use_cpp_tau2_row,
      use_cpp_b_slice        = use_cpp_b_slice
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
      use_cpp_block2         = use_cpp_block2,
      progbar                = progbar_use,
      progbar_prefix         = prefix_b2,
      progbar_finish_newline = (m == inner_sweeps)
    )
    if (m == 1L) {
      cat("--- sweep 1 chain means (R batch) ---\n")
      for (k in re_names) {
        cat(sprintf("  fixef %s:", k), colMeans(batch$fixef[[k]]), "\n")
      }
      for (k in re_names) {
        cat(sprintf("  ranef %s:", k), mean(batch$b[, k, ]), "\n")
      }
      cat("\n")
    }
    # .two_block_print_sweep_boundary(
    #   stage_label  = stage_label,
    #   sweep        = m,
    #   inner_sweeps = inner_sweeps,
    #   phase        = "Block2",
    #   boundary     = "exit"
    # )
    sweep_stats[[m]] <- .two_block_snapshot_fixef_stats(batch, re_names)
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

  out <- .two_block_pack_batch_draws(
    batch          = batch,
    design         = design,
    collect_block1 = collect_block1
  )
  out$sweep_history <- .two_block_build_sweep_history(
    stage_label = stage_label,
    sweep_stats = sweep_stats,
    fixef_mode  = fixef_mode,
    re_names    = re_names
  )
  if (isTRUE(diag_sweeps)) {
    print(out$sweep_history)
  }
  invisible(out)
}
