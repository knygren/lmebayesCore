## DEV ONLY -- does not migrate to glmbayesCore.
##
## Replay (tape), input/draw log (record), and groups prior-mean tracking for
## the indexed two-block sweep. All of it attaches to the production loops
## (.two_block_sweep_r and glmbayes::sim::twoblock::two_block_sweep) through
## their observer hooks, so the loops being checked are the loops that run.
## C++ twin: src/dev_two_block_replay.cpp (lmebayes::dev::twoblock).

#' Empty input and draw log for an indexed sweep
#'
#' A filled log is also a valid \code{tape}: replaying it reproduces the run.
#' @return List with \code{mu_D} (\code{ncol(D) x J x n x M}), \code{block1}
#'   (\code{J x ncol(D) x n x M}), \code{iters1} (\code{n x M}), and lists
#'   \code{y} (\code{nrow(W_k) x n x M}) and \code{block2}
#'   (\code{ncol(W_k) x n x M}) over components.
#' @noRd
.dev_two_block_log_new <- function(spec, n, M) {
  l1 <- ncol(spec$D)
  J <- length(spec$obs_rows)
  list(
    mu_D   = array(0, c(l1, J, n, M)),
    block1 = array(0, c(J, l1, n, M)),
    iters1 = matrix(0, n, M),
    y      = lapply(spec$W, function(Wk) array(0, c(nrow(Wk), n, M))),
    block2 = lapply(spec$W, function(Wk) array(0, c(ncol(Wk), n, M)))
  )
}

#' Observer for .two_block_sweep_r: replay, record, and prior sums
#'
#' R twin of \code{lmebayes::dev::twoblock::ReplayObserver}.
#' @param tape \code{NULL}, or a log in \code{.dev_two_block_log_new()}
#'   layout; its \code{block1}, \code{iters1}, and \code{block2} entries
#'   replace the samplers.
#' @return Observer list for \code{.two_block_sweep_r(observer = )}, plus
#'   \code{env} holding \code{log}, \code{prior_sum}, and \code{prior_sumsq}.
#' @noRd
.dev_two_block_observer <- function(spec, n, M, tape = NULL, record = FALSE) {
  l1 <- ncol(spec$D)
  J <- length(spec$obs_rows)
  env <- new.env(parent = emptyenv())
  env$log <- if (isTRUE(record)) .dev_two_block_log_new(spec, n, M) else NULL
  env$prior_sum <- array(0, c(l1, J, M))
  env$prior_sumsq <- array(0, c(l1, J, M))
  list(
    groups_source = function(m, i) {
      if (is.null(tape)) {
        return(NULL)
      }
      list(coef  = matrix(tape$block1[, , i, m], nrow = J, ncol = l1),
           iters = tape$iters1[i, m])
    },
    population_source = function(m, i, k) {
      if (is.null(tape)) NULL else tape$block2[[k]][, i, m]
    },
    groups_seen = function(m, i, mu_D, coef, iters) {
      env$prior_sum[, , m] <- env$prior_sum[, , m] + mu_D
      env$prior_sumsq[, , m] <- env$prior_sumsq[, , m] + mu_D * mu_D
      if (!is.null(env$log)) {
        env$log$mu_D[, , i, m] <- mu_D
        env$log$block1[, , i, m] <- coef
        env$log$iters1[i, m] <- iters
      }
    },
    population_seen = function(m, i, k, y_k, coef) {
      if (!is.null(env$log)) {
        env$log$y[[k]][, i, m] <- y_k
        env$log$block2[[k]][, i, m] <- coef
      }
    },
    env = env
  )
}

#' Indexed sweep with replay, record, and prior sums
#'
#' @param engine \code{"r"} (\code{.two_block_sweep_r}) or \code{"cpp"}
#'   (\code{two_block_sweep} through \code{.dev_two_block_sweep_cpp}).
#' @return State list with \code{b}, \code{fixef}, \code{iters},
#'   \code{iters_ranef}, \code{log} (\code{NULL} unless \code{record}), and
#'   \code{prior_sum}, \code{prior_sumsq}: \code{ncol(D) x J x inner_sweeps}
#'   sums over chains of the groups prior mean and of its square.
#' @noRd
.dev_two_block_sweep_replay <- function(spec, state, inner_sweeps,
                                        use_mode = FALSE, progbar = FALSE,
                                        stage_label = "", tape = NULL,
                                        record = FALSE,
                                        engine = c("r", "cpp")) {
  engine <- match.arg(engine)
  inner_sweeps <- as.integer(inner_sweeps)
  stage_label <- if (is.null(stage_label)) "" else as.character(stage_label)[1L]
  if (identical(engine, "cpp")) {
    return(.dev_two_block_sweep_cpp(
      spec, state$b, state$fixef, state$iters, state$iters_ranef,
      inner_sweeps, isTRUE(use_mode), isTRUE(progbar), stage_label,
      tape, isTRUE(record)
    ))
  }
  obs <- .dev_two_block_observer(spec, dim(state$b)[3L], inner_sweeps,
                                 tape = tape, record = record)
  out <- .two_block_sweep_r(spec, state, inner_sweeps, use_mode = use_mode,
                            progbar = progbar, stage_label = stage_label,
                            observer = obs)
  out["log"] <- list(obs$env$log)
  out$prior_sum <- obs$env$prior_sum
  out$prior_sumsq <- obs$env$prior_sumsq
  out
}

#' Groups conditional prior fed to each group, summarized over chains
#'
#' @param state Result of \code{.dev_two_block_sweep_replay()}.
#' @return List with \code{prior_mean} and \code{prior_mean_var}
#'   (\code{p_re x J x inner_sweeps}: mean and variance over chains of the
#'   prior mean fed to each group in each sweep), \code{prior_precision}
#'   (\code{p_re x p_re}, the same for every group, chain, and sweep), and
#'   \code{prior_var} (named vector \code{diag(solve(prior_precision))}).
#' @noRd
.dev_two_block_prior_track <- function(state, spec, re_names, group_levels) {
  n <- dim(state$b)[3L]
  s <- state$prior_sum[spec$d_col, , , drop = FALSE]
  q <- state$prior_sumsq[spec$d_col, , , drop = FALSE]
  prior_mean <- s / n
  prior_mean_var <- if (n > 1L) {
    (q - n * prior_mean * prior_mean) / (n - 1L)
  } else {
    array(NA_real_, dim(s))
  }
  dimnames(prior_mean) <- dimnames(prior_mean_var) <-
    list(re_names, group_levels, NULL)
  P <- spec$P_D[spec$d_col, spec$d_col, drop = FALSE]
  dimnames(P) <- list(re_names, re_names)
  list(
    prior_mean      = prior_mean,
    prior_mean_var  = prior_mean_var,
    prior_precision = P,
    prior_var       = stats::setNames(diag(solve(P)), re_names)
  )
}

#' Shared body of the tracked indexed sweeps
#' @return \code{.two_block_finish()} output plus \code{prior_track}.
#' @noRd
.dev_rGLMM_sweep_tracked <- function(engine, n_chains, start_fixef,
                                     inner_sweeps, design, block1_prior,
                                     pfamily_list, family, re_names,
                                     group_levels, collect_block1, progbar,
                                     stage_label, b_mode, b_start,
                                     tau2_start, index_map,
                                     fixef_mode = NULL, diag_sweeps = FALSE) {
  prep <- .two_block_prepare(
    n_chains, start_fixef, design, block1_prior, pfamily_list, family,
    re_names, group_levels, b_start, b_mode, tau2_start, index_map
  )
  state <- .dev_two_block_sweep_replay(
    spec         = prep$spec,
    state        = prep$state,
    inner_sweeps = inner_sweeps,
    progbar      = progbar,
    stage_label  = stage_label,
    engine       = engine
  )
  out <- .two_block_finish(
    state, prep$tau2, index_map, design, re_names, collect_block1,
    stage_label = stage_label, fixef_mode = fixef_mode,
    diag_sweeps = diag_sweeps
  )
  out$prior_track <- .dev_two_block_prior_track(
    state, prep$spec, re_names, index_map$group_levels
  )
  out
}

#' \code{.rGLMM_sweep_indexed} with groups prior tracking
#'
#' Same arguments and draws as \code{.rGLMM_sweep_indexed}; the output also
#' has \code{prior_track} (see \code{.dev_two_block_prior_track()}).
#' @noRd
.dev_rGLMM_sweep_indexed_tracked <- function(
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
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE,
    use_cpp_block2 = TRUE,
    rwt = 0.01,
    index_map
) {
  .dev_rGLMM_sweep_tracked(
    "r", n_chains, start_fixef, inner_sweeps, design, block1_prior,
    pfamily_list, family, re_names, group_levels, collect_block1, progbar,
    stage_label, b_mode, b_start, tau2_start, index_map,
    fixef_mode = fixef_mode, diag_sweeps = diag_sweeps
  )
}

#' \code{.rGLMM_sweep_cpp} with groups prior tracking
#'
#' Same arguments and draws as \code{.rGLMM_sweep_cpp}; the output also has
#' \code{prior_track} (see \code{.dev_two_block_prior_track()}).
#' @noRd
.dev_rGLMM_sweep_cpp_tracked <- function(
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
    use_cpp_tau2_row = TRUE,
    use_cpp_b_slice = TRUE,
    use_cpp_iters_ranef_add = TRUE,
    use_cpp_block2 = TRUE,
    rwt = 0.01,
    index_map
) {
  .dev_rGLMM_sweep_tracked(
    "cpp", n_chains, start_fixef, inner_sweeps, design, block1_prior,
    pfamily_list, family, re_names, group_levels, collect_block1, progbar,
    stage_label, b_mode, b_start, tau2_start, index_map,
    fixef_mode = fixef_mode, diag_sweeps = diag_sweeps
  )
}
