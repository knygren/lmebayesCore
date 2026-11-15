## Indexed two-block sweep: the R reference loop and the R side of the C++
## loop in src/two_block_sweep.cpp (glmbayes::sim::twoblock::two_block_sweep).
##
## Names are resolved once, in .two_block_spec() and .two_block_init().
## After that every step works on plain numeric matrices, numeric vectors,
## and 1-based integer index vectors, and each step has a C++ twin with the
## same inputs and the same output shape:
##
##   R                                  C++ (glmbayes::sim::twoblock)
##   .two_block_groups_prior_mean       groups_prior_mean
##   .two_block_groups_draw             groups_draw
##   .two_block_groups_place            groups_place
##   .two_block_population_response     population_response
##   .two_block_population_draw         population_draw
##   .two_block_population_write        population_write
##   .two_block_sweep_r                 two_block_sweep
##
## The draws in both languages come from the same compiled samplers:
## rNormalGLMGroups (groups) and rNormalReg (population). Names are restored
## afterwards by .rGLMM_sweep_save_indexed().
##
## State layout, identical in both languages:
##   b            numeric array  J x p_re x n
##   fixef        list of p_re numeric matrices, n x p_k, stored column order
##   iters        numeric matrix n x p_re
##   iters_ranef  numeric vector n

#' Plain-number inputs shared by every sweep step
#'
#' @param design List with \code{y}, \code{D}, and \code{W}.
#' @param block1_prior Groups prior with \code{P} (\code{p_re x p_re},
#'   \code{re_names} order) and optional \code{dispersion}.
#' @param pfamily_list Named list of \code{dNormal} pfamily objects.
#' @param family GLM family. Gaussian is not supported.
#' @param re_names Random-effect names, in stored \code{b} column order.
#' @param index_map Result of \code{.two_block_index_map()}.
#' @return A list of numeric matrices, numeric vectors, 1-based integer
#'   vectors, character scalars, and the four family functions.
#' @noRd
.two_block_spec <- function(design, block1_prior, pfamily_list, family,
                            re_names, index_map) {
  if (is.null(family$family) || is.null(family$link)) {
    stop("'family' must have family and link names.", call. = FALSE)
  }
  if (identical(family$family, "gaussian")) {
    stop(
      "The indexed sweep draws the groups block with rNormalGLMGroups. Gaussian models stay on rGLMM_sweep.",
      call. = FALSE
    )
  }
  ptypes <- vapply(pfamily_list[re_names], function(pf) pf$pfamily, character(1))
  if (any(ptypes != "dNormal")) {
    stop(
      "The indexed sweep draws the population block with rNormalReg. A component that is not dNormal stays on rGLMM_sweep.",
      call. = FALSE
    )
  }

  D <- as.matrix(design$D)
  dimnames(D) <- NULL
  l1 <- ncol(D)
  l2 <- nrow(D)
  p_re <- length(re_names)
  d_col <- as.integer(index_map$d_col)
  if (!identical(sort(d_col), seq_len(l1))) {
    stop("colnames(D) must be exactly the random-effect names.", call. = FALSE)
  }

  P_re <- as.matrix(block1_prior$P)
  if (!all(dim(P_re) == c(p_re, p_re))) {
    stop("block1_prior$P must be p_re x p_re.", call. = FALSE)
  }
  P_D <- matrix(0, l1, l1)
  P_D[d_col, d_col] <- P_re

  disp1 <- block1_prior$dispersion
  if (is.null(disp1)) {
    disp1 <- 1
  }
  disp1 <- as.numeric(disp1)
  if (length(disp1) != 1L || !is.finite(disp1)) {
    stop("Groups dispersion must be a finite scalar.", call. = FALSE)
  }

  W <- vector("list", p_re)
  mu2 <- vector("list", p_re)
  P2 <- vector("list", p_re)
  disp2 <- numeric(p_re)
  for (k in seq_len(p_re)) {
    nm <- re_names[k]
    Wk <- as.matrix(design$W[[nm]])
    dimnames(Wk) <- NULL
    pl <- pfamily_list[[nm]]$prior_list
    mu <- as.numeric(pl$mu)
    if (length(mu) != ncol(Wk)) {
      stop("length(mu) for ", nm, " must equal ncol(W).", call. = FALSE)
    }
    if (!is.null(pl$P)) {
      P <- as.matrix(pl$P)
    } else {
      R <- chol(pl$Sigma)
      Pinv <- chol2inv(R)
      P <- 0.5 * (Pinv + t(Pinv))
    }
    dimnames(P) <- NULL
    if (!all(dim(P) == c(ncol(Wk), ncol(Wk)))) {
      stop("Population precision for ", nm, " has the wrong dimension.", call. = FALSE)
    }
    disp <- pl$dispersion
    if (is.null(disp) || length(disp) != 1L || !is.finite(disp) || disp <= 0) {
      stop("dNormal component ", nm, " needs a positive scalar dispersion.",
           call. = FALSE)
    }
    W[[k]] <- Wk
    mu2[[k]] <- mu
    P2[[k]] <- P
    disp2[k] <- as.numeric(disp)
  }

  fam <- .two_block_block1_glmbfamfunc(family)
  list(
    y           = as.numeric(design$y),
    D           = D,
    offset1     = numeric(l2),
    wt1         = rep(1, l2),
    obs_rows    = unname(lapply(index_map$obs_rows, as.integer)),
    d_col       = d_col,
    P_D         = P_D,
    disp1       = disp1,
    W           = W,
    w_for_group = unname(lapply(index_map$w_for_group[re_names], as.integer)),
    b_for_wrow  = unname(lapply(index_map$b_for_wrow[re_names], as.integer)),
    w_col       = unname(lapply(index_map$w_col[re_names], as.integer)),
    mu2         = mu2,
    P2          = P2,
    disp2       = disp2,
    family      = as.character(family$family),
    link        = as.character(family$link),
    gridtype    = 2L,
    n_envopt    = 1L,
    f2          = fam$f2,
    f3          = fam$f3,
    f2_gauss    = fam$f2_gauss,
    f3_gauss    = fam$f3_gauss
  )
}

#' Initial sweep state in the indexed layout
#'
#' @param start_fixef Named list of start vectors. Names of each vector must
#'   be the stored coefficient order implied by \code{index_map$w_col}.
#' @param b_start \code{J x p_re} start for every chain. Rows are aligned to
#'   \code{group_levels} and columns to \code{re_names} by name when present.
#' @return List with \code{b}, \code{fixef}, \code{iters}, \code{iters_ranef}.
#' @noRd
.two_block_init <- function(n_chains, start_fixef, b_start, spec,
                            design, re_names, group_levels) {
  n_chains <- as.integer(n_chains)
  p_re <- length(re_names)
  J <- length(group_levels)

  fixef <- vector("list", p_re)
  for (k in seq_len(p_re)) {
    nm <- re_names[k]
    cn <- colnames(as.matrix(design$W[[nm]]))
    stored <- character(length(cn))
    stored[spec$w_col[[k]]] <- cn
    g <- start_fixef[[nm]]
    if (!identical(names(g), stored)) {
      stop("names(start_fixef[[", nm, "]]) must be the stored coefficient order.",
           call. = FALSE)
    }
    fixef[[k]] <- matrix(as.numeric(g), nrow = n_chains, ncol = length(g),
                         byrow = TRUE)
  }

  b0 <- as.matrix(b_start)
  if (!is.null(rownames(b0))) {
    b0 <- b0[group_levels, , drop = FALSE]
  }
  if (!is.null(colnames(b0))) {
    b0 <- b0[, re_names, drop = FALSE]
  }
  if (!all(dim(b0) == c(J, p_re))) {
    stop("b_start must be J x p_re.", call. = FALSE)
  }
  b <- array(rep(as.numeric(b0), n_chains), dim = c(J, p_re, n_chains))

  list(
    b           = b,
    fixef       = fixef,
    iters       = matrix(0, nrow = n_chains, ncol = p_re),
    iters_ranef = numeric(n_chains)
  )
}

#' Sweep loop in R
#'
#' @param use_mode If \code{TRUE}, each block keeps its posterior mode
#'   instead of its draw. Used only to compare the R and C++ loops exactly.
#' @param observer \code{NULL}, or a list of functions
#'   \code{groups_source(m, i)}, \code{population_source(m, i, k)} (return
#'   \code{NULL} to draw, or a replacement: \code{list(coef, iters)} for
#'   groups, a coefficient vector for population),
#'   \code{groups_seen(m, i, mu_D, coef, iters)}, and
#'   \code{population_seen(m, i, k, y_k, coef)}. Same hooks as the C++
#'   \code{SweepObserver}; indices are 1-based.
#' @return State list with \code{b}, \code{fixef}, \code{iters},
#'   \code{iters_ranef}, and the per-sweep cross-chain summary of the stacked
#'   population coefficients (components in order, stored column order):
#'   \code{sweep_mean} (\code{p_tot x inner_sweeps}) and \code{sweep_cov}
#'   (\code{p_tot x p_tot x inner_sweeps}; \code{NA} when there is one chain).
#' @noRd
.two_block_sweep_r <- function(spec, state, inner_sweeps, use_mode = FALSE,
                               progbar = FALSE, stage_label = "",
                               observer = NULL) {
  b <- state$b
  fixef <- state$fixef
  iters <- state$iters
  iters_ranef <- state$iters_ranef
  J <- dim(b)[1L]
  p_re <- dim(b)[2L]
  n <- dim(b)[3L]
  inner_sweeps <- as.integer(inner_sweeps)
  progbar <- isTRUE(progbar)
  stage_label <- if (is.null(stage_label)) "" else as.character(stage_label)[1L]
  observe <- !is.null(observer)
  p_tot <- sum(vapply(fixef, ncol, integer(1)))
  sweep_mean <- matrix(0, p_tot, inner_sweeps)
  sweep_cov <- array(0, c(p_tot, p_tot, inner_sweeps))

  for (m in seq_len(inner_sweeps)) {
    prefix_groups <- .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block1")
    prefix_population <- .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block2")
    for (i in seq_len(n)) {
      mu_D <- .two_block_groups_prior_mean(spec, fixef, i)
      src <- if (observe) observer$groups_source(m, i) else NULL
      if (is.null(src)) {
        d1 <- .two_block_groups_draw(spec, mu_D)
        coef <- if (isTRUE(use_mode)) d1$coef_mode else d1$coef
        it1 <- d1$iters
      } else {
        coef <- src$coef
        it1 <- src$iters
      }
      if (observe) observer$groups_seen(m, i, mu_D, coef, it1)
      b[, , i] <- .two_block_groups_place(spec, coef)
      iters_ranef[i] <- iters_ranef[i] + it1
      if (progbar) .two_block_progress_bar(i, n, prefix = prefix_groups)
    }
    if (progbar) .two_block_progress_bar_finish(newline = FALSE)
    for (i in seq_len(n)) {
      b_i <- matrix(b[, , i], nrow = J, ncol = p_re)
      for (k in seq_len(p_re)) {
        y_k <- .two_block_population_response(spec, b_i, k)
        coef <- if (observe) observer$population_source(m, i, k) else NULL
        if (is.null(coef)) {
          d2 <- .two_block_population_draw(spec, y_k, k)
          coef <- if (isTRUE(use_mode)) d2$coef_mode else d2$coef
        }
        if (observe) observer$population_seen(m, i, k, y_k, coef)
        fixef[[k]] <- .two_block_population_write(spec, fixef[[k]], i, k, coef)
        iters[i, k] <- iters[i, k] + 1
      }
      if (progbar) .two_block_progress_bar(i, n, prefix = prefix_population)
    }
    if (progbar) .two_block_progress_bar_finish(newline = (m == inner_sweeps))
    stacked <- do.call(cbind, fixef)
    sweep_mean[, m] <- colMeans(stacked)
    sweep_cov[, , m] <- if (n > 1L) stats::cov(stacked) else NA_real_
    if (progbar && n <= 1L) {
      prefix_sweep <- if (nzchar(stage_label)) {
        sprintf("[%s] sweep %d/%d: ", stage_label, m, inner_sweeps)
      } else {
        sprintf("sweep %d/%d: ", m, inner_sweeps)
      }
      .two_block_progress_bar(m, inner_sweeps, prefix = prefix_sweep)
      .two_block_progress_bar_finish(newline = (m == inner_sweeps))
    }
  }
  list(b = b, fixef = fixef, iters = iters, iters_ranef = iters_ranef,
       sweep_mean = sweep_mean, sweep_cov = sweep_cov)
}

#' Shared argument handling for the R and C++ indexed sweeps
#' @return List with \code{spec}, \code{state}, and \code{tau2}.
#' @noRd
.two_block_prepare <- function(n_chains, start_fixef, design, block1_prior,
                               pfamily_list, family, re_names,
                               group_levels, b_start, b_mode, tau2_start,
                               index_map) {
  if (!identical(index_map$group_levels, as.character(group_levels)) ||
      !identical(index_map$re_names, as.character(re_names))) {
    stop("'index_map' does not match re_names and group_levels.", call. = FALSE)
  }
  if (is.null(b_start)) {
    if (is.null(b_mode)) {
      stop("'b_start' or 'b_mode' required for batch init.", call. = FALSE)
    }
    b_start <- b_mode
  }
  if (is.null(tau2_start)) {
    tau2_start <- .two_block_tau2_start_from_pfamily(pfamily_list, re_names)
  } else {
    if (is.null(names(tau2_start)) || !setequal(names(tau2_start), re_names)) {
      stop("'tau2_start' must be a named vector with names(re_names).",
           call. = FALSE)
    }
    tau2_start <- tau2_start[re_names]
  }
  spec <- .two_block_spec(
    design, block1_prior, pfamily_list, family, re_names, index_map
  )
  state <- .two_block_init(
    n_chains, start_fixef, b_start, spec, design, re_names, group_levels
  )
  tau2 <- matrix(as.numeric(tau2_start), nrow = n_chains,
                 ncol = length(re_names), byrow = TRUE)
  list(spec = spec, state = state, tau2 = tau2)
}

#' Sweep history from the per-sweep summary of an indexed sweep
#'
#' Same object \code{rGLMM_sweep} stores: per-sweep chain mean and SD per
#' coefficient, and the stacked covariance with its \code{coef_index}.
#' \code{dNormal} components keep a fixed dispersion, so there are no
#' dispersion statistics.
#' @param state Sweep state with \code{sweep_mean} and \code{sweep_cov}.
#' @return Object of class \code{"two_block_sweep_history"}.
#' @noRd
.two_block_sweep_history_indexed <- function(state, index_map, design,
                                             re_names, stage_label,
                                             fixef_mode) {
  stored <- lapply(re_names, function(k) {
    cn <- colnames(as.matrix(design$W[[k]]))
    out <- character(length(cn))
    out[index_map$w_col[[k]]] <- cn
    out
  })
  names(stored) <- re_names
  coef_index <- data.frame(
    re_component = rep(re_names, lengths(stored)),
    covariate    = unlist(stored, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  lbl <- paste(coef_index$re_component, coef_index$covariate, sep = " | ")
  idx <- split(seq_len(nrow(coef_index)),
               factor(coef_index$re_component, levels = re_names))

  M <- ncol(state$sweep_mean)
  sweep_stats <- vector("list", M)
  sweep_cov <- vector("list", M)
  for (m in seq_len(M)) {
    mu <- state$sweep_mean[, m]
    V <- matrix(state$sweep_cov[, , m], nrow(coef_index), nrow(coef_index))
    sdv <- sqrt(diag(V))
    sweep_stats[[m]] <- lapply(re_names, function(k) {
      list(mean = stats::setNames(mu[idx[[k]]], stored[[k]]),
           sd   = stats::setNames(sdv[idx[[k]]], stored[[k]]))
    })
    names(sweep_stats[[m]]) <- re_names
    dimnames(V) <- list(lbl, lbl)
    sweep_cov[[m]] <- list(cov = V, coef_index = coef_index)
  }
  .two_block_build_sweep_history(
    stage_label = stage_label,
    sweep_stats = sweep_stats,
    fixef_mode  = fixef_mode,
    re_names    = re_names,
    sweep_cov   = sweep_cov,
    disp_stats  = NULL
  )
}

#' Restore names on an indexed sweep state and pack the output
#' @param stage_label,fixef_mode Stored on \code{sweep_history}.
#' @param diag_sweeps If \code{TRUE}, print the sweep history, as
#'   \code{rGLMM_sweep} does.
#' @return The list from \code{.rGLMM_sweep_save_indexed()} plus
#'   \code{sweep_history}.
#' @noRd
.two_block_finish <- function(state, tau2, index_map, design, re_names,
                              collect_block1, stage_label = "",
                              fixef_mode = NULL, diag_sweeps = FALSE) {
  fixef <- state$fixef
  names(fixef) <- re_names
  iters <- state$iters
  colnames(iters) <- re_names
  out <- .rGLMM_sweep_save_indexed(
    n              = dim(state$b)[3L],
    fixef          = fixef,
    tau2           = tau2,
    b              = state$b,
    iters          = iters,
    iters_ranef    = state$iters_ranef,
    map            = index_map,
    design         = design,
    collect_block1 = collect_block1
  )
  out$sweep_history <- .two_block_sweep_history_indexed(
    state, index_map, design, re_names, stage_label, fixef_mode
  )
  if (isTRUE(diag_sweeps)) {
    print(out$sweep_history)
  }
  out
}

#' Pack indexed-sweep output and restore the group factor
#' @noRd
.rGLMM_sweep_save_indexed <- function(
    n,
    fixef,
    tau2,
    b,
    iters,
    iters_ranef,
    map,
    design,
    collect_block1 = TRUE
) {
  re_names <- map$re_names
  group_levels <- map$group_levels
  J <- length(group_levels)
  p_re <- length(re_names)

  dimnames(b) <- list(group_levels, re_names, NULL)
  colnames(tau2) <- re_names
  for (k in re_names) {
    cn <- colnames(as.matrix(design$W[[k]]))
    stored <- character(length(cn))
    stored[map$w_col[[k]]] <- cn
    colnames(fixef[[k]]) <- stored
  }

  fixef_draws <- lapply(fixef, function(mat) {
    out <- mat
    rownames(out) <- NULL
    out
  })
  names(fixef_draws) <- re_names

  if (isTRUE(collect_block1)) {
    grp_col <- design$group_name
    if (is.null(grp_col) || !nzchar(grp_col)) grp_col <- "group"
    lev_fac <- factor(group_levels, levels = group_levels)
    coef_rows <- vector("list", n)
    for (i in seq_len(n)) {
      draw_df <- data.frame(
        draw = rep(i, J),
        stringsAsFactors = FALSE
      )
      draw_df[[grp_col]] <- lev_fac
      for (k in re_names) {
        draw_df[[k]] <- b[, match(k, re_names), i]
      }
      coef_rows[[i]] <- draw_df
    }
    coefficients <- do.call(rbind, coef_rows)
    rownames(coefficients) <- NULL
  } else {
    coefficients <- NULL
  }

  fixef_mean <- lapply(fixef, function(mat) {
    v <- colMeans(mat)
    names(v) <- colnames(mat)
    v
  })
  mu_all_last <- .two_block_index_mu_all(map, fixef_mean, design)

  list(
    fixef_draws            = fixef_draws,
    dispersion_fixef_draws = tau2,
    iters_fixef_draws      = iters,
    iters_ranef_draws      = iters_ranef,
    coefficients           = coefficients,
    mu_all_last            = mu_all_last
  )
}

#' Indexed two-block sweep, R loop
#'
#' Same arguments as \code{rGLMM_sweep}, plus \code{index_map}. Non-Gaussian
#' families with \code{dNormal} components only. \code{rwt} is accepted and
#' is not used inside the draw.
#'
#' @param index_map Result of \code{.two_block_index_map()}.
#' @noRd
.rGLMM_sweep_indexed <- function(
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
  force(rwt)
  force(diag_sweeps)
  force(fixef_mode)
  force(ptypes)
  force(use_cpp_tau2_row)
  force(use_cpp_b_slice)
  force(use_cpp_iters_ranef_add)
  force(use_cpp_block2)

  prep <- .two_block_prepare(
    n_chains, start_fixef, design, block1_prior, pfamily_list, family,
    re_names, group_levels, b_start, b_mode, tau2_start, index_map
  )
  state <- .two_block_sweep_r(
    spec         = prep$spec,
    state        = prep$state,
    inner_sweeps = inner_sweeps,
    use_mode     = FALSE,
    progbar      = progbar,
    stage_label  = stage_label
  )
  .two_block_finish(
    state, prep$tau2, index_map, design, re_names, collect_block1,
    stage_label = stage_label, fixef_mode = fixef_mode,
    diag_sweeps = diag_sweeps
  )
}

#' Indexed two-block sweep, C++ loop
#'
#' Same arguments and output as \code{.rGLMM_sweep_indexed}. Setup and the
#' final naming are the shared R steps; the sweep loop is
#' \code{.two_block_sweep_cpp}, the C++ twin of \code{.two_block_sweep_r}.
#' \code{rwt} is accepted and is not used inside the draw.
#' @param index_map Result of \code{.two_block_index_map()}.
#' @noRd
.rGLMM_sweep_cpp <- function(
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
  force(rwt)
  force(diag_sweeps)
  force(fixef_mode)
  force(ptypes)
  force(use_cpp_tau2_row)
  force(use_cpp_b_slice)
  force(use_cpp_iters_ranef_add)
  force(use_cpp_block2)

  prep <- .two_block_prepare(
    n_chains, start_fixef, design, block1_prior, pfamily_list, family,
    re_names, group_levels, b_start, b_mode, tau2_start, index_map
  )
  state <- .two_block_sweep_cpp(
    prep$spec,
    prep$state$b,
    prep$state$fixef,
    prep$state$iters,
    prep$state$iters_ranef,
    as.integer(inner_sweeps),
    FALSE,
    isTRUE(progbar),
    if (is.null(stage_label)) "" else as.character(stage_label[1L])
  )
  .two_block_finish(
    state, prep$tau2, index_map, design, re_names, collect_block1,
    stage_label = stage_label, fixef_mode = fixef_mode,
    diag_sweeps = diag_sweeps
  )
}
