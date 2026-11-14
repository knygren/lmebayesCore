## Portable two-block sweep: the R reference for src/two_block_sweep_cpp.cpp.
##
## Names are resolved once, in .two_block_port_spec() and
## .two_block_port_init(). After that every step works on plain numeric
## matrices, numeric vectors, and 1-based integer index vectors, and each
## step has a C++ twin of the same name stem with the same inputs and the
## same output shape:
##
##   R                                    C++
##   .two_block_port_block1_prior_mean    port_block1_prior_mean
##   .two_block_port_block1_draw          port_block1_draw
##   .two_block_port_block1_place         port_block1_place
##   .two_block_port_block2_response      port_block2_response
##   .two_block_port_block2_draw          port_block2_draw
##   .two_block_port_block2_write         port_block2_write
##   .two_block_port_sweep_r              port_sweep
##
## The draws in both languages come from the same compiled samplers:
## rNormalGLMBlocks (Block 1) and rNormalReg (Block 2). Names are restored
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
#' @param block1_prior Block~1 prior with \code{P} (\code{p_re x p_re},
#'   \code{re_names} order) and optional \code{dispersion}.
#' @param pfamily_list Named list of \code{dNormal} pfamily objects.
#' @param family GLM family. Gaussian is not supported.
#' @param re_names Random-effect names, in stored \code{b} column order.
#' @param index_map Result of \code{.two_block_index_map()}.
#' @return A list of numeric matrices, numeric vectors, 1-based integer
#'   vectors, character scalars, and the four family functions.
#' @noRd
.two_block_port_spec <- function(design, block1_prior, pfamily_list, family,
                                 re_names, index_map) {
  if (is.null(family$family) || is.null(family$link)) {
    stop("'family' must have family and link names.", call. = FALSE)
  }
  if (identical(family$family, "gaussian")) {
    stop(
      "The portable sweep draws Block 1 with rNormalGLMBlocks. Gaussian models stay on rGLMM_sweep.",
      call. = FALSE
    )
  }
  ptypes <- vapply(pfamily_list[re_names], function(pf) pf$pfamily, character(1))
  if (any(ptypes != "dNormal")) {
    stop(
      "The portable sweep draws Block 2 with rNormalReg. A component that is not dNormal stays on rGLMM_sweep.",
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
    stop("Block 1 dispersion must be a finite scalar.", call. = FALSE)
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
      stop("Block 2 precision for ", nm, " has the wrong dimension.", call. = FALSE)
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

#' Initial sweep state in the portable layout
#'
#' @param start_fixef Named list of start vectors. Names of each vector must
#'   be the stored coefficient order implied by \code{index_map$w_col}.
#' @param b_start \code{J x p_re} start for every chain. Rows are aligned to
#'   \code{group_levels} and columns to \code{re_names} by name when present.
#' @return List with \code{b}, \code{fixef}, \code{iters}, \code{iters_ranef}.
#' @noRd
.two_block_port_init <- function(n_chains, start_fixef, b_start, spec,
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

#' Block~1 prior mean for chain \code{i}
#'
#' Row \code{d_col[k]}, column \code{j} is
#' \code{sum_c W_k[w_for_group_k[j], c] * fixef_k[i, w_col_k[c]]}, summed over
#' \code{c} in increasing order.
#' @return Numeric matrix \code{ncol(D) x J}.
#' @noRd
.two_block_port_block1_prior_mean <- function(spec, fixef, i) {
  l1 <- ncol(spec$D)
  J <- length(spec$obs_rows)
  mu_D <- matrix(0, l1, J)
  for (k in seq_along(spec$d_col)) {
    Wk <- spec$W[[k]]
    wg <- spec$w_for_group[[k]]
    wc <- spec$w_col[[k]]
    s <- numeric(J)
    for (c in seq_len(ncol(Wk))) {
      s <- s + Wk[wg, c] * fixef[[k]][i, wc[c]]
    }
    mu_D[spec$d_col[k], ] <- s
  }
  mu_D
}

#' Block~1 draw for one chain
#' @return List with \code{coef} and \code{coef_mode} (both \code{J x ncol(D)})
#'   and \code{iters} (mean envelope candidates per group).
#' @noRd
.two_block_port_block1_draw <- function(spec, mu_D) {
  out <- rNormalGLMBlocks_cpp_export(
    n              = 1L,
    y              = spec$y,
    x              = spec$D,
    offset         = spec$offset1,
    wt             = spec$wt1,
    dispersion     = spec$disp1,
    mu             = mu_D,
    P_blocks       = list(spec$P_D),
    prior_by_block = TRUE,
    row_blocks     = spec$obs_rows,
    f2             = spec$f2,
    f3             = spec$f3,
    family         = spec$family,
    link           = spec$link,
    Gridtype       = spec$gridtype,
    n_envopt       = spec$n_envopt,
    use_parallel   = FALSE,
    use_opencl     = FALSE,
    verbose        = FALSE
  )
  it <- vapply(out$group_results, function(g) as.numeric(g$iters)[1L], numeric(1))
  list(
    coef      = out$coefficients,
    coef_mode = out$coef.mode,
    iters     = sum(it) / length(it)
  )
}

#' Columns of a Block~1 draw in stored \code{b} order
#' @return Numeric matrix \code{J x p_re}.
#' @noRd
.two_block_port_block1_place <- function(spec, coef) {
  out <- coef[, spec$d_col, drop = FALSE]
  dimnames(out) <- NULL
  out
}

#' Block~2 response for component \code{k}, in the row order of \code{W_k}
#' @return Numeric vector of length \code{nrow(W_k)}.
#' @noRd
.two_block_port_block2_response <- function(spec, b_i, k) {
  b_i[spec$b_for_wrow[[k]], k]
}

#' Block~2 draw for component \code{k}
#' @return List with \code{coef} and \code{coef_mode}, numeric vectors in
#'   \code{W_k} column order.
#' @noRd
.two_block_port_block2_draw <- function(spec, y_k, k) {
  n_w <- length(y_k)
  out <- rNormalReg_cpp_export(
    n          = 1L,
    y          = y_k,
    x          = spec$W[[k]],
    mu         = spec$mu2[[k]],
    P          = spec$P2[[k]],
    offset     = numeric(n_w),
    wt         = rep(1, n_w),
    dispersion = spec$disp2[k],
    f2         = spec$f2_gauss,
    f3         = spec$f3_gauss,
    start      = spec$mu2[[k]],
    family     = "gaussian",
    link       = "identity",
    Gridtype   = spec$gridtype
  )
  list(
    coef      = as.numeric(out$coefficients[1L, ]),
    coef_mode = as.numeric(out$coef.mode)
  )
}

#' Write a Block~2 draw into row \code{i} of the stored fixef matrix
#' @return \code{fixef_k} with \code{fixef_k[i, w_col_k] <- coef}.
#' @noRd
.two_block_port_block2_write <- function(spec, fixef_k, i, k, coef) {
  fixef_k[i, spec$w_col[[k]]] <- coef
  fixef_k
}

#' Empty input and draw log for a portable sweep
#'
#' A filled log is also a valid \code{tape}: replaying it reproduces the run.
#' @return List with \code{mu_D} (\code{ncol(D) x J x n x M}), \code{block1}
#'   (\code{J x ncol(D) x n x M}), \code{iters1} (\code{n x M}), and lists
#'   \code{y} (\code{nrow(W_k) x n x M}) and \code{block2}
#'   (\code{ncol(W_k) x n x M}) over components.
#' @noRd
.two_block_port_log_new <- function(spec, n, M) {
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

#' Sweep loop in R
#'
#' @param use_mode If \code{TRUE}, each block keeps its posterior mode
#'   instead of its draw. Used only to compare the R and C++ loops exactly.
#' @param tape \code{NULL}, or a log from \code{.two_block_port_log_new()}
#'   layout; its \code{block1}, \code{iters1}, and \code{block2} entries
#'   replace the samplers.
#' @param record If \code{TRUE}, return the sampler inputs and draws as
#'   \code{log}.
#' @return State list with \code{b}, \code{fixef}, \code{iters},
#'   \code{iters_ranef}, \code{log} (\code{NULL} unless \code{record}), and
#'   \code{prior_sum}, \code{prior_sumsq}: \code{ncol(D) x J x inner_sweeps}
#'   sums over chains of the Block~1 prior mean and of its square.
#' @noRd
.two_block_port_sweep_r <- function(spec, state, inner_sweeps, use_mode = FALSE,
                                    progbar = FALSE, stage_label = "",
                                    tape = NULL, record = FALSE) {
  b <- state$b
  fixef <- state$fixef
  iters <- state$iters
  iters_ranef <- state$iters_ranef
  J <- dim(b)[1L]
  p_re <- dim(b)[2L]
  n <- dim(b)[3L]
  l1 <- ncol(spec$D)
  inner_sweeps <- as.integer(inner_sweeps)
  progbar <- isTRUE(progbar)
  stage_label <- if (is.null(stage_label)) "" else as.character(stage_label)[1L]
  log <- if (isTRUE(record)) .two_block_port_log_new(spec, n, inner_sweeps) else NULL
  prior_sum <- array(0, c(l1, J, inner_sweeps))
  prior_sumsq <- array(0, c(l1, J, inner_sweeps))

  for (m in seq_len(inner_sweeps)) {
    prefix_b1 <- .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block1")
    prefix_b2 <- .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block2")
    for (i in seq_len(n)) {
      mu_D <- .two_block_port_block1_prior_mean(spec, fixef, i)
      prior_sum[, , m] <- prior_sum[, , m] + mu_D
      prior_sumsq[, , m] <- prior_sumsq[, , m] + mu_D * mu_D
      if (is.null(tape)) {
        d1 <- .two_block_port_block1_draw(spec, mu_D)
        coef <- if (isTRUE(use_mode)) d1$coef_mode else d1$coef
        it1 <- d1$iters
      } else {
        coef <- matrix(tape$block1[, , i, m], nrow = J, ncol = l1)
        it1 <- tape$iters1[i, m]
      }
      if (!is.null(log)) {
        log$mu_D[, , i, m] <- mu_D
        log$block1[, , i, m] <- coef
        log$iters1[i, m] <- it1
      }
      b[, , i] <- .two_block_port_block1_place(spec, coef)
      iters_ranef[i] <- iters_ranef[i] + it1
      if (progbar) .two_block_progress_bar(i, n, prefix = prefix_b1)
    }
    if (progbar) .two_block_progress_bar_finish(newline = FALSE)
    for (i in seq_len(n)) {
      b_i <- matrix(b[, , i], nrow = J, ncol = p_re)
      for (k in seq_len(p_re)) {
        y_k <- .two_block_port_block2_response(spec, b_i, k)
        if (is.null(tape)) {
          d2 <- .two_block_port_block2_draw(spec, y_k, k)
          coef <- if (isTRUE(use_mode)) d2$coef_mode else d2$coef
        } else {
          coef <- tape$block2[[k]][, i, m]
        }
        if (!is.null(log)) {
          log$y[[k]][, i, m] <- y_k
          log$block2[[k]][, i, m] <- coef
        }
        fixef[[k]] <- .two_block_port_block2_write(spec, fixef[[k]], i, k, coef)
        iters[i, k] <- iters[i, k] + 1
      }
      if (progbar) .two_block_progress_bar(i, n, prefix = prefix_b2)
    }
    if (progbar) .two_block_progress_bar_finish(newline = (m == inner_sweeps))
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
       log = log, prior_sum = prior_sum, prior_sumsq = prior_sumsq)
}

#' Block~1 conditional prior fed to each group, summarized over chains
#'
#' @return List with \code{prior_mean} and \code{prior_mean_var}
#'   (\code{p_re x J x inner_sweeps}: mean and variance over chains of the
#'   prior mean fed to each group in each sweep), \code{prior_precision}
#'   (\code{p_re x p_re}, the same for every group, chain, and sweep), and
#'   \code{prior_var} (named vector \code{diag(solve(prior_precision))}).
#' @noRd
.two_block_port_prior_track <- function(state, spec, re_names, group_levels) {
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

#' Shared argument handling for the R and C++ portable sweeps
#' @return List with \code{spec}, \code{state}, and \code{tau2}.
#' @noRd
.two_block_port_prepare <- function(n_chains, start_fixef, design, block1_prior,
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
  spec <- .two_block_port_spec(
    design, block1_prior, pfamily_list, family, re_names, index_map
  )
  state <- .two_block_port_init(
    n_chains, start_fixef, b_start, spec, design, re_names, group_levels
  )
  tau2 <- matrix(as.numeric(tau2_start), nrow = n_chains,
                 ncol = length(re_names), byrow = TRUE)
  list(spec = spec, state = state, tau2 = tau2)
}

#' Restore names on a portable sweep state and pack the output
#' @noRd
.two_block_port_finish <- function(state, tau2, index_map, design, re_names,
                                   collect_block1, spec) {
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
  out$prior_track <- .two_block_port_prior_track(
    state, spec, re_names, index_map$group_levels
  )
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

#' Portable two-block sweep, R loop
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

  prep <- .two_block_port_prepare(
    n_chains, start_fixef, design, block1_prior, pfamily_list, family,
    re_names, group_levels, b_start, b_mode, tau2_start, index_map
  )
  state <- .two_block_port_sweep_r(
    spec         = prep$spec,
    state        = prep$state,
    inner_sweeps = inner_sweeps,
    use_mode     = FALSE,
    progbar      = progbar,
    stage_label  = stage_label
  )
  .two_block_port_finish(
    state, prep$tau2, index_map, design, re_names, collect_block1, prep$spec
  )
}
