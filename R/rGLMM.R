#' Replicate-chain Gibbs sampling for Bayesian GLMMs
#'
#' Runs optional pilot replicate chains, an optional post-pilot eigenvalue
#' upper-bound calibration, then main replicate chains for two-block mixed
#' models. Block~2 priors follow the \code{\link{two_block_rNormal_reg_v2}}
#' contract (\code{pfamily_list}).
#'
#' When \code{n_pilot = 0}, this reduces to a single main-stage v2 draw (via
#' the staged C++ driver). When \code{n_pilot > 0} and \code{tv_tol} is
#' \code{NULL}, pilot sampling, the Hotelling chi-squared test vs \code{start},
#' and main sampling run in one C++ call. When \code{n_pilot > 0} and
#' \code{tv_tol} is set, the pilot runs in C++, eigenvalue upper bounds are
#' computed in R (\code{\link{two_block_mode_weights}},
#' \code{\link{two_block_rate_v2}}, \code{\link{two_block_l_for_tv}}), and the
#' main stage is rerun with the updated \code{m_convergence}.
#'
#' @param n Number of stored main-stage draws. If \code{length(n) > 1}, the
#'   length is taken to be the number required.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Level-1 design matrix \code{Z} (\code{l2 x p_re}).
#' @param block Grouping factor or block partition of length \code{l2}.
#' @param x_hyper Named list of group-level design matrices \code{X_k}
#'   (\code{J x q_k}), one per column of \code{x}.
#' @param prior_list Prior for Block~1: \code{P} or \code{Sigma},
#'   \code{dispersion} (required for \code{gaussian()}), optional \code{ddef}.
#' @param pfamily_list Named list of \code{pfamily} objects for Block~2.
#' @param start Named list of Block~2 hyper-parameter vectors at which each
#'   inner chain is initialised (pilot stage uses \code{start}; main stage uses
#'   the pilot mean when \code{n_pilot > 0}).
#' @param offset,weights,family Passed to Block~1 (length \code{l2} or recycled).
#' @param m_convergence Inner Gibbs steps per main-stage stored draw.
#' @param re_coef_names Character vector naming columns of \code{x}.
#' @param group_levels Character vector defining row order of Block~1 draws.
#' @param group_name Name for the grouping column in \code{coefficients}.
#' @param n_pilot Number of pilot replicate chains (\code{0L} skips the pilot).
#' @param m_convergence_pilot Inner Gibbs steps per pilot stored draw. Defaults
#'   to \code{m_convergence} when \code{n_pilot > 0}.
#' @param tv_tol When non-\code{NULL} and \code{n_pilot > 0}, run the
#'   post-pilot eigenvalue upper-bound calibration before the main stage.
#' @param mode_gap_max Reserved for pilot sweep calibration when
#'   \code{m_convergence_pilot} is \code{NULL} (not yet implemented).
#' @param Gridtype,n_envopt,use_parallel,use_opencl,verbose,progbar Passed to
#'   Block~1 when \code{family} is not Gaussian. \code{n_envopt} defaults to
#'   \code{n}. \code{progbar} controls progress bars for both pilot and main
#'   stages when a pilot runs.
#' @return Object of class \code{c("rGLMM", "two_block_rNormal_reg_v2",
#'   "two_block_rNormal_reg")}.  Main-stage Block~2 fields use the
#'   \code{fixef.*} namespace (\code{fixef}, \code{fixef.last}, \code{fixef.mu},
#'   \code{fixef.dispersion}, \code{fixef.iters}, \code{fixef.mode},
#'   \code{fixef.init}); Block~1 draws remain in \code{coefficients} and
#'   \code{coef.last}; layout metadata uses \code{coef.names}.  Also includes
#'   \code{m_convergence} (main stage), \code{n_pilot},
#'   \code{m_convergence_pilot} when a pilot runs, and, when applicable,
#'   \code{pilot}, \code{pilot_chisq}, \code{pilot_ub}.
#' @family simfuncs
#' @seealso \code{\link{two_block_rNormal_reg_v2}}
#' @export
rGLMM <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list,
    pfamily_list,
    start,
    offset              = NULL,
    weights             = 1,
    family              = gaussian(),
    m_convergence       = 10L,
    re_coef_names       = colnames(x),
    group_levels        = levels(block),
    group_name          = NULL,
    n_pilot             = 0L,
    m_convergence_pilot = NULL,
    tv_tol              = NULL,
    mode_gap_max        = 1.0,
    Gridtype            = 2,
    n_envopt            = NULL,
    use_parallel        = TRUE,
    use_opencl          = FALSE,
    verbose             = FALSE,
    progbar             = FALSE) {

  cl <- match.call()

  family <- .two_block_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")

  if (length(n) > 1L) {
    n <- length(n)
  }
  n <- as.integer(n)
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  n_pilot <- as.integer(n_pilot[1L])
  if (n_pilot < 0L) {
    stop("'n_pilot' must be non-negative.", call. = FALSE)
  }

  m_convergence <- as.integer(m_convergence[1L])
  if (m_convergence < 1L) {
    stop("'m_convergence' must be at least 1.", call. = FALSE)
  }

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be a single positive finite number.",
           call. = FALSE)
    }
  }

  run_pilot <- n_pilot > 0L
  run_ub <- run_pilot && !is.null(tv_tol)

  m_convergence_pilot <- if (run_pilot) {
    if (is.null(m_convergence_pilot)) {
      m_convergence
    } else {
      as.integer(m_convergence_pilot[1L])
    }
  } else {
    NA_integer_
  }
  if (run_pilot && m_convergence_pilot < 1L) {
    stop("'m_convergence_pilot' must be at least 1 when n_pilot > 0.",
         call. = FALSE)
  }

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  if (is.null(re_coef_names) || length(re_coef_names) != ncol(x)) {
    re_coef_names <- if (ncol(x) >= 1L) {
      cn <- colnames(x)
      if (is.null(cn) || length(cn) != ncol(x)) {
        paste0("RE", seq_len(ncol(x)))
      } else {
        cn
      }
    } else {
      stop("'x' must have at least one column.", call. = FALSE)
    }
  }
  colnames(x) <- re_coef_names
  re_names <- re_coef_names

  group_levels <- as.character(group_levels)
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    group_name <- tryCatch(
      deparse(substitute(block))[1L],
      error = function(e) "group"
    )
    if (!nzchar(group_name)) group_name <- "group"
  }

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop(
      "length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
      call. = FALSE
    )
  }
  if (!setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, re_names,
    J = length(group_levels)
  )

  if (!is.list(start) || is.null(names(start))) {
    stop("'start' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(start), re_names)) {
    stop("names(start) must match re_coef_names.", call. = FALSE)
  }
  start <- start[re_names]
  fixef_mode <- start

  block1_prior_meta <- .two_block_validate_block1_prior(
    prior_list,
    family = family
  )

  offset2 <- offset
  wt <- weights
  if (is.null(offset2)) {
    offset2 <- rep(0, l2)
  } else {
    offset2 <- as.numeric(offset2)
    if (length(offset2) == 1L) offset2 <- rep(offset2, l2)
    if (length(offset2) != l2) {
      stop("length(offset) must be 1 or length(y).", call. = FALSE)
    }
  }
  if (length(wt) == 1L) wt <- rep(wt, l2)
  if (length(wt) != l2) {
    stop("length(weights) must be 1 or length(y).", call. = FALSE)
  }

  famfunc_block1 <- glmbfamfunc(if (is_gaussian) gaussian() else family)
  famfunc_gauss <- glmbfamfunc(gaussian())
  n_envopt_use <- if (is.null(n_envopt)) n else as.integer(n_envopt)
  x_hyper_mats <- lapply(x_hyper, as.matrix)

  J <- length(group_levels)
  p_re <- length(re_names)
  progbar_flag <- isTRUE(progbar)

  cpp_common <- list(
    y                 = y,
    x                 = x,
    block             = block,
    x_hyper           = x_hyper_mats,
    prior_list_block1 = prior_list,
    dispersion_block1 = block1_prior_meta$dispersion,
    ddef_block1       = block1_prior_meta$ddef,
    pfamily_list      = pfamily_list,
    group_levels      = group_levels,
    family            = family$family,
    link              = family$link,
    f2                = famfunc_block1$f2,
    f3                = famfunc_block1$f3,
    f2_gauss          = famfunc_gauss$f2,
    f3_gauss          = famfunc_gauss$f3,
    offset            = offset2,
    wt                = wt,
    Gridtype          = as.integer(Gridtype),
    n_envopt          = n_envopt_use,
    use_parallel      = use_parallel,
    use_opencl        = use_opencl,
    verbose           = verbose
  )

  pilot_res     <- NULL
  pilot_chisq   <- NULL
  pilot_ub      <- NULL
  fixef_init <- fixef_mode
  m_convergence_used <- m_convergence

  if (!run_pilot) {
    cpp_out <- do.call(
      .two_block_rNormal_reg_staged_cpp,
      c(
        list(
          n_main = n,
          m_convergence_main = m_convergence,
          n_pilot = 0L,
          m_convergence_pilot = 1L,
          fixef_start = fixef_mode,
          progbar_main = progbar_flag,
          progbar_pilot = FALSE
        ),
        cpp_common
      )
    )
    main_res <- .two_block_format_v2_cpp_out(
      cpp_out         = cpp_out,
      n               = n,
      re_names        = re_names,
      fixef_start     = fixef_init,
      group_levels    = group_levels,
      group_name      = group_name,
      pfamily_list    = pfamily_list,
      family          = family,
      m_convergence   = m_convergence_used,
      sampling        = "replicate",
      cl              = cl
    )
  } else if (!run_ub) {
    cpp_out <- do.call(
      .two_block_rNormal_reg_staged_cpp,
      c(
        list(
          n_main = n,
          m_convergence_main = m_convergence,
          n_pilot = n_pilot,
          m_convergence_pilot = m_convergence_pilot,
          fixef_start = fixef_mode,
          progbar_main = progbar_flag,
          progbar_pilot = progbar_flag
        ),
        cpp_common
      )
    )
    fixef_init <- cpp_out$fixef_main_start
    m_convergence_used <- cpp_out$m_convergence_used
    pilot_chisq <- cpp_out$pilot_chisq
    pilot_res <- .two_block_format_v2_cpp_out(
      cpp_out         = cpp_out$pilot,
      n               = n_pilot,
      re_names        = re_names,
      fixef_start     = fixef_mode,
      group_levels    = group_levels,
      group_name      = group_name,
      pfamily_list    = pfamily_list,
      family          = family,
      m_convergence   = m_convergence_pilot,
      sampling        = "replicate",
      cl              = cl
    )
    main_res <- .two_block_format_v2_cpp_out(
      cpp_out         = cpp_out,
      n               = n,
      re_names        = re_names,
      fixef_start     = fixef_init,
      group_levels    = group_levels,
      group_name      = group_name,
      pfamily_list    = pfamily_list,
      family          = family,
      m_convergence   = m_convergence_used,
      sampling        = "replicate",
      cl              = cl
    )
  } else {
    pilot_cpp <- do.call(
      .two_block_rNormal_reg_v2_cpp,
      c(
        list(
          n = n_pilot,
          m_convergence = m_convergence_pilot,
          fixef_start = fixef_mode,
          progbar = progbar_flag
        ),
        cpp_common
      )
    )
    pilot_res <- .two_block_format_v2_cpp_out(
      cpp_out         = pilot_cpp,
      n               = n_pilot,
      re_names        = re_names,
      fixef_start     = fixef_mode,
      group_levels    = group_levels,
      group_name      = group_name,
      pfamily_list    = pfamily_list,
      family          = family,
      m_convergence   = m_convergence_pilot,
      sampling        = "replicate",
      cl              = cl
    )
    pilot_chisq <- .two_block_pilot_chisq_test(
      fixef_draws = pilot_res$fixef_draws,
      re_names    = re_names,
      fixef_mode  = fixef_mode,
      n_pilot     = n_pilot
    )
    fixef_init <- .two_block_fixef_colmeans(
      pilot_res$fixef_draws, re_names, fixef_mode
    )
    pilot_ub <- .two_block_pilot_eigenvalue_ub(
      b_draws           = pilot_cpp$b_draws,
      n_pilot           = n_pilot,
      J                 = J,
      p_re              = p_re,
      re_names          = re_names,
      group_levels      = group_levels,
      x                 = x,
      block             = block,
      x_hyper           = x_hyper_mats,
      prior_list        = prior_list,
      pfamily_list      = pfamily_list,
      family            = family,
      dispersion_block1 = block1_prior_meta$dispersion,
      tv_tol            = tv_tol
    )
    if (pilot_ub$m_min_upper > m_convergence_used) {
      m_convergence_used <- pilot_ub$m_min_upper
    }
    main_cpp <- do.call(
      .two_block_rNormal_reg_v2_cpp,
      c(
        list(
          n = n,
          m_convergence = m_convergence_used,
          fixef_start = fixef_init,
          progbar = progbar_flag
        ),
        cpp_common
      )
    )
    main_res <- .two_block_format_v2_cpp_out(
      cpp_out         = main_cpp,
      n               = n,
      re_names        = re_names,
      fixef_start     = fixef_init,
      group_levels    = group_levels,
      group_name      = group_name,
      pfamily_list    = pfamily_list,
      family          = family,
      m_convergence   = m_convergence_used,
      sampling        = "replicate",
      cl              = cl
    )
  }

  if (run_pilot && !is.null(pilot_res)) {
    pilot_res <- .two_block_as_staged_names(
      pilot_res,
      fixef_mode = fixef_mode,
      fixef_init = fixef_mode
    )
  }
  main_res <- .two_block_as_staged_names(
    main_res,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
  main_res$n_pilot <- n_pilot
  main_res$m_convergence <- m_convergence_used
  main_res$m_convergence_pilot <- if (run_pilot) m_convergence_pilot else NULL
  if (run_pilot) {
    main_res$pilot <- pilot_res
    main_res$pilot_chisq <- pilot_chisq
  }
  if (run_ub) {
    main_res$pilot_ub <- pilot_ub
    main_res$tv_tol <- tv_tol
  }
  class(main_res) <- c(
    "rGLMM",
    "two_block_rNormal_reg_v2",
    "two_block_rNormal_reg"
  )
  main_res
}

#' Rename v2 result fields to the rGLMM \code{fixef.*} namespace
#' @noRd
.two_block_as_staged_names <- function(x, fixef_mode, fixef_init) {
  renames <- c(
    fixef_draws            = "fixef",
    fixef_last             = "fixef.last",
    b_last                 = "coef.last",
    mu_all_last            = "fixef.mu",
    dispersion_fixef_draws = "fixef.dispersion",
    iters_fixef_draws      = "fixef.iters",
    re_coef_names          = "coef.names"
  )
  for (old_nm in names(renames)) {
    if (!is.null(x[[old_nm]])) {
      x[[renames[[old_nm]]]] <- x[[old_nm]]
      x[[old_nm]] <- NULL
    }
  }
  x$fixef_start <- NULL
  x$sampling <- NULL
  x$m_convergence <- NULL
  x$fixef.mode <- fixef_mode
  x$fixef.init <- fixef_init
  x
}

#' Hotelling chi-squared test: pilot fixef mean vs start
#' @noRd
.two_block_pilot_chisq_test <- function(fixef_draws, re_names, fixef_mode,
                                          n_pilot) {
  X_pilot <- do.call(cbind, lapply(re_names, function(k) fixef_draws[[k]]))
  pnames <- unlist(lapply(re_names, function(k) {
    paste0(k, "::", colnames(fixef_draws[[k]]))
  }))
  colnames(X_pilot) <- pnames
  mu_pilot <- colMeans(X_pilot)
  mode_vec <- unlist(lapply(re_names, function(k) fixef_mode[[k]]))
  names(mode_vec) <- pnames
  d_pm <- mu_pilot - mode_vec
  S_pilot <- stats::cov(X_pilot)
  p_dim2 <- ncol(X_pilot)
  S_inv <- tryCatch(
    solve(S_pilot),
    error = function(e) {
      ridge <- 1e-8 * mean(diag(S_pilot))
      solve(S_pilot + diag(ridge, p_dim2))
    }
  )
  Q_pm <- as.numeric(n_pilot * t(d_pm) %*% S_inv %*% d_pm)
  p_pm <- stats::pchisq(Q_pm, df = p_dim2, lower.tail = FALSE)
  list(
    Q       = Q_pm,
    df      = p_dim2,
    p_value = p_pm,
    n_pilot = n_pilot
  )
}

#' Column means of fixef draws with names copied from fixef mode
#' @noRd
.two_block_fixef_colmeans <- function(fixef_draws, re_names, fixef_mode) {
  out <- lapply(re_names, function(k) {
    v <- colMeans(fixef_draws[[k]])
    names(v) <- names(fixef_mode[[k]])
    v
  })
  stats::setNames(out, re_names)
}

#' Post-pilot eigenvalue upper bounds (per-draw rate maxima)
#' @noRd
.two_block_pilot_eigenvalue_ub <- function(
    b_draws,
    n_pilot,
    J,
    p_re,
    re_names,
    group_levels,
    x,
    block,
    x_hyper,
    prior_list,
    pfamily_list,
    family,
    dispersion_block1,
    tv_tol) {

  b_arr <- array(as.numeric(b_draws), dim = c(J, p_re, n_pilot))
  lambda_star_vec <- numeric(n_pilot)
  max_eigenvalues <- NULL
  rate_upper <- NULL
  lambda_star_best <- -Inf
  i_max_rate <- NA_integer_

  for (i in seq_len(n_pilot)) {
    b_i <- b_arr[, , i, drop = FALSE]
    dim(b_i) <- c(J, p_re)
    dimnames(b_i) <- list(group_levels, re_names)
    mode_w_i <- two_block_mode_weights(
      x            = x,
      block        = block,
      b_mode       = b_i,
      family       = family,
      dispersion   = dispersion_block1,
      group_levels = group_levels
    )
    rate_i <- two_block_rate_v2(
      x                 = x,
      block             = block,
      x_hyper           = x_hyper,
      prior_list_block1 = prior_list,
      pfamily_list      = pfamily_list,
      weights           = mode_w_i$weights,
      family            = family,
      group_levels      = group_levels
    )
    if (is.null(max_eigenvalues)) {
      max_eigenvalues <- rep(-Inf, length(rate_i$eigenvalues))
    }
    lambda_star_vec[i] <- rate_i$lambda_star
    max_eigenvalues <- pmax(max_eigenvalues, rate_i$eigenvalues)
    if (rate_i$lambda_star > lambda_star_best) {
      lambda_star_best <- rate_i$lambda_star
      rate_upper <- rate_i
      i_max_rate <- i
    }
  }

  rate_upper_eig <- rate_upper
  rate_upper_eig$eigenvalues <- max_eigenvalues
  rate_upper_eig$lambda_star <- max_eigenvalues[1L]

  m_min_upper <- two_block_l_for_tv(
    rate_upper_eig, tv_tol, method = "theorem3"
  ) + 1L

  list(
    rate_upper        = rate_upper_eig,
    m_min_upper       = m_min_upper,
    lambda_star_vec   = lambda_star_vec,
    i_max_rate        = i_max_rate,
    max_eigenvalues   = max_eigenvalues
  )
}
