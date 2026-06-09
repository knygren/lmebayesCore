#' Two-block Gaussian Gibbs sampler (Gaussian-only backup)
#'
#' Frozen copy of \code{\link{two_block_rNormal_reg}} before the
#' \code{family}-aware Block~1 path was added.  Block~1 always uses
#' \code{\link{block_rNormalReg}}; Block~2 uses \code{\link{multi_rNormal_reg}}
#' with \code{family = gaussian()}.
#'
#' @inheritParams two_block_rNormal_reg
#' @return Object of class \code{"two_block_rNormal_reg_v2"}.
#' @family simfuncs
#' @seealso \code{\link{two_block_rNormal_reg}}
#' @export
two_block_rNormal_reg_v2 <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list_block1,
    prior_list_block2,
    fixef_start,
    re_coef_names = colnames(x),
    group_levels = levels(block),
    group_name = NULL,
    m_convergence = 10L,
    sampling = c("replicate", "chain"),
    seed = NULL,
    progbar = TRUE) {

  cl <- match.call()
  sampling <- match.arg(sampling)
  if (!identical(sampling, "replicate")) {
    stop("Only sampling = \"replicate\" is implemented.", call. = FALSE)
  }

  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  m_convergence <- as.integer(m_convergence[1L])
  if (m_convergence < 1L) {
    stop("'m_convergence' must be at least 1.", call. = FALSE)
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
      if (is.null(cn) || length(cn) != ncol(x)) paste0("RE", seq_len(ncol(x))) else cn
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

  if (!is.list(prior_list_block2)) {
    stop("'prior_list_block2' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(prior_list_block2), re_names)) {
    stop(
      "names(prior_list_block2) must match re_coef_names.",
      call. = FALSE
    )
  }

  if (!is.list(fixef_start) || is.null(names(fixef_start))) {
    stop("'fixef_start' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(fixef_start), re_names)) {
    stop("names(fixef_start) must match re_coef_names.", call. = FALSE)
  }
  fixef_start <- fixef_start[re_names]

  if (!is.list(prior_list_block1)) {
    stop("'prior_list_block1' must be a list.", call. = FALSE)
  }
  P <- prior_list_block1$P
  if (is.null(P)) {
    Sigma_b <- prior_list_block1$Sigma
    if (is.null(Sigma_b)) {
      stop("prior_list_block1 must contain 'P' or 'Sigma'.", call. = FALSE)
    }
    P <- solve(Sigma_b)
  }
  dispersion <- prior_list_block1$dispersion
  if (is.null(dispersion)) {
    stop("prior_list_block1 must contain 'dispersion'.", call. = FALSE)
  }
  ddef <- if (is.null(prior_list_block1$ddef)) FALSE else prior_list_block1$ddef

  if (!is.null(seed)) {
    set.seed(seed)
  }

  fixef <- fixef_start
  mu_all <- .two_block_mu_all_v2(fixef, x_hyper, re_names, group_levels)
  block1_args <- list(
    n          = 1L,
    y          = y,
    x          = x,
    block      = block,
    prior_list = list(
      mu         = mu_all,
      P          = P,
      dispersion = dispersion,
      ddef       = ddef
    )
  )

  coef_cols <- c("draw", group_name, re_names)
  draw_rows <- vector("list", n)

  fixef_draws <- stats::setNames(
    lapply(re_names, function(k) {
      q_k <- length(fixef_start[[k]])
      matrix(NA_real_, nrow = n, ncol = q_k,
             dimnames = list(NULL, names(fixef_start[[k]])))
    }),
    re_names
  )

  if (isTRUE(progbar)) {
    pb <- utils::txtProgressBar(min = 0L, max = n, style = 3L)
    on.exit(close(pb), add = TRUE)
  }

  b_i <- NULL

  for (i in seq_len(n)) {
    if (isTRUE(progbar)) {
      utils::setTxtProgressBar(pb, i)
    }

    fixef <- fixef_start

    for (m in seq_len(m_convergence)) {

      mu_all <- .two_block_mu_all_v2(fixef, x_hyper, re_names, group_levels)
      block1_args$prior_list$mu <- mu_all
      block_i <- do.call(block_rNormalReg, block1_args)
      b_i <- block_i$coefficients
      if (is.null(rownames(b_i))) {
        rownames(b_i) <- block_i$block_info$ids
      }
      colnames(b_i) <- re_names

      fixef_draw <- multi_rNormal_reg(
        n          = 1L,
        y          = b_i,
        x          = x_hyper,
        prior_list = prior_list_block2,
        progbar    = FALSE
      )
      fixef <- stats::setNames(
        lapply(re_names, function(k) fixef_draw[[k]]$coefficients[1L, ]),
        re_names
      )
    }

    for (k in re_names) {
      fixef_draws[[k]][i, ] <- fixef[[k]]
    }

    J_i <- nrow(b_i)
    draw_df <- data.frame(
      draw = rep(i, J_i),
      stringsAsFactors = FALSE
    )
    draw_df[[group_name]] <- rownames(b_i)
    for (nm in re_names) {
      draw_df[[nm]] <- b_i[, nm]
    }
    draw_rows[[i]] <- draw_df
  }

  coefficients <- do.call(rbind, draw_rows)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, coef_cols, drop = FALSE]

  structure(
    list(
      fixef_draws   = fixef_draws,
      coefficients  = coefficients,
      fixef_last    = fixef,
      b_last        = b_i,
      mu_all_last   = mu_all,
      n             = n,
      m_convergence = m_convergence,
      sampling      = sampling,
      fixef_start   = fixef_start,
      re_coef_names = re_names,
      group_levels  = group_levels,
      group_name    = group_name,
      call          = cl
    ),
    class = "two_block_rNormal_reg_v2"
  )
}

#' @noRd
.two_block_mu_all_v2 <- function(fixef, x_hyper, re_names, group_levels) {
  p_re <- length(re_names)
  J    <- length(group_levels)
  mu_all <- matrix(NA_real_, nrow = p_re, ncol = J,
                   dimnames = list(re_names, group_levels))
  for (i in seq_len(p_re)) {
    k       <- re_names[i]
    gamma_k <- fixef[[k]]
    X_k     <- as.matrix(x_hyper[[k]])
    rn      <- rownames(X_k)
    if (is.null(rn)) {
      if (nrow(X_k) != J) {
        stop(
          "nrow(x_hyper[[", k, "]]) must equal length(group_levels).",
          call. = FALSE
        )
      }
      for (j in seq_len(J)) {
        mu_all[i, j] <- sum(X_k[j, , drop = TRUE] * gamma_k)
      }
    } else {
      for (j in seq_len(J)) {
        mu_all[i, j] <- sum(X_k[group_levels[j], , drop = TRUE] * gamma_k)
      }
    }
  }
  mu_all
}
