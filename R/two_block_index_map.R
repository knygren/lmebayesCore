#' Integer index map from factor levels and design names
#'
#' Built once, before the pilot. Row and column positions used inside the
#' sweep are the result of these name subscripts.
#'
#' @param design List with \code{D}, \code{group}, and \code{W}.
#' @param re_names Character vector of random-effect names, in the column
#'   order of the stored \code{b} array.
#' @param group_levels Character vector of group levels, in the row order of
#'   the stored \code{b} array.
#' @param fixef_start Named list of hyper-parameter vectors. Names of each
#'   vector are the column order of the stored \code{fixef} matrix.
#' @return A list of character keys and integer index vectors.
#' @noRd
.two_block_index_map <- function(design, re_names, group_levels, fixef_start) {
  if (!is.factor(design$group)) {
    stop("'design$group' must be a factor.", call. = FALSE)
  }
  group_levels <- as.character(group_levels)
  re_names <- as.character(re_names)
  if (length(group_levels) < 1L || anyDuplicated(group_levels)) {
    stop("'group_levels' must be a non-empty vector of unique names.", call. = FALSE)
  }
  if (length(re_names) < 1L || anyDuplicated(re_names)) {
    stop("'re_names' must be a non-empty vector of unique names.", call. = FALSE)
  }
  if (!identical(group_levels, levels(design$group))) {
    stop(
      "'group_levels' must be identical to levels(design$group).",
      call. = FALSE
    )
  }

  D <- as.matrix(design$D)
  d_names <- colnames(D)
  if (is.null(d_names) || anyDuplicated(d_names)) {
    stop("'D' must have unique column names.", call. = FALSE)
  }
  d_col <- match(re_names, d_names)
  if (anyNA(d_col)) {
    stop(
      "Random-effect name(s) missing from colnames(D): ",
      paste(re_names[is.na(d_col)], collapse = ", "),
      call. = FALSE
    )
  }

  g_chr <- as.character(design$group)
  if (length(g_chr) != nrow(D)) {
    stop("length(design$group) must equal nrow(D).", call. = FALSE)
  }
  obs_rows <- lapply(group_levels, function(g) which(g_chr == g))
  names(obs_rows) <- group_levels
  group_code <- integer(length(g_chr))
  for (j in seq_along(group_levels)) {
    group_code[obs_rows[[j]]] <- j
  }
  if (any(group_code < 1L)) {
    stop("An observation has a group label outside group_levels.", call. = FALSE)
  }
  if (!identical(group_levels[group_code], g_chr)) {
    stop("group_code does not rebuild as.character(design$group).", call. = FALSE)
  }
  covered <- sort(unlist(obs_rows, use.names = FALSE))
  if (!identical(covered, seq_len(nrow(D)))) {
    stop("obs_rows does not cover each row of D exactly once.", call. = FALSE)
  }

  if (!is.list(fixef_start) || is.null(names(fixef_start))) {
    stop("'fixef_start' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(fixef_start), re_names)) {
    stop("names(fixef_start) must match re_names.", call. = FALSE)
  }

  w_for_group <- vector("list", length(re_names))
  b_for_wrow <- vector("list", length(re_names))
  w_col <- vector("list", length(re_names))
  names(w_for_group) <- names(b_for_wrow) <- names(w_col) <- re_names

  for (k in re_names) {
    Xk <- as.matrix(design$W[[k]])
    rn <- rownames(Xk)
    cn <- colnames(Xk)
    if (is.null(rn) || anyDuplicated(rn)) {
      stop("W[[", k, "]] must have unique rownames.", call. = FALSE)
    }
    if (is.null(cn) || anyDuplicated(cn)) {
      stop("W[[", k, "]] must have unique colnames.", call. = FALSE)
    }
    miss_g <- setdiff(group_levels, rn)
    miss_w <- setdiff(rn, group_levels)
    if (length(miss_g) > 0L || length(miss_w) > 0L) {
      stop(
        "rownames(W[[", k, "]]) and group_levels are not the same set.",
        call. = FALSE
      )
    }
    w_for_group[[k]] <- match(group_levels, rn)
    b_for_wrow[[k]] <- match(rn, group_levels)
    gamma_names <- names(fixef_start[[k]])
    if (is.null(gamma_names) || anyDuplicated(gamma_names)) {
      stop("fixef_start[[", k, "]] must have unique names.", call. = FALSE)
    }
    w_col[[k]] <- match(cn, gamma_names)
    if (anyNA(w_col[[k]]) || anyNA(w_for_group[[k]]) || anyNA(b_for_wrow[[k]])) {
      stop("Name match failed for W[[", k, "]].", call. = FALSE)
    }
  }

  list(
    group_levels = group_levels,
    re_names     = re_names,
    obs_rows     = obs_rows,
    group_code   = group_code,
    d_col        = d_col,
    w_for_group  = w_for_group,
    b_for_wrow   = b_for_wrow,
    w_col        = w_col
  )
}

#' Factor whose codes are the stored row positions of \code{b}
#'
#' Code \code{j} is \code{group_levels[j]}. This is the factor passed to the
#' existing group sampler.
#'
#' @param map Result of \code{.two_block_index_map()}.
#' @return A factor of length \code{nrow(D)}.
#' @noRd
.two_block_index_group_factor <- function(map) {
  factor(
    map$group_code,
    levels = seq_along(map$group_levels),
    labels = map$group_levels
  )
}

#' Prior-mean matrix from the integer row map
#'
#' Column \code{j} is the prior mean for \code{group_levels[j]}.
#'
#' @param map Result of \code{.two_block_index_map()}.
#' @param fixef Named list of hyper-parameter vectors.
#' @param design List with \code{W}.
#' @return Numeric matrix \code{p_re x J}.
#' @noRd
.two_block_index_mu_all <- function(map, fixef, design) {
  p_re <- length(map$re_names)
  J <- length(map$group_levels)
  mu <- matrix(
    NA_real_,
    nrow = p_re,
    ncol = J,
    dimnames = list(map$re_names, map$group_levels)
  )
  for (i in seq_along(map$re_names)) {
    k <- map$re_names[i]
    Xk <- as.matrix(design$W[[k]])
    gamma <- fixef[[k]]
    if (is.matrix(gamma)) {
      stop("fixef[[", k, "]] must be a vector.", call. = FALSE)
    }
    g_w <- as.numeric(gamma)[map$w_col[[k]]]
    rows <- map$w_for_group[[k]]
    for (j in seq_len(J)) {
      mu[i, j] <- sum(Xk[rows[j], , drop = TRUE] * g_w)
    }
  }
  mu
}
