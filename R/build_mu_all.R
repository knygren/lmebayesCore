#' Build per-group random-effect prior means for Block 1 sampling
#'
#' Forms the \code{mu_all} matrix passed to
#' \code{\link{block_rNormalReg_update}}: for each grouping level
#' \eqn{j} and each random-effect column \eqn{k} of the level-1 design
#' \code{Z},
#' \deqn{\mu\_\text{all}[k, j] = X_{\text{hyper},k}[j,]^\top \gamma_k,}
#' where \eqn{\gamma_k} is the current hyper-parameter vector for RE \eqn{k}
#' (Block 2 state).
#'
#' @param design List with components \code{X_hyper}, \code{re_coef_names},
#'   and \code{groups} (typically supplied by a downstream mixed-effects
#'   model setup step).
#' @param fixef Named list of hyper-parameter vectors, one entry per RE column
#'   of \code{Z}. Names must match \code{design$re_coef_names}. Each
#'   \code{fixef[[k]]} is a numeric vector of length \code{ncol(X_hyper[[k]])}
#'   with names matching \code{colnames(X_hyper[[k]])}.
#' @param group_levels Character vector of grouping levels defining the
#'   \emph{column order} of \code{mu_all}.  Defaults to
#'   \code{levels(design$groups)}, which is the canonical ordering used in
#'   two-block mixed models and consistent with \code{lmer} and
#'   \code{\link{block_rNormalReg}} (which preserves the input factor's level
#'   order).
#' @param use_cpp Logical; if \code{TRUE} (default), compute \code{mu_all} in
#'   C++ via \code{.two_block_build_mu_all_cpp}; if \code{FALSE}, use
#'   the R reference implementation \code{build_mu_all_r()}.
#' @return A list with:
#'   \describe{
#'     \item{\code{mu_all}}{Numeric matrix \code{p_re x J} (\code{p_re} =
#'       number of RE columns, \code{J} = number of groups). Row names are
#'       \code{design$re_coef_names}; column names are \code{group_levels} in
#'       the order supplied.}
#'     \item{\code{re_coef_names}}{Copy of \code{design$re_coef_names}.}
#'     \item{\code{group_levels}}{Grouping levels used for columns.}
#'   }
#' @seealso \code{\link{lmerb_posterior_mean}},
#'   \code{\link{two_block_rNormal_reg_v2}},
#'   \code{\link{block_rNormalReg_update}}
#' @export
build_mu_all <- function(design, fixef, group_levels = NULL, use_cpp = TRUE) {
  if (isTRUE(use_cpp)) {
    .lmerb_validate_design(design)
    if (is.null(group_levels)) {
      group_levels <- levels(design$groups)
    } else {
      group_levels <- as.character(group_levels)
    }
    x_hyper <- lapply(design$X_hyper, as.matrix)
    mu_all <- .two_block_build_mu_all_cpp(
      x_hyper, fixef, design$re_coef_names, group_levels
    )
    return(list(
      mu_all        = mu_all,
      re_coef_names = design$re_coef_names,
      group_levels  = group_levels
    ))
  }
  build_mu_all_r(design, fixef, group_levels)
}

#' Build per-group random-effect prior means (R reference implementation)
#' @noRd
build_mu_all_r <- function(design, fixef, group_levels = NULL) {

  .lmerb_validate_design(design)

  if (!is.list(fixef) || is.null(names(fixef))) {
    stop("'fixef' must be a named list with one element per RE column.", call. = FALSE)
  }

  re <- design$re_coef_names
  Xh <- design$X_hyper

  if (length(re) < 1L) {
    stop("'design' has no random-effect columns.", call. = FALSE)
  }
  if (length(fixef) != length(re)) {
    stop(
      "length(fixef) (", length(fixef), ") must equal number of RE columns (",
      length(re), ").", call. = FALSE
    )
  }
  if (!setequal(names(fixef), re)) {
    stop(
      "names(fixef) must match design$re_coef_names: ",
      paste(re, collapse = ", "),
      call. = FALSE
    )
  }
  fixef <- fixef[re]

  if (is.null(group_levels)) {
    group_levels <- levels(design$groups)
  } else {
    group_levels <- as.character(group_levels)
  }
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  p_re <- length(re)
  J    <- length(group_levels)
  mu_all <- matrix(NA_real_, nrow = p_re, ncol = J,
                   dimnames = list(re, group_levels))

  for (i in seq_len(p_re)) {
    k       <- re[i]
    gamma_k <- fixef[[k]]
    X_k     <- Xh[[k]]

    if (is.null(X_k)) {
      stop("design$X_hyper[[", k, "]] is missing.", call. = FALSE)
    }
    if (!is.numeric(gamma_k) || length(gamma_k) < 1L) {
      stop("fixef[[", k, "]] must be a numeric vector.", call. = FALSE)
    }
    p_k <- ncol(X_k)
    if (length(gamma_k) != p_k) {
      stop(
        "length(fixef[[", k, "]]) (", length(gamma_k),
        ") must equal ncol(X_hyper[[", k, "]]) (", p_k, ").",
        call. = FALSE
      )
    }
    cn <- colnames(X_k)
    if (!is.null(cn)) {
      if (is.null(names(gamma_k)) || any(names(gamma_k) == "")) {
        stop("fixef[[", k, "]] must be a named vector.", call. = FALSE)
      }
      if (!identical(names(gamma_k), cn)) {
        stop(
          "names(fixef[[", k, "]]) must match colnames(X_hyper[[", k, "]]).",
          call. = FALSE
        )
      }
    }

    rn <- rownames(X_k)
    if (is.null(rn)) {
      if (nrow(X_k) != J) {
        stop(
          "nrow(X_hyper[[", k, "]]) (", nrow(X_k),
          ") must equal length(group_levels) (", J, ").",
          call. = FALSE
        )
      }
      for (j in seq_len(J)) {
        mu_all[i, j] <- sum(X_k[j, , drop = TRUE] * gamma_k)
      }
    } else {
      miss <- setdiff(group_levels, rn)
      if (length(miss) > 0L) {
        stop(
          "group level(s) not found in rownames(X_hyper[[", k, "]]): ",
          paste(miss, collapse = ", "),
          call. = FALSE
        )
      }
      for (j in seq_len(J)) {
        mu_all[i, j] <- sum(X_k[group_levels[j], , drop = TRUE] * gamma_k)
      }
    }
  }

  list(
    mu_all        = mu_all,
    re_coef_names = re,
    group_levels  = group_levels
  )
}

#' @noRd
.lmerb_validate_design <- function(design) {
  if (!is.list(design)) {
    stop("'design' must be a list.", call. = FALSE)
  }
  for (nm in c("X_hyper", "re_coef_names", "groups")) {
    if (is.null(design[[nm]])) {
      stop("'design' must contain '", nm, "'.", call. = FALSE)
    }
  }
  invisible(design)
}
