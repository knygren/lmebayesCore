## Floor coupling eigenvalues kappa_i^LB from floor blocks (category 6).

#' @noRd
.c05_floor_coupling_S <- function(p11, Gamma_lb, P_b = NULL) {
  if (is.null(p11$P11) || is.null(p11$H_list)) {
    stop("'p11' must contain 'P11' and 'H_list'.", call. = FALSE)
  }
  if (is.null(Gamma_lb) || !is.list(Gamma_lb)) {
    stop("'Gamma_lb' must be a named list of per-group matrices.", call. = FALSE)
  }
  if (is.null(P_b)) {
    P_b <- p11$P_b
  }
  if (is.null(P_b)) {
    stop("'P_b' is required.", call. = FALSE)
  }

  q <- nrow(p11$P11)
  S <- matrix(0, q, q)
  group_levels <- names(p11$H_list)
  if (is.null(group_levels)) {
    group_levels <- names(Gamma_lb)
  }

  for (lev in group_levels) {
    H_j <- p11$H_list[[lev]]
    Gamma_j <- Gamma_lb[[lev]]
    if (is.null(H_j) || is.null(Gamma_j)) {
      stop("Missing H_j or Gamma_lb entry for group '", lev, "'.", call. = FALSE)
    }
    B_j <- Gamma_j + P_b
    B_j <- 0.5 * (B_j + t(B_j))
    C_j <- P_b %*% solve(B_j, P_b)
    C_j <- 0.5 * (C_j + t(C_j))

    p_re <- nrow(H_j)
    x_j <- lapply(seq_len(p_re), function(k) {
      cols <- p11$gamma_cols[[k]]
      H_j[k, cols]
    })

    for (i in seq_len(p_re)) {
      for (k in i:p_re) {
        out_ik <- outer(x_j[[i]], x_j[[k]])
        S[p11$gamma_cols[[i]], p11$gamma_cols[[k]]] <-
          S[p11$gamma_cols[[i]], p11$gamma_cols[[k]]] + C_j[i, k] * out_ik
        if (k > i) {
          S[p11$gamma_cols[[k]], p11$gamma_cols[[i]]] <-
            t(S[p11$gamma_cols[[i]], p11$gamma_cols[[k]]])
        }
      }
    }
  }

  0.5 * (S + t(S))
}

#' Floor coupling eigenvalues from data-precision lower bounds.
#'
#' Forms \deqn{S^{\mathrm{LB}} = \sum_j H_j^\top P_b (P_{22,j}^{\mathrm{LB}})^{-1} P_b H_j}
#' and returns \eqn{\kappa_i^{\mathrm{LB}} = \mathrm{eig}(P_{11}^{-1/2} S^{\mathrm{LB}} P_{11}^{-1/2})}.
#' Rosenthal drift constants are computed separately by
#' \code{\link{rosenthal_drift_constants}}.
#'
#' @param mode A \code{\link{population_mode}} result with \code{p11}.
#' @param floor A \code{\link{group_precision_floor}} result with \code{Gamma_lb}.
#' @return A list with \code{kappa_lb}, \code{kappa_max_lb}, \code{S_lb}, and
#'   \code{q}.
#' @seealso \code{\link{rosenthal_drift_constants}}, \code{\link{rosenthal_tv_bound}}
#' @export
floor_coupling_eigenvalues <- function(mode, floor) {
  if (is.null(mode$p11)) {
    stop("'mode' must be a population_mode() result with 'p11'.", call. = FALSE)
  }
  if (is.null(floor$Gamma_lb)) {
    stop("'floor' must be a group_precision_floor() result with 'Gamma_lb'.",
         call. = FALSE)
  }

  p11 <- mode$p11
  P_b <- floor$P_b
  if (is.null(P_b)) {
    P_b <- p11$P_b
  }
  S_lb <- .c05_floor_coupling_S(p11, floor$Gamma_lb, P_b = P_b)
  kappa_lb <- .two_block_gen_eigen(S_lb, p11$P11, strict = FALSE)
  kappa_max_lb <- max(kappa_lb)

  list(
    kappa_lb = kappa_lb,
    kappa_max_lb = kappa_max_lb,
    S_lb = S_lb,
    q = length(kappa_lb),
    method = "floor_coupling_eigenvalues"
  )
}

#' Backward-compatible wrapper combining eigenvalues and drift constants.
#'
#' @param mode A \code{\link{population_mode}} result.
#' @param beta_set A \code{\link{group_precision_floor}} or legacy object with
#'   \code{Gamma_lb}.
#' @return List from \code{\link{floor_coupling_eigenvalues}} plus Rosenthal
#'   drift fields from \code{\link{rosenthal_drift_constants}}.
#' @export
floor_coupling_spectrum <- function(mode, beta_set) {
  ev <- floor_coupling_eigenvalues(mode, beta_set)
  drift <- rosenthal_drift_constants(ev)
  c(ev, drift[setdiff(names(drift), names(ev))])
}
