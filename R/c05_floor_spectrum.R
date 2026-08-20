## Floor coupling spectrum kappa_i^LB(delta_2) on beta safe set B-tilde.

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

#' Floor coupling eigenvalue spectrum on a beta safe set.
#'
#' Builds block-diagonal floor precisions
#' \eqn{P_{22,j}^{\mathrm{LB}} = P_b + \underline{\mathcal P}_{j,\mathrm{data}}}
#' from \code{\link{beta_marginal_safe_set}} (or compatible floor output), forms
#' \deqn{S^{\mathrm{LB}} = \sum_j H_j^\top P_b (P_{22,j}^{\mathrm{LB}})^{-1} P_b H_j,}
#' and returns \eqn{\kappa_i^{\mathrm{LB}} = \mathrm{eig}(P_{11}^{-1/2} S^{\mathrm{LB}} P_{11}^{-1/2})}
#' with Foster/Rosenthal drift constants.
#'
#' @param mode A \code{\link{population_mode}} result with \code{p11}.
#' @param beta_set A \code{\link{beta_marginal_safe_set}} or
#'   \code{\link{group_precision_floor}} result with \code{Gamma_lb} and
#'   \code{P_b}.
#' @return A list with \code{kappa_lb}, \code{kappa_max_lb}, \code{lambda_lb},
#'   \code{C_beta_plus}, \code{b_drift}, \code{S_lb}, and \code{q}.
#' @seealso \code{\link{deficiency_spectrum}}, \code{\link{rosenthal_tv_bound}}
#' @export
floor_coupling_spectrum <- function(mode, beta_set) {
  if (is.null(mode$p11)) {
    stop("'mode' must be a population_mode() result with 'p11'.", call. = FALSE)
  }
  if (is.null(beta_set$Gamma_lb)) {
    stop("'beta_set' must contain 'Gamma_lb'.", call. = FALSE)
  }

  p11 <- mode$p11
  P_b <- beta_set$P_b
  if (is.null(P_b)) {
    P_b <- p11$P_b
  }
  S_lb <- .c05_floor_coupling_S(p11, beta_set$Gamma_lb, P_b = P_b)
  kappa_lb <- .two_block_gen_eigen(S_lb, p11$P11, strict = FALSE)
  kappa_max_lb <- max(kappa_lb)
  lambda_lb <- kappa_max_lb^2
  q <- length(kappa_lb)
  C_beta_plus <- 0.5 * sum(kappa_lb)
  b_drift <- 1 - lambda_lb + q / 2 + C_beta_plus

  list(
    kappa_lb = kappa_lb,
    kappa_max_lb = kappa_max_lb,
    lambda_lb = lambda_lb,
    C_beta_plus = C_beta_plus,
    b_drift = b_drift,
    S_lb = S_lb,
    q = q,
    method = "floor_coupling"
  )
}
