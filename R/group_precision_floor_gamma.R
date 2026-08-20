## Step 1: integrate population gamma -> Lambda_beta, mu_beta on stacked beta.

#' @noRd
.group_floor_idx <- function(j, p_re) {
  ((j - 1L) * p_re + 1L):(j * p_re)
}

#' @noRd
.group_floor_stack_beta <- function(b_mat) {
  J <- nrow(b_mat)
  p_re <- ncol(b_mat)
  out <- numeric(J * p_re)
  for (j in seq_len(J)) {
    out[.group_floor_idx(j, p_re)] <- b_mat[j, ]
  }
  out
}

#' @noRd
.group_floor_unstack_beta <- function(b, J, p_re) {
  matrix(b, nrow = p_re, ncol = J, byrow = FALSE)
}

#' Integrate gamma out of the Gaussian prior; return Lambda_beta and mu_beta.
#' @noRd
.group_floor_integrate_gamma <- function(design,
                                        prior_pack,
                                        group_levels) {
  re_names <- design$groupef.names
  p_re <- length(re_names)
  J <- length(group_levels)
  n_dim <- J * p_re
  pop <- prior_pack$pop.prior_list

  P_b <- solve(as.matrix(prior_pack$group.Sigma))
  P_b <- 0.5 * (P_b + t(P_b))

  X_hyper <- stats::setNames(
    lapply(re_names, function(k) {
      X_k <- as.matrix(design$W[[k]])
      rn <- rownames(X_k)
      if (!is.null(rn)) {
        if (!all(group_levels %in% rn)) {
          stop("design$W[[\"", k, "\"]] is missing rows for some group levels.",
               call. = FALSE)
        }
        X_k <- X_k[group_levels, , drop = FALSE]
      } else if (nrow(X_k) != J) {
        stop("nrow(design$W[[\"", k, "\"]]) must equal the number of groups.",
             call. = FALSE)
      }
      X_k
    }),
    re_names
  )

  q_k <- vapply(X_hyper, ncol, integer(1L))
  q <- sum(q_k)
  gamma_cols <- split(seq_len(q), rep(seq_len(p_re), q_k))

  P11 <- matrix(0, q, q)
  for (k in seq_len(p_re)) {
    Sigma_k <- as.matrix(pop[[k]]$Sigma)
    Vk <- chol2inv(chol(Sigma_k))
    P11[gamma_cols[[k]], gamma_cols[[k]]] <- 0.5 * (Vk + t(Vk))
  }

  P22 <- matrix(0, n_dim, n_dim)
  P12 <- matrix(0, q, n_dim)

  for (j in seq_len(J)) {
    bc <- .group_floor_idx(j, p_re)
    P22[bc, bc] <- P_b

    H_j <- matrix(0, p_re, q)
    for (k in seq_len(p_re)) {
      H_j[k, gamma_cols[[k]]] <- X_hyper[[k]][j, ]
    }
    P12[, bc] <- -t(H_j) %*% P_b
    P11 <- P11 + t(H_j) %*% P_b %*% H_j
  }
  P11 <- 0.5 * (P11 + t(P11))
  P22 <- 0.5 * (P22 + t(P22))

  P21 <- t(P12)
  Lambda_beta <- P22 - P21 %*% solve(P11, P12)
  Lambda_beta <- 0.5 * (Lambda_beta + t(Lambda_beta))

  mu_beta <- numeric(n_dim)
  for (j in seq_len(J)) {
    mu_bj <- numeric(p_re)
    for (k in seq_len(p_re)) {
      mu_bj[k] <- sum(X_hyper[[k]][j, ] * pop[[k]]$mu)
    }
    mu_beta[.group_floor_idx(j, p_re)] <- mu_bj
  }

  list(
    Lambda_beta = Lambda_beta,
    mu_beta = mu_beta,
    P_b = P_b,
    P11 = P11,
    P12 = P12,
    X_hyper = X_hyper,
    gamma_cols = gamma_cols,
    re_names = re_names,
    p_re = p_re,
    J = J,
    n_dim = n_dim,
    q = q
  )
}
