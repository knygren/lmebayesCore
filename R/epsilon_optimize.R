## Stage 2: numerical epsilon(gamma*) via the convex profile (Definition 5).



#' Extract one posterior draw from an \code{rglmb()} fit.
#' @noRd
.minorization_rglmb_draw <- function(fit) {
  if (!is.null(fit$coefficients)) {
    as.numeric(fit$coefficients[1L, , drop = TRUE])
  } else if (!is.null(fit$coef)) {
    as.numeric(fit$coef)
  } else {
    stop("rglmb() result is missing coefficient draws.", call. = FALSE)
  }
}

#' @noRd
.minorization_rglmb_mode <- function(fit) {
  if (!is.null(fit$coef.mode)) {
    as.numeric(fit$coef.mode)
  } else if (!is.null(fit$fit) && !is.null(fit$fit$coefficients)) {
    as.numeric(fit$fit$coefficients)
  } else {
    stop("rglmb() result is missing coef.mode.", call. = FALSE)
  }
}



#' @noRd

.minorization_m_beta <- function(b_mat, p11) {

  rhs <- p11$Lambda_gamma %*% p11$mu_0

  for (j in seq_len(p11$J)) {

    lev <- p11$group_levels[[j]]

    H_j <- p11$H_list[[lev]]

    b_j <- b_mat[j, ]

    rhs <- rhs + t(H_j) %*% p11$P_b %*% b_j

  }

  as.vector(solve(p11$P11, rhs))

}



#' @noRd

.minorization_log_phi_q <- function(gamma_prime, mean_gamma, p11) {

  d <- gamma_prime - mean_gamma

  q <- p11$q

  mah <- as.numeric(crossprod(d, p11$P11 %*% d))

  -(q / 2) * log(2 * pi) +

    0.5 * as.numeric(determinant(p11$P11, logarithm = TRUE)$modulus) -

    0.5 * mah

}



#' @noRd

.minorization_log_q_Q <- function(gamma_prime, gamma_star, p11) {

  d <- gamma_prime - gamma_star

  q <- p11$q

  mah <- as.numeric(crossprod(d, p11$P11 %*% d))

  -(q / 2) * log(2 * pi) +

    0.5 * as.numeric(determinant(p11$P11, logarithm = TRUE)$modulus) -

    0.5 * mah

}



#' @noRd

.minorization_logsumexp <- function(x) {

  x <- x[is.finite(x)]

  if (!length(x)) {

    return(-Inf)

  }

  m <- max(x)

  m + log(sum(exp(x - m)))

}



#' @noRd

.minorization_draw_beta <- function(design,

                                    fixef,

                                    p11,

                                    family,

                                    measurement_prior_list,

                                    n,

                                    mc_seed = NULL) {

  if (!is.null(mc_seed)) {

    set.seed(mc_seed)

  }



  group_levels <- p11$group_levels

  J <- p11$J

  p_re <- p11$p_re

  re_names <- p11$re_names

  g_chr <- as.character(design$group)

  Sigma_b <- measurement_prior_list$group.Sigma

  sigma2 <- measurement_prior_list$group.dispersion

  is_gaussian <- identical(family$family, "gaussian")

  mu_all <- as.matrix(

    build_mu_all(design, fixef, group_levels = group_levels)$mu_all

  )



  draws <- vector("list", n)

  if (is_gaussian) {

    if (is.null(sigma2)) {

      stop(

        "'measurement_prior_list$group.dispersion' is required for gaussian().",

        call. = FALSE

      )

    }

    sigma2 <- as.numeric(sigma2)

    if (!(length(sigma2) %in% c(1L, J))) {

      stop(

        "'measurement_prior_list$group.dispersion' must have length 1 or J.",

        call. = FALSE

      )

    }

    P_b <- p11$P_b

    for (m in seq_len(n)) {

      b_mat <- matrix(

        0, nrow = J, ncol = p_re,

        dimnames = list(group_levels, re_names)

      )

      for (jj in seq_len(J)) {

        lev <- group_levels[[jj]]

        rows <- which(g_chr == lev)

        Z_j <- design$D[rows, , drop = FALSE]

        y_j <- design$y[rows]

        sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2

        mu_j <- mu_all[, jj]

        post_P_j <- crossprod(Z_j) / sigma2_j + P_b

        post_v_j <- crossprod(Z_j, y_j) / sigma2_j + P_b %*% mu_j

        post_mean_j <- solve(post_P_j, post_v_j)

        post_cov_j <- solve(post_P_j)

        b_mat[jj, ] <- MASS::mvrnorm(1L, post_mean_j, post_cov_j)

      }

      draws[[m]] <- b_mat

    }

    return(draws)

  }

  for (m in seq_len(n)) {

    b_mat <- matrix(

      0, nrow = J, ncol = p_re,

      dimnames = list(group_levels, re_names)

    )

    for (jj in seq_len(J)) {

      lev <- group_levels[[jj]]

      rows <- which(g_chr == lev)

      y_j <- design$y[rows]

      Z_j <- design$D[rows, , drop = FALSE]

      mu_j <- mu_all[, jj]

      pf_j <- if (is_gaussian) {

        sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2

        glmbayesCore::dNormal(

          mu = mu_j,

          Sigma = Sigma_b,

          dispersion = sigma2_j

        )

      } else {

        glmbayesCore::dNormal(mu = mu_j, Sigma = Sigma_b)

      }

      fit_j <- glmbayesCore::rglmb(

        n = 1L,

        y = y_j,

        x = Z_j,

        family = family,

        pfamily = pf_j,

        verbose = FALSE

      )

      b_mat[jj, ] <- .minorization_rglmb_draw(fit_j)

    }

    draws[[m]] <- b_mat

  }

  draws

}



#' @noRd

.minorization_log_kernel_mc <- function(gamma_prime, beta_draws, p11) {

  log_d <- vapply(

    beta_draws,

    function(b_mat) {

      m_beta <- .minorization_m_beta(b_mat, p11)

      .minorization_log_phi_q(gamma_prime, m_beta, p11)

    },

    numeric(1L)

  )

  .minorization_logsumexp(log_d) - log(length(log_d))

}



#' @noRd

.minorization_g_mc <- function(gamma_prime,

                               gamma_star,

                               beta_draws,

                               p11) {

  log_q <- .minorization_log_kernel_mc(gamma_prime, beta_draws, p11)

  log_q_Q <- .minorization_log_q_Q(gamma_prime, gamma_star, p11)

  log_q - log_q_Q

}



#' @noRd

.minorization_log_gaussian <- function(gamma_prime, mean_gamma, cov) {

  q <- length(gamma_prime)

  d <- gamma_prime - mean_gamma

  chol_S <- chol(0.5 * (cov + t(cov)))

  log_det <- 2 * sum(log(diag(chol_S)))

  u <- backsolve(chol_S, d, transpose = TRUE)

  mah <- sum(u * u)

  -(q / 2) * log(2 * pi) - 0.5 * log_det - 0.5 * mah

}



#' @noRd

.minorization_cov_m_beta <- function(p11, V_list) {

  P11_inv <- chol2inv(p11$chol_P11)

  cov_m <- matrix(0, p11$q, p11$q)

  for (j in seq_len(p11$J)) {

    lev <- p11$group_levels[[j]]

    H_j <- p11$H_list[[lev]]

    V_j <- V_list[[lev]]

    cov_m <- cov_m +

      P11_inv %*% t(H_j) %*% p11$P_b %*% V_j %*% p11$P_b %*% H_j %*% P11_inv

  }

  0.5 * (cov_m + t(cov_m))

}



#' @noRd

.minorization_g_gaussian <- function(gamma_prime, gamma_star, p11, V_list) {

  cov_kernel <- p11$Sigma_star + .minorization_cov_m_beta(p11, V_list)

  log_q <- .minorization_log_gaussian(gamma_prime, gamma_star, cov_kernel)

  log_q_Q <- .minorization_log_q_Q(gamma_prime, gamma_star, p11)

  log_q - log_q_Q

}



#' @noRd

.minorization_validate_mode_for_optimize <- function(mode, fn_name = "epsilon_optimize") {

  if (!is.list(mode) || is.null(mode$fixef)) {

    stop("'", fn_name, "()' requires 'mode' from population_mode().",

         call. = FALSE)

  }

  req <- c(

    "gamma_star", "p11", "design", "family", "measurement_prior_list"

  )

  missing <- req[!vapply(req, function(nm) !is.null(mode[[nm]]), logical(1L))]

  if (length(missing)) {

    stop(

      "'", fn_name, "()' requires a fresh population_mode() result with ",

      paste(missing, collapse = ", "), ".",

      call. = FALSE

    )

  }

  if (identical(mode$family$family, "gaussian") && is.null(mode$V_list)) {

    stop("'", fn_name, "()' requires 'mode$V_list' for the Gaussian kernel.",

         call. = FALSE)

  }

  invisible(mode)

}



#' Numerical epsilon(gamma*) from the minorization profile.

#'

#' Solves the convex program of Chapter C05 Definition 5 at the EM fixed

#' point: minimize \eqn{g(\gamma' \mid \gamma^\star) = \log q(\gamma' \mid

#' \gamma^\star) - \log q_Q(\gamma')} over destination \eqn{\gamma'}, with

#' the kernel integral over \eqn{\beta} uses an exact Gaussian convolution when

#' \code{family = gaussian()}, otherwise Monte Carlo draws from

#' \eqn{\pi(\beta \mid \gamma^\star, y)}.

#'

#' @param mode Object returned by \code{\link{population_mode}}.

#' @param n Number of \eqn{\beta} draws for the MC inner integral when
#'   the likelihood is not Gaussian. Defaults to \code{mode$n} from
#'   \code{\link{population_mode}}, or \code{10000L}.

#' @param mc_seed Optional seed for the cached \eqn{\beta} draws (fixed for the

#'   whole optimization).

#' @param gamma_prime_start Initial destination \eqn{\gamma'}; defaults to

#'   \code{mode$gamma_star}.

#' @param control A list passed to \code{\link[stats]{optim}} (\code{method =

#'   "BFGS"}).

#' @param ... Ignored.

#' @return A list with \code{eps_star}, \code{method}, \code{certified},

#'   \code{g_opt}, \code{gamma_prime}, \code{optim}, and \code{n}.

#' @seealso \code{\link{population_mode}}, \code{\link{epsilon_star}}

#' @export

epsilon_optimize <- function(mode,

                             n = NULL,

                             mc_seed = NULL,

                             gamma_prime_start = NULL,

                             control = list(maxit = 200L),

                             ...) {

  mode <- .minorization_validate_mode_for_optimize(mode)

  if (is.null(n)) {
    n <- mode$n
    if (is.null(n) || is.na(n)) {
      n <- 10000L
    }
  }

  n <- as.integer(n)

  if (!(n >= 1L)) {

    stop("'n' must be at least 1.", call. = FALSE)

  }



  p11 <- mode$p11

  gamma_star <- mode$gamma_star

  if (is.null(gamma_prime_start)) {

    gamma_prime_start <- gamma_star

  } else {

    gamma_prime_start <- as.numeric(gamma_prime_start)

    if (length(gamma_prime_start) != p11$q) {

      stop(

        "'gamma_prime_start' must have length ", p11$q, ".",

        call. = FALSE

      )

    }

  }



  is_gaussian <- identical(mode$family$family, "gaussian")
  use_mc <- identical(mode$estep, "mc")

  if (is_gaussian && !use_mc) {

    obj <- function(gamma_prime) {

      g <- .minorization_g_gaussian(

        gamma_prime = gamma_prime,

        gamma_star = gamma_star,

        p11 = p11,

        V_list = mode$V_list

      )

      if (!is.finite(g)) {

        1e100

      } else {

        g

      }

    }

    n_used <- NA_integer_

  } else {

    beta_draws <- .minorization_draw_beta(

      design = mode$design,

      fixef = mode$fixef,

      p11 = p11,

      family = mode$family,

      measurement_prior_list = mode$measurement_prior_list,

      n = n,

      mc_seed = mc_seed

    )



    obj <- function(gamma_prime) {

      g <- .minorization_g_mc(

        gamma_prime = gamma_prime,

        gamma_star = gamma_star,

        beta_draws = beta_draws,

        p11 = p11

      )

      if (!is.finite(g)) {

        1e100

      } else {

        g

      }

    }

    n_used <- n

  }



  opt <- stats::optim(

    par = gamma_prime_start,

    fn = obj,

    method = "BFGS",

    control = control

  )



  g_opt <- opt$value

  if (!is.finite(g_opt)) {

    stop("Optimization of g(gamma' | gamma*) failed to converge.", call. = FALSE)

  }



  certified <- FALSE

  if (!is.null(mode$eps_star_closure)) {

    closure_gap <- abs(exp(g_opt) - mode$eps_star_closure)

    if (is_gaussian && identical(mode$estep, "exact")) {

      certified <- closure_gap < 1e-6 * max(1, mode$eps_star_closure)

    } else if (use_mc && !is.null(mode$mc_delta_floor)) {

      certified <- closure_gap <= max(mode$mc_delta_floor, 1e-6 * max(1, mode$eps_star_closure))

    }

  }



  list(

    eps_star = exp(g_opt),

    method = "optimize",

    certified = certified,

    g_opt = g_opt,

    gamma_prime = opt$par,

    optim = opt,

    n = n_used

  )

}


