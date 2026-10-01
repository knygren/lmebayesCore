## Chapter C05 restricted Gibbs: EM fixed point gamma* and Jacobian.

#' Restricted E-step: aggregate joint draws from safe Block~1 simfuncs.
#' @noRd
.c05_estep_restricted <- function(design,
                                  fixef,
                                  p11,
                                  measurement_prior_list,
                                  family,
                                  beta_set,
                                  n = 10000L,
                                  mc_seed = NULL,
                                  max_tries = NULL) {
  safe_out <- .c05_safe_group_draw(
    design = design,
    fixef = fixef,
    p11 = p11,
    measurement_prior_list = measurement_prior_list,
    family = family,
    beta_set = beta_set,
    n = n,
    max_tries = max_tries,
    mc_seed = mc_seed
  )
  out <- .c05_safe_group_draw_to_estep(safe_out, p11$group_levels)
  out$safe_out <- safe_out
  out
}

#' Conditional-mean E-step for the population mean map.
#' @param beta_set Optional \code{beta_marginal_safe_set}; when supplied, uses
#'   \code{\link{rNormal_reg_group_safe}} or
#'   \code{\link{rNormalGLM_reg_group_safe}} by \code{family}.
#' @noRd
.c05_estep <- function(design,
                                    fixef,
                                    p11,
                                    measurement_prior_list,
                                    family,
                                    estep = c("exact", "aghq", "mc"),
                                    n = 10000L,
                                    mc_seed = NULL,
                                    beta_set = NULL) {
  estep <- match.arg(estep)
  if (identical(estep, "aghq")) {
    stop("estep = \"aghq\" is not implemented yet; use \"mc\".", call. = FALSE)
  }

  if (!is.null(beta_set)) {
    return(.c05_estep_restricted(
      design = design,
      fixef = fixef,
      p11 = p11,
      measurement_prior_list = measurement_prior_list,
      family = family,
      beta_set = beta_set,
      n = n,
      mc_seed = mc_seed
    ))
  }

  if (!identical(family$family, "gaussian") && identical(estep, "exact")) {
    stop(
      "estep = \"exact\" is only available for gaussian(); ",
      "use \"mc\" or \"aghq\".",
      call. = FALSE
    )
  }

  group_levels <- p11$group_levels
  J <- p11$J
  p_re <- p11$p_re
  re_names <- p11$re_names
  g_chr <- as.character(design$group)
  gamma <- .c05_gamma_from_fixef(fixef, p11)
  sigma2 <- measurement_prior_list$group.dispersion
  P_b <- p11$P_b
  Sigma_b <- measurement_prior_list$group.Sigma
  is_gaussian <- identical(family$family, "gaussian")

  b_mean <- matrix(
    0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )
  b_mc_se <- NULL
  V_list <- stats::setNames(vector("list", J), group_levels)

  if (is_gaussian && identical(estep, "exact")) {
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

    if (!is.null(p11$lmerb_system)) {
      b_out <- .lmerb_posterior_b_given_gamma(p11$lmerb_system, design, fixef)
      b_mean <- b_out$b
      for (lev in group_levels) {
        V_list[[lev]] <- solve(p11$lmerb_system$post_P_j_list[[lev]])
      }
    } else {
      for (jj in seq_len(J)) {
        lev <- group_levels[jj]
        rows <- which(g_chr == lev)
        Z_j <- design$D[rows, , drop = FALSE]
        y_j <- design$y[rows]
        sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
        H_j <- p11$H_list[[lev]]
        mu_j <- as.vector(H_j %*% gamma)

        post_P_j <- crossprod(Z_j) / sigma2_j + P_b
        post_v_j <- crossprod(Z_j, y_j) / sigma2_j + P_b %*% mu_j
        b_mean[jj, ] <- solve(post_P_j, post_v_j)
        V_list[[lev]] <- solve(post_P_j)
      }
    }
  } else if (identical(estep, "mc")) {
    if (!is.null(mc_seed)) set.seed(mc_seed)
    mu_all <- as.matrix(build_mu_all(design, fixef, group_levels = group_levels)$mu_all)
    if (is_gaussian) {
      if (is.null(sigma2)) {
        stop(
          "'measurement_prior_list$group.dispersion' is required for gaussian().",
          call. = FALSE
        )
      }
      sigma2 <- as.numeric(sigma2)
    }

    for (jj in seq_len(J)) {
      lev <- group_levels[jj]
      rows <- which(g_chr == lev)
      y_j <- design$y[rows]
      Z_j <- design$D[rows, , drop = FALSE]
      mu_j <- mu_all[, jj]

      pf_j <- if (is_gaussian) {
        sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
        glmbayesCore::dNormal(
          mu = mu_j, Sigma = Sigma_b, dispersion = sigma2_j
        )
      } else {
        glmbayesCore::dNormal(mu = mu_j, Sigma = Sigma_b)
      }
      draws <- matrix(0, n, p_re)
      for (m in seq_len(n)) {
        fit_j <- glmbayesCore::rglmb(
          n = 1L,
          y = y_j,
          x = Z_j,
          family = family,
          pfamily = pf_j,
          verbose = FALSE
        )
        draws[m, ] <- .minorization_rglmb_draw(fit_j)
      }
      b_mean[jj, ] <- colMeans(draws)
      V_list[[lev]] <- stats::cov(draws)
      if (is.null(b_mc_se)) {
        b_mc_se <- matrix(
          0, nrow = J, ncol = p_re,
          dimnames = list(group_levels, re_names)
        )
      }
      b_mc_se[jj, ] <- apply(draws, 2L, stats::sd) / sqrt(n)
    }
  }

  list(
    b_mean = b_mean,
    b_mc_se = b_mc_se,
    V_list = V_list,
    estep = estep,
    n = if (identical(estep, "mc")) n else NA_integer_,
    n_target = NA_integer_,
    n_tried = NA_integer_,
    accept_rate = NA_real_,
    beta_set = beta_set,
    restricted = FALSE
  )
}

#' One C05 mean-map update M(gamma) given conditional means b_j.
#' @noRd
.c05_mean_map <- function(b_mean, p11) {
  if (!is.null(p11$lmerb_system) && is.null(b_mean)) {
    return(.c05_mean_map_lmerb(b_mean, p11))
  }

  rhs <- p11$Lambda_gamma %*% p11$mu_0
  for (j in seq_len(p11$J)) {
    b_j <- b_mean[j, ]
    H_j <- p11$H_list[[p11$group_levels[[j]]]]
    rhs <- rhs + t(H_j) %*% p11$P_b %*% b_j
  }
  gamma <- solve(p11$P11, rhs)
  .c05_fixef_from_gamma(gamma, p11)
}

#' Block~2 update using the exact Gaussian Schur system (\code{lmerb}).
#' @noRd
.c05_mean_map_lmerb <- function(b_mean, p11) {
  system <- p11$lmerb_system
  gamma <- solve(system$M, system$v)
  .c05_fixef_from_gamma(gamma, p11)
}

#' Conditional-mode Block~1 update for ICM initialization.
#' @noRd
.c05_block1_modes <- function(design,
                              fixef,
                              p11,
                              measurement_prior_list,
                              family) {
  group_levels <- p11$group_levels
  J <- p11$J
  p_re <- p11$p_re
  re_names <- p11$re_names
  g_chr <- as.character(design$group)
  gamma <- .c05_gamma_from_fixef(fixef, p11)
  sigma2 <- measurement_prior_list$group.dispersion
  P_b <- p11$P_b
  Sigma_b <- measurement_prior_list$group.Sigma
  is_gaussian <- identical(family$family, "gaussian")

  b_mode <- matrix(
    0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )

  if (is_gaussian) {
    sigma2 <- as.numeric(sigma2)
    for (jj in seq_len(J)) {
      lev <- group_levels[jj]
      rows <- which(g_chr == lev)
      Z_j <- design$D[rows, , drop = FALSE]
      y_j <- design$y[rows]
      sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
      H_j <- p11$H_list[[lev]]
      mu_j <- as.vector(H_j %*% gamma)
      post_P_j <- crossprod(Z_j) / sigma2_j + P_b
      post_v_j <- crossprod(Z_j, y_j) / sigma2_j + P_b %*% mu_j
      b_mode[jj, ] <- solve(post_P_j, post_v_j)
    }
  } else {
    mu_all <- as.matrix(
      build_mu_all(design, fixef, group_levels = group_levels)$mu_all
    )
    for (jj in seq_len(J)) {
      lev <- group_levels[jj]
      rows <- which(g_chr == lev)
      y_j <- design$y[rows]
      Z_j <- design$D[rows, , drop = FALSE]
      mu_j <- mu_all[, jj]
      pf_j <- glmbayesCore::dNormal(mu = mu_j, Sigma = Sigma_b)
      fit_j <- glmbayesCore::rglmb(
        n = 1L,
        y = y_j,
        x = Z_j,
        family = family,
        pfamily = pf_j,
        verbose = FALSE
      )
      b_mode[jj, ] <- .minorization_rglmb_mode(fit_j)
    }
  }

  b_mode
}

#' ICM start for \code{population_mode()} (conditional modes + C05 mean map).
#' @noRd
.c05_icm_init <- function(design,
                          fixef_start,
                          p11,
                          measurement_prior_list,
                          family,
                          tol = 1e-8,
                          maxit = 200L) {
  if (identical(family$family, "gaussian") &&
      !is.null(p11$lmerb_system)) {
    pm <- lmerb_posterior_mean(design, measurement_prior_list, tol, maxit)
    return(list(
      fixef = pm$fixef,
      b_mode = pm$b_mean,
      converged = pm$converged,
      iterations = pm$iterations,
      delta = pm$delta
    ))
  }

  fixef <- fixef_start
  converged <- FALSE
  delta <- NA_real_

  for (iter in seq_len(maxit)) {
    b_mode <- .c05_block1_modes(
      design = design,
      fixef = fixef,
      p11 = p11,
      measurement_prior_list = measurement_prior_list,
      family = family
    )
    fixef_new <- .c05_mean_map(b_mode, p11)
    delta <- .c05_gamma_delta(fixef, fixef_new, p11)
    fixef <- fixef_new
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "ICM initialization did not converge in ", maxit,
      " iterations (final delta = ", signif(delta, 3L), ").",
      call. = FALSE
    )
  }

  list(
    fixef = fixef,
    b_mode = b_mode,
    converged = converged,
    iterations = iter,
    delta = delta
  )
}

#' Mahalanobis change in \eqn{\gamma} under \eqn{P_{11}}.
#' @noRd
.c05_gamma_delta <- function(fixef_old, fixef_new, p11) {
  d_gamma <- .c05_gamma_from_fixef(fixef_new, p11) -
    .c05_gamma_from_fixef(fixef_old, p11)
  sqrt(as.numeric(crossprod(d_gamma, p11$P11 %*% d_gamma)))
}

#' Stable log-sum-exp.
#' @noRd
.c05_logsumexp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(-Inf)
  }
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Column names for the stacked \eqn{\gamma} vector.
#' @noRd
.c05_gamma_names <- function(p11) {
  unlist(
    lapply(p11$re_names, function(k) {
      cols <- colnames(p11$X_hyper[[k]])
      if (is.null(cols)) {
        cols <- seq_len(ncol(p11$X_hyper[[k]]))
      }
      paste0(k, "_", cols)
    }),
    use.names = FALSE
  )
}

#' Log prior \eqn{\log \pi(\gamma)} under the population prior.
#' @noRd
.c05_log_logprior_gamma <- function(gamma, p11) {
  diff0 <- gamma - p11$mu_0
  q <- length(gamma)
  chol_L <- chol(p11$Lambda_gamma)
  logdet <- 2 * sum(log(diag(chol_L)))
  quad <- sum((backsolve(chol_L, diff0, transpose = TRUE))^2)
  -0.5 * (q * log(2 * pi) + logdet + quad)
}

#' Group log-likelihood \eqn{\ell_j(b_j)} at fixed \eqn{b_j}.
#' @noRd
.c05_group_loglik <- function(b_j, rows, design, family, measurement_prior_list, jj) {
  Z_j <- design$D[rows, , drop = FALSE]
  eta <- as.vector(Z_j %*% b_j)
  if (!is.null(design$offset) && length(design$offset) >= max(rows)) {
    eta <- eta + design$offset[rows]
  }
  y_j <- design$y[rows]
  wt <- if (!is.null(design$weights)) design$weights[rows] else rep(1, length(rows))

  if (identical(family$family, "gaussian")) {
    sigma2 <- measurement_prior_list$group.dispersion
    sigma2 <- as.numeric(sigma2)
    sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
    sd_j <- sqrt(sigma2_j / wt)
    sum(stats::dnorm(y_j, mean = eta, sd = sd_j, log = TRUE) + log(wt))
  } else if (identical(family$family, "binomial")) {
    p_j <- 1 / (1 + exp(-eta))
    if (all(wt == 1)) {
      sum(stats::dbinom(y_j, size = 1, prob = p_j, log = TRUE))
    } else {
      sum(stats::dbinom(y_j, size = wt, prob = p_j, log = TRUE))
    }
  } else {
    stop(
      "EM log marginal is not implemented for family ",
      family$family, "().",
      call. = FALSE
    )
  }
}

#' Log \eqn{\pi(b \mid \gamma) + \ell(b)} summed over groups at a joint draw.
#' @noRd
.c05_log_cond_beta_at_bmat <- function(b_mat,
                                       gamma,
                                       design,
                                       p11,
                                       family,
                                       measurement_prior_list) {
  group_levels <- p11$group_levels
  g_chr <- as.character(design$group)
  P_b <- p11$P_b
  chol_Pb <- chol(P_b)
  logdet_Pb <- 2 * sum(log(diag(chol_Pb)))
  p_re <- p11$p_re
  ell <- 0
  log_prior <- 0
  for (jj in seq_along(group_levels)) {
    lev <- group_levels[[jj]]
    rows <- which(g_chr == lev)
    b_j <- as.numeric(b_mat[jj, , drop = FALSE])
    H_j <- p11$H_list[[lev]]
    mu_j <- as.vector(H_j %*% gamma)
    diff_b <- b_j - mu_j
    log_prior <- log_prior - 0.5 * (
      p_re * log(2 * pi) + logdet_Pb +
        sum((backsolve(chol_Pb, diff_b, transpose = TRUE))^2)
    )
    ell <- ell + .c05_group_loglik(
      b_j, rows, design, family, measurement_prior_list, jj
    )
  }
  log_prior + ell
}

#' MC estimate of \eqn{\log \pi(\gamma \mid y)} from joint safe draws.
#' @noRd
.c05_log_marginal_logpost <- function(fixef,
                                      estep_out,
                                      design,
                                      p11,
                                      family,
                                      measurement_prior_list) {
  safe_out <- estep_out$safe_out
  if (is.null(safe_out)) {
    return(list(log_post = NA_real_, method = NA_character_))
  }
  acc <- safe_out$coefficients_all
  if (is.null(acc)) {
    b1 <- matrix(safe_out$coefficients, nrow = 1L)
    acc <- array(
      b1,
      dim = c(1L, nrow(b1), ncol(b1)),
      dimnames = list(NULL, rownames(b1), colnames(b1))
    )
  }
  gamma <- .c05_gamma_from_fixef(fixef, p11)
  R <- dim(acc)[1L]
  log_cond <- vapply(seq_len(R), function(r) {
    b_mat <- acc[r, , , drop = FALSE]
    dim(b_mat) <- c(dim(acc)[2L], dim(acc)[3L])
    .c05_log_cond_beta_at_bmat(
      b_mat, gamma, design, p11, family, measurement_prior_list
    )
  }, numeric(1L))
  log_p_gamma <- .c05_log_logprior_gamma(gamma, p11)
  list(
    log_post = log_p_gamma + .c05_logsumexp(log_cond) - log(R),
    method = "mc_marginal"
  )
}

#' EM gradient \eqn{g = P_{11}(\gamma - M(\gamma))} at the current \code{fixef}.
#' @noRd
.c05_em_gradient <- function(fixef, estep_out, p11) {
  gamma <- .c05_gamma_from_fixef(fixef, p11)
  fixef_M <- .c05_mean_map(estep_out$b_mean, p11)
  gamma_M <- .c05_gamma_from_fixef(fixef_M, p11)
  g <- as.numeric(p11$P11 %*% (gamma - gamma_M))
  gnames <- .c05_gamma_names(p11)
  names(g) <- gnames
  names(gamma) <- gnames
  names(gamma_M) <- gnames
  diff_g <- gamma - gamma_M
  list(
    g = g,
    gamma = gamma,
    gamma_M = gamma_M,
    fixef_M = fixef_M,
    delta = sqrt(as.numeric(crossprod(diff_g, p11$P11 %*% diff_g)))
  )
}

#' One EM diagnostic record (gradient, log marginal, tolerances).
#' @noRd
.c05_em_trace_record <- function(iter,
                                 fixef,
                                 estep_out,
                                 design,
                                 p11,
                                 family,
                                 measurement_prior_list,
                                 tol,
                                 mc_alpha,
                                 uses_mc,
                                 prev_delta = NA_real_,
                                 prev_log_post = NA_real_) {
  gr <- .c05_em_gradient(fixef, estep_out, p11)
  lp <- .c05_log_marginal_logpost(
    fixef, estep_out, design, p11, family, measurement_prior_list
  )
  gnames <- .c05_gamma_names(p11)
  gamma <- gr$gamma
  names(gamma) <- gnames
  g <- gr$g
  names(g) <- gnames
  tol_eff <- if (uses_mc) {
    .c05_em_tol(tol, "mc", estep_out$b_mc_se, p11, mc_alpha = mc_alpha)
  } else {
    tol
  }
  d_log_post <- if (is.finite(prev_log_post) && is.finite(lp$log_post)) {
    lp$log_post - prev_log_post
  } else {
    NA_real_
  }
  slope_g <- if (is.finite(prev_delta) && prev_delta > 0) {
    gr$delta / prev_delta
  } else {
    NA_real_
  }
  list(
    iter = iter,
    gamma = gr$gamma,
    g = gr$g,
    delta = gr$delta,
    log_post = lp$log_post,
    log_post_method = lp$method,
    d_log_post = d_log_post,
    slope_g = slope_g,
    tol_eff = tol_eff,
    pass = isTRUE(gr$delta <= tol_eff)
  )
}

#' Format a named numeric vector for EM trace messages.
#' @noRd
.c05_format_named_vec <- function(x, digits = 4L) {
  nm <- names(x)
  x <- as.vector(x, mode = "double")
  if (is.null(nm)) {
    nm <- seq_along(x)
  }
  paste0(
    nm, "=",
    format(signif(x, digits), trim = TRUE, scientific = FALSE),
    collapse = ", "
  )
}

#' Print EM trace records when \code{verbose = TRUE}.
#' @noRd
.c05_verbose_em_trace <- function(record, p11, label = NULL) {
  gnames <- .c05_gamma_names(p11)
  gamma_fmt <- stats::setNames(as.numeric(record$gamma), gnames)
  g_fmt <- stats::setNames(as.numeric(record$g), gnames)
  if (!is.null(label)) {
    message(label, ": ", .c05_format_named_vec(gamma_fmt), sep = "")
  }
  lp <- if (is.finite(record$log_post)) {
    signif(record$log_post, 6)
  } else {
    "NA"
  }
  msg <- paste0(
    "  iter ", sprintf("%2d", record$iter),
    " | ||g||=", format(record$delta, digits = 4, scientific = TRUE),
    " | slope_g=",
    if (is.finite(record$slope_g)) {
      format(record$slope_g, digits = 4)
    } else {
      "NA"
    },
    " | log_post=", lp,
    " | d_log=",
    if (is.finite(record$d_log_post)) {
      paste0(
        if (record$d_log_post >= 0) "+" else "",
        format(record$d_log_post, digits = 5)
      )
    } else {
      "NA"
    },
    " | tol_eff=", format(record$tol_eff, digits = 4, scientific = TRUE),
    " | pass=", record$pass
  )
  message(msg)
  message("      g: ", .c05_format_named_vec(g_fmt))
}

#' Print a \code{population_mode_em_trace} object.
#' @param x A \code{population_mode_em_trace} object (from
#'   \code{\link{population_mode}$em_trace}).
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.population_mode_em_trace <- function(x, ...) {
  if (!inherits(x, "population_mode_em_trace")) {
    NextMethod("print")
  }
  if (!is.null(x$icm)) {
    icm <- x$icm
    message(
      "ICM init: iters=", icm$iterations,
      " converged=", icm$converged,
      " delta=", signif(icm$delta, 4)
    )
  }
  if (!is.null(x$em_start)) {
    .c05_verbose_em_trace(
      x$em_start, x$p11,
      label = "  gamma_start"
    )
  }
  if (length(x$history)) {
    message("EM history (", length(x$history), " records):")
    for (rec in x$history) {
      .c05_verbose_em_trace(rec, x$p11)
    }
  }
  invisible(x)
}

#' MC noise floor for \eqn{\|\Delta\gamma\|_{P_{11}}} from \code{b_mc_se}.
#'
#' Propagates independent elementwise MC standard errors in \code{b_mean}
#' through the C05 mean map and returns one MC standard deviation of the
#' resulting \eqn{P_{11}}-norm update (i.e. stop when \code{delta} is not
#' statistically distinguishable from MC noise at the current \code{n}).
#' @noRd
.c05_mc_delta_floor <- function(b_mc_se, p11) {
  if (is.null(b_mc_se)) {
    return(0)
  }
  se2 <- as.numeric(b_mc_se)
  if (!any(is.finite(se2) & se2 > 0)) {
    return(0)
  }

  P11_inv <- chol2inv(p11$chol_P11)
  Cov_gamma <- matrix(0, p11$q, p11$q)
  for (j in seq_len(p11$J)) {
    lev <- p11$group_levels[[j]]
    H_j <- p11$H_list[[lev]]
    A_j <- P11_inv %*% t(H_j) %*% p11$P_b
    se2_j <- b_mc_se[j, ]^2
    Cov_gamma <- Cov_gamma + A_j %*% diag(se2_j, p11$p_re) %*% t(A_j)
  }

  sq <- sum(diag(p11$P11 %*% Cov_gamma))
  sqrt(max(sq, 0))
}

#' Whether the E-step uses Monte Carlo (including widetilde-B rejection).
#' @noRd
.c05_estep_uses_mc <- function(family, estep, beta_set) {
  if (!is.null(beta_set)) {
    return(TRUE)
  }
  if (!identical(family$family, "gaussian") && identical(estep, "exact")) {
    return(TRUE)
  }
  identical(estep, "mc")
}

#' Effective EM tolerance when the E-step is Monte Carlo.
#' @noRd
.c05_em_tol <- function(tol, estep, b_mc_se, p11, mc_alpha = 0.05) {
  if (!identical(estep, "mc") && is.null(b_mc_se)) {
    return(tol)
  }
  floor <- .c05_mc_delta_floor(b_mc_se, p11)
  if (!is.finite(floor) || floor <= 0) {
    return(tol)
  }
  z <- stats::qnorm(1 - mc_alpha / 2)
  max(tol, z * floor)
}

#' MC fixed-point screen at the ICM start (Type I + Type II combined).
#' @noRd
.c05_icm_screen <- function(design,
                            fixef,
                            p11,
                            measurement_prior_list,
                            family,
                            estep,
                            n,
                            mc_seed,
                            beta_set,
                            tol,
                            mc_alpha,
                            mc_stable) {
  mc_stable <- as.integer(mc_stable)
  if (!(mc_stable >= 1L)) {
    stop("'mc_stable' must be at least 1.", call. = FALSE)
  }

  delta <- NA_real_
  mc_delta_floor <- NA_real_
  tol_eff <- tol
  estep_out <- NULL

  for (hit in seq_len(mc_stable)) {
    estep_out <- .c05_estep(
      design = design,
      fixef = fixef,
      p11 = p11,
      measurement_prior_list = measurement_prior_list,
      family = family,
      estep = estep,
      n = n,
      mc_seed = if (hit == 1L) mc_seed else NULL,
      beta_set = beta_set
    )
    fixef_new <- .c05_mean_map(estep_out$b_mean, p11)
    delta <- .c05_gamma_delta(fixef, fixef_new, p11)
    mc_delta_floor <- .c05_mc_delta_floor(estep_out$b_mc_se, p11)
    tol_eff <- .c05_em_tol(tol, "mc", estep_out$b_mc_se, p11, mc_alpha = mc_alpha)
    if (!(delta <= tol_eff)) {
      return(list(
        passed = FALSE,
        hits = hit - 1L,
        estep_out = estep_out,
        delta = delta,
        mc_delta_floor = mc_delta_floor,
        tol_eff = tol_eff
      ))
    }
  }

  list(
    passed = TRUE,
    hits = mc_stable,
    estep_out = estep_out,
    delta = delta,
    mc_delta_floor = mc_delta_floor,
    tol_eff = tol_eff
  )
}

#' Jacobian J, spectrum, and closure objects at the current E-step.
#' @noRd
.c05_jacobian <- function(p11, estep_out) {
  q <- p11$q
  P11_inv <- chol2inv(p11$chol_P11)
  J_mat <- matrix(0, q, q)

  for (j in seq_len(p11$J)) {
    lev <- p11$group_levels[[j]]
    H_j <- p11$H_list[[lev]]
    V_j <- estep_out$V_list[[lev]]
    J_mat <- J_mat + P11_inv %*% t(H_j) %*% p11$P_b %*% V_j %*% p11$P_b %*% H_j
  }

  spec <- .c05_coupling_spectrum(p11, J_mat, strict = FALSE)
  kappa <- spec$kappa_max
  rho <- spec$rho

  Iq <- diag(q)
  Sigma_pi <- tryCatch(
    solve(p11$P11 %*% (Iq - J_mat)),
    error = function(e) NULL
  )

  list(
    tilde_J = J_mat,
    kappa_spectrum = spec$kappa,
    weights = spec$weights,
    kappa = kappa,
    rho = rho,
    Sigma_pi = Sigma_pi
  )
}

#' EM fixed point for the C05 population mode gamma*.
#'
#' Finds the fixed point of the conditional-mean map (Chapter C05 Stage 1).
#' Uses conditional means for Block~1, not modes.
#'
#' @details
#' In the restricted certificate, the EM E-step uses conditional means under
#' the \eqn{\beta}-truncated target on \eqn{\widetilde B(\delta_2)}. Pass
#' Pass \code{beta_set} from \code{\link{beta_marginal_safe_set}}; routes to
#' \code{\link{rNormal_reg_group_safe}} (Gaussian) or
#' \code{\link{rNormalGLM_reg_group_safe}} (non-Gaussian).
#'
#' @param design A \code{\link{model_setup}} list.
#' @param pfamily_list Block~2 prior list from \code{\link{pfamily_list}()}.
#' @param family A \code{\link[stats]{family}} object.
#' @param dispprior_list Optional Block~1 dispersion prior for \code{gaussian()}.
#' @param beta_set Optional \code{\link{beta_marginal_safe_set}}; enables
#'   widetilde-B rejection MC (see Details).
#' @param estep E-step tier: \code{"exact"} (default Gaussian closed form when
#'   unrestricted), \code{"mc"} (simulated conditional means; used automatically
#'   when \code{beta_set} is supplied), or \code{"aghq"} (not yet implemented).
#' @param acceleration \code{"none"} or \code{"squarem"} (not yet implemented).
#' @param n Number of \eqn{\beta_j} draws per group when \code{estep = "mc"}.
#'   The sample mean \eqn{\bar b_j} has elementwise MC standard error
#'   \eqn{\mathrm{sd}(\beta_j \mid \gamma, y)/\sqrt{n}}; relative error on
#'   means is therefore \eqn{O(1/\sqrt{n})} (e.g. \code{n = 10000} targets
#'   \eqn{\approx 1\%} of the posterior scale).
#' @param mc_seed Optional seed for the MC E-step (first iteration only).
#' @param icm_init If \code{TRUE}, run iterated conditional modes (Block~1
#'   modes via \code{\link[glmbayesCore]{rglmb}}, Block~2 via the C05 mean map)
#'   before EM. Skipped on the Gaussian exact closed-form path. Gaussian
#'   Block~1 uses exact conditional modes.
#' @param icm_tol ICM convergence tolerance on \eqn{\|\Delta\gamma\|_{P_{11}}}.
#' @param icm_maxit Maximum ICM iterations.
#' @param icm_screen If \code{TRUE} (default when the E-step is Monte Carlo),
#'   after ICM run \code{mc_stable} MC mean-map fixed-point checks at the ICM
#'   \code{fixef}; skip the EM loop when all pass. If \code{NULL}, enabled
#'   whenever \code{icm_init} is \code{TRUE} and the E-step uses MC.
#' @param mc_alpha Two-sided level for the MC fixed-point screen and EM stop
#'   (\code{tol_eff = max(tol, z_{1-alpha/2} * mc_delta_floor)}).
#' @param mc_stable Number of consecutive MC fixed-point passes required at the
#'   same \code{fixef} (ICM screen and EM stopping rule).
#' @param tol Convergence tolerance on the Mahalanobis change in \code{fixef}
#'   under the C05 metric \eqn{\|\Delta\gamma\|_{P_{11}}}. When the E-step
#'   is Monte Carlo, the effective tolerance is
#'   \code{max(tol, z_{1-alpha/2} * mc_delta_floor)}.
#' @param maxit Maximum EM iterations.
#' @param verbose If \code{TRUE}, emit \code{message()} diagnostics for ICM
#'   initialization and each EM record (gradient \eqn{g}, \eqn{\|g\|}, MC
#'   \eqn{\log \pi(\gamma \mid y)} when joint safe draws are available).
#' @return A list with \code{fixef} (\eqn{\gamma^\star}), \code{gamma_star},
#'   \code{b_mean}, \code{b_mc_se} (terminal MC standard errors when
#'   \code{estep = "mc"}), \code{V_list}, \code{n}, model context
#'   (\code{design}, \code{family}, \code{measurement_prior_list}, \code{p11}),
#'   refresh objects (\code{P11}, \code{Sigma_star}), \code{tilde_J},
#'   \code{kappa_spectrum}, \code{weights}, \code{kappa}, \code{rho},
#'   \code{Sigma_pi}, \code{eps_star_closure},
#'   \code{icm} and \code{icm_screen} diagnostics, EM diagnostics
#'   (\code{tol_eff}, \code{mc_delta_floor}, \code{em_route}, \code{em_trace}),
#'   \code{beta_set},
#'   \code{restricted}, and related fields.
#' @seealso \code{\link{group_effects_conditional_mean}},
#'   \code{\link{rNormal_reg_group_safe}},
#'   \code{\link{rNormalGLM_reg_group_safe}}
#' @export
population_mode <- function(design,
                                             pfamily_list,
                                             family = gaussian(),
                                             dispprior_list = NULL,
                                             beta_set = NULL,
                                             estep = c("exact", "aghq", "mc"),
                                             acceleration = c("none", "squarem"),
                                             n = 10000L,
                                             mc_seed = NULL,
                                             icm_init = TRUE,
                                             icm_tol = 1e-8,
                                             icm_maxit = 200L,
                                             icm_screen = NULL,
                                             mc_alpha = 0.05,
                                             mc_stable = 2L,
                                             tol = 1e-10,
                                             maxit = 200L,
                                             verbose = FALSE) {
  estep <- match.arg(estep)
  acceleration <- match.arg(acceleration)
  if (!is.null(beta_set) && identical(estep, "exact")) {
    estep <- "mc"
  }
  if (!identical(acceleration, "none")) {
    stop("acceleration = \"squarem\" is not implemented yet.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!(n >= 1L)) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  if (!is.null(beta_set) && !inherits(beta_set, "beta_marginal_safe_set")) {
    stop("'beta_set' must be from beta_marginal_safe_set().", call. = FALSE)
  }

  prep <- .c05_validate(
    design, pfamily_list, family, dispprior_list
  )
  mpl <- prep$measurement_prior_list
  p11 <- .c05_p11(
    design,
    prep$prior_pack,
    prep$group_levels,
    measurement_prior_list = mpl,
    family = family
  )

  fixef <- lapply(mpl$pop.prior_list, `[[`, "mu")
  names(fixef) <- prep$re_names

  uses_mc <- .c05_estep_uses_mc(family, estep, beta_set)
  is_gauss_exact <- identical(family$family, "gaussian") &&
    identical(estep, "exact") &&
    is.null(beta_set) &&
    !is.null(p11$lmerb_system)

  if (is.null(icm_screen)) {
    icm_screen <- isTRUE(icm_init) && uses_mc
  }

  icm <- NULL
  icm_screen_out <- NULL
  em_trace <- NULL
  mc_seed_used <- FALSE
  estep_at_start <- NULL
  if (isTRUE(icm_init) && !is_gauss_exact) {
    icm <- .c05_icm_init(
      design = design,
      fixef_start = fixef,
      p11 = p11,
      measurement_prior_list = mpl,
      family = family,
      tol = icm_tol,
      maxit = icm_maxit
    )
    fixef <- icm$fixef
    if (isTRUE(verbose)) {
      message(
        "ICM init: iters=", icm$iterations,
        " converged=", icm$converged,
        " delta=", signif(icm$delta, 4)
      )
    }
  }

  if (uses_mc && !is_gauss_exact) {
    em_trace <- list(
      icm = icm,
      p11 = p11,
      em_start = NULL,
      history = list()
    )
    class(em_trace) <- "population_mode_em_trace"
  }

  converged <- FALSE
  delta <- NA_real_
  em_iterations <- 0L
  tol_eff <- tol
  mc_delta_floor <- NA_real_
  em_route <- if (is_gauss_exact) "gaussian_exact" else "em"
  estep_out <- NULL

  if (is_gauss_exact) {
    fixef <- .c05_mean_map_lmerb(NULL, p11)
    estep_out <- .c05_estep(
      design = design,
      fixef = fixef,
      p11 = p11,
      measurement_prior_list = mpl,
      family = family,
      estep = estep,
      n = n,
      beta_set = beta_set
    )
    em_iterations <- 0L
    converged <- TRUE
    delta <- 0
  } else if (isTRUE(icm_screen) && isTRUE(icm_init)) {
    icm_screen_out <- .c05_icm_screen(
      design = design,
      fixef = fixef,
      p11 = p11,
      measurement_prior_list = mpl,
      family = family,
      estep = estep,
      n = n,
      mc_seed = mc_seed,
      beta_set = beta_set,
      tol = tol,
      mc_alpha = mc_alpha,
      mc_stable = mc_stable
    )
    estep_out <- icm_screen_out$estep_out
    estep_at_start <- estep_out
    mc_seed_used <- !is.null(mc_seed)
    delta <- icm_screen_out$delta
    mc_delta_floor <- icm_screen_out$mc_delta_floor
    tol_eff <- icm_screen_out$tol_eff
    if (isTRUE(icm_screen_out$passed)) {
      converged <- TRUE
      em_route <- "icm_screen"
    }
  }

  if (!is.null(em_trace)) {
    if (!is.null(estep_at_start)) {
      estep_for_start <- estep_at_start
    } else if (!is.null(icm_screen_out)) {
      estep_for_start <- icm_screen_out$estep_out
    } else {
      estep_for_start <- .c05_estep(
        design = design,
        fixef = fixef,
        p11 = p11,
        measurement_prior_list = mpl,
        family = family,
        estep = estep,
        n = n,
        mc_seed = mc_seed,
        beta_set = beta_set
      )
      estep_at_start <- estep_for_start
      mc_seed_used <- !is.null(mc_seed)
    }
    em_trace$em_start <- .c05_em_trace_record(
      iter = 0L,
      fixef = fixef,
      estep_out = estep_for_start,
      design = design,
      p11 = p11,
      family = family,
      measurement_prior_list = mpl,
      tol = tol,
      mc_alpha = mc_alpha,
      uses_mc = uses_mc
    )
    if (isTRUE(verbose)) {
      message("EM start (after ICM):")
      message(
        "  log_post = MC log pi(gamma|y) when joint safe draws are available"
      )
      .c05_verbose_em_trace(em_trace$em_start, p11, label = "  gamma_start")
    }
  }

  if (!converged) {
    em_route <- "em"
    prev_delta <- NA_real_
    prev_tol_eff <- NA_real_
    prev_g_norm <- if (!is.null(em_trace)) em_trace$em_start$delta else NA_real_
    prev_log_post <- if (!is.null(em_trace)) em_trace$em_start$log_post else NA_real_
    mc_stable <- as.integer(mc_stable)
    for (iter in seq_len(maxit)) {
      em_iterations <- iter
      reuse_start_estep <- isTRUE(iter == 1L && !is.null(estep_at_start))
      if (reuse_start_estep) {
        estep_out <- estep_at_start
      } else {
        estep_out <- .c05_estep(
          design = design,
          fixef = fixef,
          p11 = p11,
          measurement_prior_list = mpl,
          family = family,
          estep = estep,
          n = n,
          mc_seed = if (iter == 1L && !mc_seed_used) mc_seed else NULL,
          beta_set = beta_set
        )
      }

      if (!is.null(em_trace) && !reuse_start_estep) {
        trace_rec <- .c05_em_trace_record(
          iter = iter,
          fixef = fixef,
          estep_out = estep_out,
          design = design,
          p11 = p11,
          family = family,
          measurement_prior_list = mpl,
          tol = tol,
          mc_alpha = mc_alpha,
          uses_mc = uses_mc,
          prev_delta = prev_g_norm,
          prev_log_post = prev_log_post
        )
        em_trace$history[[length(em_trace$history) + 1L]] <- trace_rec
        if (isTRUE(verbose)) {
          .c05_verbose_em_trace(trace_rec, p11)
        }
        prev_log_post <- trace_rec$log_post
        prev_g_norm <- trace_rec$delta
      } else if (!is.null(em_trace) && reuse_start_estep) {
        prev_g_norm <- em_trace$em_start$delta
        prev_log_post <- em_trace$em_start$log_post
      }

      fixef_new <- .c05_mean_map(estep_out$b_mean, p11)

      delta <- .c05_gamma_delta(fixef, fixef_new, p11)
      mc_delta_floor <- if (uses_mc) {
        .c05_mc_delta_floor(estep_out$b_mc_se, p11)
      } else {
        NA_real_
      }
      tol_eff <- if (uses_mc) {
        .c05_em_tol(tol, "mc", estep_out$b_mc_se, p11, mc_alpha = mc_alpha)
      } else {
        tol
      }

      fixef <- fixef_new
      stable_hit <- delta <= tol_eff
      if (mc_stable <= 1L) {
        if (stable_hit) {
          converged <- TRUE
          break
        }
      } else if (stable_hit &&
                 is.finite(prev_delta) &&
                 prev_delta <= prev_tol_eff) {
        converged <- TRUE
        break
      }
      prev_tol_eff <- tol_eff
      prev_delta <- delta
    }
  }

  if (!converged) {
    msg <- paste0(
      "population_mode() did not converge in ", maxit,
      " iterations (final delta = ", signif(delta, 3L)
    )
    if (uses_mc && is.finite(mc_delta_floor)) {
      msg <- paste0(
        msg,
        ", effective tol = ", signif(tol_eff, 3L),
        " from MC n = ", n,
        ", mc_delta_floor = ", signif(mc_delta_floor, 3L), ")"
      )
    } else {
      msg <- paste0(msg, ").")
    }
    warning(msg, call. = FALSE)
  }

  jac <- .c05_jacobian(p11, estep_out)
  eps_star_closure <- .c05_epsilon_closure(jac$tilde_J)

  stationarity <- {
    fixef_check <- .c05_mean_map(estep_out$b_mean, p11)
    .c05_gamma_delta(fixef, fixef_check, p11)
  }

  list(
    fixef = fixef,
    gamma_star = .c05_gamma_from_fixef(fixef, p11),
    b_mean = estep_out$b_mean,
    b_mc_se = estep_out$b_mc_se,
    V_list = estep_out$V_list,
    icm = icm,
    icm_screen = icm_screen_out,
    em_route = em_route,
    em_trace = em_trace,
    design = design,
    family = family,
    measurement_prior_list = mpl,
    group_levels = prep$group_levels,
    p11 = p11,
    P11 = p11$P11,
    P11_RE = p11$P11_RE,
    Sigma_star = p11$Sigma_star,
    chol_P11 = p11$chol_P11,
    tilde_J = jac$tilde_J,
    kappa_spectrum = jac$kappa_spectrum,
    weights = jac$weights,
    kappa = jac$kappa,
    rho = jac$rho,
    Sigma_pi = jac$Sigma_pi,
    eps_star_closure = eps_star_closure,
    stationarity = stationarity,
    converged = converged,
    iterations = em_iterations,
    delta = delta,
    tol_eff = tol_eff,
    mc_delta_floor = mc_delta_floor,
    estep = estep_out$estep,
    n = estep_out$n,
    n_target = estep_out$n_target,
    n_tried = estep_out$n_tried,
    accept_rate = estep_out$accept_rate,
    beta_set = beta_set,
    restricted = isTRUE(estep_out$restricted),
    q = p11$q,
    call = match.call()
  )
}
