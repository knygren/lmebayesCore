#' Joint posterior mean or mode of the two-block mixed model (ICM)
#'
#' @description
#' Iterated conditional means/modes (ICM) for the two-block posterior targeted
#' by \code{\link{two_block_rNormal_reg}} and \code{\link{rGLMM_reg}}.  Block~2
#' hyperparameters \eqn{\gamma} and Block~1 random effects \eqn{b} are updated
#' alternately until \code{fixef} stabilizes.
#'
#' @details
#' \strong{Shared Block~2 update.}
#' For each RE component \eqn{k}:
#' \deqn{
#'   E[\gamma_k \mid b_k] =
#'   \bigl(X_k^\top X_k / \tau^2_k + P_{\gamma_k}\bigr)^{-1}
#'   \bigl(X_k^\top b_k / \tau^2_k + P_{\gamma_k} \mu_{\gamma_k}\bigr)
#' }
#' where \eqn{b_k} is the \eqn{k}-th column of the current Block~1 matrix,
#' \eqn{X_k =} \code{design$X_hyper[[k]]},
#' \eqn{\tau^2_k =} \code{dispersion_fixef}, and
#' \eqn{P_{\gamma_k} = \Sigma_{\gamma_k}^{-1}} from \code{prior_list}.
#'
#' \strong{Block~1 update} differs by function; see \code{\link{lmerb_posterior_mean}}
#' (closed-form Gaussian) vs \code{\link{glmerb_posterior_mode}}
#' (\code{\link{rglmb}} mode for general GLMM families).
#'
#' When the response is Gaussian and variance components are fixed, the joint
#' posterior is multivariate normal, so \code{glmerb_posterior_mode()} with
#' \code{family = gaussian()} matches \code{lmerb_posterior_mean()}.
#'
#' @param design Design list with \code{y}, \code{Z}, \code{groups},
#'   \code{X_hyper}, and \code{re_coef_names}.
#' @param measurement_prior_list List with \code{Sigma_ranef} and
#'   \code{prior_list}.  \code{dispersion_ranef} (\eqn{\sigma^2}) is required
#'   for \code{lmerb_posterior_mean()} and for \code{glmerb_posterior_mode()}
#'   when \code{family = gaussian()}; omit for non-Gaussian GLMM families.
#'   Each \code{prior_list[[k]]} must contain \code{mu_fixef},
#'   \code{Sigma_fixef}, and \code{dispersion_fixef}.
#' @param tol Convergence tolerance on the \eqn{\ell_\infty} change in
#'   \code{fixef} between successive iterations.  Default \code{1e-10}.
#' @param maxit Maximum number of ICM iterations.  Default \code{200L}.
#' @return A list with components \code{fixef}, \code{b_mean}, \code{converged},
#'   \code{iterations}, and \code{delta}.
#' @seealso \code{\link{build_mu_all}}, \code{\link{two_block_rNormal_reg}},
#'   \code{\link{rglmb}}
#' @name lmebayes_posterior_icm
#' @aliases lmerb_posterior_mean glmerb_posterior_mode
NULL

#' @describeIn lmebayes_posterior_icm Joint posterior \emph{mean} of the
#'   two-block Gaussian model (= joint mode when variance components are fixed).
#'   Block~1 uses the closed-form conditional mean update per group \eqn{j}:
#'   \deqn{
#'     E[b_j \mid \gamma] =
#'     \bigl(Z_j^\top Z_j / \sigma^2 + P_b\bigr)^{-1}
#'     \bigl(Z_j^\top y_j / \sigma^2 + P_b \,\mu_j(\gamma)\bigr)
#'   }
#'   with \eqn{P_b = \Sigma_b^{-1}} and \eqn{\mu_j(\gamma)} from
#'   \code{\link{build_mu_all}}.
#' @export
lmerb_posterior_mean <- function(design,
                                 measurement_prior_list,
                                 tol   = 1e-10,
                                 maxit = 200L) {

  .lmerb_validate_design(design)
  .lmerb_validate_measurement_prior_list(measurement_prior_list)

  if (is.null(design$y) || is.null(design$Z)) {
    stop("'design' must contain 'y' and 'Z'.", call. = FALSE)
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  J            <- length(group_levels)
  p_re         <- length(re_names)
  g_chr        <- as.character(design$groups)

  sigma2 <- measurement_prior_list$dispersion_ranef
  if (is.null(sigma2)) {
    stop(
      "'measurement_prior_list' must contain 'dispersion_ranef' ",
      "for lmerb_posterior_mean().",
      call. = FALSE
    )
  }
  P_b    <- solve(measurement_prior_list$Sigma_ranef)

  P_gamma  <- stats::setNames(
    lapply(re_names, function(k) {
      solve(measurement_prior_list$prior_list[[k]]$Sigma_fixef)
    }),
    re_names
  )
  mu_gamma <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$prior_list[[k]]$mu_fixef),
    re_names
  )
  tau2 <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$prior_list[[k]]$dispersion_fixef),
    re_names
  )

  ZtZ_scaled <- vector("list", J)
  Zty_scaled <- vector("list", J)
  names(ZtZ_scaled) <- names(Zty_scaled) <- group_levels

  for (lev in group_levels) {
    rows <- which(g_chr == lev)
    Z_j  <- design$Z[rows, , drop = FALSE]
    y_j  <- design$y[rows]
    ZtZ_scaled[[lev]] <- crossprod(Z_j) / sigma2
    Zty_scaled[[lev]] <- crossprod(Z_j, y_j) / sigma2
  }

  fixef <- lapply(measurement_prior_list$prior_list, `[[`, "mu_fixef")
  names(fixef) <- re_names

  b_mean <- matrix(
    0.0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )

  converged <- FALSE
  delta     <- NA_real_

  for (iter in seq_len(maxit)) {

    mu_all <- as.matrix(build_mu_all(design, fixef)$mu_all)

    for (jj in seq_len(J)) {
      lev      <- group_levels[jj]
      mu_j     <- mu_all[, jj]
      post_P_j <- ZtZ_scaled[[lev]] + P_b
      post_v_j <- Zty_scaled[[lev]] + P_b %*% mu_j
      b_mean[jj, ] <- solve(post_P_j, post_v_j)
    }

    fixef_new <- vector("list", p_re)
    names(fixef_new) <- re_names

    for (k in re_names) {
      X_k      <- design$X_hyper[[k]]
      b_k      <- b_mean[, k]
      tau2_k   <- tau2[[k]]
      P_gam_k  <- P_gamma[[k]]
      mu_gam_k <- mu_gamma[[k]]

      post_P_k  <- crossprod(X_k) / tau2_k + P_gam_k
      post_v_k  <- crossprod(X_k, b_k) / tau2_k + P_gam_k %*% mu_gam_k
      gam_k <- as.vector(solve(post_P_k, post_v_k))
      names(gam_k) <- colnames(X_k)
      fixef_new[[k]] <- gam_k
    }

    delta <- max(vapply(re_names, function(k) {
      max(abs(fixef_new[[k]] - fixef[[k]]))
    }, numeric(1L)))

    fixef <- fixef_new

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "lmerb_posterior_mean() did not converge in ", maxit, " iterations ",
      "(final delta = ", signif(delta, 3L), "). ",
      "Consider increasing 'maxit' or checking model identifiability.",
      call. = FALSE
    )
  }

  list(
    fixef      = fixef,
    b_mean     = b_mean,
    converged  = converged,
    iterations = iter,
    delta      = delta
  )
}

#' @describeIn lmebayes_posterior_icm Joint posterior \emph{mode} of the
#'   two-block GLMM.  Block~1 uses \code{\link{rglmb}} with \code{n = 1L} and a
#'   \code{\link{dNormal}} prior per group; the mode is read from
#'   \code{coef.mode}.  For \code{family = gaussian()}, this matches the
#'   closed-form update in \code{\link{lmerb_posterior_mean}}.
#' @param family A \code{\link[stats]{family}} object. Defaults to
#'   \code{gaussian()}.
#' @export
glmerb_posterior_mode <- function(design,
                                  family = gaussian(),
                                  measurement_prior_list,
                                  tol   = 1e-10,
                                  maxit = 200L) {

  .lmerb_validate_design(design)
  .lmerb_validate_measurement_prior_list(measurement_prior_list)

  if (is.null(design$y) || is.null(design$Z)) {
    stop("'design' must contain 'y' and 'Z'.", call. = FALSE)
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  J            <- length(group_levels)
  p_re         <- length(re_names)
  g_chr        <- as.character(design$groups)

  sigma2   <- measurement_prior_list$dispersion_ranef
  Sigma_b  <- measurement_prior_list$Sigma_ranef
  is_gaussian <- identical(family$family, "gaussian")
  if (is_gaussian && is.null(sigma2)) {
    stop(
      "'measurement_prior_list' must contain 'dispersion_ranef' ",
      "when family = gaussian().",
      call. = FALSE
    )
  }
  if (!is_gaussian && !is.null(sigma2)) {
    stop(
      "'measurement_prior_list$dispersion_ranef' must be NULL ",
      "for non-Gaussian families.",
      call. = FALSE
    )
  }

  P_gamma  <- stats::setNames(
    lapply(re_names, function(k) {
      solve(measurement_prior_list$prior_list[[k]]$Sigma_fixef)
    }),
    re_names
  )
  mu_gamma <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$prior_list[[k]]$mu_fixef),
    re_names
  )
  tau2 <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$prior_list[[k]]$dispersion_fixef),
    re_names
  )

  fixef <- lapply(measurement_prior_list$prior_list, `[[`, "mu_fixef")
  names(fixef) <- re_names

  b_mean <- matrix(
    0.0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )

  converged <- FALSE
  delta     <- NA_real_

  for (iter in seq_len(maxit)) {

    mu_all <- as.matrix(build_mu_all(design, fixef)$mu_all)

    for (jj in seq_len(J)) {
      lev  <- group_levels[jj]
      rows <- which(g_chr == lev)
      y_j  <- design$y[rows]
      Z_j  <- design$Z[rows, , drop = FALSE]
      mu_j <- mu_all[, jj]

      pf_j <- if (is.null(sigma2)) {
        dNormal(mu = mu_j, Sigma = Sigma_b)
      } else {
        dNormal(mu = mu_j, Sigma = Sigma_b, dispersion = sigma2)
      }
      fit_j <- rglmb(
        n       = 1L,
        y       = y_j,
        x       = Z_j,
        family  = family,
        pfamily = pf_j,
        verbose = FALSE
      )
      b_mean[jj, ] <- fit_j$coef.mode
    }

    fixef_new <- vector("list", p_re)
    names(fixef_new) <- re_names

    for (k in re_names) {
      X_k      <- design$X_hyper[[k]]
      b_k      <- b_mean[, k]
      tau2_k   <- tau2[[k]]
      P_gam_k  <- P_gamma[[k]]
      mu_gam_k <- mu_gamma[[k]]

      post_P_k  <- crossprod(X_k) / tau2_k + P_gam_k
      post_v_k  <- crossprod(X_k, b_k) / tau2_k + P_gam_k %*% mu_gam_k
      gam_k <- as.vector(solve(post_P_k, post_v_k))
      names(gam_k) <- colnames(X_k)
      fixef_new[[k]] <- gam_k
    }

    delta <- max(vapply(re_names, function(k) {
      max(abs(fixef_new[[k]] - fixef[[k]]))
    }, numeric(1L)))

    fixef <- fixef_new

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "glmerb_posterior_mode() did not converge in ", maxit, " iterations ",
      "(final delta = ", signif(delta, 3L), "). ",
      "Consider increasing 'maxit' or checking model identifiability.",
      call. = FALSE
    )
  }

  list(
    fixef      = fixef,
    b_mean     = b_mean,
    converged  = converged,
    iterations = iter,
    delta      = delta
  )
}

#' @noRd
.lmerb_validate_measurement_prior_list <- function(mpl) {
  if (!is.list(mpl)) {
    stop("'measurement_prior_list' must be a list.", call. = FALSE)
  }
  for (nm in c("Sigma_ranef", "prior_list")) {
    if (is.null(mpl[[nm]])) {
      stop(
        "'measurement_prior_list' must contain '", nm, "'.",
        call. = FALSE
      )
    }
  }
  invisible(mpl)
}
