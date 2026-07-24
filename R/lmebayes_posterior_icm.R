#' Joint posterior mean or mode of the two-block mixed model (ICM)
#'
#' @description
#' Joint Block~1/Block~2 posterior location for the two-block posterior
#' targeted by \code{\link{two_block_rNormal_reg}} and \code{\link{rGLMM_reg}}.
#' \code{lmerb_posterior_mean()} solves for the \emph{exact} joint Gaussian
#' mean in closed form (no iteration).  \code{glmerb_posterior_mode()} uses
#' iterated conditional modes (ICM): Block~2 hyperparameters \eqn{\gamma} and
#' Block~1 random effects \eqn{b} are updated alternately until \code{fixef}
#' stabilizes.
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
#' (exact closed-form Gaussian solve) vs \code{\link{glmerb_posterior_mode}}
#' (ICM using the \code{\link[glmbayesCore]{rglmb}} mode for general GLMM families).
#'
#' When the response is Gaussian and variance components are fixed, the joint
#' posterior is multivariate normal, so \code{glmerb_posterior_mode()} with
#' \code{family = gaussian()} targets the same mean as
#' \code{lmerb_posterior_mean()} (via ICM rather than the closed form).
#'
#' @param design Design list with \code{y}, \code{Z}, \code{groups},
#'   \code{X_hyper}, and \code{re_coef_names}.
#' @param measurement_prior_list List with \code{Sigma_ranef} and
#'   \code{prior_list}.  \code{dispersion_ranef} (\eqn{\sigma^2}) is required
#'   for \code{lmerb_posterior_mean()} and for \code{glmerb_posterior_mode()}
#'   when \code{family = gaussian()}; omit for non-Gaussian GLMM families.
#'   Each \code{prior_list[[k]]} must contain \code{mu_fixef},
#'   \code{Sigma_fixef}, and \code{dispersion_fixef}.
#' @param tol Convergence tolerance for \code{glmerb_posterior_mode()}'s ICM
#'   loop, on the change in \code{fixef} between successive iterations,
#'   measured as a Mahalanobis distance in each RE component's own
#'   posterior-precision metric (\eqn{\sqrt{(\gamma_k^{new} - \gamma_k)^\top
#'   P_{\gamma_k}^{\mathrm{post}} (\gamma_k^{new} - \gamma_k)}}, maximized
#'   over components \eqn{k}), not a raw coordinate-wise \eqn{\ell_\infty}
#'   change.  This makes convergence invariant to rescaling or whitening any
#'   \code{X_hyper[[k]]} column.  Default \code{1e-10}.  Accepted but unused
#'   by \code{lmerb_posterior_mean()}, which solves exactly (see Details) and
#'   always returns \code{converged = TRUE}, \code{iterations = 1L},
#'   \code{delta = 0}; kept for interface parity with
#'   \code{glmerb_posterior_mode()}.
#' @param maxit Maximum number of ICM iterations for
#'   \code{glmerb_posterior_mode()}.  Default \code{200L}.  Accepted but
#'   unused by \code{lmerb_posterior_mean()} (see \code{tol}).
#' @return A list with components \code{fixef}, \code{b_mean}, \code{converged},
#'   \code{iterations}, and \code{delta} (the Mahalanobis-distance stopping
#'   statistic described under \code{tol}; always \code{0} for
#'   \code{lmerb_posterior_mean()}).
#' @seealso \code{\link{build_mu_all}}, \code{\link{two_block_rNormal_reg}},
#'   \code{\link[glmbayesCore]{rglmb}}
#' @name lmebayes_posterior_icm
#' @aliases lmerb_posterior_mean glmerb_posterior_mode
NULL

#' @describeIn lmebayes_posterior_icm Joint posterior \emph{mean} of the
#'   two-block Gaussian model (= joint mode when variance components are
#'   fixed), computed \emph{exactly} in closed form -- no ICM iteration.
#'   Block~1's conditional mean per group \eqn{j},
#'   \deqn{
#'     E[b_j \mid \gamma] =
#'     \bigl(Z_j^\top Z_j / \sigma^2 + P_b\bigr)^{-1}
#'     \bigl(Z_j^\top y_j / \sigma^2 + P_b \,\mu_j(\gamma)\bigr),
#'   }
#'   is affine in \eqn{\gamma} (via \eqn{\mu_j(\gamma)} from
#'   \code{\link{build_mu_all}}), and \eqn{b_j} couples only to the shared
#'   \eqn{\gamma} -- never to another group's \eqn{b_{j'}}.  Substituting this
#'   affine relationship into the Block~2 update (see Details) eliminates all
#'   \eqn{b_j} algebraically, leaving one linear system in \eqn{\gamma} alone
#'   (dimension = total hyperparameter count, independent of the number of
#'   groups \eqn{J}); solving it once gives the exact joint mean, and one
#'   back-substitution pass gives \code{b_mean}.  Forming this system costs
#'   \eqn{O(J)} (one \eqn{p_{re} \times p_{re}} solve per group plus a small
#'   accumulation), never a \eqn{J \times J} or \eqn{J p_{re}}-dimensional
#'   matrix, so it scales to large numbers of groups.
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

  system <- .lmerb_posterior_normal_system(design, measurement_prior_list)
  gamma_full <- as.vector(solve(system$M, system$v))

  fixef <- stats::setNames(
    lapply(system$re_names, function(k) {
      g <- gamma_full[system$idx[[k]]]
      names(g) <- colnames(design$X_hyper[[k]])
      g
    }),
    system$re_names
  )

  b <- .lmerb_posterior_b_given_gamma(system, design, fixef)

  list(
    fixef      = fixef,
    b_mean     = b$b,
    converged  = TRUE,
    iterations = 1L,
    delta      = 0
  )
}

#' Exact posterior covariance of the stacked Block~2 hyperparameter vector
#'
#' @description
#' Companion to \code{\link{lmerb_posterior_mean}}: instead of solving for the
#' posterior mean of the stacked Block~2 hyperparameter vector
#' \code{gamma_full}, returns its exact posterior \emph{covariance}
#' \eqn{\Sigma_{11} = M^{-1}}, where \code{M} (the Schur-complement system's
#' posterior precision of \code{gamma_full}) is built by the package-internal
#' \code{.lmerb_posterior_normal_system()}. This is the \eqn{\Sigma_{11}}
#' of Claim~3 in the two-block Gibbs ergodicity reference (see
#' \code{\link{plot_sweep_history_var_ratio}}): the exact target covariance
#' that a two-block Gibbs sampler's cross-chain covariance converges to as the
#' number of inner sweeps grows.
#'
#' Only meaningful when \code{measurement_prior_list$dispersion_ranef} and
#' \code{measurement_prior_list$Sigma_ranef} are both \emph{fixed} (not
#' sampled) -- i.e. the same restriction as
#' \code{.lmerb_posterior_normal_system()} itself. For models with
#' estimated dispersion or estimated random-effect variance components, no
#' single exact \eqn{\Sigma_{11}} exists (it would vary by posterior draw);
#' use the empirical cross-chain covariance instead (see
#' \code{\link{plot_sweep_history_var_ratio}}'s fallback).
#'
#' @inheritParams lmerb_posterior_mean
#' @return A \code{P_total x P_total} covariance matrix (\code{P_total} = the
#'   total Block~2 hyperparameter count, summed over RE components), dimnamed
#'   \code{"re_component | covariate"} in the same stacking order as the
#'   package-internal \code{.two_block_snapshot_fixef_cov()}'s \code{coef_index}.
#' @seealso \code{\link{lmerb_posterior_mean}}
#' @export
lmerb_posterior_covariance <- function(design, measurement_prior_list) {
  .lmerb_validate_design(design)
  .lmerb_validate_measurement_prior_list(measurement_prior_list)

  system <- .lmerb_posterior_normal_system(design, measurement_prior_list)
  Sigma  <- solve(system$M)

  lbl <- unlist(lapply(system$re_names, function(k) {
    paste(k, colnames(design$X_hyper[[k]]), sep = " | ")
  }), use.names = FALSE)
  dimnames(Sigma) <- list(lbl, lbl)
  Sigma
}

#' Shared Schur-complement linear system for the exact two-block Gaussian posterior
#'
#' Builds the pieces of the exact joint Gaussian posterior over
#' \code{(gamma, b_1, ..., b_J)} that do not depend on any particular value
#' of \code{gamma}: the \eqn{P_{total} \times P_{total}} linear system
#' \code{M \%*\% gamma_full = v} (\code{M} is the posterior \emph{precision}
#' of the stacked Block~2 hyperparameter vector \code{gamma_full}, \code{v}
#' its precision-weighted mean contribution) and, per group \eqn{j}, the
#' Block~1 conditional precision \code{post_P_j} and \code{Zty_scaled}
#' needed to recover \eqn{E[b_j \mid \gamma]} for any \code{gamma} via
#' \code{.lmerb_posterior_b_given_gamma()}. Derived by eliminating
#' \eqn{b_1, ..., b_J} algebraically from the joint posterior (see
#' \code{\link{lmerb_posterior_mean}} Details) -- cost \code{O(J)}, no
#' \eqn{J \times J} or \eqn{J p_{re}}-dimensional matrix is ever formed.
#'
#' Shared by \code{lmerb_posterior_mean()} (solves \code{M \%*\% gamma = v}
#' once for the exact posterior mean) and
#' \code{\link{rLMMNormal_joint_iid}} (Cholesky-factors \code{M} and each
#' \code{post_P_j} once, then draws \code{n} iid samples from the same
#' system -- no Gibbs sweeps, no burn-in, since the target is exactly
#' Gaussian).
#' @noRd
.lmerb_posterior_normal_system <- function(design, measurement_prior_list) {
  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  J            <- length(group_levels)
  p_re         <- length(re_names)
  g_chr        <- as.character(design$groups)

  sigma2 <- measurement_prior_list$dispersion_ranef
  if (is.null(sigma2)) {
    stop(
      "'measurement_prior_list' must contain 'dispersion_ranef' ",
      "for the exact two-block Gaussian posterior.",
      call. = FALSE
    )
  }
  sigma2 <- as.numeric(sigma2)
  if (!(length(sigma2) %in% c(1L, J))) {
    stop(
      "'measurement_prior_list$dispersion_ranef' must have length 1 or the ",
      "number of groups (", J, ").",
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

  for (jj in seq_len(J)) {
    lev  <- group_levels[jj]
    rows <- which(g_chr == lev)
    Z_j  <- design$Z[rows, , drop = FALSE]
    y_j  <- design$y[rows]
    ## Per-group fixed dispersion vector (mode = "fixed_vector") is already
    ## ordered to match group_levels by the resolver in mixed_rmerb_helpers.R;
    ## the scalar "fixed" mode broadcasts the same value to every group.
    sigma2_j <- if (length(sigma2) > 1L) sigma2[[jj]] else sigma2
    ZtZ_scaled[[lev]] <- crossprod(Z_j) / sigma2_j
    Zty_scaled[[lev]] <- crossprod(Z_j, y_j) / sigma2_j
  }

  ## --- Exact joint posterior via block (Schur-complement) elimination.
  ##
  ## The target is exactly jointly Gaussian, so there is no need to alternate
  ## Block~1/Block~2 conditional means to a tolerance: eliminate Block~1
  ## algebraically and solve one small linear system for gamma_full instead.
  ##
  ## Block~1 (per group j): b_j = a_j + D_j %*% mu_j(gamma), with
  ##   post_P_j = ZtZ_scaled[[j]] + P_b  (independent of gamma),
  ##   a_j      = solve(post_P_j, Zty_scaled[[j]]),
  ##   D_j      = solve(post_P_j, P_b).
  ## mu_j(gamma) is linear in the stacked hyperparameter vector gamma_full,
  ## so b_j is affine in gamma_full.  Substituting into the Block~2 update
  ## for each RE component k (b_j only couples to the *shared* gamma_full,
  ## never to another group's b_j', so this elimination is exact and costs
  ## O(J) -- no J x J or J*p_re-dimensional matrix is ever formed):
  ##
  ##   post_P_gamma_k %*% gamma_k - sum_k2 (T_{k,k2} / tau2_k) %*% gamma_k2
  ##     = s_k / tau2_k + P_gamma_k %*% mu_gamma_k
  ##
  ## with post_P_gamma_k = crossprod(X_k) / tau2_k + P_gamma_k,
  ##      s_k            = crossprod(X_k, a_k)          (a_k = a_j[k] over j),
  ##      T_{k,k2}       = crossprod(X_k, d_{k,k2} * X_k2), d_{k,k2}[j] = D_j[k, k2].
  ##
  ## Stacking over k gives one P_total x P_total system (P_total = sum of
  ## hyperparameter counts, independent of J); solving it once yields the
  ## *exact* gamma_full, then one back-substitution pass gives b_mean.  No
  ## iteration, tolerance, or non-convergence is possible for this model.
  ## M is exactly the posterior *precision* of gamma_full, and each
  ## post_P_j the posterior precision of b_j | gamma -- this is what makes
  ## the system directly re-usable for iid sampling, not just the mean.
  post_P_j_list <- stats::setNames(vector("list", J), group_levels)
  a_list        <- stats::setNames(vector("list", J), group_levels)
  D_list        <- stats::setNames(vector("list", J), group_levels)

  for (lev in group_levels) {
    post_P_j <- ZtZ_scaled[[lev]] + P_b
    post_P_j_list[[lev]] <- post_P_j
    a_list[[lev]] <- as.vector(solve(post_P_j, Zty_scaled[[lev]]))
    D_list[[lev]] <- solve(post_P_j, P_b)
  }

  p_k     <- vapply(re_names, function(k) ncol(design$X_hyper[[k]]), integer(1L))
  P_total <- sum(p_k)
  idx     <- stats::setNames(vector("list", p_re), re_names)
  cursor  <- 0L
  for (k in re_names) {
    idx[[k]] <- seq.int(cursor + 1L, cursor + p_k[[k]])
    cursor   <- cursor + p_k[[k]]
  }

  post_P_gamma_list <- stats::setNames(
    lapply(re_names, function(k) {
      crossprod(design$X_hyper[[k]]) / tau2[[k]] + P_gamma[[k]]
    }),
    re_names
  )

  a_mat <- vapply(group_levels, function(lev) a_list[[lev]], numeric(p_re))
  if (p_re == 1L) a_mat <- matrix(a_mat, nrow = 1L)

  M <- matrix(0.0, P_total, P_total)
  v <- numeric(P_total)

  for (ik in seq_len(p_re)) {
    k    <- re_names[[ik]]
    X_k  <- design$X_hyper[[k]]
    a_k  <- a_mat[ik, ]
    s_k  <- as.vector(crossprod(X_k, a_k))

    v[idx[[k]]] <- s_k / tau2[[k]] + as.vector(P_gamma[[k]] %*% mu_gamma[[k]])
    M[idx[[k]], idx[[k]]] <- M[idx[[k]], idx[[k]]] + post_P_gamma_list[[k]]

    for (ik2 in seq_len(p_re)) {
      k2    <- re_names[[ik2]]
      X_k2  <- design$X_hyper[[k2]]
      d_vec <- vapply(group_levels, function(lev) D_list[[lev]][ik, ik2], numeric(1L))
      T_kk2 <- crossprod(X_k, X_k2 * d_vec)
      M[idx[[k]], idx[[k2]]] <- M[idx[[k]], idx[[k2]]] - T_kk2 / tau2[[k]]
    }
  }

  list(
    M             = M,
    v             = v,
    idx           = idx,
    p_k           = p_k,
    P_total       = P_total,
    post_P_j_list = post_P_j_list,
    Zty_scaled    = Zty_scaled,
    P_b           = P_b,
    re_names      = re_names,
    group_levels  = group_levels,
    J             = J,
    p_re          = p_re
  )
}

#' Block~1 conditional mean given a Block~2 hyperparameter vector
#'
#' Back-substitution step of the exact two-block Gaussian posterior: given
#' any \code{fixef} (a named list of Block~2 hyperparameter vectors, shaped
#' like \code{lmerb_posterior_mean()}'s own \code{fixef} return value --
#' either the exact posterior mean, or one iid draw from
#' \code{\link{rLMMNormal_joint_iid}}), returns \eqn{E[b_j \mid \gamma]} for
#' every group \eqn{j}, using the \code{post_P_j} / \code{Zty_scaled} /
#' \code{P_b} already built by \code{.lmerb_posterior_normal_system()}.
#' @noRd
.lmerb_posterior_b_given_gamma <- function(system, design, fixef) {
  group_levels <- system$group_levels
  J    <- system$J
  p_re <- system$p_re

  mu_all <- as.matrix(
    build_mu_all(design, fixef, group_levels = group_levels)$mu_all
  )
  b <- matrix(
    0.0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, system$re_names)
  )
  for (jj in seq_len(J)) {
    lev      <- group_levels[jj]
    mu_j     <- mu_all[, jj]
    post_v_j <- system$Zty_scaled[[lev]] + system$P_b %*% mu_j
    b[jj, ]  <- solve(system$post_P_j_list[[lev]], post_v_j)
  }

  list(b = b, mu_all = mu_all)
}

#' Cholesky factors of the exact two-block Gaussian posterior, for iid sampling
#'
#' \code{M} (posterior precision of \code{gamma_full}) and each
#' \code{post_P_j} (posterior precision of \eqn{b_j \mid \gamma}) are always
#' \emph{mathematically} symmetric -- but \code{.lmerb_posterior_normal_system()}
#' assembles \code{M} one Block~2 component-row at a time using each
#' component's own \eqn{\tau^2_k} (see that function's Details), which is
#' only guaranteed to reproduce a symmetric \code{M} when
#' \code{measurement_prior_list$Sigma_ranef} is diagonal \emph{and} its
#' \eqn{k}-th diagonal entry equals component \eqn{k}'s \code{dispersion_fixef}
#' (\eqn{\tau^2_k}) -- i.e. \code{Sigma_ranef = diag(tau2_k)}, exactly how
#' \code{lmerb()}/\code{glmerb()} always construct it (see
#' \code{lmebayes::lmerb} Details), but not a precondition enforced at this
#' matrix level. \code{solve(M, v)} (the posterior \emph{mean}, used by
#' \code{\link{lmerb_posterior_mean}}) is unaffected by this -- it is a
#' correct linear-system solve regardless of \code{M}'s symmetry. Drawing
#' \emph{samples} via \code{chol(M)}, however, is only valid when \code{M}
#' actually is (the precision of) a symmetric covariance, so this check is
#' required before any iid sampling (\code{\link{rLMMNormal_joint_iid}}).
#'
#' @param tol Maximum tolerated relative asymmetry
#'   \code{max(abs(M - t(M))) / max(abs(M))} before erroring. Default
#'   \code{1e-6} (well above floating-point roundoff, which is \code{~1e-16},
#'   but far below the asymmetry produced by a genuinely inconsistent
#'   \code{Sigma_ranef}/\code{pfamily_list} pairing).
#' @return List with \code{R_M} (\code{chol(M)}, upper triangular) and
#'   \code{R_j_list} (named list of \code{chol(post_P_j)} per group level).
#' @noRd
.lmerb_posterior_system_cholesky <- function(system, tol = 1e-6) {
  M <- system$M
  asym <- max(abs(M - t(M)))
  scale <- max(abs(M))
  if (scale > 0 && asym / scale > tol) {
    stop(
      "rLMMNormal_joint_iid(): the Block~2 posterior precision is not ",
      "symmetric (relative asymmetry ", signif(asym / scale, 3), "). ",
      "Exact iid sampling requires 'Sigma_ranef' (from 'prior_list_block1') ",
      "to be diagonal with its k-th entry exactly equal to component k's ",
      "dNormal() dispersion (tau2_k) in 'pfamily_list' -- i.e. ",
      "Sigma_ranef = diag(tau2_k), as lmerb()/glmerb() always construct it. ",
      "Use sim_method = \"TWO_BLOCK_GIBBS\" (or call ",
      "rLMMNormal_reg_known_vcov_two_bg()/two_block_rNormal_reg() directly) ",
      "for a 'Sigma_ranef' that is correlated or inconsistent with 'pfamily_list'.",
      call. = FALSE
    )
  }
  M_sym <- 0.5 * (M + t(M))

  R_M <- chol(M_sym)
  R_j_list <- stats::setNames(
    lapply(system$group_levels, function(lev) chol(system$post_P_j_list[[lev]])),
    system$group_levels
  )

  list(R_M = R_M, R_j_list = R_j_list)
}

#' @describeIn lmebayes_posterior_icm Joint posterior \emph{mode} of the
#'   two-block GLMM.  Block~1 uses \code{\link[glmbayesCore]{rglmb}} with \code{n = 1L} and a
#'   \code{\link[glmbayesCore]{dNormal}} prior per group; the mode is read from
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
  if (!is.null(sigma2)) {
    sigma2 <- as.numeric(sigma2)
    if (!(length(sigma2) %in% c(1L, J))) {
      stop(
        "'measurement_prior_list$dispersion_ranef' must have length 1 or ",
        "the number of groups (", J, ").",
        call. = FALSE
      )
    }
  }
  ## Per-group fixed dispersion vector (mode = "fixed_vector") is already
  ## ordered to match group_levels by the resolver in mixed_rmerb_helpers.R;
  ## the scalar "fixed"/dGamma() modes broadcast the same value to every
  ## group. Hoisted once since sigma2 is fixed/known across ICM iterations.
  sigma2_per_group <- if (is.null(sigma2)) {
    rep(list(NULL), J)
  } else if (length(sigma2) > 1L) {
    as.list(sigma2)
  } else {
    rep(list(sigma2), J)
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

  ## Block-2 posterior precisions are fixed across ICM iterations, so hoist
  ## them once and reuse both for the linear solve and for the scale-invariant
  ## convergence check below (see lmerb_posterior_mean()).
  post_P_list <- stats::setNames(
    lapply(re_names, function(k) {
      crossprod(design$X_hyper[[k]]) / tau2[[k]] + P_gamma[[k]]
    }),
    re_names
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

      sigma2_j <- sigma2_per_group[[jj]]
      pf_j <- if (is.null(sigma2_j)) {
        glmbayesCore::dNormal(mu = mu_j, Sigma = Sigma_b)
      } else {
        glmbayesCore::dNormal(mu = mu_j, Sigma = Sigma_b, dispersion = sigma2_j)
      }
      fit_j <- glmbayesCore::rglmb(
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

      post_v_k <- crossprod(X_k, b_k) / tau2_k + P_gam_k %*% mu_gam_k
      gam_k    <- as.vector(solve(post_P_list[[k]], post_v_k))
      names(gam_k) <- colnames(X_k)
      fixef_new[[k]] <- gam_k
    }

    ## Scale-invariant stopping rule (Mahalanobis distance); see
    ## lmerb_posterior_mean() for the rationale.
    delta <- sqrt(max(vapply(re_names, function(k) {
      d <- fixef_new[[k]] - fixef[[k]]
      as.numeric(crossprod(d, post_P_list[[k]] %*% d))
    }, numeric(1L))))

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
