## Extended (lambda, Omega)-aware local rate diagnostic ----------------------
##
## Implements inst/BLOCK_GIBBS_ERGODICITY_ING.md's extension of
## .two_block_S_P11()/two_block_rate() to a joint (gamma, beta, lambda, Omega)
## Hessian: Sections 5-6 for a sampled random-effects precision lambda_p
## (dIndependent_Normal_Gamma() on tau2_p), Section 14 for a sampled
## measurement precision Omega (pooled dGamma()) or Omega_j (per-group
## dGamma_list()). Diagnostic only -- see two_block_rate_ing()'s own
## documentation and Sections 7/9/15 of the theory note for the "no
## certified TV bound" caveat: this reports a *local* rate at a single
## reference state, not a global contraction guarantee the way
## two_block_rate()'s lambda* is (Nygren 2020, Remark 8).

#' Extend \code{.two_block_S_P11()}'s blocks with lambda-/Omega- rows and
#' diagonal entries (\code{inst/BLOCK_GIBBS_ERGODICITY_ING.md} Sections 6, 14)
#'
#' @param inp Output of \code{.two_block_rate_inputs()}.
#' @param blocks Output of \code{.two_block_S_P11(inp)}.
#' @param lambda_ing \code{NULL}, or a named list (names a subset of
#'   \code{inp$re_names}, one entry per ING random-effect component) with
#'   elements \code{lambda} (reference precision plug-in, scalar),
#'   \code{shape} (Gamma prior shape \eqn{a_p^0}), and \code{u}
#'   (numeric vector of RE-level residuals \eqn{u_{jp}}, named by group
#'   level, one entry per \code{inp$group_levels}).
#' @param omega_ing \code{NULL}, or a list with \code{scope} (\code{"pooled"}
#'   or \code{"per_group"}), \code{omega} (reference precision: scalar for
#'   \code{"pooled"}, named-by-group-level vector for \code{"per_group"}),
#'   \code{shape} (Gamma prior shape: scalar or named-by-group-level vector),
#'   \code{n} (observation count: scalar total for \code{"pooled"},
#'   named-by-group-level vector of per-group counts for \code{"per_group"}),
#'   and \code{e} (named-by-group-level list of data-level residual vectors
#'   \eqn{e_j = y_j - D_j\beta_j}, always per-group regardless of scope).
#' @return List with the extended \code{S}, \code{P11}, \code{P12}, \code{P21}
#'   (dimension \code{q_ext = q + n_lambda + n_omega}), \code{P22} (unchanged
#'   from \code{blocks$P22}), \code{ext_names}, \code{n_lambda}, \code{n_omega},
#'   \code{lambda_names}, \code{omega_scope}.
#' @keywords internal
.two_block_S_P11_ing <- function(inp, blocks, lambda_ing = NULL, omega_ing = NULL) {
  p_re <- inp$dims$p_re
  J    <- inp$dims$J
  q    <- inp$q

  lambda_names <- names(lambda_ing)
  n_lambda <- length(lambda_names)
  omega_scope <- if (is.null(omega_ing)) NULL else omega_ing$scope
  n_omega <- if (is.null(omega_ing)) {
    0L
  } else if (identical(omega_scope, "pooled")) {
    1L
  } else if (identical(omega_scope, "per_group")) {
    J
  } else {
    stop("omega_ing$scope must be \"pooled\" or \"per_group\".", call. = FALSE)
  }
  q_ext <- q + n_lambda + n_omega
  if (n_lambda == 0L && n_omega == 0L) {
    stop(
      "At least one of 'lambda_ing'/'omega_ing' must supply an ING component.",
      call. = FALSE
    )
  }

  ext_names <- c(
    inp$gamma_names,
    if (n_lambda > 0L) paste0("lambda::", lambda_names) else character(0L),
    if (identical(omega_scope, "pooled")) {
      "Omega"
    } else if (identical(omega_scope, "per_group")) {
      paste0("Omega::", inp$group_levels)
    } else {
      character(0L)
    }
  )
  b_names <- colnames(blocks$P12)

  P11 <- matrix(0, q_ext, q_ext)
  P11[seq_len(q), seq_len(q)] <- blocks$P11
  P12 <- matrix(0, q_ext, J * p_re)
  P12[seq_len(q), ] <- blocks$P12
  dimnames(P11) <- list(ext_names, ext_names)
  dimnames(P12) <- list(ext_names, b_names)

  ## --- lambda rows/diagonal (Section 6) --------------------------------------
  for (i in seq_along(lambda_names)) {
    k  <- lambda_names[i]
    li <- lambda_ing[[k]]
    row_i <- q + i
    k_idx <- match(k, inp$re_names)
    if (is.na(k_idx)) {
      stop("lambda_ing name '", k, "' is not one of inp$re_names.", call. = FALSE)
    }
    u_k <- li$u[inp$group_levels]
    if (anyNA(u_k)) {
      stop(
        "lambda_ing[[\"", k, "\"]]$u is missing an entry for some group level.",
        call. = FALSE
      )
    }
    P11[row_i, row_i] <- (J / 2 + li$shape - 1) / li$lambda^2

    Wk <- inp$X_hyper[[k]]
    cross <- -as.vector(crossprod(Wk, u_k))
    gk_cols <- inp$gamma_cols[[k_idx]]
    P11[gk_cols, row_i] <- cross
    P11[row_i, gk_cols] <- cross

    for (j in seq_len(J)) {
      bc <- (j - 1L) * p_re + k_idx
      P12[row_i, bc] <- u_k[[j]]
    }
  }

  ## --- Omega rows/diagonal (Section 14) --------------------------------------
  if (identical(omega_scope, "pooled")) {
    row_i <- q + n_lambda + 1L
    P11[row_i, row_i] <- (omega_ing$n / 2 + omega_ing$shape - 1) / omega_ing$omega^2
    for (j in seq_len(J)) {
      lev  <- inp$group_levels[[j]]
      rows <- inp$row_idx[[j]]
      D_j  <- inp$x[rows, , drop = FALSE]
      e_j  <- omega_ing$e[[lev]]
      if (length(e_j) != length(rows)) {
        stop(
          "omega_ing$e[[\"", lev, "\"]] must have length ", length(rows), ".",
          call. = FALSE
        )
      }
      bc <- (j - 1L) * p_re + seq_len(p_re)
      P12[row_i, bc] <- -as.vector(crossprod(D_j, e_j))
    }
  } else if (identical(omega_scope, "per_group")) {
    for (j in seq_len(J)) {
      lev   <- inp$group_levels[[j]]
      row_i <- q + n_lambda + j
      n_j     <- omega_ing$n[[lev]]
      shape_j <- if (length(omega_ing$shape) > 1L) omega_ing$shape[[lev]] else omega_ing$shape
      omega_j <- omega_ing$omega[[lev]]
      P11[row_i, row_i] <- (n_j / 2 + shape_j - 1) / omega_j^2

      rows <- inp$row_idx[[j]]
      D_j  <- inp$x[rows, , drop = FALSE]
      e_j  <- omega_ing$e[[lev]]
      if (length(e_j) != length(rows)) {
        stop(
          "omega_ing$e[[\"", lev, "\"]] must have length ", length(rows), ".",
          call. = FALSE
        )
      }
      bc <- (j - 1L) * p_re + seq_len(p_re)
      P12[row_i, bc] <- -as.vector(crossprod(D_j, e_j))
    }
  }

  ## P22 is block-diagonal (J * p_re small at demo scale); a direct solve()
  ## is simple and adequate for a once-per-fit diagnostic -- unlike
  ## .two_block_S_P11()'s own per-group accumulation, this is not on any
  ## sampler hot path.
  P11 <- 0.5 * (P11 + t(P11))
  P21 <- t(P12)
  S <- P12 %*% solve(blocks$P22, P21)
  S <- 0.5 * (S + t(S))
  dimnames(S)   <- list(ext_names, ext_names)
  dimnames(P21) <- list(b_names, ext_names)

  list(
    S = S, P11 = P11, P12 = P12, P21 = P21, P22 = blocks$P22,
    ext_names = ext_names, n_lambda = n_lambda, n_omega = n_omega,
    lambda_names = lambda_names, omega_scope = omega_scope
  )
}

#' Extended (lambda, Omega)-aware local rate diagnostic
#'
#' Computes the Remark~8-style rate matrix eigenvalues twice from the same
#' \code{prior_list_block1}/\code{prior_list_block2} inputs
#' \code{\link{two_block_rate}} would use: once as the plain
#' (\eqn{\gamma},\eqn{\beta})-only \code{eigenvalues_base}/\code{lambda_star_base}
#' (identical to what \code{two_block_rate()} itself would return for the same
#' inputs), and once with the joint \eqn{(\gamma,\beta,\lambda_{\mathrm{ING}},\Omega_{\mathrm{ING}})}
#' Hessian extension of \code{inst/BLOCK_GIBBS_ERGODICITY_ING.md} (Sections 5-6
#' for a sampled random-effect precision, Section 14 for a sampled measurement
#' precision) as \code{eigenvalues}/\code{lambda_star}.
#'
#' \strong{This is a local, uncertified diagnostic, not a new bound.} Unlike
#' \code{\link{two_block_rate}}'s \eqn{\lambda^*} (an exact, state-independent
#' contraction rate for the Gaussian two-block kernel, Nygren 2020 Remark~8),
#' the extended rate matrix here is evaluated at a single reference state
#' (the \code{lambda}/\code{omega}/\code{u}/\code{e} values supplied via
#' \code{lambda_ing}/\code{omega_ing}) and is genuinely state-dependent once
#' \eqn{\lambda} or \eqn{\Omega} is sampled (see Section 7 of the theory note).
#' It answers "how much would the currently-ignored \eqn{\beta}-\eqn{\lambda}/
#' \eqn{\beta}-\eqn{\Omega} coupling move the local rate at this reference
#' state?", not "what is a certified worst-case rate?" (that role is filled by
#' \code{\link{two_block_rate_from_pfamily_list}}'s existing corner plug-ins;
#' see Section 10 of the theory note).
#'
#' @inheritParams two_block_rate
#' @param lambda_ing \code{NULL} (no random-effects-precision extension), or a
#'   named list (names a subset of \code{names(x_hyper)}) as documented in
#'   \code{\link{.two_block_S_P11_ing}}.
#' @param omega_ing \code{NULL} (no measurement-precision extension), or a
#'   list as documented in \code{\link{.two_block_S_P11_ing}}.
#' @return Object of class \code{"two_block_rate_ing"}: a list with
#'   \code{lambda_star}/\code{eigenvalues} (extended system),
#'   \code{lambda_star_base}/\code{eigenvalues_base} (plain \eqn{(\gamma,\beta)}
#'   system, for direct comparison), \code{S}, \code{P11}, \code{P12},
#'   \code{P21}, \code{P22} (extended blocks), \code{ext_names}, \code{n_lambda},
#'   \code{n_omega}, \code{lambda_names}, \code{omega_scope}, \code{dims},
#'   \code{re_names}, \code{gamma_names}, \code{group_levels}, \code{family},
#'   \code{weights_source}, and \code{call}.
#' @seealso \code{\link{two_block_rate}}, \code{\link{two_block_rate_from_pfamily_list}}
#' @family simfuncs
#' @export
two_block_rate_ing <- function(x,
                                block,
                                x_hyper,
                                prior_list_block1,
                                prior_list_block2,
                                lambda_ing = NULL,
                                omega_ing = NULL,
                                weights = NULL,
                                family = gaussian(),
                                group_levels = levels(block),
                                warn_slow = TRUE) {
  cl <- match.call()
  if (is.null(lambda_ing) && is.null(omega_ing)) {
    stop(
      "At least one of 'lambda_ing'/'omega_ing' must be supplied ",
      "(otherwise call two_block_rate() directly).",
      call. = FALSE
    )
  }
  inp <- .two_block_rate_inputs(
    x = x, block = block, x_hyper = x_hyper,
    prior_list_block1 = prior_list_block1,
    prior_list_block2 = prior_list_block2,
    weights = weights, family = family, group_levels = group_levels
  )
  blocks <- .two_block_S_P11(inp)
  ev_base <- .two_block_gen_eigen(blocks$S, blocks$P11, blocks = blocks, inp = inp)
  lambda_star_base <- ev_base[1L]

  ext <- .two_block_S_P11_ing(
    inp, blocks, lambda_ing = lambda_ing, omega_ing = omega_ing
  )
  ## strict = FALSE: unlike two_block_rate()'s exact base system (where
  ## Remark 8 guarantees lambda_star < 1 and hitting the ceiling means a
  ## computation bug), this extended system is a *local* Hessian at a single
  ## reference state (lambda_ing/omega_ing's plugged-in residuals) -- see
  ## this function's own "local, uncertified diagnostic" documentation above
  ## and inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 11. lambda_star_ext >= 1
  ## is a legitimate finding (e.g. a reference-fit residual far from 0 for
  ## some group), not an error, so it is reported rather than raised.
  ev_ext <- .two_block_gen_eigen(ext$S, ext$P11, strict = FALSE)
  lambda_star_ext <- ev_ext[1L]

  if (is.finite(lambda_star_ext) && lambda_star_ext >= 1) {
    if (isTRUE(warn_slow)) {
      warning(
        "lambda* (extended) = ", signif(lambda_star_ext, 4),
        " >= 1: the local (beta, lambda/Omega) coupling at the supplied ",
        "reference state (residuals u/e) is strong enough that this ",
        "local diagnostic no longer indicates a contracting rate. This is ",
        "expected when some group's reference residual is large (see ",
        "inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 11); it does not, by ",
        "itself, mean the sampler's certified corner rate (lambda_star_base ",
        "above) is wrong.",
        call. = cl
      )
    }
  } else {
    .two_block_warn_lambda_star_slow(lambda_star_ext, warn = warn_slow, call = cl)
  }

  structure(
    list(
      lambda_star = lambda_star_ext,
      eigenvalues = ev_ext,
      lambda_star_base = lambda_star_base,
      eigenvalues_base = ev_base,
      ext_names = ext$ext_names,
      n_lambda = ext$n_lambda,
      n_omega = ext$n_omega,
      lambda_names = ext$lambda_names,
      omega_scope = ext$omega_scope,
      S = ext$S, P11 = ext$P11, P12 = ext$P12, P21 = ext$P21, P22 = ext$P22,
      dims = inp$dims,
      re_names = inp$re_names,
      gamma_names = inp$gamma_names,
      group_levels = inp$group_levels,
      family = inp$family,
      weights_source = inp$weights_source,
      call = cl
    ),
    class = "two_block_rate_ing"
  )
}

#' Print method for \code{two_block_rate_ing} objects
#'
#' @param x Object of class \code{"two_block_rate_ing"}.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @method print two_block_rate_ing
#' @export
print.two_block_rate_ing <- function(x, ...) {
  d <- x$dims
  cat("Extended (lambda, Omega)-aware LOCAL rate diagnostic\n")
  cat(
    "(NOT a certified TV bound -- a state-dependent, single-reference-state\n",
    "diagnostic; see inst/BLOCK_GIBBS_ERGODICITY_ING.md Sections 7, 9, 15)\n\n",
    sep = ""
  )
  cat(sprintf(
    "  groups J = %d, p_re = %d, q = %d (+ %d lambda, %d Omega), n_obs = %d  [family: %s]\n",
    d$J, d$p_re, d$q, x$n_lambda, x$n_omega, d$l2, x$family
  ))
  if (x$n_lambda > 0L) {
    cat("  lambda (RE precision) ING component(s): ",
        paste(x$lambda_names, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$omega_scope)) {
    cat("  Omega (measurement precision) scope: ", x$omega_scope, "\n", sep = "")
  }
  cat(sprintf("\nlambda* (base, no beta-(lambda,Omega) coupling) = %.6g\n",
              x$lambda_star_base))
  cat(sprintf("lambda* (extended, with coupling)               = %.6g\n",
              x$lambda_star))
  cat(sprintf("  delta = %.6g\n", x$lambda_star - x$lambda_star_base))
  invisible(x)
}

## Residual plug-ins for two_block_rate_ing(), built from the same lmer/
## glmmTMB reference fit the demo already fit for other purposes -- see
## inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 12 and the plan this
## implements ("Residual plug-ins", confirmed via .mer_re_reference_full()
## in lmebayes reducing coef(mer) - anchor to ranef(mer)).

#' RE-level residuals \eqn{u_{jp} = \beta_{jp} - W_{pj}\gamma_p} from a
#' reference fit, via \code{coef() - fixef()} per term
#'
#' Reduces to \code{lme4::ranef(fit)}/\code{glmmTMB::ranef(fit)} exactly
#' whenever every \code{groupef.names} entry is also a plain (non-interaction)
#' fixed-effect term name -- true for every model in
#' \code{demo/Ex_10}-\code{Ex_14} -- while routing through the existing
#' \code{\link{.lmebayes_reference_coef}}/\code{\link{.lmebayes_reference_fixef}}
#' dispatchers so the same code handles both an \code{lme4::lmer} and a
#' \code{glmmTMB::glmmTMB} reference fit.
#'
#' @param fit A fitted \code{merMod} or \code{glmmTMB} reference model.
#' @param group_name Name of the random-effects grouping factor.
#' @param groupef.names Random coefficient names to extract (a subset of the
#'   fixed-effect names in \code{fit}).
#' @return A data frame (rows = group levels, columns = \code{groupef.names})
#'   of \eqn{u_{jp}} residuals.
#' @keywords internal
#' @noRd
.lmebayes_reference_u <- function(fit, group_name, groupef.names) {
  co <- .lmebayes_reference_coef(fit)[[group_name]]
  if (is.null(co)) {
    stop(
      "Could not find grouping factor '", group_name, "' in coef(fit).",
      call. = FALSE
    )
  }
  missing_coefs <- setdiff(groupef.names, colnames(co))
  if (length(missing_coefs) > 0L) {
    stop(
      "coef(fit)[[\"", group_name, "\"]] is missing column(s): ",
      paste(missing_coefs, collapse = ", "),
      call. = FALSE
    )
  }
  fe <- .lmebayes_reference_fixef(fit)
  u  <- co[, groupef.names, drop = FALSE]
  for (k in groupef.names) {
    fe_k <- if (k %in% names(fe)) unname(fe[[k]]) else 0
    u[[k]] <- u[[k]] - fe_k
  }
  u
}

#' Data-level residuals \eqn{e_j = y_j - D_j\beta_j} from a reference fit's
#' per-group coefficients
#'
#' @param fit A fitted \code{merMod} or \code{glmmTMB} reference model.
#' @param y Response vector (as in \code{design$y}).
#' @param D Level-1 design matrix (as in \code{design$D}), \code{length(y) x
#'   length(groupef.names)}.
#' @param group Grouping factor, length \code{length(y)}, aligned to
#'   \code{y}/\code{D}.
#' @param group_name Name of the random-effects grouping factor.
#' @param groupef.names Column names of \code{D} (as in
#'   \code{design$groupef.names}).
#' @return Named list (one entry per group level in \code{coef(fit)}) of
#'   numeric residual vectors \eqn{e_j}.
#' @keywords internal
#' @noRd
.lmebayes_reference_group_residuals <- function(
    fit, y, D, group, group_name, groupef.names
) {
  co <- .lmebayes_reference_coef(fit)[[group_name]]
  if (is.null(co)) {
    stop(
      "Could not find grouping factor '", group_name, "' in coef(fit).",
      call. = FALSE
    )
  }
  missing_coefs <- setdiff(groupef.names, colnames(co))
  if (length(missing_coefs) > 0L) {
    stop(
      "coef(fit)[[\"", group_name, "\"]] is missing column(s): ",
      paste(missing_coefs, collapse = ", "),
      call. = FALSE
    )
  }
  group_chr <- as.character(group)
  levs <- rownames(co)
  out <- vector("list", length(levs))
  names(out) <- levs
  for (lev in levs) {
    rows <- which(group_chr == lev)
    if (length(rows) == 0L) {
      out[[lev]] <- numeric(0L)
      next
    }
    beta_j <- as.numeric(co[lev, groupef.names])
    out[[lev]] <- as.numeric(y[rows] - D[rows, groupef.names, drop = FALSE] %*% beta_j)
  }
  out
}

## Posterior-mean plug-ins for two_block_rate_ing(), built from a *completed*
## two-block Gibbs fit's own $coefficients/$dispersion_ranef.mean draws
## instead of an external lmer/glmmTMB reference fit -- for evaluating the
## Section 14-15 diagnostic at the sampler's own (self-consistent) posterior
## point estimate rather than the certified corner-rate's least-favorable
## disp_upper_group/fit_ref combination (which, for a group whose reference
## fit poorly predicts its own data, can put e_j and Omega_j badly out of
## step with each other; see the "det(H) < 0" derivation this fix followed
## from). Both omega and e below come from the SAME fitted object, so they
## describe one mutually-consistent reference state.

#' Posterior-mean per-group coefficients from a completed \code{rLMM_reg}/
#' \code{rGLMM_reg} fit's \code{$coefficients} draws
#'
#' @param fit A fitted object with a \code{coefficients} data frame (one row
#'   per (draw, group)), as returned by the matrix-level \code{rLMM_reg}/
#'   \code{rGLMM_reg} exports.
#' @param group_name Name of the grouping-factor column in
#'   \code{fit$coefficients}.
#' @param groupef.names Random coefficient names (columns of
#'   \code{fit$coefficients} and of \code{D}).
#' @return Named list (one entry per group level in \code{fit$coefficients})
#'   of posterior-mean coefficient vectors (length \code{length(groupef.names)}).
#' @keywords internal
#' @noRd
.lmebayes_posterior_mean_group_coef <- function(fit, group_name, groupef.names) {
  co <- fit$coefficients
  if (is.null(co) || is.null(co[[group_name]])) {
    stop(
      "fit$coefficients must contain a '", group_name, "' grouping column.",
      call. = FALSE
    )
  }
  missing_coefs <- setdiff(groupef.names, colnames(co))
  if (length(missing_coefs) > 0L) {
    stop(
      "fit$coefficients is missing column(s): ",
      paste(missing_coefs, collapse = ", "),
      call. = FALSE
    )
  }
  idx_by_group <- split(seq_len(nrow(co)), co[[group_name]])
  lapply(idx_by_group, function(idx) {
    colMeans(co[idx, groupef.names, drop = FALSE])
  })
}

#' RE-level residuals \eqn{u_{jp} = \bar\beta_{jp} - W_{pj}\bar\gamma_p} from
#' a completed fit's posterior-mean group/hyper coefficients
#'
#' @inheritParams .lmebayes_posterior_mean_group_coef
#' @param x_hyper Named list of group-level design matrices (as in
#'   \code{design$W}), one per \code{groupef.names} entry, row-named by
#'   group level.
#' @return A data frame (rows = group levels, columns = \code{groupef.names})
#'   of \eqn{u_{jp}} residuals, matching \code{.lmebayes_reference_u()}'s
#'   return shape.
#' @keywords internal
#' @noRd
.lmebayes_posterior_u <- function(fit, group_name, groupef.names, x_hyper) {
  beta_bar <- .lmebayes_posterior_mean_group_coef(fit, group_name, groupef.names)
  levs <- names(beta_bar)
  beta_mat <- do.call(rbind, beta_bar)
  rownames(beta_mat) <- levs
  gamma_bar <- stats::setNames(
    lapply(groupef.names, function(k) colMeans(fit$fixef[[k]])),
    groupef.names
  )
  u <- as.data.frame(beta_mat)
  for (k in groupef.names) {
    Wk <- x_hyper[[k]][levs, , drop = FALSE]
    u[[k]] <- beta_mat[, k] - as.vector(Wk %*% gamma_bar[[k]])
  }
  u
}

#' Data-level residuals \eqn{e_j = y_j - D_j\bar\beta_j} from a completed
#' fit's posterior-mean per-group coefficients
#'
#' @inheritParams .lmebayes_posterior_mean_group_coef
#' @param y Response vector (as in \code{design$y}).
#' @param D Level-1 design matrix (as in \code{design$D}).
#' @param group Grouping factor, length \code{length(y)}, aligned to
#'   \code{y}/\code{D}.
#' @return Named list (one entry per group level in \code{fit$coefficients})
#'   of numeric residual vectors \eqn{e_j}.
#' @keywords internal
#' @noRd
.lmebayes_posterior_group_residuals <- function(
    fit, y, D, group, group_name, groupef.names
) {
  beta_bar <- .lmebayes_posterior_mean_group_coef(fit, group_name, groupef.names)
  group_chr <- as.character(group)
  levs <- names(beta_bar)
  out <- vector("list", length(levs))
  names(out) <- levs
  for (lev in levs) {
    rows <- which(group_chr == lev)
    if (length(rows) == 0L) {
      out[[lev]] <- numeric(0L)
      next
    }
    out[[lev]] <- as.numeric(
      y[rows] - D[rows, groupef.names, drop = FALSE] %*% beta_bar[[lev]]
    )
  }
  out
}

## Empirical (not certified) worst-case rate over a completed fit's own
## main-stage draws -- mirrors .two_block_pilot_ub_from_coefficients()'s
## loop-and-pmax() pattern (R/two_block_glmm_pilot_helpers.R), but (a) runs
## over the *main*-stage draws already returned by rLMM*_reg()/rGLMM_reg()
## rather than a separate pilot batch, and (b) plugs each draw's own sampled
## RE-precision (fit$fixef.dispersion) and/or per-group measurement precision
## (fit$dispersion_ranef) into the *extended* two_block_rate_ing() system,
## rather than holding those at a fixed pilot-corner value and only letting
## beta vary. Treats the n returned draws as approximate posterior samples
## and asks: what is the largest local rate actually realized across them?
## This is a diagnostic only (inst/BLOCK_GIBBS_ERGODICITY_ING.md Sections
## 7, 9, 11, 15 still apply -- no certified TV bound for the extended
## system) and never feeds back into m_convergence/sampler calibration.

#' Slice a single main-stage draw out of a completed fit's
#' \code{$coefficients}/\code{$fixef} into a minimal pseudo-fit
#'
#' The pseudo-fit has exactly one draw's worth of rows/entries, so feeding it
#' to \code{.lmebayes_posterior_mean_group_coef()}/\code{.lmebayes_posterior_u()}/
#' \code{.lmebayes_posterior_group_residuals()} (which all reduce to
#' \code{colMeans()} over whatever rows they are given) returns that draw's
#' own values unchanged -- no separate per-draw residual code is needed.
#'
#' @param fit A fitted object with \code{$coefficients} (data frame with a
#'   \code{draw} column) and \code{$fixef} (named list of \code{n x p_k}
#'   matrices), as returned by the matrix-level \code{rLMM_reg}/
#'   \code{rGLMM_reg} exports.
#' @param i Integer draw index (\code{1..fit$n}).
#' @return List with \code{coefficients} (that draw's rows only) and
#'   \code{fixef} (that draw's row from each component, as 1-row matrices).
#' @keywords internal
#' @noRd
.lmebayes_fit_at_draw <- function(fit, i) {
  co <- fit$coefficients
  if (is.null(co) || is.null(co[["draw"]])) {
    stop("fit$coefficients must contain a 'draw' column.", call. = FALSE)
  }
  list(
    coefficients = co[co[["draw"]] == i, , drop = FALSE],
    fixef        = lapply(fit$fixef, function(m) m[i, , drop = FALSE])
  )
}

#' Empirical worst-case extended \code{(lambda, Omega)}-aware local rate over
#' a completed fit's main-stage draws
#'
#' Evaluates \code{\link{two_block_rate_ing}} once per main-stage draw
#' (treating the draws as approximate posterior samples), plugging in that
#' draw's own sampled RE precision (\code{1 / fit$fixef.dispersion[i, ]}, if
#' \code{lambda_spec} is supplied) and/or per-group measurement precision
#' (\code{1 / fit$dispersion_ranef[i, ]}, if \code{omega_spec} is supplied)
#' together with that draw's own \eqn{u_{jp}}/\eqn{e_j} residuals (via
#' \code{\link{.lmebayes_fit_at_draw}}). Tracks the pointwise maximum
#' eigenvalues and the largest \code{lambda_star} across all draws -- an
#' empirical analogue of \code{.two_block_pilot_ub_from_coefficients()}'s
#' pilot-draw scan, but for the extended system and the main-stage output.
#'
#' @param fit A fitted object as returned by the matrix-level
#'   \code{rLMM_reg}/\code{rGLMM_reg} exports (needs \code{$coefficients},
#'   \code{$fixef}, and -- as required by \code{lambda_spec}/\code{omega_spec}
#'   -- \code{$fixef.dispersion}/\code{$dispersion_ranef}).
#' @param n_draws Number of main-stage draws to scan (typically \code{fit$n}).
#' @param x,block,x_hyper,prior_list_block1,prior_list_block2,family,group_levels
#'   As in \code{\link{two_block_rate_ing}}; held fixed across draws (only
#'   \code{lambda_ing}/\code{omega_ing} vary per draw).
#' @param group_name,groupef.names As in \code{\link{.lmebayes_posterior_u}}.
#' @param y,D Response vector / level-1 design matrix, required when
#'   \code{omega_spec} is supplied (as in
#'   \code{\link{.lmebayes_posterior_group_residuals}}).
#' @param lambda_spec \code{NULL}, or a named list (subset of
#'   \code{groupef.names}) of \code{shape} values (Gamma prior shape
#'   \eqn{a_p^0}) -- one ING lambda entry per name.
#' @param omega_spec \code{NULL}, or a list with \code{scope}
#'   (\code{"pooled"}/\code{"per_group"}), \code{shape}, and \code{n}, as in
#'   \code{\link{.two_block_S_P11_ing}}'s \code{omega_ing} (everything except
#'   \code{omega}/\code{e}, which are rebuilt from each draw).
#' @return List with \code{lambda_star_vec} (length \code{n_draws}),
#'   \code{lambda_star_base} (unchanged across draws), \code{max_eigenvalues}
#'   (pointwise max of the extended spectrum across draws),
#'   \code{lambda_star_max}, \code{i_max} (draw achieving it), \code{n_over_one}
#'   (count of draws with \code{lambda_star >= 1}), and \code{n_draws}.
#' @keywords internal
#' @noRd
.two_block_rate_ing_over_draws <- function(
    fit, n_draws,
    x, block, x_hyper, prior_list_block1, prior_list_block2,
    group_name, groupef.names, y = NULL, D = NULL,
    lambda_spec = NULL, omega_spec = NULL,
    family = gaussian(), group_levels = levels(block)
) {
  if (is.null(lambda_spec) && is.null(omega_spec)) {
    stop(
      "At least one of 'lambda_spec'/'omega_spec' must be supplied.",
      call. = FALSE
    )
  }
  if (!is.null(omega_spec) && (is.null(y) || is.null(D))) {
    stop("'y' and 'D' are required when 'omega_spec' is supplied.", call. = FALSE)
  }

  lambda_names <- names(lambda_spec)
  lambda_star_vec  <- numeric(n_draws)
  lambda_star_base <- NA_real_
  max_eigenvalues  <- NULL
  lambda_star_max  <- -Inf
  i_max            <- NA_integer_
  n_over_one       <- 0L

  for (i in seq_len(n_draws)) {
    fit_i <- .lmebayes_fit_at_draw(fit, i)

    lambda_ing_i <- NULL
    if (!is.null(lambda_spec)) {
      u_i <- .lmebayes_posterior_u(fit_i, group_name, groupef.names, x_hyper)
      lambda_ing_i <- stats::setNames(lapply(lambda_names, function(k) {
        list(
          lambda = 1 / fit$fixef.dispersion[i, k],
          shape  = lambda_spec[[k]]$shape,
          u      = stats::setNames(u_i[[k]], rownames(u_i))
        )
      }), lambda_names)
    }

    omega_ing_i <- NULL
    if (!is.null(omega_spec)) {
      e_i <- .lmebayes_posterior_group_residuals(
        fit_i, y = y, D = D, group = block, group_name = group_name,
        groupef.names = groupef.names
      )
      omega_i <- 1 / fit$dispersion_ranef[i, group_levels]
      names(omega_i) <- group_levels
      omega_ing_i <- list(
        scope = omega_spec$scope,
        omega = omega_i,
        shape = omega_spec$shape,
        n     = omega_spec$n,
        e     = e_i
      )
    }

    rate_i <- two_block_rate_ing(
      x = x, block = block, x_hyper = x_hyper,
      prior_list_block1 = prior_list_block1,
      prior_list_block2 = prior_list_block2,
      lambda_ing = lambda_ing_i, omega_ing = omega_ing_i,
      family = family, group_levels = group_levels,
      warn_slow = FALSE
    )

    if (is.null(max_eigenvalues)) {
      max_eigenvalues <- rep(-Inf, length(rate_i$eigenvalues))
    }
    lambda_star_base   <- rate_i$lambda_star_base
    lambda_star_vec[i] <- rate_i$lambda_star
    max_eigenvalues    <- pmax(max_eigenvalues, rate_i$eigenvalues)
    if (is.finite(rate_i$lambda_star) && rate_i$lambda_star >= 1) {
      n_over_one <- n_over_one + 1L
    }
    if (rate_i$lambda_star > lambda_star_max) {
      lambda_star_max <- rate_i$lambda_star
      i_max <- i
    }
  }

  list(
    lambda_star_vec  = lambda_star_vec,
    lambda_star_base = lambda_star_base,
    max_eigenvalues  = max_eigenvalues,
    lambda_star_max  = lambda_star_max,
    i_max            = i_max,
    n_over_one       = n_over_one,
    n_draws          = n_draws
  )
}
