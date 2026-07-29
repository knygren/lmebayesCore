## Omega_j-MARGINALIZED (Section 16) per-group Hessian and system-wide
## lambda_star_marginal, evaluated at each pilot draw's own beta_j.
##
## Formalizes data-raw/_scratch_lambda_star_marginal_over_draws.R into an
## internal helper wired into .rLMMIngNormal_reg_run_with_pilot() (known_vcov
## engine only -- see inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 16 for the
## derivation and Section 16.5 for why this differs from the joint-extended
## system's own empirical over-draws scan). P11/P12/P21 are UNCHANGED from
## .two_block_S_P11(inp): they only depend on the RE-precision prior Lambda
## and the hyper-design, never on Omega_j; only the marginalized
## P22 = block-diag_j(Lambda + H_j(beta_j)) depends on the draw.

#' Resolve per-group ING measurement Gamma shape/rate as named vectors
#'
#' Accepts either a per-group \code{shape_group}/\code{rate_group} ING
#' measurement prior (named numeric vectors, one entry per group level) or a
#' pooled \code{shape}/\code{rate} scalar prior, broadcast to every group
#' level.
#' @noRd
.two_block_shape_rate_group_lookup <- function(ing_prior_list, group_levels) {
  if (!is.null(ing_prior_list$shape_group) || !is.null(ing_prior_list$rate_group)) {
    shape_group <- ing_prior_list$shape_group[group_levels]
    rate_group  <- ing_prior_list$rate_group[group_levels]
  } else {
    shape_group <- stats::setNames(
      rep(as.numeric(ing_prior_list$shape[1L]), length(group_levels)), group_levels
    )
    rate_group <- stats::setNames(
      rep(as.numeric(ing_prior_list$rate[1L]), length(group_levels)), group_levels
    )
  }
  list(shape_group = shape_group, rate_group = rate_group)
}

#' Per-group draw-independent setup for the Omega-marginalized Hessian
#'
#' Precomputes \eqn{D_j'D_j}, the group OLS fit, the Section 16.3
#' log-concavity threshold \eqn{2 r_j^0 + \mathrm{RSS}_{OLS,j}}, and the
#' group's ING Gamma prior shape/rate \eqn{a_j^0}/\eqn{r_j^0}, reused across
#' every draw.
#' @noRd
.two_block_marginal_group_setup <- function(D, y, group, group_levels, re_coef_names,
                                             shape_group, rate_group) {
  group_chr <- as.character(group)
  stats::setNames(lapply(group_levels, function(lev) {
    idx <- which(group_chr == lev)
    D_j <- D[idx, re_coef_names, drop = FALSE]
    y_j <- y[idx]
    DtD <- crossprod(D_j)
    beta_ols <- as.vector(solve(DtD, crossprod(D_j, y_j)))
    RSS_ols  <- sum((y_j - D_j %*% beta_ols)^2)
    list(
      rows = idx, D_j = D_j, DtD = DtD, n_j = length(idx),
      beta_ols = beta_ols, RSS_ols = RSS_ols,
      a0 = shape_group[[lev]], r0 = rate_group[[lev]],
      threshold = 2 * rate_group[[lev]] + RSS_ols
    )
  }), group_levels)
}

#' System-wide Omega-marginalized lambda_star over a set of draws
#'
#' For each draw, builds the per-group Omega_j-marginalized Hessian
#' \eqn{H_j(\beta_j)} (Section 16.2), checks that \eqn{\Lambda + H_j} is
#' positive definite for every group (skipping the draw if not, and flagging
#' the failing group(s)), and assembles the marginalized
#' \eqn{S = P_{12} P_{22,\mathrm{marginal}}^{-1} P_{21}} matrix.
#'
#' In addition to each draw's top eigenvalue, this accumulates a
#' componentwise \code{pmax} envelope of the full eigenvalue spectrum across
#' every PD-clean draw (mirroring \code{.two_block_pilot_ub_from_coefficients()}'s
#' plain-rate envelope trick), since \code{\link{two_block_tv_bound}} needs a
#' full eigenvalue vector -- not just \eqn{\lambda^*} -- to certify a TV
#' bound.
#'
#' @param coefficients Stacked \code{draw} / group-name / random-effect
#'   coefficient data frame (same layout as \code{pilot_raw$coefficients} /
#'   \code{fit$coefficients}).
#' @param n_draws Number of draws stacked in \code{coefficients}.
#' @param y,group_name As used elsewhere in this engine.
#' @param group_setup Output of \code{\link{.two_block_marginal_group_setup}}.
#' @param inp,blocks Output of \code{\link{.two_block_rate_inputs}} /
#'   \code{\link{.two_block_S_P11}} for the same design/prior.
#' @noRd
.two_block_lambda_star_marginal_over_draws <- function(coefficients, n_draws, y,
                                                        group_name, group_setup,
                                                        inp, blocks) {
  p_re <- inp$dims$p_re
  J    <- inp$dims$J
  q    <- inp$q
  P_b  <- blocks$P_b
  cols <- inp$gamma_cols
  co   <- coefficients
  re_coef_names <- inp$re_names
  group_levels  <- inp$group_levels
  co_group_chr  <- as.character(co[[group_name]])

  lambda_star_vec  <- rep(NA_real_, n_draws)
  skipped          <- logical(n_draws)
  failing_groups   <- vector("list", n_draws)
  h_violates_count <- stats::setNames(integer(length(group_levels)), group_levels)
  pd_fail_count    <- stats::setNames(integer(length(group_levels)), group_levels)
  max_eigenvalues  <- NULL
  lambda_star_max  <- -Inf
  i_max            <- NA_integer_
  n_over_one       <- 0L

  for (i in seq_len(n_draws)) {
    draw_rows <- which(co[["draw"]] == i)
    S <- matrix(0, q, q)
    ok <- TRUE
    fail_this <- character(0L)

    for (j in seq_len(J)) {
      lev <- group_levels[[j]]
      gs  <- group_setup[[lev]]
      row_i  <- draw_rows[co_group_chr[draw_rows] == lev]
      beta_j <- as.numeric(unlist(co[row_i, re_coef_names]))

      e_j <- y[gs$rows] - gs$D_j %*% beta_j
      ete <- sum(e_j^2)
      Dte <- as.vector(crossprod(gs$D_j, e_j))

      ## Section 16.3 raw criterion (H_j alone, ignoring Lambda) -- reported
      ## for comparison; q_j = e_j'e_j - RSS_ols_j exactly (Pythagorean OLS
      ## split), matching .tmp_rss_ellipsoid_test()'s q_j/threshold.
      q_j <- ete - gs$RSS_ols
      if (q_j > gs$threshold) {
        h_violates_count[[lev]] <- h_violates_count[[lev]] + 1L
      }

      Omega_eff <- (gs$a0 + gs$n_j / 2) / (gs$r0 + 0.5 * ete)
      H_j <- Omega_eff * gs$DtD -
        (Omega_eff^2 / (gs$a0 + gs$n_j / 2)) * outer(Dte, Dte)
      B_j <- P_b + H_j
      B_j <- 0.5 * (B_j + t(B_j))

      ch <- tryCatch(chol(B_j), error = function(e) NULL)
      if (is.null(ch)) {
        ok <- FALSE
        fail_this <- c(fail_this, lev)
        pd_fail_count[[lev]] <- pd_fail_count[[lev]] + 1L
        next
      }
      Binv <- chol2inv(ch)
      C_j  <- P_b %*% Binv %*% P_b

      x_j <- lapply(seq_len(p_re), function(k) inp$X_hyper[[k]][j, , drop = TRUE])
      for (a in seq_len(p_re)) {
        for (b in a:p_re) {
          out_ab <- outer(x_j[[a]], x_j[[b]])
          S[cols[[a]], cols[[b]]] <- S[cols[[a]], cols[[b]]] + C_j[a, b] * out_ab
          if (b > a) S[cols[[b]], cols[[a]]] <- t(S[cols[[a]], cols[[b]]])
        }
      }
    }

    if (!ok) {
      skipped[i] <- TRUE
      failing_groups[[i]] <- fail_this
      next
    }
    S  <- 0.5 * (S + t(S))
    ev <- .two_block_gen_eigen(S, blocks$P11, strict = FALSE)
    lambda_star_vec[i] <- ev[1L]
    if (is.null(max_eigenvalues)) {
      max_eigenvalues <- rep(-Inf, length(ev))
    }
    max_eigenvalues <- pmax(max_eigenvalues, ev)
    if (is.finite(ev[1L]) && ev[1L] >= 1) n_over_one <- n_over_one + 1L
    if (ev[1L] > lambda_star_max) {
      lambda_star_max <- ev[1L]
      i_max <- i
    }
  }

  list(
    lambda_star_vec  = lambda_star_vec,
    skipped          = skipped,
    n_skipped        = sum(skipped),
    failing_groups   = failing_groups,
    h_violates_count = h_violates_count,
    pd_fail_count    = pd_fail_count,
    lambda_star_max  = lambda_star_max,
    max_eigenvalues  = max_eigenvalues,
    i_max            = i_max,
    n_over_one       = n_over_one,
    n_draws          = n_draws
  )
}

#' Post-pilot Omega-marginalized \code{lambda_star} safeguard
#'
#' Computes the Section 16 Omega_j-marginalized system-wide
#' \code{lambda_star} at every pilot draw's own random-effect coefficients,
#' takes the componentwise \code{pmax} eigenvalue envelope across all
#' PD-clean draws, and reports whether that envelope is safe to use in place
#' of the conservative disp_upper plug-in \code{lambda_star}: valid only when
#' every pilot draw's \eqn{\Lambda + H_j(\beta_j)} block was positive
#' definite for every group AND the resulting \code{lambda_star_marginal}
#' stays below \code{cutoff}.
#'
#' This diagnostic assumes the random-effect precision \eqn{\Lambda} is
#' fixed/known, so it is only meaningful for
#' \code{rLMMindepNormalGamma_reg_known_vcov()} (\code{any_non_normal = FALSE}
#' in \code{.rLMMIngNormal_reg_run_with_pilot()}); see
#' \code{inst/BLOCK_GIBBS_ERGODICITY_ING.md} Section 16.
#'
#' @param D0 Squared standardized start distance for the main-stage TV
#'   certificate -- MUST match the \code{D0} used to derive the plain-rate
#'   main \code{m_convergence} (i.e.
#'   \code{two_block_d0_pilot_start(n_pilot, p_dim, pilot_start_tol)}), since
#'   the main-stage chains start from the pilot mean, not the exact mode.
#'   Passing \code{D0 = 0} would certify only a mode start and silently
#'   understate the sweeps needed. Default \code{0} is provided only for
#'   standalone/diagnostic use where no pilot-start recentering applies.
#' @noRd
.two_block_pilot_marginal_ub_from_coefficients <- function(
    pilot_coefficients,
    n_pilot,
    y, D, group, x_hyper,
    prior_list_block1,
    pfamily_list,
    tv_tol,
    group_name,
    group_levels,
    re_names,
    ing_prior_list,
    D0 = 0,
    cutoff = .two_block_lambda_star_slow_threshold()
) {
  sr <- .two_block_shape_rate_group_lookup(ing_prior_list, group_levels)

  re_names_x   <- names(x_hyper)
  pfamily_list <- .two_block_validate_pfamily_list(pfamily_list, re_names_x)
  prior_list_block2 <- lapply(pfamily_list, function(pf) {
    pl <- pf$prior_list
    list(
      mu    = pl$mu,
      Sigma = pl$Sigma,
      dispersion = if (identical(pf$pfamily, "dNormal")) pl$dispersion else pl$disp_lower
    )
  })

  inp <- .two_block_rate_inputs(
    x = D, block = group, x_hyper = x_hyper,
    prior_list_block1 = prior_list_block1,
    prior_list_block2 = prior_list_block2,
    family = gaussian(), group_levels = group_levels
  )
  blocks <- .two_block_S_P11(inp)

  group_setup <- .two_block_marginal_group_setup(
    D = D, y = y, group = group, group_levels = group_levels,
    re_coef_names = re_names,
    shape_group = sr$shape_group, rate_group = sr$rate_group
  )

  res <- .two_block_lambda_star_marginal_over_draws(
    coefficients = pilot_coefficients, n_draws = n_pilot, y = y,
    group_name = group_name, group_setup = group_setup, inp = inp, blocks = blocks
  )

  failing_groups <- sort(names(res$pd_fail_count)[res$pd_fail_count > 0L])
  valid <- res$n_skipped == 0L && is.finite(res$lambda_star_max) &&
    res$lambda_star_max < cutoff

  m_min_marginal   <- NA_integer_
  lambda_star_marg <- res$lambda_star_max
  rate_marginal    <- NULL

  if (valid) {
    ev_marginal <- res$max_eigenvalues
    rate_marginal <- structure(
      list(
        lambda_star  = ev_marginal[1L],
        eigenvalues  = ev_marginal,
        m_for_tol    = function(tol) {
          if (!is.numeric(tol) || length(tol) != 1L || tol <= 0 || tol >= 1) {
            stop("'tol' must be a single value in (0, 1).", call. = FALSE)
          }
          if (ev_marginal[1L] <= 0) return(1L)
          as.integer(ceiling(log(tol) / log(ev_marginal[1L])))
        },
        dims           = inp$dims,
        re_names       = inp$re_names,
        gamma_names    = inp$gamma_names,
        group_levels   = inp$group_levels,
        family         = inp$family,
        weights_source = "marginal"
      ),
      class = "two_block_rate"
    )
    m_min_marginal <- .two_block_cap_inner_sweeps(
      two_block_l_for_tv(
        rate_marginal, tv_tol, method = "theorem3", D0 = D0, warn = FALSE
      ) + 1L
    )
    lambda_star_marg <- rate_marginal$lambda_star
  }

  list(
    valid                = valid,
    lambda_star_marginal = lambda_star_marg,
    m_min_marginal       = m_min_marginal,
    rate_marginal        = rate_marginal,
    D0                   = D0,
    n_pilot              = n_pilot,
    n_skipped            = res$n_skipped,
    failing_groups       = failing_groups,
    pd_fail_count        = res$pd_fail_count,
    cutoff               = cutoff
  )
}
