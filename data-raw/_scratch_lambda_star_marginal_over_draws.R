## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## "Gradual development" step 1: a system-wide lambda_star built from the
## exact Omega_j-MARGINALIZED (Section 16.2) per-group Hessian, instead of
## the fixed-Omega_j plug-in `two_block_rate()` uses (the "base" lambda_star)
## or the joint-state extension `two_block_rate_ing()` uses (the "extended"
## lambda_star, which differentiates but never integrates Omega_j out):
##
##   H_j(beta_j) = Omega_j^eff * D_j'D_j
##                 - (Omega_j^eff)^2 / (a_j^0 + n_j/2) * (D_j'e_j)(D_j'e_j)'
##   Omega_j^eff(beta_j) = (a_j^0 + n_j/2) / (r_j^0 + 0.5 * e_j'e_j)
##   e_j = y_j - D_j * beta_j
##
## a_j^0/r_j^0 are the group's Block~1 ING Gamma PRIOR shape/rate (the same
## shape_group/rate_group already used by _scratch_rss_ellipsoid_test.R's
## threshold = 2*r0_j + RSS_ols_j).
##
## See inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 16 for the derivation and
## Section 16.5 for why this is a genuinely different diagnostic from the
## *joint* extended system's own empirical over-draws scan
## (.two_block_rate_ing_over_draws(), used by demo("Ex_13_...") etc.).
##
## Evaluated at EVERY main-stage draw's own beta_j (like
## .tmp_rss_ellipsoid_test()'s per-draw q_j check), not a single reference
## state:
##
##   P11, P12, P21      -- UNCHANGED from .two_block_S_P11(inp): they only
##                          depend on the RE-precision prior Lambda (= P_b)
##                          and the hyper-design, never on Omega_j.
##   P22_marginal^(i)    -- block-diag_j( Lambda + H_j(beta_j^(i)) ), the one
##                          thing that changes per draw.
##
## "Gradual development" plan (do NOT force this further yet): for each draw
## AND group, first check whether Lambda + H_j(beta_j^(i)) is itself
## positive-definite via chol() -- a strictly MILDER condition than Section
## 16.3's raw "H_j alone" PSD criterion (reported alongside, for comparison,
## as h_violates_count), since the RE-precision prior Lambda can absorb a
## mild H_j violation. If ANY group's Lambda + H_j block fails chol() for a
## draw, that draw's lambda_star_marginal is SKIPPED (NA) and the failing
## group(s) are flagged in pd_fail_count / failing_groups -- no attempt is
## made to force an eigendecomposition through a non-PD block. Only draws
## where every group's block is PD get a real lambda_star_marginal.

## ---------------------------------------------------------------------------
## Per-group, draw-INDEPENDENT setup (D_j'D_j, OLS fit, Section 16.3
## threshold, and the ING prior a_j^0/r_j^0) -- computed once, reused across
## all n_draws. Mirrors .tmp_rss_ellipsoid_test()'s per-group block exactly
## (data-raw/_scratch_rss_ellipsoid_test.R), so h_violates_count below should
## reproduce that script's pct_draws_outside per group as a sanity check.
## ---------------------------------------------------------------------------
.tmp_marginal_group_setup <- function(D, y, group, group_levels, re_coef_names,
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

## ---------------------------------------------------------------------------
## Main loop: one system-wide lambda_star_marginal per main-stage draw.
##
## `inp`/`blocks` are the SAME objects two_block_rate_ing()/two_block_rate()
## build internally -- pass in the output of
##   inp    <- lmebayesCore:::.two_block_rate_inputs(x = D, block = group,
##               x_hyper = W, prior_list_block1 = prior_list_block1_rate,
##               prior_list_block2 = prior_list_block2_rate)
##   blocks <- lmebayesCore:::.two_block_S_P11(inp)
## (identical inputs to the demo's own `rate_ext <- two_block_rate_ing(...)`
## call, minus lambda_ing/omega_ing) so P11/P12/P21/P_b are guaranteed
## consistent with the "base"/"extended" numbers already reported.
## ---------------------------------------------------------------------------
.tmp_lambda_star_marginal_over_draws <- function(fit, n_draws, y, group_name,
                                                   group_setup, inp, blocks) {
  p_re <- inp$dims$p_re
  J    <- inp$dims$J
  q    <- inp$q
  P_b  <- blocks$P_b
  cols <- inp$gamma_cols
  co   <- fit$coefficients
  re_coef_names <- inp$re_names
  group_levels  <- inp$group_levels
  co_group_chr  <- as.character(co[[group_name]])

  lambda_star_vec  <- rep(NA_real_, n_draws)
  skipped          <- logical(n_draws)
  failing_groups   <- vector("list", n_draws)
  h_violates_count <- stats::setNames(integer(length(group_levels)), group_levels)
  pd_fail_count    <- stats::setNames(integer(length(group_levels)), group_levels)
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
      row_i <- draw_rows[co_group_chr[draw_rows] == lev]
      beta_j <- as.numeric(unlist(co[row_i, re_coef_names]))

      e_j  <- y[gs$rows] - gs$D_j %*% beta_j
      ete  <- sum(e_j^2)
      Dte  <- as.vector(crossprod(gs$D_j, e_j))

      ## Section 16.3 raw criterion (H_j alone, ignoring Lambda) -- reported
      ## for comparison against .tmp_rss_ellipsoid_test()'s q_j/threshold
      ## (q_j = e_j'e_j - RSS_ols_j exactly, via the Pythagorean OLS split).
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
      C_j <- P_b %*% Binv %*% P_b

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
    S <- 0.5 * (S + t(S))
    ev <- lmebayesCore:::.two_block_gen_eigen(S, blocks$P11, strict = FALSE)
    lambda_star_vec[i] <- ev[1L]
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
    failing_groups    = failing_groups,
    h_violates_count = h_violates_count,
    pd_fail_count    = pd_fail_count,
    lambda_star_max  = lambda_star_max,
    i_max            = i_max,
    n_over_one       = n_over_one,
    n_draws          = n_draws
  )
}

## ---------------------------------------------------------------------------
## Print a summary in the same spirit as the ellipsoid test's cat() report.
## ---------------------------------------------------------------------------
.tmp_print_marginal_over_draws_summary <- function(res) {
  n_ok <- res$n_draws - res$n_skipped
  cat(sprintf(
    "\n=== Omega-MARGINALIZED lambda_star over %d main-stage draws ===\n\n",
    res$n_draws
  ))
  cat(sprintf(
    "  draws skipped (some group's Lambda + H_j not PD): %d / %d (%.1f%%)\n",
    res$n_skipped, res$n_draws, 100 * res$n_skipped / res$n_draws
  ))
  if (res$n_skipped > 0L) {
    pd_tab <- res$pd_fail_count[res$pd_fail_count > 0L]
    pd_tab <- sort(pd_tab, decreasing = TRUE)
    cat("  groups causing a non-PD Lambda + H_j block (count of draws):\n")
    for (nm in names(pd_tab)) {
      cat(sprintf("    %-6s  %d\n", nm, pd_tab[[nm]]))
    }
  }
  h_tab <- res$h_violates_count[res$h_violates_count > 0L]
  h_tab <- sort(h_tab, decreasing = TRUE)
  cat(sprintf(
    "\n  Section 16.3 raw H_j (Lambda-free) violations, %d groups nonzero:\n",
    length(h_tab)
  ))
  for (nm in names(h_tab)) {
    cat(sprintf("    %-6s  %d / %d draws (%.1f%%)\n",
                nm, h_tab[[nm]], res$n_draws, 100 * h_tab[[nm]] / res$n_draws))
  }
  if (n_ok > 0L) {
    lv <- res$lambda_star_vec[!res$skipped]
    cat(sprintf(
      paste0(
        "\n  Among the %d draws with every group PD:\n",
        "    lambda_star_marginal: mean = %.6f, median = %.6f, max = %.6f (draw #%d)\n",
        "    draws with lambda_star_marginal >= 1: %d / %d\n"
      ),
      n_ok, mean(lv), stats::median(lv), res$lambda_star_max, res$i_max,
      res$n_over_one, n_ok
    ))
  } else {
    cat("\n  No draws had every group PD -- no lambda_star_marginal computed.\n")
  }
  invisible(res)
}

## ---------------------------------------------------------------------------
## Usage: with an existing fit already in your session (Ex_13/Ex_13b/Ex_13c's
## variable names: fit, design, grp, re_names, shape_group, rate_group,
## prior_list_block1_rate, prior_list_block2_rate, n_draws), just call:
##
##   inp    <- lmebayesCore:::.two_block_rate_inputs(
##     x = design$D, block = grp, x_hyper = design$W,
##     prior_list_block1 = prior_list_block1_rate,
##     prior_list_block2 = prior_list_block2_rate
##   )
##   blocks <- lmebayesCore:::.two_block_S_P11(inp)
##   group_setup <- .tmp_marginal_group_setup(
##     D = design$D, y = design$y, group = grp, group_levels = group_levels,
##     re_coef_names = re_names, shape_group = shape_group, rate_group = rate_group
##   )
##   res <- .tmp_lambda_star_marginal_over_draws(
##     fit = fit, n_draws = n_draws, y = design$y,
##     group_name = design$group_name, group_setup = group_setup,
##     inp = inp, blocks = blocks
##   )
##   .tmp_print_marginal_over_draws_summary(res)
##
## No demo re-run needed -- this file only *defines* the functions above;
## sourcing it does not execute anything else.
## ---------------------------------------------------------------------------
