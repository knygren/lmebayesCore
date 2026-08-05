test_that("Prior_Setup_lmebayes: group.dispersion overrides pooled dispersion_ranef", {
  dat <- lme4::sleepstudy

  ps_default <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01
  )

  override <- as.numeric(ps_default$dispersion_ranef) * 1.25
  ps <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01,
    group.dispersion = override
  )

  expect_equal(ps$dispersion_ranef, override, ignore_attr = TRUE)
  expect_equal(ps$group.dispersion, override)
  expect_identical(attr(ps$dispersion_ranef, "source"), "user group.dispersion")
  expect_equal(ps$group.ing_prior$sigma2_hat, override)
  p_re <- length(ps$design$groupef.names)
  n_prior <- ps$group.ing_prior$n_prior
  expect_equal(
    ps$group.ing_prior$rate,
    override * (n_prior + p_re - 1) / 2
  )
  expect_identical(
    ps$pop.prior_list,
    ps_default$pop.prior_list
  )
})

test_that("Prior_Setup_lmebayes: group.dispersion vector overrides per-group dispersion", {
  skip_if_not_installed("glmmTMB")
  dat <- lme4::sleepstudy
  J <- nlevels(dat$Subject)

  ps_default <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01,
    group.dispersion.pwt = 0.1,
    dispformula = ~Subject
  )

  override <- ps_default$sigma2_group * 1.1
  ps <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01,
    group.dispersion.pwt = 0.1,
    dispformula = ~Subject,
    group.dispersion = override
  )

  expect_length(ps$dispersion_ranef, J)
  expect_named(ps$dispersion_ranef, levels(dat$Subject))
  expect_equal(
    unname(ps$dispersion_ranef),
    unname(override),
    ignore_attr = TRUE
  )
  expect_equal(ps$group.dispersion, override)
  expect_identical(attr(ps$dispersion_ranef, "source"), "user group.dispersion")
  expect_identical(ps$sigma2_group, ps_default$sigma2_group)

  p_re <- length(ps$design$groupef.names)
  for (lev in levels(dat$Subject)) {
    ing_j <- ps$group.ing_prior[[lev]]
    expect_equal(ing_j$sigma2_hat, unname(override[[lev]]))
    expect_equal(
      ing_j$rate_gamma,
      unname(override[[lev]]) * (ing_j$n_prior + p_re - 1) / 2
    )
  }
})

test_that("Prior_Setup_lmebayes: group.dispersion validation", {
  dat <- lme4::sleepstudy

  expect_error(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      group.dispersion = 1,
      dispformula = ~Subject
    ),
    "got a scalar"
  )

  expect_error(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      family = poisson(),
      group.dispersion = 1
    ),
    "only supported for family = gaussian"
  )
})

test_that("Prior_Setup_lmebayes: dispformula = ~1 keeps the pooled lme4 reference", {
  dat <- lme4::sleepstudy

  ps <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt  = 0.01
  )

  expect_identical(ps$calibration_source, "lme4")
  expect_s4_class(ps$fit_ref, "merMod")
  expect_identical(ps$fit_ref, ps$mer_fit)
  expect_identical(ps$fit_ref, ps$design$lmer)
  expect_null(ps$dispersion_fit)
  expect_null(ps$sigma2_group)
  expect_false(is.null(ps$group.ing_prior))
  expect_false(lmebayesCore:::.lmebayes_ing_prior_is_grouped(ps$group.ing_prior))
})

test_that("Prior_Setup_lmebayes: dispformula = ~group routes calibration through glmmTMB", {
  skip_if_not_installed("glmmTMB")
  dat <- lme4::sleepstudy

  ps <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data            = dat,
    pop.pwt         = 0.01,
    group.dispersion.pwt       = 0.1,
    dispformula     = ~Subject
  )

  expect_identical(ps$calibration_source, "glmmTMB")
  expect_s3_class(ps$fit_ref, "glmmTMB")
  expect_s4_class(ps$mer_fit, "merMod")
  expect_identical(ps$dispersion_fit, ps$fit_ref)
  expect_identical(ps$mer_fit, ps$design$lmer)

  ## Pooled dispersion_ranef stays lme4-derived (design$dispersion from
  ## mer_fit) regardless of dispformula.
  expect_equal(
    unname(ps$dispersion_ranef),
    unname(ps$design$dispersion),
    ignore_attr = TRUE
  )

  ## group.ing_prior holds either the pooled (flat) or per-group (named
  ## list of lists) calibration, selected by dispformula: only the
  ## per-group shape is calibrated here.
  expect_false(is.null(ps$group.ing_prior))
  expect_true(lmebayesCore:::.lmebayes_ing_prior_is_grouped(ps$group.ing_prior))
  expect_length(ps$group.ing_prior, nlevels(dat$Subject))

  expect_type(ps$sigma2_group, "double")
  expect_named(ps$sigma2_group, levels(dat$Subject))

  disp_pf <- dGamma_list(ps, warn_asymmetric = FALSE)
  expect_length(disp_pf, nlevels(dat$Subject))
  expect_true(all(vapply(disp_pf, inherits, logical(1L), "pfamily")))
  expect_identical(attr(disp_pf, "dispersion_fit"), ps$dispersion_fit)
  expect_identical(attr(disp_pf, "calibration_source"), "glmmTMB")
})

test_that("Prior_Setup_lmebayes: pop.mu supports partial (NULL-entry) overrides", {
  dat <- lme4::sleepstudy

  ps_default <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )

  ps_mu <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.mu = list("(Intercept)" = c("(Intercept)" = 250), Days = NULL)
  )

  ## Overridden component: mu_fixef is exactly the supplied value.
  expect_equal(
    unname(ps_mu$pop.prior_list[["(Intercept)"]]$mu_fixef),
    250
  )

  ## Non-overridden component (NULL entry): unchanged from the usual
  ## pop.intercept_source/pop.effects_source-derived default.
  expect_identical(
    ps_mu$pop.prior_list[["Days"]]$mu_fixef,
    ps_default$pop.prior_list[["Days"]]$mu_fixef
  )

  ## Sigma_fixef is untouched by pop.mu.
  expect_identical(
    ps_mu$pop.prior_list[["(Intercept)"]]$Sigma_fixef,
    ps_default$pop.prior_list[["(Intercept)"]]$Sigma_fixef
  )
})

test_that("Prior_Setup_lmebayes: pop.sd converts to the expected Sigma_fixef scaling", {
  dat <- lme4::sleepstudy

  sd_val <- 50
  ps_sd <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.sd = list("(Intercept)" = sd_val, Days = NULL)
  )
  ps_default <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )

  ## For a single-column component, Sigma_fixef = V_ii * (1-w)/w with
  ## w = V_ii / (V_ii + sd^2) collapses to exactly sd^2, independent of V_ii.
  expect_equal(
    ps_sd$pop.prior_list[["(Intercept)"]]$Sigma_fixef[1, 1],
    sd_val^2
  )

  ## Non-overridden component (NULL entry): unchanged pop.pwt-derived scaling.
  expect_identical(
    ps_sd$pop.prior_list[["Days"]]$Sigma_fixef,
    ps_default$pop.prior_list[["Days"]]$Sigma_fixef
  )

  ## mu_fixef is untouched by pop.sd.
  expect_identical(
    ps_sd$pop.prior_list[["(Intercept)"]]$mu_fixef,
    ps_default$pop.prior_list[["(Intercept)"]]$mu_fixef
  )
})

test_that("Prior_Setup_lmebayes: pop.pwt, pop.sd, and pop.nprior are mutually exclusive", {
  dat <- lme4::sleepstudy

  expect_error(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.pwt = 0.02,
      pop.sd  = list("(Intercept)" = 50, Days = NULL)
    ),
    "at most one of 'pop.pwt', 'pop.sd', and 'pop.nprior'"
  )

  expect_error(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.pwt    = 0.02,
      pop.nprior = list("(Intercept)" = 5, Days = NULL)
    ),
    "at most one of 'pop.pwt', 'pop.sd', and 'pop.nprior'"
  )

  expect_error(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.sd     = list("(Intercept)" = 50, Days = NULL),
      pop.nprior = list("(Intercept)" = 5, Days = NULL)
    ),
    "at most one of 'pop.pwt', 'pop.sd', and 'pop.nprior'"
  )

  ## Relying on the pop.pwt default (not explicitly supplied) does not
  ## conflict with an explicit pop.sd/pop.nprior.
  expect_no_error(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.sd = list("(Intercept)" = 50, Days = NULL)
    )
  )
  expect_no_error(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.nprior = list("(Intercept)" = 5, Days = NULL)
    )
  )
})

test_that("Prior_Setup_lmebayes: pop.nprior converts to the expected Sigma_fixef scaling", {
  dat <- lme4::sleepstudy
  J   <- nlevels(dat$Subject)

  n_val <- 5
  ps_nprior <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.nprior = list("(Intercept)" = n_val, Days = NULL)
  )
  ps_default <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )

  ## n_i = J*w_i/(1-w_i) <=> w_i = n_i/(n_i+J): an equivalent pop.pwt list
  ## must reproduce exactly the same Sigma_fixef for the overridden component.
  w_equiv <- n_val / (n_val + J)
  ps_pwt_equiv <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = list("(Intercept)" = w_equiv, Days = 0.01)
  )
  expect_equal(
    ps_nprior$pop.prior_list[["(Intercept)"]]$Sigma_fixef,
    ps_pwt_equiv$pop.prior_list[["(Intercept)"]]$Sigma_fixef
  )

  ## Non-overridden component (NULL entry): unchanged pop.pwt-derived scaling.
  expect_identical(
    ps_nprior$pop.prior_list[["Days"]]$Sigma_fixef,
    ps_default$pop.prior_list[["Days"]]$Sigma_fixef
  )

  ## mu_fixef is untouched by pop.nprior.
  expect_identical(
    ps_nprior$pop.prior_list[["(Intercept)"]]$mu_fixef,
    ps_default$pop.prior_list[["(Intercept)"]]$mu_fixef
  )

  ## pop.nprior output reproduces the supplied value for the overridden
  ## component, and the pop.pwt-derived default for the other.
  expect_equal(unname(ps_nprior$pop.nprior[["(Intercept)"]]), n_val)
  expect_equal(ps_nprior$pop.nprior[["Days"]], ps_default$pop.nprior[["Days"]])
})

test_that("Prior_Setup_lmebayes: pop.mu/pop.sd/pop.nprior outputs are always populated and round-trip", {
  dat <- lme4::sleepstudy
  J   <- nlevels(dat$Subject)

  ps_default <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )

  ## Always present (derived defaults), one named vector per RE component,
  ## matching pop.prior_list's mu_fixef / sqrt(diag(Sigma_fixef)).
  expect_named(ps_default$pop.mu, c("(Intercept)", "Days"))
  expect_named(ps_default$pop.sd, c("(Intercept)", "Days"))
  expect_named(ps_default$pop.nprior, c("(Intercept)", "Days"))
  for (k in c("(Intercept)", "Days")) {
    expect_identical(ps_default$pop.mu[[k]], ps_default$pop.prior_list[[k]]$mu_fixef)
    sd_k <- sqrt(diag(ps_default$pop.prior_list[[k]]$Sigma_fixef, names = FALSE))
    expect_equal(ps_default$pop.sd[[k]], sd_k, ignore_attr = TRUE)
    ## n_i = J*w_i/(1-w_i) with w_i from the low-d adaptive default (0.01
    ## for sleepstudy: both components have p_k < 14).
    expect_equal(
      ps_default$pop.nprior[[k]],
      stats::setNames(J * 0.01 / (1 - 0.01), names(ps_default$pop.mu[[k]])),
      ignore_attr = TRUE
    )
  }

  ## Round-trip: feeding the default pop.mu/pop.sd straight back in as
  ## overrides reproduces the same prior exactly (calling once without
  ## pop.mu/pop.sd to inspect the defaults, then supplying them back in).
  ps_roundtrip <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.mu = ps_default$pop.mu,
    pop.sd = ps_default$pop.sd
  )
  expect_identical(ps_roundtrip$pop.prior_list, ps_default$pop.prior_list)
  expect_identical(ps_roundtrip$pop.mu, ps_default$pop.mu)
  expect_identical(ps_roundtrip$pop.sd, ps_default$pop.sd)
  ## pop.nprior is re-derived from the pop.sd round-trip (sqrt()/division),
  ## so it matches only up to floating-point noise, not bit-for-bit.
  expect_equal(ps_roundtrip$pop.nprior, ps_default$pop.nprior)

  ## Same round-trip via pop.nprior instead of pop.sd.
  ps_roundtrip_n <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.mu     = ps_default$pop.mu,
    pop.nprior = ps_default$pop.nprior
  )
  expect_equal(ps_roundtrip_n$pop.prior_list, ps_default$pop.prior_list)
  expect_identical(ps_roundtrip_n$pop.nprior, ps_default$pop.nprior)

  ## Editing a copy of the default pop.mu before feeding it back in changes
  ## only the edited component/entry.
  mu_edit <- ps_default$pop.mu
  mu_edit[["(Intercept)"]]["(Intercept)"] <- 999
  ps_edited <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.mu = mu_edit
  )
  expect_equal(unname(ps_edited$pop.mu[["(Intercept)"]]), 999)
  expect_identical(ps_edited$pop.mu[["Days"]], ps_default$pop.mu[["Days"]])
})

## Minimal fixture: (Intercept) Block~2 model has p_k = 14 group-level
## predictors (intercept + n_cov covariates); optional random slope on
## 'x1' adds a second RE component with p_k = 1 (low-d by default).
high_d_prior_setup_fixture <- function(n_groups = 60L, n_cov = 13L, slope = "x1") {
  g <- factor(rep(seq_len(n_groups), each = 5L))
  n <- length(g)
  dat <- data.frame(g = g)
  for (j in seq_len(n_cov)) {
    nm <- paste0("x", j)
    dat[[nm]] <- rep(stats::rnorm(n_groups), each = 5L)
  }
  b0 <- stats::rnorm(1)
  b <- stats::rnorm(n_cov)
  X <- as.matrix(dat[, paste0("x", seq_len(n_cov)), drop = FALSE])
  dat$y <- as.numeric(X %*% b) + b0 +
    stats::rnorm(n_groups)[as.integer(g)] + stats::rnorm(n, sd = 0.5)
  covars <- paste(paste0("x", seq_len(n_cov)), collapse = " + ")
  re_part <- if (is.null(slope)) {
    "(1 | g)"
  } else {
    paste0("(1 + ", slope, " || g)")
  }
  form <- stats::as.formula(paste("y ~", covars, "+", re_part))
  list(data = dat, formula = form)
}

test_that("Prior_Setup_lmebayes: dimension-adaptive pop.pwt defaults per RE component", {
  dat <- lme4::sleepstudy

  ps_low <- suppressMessages(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat
    )
  )
  expect_true(all(unlist(ps_low$pop.pwt) == 0.01))

  fix <- high_d_prior_setup_fixture()
  ps_high <- suppressMessages(
    Prior_Setup_lmebayes(fix$formula, data = fix$data)
  )
  expect_equal(unique(unname(ps_high$pop.pwt[["(Intercept)"]])), 0.05)
  expect_equal(unique(unname(ps_high$pop.pwt[["x1"]])), 0.01)

  ps_explicit <- Prior_Setup_lmebayes(
    fix$formula,
    data = fix$data,
    pop.pwt = 0.01
  )
  expect_true(all(unlist(ps_explicit$pop.pwt) == 0.01))

  ps_custom <- suppressMessages(
    Prior_Setup_lmebayes(
      fix$formula,
      data = fix$data,
      pop.pwt_default_low  = 0.02,
      pop.pwt_default_high = 0.08
    )
  )
  expect_equal(unique(unname(ps_custom$pop.pwt[["(Intercept)"]])), 0.08)
  expect_equal(unique(unname(ps_custom$pop.pwt[["x1"]])), 0.02)

  for (k in names(ps_high$pop.pwt)) {
    expect_equal(
      ps_high$pop.dispersion.pwt[[k]],
      mean(ps_high$pop.pwt[[k]])
    )
  }

  expect_no_error(
    Prior_Setup_lmebayes(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.sd = list("(Intercept)" = 50, Days = NULL)
    )
  )
})
