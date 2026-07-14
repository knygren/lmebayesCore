test_that("Prior_Setup_lmebayes: dispformula = ~1 keeps the pooled lme4 reference", {
  dat <- lme4::sleepstudy

  ps <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pwt  = 0.01
  )

  expect_identical(ps$calibration_source, "lme4")
  expect_s4_class(ps$fit_ref, "merMod")
  expect_identical(ps$fit_ref, ps$mer_fit)
  expect_identical(ps$fit_ref, ps$design$lmer_fit)
  expect_null(ps$dispersion_fit)
  expect_null(ps$sigma2_group)
  expect_null(ps$ing_prior_measurement_group)
})

test_that("Prior_Setup_lmebayes: dispformula = ~group routes calibration through glmmTMB", {
  skip_if_not_installed("glmmTMB")
  dat <- lme4::sleepstudy

  ps <- Prior_Setup_lmebayes(
    Reaction ~ Days + (Days || Subject),
    data            = dat,
    pwt             = 0.01,
    pwt_measurement = 0.1,
    dispformula     = ~Subject
  )

  expect_identical(ps$calibration_source, "glmmTMB")
  expect_s3_class(ps$fit_ref, "glmmTMB")
  expect_s4_class(ps$mer_fit, "merMod")
  expect_identical(ps$dispersion_fit, ps$fit_ref)
  expect_identical(ps$mer_fit, ps$design$lmer_fit)

  ## Pooled dispersion_ranef stays lme4-derived (design$residual_var from
  ## mer_fit) regardless of dispformula.
  expect_identical(ps$dispersion_ranef, ps$design$residual_var)

  expect_false(is.null(ps$ing_prior_measurement_group))
  expect_length(ps$ing_prior_measurement_group, nlevels(dat$Subject))

  expect_type(ps$sigma2_group, "double")
  expect_named(ps$sigma2_group, levels(dat$Subject))

  disp_pf <- dGamma_list(ps, warn_asymmetric = FALSE)
  expect_length(disp_pf, nlevels(dat$Subject))
  expect_true(all(vapply(disp_pf, inherits, logical(1L), "pfamily")))
  expect_identical(attr(disp_pf, "dispersion_fit"), ps$dispersion_fit)
  expect_identical(attr(disp_pf, "calibration_source"), "glmmTMB")
})
