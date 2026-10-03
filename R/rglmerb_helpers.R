## Internal helpers used only by rglmerb() (glmerb only).

#' Shared matrix-level arguments for GLMM reg routes
#' @noRd
.lmebayes_matrix_args_glmm <- function(
    n,
    design,
    prior,
    family,
    gap_tol       = 0.0196,
    tv_tol        = 0.01,
    mode_gap_max  = 1.0,
    verbose       = FALSE,
    progbar       = FALSE,
    collect_block1 = TRUE,
    offset        = NULL,
    weights       = 1,
    offset_missing = FALSE,
    weights_missing = FALSE
) {
  block1_prior <- .lmebayes_block1_prior_list(prior, group.dispersion = NULL)

  ## The routed export has no 'group_name' formal: attach it to 'group'
  ## itself (design$group is never a bare variable here, so the export's
  ## substitute()-based fallback could not resolve it anyway).
  grp <- design$group
  attr(grp, "group_name") <- design$group_name
  ow <- .lmebayes_resolve_offset_weights(
    offset, weights, design,
    offset_missing = offset_missing,
    weights_missing = weights_missing
  )

  list(
    n               = n,
    y               = design$y,
    D               = design$D,
    group           = grp,
    W               = design$W,
    pfamily_list    = prior$pfamily_list,
    dispprior_list  = block1_prior,
    offset          = ow$offset,
    weights         = ow$weights,
    family          = family,
    gap_tol         = gap_tol,
    tv_tol          = tv_tol,
    mode_gap_max    = mode_gap_max,
    verbose         = verbose,
    progbar         = progbar,
    stage_verbose   = verbose,
    collect_block1  = collect_block1
  )
}

#' @noRd
.lmebayes_run_glmm_engine <- function(
    n,
    design,
    prior,
    family,
    gap_tol       = 0.0196,
    tv_tol        = 0.01,
    mode_gap_max  = 1.0,
    verbose       = FALSE,
    progbar       = FALSE,
    collect_block1 = TRUE,
    offset        = NULL,
    weights       = 1,
    offset_missing = FALSE,
    weights_missing = FALSE
) {
  route_key <- .lmebayes_reg_route_key(
    family         = family,
    disp_mode      = "none",
    any_non_normal = prior$any_non_normal
  )
  route <- .lmebayes_reg_route_fn(route_key)
  args  <- .lmebayes_matrix_args_glmm(
    n              = n,
    design         = design,
    prior          = prior,
    family         = family,
    gap_tol        = gap_tol,
    tv_tol         = tv_tol,
    mode_gap_max   = mode_gap_max,
    verbose        = verbose,
    progbar        = progbar,
    collect_block1 = collect_block1,
    offset         = offset,
    weights        = weights,
    offset_missing = offset_missing,
    weights_missing = weights_missing
  )
  out <- do.call(route$export_fn, args)
  disp_none <- list(mode = "none")
  .lmebayes_attach_group_dispersion(out, disp_none)
}

#' @noRd
.lmebayes_print_ranef_mode_reference <- function(
    ranef_mode,
    re_names,
    group_levels,
    verbose
) {
  invisible(NULL)
}
