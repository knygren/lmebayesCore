## Shared data, design and prior setup for the Gaussian C-series chapters.
##
## Sourced by data-raw/make_Chapter-C0*.R so that every Gaussian chapter is
## demonstrably fitted to the same schools with the same priors. The model
## and the full-rank filter follow
## demo("Ex_13b_rLMM_estimated_dispersion_known_vcov_BigWordClub_PartVI").

chapter_C_formula <- function() {
  score_ppvt ~
    private_school + title1 + free_reduced_lunch +
    distracted_ppvt + distracted_a1 +
    free_reduced_lunch:distracted_a1 +
    (1 + distracted_ppvt + distracted_a1 || school_id)
}

## Returns the analysis data, the design, the calibrated prior setup, the
## Block 2 pfamily list, and the dimensions the chapters quote. The dataset
## is taken as given from data-raw/make_bwc_full_rank.R.
chapter_C_setup <- function(pop.pwt = 0.01,
                            data_path = file.path("inst", "extdata",
                                                  "bwc_full_rank.rds")) {

  if (!file.exists(data_path)) {
    stop("Analysis dataset not found at '", data_path,
         "'. Run data-raw/make_bwc_full_rank.R first.", call. = FALSE)
  }
  dat  <- readRDS(data_path)
  prov <- attr(dat, "provenance")

  form <- chapter_C_formula()

  design <- model_setup(form, data = dat)
  ps     <- Prior_Setup_GLMM(form, data = dat, pop.pwt = pop.pwt)
  pf     <- pfamily_list(ps)

  list(
    data   = dat,
    form   = form,
    design = design,
    ps     = ps,
    pf     = pf,
    provenance = prov,
    dims   = list(
      n_obs          = nrow(dat),
      J              = nlevels(dat$school_id),
      sizes_kept     = table(dat$school_id),
      groupef.names  = design$groupef.names,
      group_name     = design$group_name,
      formula        = paste(deparse(form), collapse = " ")
    )
  )
}
