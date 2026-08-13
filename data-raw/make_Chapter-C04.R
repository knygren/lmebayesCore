## Precompute the Chapter-C04 artifact: the two-block Gibbs sampler on the
## same design, priors and fixed variance components as Chapter-C02 (see
## data-raw/chapter_C_gaussian_setup.R).
##
## Identical target to Chapter-C02 -- known variance components (all Block 2
## dNormal) and a fixed pooled residual variance -- but forced onto the
## iterative route with sim_method = "TWO_BLOCK_GIBBS", so the run produces a
## sweep history that the exact route has nothing to report.
##
## Run from package root:
##   Rscript data-raw/make_Chapter-C04.R
##
## The draw count is the run time here. Unlike Chapter-C02's exact route
## (~9s at n = 10000), every draw is an independent chain carried through
## m_convergence sweeps, so cost scales with n * sweeps. Override with:
##
##   $env:CHAPTER_C_N = "1000";  Rscript data-raw/make_Chapter-C04.R   # PowerShell
##   CHAPTER_C_N=1000 Rscript data-raw/make_Chapter-C04.R              # bash
##
## Start small (1000) to time a run before committing to the full 10000.
## The progress bar tracks chains within the current sweep.
##
## Writes: inst/extdata/Chapter-C04-bwc-twoblock.rds

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

## Matches Chapter-C02 so the two chapters' tables are directly comparable.
N_DRAWS <- as.integer(Sys.getenv("CHAPTER_C_N", "10000"))
SEED    <- 1L

suppressMessages(library(lmebayesCore))
source(file.path(root, "data-raw", "chapter_C_gaussian_setup.R"))

message("Chapter-C04: n = ", N_DRAWS)

t_setup <- system.time(S <- chapter_C_setup())[["elapsed"]]
design <- S$design
ps     <- S$ps
pf     <- S$pf
dims   <- S$dims

sizes <- dims$sizes_kept
show_groups <- names(sizes)[c(which.max(sizes), which.min(sizes))]

set.seed(SEED)
t_fit <- system.time({
  fit <- rlmerb(
    n              = N_DRAWS,
    design         = design,
    pfamily_list   = pf,
    dispprior_list = list(dispersion = ps$group.dispersion),
    sim_method     = "TWO_BLOCK_GIBBS",
    ## Gibbs runs many sweeps per chain and takes minutes at this n, unlike
    ## the exact route in Chapter-C02; leave the progress bar on so the build
    ## is watchable.
    progbar = TRUE, verbose = FALSE, print_icm_table = FALSE
  )
})[["elapsed"]]

message("  setup ", round(t_setup, 1), "s; sampler ", round(t_fit, 1), "s")
message("  sweeps to convergence: ", fit$m_convergence)

smry <- summary(fit, groups = show_groups)

## Stripped design + measurement_prior_list let plot_var_convergence() resolve
## the exact posterior covariance (lmerb_posterior_covariance) as Var_ref.
## Drop lmer and the data frame; keep only what the Schur system needs.
##
## IMPORTANT: use the same Block~1 precision the sampler uses
## (.rLMM_P_from_pfamily_list(pf)), NOT ps$group.Sigma -- that object carries
## the REML reference variances, while dNormal pfamily_list fixes tau2_k at the
## calibrated plug-ins (see tau2_prior_overview E[tau2] vs lmer).
design_plot <- design[c("W", "D", "y", "group", "groupef.names")]
dispprior_ref <- list(dispersion = ps$group.dispersion)
dispprior_ref$P <- lmebayesCore:::.rLMM_P_from_pfamily_list(pf, design$groupef.names)
measurement_prior_list <- lmebayesCore:::.two_block_measurement_prior_list(
  prior_list_block1 = dispprior_ref,
  pfamily_list      = pf,
  re_names          = design$groupef.names,
  x_hyper           = design$W,
  family            = gaussian()
)

## The Nygren (2020) rate object behind m_convergence. Cheap to recompute
## (one generalized eigenproblem), but storing it lets the chapter show the
## Theorem 3 / Corollary 1 sweep budget without rebuilding the design.
rate <- two_block_rate_from_pfamily_list(
  x                 = design$D,
  group             = design$group,
  x_hyper           = design$W,
  prior_list_block1 = dispprior_ref,
  pfamily_list      = pf,
  family            = gaussian(),
  group_levels      = levels(design$group)
)

out <- list(
  meta = list(
    dataset     = "bayesrules::big_word_club",
    n_draws     = N_DRAWS,
    seed        = SEED,
    n_obs       = dims$n_obs,
    J           = dims$J,
    groupef_names = dims$groupef.names,
    show_groups = show_groups,
    ## plot_var_convergence() needs the chain count, which is not carried on
    ## the history object itself.
    n_chains    = nrow(fit$popef[[1L]]),
    secs_setup  = as.numeric(t_setup),
    secs_fit    = as.numeric(t_fit),
    built       = as.character(Sys.Date()),
    R_version   = R.version.string
  ),
  summary = smry,
  ## Stored whole and plotted as-is by the chapter: plot_var_convergence()
  ## takes this object directly.
  sweep_history = fit$sweep_history,
  ## Exact Var_ref for plot_var_convergence(..., whitened = FALSE).
  design = design_plot,
  measurement_prior_list = measurement_prior_list,
  rate = rate,
  convergence = list(
    m_convergence    = fit$m_convergence,
    lambda_star      = fit$convergence_info$lambda_star,
    convergence_info = fit$convergence_info
  )
)

dir.create(file.path(root, "inst", "extdata"),
           recursive = TRUE, showWarnings = FALSE)
dest <- file.path(root, "inst", "extdata", "Chapter-C04-bwc-twoblock.rds")
saveRDS(out, dest, compress = "xz")
message("Wrote ", normalizePath(dest, winslash = "/"),
        " (", round(file.size(dest) / 1024, 1), " KB)")
