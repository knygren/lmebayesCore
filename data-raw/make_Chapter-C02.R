## Precompute the Chapter-C02 artifact: the exact Gaussian iid sampler on the
## shared C-series design (see data-raw/chapter_C_gaussian_setup.R).
##
## Known variance components (all Block 2 dNormal) and a fixed pooled
## residual variance, so rlmerb(sim_method = "DEFAULT") draws directly from
## the closed-form joint Gaussian posterior -- no sweeps, no burn-in.
##
## Run from package root:
##   Rscript data-raw/make_Chapter-C02.R
##
## Writes: inst/extdata/Chapter-C02-bwc-iid.rds

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

## Draws. 10000 keeps the Monte Carlo error small enough that the shrinkage
## toward the prior means is visible on the weakly identified coefficients.
N_DRAWS <- as.integer(Sys.getenv("CHAPTER_C_N", "10000"))
SEED    <- 1L

suppressMessages(library(lmebayesCore))
source(file.path(root, "data-raw", "chapter_C_gaussian_setup.R"))

message("Chapter-C02: n = ", N_DRAWS)

t_setup <- system.time(S <- chapter_C_setup())[["elapsed"]]
design <- S$design
ps     <- S$ps
pf     <- S$pf
dims   <- S$dims

## Two schools for the per-group detail table: the largest and the smallest,
## so the chapter can show how much the group size moves the posterior.
sizes <- dims$sizes_kept
show_groups <- names(sizes)[c(which.max(sizes), which.min(sizes))]

set.seed(SEED)
t_fit <- system.time({
  fit <- rlmerb(
    n              = N_DRAWS,
    design         = design,
    pfamily_list   = pf,
    dispprior_list = list(dispersion = ps$group.dispersion),
    sim_method     = "DEFAULT",
    progbar = FALSE, verbose = FALSE, print_icm_table = FALSE
  )
})[["elapsed"]]

message("  setup ", round(t_setup, 1), "s; sampler ", round(t_fit, 1), "s")

## The summary holds tables rather than draws, so it stores compactly and
## does not grow with n. The fit itself does: fit$pfamily_list carries ~2.4 MB
## of closures and fit$groupef is J x p_re x n rows.
smry <- summary(fit, groups = show_groups)

re_names <- design$groupef.names
gd <- fit$groupef
group_summary <- do.call(rbind, lapply(re_names, function(k) {
  agg_m <- tapply(gd[[k]], gd[[dims$group_name]], mean)
  agg_s <- tapply(gd[[k]], gd[[dims$group_name]], stats::sd)
  data.frame(
    group     = names(agg_m),
    component = k,
    mean      = as.numeric(agg_m),
    sd        = as.numeric(agg_s),
    n_obs     = as.integer(sizes[names(agg_m)]),
    stringsAsFactors = FALSE, row.names = NULL
  )
}))

## Reference REML fit for the appendix. Prior_Setup_GLMM() returns it as
## mer_fit (the same object model_setup() stores in design$lmer). Stored as
## printed text: the merMod itself carries the model frame and would dwarf
## everything else in this artifact.
##
## model_setup() builds the fit with do.call(), which inlines the data frame
## and the lmerControl object as literals into the call, so print() would dump
## all 400 rows. Restore them as symbols to get an ordinary lmer summary.
mer <- ps$mer_fit
mer@call$data    <- quote(bwc)
mer@call$control <- quote(lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5)))
lmer_summary <- capture.output(print(summary(mer)))

out <- list(
  meta = list(
    dataset     = "bayesrules::big_word_club",
    n_draws     = N_DRAWS,
    seed        = SEED,
    n_obs       = dims$n_obs,
    J           = dims$J,
    groupef_names = dims$groupef.names,
    show_groups = show_groups,
    secs_setup  = as.numeric(t_setup),
    secs_fit    = as.numeric(t_fit),
    built       = as.character(Sys.Date()),
    R_version   = R.version.string
  ),
  summary       = smry,
  lmer_summary  = lmer_summary,
  group_summary = group_summary,
  convergence = list(
    m_convergence    = fit$m_convergence,
    convergence_info = fit$convergence_info,
    sweep_history    = fit$sweep_history
  )
)

dir.create(file.path(root, "inst", "extdata"),
           recursive = TRUE, showWarnings = FALSE)
dest <- file.path(root, "inst", "extdata", "Chapter-C02-bwc-iid.rds")
saveRDS(out, dest, compress = "xz")
message("Wrote ", normalizePath(dest, winslash = "/"),
        " (", round(file.size(dest) / 1024, 1), " KB)")
