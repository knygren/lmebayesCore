## Precompute the Chapter-C01 artifact: the theoretical model, the dataset
## description, and the estimated model shared by every Gaussian chapter in
## the C series. No sampling happens here -- C01 describes the model, C02
## onwards draw from it.
##
## Run from package root:
##   Rscript data-raw/make_Chapter-C01.R
##
## Writes: inst/extdata/Chapter-C01-bwc-setup.rds

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

suppressMessages(library(lmebayesCore))
source(file.path(root, "data-raw", "chapter_C_gaussian_setup.R"))

message("Chapter-C01: model, data and design")

t_setup <- system.time(S <- chapter_C_setup())[["elapsed"]]
dat    <- S$data
design <- S$design
ps     <- S$ps
pf     <- S$pf
dims   <- S$dims

re_names <- design$groupef.names

## --- dataset description -------------------------------------------------
## Variable definitions are quoted from bayesrules::big_word_club.
var_desc <- rbind(
  c("school_id", "school identifier (the grouping factor)"),
  c("score_ppvt", "student score on the standardized vocabulary assessment"),
  c("private_school", "whether the school is private"),
  c("title1", "whether the school has Title 1 status"),
  c("free_reduced_lunch", "percent of school receiving free / reduced lunch"),
  c("distracted_a1", "student distraction during assessment 1 (0-3)"),
  c("distracted_ppvt", "student distraction during the standardized assessment (0-3)")
)
model_vars <- data.frame(
  variable = var_desc[, 1],
  meaning  = var_desc[, 2],
  stringsAsFactors = FALSE
)
model_vars$level <- ifelse(
  model_vars$variable == "school_id", "grouping",
  ifelse(
    model_vars$variable %in% c("private_school", "title1",
                               "free_reduced_lunch"),
    "school", "pupil"
  )
)
model_vars$values <- vapply(model_vars$variable, function(nm) {
  v <- dat[[nm]]
  if (is.factor(v)) {
    sprintf("%d levels", nlevels(v))
  } else if (length(unique(v)) == 2L) {
    sprintf("%.0f%% = 1", 100 * mean(v))
  } else {
    sprintf("%g to %g (median %g)", min(v), max(v), stats::median(v))
  }
}, character(1))
rownames(model_vars) <- NULL

## --- Block 2 prior -------------------------------------------------------
pop_prior <- do.call(rbind, lapply(re_names, function(k) {
  pl <- ps$pop.prior_list[[k]]
  data.frame(
    component  = k,
    parameter  = rownames(pl$Sigma),
    prior_mean = as.numeric(pl$mu),
    prior_sd   = sqrt(diag(as.matrix(pl$Sigma))),
    stringsAsFactors = FALSE, row.names = NULL
  )
}))

## --- what is held fixed, against the REML quantities ---------------------
fixed_vs_reml <- data.frame(
  component = c(re_names, "Residual"),
  calibrated = c(
    vapply(re_names, function(k) as.numeric(pf[[k]]$prior_list$dispersion),
           numeric(1)),
    as.numeric(ps$group.dispersion)
  ),
  reml = c(
    vapply(re_names, function(k) unname(design$Psi[[k]]), numeric(1)),
    design$dispersion
  ),
  stringsAsFactors = FALSE, row.names = NULL
)

out <- list(
  meta = list(
    dataset       = "bayesrules::big_word_club",
    formula       = dims$formula,
    n_obs         = dims$n_obs,
    J             = dims$J,
    groupef_names = dims$groupef.names,
    group_name    = dims$group_name,
    p_re          = length(dims$groupef.names),
    n_cols        = ncol(dat),
    secs_setup    = as.numeric(t_setup),
    built         = as.character(Sys.Date()),
    R_version     = R.version.string
  ),
  provenance  = S$provenance,
  all_columns = names(dat),
  model_vars  = model_vars,
  sizes       = as.integer(dims$sizes_kept),
  popef_names = lapply(pf, function(p) rownames(p$prior_list$Sigma)),
  pop_prior   = pop_prior,
  fixed_vs_reml = fixed_vs_reml,
  dispersion_source    = attr(ps$group.dispersion, "source"),
  dispersion_classical = attr(ps$group.dispersion, "classical")
)

dir.create(file.path(root, "inst", "extdata"),
           recursive = TRUE, showWarnings = FALSE)
dest <- file.path(root, "inst", "extdata", "Chapter-C01-bwc-setup.rds")
saveRDS(out, dest, compress = "xz")
message("Wrote ", normalizePath(dest, winslash = "/"),
        " (", round(file.size(dest) / 1024, 1), " KB)")
