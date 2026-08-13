## Build the pre-computed Big Word Club analysis dataset used by the Gaussian
## C-series chapters.
##
## Two restrictions are applied here, once, so that the chapters can take the
## dataset as given:
##   1. complete records with a valid standardized-assessment score;
##   2. schools that are full column rank for the three school-level
##      coefficients the model estimates (the per-group routes in the later
##      chapters have no fallback for a rank-deficient group).
##
## Run from package root:
##   Rscript data-raw/make_bwc_full_rank.R
##
## Writes: inst/extdata/bwc_full_rank.rds

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

suppressMessages(library(lmebayesCore))
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

data(big_word_club, package = "bayesrules")

model_vars <- c("score_ppvt", "distracted_a1", "distracted_ppvt",
                "private_school", "title1", "free_reduced_lunch",
                "school_id")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
n_raw <- nrow(dat)
J_raw <- nlevels(dat$school_id)

dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    stats::complete.cases(dat[, model_vars])
)
n_complete <- nrow(dat)

form <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design_all <- model_setup(form, data = dat)
full_rank  <- names(design_all$groupef.rank)[design_all$groupef.rank]

dat <- subset(dat, school_id %in% full_rank)
dat$school_id <- droplevels(dat$school_id)
rownames(dat) <- NULL

attr(dat, "provenance") <- list(
  source        = "bayesrules::big_word_club",
  n_raw         = n_raw,
  J_raw         = J_raw,
  n_complete    = n_complete,
  n_final       = nrow(dat),
  J_final       = nlevels(dat$school_id),
  model_vars    = model_vars,
  formula       = paste(deparse(form), collapse = " "),
  built         = as.character(Sys.Date())
)

dir.create(file.path(root, "inst", "extdata"),
           recursive = TRUE, showWarnings = FALSE)
dest <- file.path(root, "inst", "extdata", "bwc_full_rank.rds")
saveRDS(dat, dest, compress = "xz")

message("Big Word Club analysis dataset")
message("  raw:        ", n_raw, " pupils, ", J_raw, " schools")
message("  complete:   ", n_complete, " pupils")
message("  final:      ", nrow(dat), " pupils, ",
        nlevels(dat$school_id), " schools")
message("Wrote ", normalizePath(dest, winslash = "/"),
        " (", round(file.size(dest) / 1024, 1), " KB)")
