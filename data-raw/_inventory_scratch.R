# Generate inst/R_*.md inventory (run from glmbayesCore root)
root <- normalizePath("..")
core_dir <- file.path(root, "glmbayesCore")
r_dir <- file.path(core_dir, "R")

parse_exports <- function(ns_file) {
  lines <- readLines(ns_file, warn = FALSE)
  exp <- grep("^export\\(", lines, value = TRUE)
  sort(gsub("^export\\((.+)\\)$", "\\1", exp))
}

parse_importfrom <- function(ns_file, pkg = "glmbayesCore") {
  lines <- readLines(ns_file, warn = FALSE)
  pat <- paste0("^importFrom\\(", pkg, ",")
  hits <- grep(pat, lines, value = TRUE)
  sort(gsub(paste0("^importFrom\\(", pkg, ",(.+)\\)$"), "\\1", hits))
}

core <- parse_exports(file.path(core_dir, "NAMESPACE"))
glmb <- parse_exports(file.path(root, "glmbayes", "NAMESPACE"))
lme  <- parse_exports(file.path(root, "lmebayes", "NAMESPACE"))
lme_imp <- parse_importfrom(file.path(root, "lmebayes", "NAMESPACE"))
lme_reexport <- c(
  "Prior_Setup", "dNormal", "dNormal_Gamma", "dIndependent_Normal_Gamma",
  "dGamma", "pfamily_list"
)
lme_present <- sort(unique(c(intersect(core, lme), lme_reexport, lme_imp)))

scan_calls <- function(pkg, prefix = "glmbayesCore::") {
  r_files <- list.files(file.path(root, pkg, "R"), pattern = "\\.R$", full.names = TRUE)
  if (!length(r_files)) return(character())
  code <- paste(unlist(lapply(r_files, readLines, warn = FALSE)), collapse = "\n")
  hits <- grep(prefix, strsplit(code, "\n")[[1]], value = TRUE)
  sort(unique(sub(paste0(".*", prefix, "([.:A-Za-z0-9_]+).*"), "\\1", hits)))
}

lme_calls <- scan_calls("lmebayes")
lme_r_files <- list.files(file.path(root, "lmebayes", "R"), pattern = "\\.R$", full.names = TRUE)
lme_all_lines <- unlist(lapply(lme_r_files, readLines, warn = FALSE))
lme_get <- grep("getFromNamespace\\(.+glmbayesCore", lme_all_lines, value = TRUE)
lme_int <- sort(unique(c(
  grep("^\\.", scan_calls("lmebayes", "glmbayesCore:::"), value = TRUE),
  gsub('.*getFromNamespace\\(\\s*"([^"]+)"\\s*,\\s*"glmbayesCore".*', "\\1",
       lme_get[nchar(gsub('.*getFromNamespace\\(\\s*"([^"]+)"\\s*,\\s*"glmbayesCore".*', "\\1", lme_get)) > 0])
)))

cat_a <- intersect(core, glmb)
cat_b <- setdiff(intersect(core, lme_present), cat_a)
cat_c <- setdiff(intersect(lme_calls, core), Reduce(union, list(cat_a, cat_b)))
cat_d <- setdiff(core, Reduce(union, list(cat_a, cat_b, cat_c)))

r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
r_basenames <- basename(r_files)

find_def_file <- function(fn) {
  esc <- gsub("\\.", "\\\\.", fn)
  pat <- paste0("^", esc, "\\s*<-\\s*function")
  for (f in r_files) {
    if (any(grepl(pat, readLines(f, warn = FALSE)))) return(basename(f))
  }
  NA_character_
}

extract_title <- function(fn, def_file) {
  if (is.na(def_file)) return("—")
  path <- file.path(r_dir, def_file)
  lines <- readLines(path, warn = FALSE)
  esc <- gsub("\\.", "\\\\.", fn)
  idx <- grep(paste0("^", esc, "\\s*<-"), lines)[1]
  if (is.na(idx)) return("—")
  i <- idx - 1
  while (i >= 1 && grepl("^#'", lines[i])) {
    if (grepl("^#' @title ", lines[i])) {
      return(sub("^#' @title ", "", lines[i]))
    }
    if (!grepl("^#' @", lines[i])) {
      line <- trimws(sub("^#'\\s*", "", lines[i]))
      if (nzchar(line)) return(line)
    }
    i <- i - 1
  }
  "—"
}

find_callers <- function(sym) {
  esc <- gsub("\\.", "\\\\.", sym)
  pat <- paste0("(^|[^[:alnum:]._])", esc, "\\s*\\(")
  def <- find_def_file(sym)
  out <- character()
  for (f in r_files) {
    bn <- basename(f)
    if (!is.na(def) && bn == def) next
    if (any(grepl(pat, readLines(f, warn = FALSE)))) out <- c(out, bn)
  }
  sort(unique(out))
}

# @noRd / internal helpers (not in NAMESPACE)
is_exported <- function(fn) fn %in% core

no_rd <- list()
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  for (i in seq_along(lines)) {
    if (!grepl("@noRd", lines[i])) next
    j <- i
    while (j <= length(lines) && !grepl("<-\\s*function", lines[j])) j <- j + 1
    if (j > length(lines)) next
    fn <- sub("\\s*<-.*", "", trimws(lines[j]))
    if (is_exported(fn)) next
    role <- ""
    k <- i
    while (k >= 1 && grepl("^#'", lines[k])) {
      if (grepl("^#'\\s*@", lines[k])) break
      if (nzchar(role <- sub("^#'\\s*", "", lines[k]))) break
      k <- k - 1
    }
    no_rd[[fn]] <- list(file = basename(f), role = role)
  }
}

# also @keywords internal without export
for (f in r_files) {
  lines <- readLines(f, warn = FALSE)
  for (i in seq_along(lines)) {
    if (!grepl("@keywords internal", lines[i])) next
    j <- i
    while (j <= length(lines) && !grepl("<-\\s*function", lines[j])) j <- j + 1
    if (j > length(lines)) next
    fn <- sub("\\s*<-.*", "", trimws(lines[j]))
    if (is_exported(fn) || fn %in% names(no_rd)) next
    role <- ""
    k <- i
    while (k >= 1 && grepl("^#'", lines[k])) {
      if (grepl("^#'\\s*@", lines[k])) break
      if (nzchar(role <- sub("^#'\\s*", "", lines[k]))) break
      k <- k - 1
    }
    no_rd[[fn]] <- list(file = basename(f), role = role)
  }
}

md_row <- function(fn) {
  ff <- find_def_file(fn)
  if (is.na(ff)) ff <- "—"
  role <- extract_title(fn, if (ff == "—") NA_character_ else ff)
  sprintf("| `%s` | `%s` | %s |", fn, ff, role)
}

write_section <- function(title, desc, fns) {
  c(
    sprintf("## %s", title),
    "",
    desc,
    "",
    "| Function | File | Role |",
    "|----------|------|------|",
    if (length(fns)) unlist(lapply(fns, md_row)) else "| *(none)* | — | — |",
    "",
    "---",
    ""
  )
}

out_exp <- c(
  "# `R/` — exported and documented functions",
  "",
  "Symbols defined in **`R/`** that are **exported** (`NAMESPACE`) or have a **help",
  "page** (`man/*.Rd`). Use this list when reviewing the public API, `\\usage`",
  "blocks, and README coverage.",
  "",
  "Companion: [R_INTERNAL_HELPERS.md](R_INTERNAL_HELPERS.md) (`@noRd` and other",
  "undocumented helpers).",
  "",
  "Exports are grouped by overlap with **glmbayes** and **lmebayes** (non-overlapping",
  "priority: glmbayes → lmebayes → lmebayes-only callers → glmbayesCore-only).",
  "",
  "---",
  ""
)

out_exp <- c(out_exp, write_section(
  "Also exported from **glmbayes** (shared API)",
  paste0(
    "Present in **`NAMESPACE`** of both **glmbayesCore** and **glmbayes** (",
    length(cat_a), " symbols). Signatures should stay aligned per package policy."
  ),
  cat_a
))

out_exp <- c(out_exp, write_section(
  "Present in **lmebayes** (re-export or `importFrom`)",
  paste0(
    "Exported or re-exported from **lmebayes** (`NAMESPACE` export, ",
    "`reexports_glmbayesCore.R`, or `importFrom(glmbayesCore, …)`), but **not** ",
    "in **glmbayes** (", length(cat_b), " symbols)."
  ),
  cat_b
))

out_exp <- c(out_exp, write_section(
  "Called from **lmebayes** / **glmbayes** without being exported there",
  paste0(
    "Referenced from **lmebayes** `R/` as `glmbayesCore::…` (or ",
    "`getFromNamespace`), exported from **glmbayesCore** only (", length(cat_c),
    " symbols). **glmbayes** does not depend on **glmbayesCore**."
  ),
  cat_c
))

out_exp <- c(out_exp, write_section(
  "**glmbayesCore**-only exports",
  paste0(
    "Not in **glmbayes** or **lmebayes** export surfaces and not called directly ",
    "from **lmebayes** `R/` (", length(cat_d), " symbols): two-block engines, ",
    "multi-response samplers, block updates, OpenCL probe, etc."
  ),
  cat_d
))

# S3 methods in NAMESPACE (not plain exports)
ns_lines <- readLines(file.path(core_dir, "NAMESPACE"), warn = FALSE)
s3_lines <- grep("^S3method\\(", ns_lines, value = TRUE)
s3_pairs <- lapply(s3_lines, function(l) {
  m <- regmatches(l, regexec("^S3method\\(([^,]+),(.+)\\)$", l))[[1]]
  c(generic = m[2], class = m[3])
})

write_s3_section <- function() {
  rows <- vapply(s3_pairs, function(p) {
    fn <- paste0(p["generic"], ".", p["class"])
    ff <- find_def_file(fn)
    if (is.na(ff)) ff <- "—"
    role <- extract_title(fn, if (ff == "—") NA_character_ else ff)
    sprintf("| `%s()` | `%s` | %s |", fn, ff, role)
  }, character(1))
  c(
    "## S3 methods (`NAMESPACE` → `S3method`)",
    "",
    "Registered methods with help pages; not listed in the export groups above.",
    "",
    "| Method | File | Role |",
    "|--------|------|------|",
    rows,
    "",
    "---",
    ""
  )
}

out_exp <- c(out_exp, write_s3_section(),
  "## Documentation topics (no function body in `R/`)",
  "",
  "| Topic / file | Contents |",
  "|--------------|----------|",
  "| `glmbayesCore-package.R` | Package meta, imports (`\"_PACKAGE\"`). |",
  "| `gpu_diagnostics.R` | `diagnose_glmbayes()`, `glmbayesCore_has_opencl()`. |",
  "| `data-*.R` | Lazy data docs (`Boston_centered`, `BikeSharing`, `carinsca`, etc.). |",
  "",
  "---",
  "",
  "## Review checklist (exports / docs)",
  "",
  "| Priority | Item |",
  "|----------|------|",
  "| 1 | Keep **glmbayes**-shared exports signature-aligned when touching `R/`. |",
  "| 2 | Document **lmebayes** dependency surface (`cat_c`) in `inst/ARCHITECTURE_glmerb.md`. |",
  "| 3 | Run `devtools::document()` after any `@export` or `\\usage` change. |",
  ""
)

writeLines(out_exp, file.path(core_dir, "inst", "R_EXPORTED_AND_DOCUMENTED.md.generated"))
message("Note: curated inst/R_EXPORTED_AND_DOCUMENTED.md is maintained by hand; ",
        "regenerate from .generated if needed.")

helper_row <- function(fn) {
  info <- no_rd[[fn]]
  cl <- find_callers(fn)
  cf <- if (length(cl)) paste(cl, collapse = ", ") else "*(unused)*"
  role <- info$role
  if (!nzchar(role) || grepl("^@", role) || grepl("[\\\\{}]", role)) {
    role <- extract_title(fn, info$file)
  }
  if (!nzchar(role) || grepl("[\\\\{}]", role)) role <- "—"
  sprintf("| `%s` | `%s` | %s | %s |", fn, info$file, role, cf)
}

write_helper_section <- function(title, fns) {
  fns <- intersect(fns, names(no_rd))
  if (!length(fns)) return(character())
  c(
    sprintf("## %s", title),
    "",
    "| Function | File | Role | Called from |",
    "|----------|------|------|-------------|",
    unlist(lapply(fns, helper_row)),
    "",
    "---",
    ""
  )
}

# Helpers grouped by theme (exclusive — specific groups before broad patterns)
helper_groups <- list(
  "Two-block rate / TV (`two_block_ergodicity.R`)" = grep("^\\.two_block_(rate|S_P|gen_eigen|erfn|tv_bound)|^two_block_mode_weights$", names(no_rd), value = TRUE),
  "Two-block pilot / GLMM (`two_block_glmm_pilot_helpers.R`, `two_block_pilot_cost.R`)" = grep("^\\.two_block_(print_pilot|as_staged|pilot|fixef_col|resolve_n|resolve_pilot|pilot_will)", names(no_rd), value = TRUE),
  "Two-block sweep history (`two_block_sweep_history.R`)" = grep("^\\.two_block_(build_sweep|filter_sweep|sweep_history|print_sweep)", names(no_rd), value = TRUE),
  "Two-block measurement / tau2 (`two_block_measurement_prior.R`, `two_block_tau2_ref.R`)" = grep("^\\.two_block_(measurement|icm_at|validate_gap|tau2)", names(no_rd), value = TRUE),
  "Two-block drivers (`two_block_rNormal_reg*.R`)" = grep("^\\.two_block_(normalize|validate|block1|mu_all|format_v)", names(no_rd), value = TRUE),
  "Two-block batch Gibbs (`two_block_batch_gibbs.R`)" = names(no_rd)[vapply(names(no_rd), function(x) identical(no_rd[[x]]$file, "two_block_batch_gibbs.R"), logical(1))],
  "LMM engines (`rLMMNormal_reg.R`, `two_block_lmm_staged_sweep_outer.R`)" = grep("^\\.rLMM", names(no_rd), value = TRUE),
  "GLMM sweep (`rGLMM.R`, `two_block_batch_gibbs.R`)" = grep("^\\.rGLMM", names(no_rd), value = TRUE),
  "Multi-response / pfamily validation (`multi_rlmb.R`, `multi_rNormal_reg.R`)" = grep("^\\.mrglmb|^\\.validate_|^\\.check_symmetric", names(no_rd), value = TRUE),
  "Block / simfunction utils (`simfunction_block_utils.R`)" = grep("^\\.prior_|^normalize_prior_for_blocks|^\\.check_P", names(no_rd), value = TRUE),
  "lmerb / build_mu (`build_mu_all.R`, `lmerb_posterior_mean.R`)" = grep("^\\.lmerb|^build_mu_all_r", names(no_rd), value = TRUE),
  "ING guard (`ing_prior_guard.R`)" = grep("^\\.ing_", names(no_rd), value = TRUE),
  "C++ R wrappers (`rcpp_wrappers.R`)" = names(no_rd)[vapply(names(no_rd), function(x) identical(no_rd[[x]]$file, "rcpp_wrappers.R"), logical(1))],
  "Build, attach, misc." = intersect(
    c(".opencl_startup_quiet", ".opencl_runtime_sniff", ".opencl_startup_message",
      ".onAttach", "use_RcppParallel", "DIC_Info", "dpois2", "simfunction.default"),
    names(no_rd)
  )
)

out_int <- c(
  "# `R/` — internal helpers (undocumented in `man/`)",
  "",
  "Functions and symbols in **`R/`** with **`@noRd`**, **`@keywords internal`**, or",
  "no roxygen, without a dedicated help page. Intended for **`glmbayesCore:::`**",
  "or in-package use — not part of the exported API unless promoted.",
  "",
  "**Columns:** *File* is the defining source; *Called from* lists direct callers",
  "in `R/` (comma-separated). Helpers with no callers are marked *(unused)*.",
  "",
  "Companion: [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md).",
  "",
  sprintf("**lmebayes** resolves these via `getFromNamespace` / `:::`: %s.",
          paste(sprintf("`%s`", c(lme_int, ".mrglmb_normalize_pfamily_lists", ".validate_pfamily_for_rlmb")), collapse = ", ")),
  "",
  "---",
  ""
)

assigned <- character()
for (nm in names(helper_groups)) {
  fns <- intersect(setdiff(helper_groups[[nm]], assigned), names(no_rd))
  assigned <- c(assigned, fns)
  out_int <- c(out_int, write_helper_section(nm, fns))
}
rest <- setdiff(names(no_rd), assigned)
if (length(rest)) {
  out_int <- c(out_int, write_helper_section("Other internals", rest))
}

out_int <- c(out_int,
  "## Review checklist (internals)",
  "",
  "| Priority | Item |",
  "|----------|------|",
  "| 1 | Avoid new `@noRd` helpers unless tied to exported behavior. |",
  "| 2 | Keep `.mrglmb_normalize_pfamily_lists` / `.validate_pfamily_for_rlmb` stable for **lmebayes** `block_core_pfamily.R`. |",
  "| 3 | Remove or wire up *(unused)* helpers when touching related code. |",
  ""
)

writeLines(out_int, file.path(core_dir, "inst", "R_INTERNAL_HELPERS.md"))

out_idx <- c(
  "# `R/` function inventory (index)",
  "",
  "Maintainer index for symbols defined under **`R/`**. Split into two lists:",
  "",
  "| Document | Contents |",
  "|----------|----------|",
  "| **[R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md)** | `NAMESPACE` exports grouped by **glmbayes** / **lmebayes** overlap, S3 methods, doc topics. |",
  "| **[R_INTERNAL_HELPERS.md](R_INTERNAL_HELPERS.md)** | `@noRd` / `@keywords internal` helpers, C++ glue, attach hooks. |",
  "",
  "Scratch checks and one-off scripts live in `data-raw/` (not run by `test_check`).",
  ""
)
writeLines(out_idx, file.path(core_dir, "inst", "R_FUNCTION_INVENTORY.md"))

cat("Wrote inst/R_*.md\n")
cat("Exports:", length(core), " A:", length(cat_a), " B:", length(cat_b),
    " C:", length(cat_c), " D:", length(cat_d), " Helpers:", length(no_rd), "\n")
