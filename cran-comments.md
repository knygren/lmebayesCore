# CRAN resubmission comments — glmbayes 0.9.2

Thank you for the review. This is a resubmission addressing the
requested changes.

Each flagged item was tracked as a GitHub issue and addressed with a
dedicated GitHub push. The issues contain additional implementation
details for each update.

## Response to reviewer comments

- Expanded the DESCRIPTION Description field into a fuller paragraph
  describing the package purpose, supported model families, iid posterior
  sampling approach, diagnostics, simulation tools, vignettes, and
  optional OpenCL acceleration.

- Added the method reference in DESCRIPTION using CRAN-preferred
  auto-linking format:
  Nygren and Nygren (2006) <doi:10.1198/016214506000000357>.

- Added the missing \value{} documentation for summary.rgamma_reg.Rd,
  including the output class/structure and print method behavior.

- Removed commented-out executable code from examples, including the
  examples underlying glmb.influence.measures.Rd and influence.glmb.Rd.
  Related commented executable fragments found during the package scan
  were also removed or rewritten as prose comments.

- Replaced \dontrun{} examples with \donttest{}. 

- Reviewed examples, vignettes, and demos for changes to par() and
  options(). Added save/restore blocks for graphical parameters and
  options in the affected files.

- Updated Authors@R to include contributors and copyright holders for
  code derived from R Mathlib and R stats, including The R Core Team,
  The R Foundation, Ross Ihaka, Robert Gentleman, Simon Davies, Morten
  Welinder, and Martin Maechler, with appropriate ctb/cph roles.

- Added inst/COPYRIGHTS documenting bundled/adapted R Mathlib code,
  OpenCL ports/adaptations of selected R Mathlib routines and support
  code, and R stats-derived lm/glm modeling conventions. Added
  source-level notices in affected R and OpenCL files pointing to this
  documentation.

## Test environments

- Local check:
  devtools::check(vignettes = TRUE, args = "--as-cran", remote = TRUE,
                  manual = TRUE)
  - Duration: 9m 33.4s
  - Result: 0 errors, 0 warnings, 1 note

- macOS release builder, check_mac_release:
  - Result: Status OK; clean install.
  - mac-builder reports final status rather than the standard
    errors/warnings/notes summary.

- macOS devel builder, check_mac_devel:
  - Result: Status OK; clean install.
  - mac-builder reports final status rather than the standard
    errors/warnings/notes summary.

- Windows release builder, win-builder:
  - R version: 4.6.0 (2026-04-24 ucrt)
  - Installation time: 550 seconds
  - Check time: 297 seconds
  - Result: Status 1 NOTE

- Windows old-release builder, win-builder:
  - Result: Status 1 NOTE

- rhub checks:
  - 27 platforms checked; all jobs succeeded.
  - 22 platforms reported OK.
  - 5 clang-based platforms reported 1 NOTE.
  - The clang-platform notes match the prior toolchain-note pattern
    observed for this package and are not from package Makevars.

## Primary R CMD check result

0 errors | 0 warnings | 1 note

The remaining note in the local, win-builder release, and win-builder
old-release checks is the expected CRAN incoming feasibility note for
new submissions:

    Maintainer: 'Kjell Nygren <kjell.a.nygren@gmail.com>'

    New submission

---
_This file is listed in `.Rbuildignore` and is not included in the
built source tarball. When submitting, paste the content above into the
"Optional comments" field on the CRAN submission form at
https://cran.r-project.org/submit.html_
