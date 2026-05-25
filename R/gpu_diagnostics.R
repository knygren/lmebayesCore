#' GPU and OpenCL Diagnostics for glmbayes
#'
#' @description
#' A collection of tools for detecting GPU hardware, verifying OpenCL
#' availability, checking driver installation, validating environment
#' configuration, and diagnosing whether \pkg{glmbayes} can use GPU
#' acceleration. These functions provide both high-level diagnostic
#' summaries and low-level checks of system components such as PATH,
#' library directories, OpenCL headers, and the ICD loader.
#'
#' The diagnostic workflow is centered around
#' \code{diagnose_glmbayes()}, which orchestrates all other checks and
#' prints a detailed, human-readable report. Lower-level helpers can be
#' called individually for programmatic inspection or automated testing.
#' @param info A list returned by `detect_environment_and_gpus()`. The list must
#'   contain the following elements:
#'   \describe{
#'     \item{environment}{One of "windows", "msys2", "linux", "wsl", or "unknown".}
#'     \item{nvidia}{A list with elements `present` (logical) and `names` (character).}
#'     \item{amd}{A list with elements `present` (logical) and `names` (character).}
#'     \item{intel}{A list with elements `present` (logical) and `names` (character).}
#'   }
#' @param lib_dirs A list of OpenCL directories
#' @param runtime_info The structured list returned by \code{detect_compute_runtimes()}.
#' @return A structured list with diagnostics for each runtime, including:
#'   \itemize{
#'     \item installed: whether the runtime was detected
#'     \item found_path_dirs: directories already present in PATH
#'     \item missing_path_dirs: directories that should be added to PATH
#'     \item found_lib_dirs: directories already present in LD_LIBRARY_PATH
#'     \item missing_lib_dirs: directories that should be added to LD_LIBRARY_PATH
#'     \item include_dirs: include directories detected for headers
#'   }
#' @section High-level diagnostic:
#' \itemize{
#'   \item \code{diagnose_glmbayes()} — full GPU/OpenCL diagnostic report.
#' }
#'
#' @section Environment and hardware detection:
#' \itemize{
#'   \item \code{detect_environment_and_gpus()} — detect OS and GPU vendor.
#'   \item \code{gpu_names()} — enumerate available GPU device names.
#'   \item \code{detect_compute_runtimes()} — detect CUDA/OpenCL runtimes.
#' }
#'
#' @section OpenCL availability and runtime checks:
#' \itemize{
#'   \item \code{has_opencl()} — quick check for OpenCL support.
#'   \item \code{verify_opencl_runtime()} — probe OpenCL platform/device availability.
#'   \item \code{check_runtime_env()} — validate PATH and library directories.
#' }
#'
#' @section Driver installation helpers:
#' \itemize{
#'   \item \code{detect_or_install_gpu_drivers()} — detect driver presence and issues.
#' }
#'
#' @section PATH and library path utilities:
#' These are optional helpers used by the diagnostic pipeline.
#' \itemize{
#'   \item \code{add_to_path_windows()}
#'   \item \code{add_to_path_linux()}
#'   \item \code{add_to_libpath_linux()}
#' }
#'
#' @details
#' GPU acceleration speeds up **envelope construction and grid evaluation**
#' (e.g. large \eqn{3^p} grids or many tangency evaluations) when you pass
#' \code{use_opencl = TRUE} in modeling and envelope functions such as
#' \code{\link{glmb}} and \code{\link{rglmb}}. OpenCL is **vendor-neutral**
#' (NVIDIA, AMD, Intel); CPU-only builds remain valid and are often used when
#' no OpenCL stack is present.
#'
#' **Practical setup (summary).** The CRAN release and prebuilt R-Universe
#' binaries are built **without** OpenCL GPU support; enabling the GPU path
#' usually requires installing the development \pkg{glmbayes} **from source**
#' on a machine
#' with OpenCL **headers**, a linkable **OpenCL library / ICD loader**, and
#' a working **vendor runtime** (GPU driver). You need a normal C/C++
#' toolchain (e.g. Rtools on Windows, \code{build-essential} and
#' \code{r-base-dev} on Linux, Xcode CLT plus GCC on macOS for source installs).
#' Vendor-specific notes (CUDA Toolkit vs Intel SDK vs Khronos headers on
#' Windows, \code{opencl-headers} and \code{ocl-icd} packages on Linux, etc.)
#' are spelled out in \insertCite{glmbayesChapter12}{glmbayes}.
#'
#' **What this help page checks.** A usable OpenCL environment requires:
#' \enumerate{
#'   \item OpenCL headers (e.g., \code{CL/cl.h}) at compile time,
#'   \item the OpenCL ICD loader (e.g., \code{libOpenCL.so.1}) at runtime,
#'   \item correct PATH and library search paths (especially on Linux/WSL),
#'   \item a functional OpenCL platform and device (driver installed).
#' }
#' The functions here inspect these pieces. On Linux and WSL,
#' \code{verify_opencl_runtime()} tries to create a platform, device, context,
#' queue, and compile a minimal kernel. On Windows, that probe is skipped
#' because platform-creation failures are often uninformative; rely on
#' \code{diagnose_glmbayes()} and driver/runtime detection instead.
#'
#' Start with \code{\link{diagnose_glmbayes}()} for a single readable report;
#' use \code{\link{has_opencl}()} for a quick boolean when scripting.
#'
#' @return
#' Most functions return structured lists describing detected hardware,
#' drivers, runtimes, or environment issues. \code{diagnose_glmbayes()}
#' prints a formatted report and invisibly returns a named list containing
#' all intermediate diagnostic results.
#'
#' @seealso
#' \code{\link{diagnose_glmbayes}},
#' \code{\link{detect_environment_and_gpus}},
#' \code{\link{detect_compute_runtimes}},
#' \code{\link{verify_opencl_runtime}},
#' \code{\link{has_opencl}}.
#'
#' Modeling with \code{use_opencl}: \code{\link{glmb}}, \code{\link{rglmb}}.
#' Envelope helpers: \code{\link{EnvelopeBuild}}, \code{\link{EnvelopeEval}}.
#'
#' Full install and troubleshooting: \code{vignette("Chapter-12", package = "glmbayes")}
#' (\insertCite{glmbayesChapter12}{glmbayes}); implementation notes:
#' \insertCite{glmbayesChapterA10}{glmbayes}.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @keywords diagnostics gpu opencl environment
#' @name gpu_diagnostics
NULL


#' @export
#' @rdname gpu_diagnostics
#' @order 1
diagnose_glmbayes <- function() {
  cat("=== glmbayes Diagnostic Report ===\n")

  # Step 1: Environment + GPU detection (opencltools)
  info     <- detect_environment_and_gpus()
  drivers  <- detect_or_install_gpu_drivers(info)
  runtimes <- detect_compute_runtimes(info)
  env_diag <- check_runtime_env(runtimes)

  cat("Environment:", info$environment, "\n\n")

  # Step 2: Preference order (NVIDIA > AMD > Intel)
  gpu_vendor <- if (info$nvidia$present) "nvidia"
  else if (info$amd$present) "amd"
  else if (info$intel$present) "intel"
  else NULL

  diag <- NULL
  runtime_ok <- NA

  if (!is.null(gpu_vendor)) {
    cat("GPU:", toupper(gpu_vendor), "\n")
    drv  <- drivers$drivers[[gpu_vendor]]
    rt   <- runtimes$runtimes[[gpu_vendor]]
    diag <- env_diag$diagnostics[[gpu_vendor]]

    # 1. Driver check
    if (drv$installed) {
      cat("  [OK] Driver installed\n")
    } else {
      cat("  [FAIL] Driver not installed\n")
      if (length(drv$issues) > 0)
        cat("    Issues:", paste(drv$issues, collapse=", "), "\n")
    }

    # 2. OpenCL header/runtime presence
    hdr <- rt$opencl$headers_present
    rtm <- rt$opencl$runtime_present
    inst <- rt$opencl$installed

    if (hdr) {
      cat("  [OK] OpenCL headers found (CL/cl.h)\n")
    } else {
      cat("  [FAIL] OpenCL headers not found (CL/cl.h missing)\n")
    }

    if (rtm) {
      cat("  [OK] OpenCL runtime found (OpenCL.dll / ICD)\n")
    } else {
      cat("  [FAIL] OpenCL runtime not found\n")
    }

    if (inst) {
      cat("  [OK] OpenCL fully available (headers + runtime)\n")
    } else {
      cat("  [FAIL] OpenCL incomplete (missing headers or runtime)\n")
    }

    # 3. PATH/lib environment validation
    paths_ok <- (length(diag$opencl$missing_path_dirs) == 0 &&
                   length(diag$opencl$missing_lib_dirs) == 0)

    if (paths_ok) {
      cat("  [OK] Required PATH and library dirs present\n")
    } else {
      if (length(diag$opencl$missing_path_dirs) > 0)
        cat("  [WARN] Missing PATH entries:",
            paste(diag$opencl$missing_path_dirs, collapse=", "), "\n")
      if (length(diag$opencl$missing_lib_dirs) > 0)
        cat("  [WARN] Missing library dirs:",
            paste(diag$opencl$missing_lib_dirs, collapse=", "), "\n")
    }

    # 4. Runtime probe (Linux/WSL only)
    if (paths_ok && tolower(info$environment) %in% c("linux", "wsl")) {
      runtime_ok <- verify_opencl_runtime(rt$opencl$lib_dirs)
      if (runtime_ok) {
        cat("  [OK] OpenCL runtime probe succeeded (platform available)\n")
      } else {
        cat("  [FAIL] OpenCL runtime probe failed (no usable platform)\n")
      }
    } else if (!paths_ok) {
      cat("  [SKIP] Runtime probe skipped (missing PATH/lib dirs)\n")
    } else {
      cat("  [SKIP] Runtime probe skipped on Windows\n")
    }

  } else {
    cat("[FAIL] No supported GPU detected. glmbayes will run in CPU-only mode.\n")
  }

  # Step 3: Report compile-time OpenCL status (glmbayes build)
  opencl_enabled <- has_opencl()
  if (opencl_enabled) {
    cat("\n[OK] glmbayes was compiled with OpenCL support.\n")
  } else {
    cat("\n[FAIL] glmbayes was compiled without OpenCL support.\n")
  }

  # Step 4: Interactive PATH/lib fixes
  missing_items <- !is.null(diag) &&
    (length(diag$opencl$missing_path_dirs) > 0 ||
       length(diag$opencl$missing_lib_dirs) > 0)

  if (missing_items && !isTRUE(opencl_enabled)) {
    cat("\n[INFO] Missing PATH/lib entries detected and OpenCL is not enabled.\n")

    if (length(diag$opencl$missing_path_dirs) > 0) {
      cat("  Missing PATH entries:\n")
      cat("   -", paste(diag$opencl$missing_path_dirs, collapse="\n   - "), "\n")
      ans <- readline("Would you like to permanently add missing PATH dirs? [y/N]: ")
      if (tolower(ans) == "y") {
        if (tolower(info$environment) == "windows") {
          add_to_path_windows(diag$opencl$missing_path_dirs)
        } else {
          add_to_path_linux(diag$opencl$missing_path_dirs)
        }
      }
    }

    if (length(diag$opencl$missing_lib_dirs) > 0 &&
        tolower(info$environment) %in% c("linux", "wsl")) {
      cat("  Missing library dirs:\n")
      cat("   -", paste(diag$opencl$missing_lib_dirs, collapse="\n   - "), "\n")
      ans_lib <- readline("Would you like to permanently add missing library dirs to LD_LIBRARY_PATH? [y/N]: ")
      if (tolower(ans_lib) == "y") {
        add_to_libpath_linux(diag$opencl$missing_lib_dirs)
      }
    }
  }

  cat("\n=== End of Diagnostic Report ===\n")

  invisible(list(
    environment_info      = info,
    driver_status         = drivers,
    runtime_status        = runtimes,
    env_diag              = env_diag,
    opencl_runtime_probe  = runtime_ok,
    opencl_enabled        = opencl_enabled
  ))
}


#' @export
#' @rdname gpu_diagnostics
#' @order 6
has_opencl <- function() {
  .has_opencl_cpp()
}
