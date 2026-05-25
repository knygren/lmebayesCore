# Re-exports from opencltools (Tier 3: runtime/system diagnostics).

# ---- Tier 3: gpu_diagnostics (system / runtime; not package compile-flag probes) ----

#' @export
#' @rdname gpu_diagnostics
#' @order 2
detect_environment_and_gpus <- opencltools::detect_environment_and_gpus

#' @export
#' @rdname gpu_diagnostics
#' @order 3
gpu_names <- opencltools::gpu_names

#' @export
#' @rdname gpu_diagnostics
#' @order 4
detect_or_install_gpu_drivers <- opencltools::detect_or_install_gpu_drivers

#' @export
#' @rdname gpu_diagnostics
#' @order 5
detect_compute_runtimes <- opencltools::detect_compute_runtimes

#' @export
#' @rdname gpu_diagnostics
#' @order 7
verify_opencl_runtime <- opencltools::verify_opencl_runtime

#' @export
#' @rdname gpu_diagnostics
#' @order 8
check_runtime_env <- opencltools::check_runtime_env

# ---- Tier 3: PATH helpers ----

#' @title Add Directories to PATH or LD_LIBRARY_PATH
#' @description
#' These helper functions allow you to add missing directories to the PATH
#' or library search environment variables in a permanent way, minimizing
#' manual editing.
#' @details
#' - On **Windows**, updates the user-level PATH via PowerShell.
#' - On **Linux/WSL**, appends export lines to ~/.bashrc for PATH or LD_LIBRARY_PATH.
#' @param dirs Character vector of directories to add.
#' @return No return value; called for side effects.
#' @seealso [Sys.getenv], [Sys.setenv]
#' @name add_to_path
#' @rdname add_to_path
NULL

#' @export
#' @rdname add_to_path
add_to_path_windows <- opencltools::add_to_path_windows

#' @export
#' @rdname add_to_path
add_to_path_linux <- opencltools::add_to_path_linux

#' @export
#' @rdname add_to_path
add_to_libpath_linux <- opencltools::add_to_libpath_linux
