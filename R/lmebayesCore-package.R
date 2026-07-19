#' @aliases lmebayesCore
#'
#' @title lmebayesCore: Core C++ Sampling Engine for lmebayes
#'
#' @description
#' Core C++ engine for envelope-based iid samplers, two-block Gibbs
#' mixed-model engines, and optional OpenCL acceleration. Full-featured
#' developer backend for \pkg{lmebayes} and related extensions. End users
#' should install \pkg{lmebayes} for lmer/glmer-style mixed-model workflows.
#'
#' @details
#' Low-level entry points include envelope construction
#' (\code{\link[glmbayesCore]{EnvelopeBuild}}, \code{\link[glmbayesCore]{EnvelopeOrchestrator}}),
#' registered simulation pipelines (\code{\link[glmbayesCore]{rNormalGLM_std}},
#' \code{\link[glmbayesCore]{rIndepNormalGammaReg_std}}), matrix-input samplers
#' (\code{\link[glmbayesCore]{rglmb}}, \code{\link[glmbayesCore]{rlmb}}), the two-block Gibbs
#' mixed-model engines (\code{\link{rlmerb}}, \code{\link{rglmerb}}), and
#' OpenCL kernel loaders. Formula interfaces \code{lmerb} and \code{glmerb}
#' live in \pkg{lmebayes}; \pkg{lmebayesCore} supplies the sampling engine.
#'
#' @section OpenCL startup checks:
#' In interactive sessions, attaching the package with \code{library(lmebayesCore)}
#' may emit a short \code{\link{packageStartupMessage}} when \code{glmbayesCore_has_opencl()}
#' is \code{FALSE} but a GPU or OpenCL stack appears available on the host.
#' Set \code{options(glmbayes.quiet_opencl_startup = TRUE)} to suppress attach
#' notes (recommended for CI and \command{R CMD check}).
#'
#' @seealso \pkg{lmebayes} for the end-user modelling package.
#'
#' @references
#' \insertAllCited{}
#'
#' @import stats Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rdpack reprompt
#' @importFrom utils flush.console
#' @import opencltools
#' @import nmathopencl
#' @useDynLib lmebayesCore, .registration = TRUE
"_PACKAGE"
