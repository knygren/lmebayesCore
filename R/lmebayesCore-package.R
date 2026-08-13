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
#' @section Vignette index (the B series):
#' The \code{Chapter-B*} vignettes document the engine itself: the theory it
#' implements, the routes it selects between, and the diagnostics it emits.
#' They assume the reader has met the models already, for which see the
#' \code{Chapter-00} to \code{Chapter-13} curriculum in \pkg{lmebayes}.
#' \describe{
#'   \item{\code{vignette("Chapter-B00")}}{Roadmap, notation, and the
#'     correspondence to \pkg{lme4}.}
#'   \item{\code{vignette("Chapter-B01")}}{The two-block architecture:
#'     replicate chains, the sweep loop, and the call stack.}
#'   \item{\code{vignette("Chapter-B02")}}{Ergodicity and total-variation
#'     convergence; the rate \eqn{\lambda^*} and the rank conditions.}
#'   \item{\code{vignette("Chapter-B03")}}{Calibrating \code{m_convergence}
#'     from \code{tv_tol}, and the pilot-versus-main cost trade-off.}
#'   \item{\code{vignette("Chapter-B04")}}{The known-vcov Gaussian route and
#'     the exact iid sampler.}
#'   \item{\code{vignette("Chapter-B05")}}{Estimated variance components
#'     (\code{dIndependent_Normal_Gamma} Block 2).}
#'   \item{\code{vignette("Chapter-B06")}}{The extended local rate
#'     diagnostic.}
#'   \item{\code{vignette("Chapter-B07")}}{Per-group measurement dispersion.}
#'   \item{\code{vignette("Chapter-B08")}}{Design, identifiability and rank
#'     conditions.}
#'   \item{\code{vignette("Chapter-B09")}}{Row-block engines.}
#'   \item{\code{vignette("Chapter-B10")}}{Sweep history and convergence
#'     plots.}
#' }
#'
#' @section Vignette index (the C series):
#' The \code{Chapter-C*} vignettes take one model class each and work through
#' its posterior algebra, the eigenvalue calculation behind the convergence
#' rate \eqn{\lambda^*}{lambda*}, and simulation results with convergence
#' plots. Where the B series describes what the engine does, the C series
#' derives why it is correct and then shows the numbers.
#' \describe{
#'   \item{\code{vignette("Chapter-C01")}}{The Gaussian running example: data,
#'     model, design, rank restriction and calibrated priors, shared by every
#'     Gaussian chapter in the series.}
#'   \item{\code{vignette("Chapter-C02")}}{The exact Gaussian iid sampler:
#'     known \eqn{\Psi}{Psi} and \eqn{\sigma^2}{sigma^2}, closed-form
#'     posterior, no Markov chain.}
#'   \item{\code{vignette("Chapter-C03")}}{Total variation distances between
#'     multivariate normal densities: Nygren (2020) reproduced in full ---
#'     the generalized error functions, Lemmas 1--2, Theorems 1--3,
#'     Claims 1--4, Corollary 1 and their proofs. Model-free.}
#'   \item{\code{vignette("Chapter-C04")}}{The same Gaussian model via the
#'     two-block Gibbs sampler: the Chapter-C03 rate and sweep budget
#'     instantiated on the running design, sweep history, and
#'     variance-convergence diagnostics.}
#' }
#'
#' The \code{Chapter-A*} vignettes shipped in the sources are inherited iid
#' GLM reference material from \pkg{glmbayes} and are not built into the
#' package tarball.
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
