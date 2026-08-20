## Chapter C05 Stages 3--4: spectrum-calibrated set sizing and epsilon multiplier.

#' MVN reference deficiency d from tail budget delta (legacy equal-weight route).
#' @noRd
.c05_mvn_d <- function(q, delta) {
  if (!(delta > 0 && delta < 1)) {
    stop("'delta' must lie in (0, 1).", call. = FALSE)
  }
  if (q < 1L) {
    stop("'q' must be a positive integer.", call. = FALSE)
  }
  r2 <- stats::qchisq(1 - delta, df = q)
  d <- 0.5 * r2
  list(
    delta = delta,
    q = as.integer(q),
    r2 = r2,
    r = sqrt(r2),
    d = d,
    set_cut = exp(-d)
  )
}

#' Kernel floor and Theorem~2 multiplier from eps_star, d, and spectrum.
#' @noRd
.c05_minorization_constants <- function(eps_star, d, spectrum) {
  Q_mass_lb <- .c05_weighted_chisq_cdf(2 * d, spectrum$kappa)
  eps_d <- exp(-d) * eps_star
  eps <- eps_d * Q_mass_lb
  list(
    eps_star = eps_star,
    d = d,
    eps_d = eps_d,
    Q_mass_lb = Q_mass_lb,
    eps = eps,
    spectrum = spectrum,
    kappa = spectrum$kappa,
    weights = spectrum$weights
  )
}

#' Legacy MVN-calibrated reference radius (equal-weight chi-square shortcut).
#'
#' Computes \code{d = 0.5 * qchisq(1 - delta, q)}. This equal-weight route is
#' retained for reference only; \code{\link{epsilon}} and
#' \code{\link{certificate}} use \code{\link{deficiency_calibrate}} on the
#' coupling spectrum from \code{\link{deficiency_spectrum}} instead.
#'
#' @param q Hyperparameter dimension.
#' @param delta Tail / escape probability budget in \code{(0, 1)}.
#' @return A list with \code{r2}, \code{r}, \code{d}, and \code{set_cut = exp(-d)}.
#' @seealso \code{\link{deficiency_calibrate}}, \code{\link{epsilon}}
#' @export
mvn_calibrate <- function(q, delta) {
  .c05_mvn_d(q, delta)
}

#' Restricted-chain Doeblin multiplier epsilon from spectrum-calibrated sizing.
#'
#' Combines the mode profile \code{eps_star} with the reference radius
#' \code{r(delta)} and deficiency \code{d = r^2/2} from the \eqn{B}-spectrum
#' weights \eqn{w_i = \kappa_i/(1-\kappa_i)} at \code{mode}. Produces
#' \code{eps_d = exp(-d) eps_star} and \code{eps = eps_d Q(C_tilde_d)} with
#' \code{Q(C_tilde_d)} at \code{Pr(sum kappa_i Z_i^2 <= r^2)}.
#'
#' @param eps_star Mode profile \eqn{\varepsilon(\gamma^\star)} from
#'   \code{\link{epsilon_star}}.
#' @param delta Tail budget used to calibrate \code{d}.
#' @param mode A \code{\link{population_mode}} result with \code{p11} and
#'   \code{tilde_J}. The spectrum uses the terminal conditional-mean E-step.
#' @details
#' Weighted-\eqn{\chi^2} tails for general coupling spectra require suggested
#' package \pkg{CompQuadForm} (preferred) or \pkg{mgcv}.
#' @return A list with \code{eps_star}, \code{d}, \code{r}, \code{r2},
#'   \code{eps_d}, \code{Q_mass_lb}, \code{eps}, \code{spectrum}, \code{kappa},
#'   \code{cal}, and \code{route}.
#' @export
epsilon <- function(eps_star, delta, mode) {
  if (!(eps_star > 0 && eps_star <= 1)) {
    stop("'eps_star' must lie in (0, 1].", call. = FALSE)
  }
  if (is.null(mode$p11) || is.null(mode$tilde_J)) {
    stop(
      "'mode' must be a population_mode() result with 'p11' and 'tilde_J'.",
      call. = FALSE
    )
  }

  spectrum <- deficiency_spectrum(mode)
  cal <- deficiency_calibrate(delta, spectrum)
  const <- .c05_minorization_constants(eps_star, cal$d, spectrum)
  const$r <- cal$r
  const$r2 <- cal$r2
  const$cal <- cal
  const$route <- "spectrum_calibrated"
  const
}
