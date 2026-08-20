## Coupling eigenvalue spectrum for C05 deficiency calibration.

#' Symmetrized coupling spectrum from refresh precision and Jacobian.
#'
#' With \eqn{S = P_{11}\tilde J}, \eqn{A = P_{11}^{-1/2} S P_{11}^{-1/2}}, and
#' \eqn{B = A(I-A)^{-1}}, the returned \code{kappa} are \eqn{\mathrm{eig}(A)}
#' and \code{weights} are \eqn{\mathrm{eig}(B) = \kappa_i/(1-\kappa_i)} (spectral
#' identity; no \eqn{P_{11}}-commute assumption).
#' @noRd
.c05_coupling_spectrum <- function(p11, tilde_J, strict = FALSE) {
  if (is.null(p11$P11) || is.null(p11$chol_P11)) {
    stop("'p11' must contain 'P11' and 'chol_P11'.", call. = FALSE)
  }
  if (!is.matrix(tilde_J)) {
    stop("'tilde_J' must be a matrix.", call. = FALSE)
  }
  q <- nrow(tilde_J)
  if (ncol(tilde_J) != q || q != nrow(p11$P11)) {
    stop("'tilde_J' and 'p11$P11' must have matching dimensions.", call. = FALSE)
  }

  S <- p11$P11 %*% tilde_J
  kappa <- .two_block_gen_eigen(S, p11$P11, strict = strict)
  weights <- kappa / (1 - kappa)
  kappa_max <- max(kappa)

  list(
    kappa = kappa,
    weights = weights,
    kappa_max = kappa_max,
    rho = if (kappa_max < 1) 1 / (1 - kappa_max) else Inf,
    q = length(kappa),
    method = "symmetrized_coupling"
  )
}

#' Coupling eigenvalue spectrum at the population mode.
#'
#' Returns symmetrized coupling eigenvalues \eqn{\kappa_i = \mathrm{eig}(A)} and
#' deficiency weights \eqn{w_i = \kappa_i/(1-\kappa_i) = \mathrm{eig}(B)} with
#' \eqn{A = P_{11}^{-1/2} S P_{11}^{-1/2}}, \eqn{S = P_{11}\tilde J},
#' \eqn{B = A(I-A)^{-1}} (spectrally; matrix product only when \eqn{A} and
#' \eqn{P_{11}} share eigenvectors). Computed at the terminal conditional-mean
#' E-step of \code{\link{population_mode}}.
#'
#' Limits: \eqn{\kappa_i \to 0 \Rightarrow w_i \to 0 \Rightarrow d \to 0};
#' \eqn{\kappa_i \to 1 \Rightarrow w_i \to \infty \Rightarrow d \to \infty}.
#'
#' @param mode A \code{\link{population_mode}} result with \code{p11} and
#'   \code{tilde_J}.
#' @return A list with \code{kappa}, \code{weights}, \code{kappa_max},
#'   \code{rho}, \code{q}, and \code{method}.
#' @seealso \code{\link{deficiency_calibrate}}, \code{\link{epsilon}}
#' @export
deficiency_spectrum <- function(mode) {
  if (is.null(mode$p11) || is.null(mode$tilde_J)) {
    stop(
      "'mode' must be a population_mode() result with 'p11' and 'tilde_J'.",
      call. = FALSE
    )
  }
  .c05_coupling_spectrum(mode$p11, mode$tilde_J, strict = FALSE)
}

#' Calibrate reference radius r(delta) and deficiency d = r^2/2.
#'
#' Inverts the Gaussian-reference escape tail on \eqn{B}-weights,
#' \deqn{\Pr\Bigl(\sum_i w_i Z_i^2 > r^2\Bigr) = \delta,
#' \qquad w_i = \frac{\kappa_i}{1-\kappa_i},}
#' then sets \code{d = r^2/2} (equivalently \code{2d = r^2}) and
#' \code{set_cut = exp(-d)}.
#'
#' @param delta Tail / escape probability budget in \code{(0, 1)}.
#' @param spectrum A list from \code{\link{deficiency_spectrum}} with
#'   \code{weights}.
#' @details
#' Weighted-\eqn{\chi^2} tails for general spectra use
#' \code{CompQuadForm::imhof()} when the suggested package \pkg{CompQuadForm}
#' is installed, otherwise \code{mgcv::psum.chisq()} if \pkg{mgcv} is available.
#' @return A list with \code{r2}, \code{r}, \code{d}, \code{delta},
#'   \code{escape_mass}, \code{set_cut}, and \code{spectrum}.
#' @seealso \code{\link{epsilon}}, \code{\link{deficiency_spectrum}}
#' @export
deficiency_calibrate <- function(delta, spectrum) {
  .c05_validate_delta(delta)
  if (is.null(spectrum$weights)) {
    stop("'spectrum' must contain 'weights' (from deficiency_spectrum()).",
         call. = FALSE)
  }

  weights <- .c05_validate_weights(
    spectrum$weights,
    fn_name = "deficiency_calibrate"
  )
  rad <- .c05_r_from_delta(delta, weights)
  escape_mass <- .c05_weighted_chisq_tail(rad$r2, weights)

  list(
    delta = delta,
    r2 = rad$r2,
    r = rad$r,
    d = rad$d,
    escape_mass = escape_mass,
    set_cut = exp(-rad$d),
    spectrum = spectrum
  )
}
