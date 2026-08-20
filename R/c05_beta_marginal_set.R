## Beta-marginal level-set helpers (BETA_MARGINAL_MODE_LEVELSET.md).
## @noRd helpers only; no exports until beta_marginal_safe_set() is approved.

#' Gaussian-reference deficiency level r_Gauss(n, delta_2).
#'
#' Returns r = 0.5 * qchisq(1 - delta_2, df = n), the level at which a
#' N(beta_dagger, H^{-1}) Laplace reference has tail mass delta_2 on
#' {Xi_Lap > r}. Same convention as \code{.c05_mvn_d()}.
#'
#' @param n State dimension (typically J * p_re).
#' @param delta_2 Tail mass budget in (0, 1).
#' @return List with \code{r}, \code{r2 = 2*r}, \code{n}, \code{delta_2},
#'   \code{s = sqrt(2*r)}, and \code{laplace_tail_mass} (= \code{delta_2}).
#' @noRd
.c05_beta_r_gauss_level <- function(n, delta_2) {
  n <- as.integer(n)
  if (n < 1L) {
    stop("'n' must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(delta_2) || delta_2 <= 0 || delta_2 >= 1) {
    stop("'delta_2' must lie in (0, 1).", call. = FALSE)
  }
  r2 <- stats::qchisq(1 - delta_2, df = n)
  r <- 0.5 * r2
  list(
    n = n,
    delta_2 = delta_2,
    r2 = r2,
    r = r,
    s = sqrt(r2),
    laplace_tail_mass = delta_2
  )
}

#' Laplace-reference tail mass P(Xi_Lap > r) for N(0, I_n) deficiency.
#'
#' For exact Gaussian with Xi = ||beta||^2/2 ~ Gamma(n/2, 1), this equals
#' pgamma(r, n/2, lower.tail = FALSE). Used to compare r_Gauss vs Prop (P2).
#'
#' @param r Deficiency level.
#' @param n State dimension.
#' @noRd
.c05_beta_laplace_tail_mass <- function(r, n) {
  stats::pgamma(r, shape = n / 2, rate = 1, lower.tail = FALSE)
}

#' Proposition (P2) tail bound ratio Gamma(n,r)/gamma(n,r) at level r.
#'
#' Legacy convex worst-case mass upper bound; not used for r_Gauss calibration.
#'
#' @param r Deficiency level.
#' @param n State dimension.
#' @noRd
.c05_beta_prop2_tail_bound <- function(r, n) {
  stats::pgamma(r, n, lower.tail = FALSE) /
    stats::pgamma(r, n, lower.tail = TRUE)
}

#' Invert Proposition (P2) for a per-dimension budget (legacy group_precision_floor).
#'
#' @param d Dimension passed to pgamma shape (Prop 2 convention).
#' @param epsilon Budget in (0, 1).
#' @noRd
.c05_beta_prop2_level <- function(d, epsilon) {
  if (!is.finite(epsilon) || epsilon <= 0 || epsilon >= 1) {
    stop("'epsilon' must be in (0, 1).", call. = FALSE)
  }
  d <- as.numeric(d)
  if (d < 1) stop("'d' must be at least 1.", call. = FALSE)
  uniroot(
    function(r) .c05_beta_prop2_tail_bound(r, d) - epsilon,
    c(1e-8, 200 * d + 200),
    tol = 1e-10
  )$root
}

#' Find marginal beta mode by Newton on the gamma-integrated posterior.
#'
#' Uses \code{.group_floor_engine()} with \code{b_dag = NULL} (Newton from mu_beta).
#'
#' @param prior_int Output of \code{.group_floor_integrate_gamma()}.
#' @param group_data Output of \code{.group_floor_prepare_group_data()}.
#' @param family_hook Output of \code{.group_floor_family()}.
#' @return List with \code{beta}, \code{hessian}, \code{f_mode}, \code{engine}.
#' @noRd
.c05_beta_marginal_mode <- function(prior_int, group_data, family_hook) {
  engine <- .group_floor_engine(
    prior_int = prior_int,
    group_data = group_data,
    family_hook = family_hook,
    b_dag = NULL
  )
  list(
    beta = engine$b_dag,
    hessian = engine$He_f(engine$b_dag),
    f_mode = engine$f_dag,
    engine = engine
  )
}
