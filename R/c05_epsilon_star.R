## Stage 2: epsilon(gamma*) profile (closure and optimize routes).



#' Closure value det(I + tilde J)^{-1/2}.

#' @noRd

.c05_epsilon_closure <- function(tilde_J) {

  q <- nrow(tilde_J)

  exp(-0.5 * as.numeric(determinant(diag(q) + tilde_J, logarithm = TRUE)$modulus))

}



#' Mode profile epsilon(gamma*) for the restricted Gibbs certificate.

#'

#' @param mode Object returned by \code{\link{population_mode}}.

#' @param method \code{"closure"} uses \eqn{\det(I+\tilde J)^{-1/2}} when

#'   available (Gaussian reference); \code{"optimize"} calls

#'   \code{\link{epsilon_optimize}}.

#' @param ... Passed to \code{\link{epsilon_optimize}} when

#'   \code{method = "optimize"}.

#' @return A list with \code{eps_star}, \code{method}, and \code{certified}.

#' @export

epsilon_star <- function(mode,

                         method = c("closure", "optimize"),

                         ...) {

  method <- match.arg(method)

  if (!is.list(mode) || is.null(mode$fixef)) {

    stop("'mode' must be the result of population_mode().",

         call. = FALSE)

  }



  if (identical(method, "optimize")) {

    return(epsilon_optimize(mode, ...))

  }



  if (is.null(mode$eps_star_closure)) {

    stop(

      "closure epsilon(gamma*) is unavailable; use method = \"optimize\".",

      call. = FALSE

    )

  }



  is_gaussian <- !is.null(mode$family) &&

    identical(mode$family$family, "gaussian")

  is_exact <- identical(mode$estep, "exact")



  list(

    eps_star = mode$eps_star_closure,

    method = "closure",

    certified = isTRUE(is_gaussian && is_exact),

    g_opt = log(mode$eps_star_closure),

    gamma_prime = NULL

  )

}


