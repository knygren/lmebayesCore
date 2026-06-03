#' Prior setup for multiple Gaussian responses
#'
#' @inheritParams Prior_Setup
#' @return A named list of class \code{"multi_PriorSetup"}. Each element is a
#'   \code{\link{Prior_Setup}} result for one column of the response (names from
#'   \code{colnames(y)} or \code{Y1}, \code{Y2}, \ldots).
#' @family prior
#' @export
multi_prior_setup <- function(
    formula,
    family = gaussian(),
    data = NULL,
    weights = NULL,
    subset = NULL,
    na.action = na.fail,
    offset = NULL,
    contrasts = NULL,
    pwt = NULL,
    pwt_default_low = 0.01,
    pwt_default_high = 0.05,
    n_prior = NULL,
    sd = NULL,
    dispersion = NULL,
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model"),
    mu = NULL,
    k = 1,
    ...
) {
  call <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family) || family$family != "gaussian") {
    stop(
      "multi_prior_setup() currently supports family = gaussian() only.",
      call. = FALSE
    )
  }

  if (missing(data)) {
    data <- environment(formula)
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(
    c("formula", "data", "subset", "weights", "na.action", "offset"),
    names(mf),
    0L
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- as.matrix(model.response(mf, "any"))
  l1 <- ncol(Y)
  if (l1 < 1L) {
    stop("formula must specify at least one response column.", call. = FALSE)
  }
  coef_names <- colnames(Y)
  if (is.null(coef_names) || length(coef_names) != l1) {
    coef_names <- paste0("Y", seq_len(l1))
  }

  termlabels <- attr(mt, "term.labels")
  ps_args <- list(
    family = gaussian(),
    data = data,
    weights = weights,
    subset = subset,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    pwt = pwt,
    pwt_default_low = pwt_default_low,
    pwt_default_high = pwt_default_high,
    n_prior = n_prior,
    sd = sd,
    dispersion = dispersion,
    intercept_source = intercept_source,
    effects_source = effects_source,
    mu = mu,
    k = k
  )

  setups <- setNames(vector("list", l1), coef_names)
  for (j in seq_len(l1)) {
    f_j <- stats::reformulate(termlabels, response = coef_names[j])
    setups[[j]] <- do.call(
      Prior_Setup,
      c(list(formula = f_j), ps_args, list(...))
    )
  }

  attr(setups, "call") <- call
  attr(setups, "formula") <- formula
  class(setups) <- c("multi_PriorSetup", "list")
  setups
}
