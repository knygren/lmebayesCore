#' @noRd
.two_block_normalize_family <- function(family) {
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    stop("'family' not recognized.", call. = FALSE)
  }

  okfamilies <- c(
    "gaussian", "poisson", "binomial",
    "quasipoisson", "quasibinomial", "Gamma"
  )
  if (!family$family %in% okfamilies) {
    stop(
      "family \"", family$family, "\" is not supported by two-block samplers.",
      call. = FALSE
    )
  }

  oklinks <- switch(
    family$family,
    gaussian = "identity",
    poisson = "log",
    quasipoisson = "log",
    binomial = c("logit", "probit", "cloglog"),
    quasibinomial = c("logit", "probit", "cloglog"),
    Gamma = "log",
    character(0)
  )
  if (!family$link %in% oklinks) {
    stop(
      "link \"", family$link, "\" not available for family \"",
      family$family, "\".",
      call. = FALSE
    )
  }

  family
}

#' @noRd
.two_block_validate_block1_prior <- function(prior_list_block1, family) {
  if (!is.list(prior_list_block1)) {
    stop("'prior_list_block1' must be a list.", call. = FALSE)
  }
  if (is.null(prior_list_block1$P) && is.null(prior_list_block1$Sigma)) {
    stop("prior_list_block1 must contain 'P' or 'Sigma'.", call. = FALSE)
  }

  ddef <- if ("ddef" %in% names(prior_list_block1)) {
    prior_list_block1$ddef
  } else {
    is.null(prior_list_block1$dispersion)
  }

  dispersion <- prior_list_block1$dispersion

  if (identical(family$family, "gaussian")) {
    if (is.null(dispersion)) {
      stop(
        "prior_list_block1 must contain 'dispersion' for gaussian() Block~1.",
        call. = FALSE
      )
    }
    if (isTRUE(ddef)) {
      stop(
        "For gaussian() Block~1, dNormal() requires an explicit dispersion.",
        call. = FALSE
      )
    }
  }

  if (family$family %in% c("gaussian", "Gamma") && isTRUE(ddef)) {
    stop(
      "For gaussian() and Gamma() models, dNormal() requires an explicit dispersion.",
      call. = FALSE
    )
  }

  list(ddef = ddef, dispersion = dispersion)
}
