## Chapter C05 restricted Gibbs minorization: shared assembly and validation.

#' @noRd
.c05_validate <- function(design,
                                        pfamily_list,
                                        family,
                                        dispprior_list = NULL,
                                        fn_name = "population_mode") {
  .lmerb_validate_design(design)
  if (is.null(design$y) || is.null(design$D)) {
    stop("'", fn_name, "()' requires 'design$y' and 'design$D'.", call. = FALSE)
  }

  group_levels <- levels(design$group)
  J <- length(group_levels)
  re_names <- design$groupef.names
  p_re <- length(re_names)

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, re_names, J = J
  )

  prior_pack <- .priors_from_pfamily_list(
    pfamily_list = pfamily_list,
    group.dispersion = .lmebayes_dispprior_list_as_group_dispersion(dispprior_list),
    design = design,
    family = family,
    fn_name = fn_name
  )

  if (isTRUE(prior_pack$any_non_normal)) {
    stop(
      fn_name, "() requires fixed dNormal tau^2 in 'pfamily_list' ",
      "(dIndependent_Normal_Gamma is not supported in v1).",
      call. = FALSE
    )
  }

  measurement_prior_list <- list(
    group.Sigma = prior_pack$group.Sigma,
    pop.prior_list = prior_pack$pop.prior_list
  )
  if (!is.null(prior_pack$group.dispersion)) {
    measurement_prior_list$group.dispersion <- prior_pack$group.dispersion
  }

  list(
    prior_pack = prior_pack,
    measurement_prior_list = measurement_prior_list,
    group_levels = group_levels,
    re_names = re_names,
    J = J,
    p_re = p_re
  )
}

#' Build C05 refresh precision P11 and hyper-design maps H_j.
#' @noRd
.c05_p11 <- function(design,
                     prior_pack,
                     group_levels,
                     measurement_prior_list = NULL,
                     family = gaussian()) {
  re_names <- design$groupef.names
  p_re <- length(re_names)
  J <- length(group_levels)
  pop <- prior_pack$pop.prior_list

  P_b <- solve(as.matrix(prior_pack$group.Sigma))
  P_b <- 0.5 * (P_b + t(P_b))

  X_hyper <- stats::setNames(
    lapply(re_names, function(k) {
      X_k <- as.matrix(design$W[[k]])
      rn <- rownames(X_k)
      if (!is.null(rn)) {
        if (!all(group_levels %in% rn)) {
          stop("design$W[[\"", k, "\"]] is missing rows for some group levels.",
               call. = FALSE)
        }
        X_k <- X_k[group_levels, , drop = FALSE]
      } else if (nrow(X_k) != J) {
        stop("nrow(design$W[[\"", k, "\"]]) must equal the number of groups.",
             call. = FALSE)
      }
      X_k
    }),
    re_names
  )

  q_k <- vapply(X_hyper, ncol, integer(1L))
  q <- sum(q_k)
  gamma_cols <- split(seq_len(q), rep(seq_len(p_re), q_k))

  P_gamma <- stats::setNames(
    lapply(re_names, function(k) {
      solve(as.matrix(pop[[k]]$Sigma))
    }),
    re_names
  )
  mu_gamma <- stats::setNames(
    lapply(re_names, function(k) pop[[k]]$mu),
    re_names
  )
  tau2 <- stats::setNames(
    lapply(re_names, function(k) pop[[k]]$dispersion),
    re_names
  )

  Lambda_gamma <- matrix(0, q, q)
  mu_0 <- numeric(q)
  for (k in seq_len(p_re)) {
    Lambda_gamma[gamma_cols[[k]], gamma_cols[[k]]] <- P_gamma[[re_names[[k]]]]
    mu_0[gamma_cols[[k]]] <- mu_gamma[[re_names[[k]]]]
  }

  P11_RE <- matrix(0, q, q)
  H_list <- vector("list", J)
  for (j in seq_len(J)) {
    H_j <- matrix(0, p_re, q)
    for (k in seq_len(p_re)) {
      H_j[k, gamma_cols[[k]]] <- X_hyper[[re_names[[k]]]][j, ]
    }
    H_list[[j]] <- H_j
    P11_RE <- P11_RE + t(H_j) %*% P_b %*% H_j
  }

  lmerb_system <- NULL
  if (identical(family$family, "gaussian")) {
    if (is.null(measurement_prior_list)) {
      stop(
        "measurement_prior_list is required to build the Gaussian ",
        "posterior system for P11.",
        call. = FALSE
      )
    }
    lmerb_system <- .lmerb_posterior_normal_system(
      design, measurement_prior_list
    )
    Lambda_gamma <- 0.5 * (Lambda_gamma + t(Lambda_gamma))
    P11_RE <- 0.5 * (P11_RE + t(P11_RE))
    P11 <- Lambda_gamma + P11_RE
    P11 <- 0.5 * (P11 + t(P11))
  } else {
    P11 <- Lambda_gamma + P11_RE
    P11 <- 0.5 * (P11 + t(P11))
    P11_RE <- 0.5 * (P11_RE + t(P11_RE))
    Lambda_gamma <- 0.5 * (Lambda_gamma + t(Lambda_gamma))
  }

  chol_P11 <- tryCatch(chol(P11), error = function(e) {
    stop(
      "C05 hyper-design rank deficiency: P11 is not positive definite ",
      "(check design$W and population priors).",
      call. = FALSE
    )
  })

  list(
    P_b = P_b,
    P11 = P11,
    P11_RE = P11_RE,
    Lambda_gamma = Lambda_gamma,
    mu_0 = mu_0,
    P_gamma = P_gamma,
    mu_gamma = mu_gamma,
    tau2 = tau2,
    H_list = stats::setNames(H_list, group_levels),
    X_hyper = X_hyper,
    gamma_cols = gamma_cols,
    q_k = q_k,
    q = q,
    re_names = re_names,
    p_re = p_re,
    J = J,
    group_levels = group_levels,
    lmerb_system = lmerb_system,
    chol_P11 = chol_P11,
    Sigma_star = chol2inv(chol_P11)
  )
}

#' @noRd
.c05_gamma_from_fixef <- function(fixef, p11) {
  gamma <- numeric(p11$q)
  for (k in seq_along(p11$re_names)) {
    nm <- p11$re_names[[k]]
    gamma[p11$gamma_cols[[k]]] <- fixef[[nm]]
  }
  gamma
}

#' @noRd
.c05_fixef_from_gamma <- function(gamma, p11) {
  stats::setNames(
    lapply(p11$re_names, function(k) {
      g <- gamma[p11$gamma_cols[[which(p11$re_names == k)]]]
      names(g) <- colnames(p11$X_hyper[[k]])
      g
    }),
    p11$re_names
  )
}

#' @noRd
.c05_post_P_gamma_list <- function(design, measurement_prior_list) {
  re_names <- design$groupef.names
  stats::setNames(
    lapply(re_names, function(k) {
      crossprod(design$W[[k]]) /
        measurement_prior_list$pop.prior_list[[k]]$dispersion +
        solve(measurement_prior_list$pop.prior_list[[k]]$Sigma)
    }),
    re_names
  )
}
