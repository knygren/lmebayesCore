## Safe-set restricted Block~1 group draws (widetilde B rejection).

#' @noRd
.c05_prior_list_from_fixef <- function(design,
                                       fixef,
                                       p11,
                                       measurement_prior_list,
                                       group_levels) {
  mu_all <- as.matrix(
    build_mu_all(design, fixef, group_levels = group_levels)$mu_all
  )
  list(
    mu = mu_all,
    P = p11$P_b,
    dispersion = measurement_prior_list$group.dispersion,
    ddef = is.null(measurement_prior_list$group.dispersion)
  )
}

#' Adaptive rejection batch for joint safe-set draws.
#' @noRd
.simfunc_group_safe_rejection <- function(draw_once,
                                          beta_set,
                                          n_target,
                                          max_tries = NULL,
                                          warn_prefix = "simfunc_group_safe") {
  if (is.null(max_tries)) {
    max_tries <- max(1000000L, n_target * 10000L)
  }
  max_tries <- as.integer(max_tries)

  accepted <- vector("list", n_target)
  n_accepted <- 0L
  n_tried <- 0L
  accept_rate <- NA_real_
  rate_floor <- 0.01
  last_out <- NULL

  while (n_accepted < n_target && n_tried < max_tries) {
    need <- n_target - n_accepted
    if (!is.finite(accept_rate) || accept_rate <= 0) {
      batch <- max(need, 100L)
    } else {
      batch <- max(
        need,
        as.integer(ceiling(need / max(accept_rate, rate_floor))),
        10L
      )
    }
    batch <- min(batch, max_tries - n_tried)

    for (m in seq_len(batch)) {
      n_tried <- n_tried + 1L
      out <- draw_once()
      b_draw <- out$coefficients
      if (beta_in_marginal_safe_set(b_draw, beta_set)) {
        n_accepted <- n_accepted + 1L
        accepted[[n_accepted]] <- b_draw
        last_out <- out
        if (n_accepted >= n_target) {
          break
        }
      }
    }

    accept_rate <- n_accepted / n_tried
  }

  if (n_accepted < 1L) {
    stop(
      "No joint draws accepted in widetilde B after ", n_tried,
      " attempts; check beta_set calibration or increase max_tries.",
      call. = FALSE
    )
  }

  if (n_accepted < n_target) {
    warning(
      warn_prefix, "(): only ", n_accepted, " of ", n_target,
      " target draws accepted (", n_tried, " attempts, rate = ",
      signif(accept_rate, 3L), ").",
      call. = FALSE
    )
  }

  coef_mat <- accepted[[1L]]
  J <- nrow(coef_mat)
  p_re <- ncol(coef_mat)
  group_levels <- rownames(coef_mat)
  re_names <- colnames(coef_mat)
  acc_arr <- array(
    0,
    dim = c(n_accepted, J, p_re),
    dimnames = list(NULL, group_levels, re_names)
  )
  for (k in seq_len(n_accepted)) {
    acc_arr[k, , ] <- accepted[[k]]
  }

  list(
    last_out = last_out,
    coefficients_all = acc_arr,
    b_mean = apply(acc_arr, c(2L, 3L), mean),
    n_accepted = n_accepted,
    n_target = n_target,
    n_tried = n_tried,
    accept_rate = accept_rate
  )
}

#' @noRd
.simfunc_group_safe_finalize <- function(out,
                                         beta_set,
                                         rejection = NULL,
                                         safe_class) {
  if (!is.null(rejection)) {
    out$coefficients <- rejection$last_out$coefficients
    out$coef.mode <- rejection$last_out$coef.mode
    out$coefficients_all <- rejection$coefficients_all
    out$b_mean <- rejection$b_mean
    out$n_accepted <- rejection$n_accepted
    out$n_target <- rejection$n_target
    out$n_tried <- rejection$n_tried
    out$accept_rate <- rejection$accept_rate
    out$restricted <- TRUE
  } else {
    out$n_accepted <- out$n
    out$n_target <- out$n
    out$n_tried <- out$n
    out$accept_rate <- NA_real_
    out$restricted <- FALSE
    out$coefficients_all <- NULL
    out$b_mean <- out$coefficients
  }
  out$beta_set <- beta_set
  class(out) <- c(safe_class, class(out))
  out
}

#' Test whether group random effects lie in widetilde B(delta_2).
#'
#' @param b Group effects as a \code{J x p_re} matrix (rows = groups) or a
#'   stacked numeric vector in \code{\link{group_precision_floor}} order.
#' @param beta_set A \code{\link{beta_marginal_safe_set}} object.
#' @return Logical scalar: \code{TRUE} if \eqn{\Xi(\beta) \le r_{\mathrm{Gauss}}}.
#' @seealso \code{\link{rNormal_reg_group_safe}},
#'   \code{\link{rNormalGLM_reg_group_safe}}
#' @export
beta_in_marginal_safe_set <- function(b, beta_set) {
  if (!inherits(beta_set, "beta_marginal_safe_set")) {
    stop("'beta_set' must be from beta_marginal_safe_set().", call. = FALSE)
  }
  prep <- beta_set$beta_mode$prep
  if (is.null(prep)) {
    stop(
      "'beta_set$beta_mode' is missing 'prep'; recompute beta_marginal_mode().",
      call. = FALSE
    )
  }

  if (is.matrix(b)) {
    b_mat <- b
    if (nrow(b_mat) != prep$J || ncol(b_mat) != prep$p_re) {
      stop(
        "'b' matrix must be ", prep$J, " x ", prep$p_re, ".",
        call. = FALSE
      )
    }
  } else {
    b <- as.numeric(b)
    if (length(b) != prep$J * prep$p_re) {
      stop(
        "Stacked 'b' must have length ", prep$J * prep$p_re, ".",
        call. = FALSE
      )
    }
    b_mat <- matrix(
      b,
      nrow = prep$J,
      ncol = prep$p_re,
      byrow = TRUE,
      dimnames = list(prep$group_levels, prep$re_names)
    )
  }

  engine <- beta_set$beta_mode$engine
  if (is.null(engine) || !is.function(engine$Xi)) {
    stop(
      "'beta_set$beta_mode' must contain an 'engine' with Xi().",
      call. = FALSE
    )
  }
  r_level <- beta_set$level$r_gauss
  if (is.null(r_level) || !is.finite(r_level)) {
    stop("'beta_set' is missing level$r_gauss.", call. = FALSE)
  }

  beta_vec <- .group_floor_stack_beta(b_mat)
  engine$Xi(beta_vec) <= r_level
}

#' Gaussian blockwise draw with widetilde-B rejection
#'
#' Like \code{\link{rNormal_reg_group}}, but when \code{beta_set} is supplied
#' draws independent unrestricted joint samples via \code{rNormal_reg_group},
#' stacks group coefficients, and keeps only draws with
#' \code{\link{beta_in_marginal_safe_set}}. Additional batches use the
#' observed acceptance rate until \code{n} accepted joint draws are collected.
#'
#' @inheritParams rNormal_reg_group
#' @param beta_set Optional \code{\link{beta_marginal_safe_set}}; when
#'   \code{NULL}, behaves as \code{\link{rNormal_reg_group}}.
#' @param max_tries Maximum attempts when \code{beta_set} is set.
#' @param mc_seed Optional RNG seed.
#' @return Like \code{\link{rNormal_reg_group}}, plus safe-set fields
#'   \code{n_accepted}, \code{n_target}, \code{n_tried}, \code{accept_rate},
#'   \code{restricted}, \code{beta_set}, and \code{coefficients_all} when
#'   \code{n > 1} accepted draws.
#' @seealso \code{\link{rNormalGLM_reg_group_safe}},
#'   \code{\link{beta_in_marginal_safe_set}}
#' @export
rNormal_reg_group_safe <- function(n,
                                   y,
                                   x,
                                   group,
                                   prior_list = NULL,
                                   prior_lists = NULL,
                                   beta_set = NULL,
                                   offset = NULL,
                                   weights = 1,
                                   max_tries = NULL,
                                   mc_seed = NULL,
                                   Gridtype = 2L) {
  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  if (!is.null(beta_set) && !inherits(beta_set, "beta_marginal_safe_set")) {
    stop("'beta_set' must be from beta_marginal_safe_set().", call. = FALSE)
  }
  if (!is.null(mc_seed)) {
    set.seed(mc_seed)
  }

  draw_args <- list(
    n = 1L,
    y = y,
    x = x,
    group = group,
    prior_list = prior_list,
    prior_lists = prior_lists,
    offset = offset,
    weights = weights,
    Gridtype = Gridtype
  )

  if (is.null(beta_set)) {
    draw_args$n <- n
    out <- do.call(rNormal_reg_group, draw_args)
    return(.simfunc_group_safe_finalize(
      out,
      beta_set = NULL,
      rejection = NULL,
      safe_class = "rNormal_reg_group_safe"
    ))
  }

  draw_once <- function() {
    do.call(rNormal_reg_group, draw_args)
  }
  rejection <- .simfunc_group_safe_rejection(
    draw_once = draw_once,
    beta_set = beta_set,
    n_target = n,
    max_tries = max_tries,
    warn_prefix = "rNormal_reg_group_safe"
  )
  .simfunc_group_safe_finalize(
    rejection$last_out,
    beta_set = beta_set,
    rejection = rejection,
    safe_class = "rNormal_reg_group_safe"
  )
}

#' GLM blockwise draw with widetilde-B rejection
#'
#' Non-Gaussian counterpart of \code{\link{rNormal_reg_group_safe}} using
#' \code{\link{rNormalGLM_reg_group}} for each unrestricted joint attempt,
#' then joint rejection on \code{beta_set}.
#'
#' @inheritParams rNormalGLM_reg_group
#' @inheritParams rNormal_reg_group_safe
#' @return Like \code{\link{rNormalGLM_reg_group}}, plus safe-set metadata
#'   (see \code{\link{rNormal_reg_group_safe}}).
#' @seealso \code{\link{rNormal_reg_group_safe}},
#'   \code{\link{beta_in_marginal_safe_set}}
#' @export
rNormalGLM_reg_group_safe <- function(n,
                                      y,
                                      x,
                                      group,
                                      prior_list = NULL,
                                      prior_lists = NULL,
                                      family,
                                      beta_set = NULL,
                                      offset = NULL,
                                      weights = 1,
                                      max_tries = NULL,
                                      mc_seed = NULL,
                                      Gridtype = 2L,
                                      n_envopt = NULL,
                                      use_parallel = TRUE,
                                      use_opencl = FALSE,
                                      verbose = FALSE,
                                      progbar = FALSE) {
  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  if (!is.null(beta_set) && !inherits(beta_set, "beta_marginal_safe_set")) {
    stop("'beta_set' must be from beta_marginal_safe_set().", call. = FALSE)
  }
  if (!is.null(mc_seed)) {
    set.seed(mc_seed)
  }

  draw_args <- list(
    n = 1L,
    y = y,
    x = x,
    group = group,
    prior_list = prior_list,
    prior_lists = prior_lists,
    family = family,
    offset = offset,
    weights = weights,
    Gridtype = Gridtype,
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose,
    progbar = progbar
  )

  if (is.null(beta_set)) {
    draw_args$n <- n
    out <- do.call(rNormalGLM_reg_group, draw_args)
    return(.simfunc_group_safe_finalize(
      out,
      beta_set = NULL,
      rejection = NULL,
      safe_class = "rNormalGLM_reg_group_safe"
    ))
  }

  draw_once <- function() {
    do.call(rNormalGLM_reg_group, draw_args)
  }
  rejection <- .simfunc_group_safe_rejection(
    draw_once = draw_once,
    beta_set = beta_set,
    n_target = n,
    max_tries = max_tries,
    warn_prefix = "rNormalGLM_reg_group_safe"
  )
  .simfunc_group_safe_finalize(
    rejection$last_out,
    beta_set = beta_set,
    rejection = rejection,
    safe_class = "rNormalGLM_reg_group_safe"
  )
}

#' Route a safe Block~1 draw by \code{family} (Gaussian vs GLM).
#' @noRd
.c05_safe_group_draw <- function(design,
                                 fixef,
                                 p11,
                                 measurement_prior_list,
                                 family,
                                 beta_set,
                                 n = 1L,
                                 max_tries = NULL,
                                 mc_seed = NULL) {
  group_levels <- p11$group_levels
  prior_list <- .c05_prior_list_from_fixef(
    design = design,
    fixef = fixef,
    p11 = p11,
    measurement_prior_list = measurement_prior_list,
    group_levels = group_levels
  )
  offset <- if (!is.null(design$offset)) design$offset else NULL
  weights <- if (!is.null(design$weights)) design$weights else 1

  common <- list(
    n = n,
    y = design$y,
    x = design$D,
    group = design$group,
    prior_list = prior_list,
    beta_set = beta_set,
    offset = offset,
    weights = weights,
    max_tries = max_tries,
    mc_seed = mc_seed
  )

  if (identical(family$family, "gaussian")) {
    if (is.null(measurement_prior_list$group.dispersion)) {
      stop(
        "'measurement_prior_list$group.dispersion' is required for gaussian().",
        call. = FALSE
      )
    }
    do.call(rNormal_reg_group_safe, common)
  } else {
    common$family <- family
    do.call(rNormalGLM_reg_group_safe, common)
  }
}

#' @noRd
.c05_safe_group_draw_to_estep <- function(safe_out, group_levels) {
  J <- length(group_levels)
  p_re <- ncol(safe_out$coefficients)
  re_names <- colnames(safe_out$coefficients)

  if (!is.null(safe_out$coefficients_all)) {
    acc_arr <- safe_out$coefficients_all
    n_accepted <- safe_out$n_accepted
  } else {
    n_accepted <- 1L
    acc_arr <- array(
      safe_out$coefficients,
      dim = c(1L, J, p_re),
      dimnames = list(NULL, group_levels, re_names)
    )
  }

  b_mean <- safe_out$b_mean
  V_list <- stats::setNames(vector("list", J), group_levels)
  b_mc_se <- matrix(
    0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )
  for (jj in seq_len(J)) {
    lev <- group_levels[jj]
    draws_j <- acc_arr[, jj, , drop = TRUE]
    if (n_accepted >= 2L) {
      V_list[[lev]] <- stats::cov(draws_j)
      b_mc_se[jj, ] <- apply(draws_j, 2L, stats::sd) / sqrt(n_accepted)
    } else {
      V_list[[lev]] <- matrix(0, p_re, p_re)
    }
  }

  list(
    b_mean = b_mean,
    b_mc_se = b_mc_se,
    V_list = V_list,
    estep = "mc",
    n = n_accepted,
    n_target = safe_out$n_target,
    n_tried = safe_out$n_tried,
    accept_rate = safe_out$accept_rate,
    beta_set = safe_out$beta_set,
    restricted = isTRUE(safe_out$restricted)
  )
}
