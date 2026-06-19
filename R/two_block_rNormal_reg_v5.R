#' Two-block Gibbs sampler via independent short chains (development v5)
#'
#' Same short-chain Gibbs semantics as \code{\link{two_block_rNormal_reg_v4}},
#' with per-chain state buffers pre-allocated in C++ via
#' \code{.two_block_rNormal_reg_v5_cpp}.  Loop order is sweep-outer
#' (\code{m} then \code{i}) with per-sweep logging and chain progress bars.
#'
#' @inheritParams two_block_rNormal_reg_v2
#' @param seed_offset Integer added to \code{seed} for chain \code{i}
#'   (\code{seed + seed_offset + i + 1} in C++). Default \code{0L}.
#' @param collect_block1 Logical. If \code{TRUE}, row-bind Block~1
#'   (\code{coefficients}) draws from every chain.  Default \code{TRUE}.
#' @return Object of class \code{c("two_block_rNormal_reg_v5",
#'   "two_block_rNormal_reg_v2", "two_block_rNormal_reg")}.  Same fields as
#'   \code{\link{two_block_rNormal_reg_v2}}.
#' @family simfuncs
#' @seealso \code{\link{two_block_rNormal_reg_v4}}, \code{\link{rGLMM}}
#' @export
two_block_rNormal_reg_v5 <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list_block1,
    pfamily_list,
    fixef_start,
    re_coef_names = colnames(x),
    group_levels = levels(block),
    group_name = NULL,
    m_convergence = 10L,
    sampling = c("replicate", "chain"),
    family = gaussian(),
    offset = NULL,
    weights = 1,
    Gridtype = 2L,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    seed = NULL,
    seed_offset = 0L,
    collect_block1 = TRUE,
    progbar = FALSE,
    stage_label = "",
    diag_sweeps = FALSE,
    fixef_mode = NULL,
    b_mode = NULL
) {

  cl <- match.call()
  sampling <- match.arg(sampling)
  if (!identical(sampling, "replicate")) {
    stop("Only sampling = \"replicate\" is implemented.", call. = FALSE)
  }

  family <- .two_block_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")

  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  m_convergence <- as.integer(m_convergence[1L])
  if (m_convergence < 1L) {
    stop("'m_convergence' must be at least 1.", call. = FALSE)
  }

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  if (is.null(re_coef_names) || length(re_coef_names) != ncol(x)) {
    re_coef_names <- if (ncol(x) >= 1L) {
      cn <- colnames(x)
      if (is.null(cn) || length(cn) != ncol(x)) paste0("RE", seq_len(ncol(x))) else cn
    } else {
      stop("'x' must have at least one column.", call. = FALSE)
    }
  }
  colnames(x) <- re_coef_names
  re_names <- re_coef_names

  group_levels <- as.character(group_levels)
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    group_name <- tryCatch(
      deparse(substitute(block))[1L],
      error = function(e) "group"
    )
    if (!nzchar(group_name)) group_name <- "group"
  }

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop(
      "length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
      call. = FALSE
    )
  }
  if (!setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, re_names,
    J = length(group_levels)
  )

  if (!is.list(fixef_start) || is.null(names(fixef_start))) {
    stop("'fixef_start' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(fixef_start), re_names)) {
    stop("names(fixef_start) must match re_coef_names.", call. = FALSE)
  }
  fixef_start <- fixef_start[re_names]

  block1_prior_meta <- .two_block_validate_block1_prior(
    prior_list_block1,
    family = family
  )

  offset2 <- offset
  wt <- weights
  if (is.null(offset2)) {
    offset2 <- rep(0, l2)
  } else {
    offset2 <- as.numeric(offset2)
    if (length(offset2) == 1L) offset2 <- rep(offset2, l2)
    if (length(offset2) != l2) {
      stop("length(offset) must be 1 or length(y).", call. = FALSE)
    }
  }
  if (length(wt) == 1L) wt <- rep(wt, l2)
  if (length(wt) != l2) {
    stop("length(weights) must be 1 or length(y).", call. = FALSE)
  }

  famfunc_block1 <- glmbfamfunc(if (is_gaussian) gaussian() else family)
  famfunc_gauss <- glmbfamfunc(gaussian())
  n_envopt_use <- if (is.null(n_envopt)) 1L else as.integer(n_envopt)

  x_hyper_mats <- lapply(x_hyper, as.matrix)

  cpp_out <- .two_block_rNormal_reg_v5_cpp(
    n                 = n,
    m_convergence     = m_convergence,
    y                 = y,
    x                 = x,
    block             = block,
    x_hyper           = x_hyper_mats,
    prior_list_block1 = prior_list_block1,
    dispersion_block1 = block1_prior_meta$dispersion,
    ddef_block1       = block1_prior_meta$ddef,
    pfamily_list      = pfamily_list,
    fixef_start       = fixef_start,
    group_levels      = group_levels,
    family            = family$family,
    link              = family$link,
    f2                = famfunc_block1$f2,
    f3                = famfunc_block1$f3,
    f2_gauss          = famfunc_gauss$f2,
    f3_gauss          = famfunc_gauss$f3,
    offset            = offset2,
    wt                = wt,
    Gridtype          = as.integer(Gridtype),
    n_envopt          = n_envopt_use,
    use_parallel      = use_parallel,
    use_opencl        = use_opencl,
    verbose           = verbose,
    seed              = seed,
    seed_offset       = as.integer(seed_offset),
    progbar           = isTRUE(progbar),
    stage_label       = as.character(stage_label)[1L],
    diag_sweeps       = isTRUE(diag_sweeps),
    fixef_mode        = fixef_mode,
    b_mode            = b_mode
  )

  res <- .two_block_format_v5_cpp_out(
    cpp_out         = cpp_out,
    n               = n,
    re_names        = re_names,
    fixef_start     = fixef_start,
    group_levels    = group_levels,
    group_name      = group_name,
    pfamily_list    = pfamily_list,
    family          = family,
    m_convergence   = m_convergence,
    sampling        = sampling,
    cl              = cl
  )

  if (!isTRUE(collect_block1)) {
    res$coefficients <- NULL
  }

  structure(
    res,
    class = c("two_block_rNormal_reg_v5", class(res))
  )
}

#' Format raw v5 C++ output into a two_block_rNormal_reg object
#' @noRd
.two_block_format_v5_cpp_out <- function(
    cpp_out,
    n,
    re_names,
    fixef_start,
    group_levels,
    group_name,
    pfamily_list,
    family,
    m_convergence,
    sampling,
    cl
) {
  J <- length(group_levels)
  p_re <- length(re_names)
  group_ids <- as.character(cpp_out$group_ids)

  fixef_draws <- stats::setNames(cpp_out$fixef_draws, re_names)
  for (k in re_names) {
    dimnames(fixef_draws[[k]]) <- list(NULL, names(fixef_start[[k]]))
  }

  fixef <- stats::setNames(cpp_out$fixef_last, re_names)

  b_arr <- array(as.numeric(cpp_out$b_draws), dim = c(J, p_re, n))
  b_i <- matrix(b_arr[, , n], nrow = J, ncol = p_re,
                dimnames = list(group_ids, re_names))

  mu_all <- cpp_out$mu_all_last
  dimnames(mu_all) <- list(re_names, group_levels)

  dispersion_fixef_draws <- cpp_out$dispersion_fixef_draws
  dimnames(dispersion_fixef_draws) <- list(NULL, re_names)

  iters_fixef_draws <- cpp_out$iters_fixef_draws
  dimnames(iters_fixef_draws) <- list(NULL, re_names)

  coef_cols <- c("draw", group_name, re_names)
  draw_rows <- vector("list", n)
  for (i in seq_len(n)) {
    draw_df <- data.frame(
      draw = rep(i, J),
      stringsAsFactors = FALSE
    )
    draw_df[[group_name]] <- group_ids
    for (jj in seq_len(p_re)) {
      draw_df[[re_names[jj]]] <- b_arr[, jj, i]
    }
    draw_rows[[i]] <- draw_df
  }

  coefficients <- do.call(rbind, draw_rows)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, coef_cols, drop = FALSE]

  structure(
    list(
      fixef_draws            = fixef_draws,
      coefficients           = coefficients,
      fixef_last             = fixef,
      b_last                 = b_i,
      mu_all_last            = mu_all,
      dispersion_fixef_draws = dispersion_fixef_draws,
      iters_fixef_draws      = iters_fixef_draws,
      pfamily_list           = pfamily_list,
      family                 = family,
      n                      = n,
      m_convergence          = m_convergence,
      sampling               = sampling,
      fixef_start            = fixef_start,
      re_coef_names          = re_names,
      group_levels           = group_levels,
      group_name             = group_name,
      call                   = cl
    ),
    class = c("two_block_rNormal_reg_v2", "two_block_rNormal_reg")
  )
}
