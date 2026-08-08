#' Prior setup for row-block regressions
#'
#' Runs \code{\link[glmbayesCore]{Prior_Setup}} on each block subset of the data.
#' Typical callers are row-block engines such as
#' \code{\link{block_rNormalReg}} / \code{\link{block_rNormalGLM}}, or the
#' higher-level \code{lmbBlock} / \code{glmbBlock} wrappers in \pkg{lmebayes}.
#'
#' @param formula A \code{\link{formula}} with a single response.
#' @param block Block partition: \code{factor} or vector of length \code{nrow(data)}
#'   (after \code{model.frame}), a column name in \code{data}, \code{l2_blocks}
#'   counts, or a list of row index vectors (see \code{\link{normalize_block}}).
#' @param pwt,n_prior,sd,mu,dispersion,k Each of these six calibration
#'   arguments (passed to \code{\link[glmbayesCore]{Prior_Setup}} for every
#'   block) may be supplied in \strong{either} of two forms:
#'   \itemize{
#'     \item the ordinary shared form documented in
#'       \code{\link[glmbayesCore]{Prior_Setup}} (a scalar, or a
#'       per-coefficient vector for \code{pwt}/\code{sd}/\code{mu}) -- applied
#'       identically to every block (default, unchanged behavior); or
#'     \item a \strong{named list with one element per block}, keyed by the
#'       block/group level IDs (\code{names(Prior_SetupBlock(...))}, i.e.
#'       \code{normalize_block()$ids}) -- each element supplies that block's
#'       own value, in the same shape \code{Prior_Setup()} would otherwise
#'       accept. \code{names(x)} must exactly match the block IDs (extra or
#'       missing names raise an error). An element may be \code{NULL},
#'       meaning "use \code{Prior_Setup()}'s own default for this argument in
#'       this block."
#'   }
#'   The two forms may be mixed freely across the six arguments in one call
#'   (e.g. a shared \code{sd} together with a per-block \code{n_prior}).
#' @inheritParams glmbayesCore::Prior_Setup
#' @return A named list of class \code{"block_PriorSetup"}. Each element is a
#'   \code{\link[glmbayesCore]{Prior_Setup}} result for one block.
#' @seealso \code{\link[glmbayesCore]{Prior_Setup}},
#'   \code{\link[glmbayesCore]{multi_prior_setup}},
#'   \code{\link{Prior_Setup_GLMM}}, \code{\link{normalize_block}},
#'   \code{\link{block_rNormalReg}}, \code{\link{block_rNormalGLM}}
#' @example inst/examples/Ex_Prior_SetupBlock.R
#' @export
Prior_SetupBlock <- function(
    formula,
    block,
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
  fam_ok <- family$family %in% c("gaussian", "poisson")
  if (is.null(family$family) || !fam_ok) {
    stop(
      "Prior_SetupBlock() supports family = gaussian() or poisson() only.",
      call. = FALSE
    )
  }
  if (missing(data)) {
    data <- environment(formula)
  }

  meta <- .lmebayes_formula_block_meta(
    formula = formula,
    block = block,
    data = data,
    subset = if (!missing(subset)) subset else NULL,
    weights = if (!missing(weights)) weights else NULL,
    na.action = if (!missing(na.action)) na.action else NULL,
    offset = if (!missing(offset)) offset else NULL,
    contrasts = if (!missing(contrasts)) contrasts else NULL
  )
  block_ids <- meta$block_info$ids

  ## Each of these six calibration arguments may be shared (reused for every
  ## block, as before) or a named list keyed by `block_ids` (one value per
  ## block); see .lmebayes_resolve_block_calibration_arg().
  pwt_r        <- .lmebayes_resolve_block_calibration_arg(pwt, block_ids, "pwt")
  n_prior_r    <- .lmebayes_resolve_block_calibration_arg(n_prior, block_ids, "n_prior")
  sd_r         <- .lmebayes_resolve_block_calibration_arg(sd, block_ids, "sd")
  mu_r         <- .lmebayes_resolve_block_calibration_arg(mu, block_ids, "mu")
  dispersion_r <- .lmebayes_resolve_block_calibration_arg(dispersion, block_ids, "dispersion")
  k_r          <- .lmebayes_resolve_block_calibration_arg(k, block_ids, "k")

  ps_args <- list(
    family = family,
    data = data,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    pwt_default_low = pwt_default_low,
    pwt_default_high = pwt_default_high,
    intercept_source = intercept_source,
    effects_source = effects_source
  )

  setups <- vector("list", meta$block_info$k)
  for (b in seq_len(meta$block_info$k)) {
    block_id <- block_ids[[b]]
    rows_b <- .lmebayes_rows_to_data_subset(
      meta$block_info$rows[[b]], meta$mf, data
    )
    setups[[b]] <- tryCatch(
      do.call(
        glmbayesCore::Prior_Setup,
        c(
          list(formula = formula, subset = rows_b),
          ps_args,
          list(
            pwt        = pwt_r[[block_id]],
            n_prior    = n_prior_r[[block_id]],
            sd         = sd_r[[block_id]],
            dispersion = dispersion_r[[block_id]],
            mu         = mu_r[[block_id]],
            k          = k_r[[block_id]]
          ),
          list(...)
        )
      ),
      error = function(e) {
        stop(
          "Prior_SetupBlock(): block '", block_id, "': ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }
  names(setups) <- block_ids

  attr(setups, "call") <- call
  attr(setups, "formula") <- formula
  attr(setups, "block") <- block
  attr(setups, "block_info") <- meta$block_info
  class(setups) <- c("block_PriorSetup", "list")
  setups
}

#' @noRd
.lmebayes_resolve_block_calibration_arg <- function(x, block_ids, arg_name) {
  if (!is.list(x)) {
    return(stats::setNames(rep(list(x), length(block_ids)), block_ids))
  }
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop(
      "Prior_SetupBlock(): '", arg_name, "' is a list but is not fully named; ",
      "supply one named element per block/group level (",
      paste(block_ids, collapse = ", "), ").",
      call. = FALSE
    )
  }
  missing_ids <- setdiff(block_ids, nm)
  extra_ids   <- setdiff(nm, block_ids)
  if (length(missing_ids) || length(extra_ids)) {
    stop(
      "Prior_SetupBlock(): names('", arg_name, "') must exactly match the ",
      "block/group level IDs.",
      if (length(missing_ids)) {
        paste0("\n  Missing: ", paste(missing_ids, collapse = ", "))
      } else {
        ""
      },
      if (length(extra_ids)) {
        paste0("\n  Unexpected: ", paste(extra_ids, collapse = ", "))
      } else {
        ""
      },
      call. = FALSE
    )
  }
  x[block_ids]
}

#' @noRd
.lmebayes_formula_block_meta <- function(
    formula,
    block,
    data,
    subset = NULL,
    weights = NULL,
    na.action = NULL,
    offset = NULL,
    contrasts = NULL
) {
  mf_args <- list(
    formula = formula,
    data = data,
    drop.unused.levels = TRUE
  )
  if (!is.null(subset)) mf_args$subset <- subset
  if (!is.null(weights)) mf_args$weights <- weights
  if (!is.null(na.action)) mf_args$na.action <- na.action
  if (!is.null(offset)) mf_args$offset <- offset
  if (!is.null(contrasts)) mf_args$contrasts <- contrasts
  mf <- do.call(stats::model.frame, mf_args)

  l2 <- nrow(mf)
  block_vec <- .lmebayes_resolve_block(block, data, mf, l2)
  block_info <- normalize_block(block_vec, l2)

  mt <- attr(mf, "terms")
  x_mat <- stats::model.matrix(mt, mf, contrasts)
  p <- ncol(x_mat)
  pred_names <- colnames(x_mat)
  if (is.null(pred_names) || length(pred_names) != p) {
    pred_names <- paste0("X", seq_len(p))
  }

  list(
    mf = mf,
    block_info = block_info,
    p = p,
    pred_names = pred_names
  )
}

#' @noRd
.lmebayes_rows_to_data_subset <- function(rows_mf, mf, data) {
  rn <- rownames(mf)
  if (is.null(rn)) {
    return(rows_mf)
  }
  if (!is.null(rownames(data))) {
    out <- match(rn[rows_mf], rownames(data))
    if (anyNA(out)) {
      stop("Could not map model.frame rows to rownames(data).", call. = FALSE)
    }
    return(out)
  }
  as.integer(rn[rows_mf])
}

#' @noRd
.lmebayes_resolve_block <- function(block, data, mf, l2) {
  if (is.list(block)) {
    return(block)
  }
  if (is.character(block) && length(block) == 1L && block %in% names(data)) {
    rn_mf <- rownames(mf)
    col <- data[[block]]
    if (length(col) == l2 && is.null(rn_mf)) {
      return(col)
    }
    if (!is.null(rn_mf)) {
      if (!is.null(rownames(data))) {
        return(col[match(rn_mf, rownames(data))])
      }
      return(col[as.integer(rn_mf)])
    }
  }
  block <- as.vector(block)
  if (length(block) == l2) {
    return(block)
  }
  stop(
    "'block' must have length nrow(model.frame), be a list of row indices, ",
    "l2_blocks counts, or a single column name in 'data'.",
    call. = FALSE
  )
}
