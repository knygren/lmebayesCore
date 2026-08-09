#' @rdname pfamily_list
#' @order 3
#' @details
#' For a \code{\link{Prior_SetupGroup}} object, each \strong{row block}
#' (named as in \code{names(object)}) is mapped to one
#' \code{\link[glmbayesCore]{pfamily}} for use with
#' \code{lmbBlock} / \code{glmbBlock} (or any API that takes a per-block
#' \code{pfamily_list}):
#' \itemize{
#'   \item \code{"dNormal"} (default): \code{dNormal(mu, Sigma, dispersion)}
#'     from that block's \code{\link[glmbayesCore]{Prior_Setup}} result.
#'     Available for all families supported by \code{Prior_SetupGroup}.
#'   \item \code{"dNormal_Gamma"}: \code{dNormal_Gamma(mu, Sigma_0, shape, rate)}
#'     (Gaussian only; uses conjugate Normal--Gamma hyperparameters).
#'   \item \code{"dIndependent_Normal_Gamma"}:
#'     \code{dIndependent_Normal_Gamma(mu, Sigma, shape_ING, rate)}
#'     (Gaussian only; uses \code{shape_ING} from \code{Prior_Setup}).
#' }
#'
#' When \code{ptypes} is omitted (\code{NULL}), every block defaults to
#' \code{"dNormal"}.
#'
#' @export
#' @method pfamily_list Prior_SetupGroup
pfamily_list.Prior_SetupGroup <- function(object, ptypes = NULL, ...) {
  allowed <- c("dNormal", "dNormal_Gamma", "dIndependent_Normal_Gamma")
  block_ids <- names(object)
  if (is.null(block_ids) || any(!nzchar(block_ids))) {
    stop("'object' must be a named Prior_SetupGroup list.", call. = FALSE)
  }
  k <- length(block_ids)

  if (is.null(ptypes)) {
    ptypes <- "dNormal"
  }

  if (is.list(ptypes)) {
    ok <- vapply(
      ptypes,
      function(p) is.character(p) && length(p) == 1L && !is.na(p),
      logical(1L)
    )
    if (!all(ok)) {
      stop("'ptypes' list elements must each be a single string.",
           call. = FALSE)
    }
    nms <- names(ptypes)
    ptypes <- vapply(ptypes, identity, character(1L))
    names(ptypes) <- nms
  }
  if (!is.character(ptypes) || length(ptypes) < 1L || anyNA(ptypes)) {
    stop("'ptypes' must be a character vector or list of strings.",
         call. = FALSE)
  }
  bad <- setdiff(unique(ptypes), allowed)
  if (length(bad) > 0L) {
    stop(
      "Invalid 'ptypes' value(s): ", paste(bad, collapse = ", "),
      ". Allowed: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (length(ptypes) == 1L) {
    ptypes <- stats::setNames(rep(unname(ptypes), k), block_ids)
  } else {
    if (length(ptypes) != k) {
      stop(
        sprintf(
          "'ptypes' has length %d but Prior_SetupGroup has %d block(s): %s.",
          length(ptypes), k, paste(block_ids, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (!is.null(names(ptypes)) && any(nzchar(names(ptypes)))) {
      if (!setequal(names(ptypes), block_ids)) {
        stop(
          "Names of 'ptypes' must match the block IDs: ",
          paste(block_ids, collapse = ", "), ".",
          call. = FALSE
        )
      }
      ptypes <- ptypes[block_ids]
    } else {
      names(ptypes) <- block_ids
    }
  }

  out <- stats::setNames(vector("list", k), block_ids)
  for (id in block_ids) {
    ps <- object[[id]]
    out[[id]] <- switch(
      ptypes[[id]],
      dNormal = {
        disp <- ps$dispersion
        if (is.null(disp)) {
          disp <- 1
        }
        glmbayesCore::dNormal(
          mu         = ps$mu,
          Sigma      = ps$Sigma,
          dispersion = disp
        )
      },
      dNormal_Gamma = {
        if (is.null(ps$Sigma_0) || is.null(ps$shape) || is.null(ps$rate)) {
          stop(
            "Block '", id, "': ptypes = \"dNormal_Gamma\" requires ",
            "Sigma_0, shape, and rate on the Prior_Setup result ",
            "(family = gaussian()).",
            call. = FALSE
          )
        }
        glmbayesCore::dNormal_Gamma(
          mu      = ps$mu,
          Sigma_0 = ps$Sigma_0,
          shape   = ps$shape,
          rate    = ps$rate
        )
      },
      dIndependent_Normal_Gamma = {
        if (is.null(ps$Sigma) || is.null(ps$shape_ING) || is.null(ps$rate)) {
          stop(
            "Block '", id, "': ptypes = \"dIndependent_Normal_Gamma\" requires ",
            "Sigma, shape_ING, and rate on the Prior_Setup result ",
            "(family = gaussian()).",
            call. = FALSE
          )
        }
        glmbayesCore::dIndependent_Normal_Gamma(
          mu    = ps$mu,
          Sigma = ps$Sigma,
          shape = ps$shape_ING,
          rate  = ps$rate
        )
      }
    )
  }

  attr(out, "ptypes") <- ptypes
  class(out) <- c("pfamily_list", "list")
  out
}
