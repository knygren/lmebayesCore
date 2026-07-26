## Family-specific "quick reject" pre-check for classical GLM MLE existence
## (or, for Gamma, dispersion estimability), evaluated *before* attempting
## the per-group glm() fit.  Returns a single note string, or character(0)
## if the group passes.  @noRd
.lmebayes_glm_estimable_precheck <- function(fam_name, y_j, n_j, p_j) {
  switch(
    fam_name,
    binomial = if (length(unique(y_j)) < 2L) "single outcome level",
    poisson  = if (all(y_j == 0)) {
      "all-zero response (MLE diverges under log link)"
    },
    Gamma = {
      if (any(y_j <= 0)) {
        "non-positive response (invalid for Gamma family)"
      } else if (n_j <= p_j) {
        "zero residual degrees of freedom for dispersion (n_j <= p_j)"
      }
    },
    NULL
  )
}

## Families for which a per-group glm() fit is attempted (mean-model MLE
## existence, plus family-specific pre-checks above).
.lmebayes_glm_fit_families <- c("binomial", "poisson", "Gamma")

#' Per-group classical-GLM MLE / dispersion-estimability check for Block-1
#'
#' For \code{binomial}, \code{poisson}, and \code{Gamma} families, fits
#' \code{glm(y ~ Z_j - 1, family = family)} on each group with algebraically
#' full-rank \code{Z_j} and marks the group estimable when a family-specific
#' pre-check passes (see below) \emph{and} all coefficients and \code{vcov}
#' entries from the fit are finite.  For \code{gaussian}, no fit is
#' attempted -- the coefficient MLE (OLS) always exists given full column
#' rank, but the group-level residual \emph{dispersion} additionally
#' requires strictly more observations than parameters (\code{n_j > p_j});
#' this is checked directly via arithmetic on \code{n_j}/\code{p_j}. Other
#' families return \code{re_estimable = re_rank} (no check).
#'
#' Family-specific pre-checks (before any fit is attempted):
#' \itemize{
#'   \item \code{binomial}: single outcome level (complete/quasi-complete
#'     separation).
#'   \item \code{poisson}: all-zero response (log-link intercept MLE
#'     diverges to \eqn{-\infty}).
#'   \item \code{Gamma}: non-positive response (invalid domain), or
#'     \code{n_j <= p_j} (no residual degrees of freedom for the
#'     dispersion/shape parameter).
#'   \item \code{gaussian}: \code{n_j <= p_j} (no residual degrees of
#'     freedom for the variance).
#' }
#' @noRd
.lmebayes_block_glm_estimable <- function(y, groups, Z, re_rank, family) {
  g_levs <- names(re_rank)
  if (is.null(g_levs)) {
    g_levs <- levels(groups)
    names(re_rank) <- g_levs
  }

  re_estimable <- stats::setNames(rep(FALSE, length(g_levs)), g_levs)
  fam_name <- family$family

  is_fit_family <- fam_name %in% .lmebayes_glm_fit_families
  is_gaussian   <- identical(fam_name, "gaussian")

  if (!is_fit_family && !is_gaussian) {
    re_estimable[re_rank] <- TRUE
    return(list(
      re_estimable = re_estimable,
      re_glm_check = NULL
    ))
  }

  y <- as.numeric(y)
  g_chr <- as.character(groups)
  fr_levs <- g_levs[re_rank]

  rows_out <- vector("list", length(fr_levs))

  for (ii in seq_along(fr_levs)) {
    lev <- fr_levs[ii]
    rows <- which(g_chr == lev)
    y_j  <- y[rows]
    X_j  <- Z[rows, , drop = FALSE]
    n_j  <- length(y_j)
    p_j  <- ncol(X_j)

    estimable <- FALSE
    note      <- character(0)

    if (n_j < 2L) {
      note <- "fewer than 2 observations"
    } else if (is_gaussian) {
      ## OLS coefficient MLE always exists given full column rank; only the
      ## residual degrees of freedom for the dispersion is at issue here --
      ## no fit is attempted.
      if (n_j <= p_j) {
        note <- "zero residual degrees of freedom for dispersion (n_j <= p_j)"
      } else {
        estimable <- TRUE
      }
    } else {
      note <- .lmebayes_glm_estimable_precheck(fam_name, y_j, n_j, p_j)
      if (length(note) == 0L) {
        df_j <- data.frame(y = y_j, X_j, check.names = FALSE)
        fit <- tryCatch(
          suppressWarnings(
            stats::glm(
              y ~ . - 1,
              data    = df_j,
              family  = family,
              control = stats::glm.control(maxit = 50L)
            )
          ),
          error = function(e) e
        )
        if (inherits(fit, "error")) {
          note <- conditionMessage(fit)
        } else {
          cf <- stats::coef(fit)
          if (length(cf) != p_j) {
            note <- sprintf(
              "glm returned %d coefficient(s), expected %d",
              length(cf), p_j
            )
          } else if (!isTRUE(fit$rank == p_j)) {
            note <- sprintf("rank-deficient glm fit (rank %d, expected %d)",
                            fit$rank, p_j)
          } else if (any(is.na(cf))) {
            note <- "NA coefficient(s)"
          } else if (any(!is.finite(cf))) {
            note <- "non-finite coefficient(s) (possible separation)"
          } else {
            V_ok <- tryCatch({
              V <- stats::vcov(fit)
              is.matrix(V) && all(is.finite(V))
            }, error = function(e) FALSE)
            if (!isTRUE(V_ok)) {
              note <- "vcov not finite (possible separation)"
            } else {
              estimable <- TRUE
            }
          }
        }
      }
    }

    re_estimable[[lev]] <- estimable
    rows_out[[ii]] <- data.frame(
      group     = lev,
      n         = n_j,
      p         = p_j,
      re_rank   = TRUE,
      estimable = estimable,
      note      = if (length(note)) paste(note, collapse = "; ") else "",
      stringsAsFactors = FALSE
    )
  }

  re_glm_check <- if (length(rows_out)) {
    do.call(rbind, rows_out)
  } else {
    data.frame(
      group = character(0),
      n = integer(0),
      p = integer(0),
      re_rank = logical(0),
      estimable = logical(0),
      note = character(0),
      stringsAsFactors = FALSE
    )
  }

  list(re_estimable = re_estimable, re_glm_check = re_glm_check)
}

#' Check identifiability and estimability of a single-factor mixed-model design
#'
#' @description
#' Runs the "two-step identifiability assessment" used internally by
#' \code{\link{model_setup}}: a Level~1 (within-group) rank/estimability
#' check on the random-effects design \code{Z}, followed by a Level~2
#' (across-group) rank check on the level-2 hyper-design matrices
#' \code{X_hyper}, restricted to the Level-1-estimable groups. Exposed as a
#' standalone, documented function so it can be run directly on hand-built
#' \code{Z}/\code{X_hyper} matrices (e.g. matrix-level \code{rLMM_reg()}/
#' \code{rGLMM_reg()} workflows) without requiring a full \code{lmer}/
#' \code{glmer} reference fit.
#'
#' @details
#' \strong{Level 1 (within-group).} For each group \eqn{j}, the within-group
#' random-effects design submatrix \eqn{\mathbf{Z}_j} is checked for full
#' column rank (\code{re_rank}). A rank-deficient group has too few distinct
#' observations to estimate all random slopes independently. Algebraically
#' full-rank groups are then additionally checked for estimability
#' (\code{re_estimable}, with per-group detail in \code{re_glm_check}), with
#' the specific check depending on \code{family}:
#' \itemize{
#'   \item \code{binomial()}, \code{poisson()}, \code{Gamma()}: a classical
#'     \code{glm(y ~ Z_j - 1, family)} fit is attempted; a group can be full
#'     column rank yet have no finite MLE under complete or quasi-complete
#'     separation (\code{binomial}), an all-zero count response
#'     (\code{poisson}), or a non-positive response (\code{Gamma}) --
#'     \code{Gamma} groups are additionally required to have \code{n_j > p_j}
#'     so the dispersion/shape parameter is estimable.
#'   \item \code{gaussian()}: no fit is attempted (the OLS coefficient MLE
#'     always exists given full column rank), but the group must still have
#'     \code{n_j > p_j} so the residual variance/dispersion is estimable
#'     (otherwise the group achieves a perfect/saturated fit with zero
#'     residual degrees of freedom).
#'   \item other families: \code{re_estimable} is set equal to \code{re_rank}
#'     (no additional check yet).
#' }
#'
#' \strong{Level 2 (across-group).} Restricting to the \strong{estimable}
#' groups from Level 1 (\code{re_estimable}), each hyper-design matrix
#' \code{X_hyper[[k]]} is checked for full column rank (\code{hyper_rank}).
#' Rank deficiency at this level means the level-2 hyperparameters are not
#' identified by the data, even as the number of estimable groups grows. The
#' scalar \code{rank_ok} is \code{TRUE} only when every \code{hyper_rank}
#' entry is \code{TRUE}.
#'
#' \strong{Inputs must already be consistently named.} This function does
#' \emph{not} infer or coerce missing names -- it validates them up front and
#' errors immediately on a mismatch, rather than risking a silently wrong
#' match deeper in the rank/estimability computation:
#' \itemize{
#'   \item \code{groups} must be a \code{factor} (not coerced from a
#'     character vector).
#'   \item \code{X_hyper} must be a named list, with names covering
#'     \code{re_coef_names} exactly (no missing, empty, or duplicate names).
#'   \item Every \code{X_hyper[[k]]} must already have \code{rownames()} set,
#'     and those row names must equal \code{levels(groups)} (as sets) -- this
#'     is how \code{hyper_rank} restricts to estimable groups by name, not
#'     position.
#'   \item \code{group_name}, used only for error/verbose messages, must be
#'     resolvable from either the \code{group_name} argument or
#'     \code{attr(groups, "group_name")}.
#' }
#'
#' @param y Numeric response vector, length \code{nrow(Z)}.
#' @param Z Numeric (or \code{\link[Matrix]{Matrix}}) random-effects design
#'   matrix, stacked across observations (\code{design$Z} from
#'   \code{\link{extract_re_hyper_matrices}}).
#' @param groups Factor of length \code{nrow(Z)} giving each observation's
#'   group membership. Must be a \code{factor}; \code{levels(groups)} name
#'   \code{re_rank}/\code{re_estimable} in the return value.
#' @param X_hyper Named list of per-RE-coefficient level-2 design matrices
#'   (one row per group level; \code{rownames} must equal
#'   \code{levels(groups)}).
#' @param family A \code{\link[stats]{family}} object (or a string/function
#'   accepted by \code{\link[stats]{family}}); normalized internally via the
#'   same logic \code{\link{model_setup}} uses. Defaults to \code{gaussian()}.
#' @param re_coef_names Character vector of random-effects coefficient names;
#'   defaults to \code{names(X_hyper)}. If supplied explicitly, must match
#'   \code{names(X_hyper)} exactly (as sets).
#' @param group_name Display name of the grouping factor (e.g.
#'   \code{"school_id"}), used only in error and \code{verbose} messages
#'   (not in any computation). Defaults to \code{attr(groups, "group_name")}
#'   (the convention already used by \code{\link{rLMM_reg}}); must be
#'   resolvable to a non-empty string, either this way or by passing it
#'   explicitly.
#' @param verbose If \code{TRUE}, emit \code{message()}s summarizing the
#'   Level 1 and Level 2 outcomes. Default \code{FALSE} (silent), so calling
#'   this from \code{\link{model_setup}} introduces no new console output.
#' @return A list with:
#'   \describe{
#'     \item{\code{re_rank}}{Named logical, one entry per group (names
#'       \code{levels(groups)}): \code{TRUE} if \code{Z_j} is full column
#'       rank.}
#'     \item{\code{re_estimable}}{Named logical, same names as
#'       \code{re_rank}: \code{TRUE} if the group is additionally estimable
#'       (see Details).}
#'     \item{\code{re_glm_check}}{Per-group diagnostic data frame for
#'       \code{binomial()}/\code{poisson()}/\code{Gamma()}/\code{gaussian()},
#'       or \code{NULL} for other families.}
#'     \item{\code{hyper_rank}}{Named logical, one entry per RE coefficient
#'       (names \code{re_coef_names}): \code{TRUE} if \code{X_hyper[[k]]} is
#'       full column rank when restricted to \code{re_estimable} groups.}
#'     \item{\code{hyper_deficient}}{Negation of \code{hyper_rank}.}
#'     \item{\code{rank_ok}}{Scalar \code{TRUE} only when every
#'       \code{hyper_rank} entry is \code{TRUE}.}
#'     \item{\code{group_levels}}{\code{levels(groups)}, i.e. the exact
#'       names used on \code{re_rank}/\code{re_estimable}.}
#'     \item{\code{group_name}}{The resolved grouping-factor display name.}
#'   }
#' @seealso \code{\link{model_setup}}, \code{\link{extract_re_hyper_matrices}}
#' @examples
#' set.seed(1)
#' J <- 8L; n_j <- 6L
#' groups <- factor(rep(paste0("g", seq_len(J)), each = n_j))
#' x1 <- rnorm(J * n_j)
#' Z <- cbind("(Intercept)" = 1, x1 = x1)
#' y <- 1 + x1 + rnorm(J * n_j)
#' X_hyper <- list(
#'   "(Intercept)" = matrix(
#'     1, J, 1, dimnames = list(levels(groups), "(Intercept)")
#'   ),
#'   x1 = matrix(1, J, 1, dimnames = list(levels(groups), "(Intercept)"))
#' )
#' ident <- check_identifiability(
#'   y = y, Z = Z, groups = groups, X_hyper = X_hyper,
#'   family = gaussian(), group_name = "group"
#' )
#' ident$rank_ok
#' @export
check_identifiability <- function(
    y,
    Z,
    groups,
    X_hyper,
    family = gaussian(),
    re_coef_names = names(X_hyper),
    group_name = attr(groups, "group_name"),
    verbose = FALSE
) {
  family <- .lmebayes_normalize_family(family)

  if (!is.factor(groups)) {
    stop(
      "'groups' must be a factor; levels(groups) label re_rank/",
      "re_estimable/hyper_rank in the return value.",
      call. = FALSE
    )
  }
  if (!is.character(group_name) || length(group_name) != 1L ||
      is.na(group_name) || !nzchar(group_name)) {
    stop(
      "'group_name' could not be resolved to a single non-empty string; ",
      "pass it explicitly or set attr(groups, 'group_name').",
      call. = FALSE
    )
  }

  n_obs <- nrow(Z)
  if (length(y) != n_obs || length(groups) != n_obs) {
    stop(
      sprintf(
        "'y' (length %d), 'Z' (%d rows), and 'groups' (length %d) must all agree.",
        length(y), n_obs, length(groups)
      ),
      call. = FALSE
    )
  }

  if (!is.list(X_hyper) || is.null(names(X_hyper)) ||
      any(!nzchar(names(X_hyper))) || anyDuplicated(names(X_hyper)) != 0L) {
    stop(
      "'X_hyper' must be a named list with unique, non-empty names ",
      "(one entry per random-effects coefficient).",
      call. = FALSE
    )
  }
  if (!setequal(re_coef_names, names(X_hyper))) {
    stop(
      "'re_coef_names' must match names(X_hyper) exactly (as sets); got ",
      "re_coef_names = {", paste(re_coef_names, collapse = ", "), "}, ",
      "names(X_hyper) = {", paste(names(X_hyper), collapse = ", "), "}.",
      call. = FALSE
    )
  }

  g_levs <- levels(groups)
  for (nm in re_coef_names) {
    Xh <- X_hyper[[nm]]
    if (!is.matrix(Xh) && !inherits(Xh, "Matrix")) {
      stop(
        sprintf("X_hyper[['%s']] must be a matrix.", nm),
        call. = FALSE
      )
    }
    rn <- rownames(Xh)
    if (is.null(rn)) {
      stop(
        sprintf(
          "X_hyper[['%s']] has no rownames; rownames must equal levels(%s).",
          nm, group_name
        ),
        call. = FALSE
      )
    }
    if (!setequal(rn, g_levs)) {
      stop(
        sprintf(
          "X_hyper[['%s']] row names do not match levels(%s).",
          nm, group_name
        ),
        call. = FALSE
      )
    }
    if (nrow(Xh) != length(g_levs)) {
      stop(
        sprintf(
          "X_hyper[['%s']] must have exactly one row per level of %s (%d), got %d.",
          nm, group_name, length(g_levs), nrow(Xh)
        ),
        call. = FALSE
      )
    }
  }

  ## Level 1: per-group rank check -- is Z_j full column rank for each level?
  p_re  <- ncol(Z)
  g_chr <- as.character(groups)
  re_rank <- vapply(
    g_levs,
    function(lev) {
      rows <- which(g_chr == lev)
      Z_j  <- Z[rows, , drop = FALSE]
      nrow(Z_j) >= p_re &&
        Matrix::rankMatrix(Z_j, method = "qr")[1L] == p_re
    },
    logical(1L)
  )

  if (isTRUE(verbose)) {
    message(sprintf(
      "check_identifiability: Level 1 (rank) -- %d of %d %s full column rank.",
      sum(re_rank), length(g_levs), group_name
    ))
  }

  ## Level 1: estimability (classical-glm MLE / residual-df, per family) --
  ## must run before the Level 2 check below, since Level 2 restricts to
  ## *estimable* groups, not merely algebraically full-rank ones.
  glm_est <- .lmebayes_block_glm_estimable(
    y = y, groups = groups, Z = Z, re_rank = re_rank, family = family
  )
  re_estimable <- glm_est$re_estimable
  re_glm_check <- glm_est$re_glm_check

  if (!is.null(re_glm_check) && nrow(re_glm_check) > 0L &&
      !all(re_glm_check$group %in% g_levs)) {
    stop(
      "Internal error: re_glm_check$group contains values not in ",
      "levels(", group_name, ").",
      call. = FALSE
    )
  }

  if (isTRUE(verbose)) {
    message(sprintf(
      "check_identifiability: Level 1 (estimability) -- %d of %d full-rank %s are estimable.",
      sum(re_estimable[re_rank]), sum(re_rank), group_name
    ))
  }

  ## Level 2: hyper-design rank check, restricted to Level-1-estimable groups.
  estimable_levs <- names(re_estimable)[re_estimable]
  hyper_rank <- vapply(
    re_coef_names,
    function(nm) {
      Xh <- X_hyper[[nm]][estimable_levs, , drop = FALSE]
      p  <- ncol(Xh)
      nrow(Xh) >= p && Matrix::rankMatrix(Xh, method = "qr")[1L] == p
    },
    logical(1L)
  )
  hyper_deficient <- !hyper_rank
  rank_ok <- all(hyper_rank)

  if (isTRUE(verbose)) {
    if (isTRUE(rank_ok)) {
      message(sprintf(
        "check_identifiability: Level 2 -- all hyper-design matrices are full rank (%d estimable %s).",
        length(estimable_levs), group_name
      ))
    } else {
      message(sprintf(
        "check_identifiability: Level 2 -- rank-deficient hyper-design matrix(es): %s.",
        paste(names(hyper_rank)[hyper_deficient], collapse = ", ")
      ))
    }
  }

  list(
    re_rank         = re_rank,
    re_estimable    = re_estimable,
    re_glm_check    = re_glm_check,
    hyper_rank      = hyper_rank,
    hyper_deficient = hyper_deficient,
    rank_ok         = rank_ok,
    group_levels    = g_levs,
    group_name      = group_name
  )
}
