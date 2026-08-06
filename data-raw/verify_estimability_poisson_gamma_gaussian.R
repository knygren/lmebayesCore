## Verify: the estimability check (re_estimable/re_glm_check) generalized
## from binomial-only to also cover poisson()/Gamma() (glm-fit based) and
## gaussian() (new residual-degrees-of-freedom-for-dispersion check).
##
## Mirrors the earlier data-raw/verify_re_estimable_move.R pattern: exercise
## this only through the real Prior_Setup_GLMM()/model_setup() call
## chain -- no hand-rolled glm() calls or fabricated coefficient/vcov
## values.
devtools::load_all(".", quiet = TRUE)

form <- y ~ x1 + (1 + x1 || group)

## ---------------------------------------------------------------------
## 1. poisson(): one group with an all-zero count response.
##    (A within-group all-zero count is *not* a global domain violation --
##    the pooled glmer(poisson()) reference fit converges fine -- so this
##    goes through the normal Prior_Setup_GLMM() call chain.)
## ---------------------------------------------------------------------
set.seed(1)
J   <- 20L
n_j <- 20L
dat_pois <- data.frame(
  group = factor(rep(paste0("g", seq_len(J)), each = n_j)),
  x1    = stats::rnorm(J * n_j, sd = 0.3)
)
b0 <- stats::setNames(stats::rnorm(J, sd = 0.3), levels(dat_pois$group))
b1 <- stats::setNames(stats::rnorm(J, sd = 0.15), levels(dat_pois$group))
eta <- 1.5 + b0[as.character(dat_pois$group)] +
  (0.3 + b1[as.character(dat_pois$group)]) * dat_pois$x1
dat_pois$y <- stats::rpois(nrow(dat_pois), lambda = exp(eta))
dat_pois$y[dat_pois$group == "g1"] <- 0L  ## force one group all-zero

ps_pois <- Prior_Setup_GLMM(form, data = dat_pois, family = poisson(), pwt = 0.05)
d_pois  <- ps_pois$design

cat("=== poisson(): all-zero-response group ===\n")
print(d_pois$re_glm_check[d_pois$re_glm_check$group == "g1", ])
note_g1 <- d_pois$re_glm_check$note[d_pois$re_glm_check$group == "g1"]
stopifnot(
  isTRUE(d_pois$re_rank[["g1"]]),                  ## design itself is full rank
  identical(unname(d_pois$re_estimable[["g1"]]), FALSE), ## but not estimable
  grepl("all-zero", note_g1, fixed = TRUE)
)
other_pois <- setdiff(levels(dat_pois$group), "g1")
stopifnot(all(d_pois$re_estimable[other_pois] == d_pois$re_rank[other_pois]))
cat("OK\n\n")

## ---------------------------------------------------------------------
## 2. Gamma(): residual-df check (n_j == p_j, otherwise-valid response).
## ---------------------------------------------------------------------
set.seed(2)
J2     <- 12L
n_main <- 10L
grp_ids2 <- paste0("h", seq_len(J2))
dat_gamma <- do.call(rbind, lapply(grp_ids2, function(g) {
  n_g <- if (identical(g, "h1")) 2L else n_main  ## h1: exactly p_j = 2 obs
  data.frame(group = g, x1 = stats::rnorm(n_g, sd = 0.5))
}))
dat_gamma$group <- factor(dat_gamma$group, levels = grp_ids2)
b0g <- stats::setNames(stats::rnorm(J2, sd = 0.2), grp_ids2)
b1g <- stats::setNames(stats::rnorm(J2, sd = 0.1), grp_ids2)
eta_g <- 1 + b0g[as.character(dat_gamma$group)] +
  (0.2 + b1g[as.character(dat_gamma$group)]) * dat_gamma$x1
mu_g <- exp(eta_g)
dat_gamma$y <- stats::rgamma(nrow(dat_gamma), shape = 5, rate = 5 / mu_g)

ps_gamma <- Prior_Setup_GLMM(
  form, data = dat_gamma, family = Gamma(link = "log"), pwt = 0.05
)
d_gamma <- ps_gamma$design

cat("=== Gamma(): n_j == p_j group (zero residual df for dispersion) ===\n")
print(d_gamma$re_glm_check[d_gamma$re_glm_check$group == "h1", ])
note_h1 <- d_gamma$re_glm_check$note[d_gamma$re_glm_check$group == "h1"]
stopifnot(
  isTRUE(d_gamma$re_rank[["h1"]]),
  identical(unname(d_gamma$re_estimable[["h1"]]), FALSE),
  grepl("residual degrees of freedom", note_h1, fixed = TRUE)
)
other_gamma <- setdiff(grp_ids2, "h1")
stopifnot(all(d_gamma$re_estimable[other_gamma]))
cat("OK\n\n")

## ---------------------------------------------------------------------
## 3. Gamma(): domain check (non-positive response in one group).
##    A non-positive value anywhere breaks the *pooled* glmer(Gamma())
##    reference fit itself, so this is exercised via model_setup() directly
##    with fit_mer = FALSE (an existing, documented model_setup() mode used
##    by glmerb() in lmebayes) rather than Prior_Setup_GLMM().
## ---------------------------------------------------------------------
dat_gamma_bad <- dat_gamma
bad_rows <- which(dat_gamma_bad$group == "h2")
dat_gamma_bad$y[bad_rows[1L]] <- -0.5  ## inject one negative response

design_bad <- model_setup(
  form, data = dat_gamma_bad, family = Gamma(link = "log"), fit_mer = FALSE
)

cat("=== Gamma(): non-positive-response group ===\n")
print(design_bad$re_glm_check[design_bad$re_glm_check$group == "h2", ])
note_h2 <- design_bad$re_glm_check$note[design_bad$re_glm_check$group == "h2"]
stopifnot(
  isTRUE(design_bad$re_rank[["h2"]]),
  identical(unname(design_bad$re_estimable[["h2"]]), FALSE),
  grepl("non-positive", note_h2, fixed = TRUE)
)
cat("OK\n\n")

## ---------------------------------------------------------------------
## 4. gaussian(): residual-df check (n_j == p_j) -- behavior change: this
##    group previously would have been re_estimable = re_rank = TRUE, and
##    design$re_glm_check was always NULL for gaussian().
## ---------------------------------------------------------------------
set.seed(3)
J3 <- 15L
grp_ids3 <- paste0("k", seq_len(J3))
dat_g <- do.call(rbind, lapply(grp_ids3, function(g) {
  n_g <- if (identical(g, "k1")) 2L else 15L  ## k1: exactly p_j = 2 obs
  data.frame(group = g, x1 = stats::rnorm(n_g))
}))
dat_g$group <- factor(dat_g$group, levels = grp_ids3)
b0k <- stats::setNames(stats::rnorm(J3, sd = 1.5), grp_ids3)
b1k <- stats::setNames(stats::rnorm(J3, sd = 0.8), grp_ids3)
dat_g$y <- 2 + b0k[as.character(dat_g$group)] +
  (1.5 + b1k[as.character(dat_g$group)]) * dat_g$x1 +
  stats::rnorm(nrow(dat_g), sd = 1)

ps_g <- Prior_Setup_GLMM(form, data = dat_g, pwt = 0.05)  ## family = gaussian() default
d_g  <- ps_g$design

cat("=== gaussian(): n_j == p_j group (zero residual df for dispersion) ===\n")
print(d_g$re_glm_check[d_g$re_glm_check$group == "k1", ])
note_k1 <- d_g$re_glm_check$note[d_g$re_glm_check$group == "k1"]
stopifnot(
  !is.null(d_g$re_glm_check),  ## re_glm_check no longer always NULL for gaussian
  isTRUE(d_g$re_rank[["k1"]]),
  identical(unname(d_g$re_estimable[["k1"]]), FALSE),
  grepl("residual degrees of freedom", note_k1, fixed = TRUE)
)
other_g <- setdiff(grp_ids3, "k1")
stopifnot(all(d_g$re_estimable[other_g]))
cat("OK\n\n")

cat("ALL OK\n")
