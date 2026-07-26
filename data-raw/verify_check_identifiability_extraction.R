## Scratch verification for the check_identifiability() extraction, and for
## its D/group/W argument renaming to match rLMM_reg()/rGLMM_reg() naming
## conventions.
## Confirms:
##   1. check_identifiability()'s roxygen example runs standalone.
##   2. model_setup()'s re_rank/re_estimable/re_glm_check/hyper_rank/
##      hyper_deficient/rank_ok are unchanged (same values as calling
##      check_identifiability() directly on the extracted design fields,
##      keyword-remapped D = design$Z, group = design$groups,
##      W = design$X_hyper).
##   3. Input-validation errors fire as documented.

devtools::load_all(quiet = TRUE)

cat("=== 1. Standalone example ===\n")
set.seed(1)
J <- 8L; n_j <- 6L
group <- factor(rep(paste0("g", seq_len(J)), each = n_j))
x1 <- rnorm(J * n_j)
D <- cbind("(Intercept)" = 1, x1 = x1)
y <- 1 + x1 + rnorm(J * n_j)
W <- list(
  "(Intercept)" = matrix(1, J, 1, dimnames = list(NULL, "(Intercept)")),
  x1            = matrix(1, J, 1, dimnames = list(NULL, "(Intercept)"))
)
ident <- check_identifiability(
  y = y, D = D, group = group, W = W,
  family = gaussian(), group_name = "group", verbose = TRUE
)
stopifnot(isTRUE(ident$rank_ok))
cat("OK: rank_ok =", ident$rank_ok, "\n\n")

cat("=== 2. model_setup() equivalence (gaussian) ===\n")
set.seed(2)
J <- 20L; n_j <- 8L
grp <- factor(rep(paste0("g", seq_len(J)), each = n_j))
x1  <- rnorm(J * n_j)
b0  <- rnorm(J, sd = 2)[as.integer(grp)]
b1  <- rnorm(J, sd = 1)[as.integer(grp)]
y2  <- 1 + x1 + b0 + b1 * x1 + rnorm(J * n_j)
df  <- data.frame(y = y2, x1 = x1, grp = grp)

ms <- model_setup(y ~ x1 + (1 + x1 || grp), data = df, family = gaussian())

ident2 <- check_identifiability(
  y = ms$y, D = ms$Z, group = ms$groups, W = ms$X_hyper,
  family = gaussian(), group_name = ms$group_name
)

stopifnot(identical(ms$re_rank, ident2$re_rank))
stopifnot(identical(ms$re_estimable, ident2$re_estimable))
stopifnot(identical(ms$re_glm_check, ident2$re_glm_check))
stopifnot(identical(ms$hyper_rank, ident2$hyper_rank))
stopifnot(identical(ms$hyper_deficient, ident2$hyper_deficient))
stopifnot(identical(ms$rank_ok, ident2$rank_ok))
cat("OK: model_setup() fields match check_identifiability() directly (gaussian)\n\n")

cat("=== 3. model_setup() equivalence (binomial, with a rank-deficient/non-estimable group) ===\n")
set.seed(3)
J <- 15L; n_j <- 20L
grp3 <- factor(rep(paste0("g", seq_len(J)), each = n_j))
x1b  <- rnorm(J * n_j)
lin  <- -0.5 + x1b + rnorm(J, sd = 1.5)[as.integer(grp3)]
yb   <- rbinom(J * n_j, 1, plogis(lin))
## Force one group to a single outcome level (non-estimable under binomial).
yb[as.integer(grp3) == 1L] <- 0L
df3 <- data.frame(y = yb, x1 = x1b, grp = grp3)

ms3 <- model_setup(y ~ x1 + (1 + x1 || grp), data = df3, family = binomial())
ident3 <- check_identifiability(
  y = ms3$y, D = ms3$Z, group = ms3$groups, W = ms3$X_hyper,
  family = binomial(), group_name = ms3$group_name
)
stopifnot(identical(ms3$re_estimable, ident3$re_estimable))
stopifnot(identical(ms3$re_glm_check, ident3$re_glm_check))
stopifnot(identical(ms3$hyper_rank, ident3$hyper_rank))
stopifnot(identical(ms3$rank_ok, ident3$rank_ok))
stopifnot(!all(ms3$re_estimable))
cat("OK: model_setup() fields match check_identifiability() directly (binomial, non-estimable group present)\n")
cat("re_estimable:\n"); print(ms3$re_estimable)
cat("\n")

cat("=== 4. Input validation errors ===\n")
err <- function(expr) tryCatch({ expr; "NO ERROR" }, error = function(e) conditionMessage(e))

cat("- group not a factor: ",
    err(check_identifiability(y = y, D = D, group = as.character(group),
                               W = W, group_name = "group")),
    "\n")

cat("- unresolvable group_name: ",
    err(check_identifiability(y = y, D = D, group = group, W = W)),
    "\n")

D_bad_names <- D
colnames(D_bad_names) <- c("(Intercept)", "not_x1")
cat("- colnames(D)/names(W) mismatch: ",
    err(check_identifiability(y = y, D = D_bad_names, group = group, W = W,
                               group_name = "group")),
    "\n")

W_bad_nrow <- W
W_bad_nrow[["x1"]] <- W[["x1"]][seq_len(nrow(W[["x1"]]) - 1L), , drop = FALSE]
cat("- W[[k]] wrong nrow: ",
    err(check_identifiability(y = y, D = D, group = group, W = W_bad_nrow,
                               group_name = "group")),
    "\n")

cat("\nAll checks passed.\n")
