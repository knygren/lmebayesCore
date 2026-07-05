devtools::load_all(quiet = TRUE)
set.seed(1L)
J <- 4L
g <- factor(rep(seq_len(J), each = 3L))
n <- length(g)
x <- matrix(1, n, 1L)
xh <- stats::setNames(list(matrix(1, J, 1L)), "(Intercept)")
P <- matrix(1)
pf <- dNormal(c(`(Intercept)` = 0), matrix(1), dispersion = 1)
pfl <- stats::setNames(list(pf), "(Intercept)")
r_fix <- two_block_rate_from_pfamily_list(
  x, g, xh, list(P = P, dispersion = 1, ddef = FALSE), pfl
)
r_up <- two_block_rate_from_pfamily_list(
  x, g, xh, list(P = P, dispersion = 50, ddef = FALSE), pfl
)
stopifnot(r_up$lambda_star >= r_fix$lambda_star)
stopifnot(
  glmbayesCore:::.two_block_pilot_will_run(
    TRUE, NULL, 0.02, 0.05, FALSE, TRUE
  )
)
stopifnot(
  !glmbayesCore:::.two_block_pilot_will_run(
    TRUE, NULL, 0.02, 0.05, FALSE, FALSE
  )
)
message("disp_upper rate monotonicity OK: ",
        "lambda(fix)=", r_fix$lambda_star,
        " lambda(upper)=", r_up$lambda_star)
