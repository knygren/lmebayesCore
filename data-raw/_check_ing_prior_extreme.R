suppressPackageStartupMessages({
  devtools::load_all("c:/Rpackages/glmbayesCore", quiet = TRUE)
  devtools::load_all("c:/Rpackages/lmebayes", quiet = TRUE)
  library(bayesrules)
})
source("c:/Rpackages/glmbayesCore/data-raw/_check_ing_prior.R", local = TRUE)

design <- ps$design
k <- "(Intercept)"
X_k <- design$X_hyper[[k]]
for (lab in c("tiny", "huge", "mixed")) {
  b <- switch(
    lab,
    tiny  = rep(1e-6, nrow(X_k)),
    huge  = rep(1e6, nrow(X_k)),
    mixed = seq(-500, 500, length.out = nrow(X_k))
  )
  cat("\n---", lab, "b ---\n")
  tryCatch(
    {
      fit <- glmbayesCore::rglmb(
        1L, b, X_k, gaussian(), pf[[k]],
        verbose = FALSE, use_parallel = FALSE
      )
      cat("OK tau2 =", fit$dispersion, "\n")
    },
    error = function(e) cat("ERROR:", conditionMessage(e), "\n")
  )
}

## Wrong prior: conjugate shape (missing + p/2) like dNormal_Gamma
pl_wrong <- pf[["(Intercept)"]]$prior_list
shape_wrong <- (ps$n_prior_dispersion[["(Intercept)"]] + 1) / 2
rate_wrong <- pl_wrong$rate
pf_bad <- glmbayesCore::dIndependent_Normal_Gamma(
  mu = pl_wrong$mu,
  Sigma = pl_wrong$Sigma,
  shape = shape_wrong,
  rate = rate_wrong,
  disp_lower = pl_wrong$disp_lower,
  disp_upper = pl_wrong$disp_upper
)
pf_bad$prior_list$tau2_ref <- ps$prior_list[["(Intercept)"]]$dispersion_fixef
cat("\n=== wrong shape (NG not ING): shape=", shape_wrong, " vs correct=", pl_wrong$shape, "===\n")
b <- rnorm(nrow(X_k), 100, 10)
tryCatch(
  {
    fit <- glmbayesCore::rglmb(
      1L, b, X_k, gaussian(), pf_bad,
      verbose = FALSE, use_parallel = FALSE
    )
    cat("OK tau2 =", fit$dispersion, "\n")
  },
  error = function(e) cat("ERROR:", conditionMessage(e), "\n")
)
