compare_pkg <- function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  p_setup <- Prior_Setup(
    weight ~ group,
    data = data.frame(weight = weight, group = group),
    family = gaussian()
  )
  prior_list <- list(
    mu = p_setup$mu,
    Sigma = p_setup$Sigma,
    shape = p_setup$shape,
    rate = p_setup$rate,
    max_disp_perc = 0.99
  )
  set.seed(1)
  sim <- get("rindepNormalGamma_reg", asNamespace(pkg))(
    n = 2,
    y = p_setup$y,
    x = p_setup$x,
    prior_list = prior_list,
    n_envopt = 1000L,
    use_parallel = FALSE,
    progbar = FALSE
  )
  list(
    pkg = pkg,
    names = names(sim),
    class = class(sim),
    null_names = names(sim)[vapply(sim, is.null, logical(1))]
  )
}

pkg <- commandArgs(trailingOnly = TRUE)[1]
if (length(pkg) == 0L) {
  stop("usage: Rscript compare_rindep_return_structure.R glmbayes|glmbayesCore")
}
if (pkg == "glmbayesCore") {
  pkgload::load_all(quiet = TRUE)
}
out <- compare_pkg(pkg)
saveRDS(out, file = paste0("data-raw/_rindep_return_", pkg, ".rds"))
