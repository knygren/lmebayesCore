pkgload::load_all(quiet = TRUE)
suppressPackageStartupMessages(library(glmbayes))

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
sim_old <- glmbayes::rindepNormalGamma_reg(
  n = 2, y = p_setup$y, x = p_setup$x, prior_list = prior_list,
  n_envopt = 1000L, use_parallel = FALSE, progbar = FALSE
)
set.seed(1)
sim_new <- rindepNormalGamma_reg(
  n = 2, y = p_setup$y, x = p_setup$x, prior_list = prior_list,
  n_envopt = 1000L, use_parallel = FALSE, progbar = FALSE
)

cat("names match:", identical(names(sim_old), names(sim_new)), "\n")
cat("names:\n")
print(names(sim_old))
cat("class old:", paste(class(sim_old), collapse = ", "), "\n")
cat("class new:", paste(class(sim_new), collapse = ", "), "\n")
cat("\nNULL elements (old):\n")
print(names(sim_old)[vapply(sim_old, is.null, logical(1))])
cat("NULL elements (new):\n")
print(names(sim_new)[vapply(sim_new, is.null, logical(1))])
cat("\nstr glmbayes (max.level=1):\n")
str(sim_old, max.level = 1)
cat("\nstr glmbayesCore (max.level=1):\n")
str(sim_new, max.level = 1)
cat("\nExtra names in new only:\n")
print(setdiff(names(sim_new), names(sim_old)))
cat("Extra names in old only:\n")
print(setdiff(names(sim_old), names(sim_new)))
