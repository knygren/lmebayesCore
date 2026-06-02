## Multi-response lmb (iris): formula interface, mlmb class

set.seed(42)
data("iris", package = "datasets")

form_multi <- cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species

ps_multi <- multi_prior_setup(
  form_multi,
  data = iris,
  family = gaussian()
)

pfamily_list <- lapply(ps_multi, function(ps) {
  dNormal_Gamma(
    mu = ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate
  )
})

out_lmb <- lmb(
  form_multi,
  pfamily_list = pfamily_list,
  data = iris,
  n = 150L,
  use_parallel = FALSE
)

print(out_lmb)
summary(out_lmb)
