## Prior_SetupBlock — independent Prior_Setup() per Species block (iris)

data("iris", package = "datasets")

ps_block <- Prior_SetupBlock(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  block = "Species",
  data = iris,
  family = gaussian()
)

names(ps_block)
ps_block$setosa$mu
length(ps_block$setosa$mu)
