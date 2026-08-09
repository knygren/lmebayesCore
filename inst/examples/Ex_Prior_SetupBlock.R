## Prior_SetupBlock — independent Prior_Setup() per Species block (iris)

data("iris", package = "datasets")

ps_block <- Prior_SetupBlock(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  block = "Species",
  data = iris,
  family = gaussian()
)

print(ps_block)
names(ps_block)
ps_block$setosa$mu

## Default: dNormal() per block (also valid for binomial/poisson blocks)
pf <- pfamily_list(ps_block)
names(pf)
attr(pf, "ptypes")

## Gaussian-only alternatives
pf_ng <- pfamily_list(ps_block, ptypes = "dNormal_Gamma")
pf_ing <- pfamily_list(ps_block, ptypes = "dIndependent_Normal_Gamma")
print(pf_ng, components = "setosa")
