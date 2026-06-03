set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)

ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
mu <- ps$mu
V0 <- ps$Sigma
out <- rglmb(
  n = 1000,
  y = ps$y,
  x = as.matrix(ps$x),
  pfamily = dNormal(mu = mu, Sigma = V0),
  family = poisson(),
  weights = rep(1, nrow(ps$x))
)

betastar <- out$coef.mode
x <- out$x
y <- out$y
offset2 <- 0 * y
weights2 <- out$prior.weights

fit <- glmb.wfit(x, y, weights2, offset2, family = poisson(), Bbar = mu, P = solve(V0), betastar)
influence.measures(fit)

print(fit)
print(out$coef.mode)

mu1 <- 0 * mu
V1 <- 0.1 * V0
out2 <- rglmb(
  n = 1000,
  y = ps$y,
  x = as.matrix(ps$x),
  pfamily = dNormal(mu = mu1, Sigma = V1),
  family = poisson(),
  weights = rep(1, nrow(ps$x))
)

Bbar2 <- mu1
betastar2 <- out2$coef.mode
fit2 <- glmb.wfit(x, y, weights2, offset2, family = poisson(), Bbar2, P = solve(V1), betastar2)

influence.measures(fit2)

print(fit2)
print(out2$coef.mode)
