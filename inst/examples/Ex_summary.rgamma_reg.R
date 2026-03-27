## summary.rGamma_reg: dGamma prior (dispersion-only; coefficients fixed)
## All five functions (rGamma_reg, rglmb, rlmb, glmb, lmb) use summary.rGamma_reg when
## prior is dGamma.
##
## This example uses the Boston data in a two-step setup:
## 1) fit lm() to obtain fixed coefficient values (beta),
## 2) estimate dispersion under dGamma prior using default Prior_Setup().
data("Boston", package = "MASS")

predictors <- setdiff(names(Boston), "medv")
Boston_centered <- Boston
Boston_centered[predictors] <- scale(Boston[predictors], center = TRUE, scale = FALSE)

form <- medv ~
  crim + zn +
  indus + chas + nox + age + dis + rad + tax + ptratio + black + lstat + rm

lm.boston <- lm(form, data = Boston_centered, x = TRUE, y = TRUE)
ps.boston <- Prior_Setup(form, gaussian(), data = Boston_centered)

y <- lm.boston$y
x <- as.matrix(lm.boston$x)
wt <- rep(1, length(y))

## 1. rGamma_reg
out1 <- rGamma_reg(
  n = 1000,
  y = y,
  x = x,
  prior_list = list(beta = coef(lm.boston), shape = ps.boston$shape, rate = ps.boston$rate),
  offset = rep(0, length(y)),
  weights = wt,
  family = gaussian()
)
summary(out1)

## 2. rglmb
out2 <- rglmb(n = 1000, y = y, x = x,
  pfamily = dGamma(shape = ps.boston$shape, rate = ps.boston$rate, beta = coef(lm.boston)),
  weights = wt, family = gaussian())
summary(out2)

## 3. rlmb
out3 <- rlmb(n = 1000, y = y, x = x,
  pfamily = dGamma(shape = ps.boston$shape, rate = ps.boston$rate, beta = coef(lm.boston)),
  weights = wt)
summary(out3)

## 4. glmb
out4 <- glmb(form, data = Boston_centered, family = gaussian(),
  pfamily = dGamma(shape = ps.boston$shape, rate = ps.boston$rate, beta = coef(lm.boston)))
summary(out4)

## 5. lmb
out5 <- lmb(form, data = Boston_centered,
  pfamily = dGamma(shape = ps.boston$shape, rate = ps.boston$rate, beta = coef(lm.boston)))
summary(out5)
