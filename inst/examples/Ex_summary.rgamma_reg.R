## summary.rGamma_reg: dGamma prior (dispersion-only; coefficients fixed)
## All five functions (rGamma_reg, rglmb, rlmb, glmb, lmb) use summary.rGamma_reg when prior is dGamma.

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

lm.D9 <- lm(weight ~ group, x = TRUE, y = TRUE)
ps <- Prior_Setup(weight ~ group, family = gaussian())

y <- lm.D9$y
x <- as.matrix(lm.D9$x)
wt <- rep(1, length(y))

## 1. rGamma_reg
out1 <- rGamma_reg(n = 100, y = y, x = x,
  prior_list = list(beta = coef(lm.D9), shape = ps$shape, rate = ps$rate),
  offset = rep(0, length(y)), weights = wt, family = gaussian())
summary(out1)

## 2. rglmb
out2 <- rglmb(n = 100, y = y, x = x,
  pfamily = dGamma(shape = ps$shape, rate = ps$rate, beta = coef(lm.D9)),
  weights = wt, family = gaussian())
summary(out2)

## 3. rlmb
out3 <- rlmb(n = 100, y = y, x = x,
  pfamily = dGamma(shape = ps$shape, rate = ps$rate, beta = coef(lm.D9)),
  weights = wt)
summary(out3)

## 4. glmb
out4 <- glmb(weight ~ group, family = gaussian(),
  pfamily = dGamma(shape = ps$shape, rate = ps$rate, beta = coef(lm.D9)))
summary(out4)

## 5. lmb
out5 <- lmb(weight ~ group,
  pfamily = dGamma(shape = ps$shape, rate = ps$rate, beta = coef(lm.D9)))
summary(out5)
