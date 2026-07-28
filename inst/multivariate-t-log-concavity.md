# Local Log-Concavity of the Multivariate Student-*t* Distribution

## Motivation

The multivariate Student-*t* density is **not** globally log-concave — that's
precisely what produces its heavy tails. But near the mean it behaves like a
(scaled) Gaussian and *is* locally log-concave. This note derives the exact
boundary of that region.

To be precise about what's being characterized: the density $f(x)$ is
positive everywhere. What can fail is positive-definiteness of the **negative
Hessian of the log-density** (the local curvature / information matrix). That
is the object characterized below.

## Setup

For $x \in \mathbb{R}^p$, mean $\mu$, scale matrix $\Sigma$, and $\nu$ degrees
of freedom:

$$
f(x) \;\propto\; \left[1 + \tfrac{1}{\nu}(x-\mu)^\top \Sigma^{-1} (x-\mu)\right]^{-(\nu+p)/2}
$$

Define:

- $A = \Sigma^{-1}$
- $q = x - \mu$
- $d = q^\top A q$ (squared Mahalanobis distance)
- $s = 1 + d/\nu$

## Gradient and Hessian of the log-density

$$
\nabla \log f = -\frac{\nu+p}{\nu s}\, A q
$$

$$
H = \nabla^2 \log f = -\frac{\nu+p}{\nu s} A \;+\; \frac{2(\nu+p)}{\nu^2 s^2} (Aq)(Aq)^\top
$$

So the negative Hessian has the rank-one structure

$$
-H = c\left[A - k\, v v^\top\right], \qquad
c = \frac{\nu+p}{\nu s} > 0, \quad
v = Aq, \quad
k = \frac{2}{\nu s}
$$

## When is $-H$ positive definite?

This is a rank-one downdate of the positive definite matrix $A$. Diagonalizing

$$
A^{-1/2}(A - k v v^\top) A^{-1/2} = I - k\, u u^\top, \qquad u = A^{-1/2} v
$$

gives eigenvalues:

- $1$, with multiplicity $p-1$
- $1 - k\|u\|^2$, with multiplicity $1$

and $\|u\|^2 = v^\top A^{-1} v = q^\top A \Sigma A q = q^\top A q = d$.

Therefore:

$$
-H \succ 0 \iff k\,d < 1 \iff \frac{2d}{\nu s} < 1 \iff \nu(1 + d/\nu) > 2d \iff d < \nu
$$

## Result

$$
\boxed{-\nabla^2 \log f(x) \succ 0 \quad \Longleftrightarrow \quad (x-\mu)^\top \Sigma^{-1} (x-\mu) < \nu}
$$

**The log-density is strictly concave exactly inside the Mahalanobis
ellipsoid of squared radius $\nu$.** The degrees of freedom directly set the
"radius of concavity" in Mahalanobis units.

## Additional notes

- **Direction of first failure.** The vanishing eigenvalue's eigenvector
  points along $u = A^{-1/2} v \propto \Sigma^{-1/2} q$ — the *radial*
  direction from $\mu$ through $x$. The $p-1$ tangential directions retain
  eigenvalue $c$ (positive) regardless of $d$. So concavity is lost radially
  first, exactly at $d = \nu$; just past that boundary the density is
  saddle-shaped before convex behavior eventually dominates further out.

- **Gaussian limit.** As $\nu \to \infty$, $s \to 1$ and $-H \to \Sigma^{-1}$
  for all $x$, recovering global log-concavity — consistent with $d < \nu$
  holding for any finite $d$ as $\nu \to \infty$.

- **Small $\nu$.** For $\nu \lesssim 2$, the concave region shrinks
  considerably in Mahalanobis terms, consistent with the visibly "peakier,
  then heavy-tailed" shape of low-df Student-*t* densities.

- **Practical relevance.** For Laplace approximations or Newton-type
  optimization on *t*-distributed (or *t*-prior) posteriors, the "trust
  region" in which treating the mode's local curvature as informative is
  valid is bounded by $d < \nu$ — not global, as it would be under a Gaussian
  assumption. This matters for multivariate-*t* proposals or priors in MCMC
  where a local quadratic approximation is being relied on.

## Application: the ellipsoid for OLS regression coefficients

The Mahalanobis ellipsoid above specializes cleanly when $\mu$ and $\Sigma$
arise from ordinary least squares.

### Setup

With design matrix $X \in \mathbb{R}^{n\times p}$ (rank $p$) and response
$y \in \mathbb{R}^n$, assume the classical model $y = X\beta + \varepsilon$,
$\varepsilon \sim N(0,\sigma^2 I)$. Then:

$$
\hat\beta = (X^\top X)^{-1}X^\top y, \qquad
\hat\sigma^2 = \frac{\text{RSS}}{n-p}, \qquad
\text{RSS} = \|y - X\hat\beta\|^2
$$

Studentizing $\hat\beta$ by the *estimated* $\hat\sigma^2$ (instead of the
unknown $\sigma^2$) gives the classical result that $\hat\beta$ follows a
**multivariate Student-$t$** distribution with:

$$
\mu = \hat\beta, \qquad \Sigma = \hat\sigma^2 (X^\top X)^{-1}, \qquad \nu = n-p
$$

(Same construction as the univariate case: $Z/\sqrt{W/\nu}$ with
$Z \sim N(0,\Sigma)$ and $W = (n-p)\hat\sigma^2/\sigma^2 \sim \chi^2_{n-p}$,
independent of $\hat\beta$.)

### Substituting into the ellipsoid

Recall the boundary derived above: $(x-\mu)^\top \Sigma^{-1}(x-\mu) < \nu$.
Here $\Sigma^{-1} = \frac{1}{\hat\sigma^2}X^\top X$, so:

$$
d(x) = \frac{1}{\hat\sigma^2}(x-\hat\beta)^\top X^\top X (x-\hat\beta) < \nu = n-p
$$

Multiplying both sides by $\hat\sigma^2 = \text{RSS}/(n-p)$, the $(n-p)$
factor cancels cleanly:

$$
\boxed{(x-\hat\beta)^\top X^\top X\,(x-\hat\beta) \;<\; \text{RSS} \;=\; \|y - X\hat\beta\|^2}
$$

with $\hat\beta = (X^\top X)^{-1}X^\top y$.

### Writing it purely in terms of $X$ and $y$

Let $H = X(X^\top X)^{-1}X^\top$ be the hat/projection matrix. Then
$\text{RSS} = y^\top(I-H)y$, and the ellipsoid becomes:

$$
\left(x - (X^\top X)^{-1}X^\top y\right)^\top X^\top X\,\left(x - (X^\top X)^{-1}X^\top y\right) \;<\; y^\top(I-H)y
$$

### Interpretation

- The region is an ellipsoid in $\beta$-space centered at $\hat\beta$, with
  shape governed by $X^\top X$ (the design's information geometry) and
  radius set by the **residual sum of squares** itself.
- This is essentially the same ellipsoid that appears in the classical
  **Scheffé / joint confidence region** for $\beta$:
  $(\beta-\hat\beta)^\top X^\top X(\beta-\hat\beta) \le p\,\hat\sigma^2 F_{p,n-p,\alpha}$.
  The log-concavity boundary derived here is the special case where the
  F-quantile is replaced by $\frac{n-p}{p}$ — the shape is identical, just
  RSS-scaled rather than $\hat\sigma^2$-scaled. (Worth checking the exact
  correspondence carefully if a confidence-level interpretation is wanted.)
- Practically: past this ellipsoid, the multivariate-$t$ likelihood surface
  for $\beta$ (built from studentized residuals) stops being locally
  log-concave — so Newton-type steps or Laplace approximations around
  $\hat\beta$ are only trustworthy strictly inside the region where the
  quadratic form in $X^\top X$ stays below the observed RSS.
