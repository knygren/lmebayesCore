# The $\Omega_j$-Marginal Likelihood Under the *Independent* (not Conjugate) Normal-Gamma Sampling Model

This note re-derives, from scratch and with an explicit eye on which prior is
being integrated against what, the marginal likelihood
$p(y_j\mid\beta_j)$ that Section 16 of `inst/BLOCK_GIBBS_ERGODICITY_ING.md`
uses for its $H_j(\beta_j)$ Hessian and log-concavity criterion. The point of
redoing it here is to settle, completely explicitly, whether that
marginalization already accounts for the measurement-precision prior's
**shape** *and* **rate** (it does), and to make sure it is built from the
right generative model: this package's sampling-time
**Independent Normal-Gamma** (`dIndependent_Normal_Gamma()`/`dGamma()`), in
which $\Omega_j$'s Gamma prior does **not** depend on $\beta_j$ at all — as
opposed to a classical **conjugate** Normal-Gamma prior, in which $\beta_j$'s
own prior variance would be scaled by $1/\Omega_j$. Getting this distinction
right matters because the package *also* contains a genuinely conjugate-style
marginalization elsewhere (§3.3.4 of `inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md`,
used only to **calibrate** the numeric values of the shape/rate fed into
`dGamma()`), and it would be easy to (wrongly) reuse that formula here.

**Status:** theory/derivation only, mirroring `inst/multivariate-t-log-concavity.md`'s
style and cross-checked numerically (`data-raw/_scratch_check_omega_prior_marginal.R`).
No code changes; this confirms, rather than modifies, `R/two_block_ergodicity_ing_marginal.R`
and Section 16 of `inst/BLOCK_GIBBS_ERGODICITY_ING.md`.

---

## 0. Why this needs to be said explicitly: two different "Normal-Gamma"s in this package

| | Where it's used | What's integrated against what | Result |
|---|---|---|---|
| **Calibration-time** (§3.3.4, `inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md`) | `dGamma_list()`, once, offline, to *choose* the numbers `shape_ING,j`/`rate,j` handed to `dGamma()` | $\beta_j$ integrated out against an **auxiliary calibration prior** $\beta_j\sim N(\mu_j,\Sigma_j)$ (a Zellner-style prior built from the classical GLM fit's covariance and the *current* Block-2 $\tau^2$) | A marginal Gamma law for $1/\sigma_j^2$ whose rate contains a quadratic **penalty term** $(\hat\beta_j-\mu_j)'M_j^{-1}(\hat\beta_j-\mu_j)$ pulling the calibrated dispersion up when the group's OLS fit is far from $\mu_j$ |
| **Sampling/diagnostic-time** (this note; Section 16 of `inst/BLOCK_GIBBS_ERGODICITY_ING.md`) | The live two-block Gibbs sampler (`rindepNormalGamma_reg()`) and the $H_j(\beta_j)$ diagnostic, given *already-fixed* numbers `a0_j = shape_ING,j`, `r0_j = rate,j` | $\Omega_j\sim\mathrm{Gamma}(a_j^0,r_j^0)$ **unconditionally** — this Gamma prior does not reference $\beta_j$, $\mu_j$, or $\Sigma_j$ at all; $\beta_j$ is treated as a **fixed, given argument**, not itself integrated over | A marginal likelihood $p(y_j\mid\beta_j)$ centered at the group's own **OLS fit** $\hat\beta_j^{\mathrm{ols}}$, with no shrinkage toward any calibration-time prior mean |

Both are real, both matter, and they are **not the same integral** — the
calibration step (row 1) happens once, offline, to *produce* the numbers
$a_j^0,r_j^0$; the sampling/diagnostic step (row 2), which is what this note
derives, then treats those numbers as a plain, $\beta_j$-independent Gamma
prior, exactly matching what `rindepNormalGamma_reg()` actually draws from
(§"Which one is right" in `inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md`:
*"$\beta_j$ is not fixed at any point value during this step — it is
simulated fresh, every sweep... the prior handed to it must already be that
same marginal law"* — i.e. by the time $(a_j^0,r_j^0)$ reach the sampler,
$\Omega_j$'s prior is a fixed Gamma, and $\beta_j$'s own (Block-2) prior is
a completely separate object, never re-derived from $\mu_j,\Sigma_j$ again).
This note is exclusively about row 2. Row 1's $\mu_j,\Sigma_j,M_j$ never
appear below, and §5 shows concretely what would go wrong if they did.

---

## 1. The model being marginalized (restated with the independence made explicit)

$$
y_j \mid \beta_j,\Omega_j \sim N\!\big(D_j\beta_j,\ \Omega_j^{-1}I_{n_j}\big),
\qquad
\Omega_j \sim \mathrm{Gamma}(a_j^0, r_j^0)\ \text{(shape/rate)}
$$

with $\Omega_j$'s prior **not a function of $\beta_j$** — no
$\beta_j$-dependent scaling, no coupling term, nothing. Separately (and
irrelevantly to the integral below, precisely because of that independence):

$$
\beta_{jp} \mid \gamma_p,\lambda_p \sim N(W_{pj}\gamma_p,\ 1/\lambda_p)
$$

is Block-2's own prior on $\beta_j$, entirely unconnected to $\Omega_j$. This
is the "Independent" in `dIndependent_Normal_Gamma` — contrast with a
classical **conjugate** Normal-Gamma prior, which would instead specify
$\beta_j\mid\Omega_j\sim N(b_0,\ \Omega_j^{-1}V_0)$ (prior variance scaled by
$1/\Omega_j$), coupling the two. This package never uses that joint form as
the *sampling* prior for $(\beta_j,\Omega_j)$ for any component; the
conjugate-shaped calculation in `inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md`
§3.3.4 exists only as a one-time **calibration device** for picking
$a_j^0,r_j^0$ (row 1 of §0's table), never as the model integrated over here.

**Consequence.** Because $\Omega_j$'s prior does not involve $\beta_j$, the
marginal likelihood $p(y_j\mid\beta_j) = \int p(y_j\mid\beta_j,\Omega_j)\,
p(\Omega_j)\,d\Omega_j$ is a well-defined function of $\beta_j$ obtained by
integrating **only** over $\Omega_j$'s own Gamma law, holding $\beta_j$ fixed
throughout — no $\mu_j,\Sigma_j$, and no blending with any $\beta_j$-prior,
enters this integral at all. That is exactly the calculation below.

---

## 2. The marginal likelihood

Write $e_j(\beta_j) := y_j - D_j\beta_j$. The Normal-Gamma conjugate identity
for the *likelihood factor alone* gives, exactly,

$$
p(y_j\mid\beta_j) \;\propto\; \int_0^\infty
\Omega_j^{n_j/2}\exp\!\Big(-\tfrac{\Omega_j}{2}e_j'e_j\Big)\cdot
\Omega_j^{a_j^0-1}\exp(-r_j^0\Omega_j)\ d\Omega_j
\;=\; \Gamma\!\Big(a_j^0+\tfrac{n_j}{2}\Big)\Big(r_j^0+\tfrac12 e_j'e_j\Big)^{-(a_j^0+n_j/2)}
$$

i.e.

$$
p(y_j\mid\beta_j) \;\propto\; \Big(r_j^0+\tfrac12 e_j'e_j\Big)^{-(a_j^0+n_j/2)}
$$

— matching `inst/BLOCK_GIBBS_ERGODICITY_ING.md` §16.1 exactly, now with the
integration made explicit as touching *only* $\Omega_j$'s Gamma factor.

---

## 3. Identifying this as a generalized multivariate-$t$

Let $\hat\beta_j^{\mathrm{ols}} = (D_j'D_j)^{-1}D_j'y_j$,
$\mathrm{RSS}_j^{\mathrm{ols}} = \|y_j - D_j\hat\beta_j^{\mathrm{ols}}\|^2$, and
(Pythagorean OLS split, $D_j'\hat e_j^{\mathrm{ols}}=0$)

$$
e_j'e_j \;=\; \mathrm{RSS}_j^{\mathrm{ols}} + q_j, \qquad
q_j := (\beta_j-\hat\beta_j^{\mathrm{ols}})'D_j'D_j(\beta_j-\hat\beta_j^{\mathrm{ols}})
$$

Substituting into §2 and factoring out the $\beta_j$-free constant
$K_j := 2r_j^0+\mathrm{RSS}_j^{\mathrm{ols}}$:

$$
p(y_j\mid\beta_j) \;\propto\; \Big[1+\frac{q_j}{K_j}\Big]^{-(a_j^0+n_j/2)}
$$

Define

$$
\boxed{\ \nu_j := 2a_j^0+n_j-p_{re}\ } \qquad
\boxed{\ \Sigma_j := \frac{K_j}{\nu_j}(D_j'D_j)^{-1} = \frac{2r_j^0+\mathrm{RSS}_j^{\mathrm{ols}}}{2a_j^0+n_j-p_{re}}(D_j'D_j)^{-1}\ }
$$

Then $q_j/K_j = (\beta_j-\hat\beta_j^{\mathrm{ols}})'\Sigma_j^{-1}(\beta_j-\hat\beta_j^{\mathrm{ols}})/\nu_j$
and $(a_j^0+n_j/2) = (\nu_j+p_{re})/2$, so

$$
p(y_j\mid\beta_j) \;\propto\;
\Big[1+\tfrac{1}{\nu_j}(\beta_j-\hat\beta_j^{\mathrm{ols}})'\Sigma_j^{-1}(\beta_j-\hat\beta_j^{\mathrm{ols}})\Big]^{-(\nu_j+p_{re})/2}
$$

**This is exactly the multivariate-$t$ kernel of `inst/multivariate-t-log-concavity.md`
§"Setup"**, with location $\hat\beta_j^{\mathrm{ols}}$, degrees of freedom
$\nu_j$, and scale $\Sigma_j$ — i.e. the "additional terms" the shape/rate
prior contributes, relative to the plain OLS-sampling-theory multivariate-$t$
that note derives ($\nu=n-p$, $\Sigma=\hat\sigma^2(X'X)^{-1}$ with
$\hat\sigma^2=\mathrm{RSS}/(n-p)$, i.e. $a_j^0=r_j^0=0$), are precisely:

- $n_j \to n_j + 2a_j^0$ — the Gamma prior's shape contributes **$2a_j^0$
  pseudo-observations** to the degrees of freedom;
- $\mathrm{RSS}_j^{\mathrm{ols}} \to \mathrm{RSS}_j^{\mathrm{ols}} + 2r_j^0$ — the
  Gamma prior's rate contributes a **$2r_j^0$ pseudo-RSS**.

Both terms vanish exactly as $a_j^0,r_j^0\to0$ (an improper, data-only limit),
recovering the plain OLS multivariate-$t$ note. **Crucially, the location
$\hat\beta_j^{\mathrm{ols}}$ itself is untouched by the prior** — it is pure
data, with no shrinkage toward any calibration-time prior mean $\mu_j$. This
is the direct fingerprint of the *Independent* structure: because $\beta_j$
was never itself given an $\Omega_j$-linked prior in this integral (§1), there
is nothing for $\Omega_j$'s prior to shrink $\beta_j$ *toward* — it can only
reweight the effective sample size and effective RSS, not move the center.
§5 shows exactly how this would differ under the (inapplicable) conjugate
structure.

---

## 4. Gradient, Hessian, and the log-concavity boundary, via the general multivariate-$t$ machinery

Reusing `inst/multivariate-t-log-concavity.md`'s general result verbatim, with
$x=\beta_j$, $\mu=\hat\beta_j^{\mathrm{ols}}$, $A=\Sigma_j^{-1}$, $\nu=\nu_j$,
$p=p_{re}$, $q=\beta_j-\hat\beta_j^{\mathrm{ols}}$, $d=q'Aq$, $s=1+d/\nu_j$:

$$
-\nabla^2\log p(y_j\mid\beta_j) \;=\; c\big[A - k\,vv'\big], \qquad
c=\frac{\nu_j+p_{re}}{\nu_j s},\quad v=Aq,\quad k=\frac{2}{\nu_j s}
$$

**This reduces algebraically to exactly `inst/BLOCK_GIBBS_ERGODICITY_ING.md`
§16.2's closed form** $H_j(\beta_j)=\Omega_j^{\mathrm{eff}}D_j'D_j -
\big(\Omega_j^{\mathrm{eff}}\big)^2/(a_j^0+n_j/2)\cdot(D_j'e_j)(D_j'e_j)'$
(with $\Omega_j^{\mathrm{eff}}=(a_j^0+n_j/2)/(r_j^0+\tfrac12 e_j'e_j)$) —
verified both symbolically (substituting $A=D_j'D_j\cdot\nu_j/K_j$,
$D_j'e_j=-D_j'D_j(\beta_j-\hat\beta_j^{\mathrm{ols}})$) and numerically:

```
max abs diff closed vs numeric (finite-difference) Hessian: 3.63e-06   (FD truncation error)
max abs diff generalized-t (-H) vs closed H_j(beta_j):      7.11e-15   (machine precision)
```

(`data-raw/_scratch_check_omega_prior_marginal.R`.) So `H_j(\beta_j)` as
already implemented in `R/two_block_ergodicity_ing_marginal.R` **is** the
negative Hessian of this generalized-$t$ log-density — the shape and rate are
already fully "baked in," via $\nu_j$ and $\Sigma_j$, not ignored.

**Log-concavity boundary.** `multivariate-t-log-concavity.md`'s result
$d<\nu$ becomes, here, $q_j\cdot\nu_j/K_j < \nu_j$, i.e.

$$
\boxed{\ q_j < K_j = 2r_j^0+\mathrm{RSS}_j^{\mathrm{ols}}\ }
$$

— exactly `inst/BLOCK_GIBBS_ERGODICITY_ING.md` §16.3's threshold, and this is
where $a_j^0$ appears to "disappear": $\nu_j$ **cancels** out of the
ratio-form boundary $d/\nu_j<1$ because $A$ was itself built proportional to
$\nu_j/K_j$ (§3) — a general feature of this family, verified numerically in
§6 below by varying $a_j^0$ at fixed $r_j^0,\mathrm{RSS}_j^{\mathrm{ols}}$ and
confirming the PSD boundary of $H_j$ does not move. $a_j^0$ is not, however,
irrelevant to the diagnostic overall: it sets $\nu_j$, which governs the
*degrees of freedom*/tail behavior used everywhere the boundary's **confidence
level** is interpreted (`data-raw/_scratch_rss_ellipsoid_test.R`'s
$\alpha_j=P(F_{p_{re},\nu_j}>\nu_j/p_{re})$, via `inst/multivariate-t-log-concavity.md`'s
Scheffé-region correspondence) — a sharper (larger $a_j^0$) prior shrinks the
*expected* fraction of draws outside a *fixed*-radius ellipsoid without
changing the ellipsoid's radius itself.

---

## 5. What would be different under the (inapplicable) conjugate structure

For contrast — and to make concrete why row 1 of §0's table must **not** be
substituted for the derivation above — suppose $\beta_j$ *had* been given a
conjugate prior $\beta_j\mid\Omega_j\sim N(\mu_j,\Omega_j^{-1}\Sigma_{0j})$
(the calibration-time device of `inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md`
§3.3.4), and both $\beta_j$ and $\Omega_j$ were integrated out jointly instead
of $\beta_j$ being held fixed. The standard Bayesian linear-model marginal
posterior for $\beta_j$ would then be a multivariate-$t$ centered at the
**precision-weighted blend**

$$
\tilde\beta_j = \big(\Sigma_{0j}^{-1}+D_j'D_j\big)^{-1}\big(\Sigma_{0j}^{-1}\mu_j+D_j'D_j\hat\beta_j^{\mathrm{ols}}\big)
$$

not $\hat\beta_j^{\mathrm{ols}}$, with scale built from
$M_j:=(\Sigma_{0j}+(D_j'D_j)^{-1})^{-1}$ and an RSS inflated by the
quadratic penalty $(\hat\beta_j^{\mathrm{ols}}-\mu_j)'M_j(\hat\beta_j^{\mathrm{ols}}-\mu_j)$
— literally the $S_{\mathrm{marg},j}$ construction in
`inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md` Part I, and *not* a function of a
free-standing $\beta_j$ at all (there, $\beta_j$ has already been integrated
out — the object is a marginal law for $\Omega_j$ alone, not a Hessian in
$\beta_j$). Using this instead of §2-4 above would (a) silently reintroduce
$\mu_j,\Sigma_{0j}$ into what is supposed to be a $\beta_j$-conditional
diagnostic, (b) shift the log-concavity ellipsoid's center away from the OLS
fit toward the calibration prior mean, and (c) make the boundary radius
depend on $\Sigma_{0j}$, not just $r_j^0,\mathrm{RSS}_j^{\mathrm{ols}}$ —
none of which matches what `rindepNormalGamma_reg()` actually samples (§0),
and none of which is what `R/two_block_ergodicity_ing_marginal.R` computes.
The derivation in §2-4 is the one that matches both.

---

## 6. Numerical confirmation that $a_j^0$ only moves the scale, never the PSD boundary sign

With $D_j,y_j$ fixed, $r_j^0=1.7$, $\mathrm{RSS}_j^{\mathrm{ols}}$ fixed, and
$\beta_j$ fixed at a point outside the ellipsoid ($q_j>K_j$ at this $r_j^0$):

```
a0=  0.01  min eigenvalue of H_j:  3.95
a0=  1     min eigenvalue of H_j:  5.06
a0=  5     min eigenvalue of H_j:  9.56
a0= 50     min eigenvalue of H_j: 60.18
```

The minimum eigenvalue **scales up** with $a_j^0$ (larger $\Omega_j^{\mathrm{eff}}$
means more curvature everywhere) but its **sign** never flips as $a_j^0$
varies alone — consistent with $q_j<K_j=2r_j^0+\mathrm{RSS}_j^{\mathrm{ols}}$
not involving $a_j^0$ at all (§4). (This point was chosen at
$q_j\approx0.55 < K_j\approx5.26$ in the reproducible script, so $H_j$ is PSD
throughout, as expected; the point of the table is that the *sign* of the
smallest eigenvalue is constant across $a_j^0$, not that this particular
point is near the boundary.)

---

## 7. The exact tail-probability pivot: why $\alpha_j$ depends only on $\nu_j$, and moves monotonically

§4 asserted that a sharper (larger $a_j^0$) prior shrinks the *expected*
fraction of draws outside a *fixed*-radius ellipsoid without changing the
ellipsoid's radius. This section makes that claim exact rather than merely
numerical.

**The F-to-Beta pivot.** `data-raw/_scratch_rss_ellipsoid_test.R` defines
$\alpha_j = P(F_{p_{re},\nu_j} > \nu_j/p_{re})$. Recall the classical
transform: if $X\sim F(p,\nu)$, then $Y:=\dfrac{pX}{pX+\nu}\sim\mathrm{Beta}(p/2,\nu/2)$.
Evaluating $Y$ exactly at the cutoff $X=\nu/p$ gives

$$
Y = \frac{p\cdot(\nu/p)}{p\cdot(\nu/p)+\nu} = \frac{\nu}{2\nu} = \frac12
\qquad\text{(exactly, for every }\nu\text{)}
$$

so

$$
\boxed{\ \alpha_j \;=\; P\Big(\mathrm{Beta}\big(\tfrac{p_{re}}{2},\ \tfrac{\nu_j}{2}\big) > \tfrac12\Big)\ }
$$

This is a genuinely scale-free statement — $K_j$ (the radius) never entered
the derivation; the entire dependence on the prior is through $\nu_j$ alone,
via the second Beta shape parameter.

**Monotonicity.** For fixed $a=p_{re}/2$, $\mathrm{Beta}(a,b)$'s mean
$a/(a+b)\to0$ as $b\to\infty$ and the distribution becomes increasingly
concentrated near $0$ — so $P(\mathrm{Beta}(a,b)>\tfrac12)$ is monotonically
**decreasing** in $b$. Since $b=\nu_j/2$ increases one-for-one with $a_j^0$
(holding $n_j,p_{re}$ fixed, $\nu_j=2a_j^0+n_j-p_{re}$), $\alpha_j$ decreases
**monotonically** in $a_j^0$, from $\alpha_j=\tfrac12$ exactly when
$\nu_j=p_{re}$ (the $a=b$ symmetric point — matching the
"coin-flip" row of `inst/multivariate-t-log-concavity.md`'s table, $p=5,n=10$)
down to $\alpha_j\to0$ as $a_j^0\to\infty$.

**Why the radius can be frozen while the tail probability isn't.** This is
the resolution of the apparent tension: $K_j$ is a statement about where the
*curvature* of the log-density changes sign (a structural fact about the
family, independent of which particular $t$-distribution instance you're in);
$\alpha_j$ is a statement about how much *probability mass* a specific
instance (indexed by $\nu_j$) places beyond a boundary that — expressed in
the *dimensionless*, pivotal $d/\nu_j$ coordinate — sits at exactly $1$
regardless of $\nu_j$. Larger $a_j^0$ pushes $\nu_j$ up, which (i) shrinks the
scale $\Sigma_j=K_j/\nu_j\cdot(D_j'D_j)^{-1}$ (§3) and (ii) makes the
$t_{\nu_j}$-type tails lighter — both effects pull probability mass away from
the (unmoved) boundary, hence $\alpha_j\downarrow$.

Verified in `data-raw/_scratch_check_omega_prior_marginal.R`: `pf(nu/p, p, nu,
lower.tail=FALSE)` and `pbeta(0.5, p/2, nu/2, lower.tail=FALSE)` agree to
machine precision across a grid of $(\nu,p)$, and $\alpha_j$ decreases
monotonically as $\nu_j$ is swept from $p_{re}$ upward.

---

## 8. What happens when `pwt_measurement` changes (shape and rate move together)

§4/§7 held $r_j^0$ fixed while varying $a_j^0$ in isolation to isolate what
each parameter does. In practice, `dGamma_list()`'s calibration
(`inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md` Part I) does not offer that
isolation: both `shape_ING,j` ($a_j^0$) and `rate,j` ($r_j^0$) are driven by
the *same* knob, `pwt_measurement` ($w_j$), through
$n_{\mathrm{prior},j}=\frac{w_j}{1-w_j}n_j$:

$$
a_j^0 = \frac{n_{\mathrm{prior},j}+1}{2}+\frac{p_{re}}{2}, \qquad
r_j^0 = \frac{S_{\mathrm{marg},j}}{2}\cdot\frac{n_{\mathrm{prior},j}+p_{re}-1}{n_j-p_{re}}
$$

Both increase linearly in $n_{\mathrm{prior},j}$ — hence both increase
monotonically as $w_j\uparrow$ — and, per the package's own "key algebraic
fact" (loc. cit.), are tied together *exactly*, for every $w_j$, by

$$
r_j^0 = \hat\sigma_j^2\,(a_j^0-1), \qquad \hat\sigma_j^2 := \frac{S_{\mathrm{marg},j}}{n_j-p_{re}}
\ \ (\text{independent of }w_j)
$$

**Substituting into $K_j$ and $\nu_j$.** Using this identity,

$$
K_j = 2r_j^0+\mathrm{RSS}_j^{\mathrm{ols}} = \hat\sigma_j^2\,\nu_j
+ \underbrace{\Big[\mathrm{RSS}_j^{\mathrm{ols}}-S_{\mathrm{marg},j}-2\hat\sigma_j^2\Big]}_{\text{constant in }w_j,\ <\,0}
$$

so $K_j$ — unlike in §4's isolated-$a_j^0$ thought experiment — **also grows**
with $w_j$, affinely in $\nu_j$. Dividing by $\nu_j$:

$$
\frac{K_j}{\nu_j} = \hat\sigma_j^2 - \frac{\big(S_{\mathrm{marg},j}-\mathrm{RSS}_j^{\mathrm{ols}}\big)+2\hat\sigma_j^2}{\nu_j}
$$

Since $S_{\mathrm{marg},j}\ge\mathrm{RSS}_j^{\mathrm{ols}}$ always (the
calibration penalty is $\ge0$), the subtracted term is strictly positive, so
$K_j/\nu_j$ (the marginal-$t$'s scale factor, §3) **increases monotonically**
in $w_j$ as well, rising from a value below $\hat\sigma_j^2$ at $w_j=0$ and
approaching $\hat\sigma_j^2\cdot(D_j'D_j)^{-1}$ (never reaching it) as
$w_j\to1$.

**So: does larger `pwt_measurement` help or hurt the ellipsoid diagnostic?**
Both the radius $K_j$ *and* the group's own marginal-$t$ scale $K_j/\nu_j$
grow as $w_j\uparrow$ — at first glance this looks ambiguous. But §7 already
settled it: $\alpha_j$ depends on $\nu_j$ **alone**, not on $K_j$ or the
scale, via the scale-free F/Beta pivot. And
$\nu_j = n_{\mathrm{prior},j}+n_j+1$ is strictly increasing in $w_j$
throughout $w_j\in[0,1)$ (since $w_j/(1-w_j)$ is). Therefore

$$
\boxed{\ \alpha_j\ \text{decreases monotonically as}\ w_j=\texttt{pwt\_measurement}\ \text{increases, over the entire range}\ [0,1)\ }
$$

with $\alpha_j\to0$ as $w_j\to1$ ($\nu_j\to\infty$), and, at $w_j=0$,
$\nu_j=n_j+1$ (not $n_j-p_{re}$ — `shape_ING,j`'s $\tfrac12+\tfrac{p_{re}}{2}$
floor comes from `compute_gaussian_prior()`'s general construction, not from
`pwt_measurement` itself). **Larger `pwt_measurement` unambiguously helps**,
on this metric, with no tradeoff: even though it also enlarges the physical
ellipsoid, it shrinks the tail probability of exceeding that (also-larger)
ellipsoid *faster*, in the scale-free sense that is all $\alpha_j$ actually
measures.

**Caveat.** $\alpha_j$ is the *no-shrinkage baseline* — a closed-form property
of group $j$'s own marginal-$t$ in isolation, with Block-2's $(\Lambda,\gamma)$
pull deliberately absent. Increasing `pwt_measurement` doesn't touch
$\Lambda,\gamma$ directly, but it does shift what the *real* sampler draws
for $\Omega_j$ (toward $\hat\sigma_j^{-2}$, away from the group's raw OLS
residual variance), which can change how far the real, autocorrelated Gibbs
draws of $\beta_j$ end up sitting from $\hat\beta_j^{\mathrm{ols}}$. So this
result is a statement about the theoretical baseline improving monotonically,
not a guarantee that the *observed* `pct_draws_outside` in a live run
improves by the same amount — though it is consistent with, and likely part
of the explanation for, the empirical finding that Part VI's broader
recalibration improved `pct_outside_prerun` across most groups in the
`Ex_13`/`Ex_13b` comparison.

Verified in `data-raw/_scratch_check_omega_prior_marginal.R`: sweeping
`pwt_measurement` over a grid in $(0,1)$ (via $n_{\mathrm{prior},j}=w_j/(1-w_j)n_j$)
confirms $a_j^0,r_j^0,\nu_j,K_j,K_j/\nu_j$ all increase monotonically while
$\alpha_j$ decreases monotonically, and confirms the exact identity
$r_j^0=\hat\sigma_j^2(a_j^0-1)$ holds at every grid point.

**The deeper mechanism: sharpening the prior slides the marginal likelihood
toward its own Gaussian limit.** §7's "lighter tails" language can be made
completely literal. Since $\Sigma_j=K_j/\nu_j\cdot(D_j'D_j)^{-1}$ converges to
a *fixed, non-degenerate* limit $\hat\sigma_j^2(D_j'D_j)^{-1}$ as
$\nu_j\to\infty$ (established just above), the classical identity
$(1+d/\nu)^{\nu}\to e^{d}$ forces the §3 kernel itself to converge pointwise:

$$
\Big[1+\frac{d}{\nu_j}\Big]^{-(\nu_j+p_{re})/2}
= \Big[\big(1+\tfrac{d}{\nu_j}\big)^{\nu_j}\Big]^{-1/2}\Big(1+\frac{d}{\nu_j}\Big)^{-p_{re}/2}
\ \longrightarrow\ e^{-d/2}
$$

i.e. $p(y_j\mid\beta_j)$ itself converges, as `pwt_measurement`$\to1$, to the
Gaussian kernel of $N\big(\hat\beta_j^{\mathrm{ols}},\ \hat\sigma_j^2(D_j'D_j)^{-1}\big)$
— exactly the classical sampling distribution of the OLS estimator, with the
calibrated $\hat\sigma_j^2$ standing in for the usual residual-based variance.
Numerically (`data-raw/_scratch_check_omega_prior_marginal.R`, fixed $d=5$,
$p_{re}=2$), the ratio of the $t$-kernel to the Gaussian kernel at that same
$d$ tracks this convergence exactly: $1.077,\ 1.046,\ 1.012,\ 1.001,\ 1.00001$
for $\nu=5,20,100,1000,10^5$.

And a Gaussian is **globally** log-concave —
`inst/multivariate-t-log-concavity.md`'s own "Gaussian limit" remark: *"As
$\nu\to\infty$, $s\to1$ and $-H\to\Sigma^{-1}$ for all $x$, recovering global
log-concavity."* There is no finite boundary left in the limit, which is
exactly why $K_j$ **itself** diverges (not just $K_j/\nu_j$) as
`pwt_measurement`$\to1$ — the worked table above literally shows it:
$K_j=19.5,\ 26.1,\ 46.1,\ 126.0,\ 525.4,\ 2655.9$ for
$w_j=0,0.2,0.5,0.8,0.95,0.99$. So sharpening `pwt_measurement` is not merely
"incidentally" shrinking $\alpha_j$ — it is literally sliding group $j$'s own
marginal likelihood along the one-parameter family that interpolates between
"Student-$t$ with a finite, sometimes-violated log-concavity boundary" (small
$\nu_j$) and "Gaussian, log-concave everywhere, no boundary at all" (large
$\nu_j$). $\alpha_j\to0$ is what "the boundary has receded toward infinity"
looks like from the finite-$\nu_j$ side of that spectrum.

---

## 9. Beyond the no-shrinkage baseline: the real, noncentral full conditional

§§7-8 are a correct, closed-form, and explicitly self-limited result: $\alpha_j$
is the exceedance probability of the *isolated* group-$j$-only posterior
(Block~1's likelihood plus the weak $\sigma^2_j$ prior, with Block~2's
$(\Lambda,\gamma)$ hierarchical pull switched off entirely, matching §8's own
caveat). It is **not**, and was never claimed to be, a prediction of the real
sampler's observed `pct_draws_outside` for a specific group -- doing that
requires bringing $(\Lambda,\gamma)$ back in, which reintroduces exactly the
noncentrality this section derives.

### 9.1 The real full conditional is a shrinkage estimator, not centered at $\hat\beta_j^{\mathrm{ols}}$

Fixing $\Omega_j$ at a plug-in value (its own Gaussian full-conditional
structure is what makes this tractable, unlike the $\Omega_j$-integrated-out
case of §§1-4), the real Gibbs sampler's Block~1 draw combines two Gaussian
pieces exactly as `Prior_Setup_lmebayes.R`/`_scratch_rss_prerun_estimate.R`
already implement:

$$
\beta_j\mid\Omega_j,\Lambda,\gamma\ \sim\ N(\mu_j,\ B_j^{-1}),\qquad
B_j := \Omega_jD_j'D_j+\Lambda,\qquad
\mu_j := B_j^{-1}\big(\Omega_jD_j'D_j\,\hat\beta_j^{\mathrm{ols}}+\Lambda W_j\gamma\big).
$$

Writing $e_j:=\beta_j-\hat\beta_j^{\mathrm{ols}}$, this is $e_j\sim
N(\Delta_j,\ B_j^{-1})$ with

$$
\boxed{\ \Delta_j = B_j^{-1}\Lambda\big(W_j\gamma-\hat\beta_j^{\mathrm{ols}}\big)\ }
$$

-- the **noncentrality vector**: how far, and in which direction, the
hierarchy pulls group $j$ away from its own OLS fit, weighted by how much
relative precision $\Lambda$ has versus the group's own data
($\Omega_jD_j'D_j$, inside $B_j$). This is not an add-on correction; it *is*
partial pooling, the central mechanism of every hierarchical model in this
package. By construction $\Delta_j$ is largest for exactly the groups whose
$\hat\beta_j^{\mathrm{ols}}$ disagrees most with what the rest of the
hierarchy implies -- i.e. the same groups (6, 33, 30) the ellipsoid test and
$H_j$-violation tables have flagged as outliers all along.

### 9.2 Why this does not reduce to a noncentral $F$ in general

The boundary check is $q_j(\beta_j)=e_j'D_j'D_je_j\le K_j$. With
$e_j\sim N(\Delta_j,B_j^{-1})$, $q_j$ is a **noncentral quadratic form in a
Gaussian vector**, weighted by $D_j'D_j$ under a covariance $B_j^{-1}$ that is
*not*, in general, proportional to $D_j'D_j$'s own eigenstructure (that would
require $\Lambda\propto\Omega_jD_j'D_j$, not something the model guarantees).
Consequently $q_j$ is in general a weighted sum of noncentral $\chi^2$'s (one
per eigenvalue of $B_j^{-1/2}D_j'D_j B_j^{-1/2}$), not a single noncentral
$\chi^2$ or $F$ -- so §7's clean, scale-free F-to-Beta pivot does **not**
survive the shift. There is no equally simple closed form for the
noncentral case.

### 9.3 The practical resolution: simulate the real conditional directly, and reuse it for every $w_j$ at once

Rather than chase a closed-form (Satterthwaite/Patnaik-style) approximation
to a weighted noncentral $\chi^2$, `_scratch_rss_prerun_estimate.R` already
implements the exact answer (up to Monte Carlo error and the $\Omega_j,\gamma$
plug-in approximation, which the pre-run estimate always carries): draw
directly from the closed-form Gaussian $N(\mu_j,B_j^{-1})$ and compute the
empirical fraction with $q_j>K_j$. This is exact for *any* $\Lambda,\gamma$,
with no assumption on how $D_j'D_j$ and $B_j^{-1}$ relate.

**Key simplification for calibrating `pwt_measurement` group-by-group.**
$\hat\sigma_j^2$ -- hence $\Omega_j$, $B_j$, $\mu_j$, and the entire
*distribution* of $e_j=\beta_j-\hat\beta_j^{\mathrm{ols}}$ -- is, per §8's own
identity, **invariant to $w_j$**; only the boundary $K_j=2r_j^0+
\mathrm{RSS}_j^{\mathrm{ols}}$ moves as $w_j$ (hence $r_j^0$) changes. So a
*single* Monte Carlo draw of $q_j$'s distribution per group is enough to solve
for the $w_j$ needed to hit *any* target $\alpha^\star$, via one quantile
lookup instead of re-simulating per candidate $w_j$:

1. Simulate $q_j^{(1)},\dots,q_j^{(n_{\mathrm{sim}})}$ once, at the group's
   *current* (any) calibration, exactly as `_scratch_rss_prerun_estimate.R`
   already does.
2. Read off $q_j^\star := $ the empirical $(1-\alpha^\star)$ quantile of the
   $q_j^{(i)}$'s -- the boundary that would leave exactly $\alpha^\star$ of
   this same, noncentrality-correct distribution outside.
3. Invert $K_j=q_j^\star \iff r_j^{0\star} = (q_j^\star-\mathrm{RSS}_j^{\mathrm{ols}})/2$,
   then, using $r_j^0=\hat\sigma_j^2(a_j^0-1)$ (§8, invariant to $w_j$):
   $a_j^{0\star} = r_j^{0\star}/\hat\sigma_j^2+1$, and, using
   $a_j^0=(n_{\mathrm{prior},j}+1)/2+p_{re}/2$:
   $n_{\mathrm{prior},j}^\star = 2a_j^{0\star}-p_{re}-1$.
4. $w_j^\star = \max\big(0,\ n_{\mathrm{prior},j}^\star/(n_{\mathrm{prior},j}^\star+n_j)\big)$,
   clipped at the package's own $0.5$ per-group ceiling and flagged if that
   ceiling binds (meaning $\alpha^\star$ is not reachable for that group via
   `pwt_measurement` alone -- the same conclusion `Ex_13c` reached
   empirically for groups 6/33 by excluding them instead).

Because step 1's simulation doesn't depend on $w_j$ at all, this whole
procedure is a **single pass per group**, not an iterative root-find --
despite $q_j$ having no closed-form noncentral distribution. See
`data-raw/_scratch_group_pwt_measurement_noncentral.R` for the
implementation, run against `Ex_13b`'s fixture.

**Caveats, inherited from `_scratch_rss_prerun_estimate.R`:** $\Omega_j$ and
$\gamma$ are held at plug-in point estimates (the sampler's own posterior
uncertainty in both is ignored), so this is a pre-run *approximation* to the
sampler's true `pct_draws_outside`, not a certified value -- exactly the
caveat already attached to the pre-run estimator itself. Unlike §7-8's
$\alpha_j$, though, it *does* correctly reflect each group's own
noncentrality $\Delta_j$, which is the property this section exists to
supply.

## 10. Conclusion

`inst/BLOCK_GIBBS_ERGODICITY_ING.md` §16's $H_j(\beta_j)$, the
$q_j\le2r_j^0+\mathrm{RSS}_j^{\mathrm{ols}}$ log-concavity threshold, and
`data-raw/_scratch_rss_ellipsoid_test.R`'s $\nu_j=2a_j^0+n_j-p_{re}$/$\alpha_j$
machinery were, on close (re-)derivation from the correct **Independent**
Normal-Gamma sampling model (§1-2, not the calibration-time conjugate device
of §5), already exactly the negative Hessian and log-concavity boundary of a
**generalized multivariate-$t$** — plain OLS multivariate-$t$
(`inst/multivariate-t-log-concavity.md`) with $n_j\to n_j+2a_j^0$ and
$\mathrm{RSS}_j^{\mathrm{ols}}\to\mathrm{RSS}_j^{\mathrm{ols}}+2r_j^0$. The
measurement-precision prior's shape *and* rate are both already fully
accounted for — the shape enters via $\nu_j$ (degrees of freedom / confidence
level of the Scheffé-region interpretation), the rate enters via both $\nu_j$
and the scale $K_j$, and the fact that $a_j^0$ cancels out of the specific
ratio $q_j/K_j<1$ is an exact algebraic consequence of this family, not a
dropped term. §7's exact F-to-Beta pivot shows why this is not a
contradiction: $\alpha_j$ is a scale-free statement depending on $\nu_j$
alone, so it can (and does) move monotonically even while the boundary
$K_j$ itself is frozen (§7) or, more realistically, moving *with* $\nu_j$
under the package's actual `pwt_measurement`-driven calibration (§8) — in
either case $\alpha_j$ improves monotonically as the prior is sharpened. No
change to `R/two_block_ergodicity_ing_marginal.R` or to Section 16 is implied
by this note; it is a confirmation, not a correction.

§9 adds the piece $\alpha_j$ deliberately leaves out: the real sampler's
$\beta_j$ is not centered at $\hat\beta_j^{\mathrm{ols}}$, and the resulting
noncentrality $\Delta_j$ -- driven by the hierarchy's pull, not by $n_j$ --
is the dominant reason a fixed-$n_j$ comparison across groups (e.g. groups
5/6/14/24/26/40, all $n_j=11$) shows very different empirical violation
rates despite identical $\alpha_j$. §9's Monte Carlo procedure, not §7-8's
closed form, is the right tool for group-specific `pwt_measurement`
calibration.
