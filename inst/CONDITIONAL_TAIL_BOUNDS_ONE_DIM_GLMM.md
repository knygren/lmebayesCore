# Bounds on conditional tail probabilities — one-dimensional GLMMs

Per-group, drift-free **upper bounds** on
\(\pi\bigl(\beta_j \in E_j(\omega_j) \mid \gamma, y_j\bigr)\)
for the package's one-dimensional group RE models (Gaussian, logit, probit,
Poisson, cloglog, Gamma regression). These are the building blocks for
`LOGIT_STATIC_TAIL_CERTIFICATION.md` (profile \(h\), gates \(\Gamma_{\mathrm{safe}}\),
\(\varepsilon_2(\boldsymbol{\omega})\), \(B(\boldsymbol{\omega})\)) and for
Claim 1 in `LOGIT_MARGINAL_INTEGRATE_GAMMA.md`.

**Companion notes.**

| Note | Role |
|---|---|
| `notation.md` | \(D_j\), \(\mathcal{W}_j\), \(\Psi\), \(P_b=\Psi^{-1}\) |
| `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0–§6 | Likelihood specializations (source proofs) |
| `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §IV | \(A/A^c\), \(\omega_j(\beta)\), box certification |
| `LOGIT_STATIC_TAIL_CERTIFICATION.md` | Multi-group union, \(\gamma\) partition, §4.5–§4.6 |
| `LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` | Template proof (logit \(J=1\)) |

**Scope.** Scalar intercept RE per group (\(p_{\mathrm{re}}=1\)):
\(\eta_{j,i}=\mathrm{offset}_{j,i}+x_j\beta_j\),
\(\beta_j \mid \gamma \sim N(x_j\gamma,\,\lambda_b^{-1})\).
Multivariate RE (\(q>1\)): same lemmas apply group-wise once \(E_j(\omega_j)\) is
defined from \(\lambda_{\max}(\Theta_j(\beta_j))\); margin identities (§1.4) are
intercept-specific.

---

## 1. Common notation (all likelihoods)

### 1.1 Hierarchy and conditional posterior

Groups \(j=1,\ldots,J\). **Population** \(\gamma\in\mathbb R^q\) (often \(q=1\)).
**Group RE** with precision \(P_b=\lambda_b I\) in the scalar case:

\[
\beta_j \mid \gamma \;\sim\; N(x_j\gamma,\,\lambda_b^{-1}),
\qquad
\pi(\beta \mid \gamma, y)
=
\prod_{j=1}^{J}\pi(\beta_j \mid \gamma, y_j),
\]

\[
\pi(\beta_j \mid \gamma, y_j)
\propto
\exp\bigl(\ell_j(\beta_j)\bigr)\,
\phi\bigl(\beta_j;\,x_j\gamma,\,\lambda_b^{-1}\bigr).
\]

Group log-likelihood \(\ell_j(\beta_j)=\sum_{i=1}^{n_j}\ell_{j,i}(\eta_{j,i})\),
\(\eta_{j,i}=\mathrm{offset}_{j,i}+x_j\beta_j\).

*(Package mapping: level-2 design \(\mathcal{W}_j\equiv H_j\); scalar intercept
\(H_j=x_j\). Do not confuse data curvature sum \(W_j^{\mathrm{sum}}(\beta)\) below
with \(\mathcal{W}_j\) in `notation.md`.)*

### 1.2 Shrinkage and certification threshold

IRLS / expected Fisher weight \(w_{j,i}(\eta)\) (family-specific; §2–§7).
Curvature sum and shrinkage:

\[
W_j^{\mathrm{sum}}(\beta_j)
:=
\sum_{i=1}^{n_j} w_{j,i}\bigl(\eta_{j,i}(\beta_j)\bigr),
\qquad
\omega_j(\beta_j)
:=
\frac{\lambda_b}{\lambda_b + W_j^{\mathrm{sum}}(\beta_j)}.
\]

Fix a **certification vector** \(\boldsymbol{\omega}=(\omega_1,\ldots,\omega_J)\),
\(\omega_j\in(0,1]\) (threshold; not the state function \(\omega_j(\beta_j)\)).

**Per-group tail event** (primary definition):

\[
\boxed{
E_j(\omega_j)
:=
\{\,\beta_j\in\mathbb R : \omega_j(\beta_j)\ge \omega_j\,\}.
}
\]

**Conditional tail mass** (object to bound):

\[
\boxed{
g_j(\gamma;\omega_j)
:=
\pi\bigl(\beta_j\in E_j(\omega_j)\mid\gamma,y_j\bigr).
}
\]

**Target (Proposition 1 template, all GLMs):** find explicit \(G_j(\gamma;\omega_j)\) with

\[
\boxed{
g_j(\gamma;\omega_j)\;\le\;G_j(\gamma;\omega_j).
}
\]

Multi-group profile (independence given \(\gamma\)):

\[
g(\gamma;\boldsymbol{\omega})
=
1-\prod_{j=1}^{J}\bigl(1-g_j(\gamma;\omega_j)\bigr),
\qquad
h(\gamma;\boldsymbol{\omega})
=
1-\prod_{j=1}^{J}\bigl(1-G_j(\gamma;\omega_j)\bigr)
\;\ge\;
g(\gamma;\boldsymbol{\omega}).
\]

### 1.3 Shared lemmas and Proposition 1 template

Write \(\sigma_j := \lambda_b^{-1/2}\). When \(\ell_j''\le 0\) and
\(|\ell_j'(\beta_j)|\le L_j\) for **every** \(\beta_j\in E_j(\omega_j)\):

**Lemma T.** \(\ell_j(\beta_j)\ge \ell_j(x_j\gamma)-L_j|\beta_j-x_j\gamma|\).

**Lemma G (upper tail).**
\(\int_t^\infty e^{cz}\phi(z;x_j\gamma,\sigma_j^2)\,dz
= e^{cx_j\gamma+c^2\sigma_j^2/2}\Phi(-(t-x_j\gamma-c\sigma_j^2)/\sigma_j)\).

**Lemma G′ (lower tail).**
\(\int_{-\infty}^t e^{cz}\phi(z;x_j\gamma,\sigma_j^2)\,dz
= e^{cx_j\gamma+c^2\sigma_j^2/2}\Phi((t-x_j\gamma-c\sigma_j^2)/\sigma_j)\).

**Lemma Z.** With
\(Z_j(\gamma):=\int e^{\ell_j(\beta_j)}\phi(\beta_j;x_j\gamma,\lambda_b^{-1})\,d\beta_j\),

\[
Z_j(\gamma,y_j)
\ge
\exp\!\Big(\ell_j(x_j\gamma) - L_j\,\sigma_j\sqrt{2/\pi}\Big).
\]

*Proofs:* `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0.3;
`LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md`.

**Proposition 1 template.** If \(E_j(\omega_j)\) is a one- or two-sided tail and
\(L_j\) certifies Lemma T on **all** of \(E_j(\omega_j)\), then

\[
g_j(\gamma;\omega_j)
=
\frac{\int_{E_j(\omega_j)} e^{\ell_j}\phi\,d\beta_j}{Z_j(\gamma)}
\le
e^{L_j\sigma_j\sqrt{2/\pi}+L_j^2/(2\lambda_b)}
\times (\text{tail sum of }\Phi\text{-terms}).
\]

*Proof.* Lemma T with \(b_0=x_j\gamma\); majorize the numerator by tilted-Gaussian
integrals on each tail piece (Lemmas G/G′); divide by Lemma Z. \(\square\)

**One-sided bound** (integrate over \(\beta_j<t_j\) with \(E_j(\omega_j)\subseteq\{\beta_j<t_j\}\)):

\[
G_j^{\mathrm{low}}(\gamma;\omega_j)
=
\exp\!\Big(L_j\,\sigma_j\sqrt{2/\pi} + \tfrac{L_j^2}{2\lambda_b}\Big)\,
\Phi\!\Big(-\tfrac{t_j - x_j\gamma - L_j/\lambda_b}{\sigma_j}\Big),
\qquad
t_j=\beta_j^\star(\omega_j).
\]

**Two-sided bound** (logit; \(E_j=\{|\beta_j|>b_j^\star(\omega_j)\}\)):

\[
G_j^{\mathrm{2s}}(\gamma;\omega_j)
=
\exp\!\Big(L_j\,\sigma_j\sqrt{2/\pi} + \tfrac{L_j^2}{2\lambda_b}\Big)
\left[
\Phi\!\Big(-\tfrac{b_j^\star(\omega_j)-x_j\gamma-L_j/\lambda_b}{\sigma_j}\Big)
+
\Phi\!\Big(-\tfrac{b_j^\star(\omega_j)+x_j\gamma-L_j/\lambda_b}{\sigma_j}\Big)
\right].
\]

*(Likelihood-specific inputs: \(L_j(\omega_j)\), margin \(\beta_j^\star(\omega_j)\)
or \(b_j^\star(\omega_j)\), and verification that \(E_j(\omega_j)\) matches the
tail used in the \(\Phi\) formula.)*

### 1.4 Margin from \(\omega_j\)

Always start from the **primary** definition \(E_j(\omega_j)=\{\omega_j(\beta_j)\ge\omega_j\}\)
and the curvature floor

\[
\boxed{
\kappa_j(\omega_j)
:=
\frac{\lambda_b\,(1-\omega_j)}{\omega_j},
\qquad
E_j(\omega_j)
=
\{\,\beta_j : W_j^{\mathrm{sum}}(\beta_j)\le \kappa_j(\omega_j)\,\}
\quad\text{(when equivalent)}.
}
\]

#### 1.4.1 Logit (two-sided, intercept, offsets \(0\))

\(W_j^{\mathrm{sum}}(\beta_j)=N_j\,u(x_j\beta_j)\), \(u(\eta)=p(\eta)(1-p(\eta))\),
strictly decreasing in \(|\eta|\). Define \(\eta_j^\star(\omega_j)\ge 0\) by
\(N_j u(\eta_j^\star(\omega_j))=\kappa_j(\omega_j)\), then
\(b_j^\star(\omega_j):=\eta_j^\star(\omega_j)/|x_j|\) and
\(E_j(\omega_j)=\{|\beta_j|>b_j^\star(\omega_j)\}\).
*(Lemma 3: `LOGIT_STATIC_TAIL_CERTIFICATION.md` §2.)*

#### 1.4.2 Poisson (one-sided low-curvature tail)

\(w_{j,i}(\eta)=e^{\eta}\). On the **unsafe tail** \(E_j(\omega_j)\), each
\(\eta_{j,i}(\beta_j)\) is bounded above by the margin value at the boundary
\(\beta_j^\star(\omega_j)\) solving

\[
W_j^{\mathrm{sum}}(\beta_j^\star(\omega_j))=\kappa_j(\omega_j),
\qquad
\eta_{j,i}^\star(\omega_j):=\eta_{j,i}(\beta_j^\star(\omega_j)).
\]

Define the **score bound on all of \(E_j(\omega_j)\)**:

\[
\boxed{
L_j^{\mathrm{po}}(\omega_j)
:=
\sum_{i=1}^{n_j}\max\bigl(y_{j,i},\,e^{\eta_{j,i}^\star(\omega_j)}\bigr).
}
\]

(On \(E_j\), \(\eta_{j,i}\le \eta_{j,i}^\star\) when \(W_j^{\mathrm{sum}}\) is
increasing in \(\eta_{j,i}\) along the tail — intercept with common \(x_j>0\),
offsets fixed.) Then \(E_j(\omega_j)=\{\beta_j:\beta_j\le\beta_j^\star(\omega_j)\}\)
(one-sided; flip inequalities if \(x_j<0\)).

**Do not confuse** with the **safe-box floor** \(\eta\ge\log s_j\) in
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §4.3, used only to certify \(B\subseteq A\)
(high-\(\eta\) working region). Proposition Po1 integrates over the **full**
low-\(W_j^{\mathrm{sum}}\) tail \(E_j(\omega_j)\), not a truncated strip.

#### 1.4.3 Probit (two-sided low-\(W\) tail)

Per-observation Fisher weight (package convention):

\[
w_{j,i}(\eta)=\frac{n_{j,i}\,\varphi(\eta)^2}{\bigl[\Phi(\eta)\,\Phi(-\eta)\bigr]^2},
\qquad
W_j^{\mathrm{sum}}(\beta_j)=\sum_i w_{j,i}(\eta_{j,i}(\beta_j)).
\]

**Both tails.** \(w_{j,i}(\eta)\) is **symmetric** in \(\eta\) and vanishes as
\(|\eta|\to\infty\) (exponentially fast): equivalently, the probit log-likelihood
Hessian \(-\partial^2\ell/\partial\beta^2 = x^2\cdot\lambda(\eta)[\lambda(\eta)+\eta]\)
with inverse Mills ratio \(\lambda\) goes to \(0\) at \(z\to\pm\infty\). The global
**upper bound** \(w_{j,i}\le (2/\pi)\,n_{j,i}\) is a **maximum at \(\eta=0\)**, not
a lower bound at \(+\infty\).

Define \(\eta_j^\star(\omega_j)>0\) from the margin
\(W_j^{\mathrm{sum}}(\beta_j^\star)=\kappa_j(\omega_j)\) (intercept, offsets \(0\):
solve \(\sum_i n_{j,i}\,w_{j,i}(\eta_j^\star)/[\cdots]=\kappa_j\), or common-\(\eta\)
reduction). Set \(b_j^\star(\omega_j):=\eta_j^\star(\omega_j)/|x_j|\). Then

\[
E_j(\omega_j)=\{\,|\beta_j|>b_j^\star(\omega_j)\,\}
\]

(same geometry as logit §1.4.1; only the link function \(w(\eta)\) differs).
Score: \(L_j^{\mathrm{pr}}=\sum_i n_{j,i}(2/\pi)\) (valid on all of \(\mathbb R\)).

#### 1.4.4 Poisson (one-sided low-\(W\) tail)

Same pattern as §1.4.2: \(W_j^{\mathrm{sum}}\to 0\) only as \(\eta\to-\infty\);
the **upper** \(\eta\)-tail has \(w=e^\eta\to\infty\). Hence \(E_j(\omega_j)\) is
**one-sided** in \(\beta_j\) (§5).

#### 1.4.5 Cloglog (one-sided **observed** flattening; Fisher-weight trap)

Per-observation **Fisher / IRLS weight** (package convention; enters \(\omega_j\)):

\[
\mu_{j,i}(\eta)=1-\exp\bigl(-\exp(\eta)\bigr),
\qquad
w_{j,i}(\eta)=n_{j,i}\,\frac{\exp(2\eta)\,(1-\mu_{j,i})}{\mu_{j,i}}.
\]

Numerically, \(w(\eta)\to 0\) as \(\eta\to-\infty\) **and** as \(\eta\to+\infty\)
(since \(w\sim n\,\exp(2\eta-\exp(\eta))\to 0\); the old claim \(w\to n_i\) at
\(+\infty\) used the wrong link \(1-\mu\sim e^{-\eta}\)).

**Observed log-likelihood Hessian (Prop 1 geometry).** This is **not** the Mills-ratio
pattern. For a common linear predictor \(\eta\), aggregate
\(\mathcal H_{\mathrm{total}}=n_1\mathcal H_{y=1}+n_0\mathcal H_{y=0}\) with
\(n_1=\sum_i y_{j,i}\), \(n_0=\sum_i(n_{j,i}-y_{j,i})\). Per Bernoulli component,

\[
\frac{\partial^2 \ell_{y=0}}{\partial \eta^2}=-e^{\eta},
\qquad
\frac{\partial^2 \ell_{y=1}}{\partial \eta^2}\to 0
\quad\text{as }\eta\to+\infty.
\]

Hence, when **both** \(n_1>0\) and \(n_0>0\):

- \(\eta\to-\infty\): \(\mathcal H_{\mathrm{total}}\to 0\) (both components vanish);
- \(\eta\to+\infty\): \(\mathcal H_{\mathrm{total}}\sim n_0(-e^{\eta})\to -\infty\)
  (failure term dominates; **not** flat).

**Do not conflate Fisher weight with observed curvature.** \(w\to 0\) at \(+\infty\)
while \(\mathcal H_{\mathrm{total}}\to -\infty\) when failures are present; the score
also grows (\(|\partial\ell/\partial\eta|\lesssim n_0 e^{\eta}\) on the upper tail).
Prop 1 therefore uses the **one-sided lower tail** (same Poisson-style template as
§1.4.2), not the Fisher-weight outer set \(\{W_j^{\mathrm{sum}}\le\kappa\}\), which
would incorrectly include an upper outer interval.

On the **lower** tail \(E_j(\omega_j)\), with margin
\(\eta_{j,i}^\star(\omega_j)\) at \(\beta_j^\star(\omega_j)\) solving
\(W_j^{\mathrm{sum}}(\beta_j^\star)=\kappa_j(\omega_j)\) on the \(\eta\to-\infty\) branch,

\[
\boxed{
L_j^{\mathrm{cl}}(\omega_j)
:=
\sum_{i=1}^{n_j}\max\bigl(y_{j,i},\,e^{\eta_{j,i}^\star(\omega_j)}\bigr).
}
\]

*(All-success groups \(n_0=0\) are a separate edge case: observed Hessian can also
flatten on the upper tail; not covered by the mixed-data Prop 1 bound above.)*

### 1.5 Pipeline after Proposition 1

Once \(G_j(\gamma;\omega_j)\) is available:

1. Gates \(\Gamma_{\mathrm{low},j}(\omega_j)=\{\gamma:G_j(\gamma;\omega_j)\le\varepsilon_{1,j}\}\);
   invert \(G_j\) (bisection) for cutoff radii.
2. \(\varepsilon_2(\boldsymbol{\omega})=\pi(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})^c\mid y)\) via
   \(\gamma\)-marginal quadrature (`LOGIT_STATIC_TAIL_CERTIFICATION.md` §4.5).
3. \(B(\boldsymbol{\omega})=\varepsilon_2(\boldsymbol{\omega})+(1-\varepsilon_2(\boldsymbol{\omega}))\,\varepsilon_1^{\mathrm{eff}}\).

---

## 2. Gaussian (Normal regression)

### 2.1 Specialization

\(w_{j,i}=1/\phi\), \(W_j^{\mathrm{sum}}=n_j/\phi\), **constant in \(\beta_j\)**:

\[
\omega_j(\beta_j)
\equiv
\omega_j^{\mathrm{Ga}}
=
\frac{\lambda_b}{\lambda_b+n_j/\phi}.
\]

### 2.2 Proposition N1 (trivial safe region)

Because \(\omega_j(\beta_j)\) does not depend on \(\beta\),

\[
E_j(\omega_j)=\varnothing
\quad\text{if}\quad
\omega_j^{\mathrm{Ga}}<\omega_j,
\qquad
E_j(\omega_j)=\mathbb R
\quad\text{if}\quad
\omega_j^{\mathrm{Ga}}\ge \omega_j.
\]

Hence, on a **certified safe** choice \(\omega_j^{\mathrm{Ga}}<\omega_j\),

\[
\boxed{
g_j(\gamma;\omega_j)=0,
\qquad
G_j(\gamma;\omega_j)\equiv 0.
}
\]

**No Lemmas T/G/Z needed** — \(\pi(\beta_j\mid\gamma,y_j)\) is Gaussian; once
\(\boldsymbol{\omega}\) is set so \(A^c(\boldsymbol{\omega})=\varnothing\), tail
mass is **exactly zero** for every \(\gamma\).

*Reference:* `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §1.2.

---

## 3. Logit (binomial, canonical link)

### 3.1 Specialization

\(w_{j,i}=n_{j,i}\,p_{j,i}(1-p_{j,i})\),
\(L_j=\sum_i\max(y_{j,i},n_{j,i}-y_{j,i})\le N_j\),
\(\ell_j''\le 0\). **Two-sided** tails; margin \(b_j^\star(\omega_j)\) from §1.4.

### 3.2 Proposition L1 (Proposition 1)

Under §1.4.1 (\(\mathrm{offset}_{j,i}\equiv 0\), \(x_j\neq 0\)), Lemma T holds on
\(\mathbb R\) with \(L_j=\sum_i\max(y_{j,i},n_{j,i}-y_{j,i})\). On
\(E_j(\omega_j)=\{|\beta_j|>b_j^\star(\omega_j)\}\),

\[
\boxed{
g_j(\gamma;\omega_j)
\le
G_j^{\mathrm{logit}}(\gamma;\omega_j)
:=
G_j^{\mathrm{2s}}(\gamma;\omega_j)
}
\]

*Proof.* Proposition 1 template: split \(E_j\) into upper and lower tails at
\(\pm b_j^\star(\omega_j)\); apply Lemmas G/G′ with \(c=L_j\); divide by Lemma Z.
Same as `LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` Proposition 1 with
\(\theta=x_j\gamma\), \(\tau=\sigma_j\). \(\square\)

**Status.** **Complete** for intercept RE with sufficient box; necessary bound on
\(\{\bar\omega\ge\tau\}\) without a box — open (`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §2.3).

---

## 4. Probit (binomial, probit link)

### 4.1 Geometry: two-sided low-curvature tail

\[
w_{j,i}(\eta)=\frac{n_{j,i}\,\varphi(\eta)^2}{\bigl[\Phi(\eta)\,\Phi(-\eta)\bigr]^2},
\qquad
W_j^{\mathrm{sum}}(\beta_j)=\sum_i w_{j,i}(\eta_{j,i}(\beta_j)).
\]

**Hessian / weight vanish at both tails (Mills ratio).** For one observation with
\(z=\eta\), \(x=1\),

\[
-\frac{\partial^2 \ell}{\partial \beta^2}
=
x^2\,\lambda(z)\bigl[\lambda(z)+z\bigr],
\qquad
\lambda(z)=
\begin{cases}
\phi(z)/\Phi(z) & y=1,\\[4pt]
\phi(-z)/\Phi(-z) & y=0,
\end{cases}
\]

**Left tail (\(z\to-\infty\)).** \(\phi(z)\to 0\) and \(\lambda(z)\sim -z\), so
\(\lambda(z)+z\to 0\) and the Hessian \(\to 0\).

**Right tail (\(z\to+\infty\)).** \(\phi(z)\to 0\), \(\lambda(z)\to 0\), so the product
\(\lambda(z)[\lambda(z)+z]\to 0\).

With **both** successes and failures present, \(\mathcal H_{\mathrm{total}}=
\sum_i \mathcal H_{y_i}\) still vanishes at **both** tails because each component
does. The IRLS weight \(w(\eta)\) above is symmetric, maximal at \(\eta=0\)
(\(w\le (2/\pi)\,n\)), and **\(w(\eta)\to 0\) as \(|\eta|\to\infty\)** — aligned
with observed flattening on both sides.

Hence (intercept, offsets \(0\), as §1.4.3)

\[
E_j(\omega_j)
=
\{\,|\beta_j|>b_j^\star(\omega_j)\,\},
\]

the **full** low-\(W\) region on **both** sides — not a one-sided tail.

*(Earlier notes said “mostly lower tail” because complete separation is often
**observed** on the \(y=1\) / negative-\(\eta\) side, and because \(w\le 2/\pi\) was
misread as “bounded away from zero” at \(+\infty\). For **Prop 1 on
\(E_j(\omega_j)\)**, the correct geometry is two-sided.)*

### 4.2 Proposition P1 (Proposition 1)

With \(b_j^\star(\omega_j)\) from §1.4.3 and \(L_j=L_j^{\mathrm{pr}}\),

\[
\boxed{
g_j(\gamma;\omega_j)
\le
G_j^{\mathrm{pr}}(\gamma;\omega_j)
:=
G_j^{\mathrm{2s}}(\gamma;\omega_j)
}
\]

*Proof.* Same as Proposition L1 (§3.2): split \(E_j\) at \(\pm b_j^\star(\omega_j)\);
Lemma T holds globally; Lemmas G/G′; Lemma Z. \(\square\)

**Status.** **Complete** on \(E_j(\omega_j)\) for intercept RE (margin via numeric
inversion of \(w(\eta^\star)=\kappa_j/\sum n\)). Necessary bound on \(\{\bar\omega\ge\tau\}\)
without a box — open.

---

## 5. Poisson (log link)

### 5.1 Geometry: entire low-curvature tail

\(w_{j,i}(\eta)=e^{\eta}\), \(\ell_{j,i}''=-e^{\eta}\le 0\).
**Rate-unsafe region:** \(W_j^{\mathrm{sum}}\to 0\) as \(\eta_{j,i}\to-\infty\);
\[
E_j(\omega_j)
=
\{\,\beta_j:\omega_j(\beta_j)\ge\omega_j\,\}
=
\{\,\beta_j:W_j^{\mathrm{sum}}(\beta_j)\le \kappa_j(\omega_j)\,\}.
\]
This is the **full** tail where effective Hessian weight is below the certification
floor — not a subset cut off by an unrelated safe-box floor.

**Intercept, offsets \(0\), \(x_j>0\):** \(W_j^{\mathrm{sum}}(\beta_j)=\sum_i e^{x_j\beta_j}\)
is **strictly increasing** in \(\beta_j\), so
\[
E_j(\omega_j)
=
\{\,\beta_j\le \beta_j^\star(\omega_j)\,\},
\qquad
W_j^{\mathrm{sum}}(\beta_j^\star(\omega_j))=\kappa_j(\omega_j).
\]
*(If \(x_j<0\), the tail is \(\beta_j\ge\beta_j^\star(\omega_j)\); flip signs below.)*

**Closed form (intercept, offsets \(0\), \(x_j\neq 0\)).**
\(W_j^{\mathrm{sum}}(\beta_j)=\bigl(\sum_i 1\bigr)\,e^{x_j\beta_j}=n_j e^{x_j\beta_j}\) when
each observation shares the same \(\eta=x_j\beta_j\), giving
\[
\beta_j^\star(\omega_j)
=
\frac{1}{x_j}\log\!\Bigl(\frac{\kappa_j(\omega_j)}{n_j}\Bigr).
\]
With offsets, \(W_j^{\mathrm{sum}}(\beta_j)=\sum_i e^{\mathrm{offset}_{j,i}+x_j\beta_j}
=e^{x_j\beta_j}\sum_i e^{\mathrm{offset}_{j,i}}\), so
\(\beta_j^\star(\omega_j)=x_j^{-1}\bigl(\log\kappa_j(\omega_j)-\log\sum_i e^{\mathrm{offset}_{j,i}}\bigr)\).

### 5.2 Score bound on all of \(E_j(\omega_j)\)

Poisson score \(\partial\ell_{j,i}/\partial\eta=y_{j,i}-e^{\eta_{j,i}}\).
On \(E_j(\omega_j)\), each \(\eta_{j,i}(\beta_j)\le \eta_{j,i}^\star(\omega_j)\), so
\(|y_{j,i}-e^{\eta_{j,i}}|\le \max(y_{j,i},e^{\eta_{j,i}^\star(\omega_j)})\).
Hence Lemma T applies on **every** \(\beta_j\in E_j(\omega_j)\) with
\(L_j=L_j^{\mathrm{po}}(\omega_j)\) from §1.4.2.

*(Remark: the score is unbounded on all of \(\mathbb R\) as \(\eta\to+\infty\), but
that region is **outside** \(E_j(\omega_j)\) — large \(\eta\) gives large \(W\),
low shrinkage, not in the tail.)*

### 5.3 Proposition Po1 (Proposition 1)

Let \(\beta_j^\star(\omega_j)\) and \(L_j^{\mathrm{po}}(\omega_j)\) be as in §1.4.2.
For the one-sided tail \(E_j(\omega_j)\subseteq\{\beta_j\le\beta_j^\star(\omega_j)\}\),

\[
\boxed{
g_j(\gamma;\omega_j)
\le
G_j^{\mathrm{po}}(\gamma;\omega_j)
:=
G_j^{\mathrm{low}}(\gamma;\omega_j)
}
\]

with \(t_j=\beta_j^\star(\omega_j)\) and \(L_j=L_j^{\mathrm{po}}(\omega_j)\) in
\(G_j^{\mathrm{low}}\) (§1.3).

*Proof.* Proposition 1 template: on \(E_j\), Lemma T with \(b_0=x_j\gamma\);
majorize \(\int_{E_j}e^{\ell_j}\phi\,d\beta_j\) by the lower-tail tilted integral
(Lemma G′ / one-sided \(G_j^{\mathrm{low}}\)); divide by Lemma Z. \(\square\)

**Safe box (separate).** To certify \(B(\boldsymbol{\omega})\subseteq A\), one may
still use a high-\(\eta\) floor \(\beta_j\ge\log s_j\) as in
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §4.3; that construction certifies the
**working safe region**, not the domain of Proposition Po1.

**Status.** **Complete** for the full low-\(W\) tail \(E_j(\omega_j)\) when
\(\beta_j^\star(\omega_j)\) is computed from \(\kappa_j(\omega_j)\); varying
offsets require solving \(W_j^{\mathrm{sum}}(\beta_j^\star)=\kappa_j\) numerically.

---

## 6. Complementary log-log (binomial, cloglog link)

### 6.1 Geometry: one-sided observed flattening (mixed data)

**Fisher weight** \(w(\eta)\to 0\) at both \(\eta\)-tails (§1.4.5), but with both
successes and failures the **observed** Hessian vanishes only on the **lower** tail
(\(\eta\to-\infty\)) and **diverges** on the upper tail (\(\eta\to+\infty\),
\(\partial^2\ell/\partial\eta^2\sim -n_0 e^{\eta}\)). Prop 1 targets the flat,
score-bounded regime, so

\[
E_j(\omega_j)
=
\{\,\beta_j:\omega_j(\beta_j)\ge\omega_j\,\}
\quad\text{on the lower branch only}
=
\{\,\beta_j\le \beta_j^\star(\omega_j)\,\}
\]

(intercept, \(x_j>0\); flip if \(x_j<0\)), **not** a two-sided outer union from
\(\{W_j^{\mathrm{sum}}\le\kappa\}\).

*(Contrast **probit** §4: Mills ratio gives observed \(\mathcal H\to 0\) at **both**
tails even with mixed data.)*

### 6.2 Score on \(E_j(\omega_j)\)

Use \(L_j^{\mathrm{cl}}(\omega_j)\) from §1.4.5 on the lower tail margin
\(\eta_{j,i}^\star(\omega_j)\). Lemma T applies on \(E_j\) as in Poisson §5.

### 6.3 Proposition C1 (Proposition 1)

\[
\boxed{
g_j(\gamma;\omega_j)
\le
G_j^{\mathrm{cl}}(\gamma;\omega_j)
:=
G_j^{\mathrm{low}}(\gamma;\omega_j)
}
\]

with \(t_j=\beta_j^\star(\omega_j)\), \(L_j=L_j^{\mathrm{cl}}(\omega_j)\).

*Proof.* Same as Proposition Po1 (§5), lower tail only. \(\square\)

**Status.** **Complete** on the lower flat tail \(E_j(\omega_j)\) for mixed binomial
groups when \(\beta_j^\star(\omega_j)\) is computed from \(\kappa_j(\omega_j)\) on the
\(\eta\to-\infty\) branch.

---

## 7. Gamma regression (log link, fixed dispersion)

### 7.1 Specialization

Constant \(\phi\): \(w_{j,i}=1/\phi\), \(W_j^{\mathrm{sum}}=n_j/\phi\),
**\(\beta\)-independent** (same as Gaussian §2).

### 7.2 Proposition G1 (Proposition 1 analog)

On certified safe \(\omega_j\) with \(\omega_j^{\mathrm{Ga}}(\phi)<\omega_j\),

\[
\boxed{
g_j(\gamma;\omega_j)=0,
\qquad
G_j(\gamma;\omega_j)\equiv 0.
}
\]

Identical to Proposition N1; \(\lambda^{\star,(0)}(\boldsymbol{\omega})\) from
`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §6.

*Reference:* `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §6 Proposition G1.

---

## 8. Summary

| Likelihood | \(L_j(\omega_j)\) | Tail \(E_j(\omega_j)\) | \(G_j\) | Prop 1 status |
|---|---|---|---|---|
| Gaussian | — | \(\varnothing\) if certified | \(0\) | **Exact / trivial** |
| Gamma (log) | — | \(\varnothing\) if certified | \(0\) | **Exact / trivial** |
| Logit | global \(L_j\) | two-sided \(|\beta_j|>b_j^\star(\omega_j)\) | \(G_j^{\mathrm{2s}}\) | **Complete** |
| Probit | global \(L_j^{\mathrm{pr}}\) | two-sided \(|\beta_j|>b_j^\star(\omega_j)\) | \(G_j^{\mathrm{2s}}\) | **Complete** on \(E_j\) |
| Poisson | \(L_j^{\mathrm{po}}(\omega_j)\) on \(E_j\) | full low-\(W\) tail (one-sided) | \(G_j^{\mathrm{low}}\) | **Complete** on \(E_j\) |
| Cloglog | \(L_j^{\mathrm{cl}}(\omega_j)\) on \(E_j\) | lower tail only (mixed data) | \(G_j^{\mathrm{low}}\) | **Complete** on \(E_j\) |

**Next steps (likelihood-agnostic).** (i) implement \(\beta_j^\star(\omega_j)\) from
\(\kappa_j(\omega_j)\) (closed form for intercept; numeric if offsets vary);
(ii) verify monotonicity of \(G_j\) in \(\gamma\) for gate inversion;
(iii) plug into §1.5 for \(\varepsilon_2(\boldsymbol{\omega})\) and \(B(\boldsymbol{\omega})\).

**Sampler programs.** Two ways to connect these bounds to MCMC convergence are
summarized in `LOGIT_STATIC_TAIL_CERTIFICATION.md` §7: **Approach A** uses Rosenthal
**Proposition 2** (minorization on the safe space); **Approach B** uses static tail
bounds and Claims 1–3.

---

## References (in repo)

- `inst/RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0–§6
- `inst/LOGIT_STATIC_TAIL_CERTIFICATION.md` — §7 two program approaches
- `inst/LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md`
- `inst/CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §IV, §6
- `inst/notation.md`
