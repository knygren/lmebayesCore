# Logit likelihood: marginal on group coefficients after integrating population parameters

Draft companion to `inst/notation.md`, `inst/CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md`
§2 (logit), and `inst/RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0–§2.

**Goal.** Write the **unnormalized** target density on the group-level
coefficients \(\beta = (\beta_1,\ldots,\beta_J)\) after integrating out the
population (hyper-mean) parameters \(\gamma\), for **binomial logit** groups.
The logit part does not integrate in closed form; the **population integral**
is exact Gaussian and yields a **Schur-complement** quadratic on \(\beta\).
**Static tail certification** (conditional \(\pi(\beta \mid \gamma,y)\), per-group
\(G_j(\gamma)\), profile \(h(\gamma)\), two-level \(\gamma\) partition — in
`LOGIT_STATIC_TAIL_CERTIFICATION.md`. **§11 below:** Gibbs escape program (separate topic).

**Notation.** Follow `notation.md`: \(\mathcal{W}_j\) is the level-2 design
(group \(j\), dimension \(P\)), \(\gamma \in \mathbb{R}^q\) is the population
vector, \(\Psi\) is the stage-2 covariance, and \(P_b := \Psi^{-1}\) is the RE
prior precision (rate-note / Chapter-C03 convention). In
`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` the same level-2 object is
written \(H_j\); we use \(\mathcal{W}_j \equiv H_j\) when citing \(P\)-blocks.
Do **not** confuse \(\mathcal{W}_j\) with the **data** Gram matrix
\(\sum_i w_{j,i} z_{j,i} z_{j,i}^\top\) from Chapter-C03 §2.2 (there denoted
\(C_j(\beta_j)\) or \(\mathcal{W}_j^{\mathrm{data}}\) below).

---

## 1. Hierarchical model (logit stage 1)

For groups \(j = 1,\ldots,J\), observations \(i = 1,\ldots,n_j\):

**Likelihood (binomial, canonical logit link).**

\[
y_{j,i} \mid \beta_j
\;\sim\;
\mathrm{Binomial}\bigl(n_{j,i},\, p_{j,i}(\beta_j)\bigr),
\qquad
p_{j,i}(\beta_j) = \mathrm{expit}(\eta_{j,i}),
\]

\[
\eta_{j,i}
=
\mathrm{offset}_{j,i} + d_{j,i}^\top \beta_j,
\qquad
d_{j,i}^\top = \text{row $i$ of $D_j$ in \texttt{notation.md}.}
\]

Write the group log-likelihood (logit, \(\phi = 1\)):

\[
\boxed{
\ell_j(\beta_j)
=
\sum_{i=1}^{n_j}
\Bigl[
y_{j,i}\,\eta_{j,i}
-
n_{j,i}\,\log\bigl(1 + e^{\eta_{j,i}}\bigr)
\Bigr].
}
\]

**Stage 2 (centered RE prior).**

\[
\beta_j \mid \gamma,\, \Psi
\;\sim\;
N\bigl(\mathcal{W}_j \gamma,\, \Psi\bigr),
\qquad
\Psi \succ 0 \ \text{(fixed for this note)}.
\]

**Stage 3 (population / hyper-mean prior).**

\[
\gamma
\;\sim\;
N\bigl(\mu_0,\, \Lambda_\gamma^{-1}\bigr),
\qquad
\Lambda_\gamma
=
\mathrm{blockdiag}\bigl(V_1^{-1},\ldots,V_{p_{\mathrm{re}}}^{-1}\bigr)
\ \text{from \texttt{notation.md}.}
\]

Stack \(\beta = (\beta_1^\top,\ldots,\beta_J^\top)^\top \in \mathbb{R}^{JP}\),
\(\mathcal{W} = \mathrm{blockdiag}(\mathcal{W}_1,\ldots,\mathcal{W}_J)\).

---

## 2. Joint unnormalized posterior on \((\beta, \gamma)\)

\[
\boxed{
\widetilde{\pi}(\beta, \gamma \mid y)
=
\exp\Bigl(\sum_{j=1}^{J} \ell_j(\beta_j)\Bigr)
\;
\exp\Bigl(
-\tfrac{1}{2}(\beta - \mathcal{W}\gamma)^\top \Psi^{-1}(\beta - \mathcal{W}\gamma)
-\tfrac{1}{2}(\gamma - \mu_0)^\top \Lambda_\gamma (\gamma - \mu_0)
\Bigr).
}
\]

Equivalently, with \(P_b = \Psi^{-1}\),

\[
\widetilde{\pi}(\beta, \gamma \mid y)
\propto
\exp\Bigl(\sum_j \ell_j(\beta_j)\Bigr)
\;
\exp\Bigl(
-\tfrac{1}{2}(\beta - \mathcal{W}\gamma)^\top P_b^{\mathrm{stack}}(\beta - \mathcal{W}\gamma)
-\tfrac{1}{2}(\gamma - \mu_0)^\top \Lambda_\gamma (\gamma - \mu_0)
\Bigr),
\]

where \(P_b^{\mathrm{stack}} = \mathrm{blockdiag}(P_b,\ldots,P_b)\) when \(\Psi\)
is shared across groups (general block structure if \(\Psi\) varies by group).

This is **not** a Gaussian joint density because of \(\ell_j\); only the
**prior** factor is quadratic in \((\beta, \gamma)\).

---

## 3. Prior precision in \((\gamma, \beta)\) block form

Order blocks as **population first**, **RE second** (Chapter-C03 / block updated
second for \(\gamma\)). The **prior** log-density alone is
\(-\tfrac{1}{2}(\gamma,\beta)^\top \mathfrak{P} (\gamma,\beta) + \text{linear}\)
with

\[
\boxed{
\mathfrak{P}_{11}
=
\Lambda_\gamma
+
\sum_{j=1}^{J} \mathcal{W}_j^\top P_b \mathcal{W}_j,
\qquad
\mathfrak{P}_{22}
=
P_b^{\mathrm{stack}},
}
\]

\[
\boxed{
\mathfrak{P}_{12}
=
\Bigl[\,-\mathcal{W}_1^\top P_b \;\Big|\; \cdots \;\Big|\; -\mathcal{W}_J^\top P_b \,\Bigr],
\qquad
\mathfrak{P}_{21} = \mathfrak{P}_{12}^\top.
}
\]

These are the **prior-only** analogues of \(P_{11}, P_{22}, P_{12}, P_{21}\) in
Chapter-C03 §I.3 **before** adding data curvature
\(Z_j^\top W_j^{\mathrm{data}} Z_j\) to the \(\beta\)-block.

---

## 4. Integrate \(\gamma\): exact Gaussian marginal prior on \(\beta\)

Fix \(\beta\). The \(\gamma\)-integral involves only the prior quadratic part:

\[
\widetilde{\pi}_{\mathrm{prior}}(\beta)
:=
\int_{\mathbb{R}^q}
\exp\Bigl(
-\tfrac{1}{2}(\beta - \mathcal{W}\gamma)^\top P_b^{\mathrm{stack}}(\beta - \mathcal{W}\gamma)
-\tfrac{1}{2}(\gamma - \mu_0)^\top \Lambda_\gamma (\gamma - \mu_0)
\Bigr)\,d\gamma.
\]

This is a standard Lindley–Smith marginalization. Define the **Schur precision**
on \(\beta\):

\[
\boxed{
\Lambda_\beta
:=
\mathfrak{P}_{22} - \mathfrak{P}_{21}\,\mathfrak{P}_{11}^{-1}\,\mathfrak{P}_{12}
=
P_b^{\mathrm{stack}}
-
P_b^{\mathrm{stack}} \mathcal{W}\,
\mathfrak{P}_{11}^{-1}\,
\mathcal{W}^\top P_b^{\mathrm{stack}}.
}
\]

The **prior mean** of \(\beta\) before data is

\[
\boxed{
\mu_\beta := \mathcal{W}\,\mu_0.
}
\]

Then

\[
\boxed{
\widetilde{\pi}_{\mathrm{prior}}(\beta)
=
|\mathfrak{P}_{11}|^{-1/2}\,(2\pi)^{q/2}\;
\exp\Bigl(
-\tfrac{1}{2}(\beta - \mu_\beta)^\top \Lambda_\beta\,(\beta - \mu_\beta)
\Bigr),
}
\]

up to a \(\gamma\)-only constant. Equivalently \(\Lambda_\beta^{-1} = \Psi^{\mathrm{stack}}
+ \mathcal{W}\,\Lambda_\gamma^{-1}\,\mathcal{W}^\top\) (covariance form).

**Important:** integrating \(\gamma\) **couples groups** in the prior on \(\beta\):
\(\Lambda_\beta\) is generally **not** block-diagonal across \(j\), even when
\(\Psi\) is block-diagonal and \(\ell_j\) factorizes over \(j\).

---

## 5. Unnormalized marginal on group parameters (functional form)

Combining §2 and §4,

\[
\boxed{
\widetilde{\pi}(\beta \mid y)
=
\exp\Bigl(\sum_{j=1}^{J} \ell_j(\beta_j)\Bigr)
\;
\exp\Bigl(
-\tfrac{1}{2}(\beta - \mu_\beta)^\top \Lambda_\beta\,(\beta - \mu_\beta)
\Bigr).
}
\]

Absorbing constants into proportionality,

\[
\boxed{
\widetilde{\pi}(\beta \mid y)
\;\propto\;
\Biggl[\prod_{j=1}^{J}
\exp\bigl(\ell_j(\beta_j)\bigr)\Biggr]
\;
\exp\Bigl(
-\tfrac{1}{2}(\beta - \mathcal{W}\mu_0)^\top \Lambda_\beta\,(\beta - \mathcal{W}\mu_0)
\Bigr).
}
\]

**Interpretation.**

| Factor | Role |
|---|---|
| \(\exp(\sum_j \ell_j(\beta_j))\) | **Logit likelihood**; non-Gaussian; no closed \(\beta\)-integrals |
| Gaussian in \(\beta\) with \(\Lambda_\beta\) | **Integrated population prior**; exact from §4 |
| \(\mu_\beta = \mathcal{W}\mu_0\) | Center from population prior mean, not from data |

This is the **unnormalized marginal target** on \(\beta\) after integrating
\(\gamma\) out of the **prior × likelihood** integrand. The normalizing constant
\(Z(y) = \int \widetilde{\pi}(\beta \mid y)\,d\beta\) is **not** available in
closed form for logit.

**Contrast with Gibbs conditioning.** The sampler uses
\(\pi(\beta_j \mid \gamma, y_j) \propto \exp(\ell_j(\beta_j))\,
\phi(\beta_j; \mathcal{W}_j\gamma, \Psi)\) for fixed \(\gamma\)
(`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0.1). That is **not** the same object as
\(\widetilde{\pi}(\beta \mid y)\): conditioning keeps \(\gamma\); here \(\gamma\)
is integrated out and cross-group coupling enters through \(\Lambda_\beta\).

### 5.1 Partial factorization by group ( \(\beta_j\) not independent)

There is **no** rewrite of \(\log \widetilde{\pi}(\beta \mid y)\) as
\(\sum_{j=1}^{J} g_j(\beta_j)\) with each \(g_j\) depending on **only**
\(\beta_j\): integrating \(\gamma\) couples groups through \(\Lambda_\beta\) (§4).
What follows is the best **partial** grouping—per-group logit/RE pieces plus one
global coupling term in \(s(\beta)\)—not a product of independent group factors.

Write \(u_j := \beta_j - \mathcal{W}_j \mu_0\) and stack \(u = (\beta - \mu_\beta)\).
The **joint** log-target on \((\beta, \gamma)\) already separates over groups in
the logit and stage-2 factors:

\[
\boxed{
\log \widetilde{\pi}(\beta, \gamma \mid y)
=
-\tfrac{1}{2}(\gamma - \mu_0)^\top \Lambda_\gamma (\gamma - \mu_0)
+
\sum_{j=1}^{J}
\Bigl[
\ell_j(\beta_j)
-
\tfrac{1}{2}(\beta_j - \mathcal{W}_j \gamma)^\top P_b (\beta_j - \mathcal{W}_j \gamma)
\Bigr].
}
\]

(Use \(P_b = \Psi^{-1}\); replace by \(P_{b,j}\) if \(\Psi\) varies by group.)

**After integrating \(\gamma\).** Complete the square in \(\gamma\); define the
population sufficient statistic (a \(q\)-vector)

\[
\boxed{
s(\beta)
:=
\Lambda_\gamma \mu_0
+
\sum_{j=1}^{J} \mathcal{W}_j^\top P_b\, u_j
=
\Lambda_\gamma \mu_0
+
\sum_{j=1}^{J} \mathcal{W}_j^\top P_b \bigl(\beta_j - \mathcal{W}_j \mu_0\bigr).
}
\]

Then the unnormalized log-marginal splits into a **sum over \(j\)** of
\(\beta_j\)-only terms plus a **single** quadratic in \(s(\beta)\) that mixes
all groups:

\[
\boxed{
\log \widetilde{\pi}(\beta \mid y)
=
C
+
\sum_{j=1}^{J}
\Bigl[
\ell_j(\beta_j)
-
\tfrac{1}{2}\, u_j^\top P_b\, u_j
\Bigr]
+
\tfrac{1}{2}\,
s(\beta)^\top \mathfrak{P}_{11}^{-1}\, s(\beta),
}
\]

where \(C\) is a constant in \(\beta\) (depends on \((\mu_0, \Lambda_\gamma,
\mathfrak{P}_{11})\) only) and \(\mathfrak{P}_{11} = \Lambda_\gamma +
\sum_j \mathcal{W}_j^\top P_b \mathcal{W}_j\) as in §3. Equivalently,
\(\tfrac{1}{2} s^\top \mathfrak{P}_{11}^{-1} s
= \tfrac{1}{2}\|\mathfrak{P}_{11}^{-1/2} s(\beta)\|^2\).

*Proof sketch.* Expand the prior exponent as a quadratic in \((\gamma, \beta)\),
integrate the Gaussian in \(\gamma\); the \(\gamma\)-linear term produces
\(\mathfrak{P}_{11}^{-1} s(\beta)\) and the \(\gamma\)-quadratic produces
\(-\tfrac{1}{2} s^\top \mathfrak{P}_{11}^{-1} s\); the remaining \(\beta\)-only
piece is \(-\tfrac{1}{2}\sum_j u_j^\top P_b u_j\). This is the same algebra as
§4, only grouped by \(j\). \(\square\)

**Expanded square (pairwise form).** Since
\(s(\beta) = \Lambda_\gamma \mu_0 + \sum_j a_j(\beta_j)\) with
\(a_j(\beta_j) := \mathcal{W}_j^\top P_b u_j\),

\[
\tfrac{1}{2} s^\top \mathfrak{P}_{11}^{-1} s
=
\tfrac{1}{2}(\Lambda_\gamma \mu_0)^\top \mathfrak{P}_{11}^{-1} (\Lambda_\gamma \mu_0)
+
\sum_{j=1}^{J} (\Lambda_\gamma \mu_0)^\top \mathfrak{P}_{11}^{-1} a_j(\beta_j)
+
\tfrac{1}{2}
\sum_{j=1}^{J}\sum_{k=1}^{J}
a_j(\beta_j)^\top \mathfrak{P}_{11}^{-1} a_k(\beta_k).
\]

The double sum is the **cross-group coupling** from shared \(\gamma\). It matches
the block Schur form

\[
-\tfrac{1}{2} u^\top \Lambda_\beta u
=
\sum_{j=1}^{J}\Bigl(-\tfrac{1}{2} u_j^\top \Lambda_{\beta,jj} u_j\Bigr)
-
\sum_{1 \le j < k \le J} u_j^\top \Lambda_{\beta,jk} u_k,
\]

with \(\Lambda_{\beta,jj} = P_b - P_b \mathcal{W}_j \mathfrak{P}_{11}^{-1}
\mathcal{W}_j^\top P_b\) and
\(\Lambda_{\beta,jk} = - P_b \mathcal{W}_j \mathfrak{P}_{11}^{-1}
\mathcal{W}_k^\top P_b\) for \(j \neq k\).

**Exponentiated form (product over \(j\) × one coupling factor).**

\[
\widetilde{\pi}(\beta \mid y)
\;\propto\;
\Biggl[\prod_{j=1}^{J}
\exp\Bigl(
\ell_j(\beta_j)
-
\tfrac{1}{2} u_j^\top P_b u_j
\Bigr)\Biggr]
\;
\exp\Bigl(
\tfrac{1}{2} s(\beta)^\top \mathfrak{P}_{11}^{-1} s(\beta)
\Bigr).
\]

The logit and per-group RE pieces factor by \(j\); the **only** term that is not
a function of a **single** \(\beta_j\) alone is
\(\exp\bigl(\tfrac{1}{2} s(\beta)^\top \mathfrak{P}_{11}^{-1} s(\beta)\bigr)\),
which depends on \(\beta\) through the **sum** \(\sum_j \mathcal{W}_j^\top P_b
u_j\).

**Scalar intercept case (\(\mathcal{W}_j = x_j\), \(q=1\)).** Then
\(s(\beta) = \lambda_\gamma \mu_0 + \lambda_b \sum_j x_j u_j\) and

\[
\log \widetilde{\pi}(\beta \mid y)
=
C
+
\sum_{j=1}^{J}
\Bigl[
\ell_j(\beta_j)
-
\tfrac{\lambda_b}{2}(\beta_j - x_j \mu_0)^2
\Bigr]
+
\tfrac{1}{2\mathfrak{P}_{11}}
\Bigl(
\lambda_\gamma \mu_0
+
\lambda_b \sum_{j=1}^{J} x_j (\beta_j - x_j \mu_0)
\Bigr)^{\!2},
\qquad
\mathfrak{P}_{11} = \lambda_\gamma + \lambda_b \textstyle\sum_j x_j^2.
\]

When \(x_j \equiv 1\), the coupling is the single square
\(\tfrac{1}{2\mathfrak{P}_{11}}(\lambda_\gamma \mu_0
+ \lambda_b \sum_j u_j)^2\).

### 5.2 Static tail certification (moved)

Certified bounds on **\(g(\gamma)=\pi(A^c \mid \gamma,y)\)** — per-group Prop L1
tails \(G_j(\gamma)\), profile \(h(\gamma)=1-\prod_j(1-G_j(\gamma))\), and
two-level \((\varepsilon_1,\varepsilon_2)\) partition with per-group thresholds —
are in **`LOGIT_STATIC_TAIL_CERTIFICATION.md`** (conditional on \(\gamma\);
optional average over \(\pi(\gamma \mid y)\) only in §4.2 there).

---

## 6. Block Schur form (equivalent to §5.1)

Write \(u_j := \beta_j - \mathcal{W}_j \mu_0\). Stack \(u = \beta - \mu_\beta\).
The §5.1 sum-over-\(j\) identity is equivalent to

\[
\widetilde{\pi}(\beta \mid y)
\;\propto\;
\exp\Bigl(\sum_{j=1}^{J} \ell_j(\beta_j)\Bigr)
\;
\exp\Bigl(-\tfrac{1}{2}\, u^\top \Lambda_\beta\, u\Bigr).
\]

Per group, the logit factor depends only on \(\beta_j\). The integrated-prior
coupling is the off-diagonal part of the same quadratic:

\[
-\tfrac{1}{2} u^\top \Lambda_\beta u
=
\sum_{j=1}^{J}\Bigl(-\tfrac{1}{2} u_j^\top \Lambda_{\beta,jj} u_j\Bigr)
-
\sum_{1 \le j < k \le J} u_j^\top \Lambda_{\beta,jk} u_k
\qquad\text{(see §5.1)}.
\]

---

## 7. Scalar intercept RE specialization (\(P = 1\), \(\mathcal{W}_j = x_j\))

One random intercept per group, scalar \(\beta_j \in \mathbb{R}\),
\(\mathcal{W}_j = x_j \in \mathbb{R}^{1\times q}\) (often \(x_j \equiv 1\)),
\(\Psi = \lambda_b^{-1}\), population precision \(\lambda_\gamma\) on \(\gamma\)
(scalar \(q=1\) case).

**Logit.**

\[
\ell_j(\beta_j)
=
\sum_{i\in j}
\bigl[
y_{j,i}(\mathrm{offset}_{j,i} + \beta_j)
-
n_{j,i}\log(1 + e^{\mathrm{offset}_{j,i} + \beta_j})
\bigr].
\]

**Prior blocks.**

\[
\mathfrak{P}_{11}
=
\lambda_\gamma
+
\lambda_b \sum_{j=1}^{J} x_j^2,
\qquad
\mathfrak{P}_{22}
=
\lambda_b I_J,
\qquad
\mathfrak{P}_{12,j}
=
-\lambda_b x_j.
\]

**Schur precision on \(\beta = (\beta_1,\ldots,\beta_J)^\top\).**

\[
\Lambda_\beta
=
\lambda_b I_J
-
\lambda_b^2
\begin{pmatrix} x_1 \\ \vdots \\ x_J \end{pmatrix}
\bigl[\lambda_\gamma + \lambda_b \textstyle\sum_j x_j^2\bigr]^{-1}
\begin{pmatrix} x_1 & \cdots & x_J \end{pmatrix}.
\]

**Center.** \(\mu_\beta = x_j \mu_0\) per component when \(x_j \equiv 1\):
\(\mu_\beta = \mu_0 \mathbf{1}_J\).

**Unnormalized marginal.**

\[
\boxed{
\widetilde{\pi}(\beta \mid y)
\;\propto\;
\exp\Bigl(\sum_{j=1}^{J} \ell_j(\beta_j)\Bigr)
\;
\exp\Bigl(
-\tfrac{1}{2}(\beta - \mu_0 \mathbf{1}_J)^\top \Lambda_\beta\,(\beta - \mu_0 \mathbf{1}_J)
\Bigr).
}
\]

When \(x_j \equiv 1\), \(\mathfrak{P}_{11} = \lambda_\gamma + J\lambda_b\) and
\(\Lambda_\beta = \lambda_b I_J - \lambda_b^2 (J\lambda_b + \lambda_\gamma)^{-1}
\mathbf{1}\mathbf{1}^\top\).

---

## 8. Single-group marginal on \(\beta_j\) (still not closed)

Marginalizing **all other groups** as well would require integrating the coupled
Gaussian factor; for \(J>1\) there is no closed form for
\(\widetilde{\pi}(\beta_j \mid y)\) because of \(\Lambda_\beta\).

For **one group** (\(J=1\)), the integrated population prior is scalar Gaussian
and

\[
\widetilde{\pi}(\beta_1 \mid y)
\;\propto\;
\exp\bigl(\ell_1(\beta_1)\bigr)\,
\exp\Bigl(
-\tfrac{1}{2}\,\lambda_{\beta,1}\,(\beta_1 - x_1 \mu_0)^2
\Bigr),
\qquad
\lambda_{\beta,1}
=
\lambda_b
-
\frac{\lambda_b^2 x_1^2}{\lambda_\gamma + \lambda_b x_1^2}.
\]

Still **no closed-form normalizer** in \(\beta_1\) because \(\ell_1\) is logit.

---

## 9. Relation to Chapter-C03 \(P\)-matrix at the mode

At a working point \(\beta\), the **posterior** joint precision adds data
curvature to the \(\beta\)-block:

\[
P_{22,j} = \underbrace{Z_j^\top W_j^{\mathrm{data}}(\beta_j)\, Z_j}_{\text{logit Fisher}}
+ P_b,
\qquad
W_{j,i}^{\mathrm{data}} = n_{j,i}\, p_{j,i}(1-p_{j,i}).
\]

That **Laplace / IRLS** \(P\)-matrix is what enters the rate matrix
\(A = P_{11}^{-1/2} P_{12} P_{22}^{-1} P_{21} P_{11}^{-1/2}\). It is **not**
the same as the prior Schur \(\Lambda_\beta\) in §4–§5: §5 integrates \(\gamma\)
from the **Gaussian prior** only; the logit enters **multiplicatively** outside
that integral. A Laplace approximation to \(\widetilde{\pi}(\beta \mid y)\) would
replace \(\sum_j \ell_j\) by a quadratic with precision \(P_{22,j} - P_b\) at
the mode — a different object from exact marginalization.

---

## 10. Static tail certification (moved)

See **LOGIT_STATIC_TAIL_CERTIFICATION.md**. This note (§1–§9) stops at the
unnormalized marginal \(\widetilde{\pi}(\beta \mid y)\) and its structure; unsafe
mass certification is developed separately.

---

## 11. Target program: three escape claims (no false drift)

The Gibbs chain can **stick in the tail** (\(A^c\)); Foster–Lyapunov drift (rate
note §8.3) is **false** there. The viable theory splits: **equilibrium** tail mass is certified statically in
`LOGIT_STATIC_TAIL_CERTIFICATION.md`; **sampler** comparison uses the three claims
below (no false drift on \(A^c\)).

**Notation (Gibbs on \((\beta,\gamma)\)).** Let \(\pi\) be the **target**
posterior. One sweep: \(\beta_j \mid \gamma,y_j\), then \(\gamma \mid \beta\).
Write \(g(\gamma) := \pi(A^c \mid \gamma,y)\) and, for the \(\gamma\)-marginal
chain \(\{\gamma_n\}\) started at \(\gamma_0=x\),

\[
P^n(x,A^c)
\;:=\;
P\bigl(\beta_n \in A^c \,\big|\, \gamma_0 = x\bigr)
\;=\;
E\bigl[g(\gamma_{n-1}) \mid \gamma_0 = x\bigr]
\]

(stochastic \(\beta_n \mid \gamma_{n-1}\)). Stationary escape
\(\pi(A^c) := P(\beta \in A^c)\) under \(\pi\) (equivalently
\(E_\pi[g(\gamma)]\) when the chain targets \(\pi\)).

**Tradeoff (design of \(A\)).** Enlarging the **safe** set \(A\) (e.g.\
raising the certification threshold \(\tau\) so more \((\beta,\gamma)\) working
points qualify as “safe”, or allowing larger average shrinkage \(\bar\omega\)
inside \(A\)) shrinks \(A^c\) and should make \(\pi(A^c)\) **smaller**, at the
cost of **slower** within-\(A\) Gibbs contraction (\(\lambda^\star \uparrow\)
on \(A\)). Claim (1) makes that tradeoff **quantitative**.

---

### Claim 1 — Tight bound; \(\pi(A^c)\) arbitrarily small by enlarging \(A\)

**Target.**

\[
\boxed{
\pi(A^c) \;\le\; B(\text{geometry of } A;\ \text{design},\ \text{priors})
\;=\; \varepsilon(\tau,\ \lambda_b,\ \lambda_\gamma,\ \{n_j\},\ \ldots),
}
\]

with \(B\) **explicit** and **tight** enough to use in TV, and such that for
GLMs (logit first) one can make \(\varepsilon \downarrow 0\) by enlarging \(A\)
(increasing allowed shrinkage / \(\tau\)) while accepting slower mixing **on**
\(A\) (\(\sup_{A} \lambda^\star \uparrow\)).

**Status.** **Partial.** Per-\(\gamma\) certificates and the two-level partition are in
`LOGIT_STATIC_TAIL_CERTIFICATION.md` (\(g(\gamma)\), \(h(\gamma)\),
\(\varepsilon_1^{\mathrm{eff}}\)). **Open:** certifying \(\varepsilon_2\); optional
\(\beta\)-marginal average if needed for TV.

**Role.** Certify tail mass at equilibrium is negligible **without** claiming
drift exits the tail.

---

### Claim 2 — \(n\)-step escape never exceeds stationary: \(P^n(x,A^c) \le \pi(A^c)\)

**Target.**

\[
\boxed{
P^n(x,A^c) \;\le\; \pi(A^c)
\qquad \forall n \ge 1,\ \forall \text{ admissible starts } x.
}
\]

**Interpretation.** With initialization in (or near) the safe region, the chain
does **not accumulate** more one-step unsafe mass than the stationary tail
proportion. The long-run fraction in \(A^c\) is capped by \(\pi(A^c)\) from
**above** along the path — you cannot “pile up” more than equilibrium tail
density without already being at or above \(\pi(A^c)\).

**Status.** **Open (must be proved; not from drift).** This is **not** a
generic Markov-chain fact; it requires a **model-specific** argument (likely:
structure of \(g(\gamma)\), initialization \(x \in A\) or small \(g(x)\), and
tail stickiness — once unsafe, hard to leave, but also **entry** rate to \(A^c\)
controlled by stationary \(\pi(A^c)\)). **Cannot** use §8.3. Candidate tools:
coupling from a stationary copy, comparison of \(g(\gamma_n)\) under
\(\gamma\)-marginal dynamics without Lyapunov drift on \(A^c\).

**Role.** Ensures finite-\(n\) sampling does not **overshoot** tail mass relative
to the certified equilibrium bound in Claim 1.

---

### Claim 3 — Gap closes: \(\big|\pi(A^c) - P^n(x,A^c)\big|\) small (no drift)

**Target.**

\[
\boxed{
\big|\pi(A^c) - P^n(x,A^c)\big|
\;\le\;
H\bigl(n,\ x;\ \pi(A^c),\ P^n(x,A^c)\bigr)
\;\xrightarrow[n\to\infty]{}\; 0,
}
\]

with \(H\) explicit **without** Foster–Lyapunov drift on \(V(\gamma)\) or
\(d > 2C/(1-\lambda)\). Plausible routes: (i) coupling + renewal on **canonical**
\(C_d\) (minorization notes, \(\delta\)-levels at \(\mu^\star\)); (ii) triangle
inequality from Claims 1–2:
\(|\pi(A^c)-P^n| \le \pi(A^c) + P^n \le 2\pi(A^c)\) if Claim 2 holds; (iii)
Cauchy–Schwarz gap from `LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §6 **once**
\(P(T>n-1)\) is bounded by **valid** minorization, not §8.3.

**Status.** **Open.** Rate note §10.2 algebra is a template; its input is
**invalid** (false drift). **Role.** Prevents **understating** tail mass for
long \(n\): together with Claim 2, \(P^n(x,A^c)\) is sandwiched near
\(\pi(A^c)\).

---

### Why all three

| Claim | Prevents |
|---|---|
| (1) Small \(\pi(A^c)\) | Large **equilibrium** tail mass (design certificate) |
| (2) \(P^n \le \pi(A^c)\) | **Accumulation** above equilibrium tail along the path |
| (3) Small gap | **Underestimating** tail mass at finite \(n\) |

Without (1), TV escape terms are unchecked. Without (2), transient paths could
exceed \(\pi(A^c)\) even when (1) is small. Without (3), \(P^n\) could stay
**below** \(\pi(A^c)\) for many steps and understate coupling/TV. **None** of
the three may use false global drift on the tail.

**Relation to static certification.** Claim (1) uses
`LOGIT_STATIC_TAIL_CERTIFICATION.md` (Props 1–3; **Approach B** in §7). Claims (2)–(3) for
\(P^n(x,A^c)\) remain **open** Markov comparisons without rate-note §8–§11 drift.

---

## 12. Summary

| Quantity | Closed form? |
|---|---|
| \(\int \exp(\text{prior}(\beta,\gamma))\,d\gamma\) | **Yes** — Schur \(\Lambda_\beta\), mean \(\mathcal{W}\mu_0\) |
| \(\widetilde{\pi}(\beta \mid y)\) up to constant | **Functional form yes** — logit × integrated prior Gaussian |
| Normalizer \(Z(y) = \int \widetilde{\pi}(\beta \mid y)\,d\beta\) | **No** (logit) |
| \(\pi(\beta_j \mid y)\) for \(J \ge 2\) | **No** — coupling in \(\Lambda_\beta\) |
| \(\pi(A^c \mid \gamma,y)\) small on \(\Gamma_{\mathrm{safe}}\) | **`LOGIT_STATIC_TAIL_CERTIFICATION.md`** |
| \(E_{\pi(\gamma\mid y)}[g(\gamma)]\) (optional) | same, §4.2 |
| \(P^n(x,A^c) \le \pi(A^c)\) | **Open** Claim 2 §11 (no drift) |
| \(\|\pi(A^c)-P^n(x,A^c)\|\) | **Open** Claim 3 §11 (no §8.3 drift) |

**Boxed target (general \(P\), multivariate \(\beta_j\)):**

\[
\widetilde{\pi}(\beta \mid y)
\;\propto\;
\exp\Bigl(\sum_{j=1}^{J} \ell_j(\beta_j)\Bigr)
\;
\exp\Bigl(
-\tfrac{1}{2}(\beta - \mathcal{W}\mu_0)^\top
\bigl[
P_b^{\mathrm{stack}}
-
P_b^{\mathrm{stack}} \mathcal{W}\,
\bigl(\Lambda_\gamma + \sum_j \mathcal{W}_j^\top P_b \mathcal{W}_j\bigr)^{-1}
\mathcal{W}^\top P_b^{\mathrm{stack}}
\bigr]
(\beta - \mathcal{W}\mu_0)
\Bigr).
\]

---

## References (in repo)

- `inst/notation.md` — \(\beta_j\), \(\mathcal{W}_j\), \(\gamma\), \(\Psi\), \(D_j\)
- `inst/CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §2 — logit weights, \(P\)-blocks
- `inst/RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0.1, §0.3–§0.4, §2.2–§2.3 — factorization, Lemmas T/G/Z, Prop L1, \(G_j(\gamma)\) (**not** §8–§9 drift)
- `inst/CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §IV.2–IV.3 — \(\bar\omega\), box certification
- `inst/LOGIT_STATIC_TAIL_CERTIFICATION.md` — static \(\pi_{\beta \mid y}(A^c)\) bounds (group-level)
- `inst/LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` — single-group logit score bounds
