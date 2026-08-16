# Symmetric non-Gaussian population posteriors — minorization extension

Draft extension of `MINORIZATION_GAUSSIAN_REFRESH.md` to **logit, probit, and
other GLM** two-block Gibbs samplers when the **marginal posterior** on the
population parameter \(\gamma\) is **symmetric** about a known center. In that
special case, **mean = mode = symmetry center**; the Gaussian refresh pipeline
(Hypotheses 1–2, two-stage \(\varepsilon_d\), penalty ellipsoid \(C_d\)) largely
survives, but **closed Schur forms** for \(\delta(\gamma)\) do not unless one
adds a **local quadratic surrogate**.

**Status.** Draft. Notation matches `inst/notation.md`, `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md`
§10, and the Gaussian minorization note.

**Companion notes.**

- `MINORIZATION_GAUSSIAN_REFRESH.md` — Normal-likelihood closed forms
  (\(\Sigma\), \(\Sigma^\star\), \(\mathcal{L}\), \(\mu_\delta\), \(\Gamma_\gamma\)).
- `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §10.7 — precision dominance
  \(\Lambda_q \prec \Sigma_{\mathrm{upd}}^{-1}\) for general GLMs.
- `LOGIT_SINGLE_GROUP_SAFE_REGION.md` — scalar logit symmetry of
  \(p(1-p)\) at \(\eta=0\).

---

## 1. What is unchanged (two-block Gibbs)

One sweep: \(\beta \mid \gamma, y\), then \(\gamma' \mid \beta\).

**Block 2 (always Gaussian in \(\gamma'\)).**

\[
\gamma' \mid \beta \;\sim\; N\bigl(m(\beta),\, \Sigma^\star\bigr),
\qquad
\Sigma^\star = \Sigma_{\mathrm{upd}},
\]

with \(m(\beta)\) from §8.1 / `notation.md` (scalar:
\(m(\beta)=w\bar\beta+(1-w)\mu_0\); multivariate: GLS form in the Gaussian note).

**Refresh density (Hypotheses 1–2, unchanged).**

\[
Q = N(\mu^\star,\, \Sigma^\star),
\qquad
\mu^\star := E[\gamma \mid y],
\qquad
\Sigma^\star \ \text{fixed from the hierarchy (not the marginal kernel spread).}
\]

**Minorization target (unchanged).**

\[
q(\gamma', \gamma) \;\ge\; \varepsilon_d\, q_Q(\gamma')
\quad
\forall\,\gamma \in C_d,\;\gamma' \in \mathbb R^q,
\]

with \(q(\gamma', \gamma) = \int \kappa(\gamma', \beta)\,\pi(\beta \mid \gamma, y)\,d\beta\).

**Two-stage recipe (unchanged).**

\[
g(\gamma', \gamma) = \log q(\gamma', \gamma) - \log q_Q(\gamma'),
\qquad
\gamma'_\star(\gamma) \in \arg\min_{\gamma'} g(\gamma', \gamma),
\qquad
\delta(\gamma) = g\bigl(\gamma'_\star(\gamma), \gamma\bigr),
\qquad
\varepsilon_d = \exp\Bigl(\min_{\gamma \in C_d} \delta(\gamma)\Bigr).
\]

(Subscript \(d\): certified constant on the level set \(C_d\).)

---

## 2. \(\Sigma^\star\) vs \(\Sigma\) (still distinct)

| Symbol | Meaning |
|---|---|
| \(\Sigma^\star\) | Conditional spread of \(\gamma' \mid \beta\) (Block~2); **refresh** covariance |
| \(\Sigma\) | Second moment \(\mathrm{Cov}(\gamma' \mid \gamma)\) after integrating \(\beta\) |

**Always (any likelihood):**

\[
\Sigma
=
\Sigma^\star
+
\mathrm{Cov}\bigl(m(\beta) \mid \gamma, y\bigr).
\]

**Gaussian closure:** \(q(\gamma', \gamma) = \phi_q(\gamma'; M(\gamma), \Sigma)\) with
exact \(\Sigma\) above.

**Non-Gaussian:** \(q\) is **not** Normal. The same **identity for \(\Sigma\)** still
holds as a **second-moment** decomposition, but \(\Sigma\) is **not** the covariance
parameter of a Gaussian kernel. Do **not** plug this \(\Sigma\) into the Normal
Schur formula for \(\delta(\gamma)\) without justification.

For curvature certification, use the **precision** (when defined)

\[
\Lambda_q(\gamma', \gamma)
:=
-\nabla_{\gamma'}^2 \log q(\gamma', \gamma),
\]

and require \(\Lambda_q \prec \Lambda_Q := \Sigma^{\star-1}\) on \(C_d \times \mathbb R^q\)
(rate note §10.7). **Do not** identify \(\Lambda_q\) with \(\Sigma^{-1}\) where
\(\Sigma = \mathrm{Cov}(\gamma' \mid \gamma)\): that identification holds **only**
in the Gaussian closure (§2.1).

### 2.1 Precision vs variance — be careful

In the **Gaussian** note, every object is Normal and **precision and covariance are
just inverses**. You may rewrite the same formula either way:

\[
\Lambda = \Sigma^{-1},
\qquad
\|u\|_\Sigma^2 = u^\top \Sigma^{-1} u = u^\top \Lambda u,
\qquad
\Lambda_q = \Sigma^{-1},\;
\Lambda_Q = \Sigma^{\star-1}.
\]

In **symmetric non-Gaussian** models this convenience **disappears**. Several
different matrices live in the same notation unless you keep them separate.

| Object | Type | Role | Gaussian special case |
|---|---|---|---|
| \(\Sigma^\star\) | **Variance** (Block~2) | Spread of \(\gamma' \mid \beta\); defines \(Q\) | Same |
| \(\Lambda_Q = \Sigma^{\star-1}\) | **Precision** (refresh) | Hessian of \(-\log q_Q\) | \(=\Sigma^{\star-1}\) |
| \(\Sigma = \mathrm{Cov}(\gamma' \mid \gamma)\) | **Second moment** | Total variance identity §2 | **Equals** kernel cov of \(q\) |
| \(\Lambda_q(\gamma', \gamma)\) | **Precision** (kernel) | Hessian \(-\nabla_{\gamma'}^2 \log q\) | **Constant** \(=\Sigma^{-1}\) |
| \(\Gamma_\gamma\) | **Precision** (outer / \(C_d\)) | Curvature of penalty in \(\gamma\) | From §5.1.1 of Gaussian note |
| \((\Sigma + \Sigma^\star)^{-1}\) | **Precision combination** | Schur / Mahalanobis in **Gaussian** \(\delta\) | Valid **only** when \(q\) is Normal |

**Rule 1 — Inner minimization in \(\gamma'\) (always).** Work in **precisions**:

\[
\nabla_{\gamma'}^2 g(\gamma', \gamma)
=
\Lambda_Q - \Lambda_q(\gamma', \gamma),
\qquad
\Lambda_Q = \Sigma^{\star-1}.
\]

Certify \(\Lambda_q \prec \Lambda_Q\), then solve the **FOC in log-density form**

\[
\nabla_{\gamma'} \log q(\gamma', \gamma)
=
\Lambda_Q\,(\gamma' - \mu^\star)
\]

(or equivalent numeric minimization of \(g\)). Do **not** substitute
\(\Lambda_q := \Sigma(\gamma)^{-1}\) using the second-moment \(\Sigma\) from §2
unless you have already proved Gaussian closure.

**Rule 2 — Second moments are not kernel precisions.** The identity
\(\Sigma = \Sigma^\star + \mathrm{Cov}(m(\beta)\mid\gamma,y)\) is about **variances**.
It does **not** imply

\[
\Lambda_q
\stackrel{?}{=}
\Sigma^{-1},
\qquad
\Lambda_q
\stackrel{?}{=}
\Sigma^{\star-1} + \text{(something simple in $\mathrm{Cov}(m)$)}.
\]

For logit/probit, \(\Lambda_q(\gamma', \gamma)\) depends on **both** \(\gamma'\) and
\(\gamma\). Compute it from the **expectation formula** §2.2.4 (tilted covariance of
\(m(\beta)\)), not from \(\Sigma^{-1}\).

**Rule 3 — Gaussian Schur formulas use a special precision.** The closed
\(\delta(\gamma)\) in `MINORIZATION_GAUSSIAN_REFRESH.md` uses
\((\Sigma + \Sigma^\star)^{-1}\) on \((M(\gamma)-\mu^\star)\). This is **not** the
same as \((\Sigma^{-1} - \Sigma^{\star-1})\) nor as \(\Lambda_Q - \Lambda_q\) at a
single point. **Do not** port \((\Sigma + \Sigma^\star)^{-1}\) to GLMs by plugging in
moment-based \(\Sigma\).

**Rule 4 — Outer ellipsoid \(C_d\) uses a precision.** The matched small set

\[
C_d
=
\{\gamma : (\gamma - \mu^\star)^\top \Gamma_\gamma (\gamma - \mu^\star) \le d\}
\]

is defined by a **precision** \(\Gamma_\gamma\) (curvature of the penalty in
\(\gamma\)), not by \(\mathrm{Cov}(\gamma \mid y)\). Under symmetry,
\(\Gamma_\gamma\) may equal \(-\nabla_\gamma^2 \log \pi_\gamma(\mu^\star)\); its
inverse is a **local variance** surrogate, not the global posterior covariance.

**Rule 5 — Determinant ratio.** \(\tfrac{1}{2}\log(\det \Sigma^\star / \det \Sigma)\)
in the Gaussian \(\delta\) is a **Normalizing constant** identity for **Gaussian**
\(q\). For non-Gaussian \(q\), extract \(\log q - \log q_Q\) **directly** at
\((\gamma'_\star(\gamma), \gamma)\); do not replace the log-det term with moments.

**Summary.** Variances (\(\Sigma^\star\), \(\Sigma\), posterior cov) describe **spread**.
Precisions (\(\Lambda_Q\), \(\Lambda_q\), \(\Gamma_\gamma\)) describe **curvature**
and enter **convexity and FOCs**. Minimization and certification should be written
in **precision / Hessian language**; convert to variances only when the law is
Gaussian and the inverse is exact.

### 2.2 Kernel precision \(\Lambda_q\) as an expectation (tilted measure)

This section derives the **correct** precision of \(\log q\) for the two-block kernel
when Block~2 is Gaussian in \(\gamma'\). It replaces the false shortcut
\(\Lambda_q \stackrel{?}{=} \Sigma^{-1}\) from §2.1.

#### 2.2.1 Tilted RE measure

Fix \(\gamma \in C_d\) and \(\gamma' \in \mathbb R^q\). Write
\(\pi(\beta) := \pi(\beta \mid \gamma, y)\) and

\[
q(\gamma')
:=
\int \kappa(\gamma', \beta)\,\pi(\beta)\,d\beta,
\qquad
\kappa(\gamma', \beta)
=
\phi_q\bigl(\gamma';\, m(\beta),\, \Sigma^\star\bigr).
\]

Define the **tilted** (kernel-weighted) RE posterior

\[
\boxed{
\widetilde\pi(\beta \mid \gamma', \gamma)
\;:=\;
\frac{\kappa(\gamma', \beta)\,\pi(\beta \mid \gamma, y)}{q(\gamma')},
\qquad
\int \widetilde\pi\,d\beta = 1.
}
\]

Expectations under \(\widetilde\pi\) are denoted \(\widetilde E[\cdot]\).

#### 2.2.2 Score identity (gradient of \(\log q\))

Because \(\kappa\) is Gaussian in \(\gamma'\),

\[
\nabla_{\gamma'} \log \kappa(\gamma', \beta)
=
-\,\Sigma^{\star-1}\bigl(\gamma' - m(\beta)\bigr).
\]

For \(q(\gamma') = \int \kappa\,\pi\,d\beta\),

\[
\boxed{
\nabla_{\gamma'} \log q(\gamma', \gamma)
=
\widetilde E\bigl[\nabla_{\gamma'} \log \kappa(\gamma', \beta)\bigr]
=
-\,\Sigma^{\star-1}\Bigl(\gamma' - \widetilde E[m(\beta)]\Bigr).
}
\]

*Proof.* \(\nabla q = \int (\nabla \kappa)\,\pi\,d\beta\), so
\(\nabla \log q = (\int (\nabla \kappa)\pi)/q
= \int (\nabla \log \kappa)\,(\kappa\pi/q)\,d\beta\). \(\square\)

**Contrast with the kernel mean.** The unconditional
\(M(\gamma) = E[m(\beta)\mid\gamma,y]\) uses \(\pi(\beta\mid\gamma,y)\), **not**
\(\widetilde\pi\). Only when \(\gamma'\) is such that \(\widetilde E[m] = M(\gamma)\)
(e.g.\ at a special point) do the two agree.

#### 2.2.3 General Hessian identity (log of an integral)

For any positive integrand \(f(\gamma', \beta)\),

\[
\nabla^2_{\gamma'} \log q
=
\frac{\int (\nabla^2 f)\,\pi\,d\beta}{q}
\;-\;
\Bigl(\nabla_{\gamma'} \log q\Bigr)
\Bigl(\nabla_{\gamma'} \log q\Bigr)^\top.
\]

With \(f = \kappa\) and the score identity,

\[
\boxed{
\nabla_{\gamma'}^2 \log q(\gamma', \gamma)
=
\widetilde E\bigl[\nabla_{\gamma'}^2 \log \kappa\bigr]
\;+\;
\widetilde{\mathrm{Cov}}_{\gamma'}\bigl(\nabla_{\gamma'} \log \kappa\bigr),
}
\]

where \(\widetilde{\mathrm{Cov}}\) is covariance under \(\widetilde\pi\). Equivalently,
\(\nabla^2 \log q = \widetilde E[\nabla^2 \log \kappa + (\nabla \log \kappa)(\nabla \log \kappa)^\top]
- \widetilde E[\nabla \log \kappa]\,\widetilde E[\nabla \log \kappa]^\top\).

#### 2.2.4 Block-2 Gaussian specialization (main formula)

Since \(\nabla_{\gamma'}^2 \log \kappa = -\Sigma^{\star-1}\) is **constant** in \(\beta\),

\[
\widetilde E\bigl[\nabla_{\gamma'}^2 \log \kappa\bigr]
=
-\,\Sigma^{\star-1}.
\]

Also \(\nabla_{\gamma'} \log \kappa = -\Sigma^{\star-1}(\gamma' - m(\beta))\) with
\(\gamma'\) fixed inside the covariance, so

\[
\widetilde{\mathrm{Cov}}_{\gamma'}\bigl(\nabla_{\gamma'} \log \kappa\bigr)
=
\Sigma^{\star-1}\,
\widetilde{\mathrm{Cov}}_{\beta}\bigl(m(\beta)\bigr)\,
\Sigma^{\star-1}.
\]

Therefore the **kernel precision** is

\[
\boxed{
\Lambda_q(\gamma', \gamma)
\;:=\;
-\,\nabla_{\gamma'}^2 \log q(\gamma', \gamma)
\;=\;
\Sigma^{\star-1}
\;-\;
\Sigma^{\star-1}\,
\widetilde{\mathrm{Cov}}_{\beta}\bigl(m(\beta)\,\big|\,\gamma', \gamma\bigr)\,
\Sigma^{\star-1}.
}
\]

**Expectation form (single line).**

\[
\boxed{
\Lambda_q
=
\Sigma^{\star-1}
\;-\;
\Sigma^{\star-1}\,
\widetilde E\Bigl[
\bigl(m - \widetilde E[m]\bigr)
\bigl(m - \widetilde E[m]\bigr)^\top
\Bigr]\,
\Sigma^{\star-1},
\qquad
\widetilde E[\cdot]
=
\int (\cdot)\,\widetilde\pi(\beta\mid\gamma',\gamma)\,d\beta.
}
\]

This is the quantity to use in §6 — **not** \(\Sigma^{-1}\) from §2.

#### 2.2.5 Precision dominance as an expectation inequality

Refresh precision \(\Lambda_Q = \Sigma^{\star-1}\). Then

\[
\Lambda_Q - \Lambda_q
=
\Sigma^{\star-1}\,
\widetilde{\mathrm{Cov}}_{\beta}\bigl(m(\beta)\,\big|\,\gamma', \gamma\bigr)\,
\Sigma^{\star-1}.
\]

So

\[
\boxed{
\Lambda_q(\gamma', \gamma) \;\prec\; \Lambda_Q
\quad\Longleftrightarrow\quad
\widetilde{\mathrm{Cov}}_{\beta}\bigl(m(\beta)\,\big|\,\gamma', \gamma\bigr) \;\succ\; 0
}
\]

(nontrivial spread of the Block~2 mean \(m(\beta)\) under the **tilted** RE law).
Strict inner convexity \(\nabla_{\gamma'}^2 g = \Lambda_Q - \Lambda_q \succ 0\) is
**equivalent** to this tilted covariance being positive definite.

**Important.** The **unconditional** variance identity §2 uses
\(\mathrm{Cov}(m(\beta)\mid\gamma,y)\) under \(\pi(\beta\mid\gamma,y)\), **not**
\(\widetilde{\mathrm{Cov}}\). Both involve \(m(\beta)\), but the **expectation is
over different measures**.

#### 2.2.6 Inner FOC in expectation form

Write \(\widetilde m(\gamma') := \widetilde E[m(\beta)\mid\gamma', \gamma]\) (tilted mean of
the Block~2 map). From §2.2.2,

\[
\nabla_{\gamma'} \log q(\gamma', \gamma)
=
-\,\Sigma^{\star-1}\bigl(\gamma' - \widetilde m(\gamma')\bigr).
\]

The log gap \(g = \log q - \log q_Q\) satisfies

\[
\boxed{
\nabla_{\gamma'} g(\gamma', \gamma)
=
\Sigma^{\star-1}\Bigl(\widetilde m(\gamma') - \mu^\star\Bigr).
}
\]

Note \(\widetilde m(\gamma')\) **depends on \(\gamma'\)** through the tilt
\(\widetilde\pi(\cdot\mid\gamma',\gamma)\); this is **not** the unconditional
\(M(\gamma)\).

The inner critical point \(\gamma'_\star(\gamma)\) solves

\[
\boxed{
\widetilde E\bigl[m(\beta) \,\big|\, \gamma'_\star(\gamma), \gamma\bigr]
=
\mu^\star.
}
\]

At \(\gamma = \mu^\star\) under symmetry (§5), if \(\gamma'_\star(\mu^\star) = \mu^\star\)
then \(\widetilde\pi(\cdot\mid\mu^\star,\mu^\star)\) is symmetric and
\(\widetilde E[m]=\mu^\star\) holds at the center fiber.

#### 2.2.7 Gaussian closure (check)

When \(\pi(\beta\mid\gamma,y)\) and the likelihood are Normal-conjugate so that
\(q(\gamma',\gamma) = \phi_q(\gamma'; M(\gamma), \Sigma)\), one has
\(\Lambda_q = \Sigma^{-1}\) **constant** in \(\gamma'\). The expectation formula
§2.2.4 must agree:

\[
\Sigma^{-1}
=
\Sigma^{\star-1}
\;-\;
\Sigma^{\star-1}\,
\widetilde{\mathrm{Cov}}_{\beta}(m)\,
\Sigma^{\star-1}
\quad\text{at each $(\gamma',\gamma)$ with matching tilt.}
\]

Solving for the tilted covariance,

\[
\widetilde{\mathrm{Cov}}_{\beta}(m)
=
\Sigma^\star - \Sigma^\star \Sigma^{-1} \Sigma^\star
\;=\;
\Sigma^\star\bigl(\Sigma^{-1} - \Sigma^{\star-1}\bigr)^{-1}\Sigma^\star
\quad\text{(when invertible).}
\]

The **unconditional** identity §2 gives
\(\Sigma - \Sigma^\star = \mathrm{Cov}(m \mid \gamma,y)\) under \(\pi\), not under
\(\widetilde\pi\). Inverting \(\Lambda_q = \Sigma^{-1}\) is valid **only** in this
Normal closure.

#### 2.2.8 Monte Carlo evaluation

Given draws \(\beta^{(s)} \sim \pi(\beta\mid\gamma,y)\), importance weights
\(w_s \propto \kappa(\gamma', \beta^{(s)})\), normalized, estimate

\[
\widetilde E[m]
\approx
\sum_s w_s m(\beta^{(s)}),
\qquad
\widetilde{\mathrm{Cov}}(m)
\approx
\sum_s w_s \bigl(m^{(s)} - \widetilde E[m]\bigr)\bigl(m^{(s)} - \widetilde E[m]\bigr)^\top,
\]

then \(\Lambda_q \approx \Sigma^{\star-1} - \Sigma^{\star-1}\widetilde{\mathrm{Cov}}(m)\Sigma^{\star-1}\).
Alternatively differentiate \(\log q\) numerically; the expectation form is the
**variance–precision** certificate aligned with §2.1.

---

## 3. Symmetry hypotheses

Fix a center \(\mu^\star \in \mathbb R^q\) (to be identified with
\(E[\gamma \mid y]\)).

**Hypothesis S1 (posterior symmetry).** The marginal posterior density
\(\pi_\gamma(\gamma) := \pi(\gamma \mid y)\) is **symmetric about \(\mu^\star\)**:

\[
\pi_\gamma(\mu^\star + u)
=
\pi_\gamma(\mu^\star - u)
\qquad
\forall\, u \in \mathbb R^q
\]

(equivalently, \(\gamma - \mu^\star\) has a symmetric distribution).

**Consequences of S1 (when \(\pi_\gamma\) is unimodal / log-concave).**

\[
E[\gamma \mid y]
=
\operatorname{mode}(\pi_\gamma)
=
\mu^\star,
\qquad
\nabla_\gamma \log \pi_\gamma(\gamma)\big|_{\gamma=\mu^\star} = 0.
\]

**Hypothesis S2 (prior symmetry at the center).** Prior on \(\gamma\) is symmetric
about the same center (e.g.\ \(\mu_0 = \mu^\star\) or diffuse symmetric prior).

**Hypothesis S3 (data / design symmetry).** The likelihood and group structure are
invariant under the **sign flip** that makes S1 hold. Examples below (§4).

**Hypothesis S4 (kernel parity at the center).** At \(\gamma = \mu^\star\),

\[
q(\mu^\star + v,\, \mu^\star)
=
q(\mu^\star - v,\, \mu^\star)
\qquad
\forall\, v \in \mathbb R^q,
\]

and \(M(\mu^\star) := E[\gamma' \mid \gamma=\mu^\star] = \mu^\star\).

S4 follows from S1–S3 under the Block~1 / Block~2 maps in §5 when the
involution is consistent (see Proposition 1).

---

## 4. Hypothetical symmetric datasets (logit and probit)

### 4.1 Single-group logit (scalar)

Model as `LOGIT_SINGLE_GROUP_SAFE_REGION.md` §1 with population \(\theta\) and
RE \(b\), prior \(b \mid \theta \sim N(\theta, \tau^2)\), population prior
\(\theta \sim N(0, \sigma_0^2)\) (**centered**). Binomial records with
\(\eta = b\).

**Balanced counts.** For each record \(i\), choose \(y_i = n_i/2\) (requires even
\(n_i\)). Total log-likelihood \(\ell(b)\) satisfies \(\ell(b) = \ell(-b)\) because
each contribution depends on \(b\) only through \(\log(1+e^b)\) and \(\log(1+e^{-b})\)
in a symmetric way at balanced \(y/n = 1/2\).

**Posterior symmetry.** Joint \((\theta, b)\) is symmetric under
\((\theta, b) \mapsto (-\theta, -b)\); marginal \(\pi(\theta \mid y)\) is
symmetric about \(\mu^\star = 0 = E[\theta \mid y]\).

**Same construction, multiple groups.** Groups \(j=1,\ldots,J\), RE \(\beta_j\),
prior \(\beta_j \mid \gamma \sim N(\gamma, \lambda_b^{-1})\), \(\gamma \sim N(0,
\lambda_\gamma^{-1})\). In each group, set \(y_j = n_j/2\) (even \(n_j\)). Then
\(\pi(\gamma \mid y)\) is symmetric about \(0\).

### 4.2 Single-group probit

Replace logit link by probit; keep centered prior and **balanced** binomial outcomes
(\(y_i = n_i/2\)). Likelihood terms depend on \(\Phi(b)\) and \(\Phi(-b)\) symmetrically;
\(\ell(b)=\ell(-b)\). Same posterior symmetry about \(\mu^\star=0\).

### 4.3 Multivariate population (`notation.md`)

Symmetry under **coordinate reflection** about \(\mu^\star\):

\[
(\mu^\star + u,\, \beta,\, y)
\ \text{same density as}\
(\mu^\star - u,\, \beta',\, y)
\]

for an involution on \(\beta\) coupled to \(u\) (e.g.\ \(\beta_j \mapsto -\beta_j\)
with \(\mathcal{W}_j\) and \(D_j\) even in the appropriate sense). A sufficient
recipe: centered prior on \(\gamma\), identical balanced counts in every group,
identical designs across groups, and \(\mu_0 = \mu^\star\).

**Role in development.** These are **certification fixtures**: they isolate
non-Gaussian \(q\) while keeping \(\mu^\star = \text{mode}\) and \(M(\mu^\star)=\mu^\star\).

---

## 5. Algebraic structure induced by symmetry

### 5.1 Involution on the Gibbs state

Assume an involution \(\iota\) on \((\gamma, \beta)\) of the form

\[
\iota(\gamma, \beta) = (2\mu^\star - \gamma,\; R_\beta \beta),
\]

with \(R_\beta\) an orthogonal sign / block reflection on RE coefficients, such that:

1. \(\pi(\beta \mid \gamma, y)\) is invariant: \(\pi(\beta \mid \gamma, y) =
   \pi(R_\beta \beta \mid 2\mu^\star - \gamma, y)\).
2. \(m(R_\beta \beta) = 2\mu^\star - m(\beta)\) (Block~2 mean **pointwise** odd about
   \(\mu^\star\)).
3. Prior and likelihood jointly invariant under \(\iota\).

**Proposition 1 (center fixed point).** Under (1)–(3),

\[
M(\mu^\star) = \mu^\star,
\qquad
q(\mu^\star + v,\, \mu^\star) = q(\mu^\star - v,\, \mu^\star)
\quad
\forall v.
\]

*Proof sketch.* At \(\gamma = \mu^\star\), change variables \(\beta \mapsto R_\beta\beta\)
in \(q(\gamma', \mu^\star) = \int \phi(\gamma'; m(\beta), \Sigma^\star)\pi(\beta \mid
\mu^\star, y)\,d\beta\). Invariance of \(\pi\) and oddness of \(m\) give
\(q(\gamma', \mu^\star) = q(2\mu^\star - \gamma', \mu^\star)\). Take \(\gamma' =
\mu^\star + v\). For \(M(\mu^\star)\), use \(\iota\) on the expectation defining
\(E[\gamma' \mid \gamma]\). \(\square\)

**Corollary (Gaussian note offset).** In the affine notation of
`MINORIZATION_GAUSSIAN_REFRESH.md` §5.1.1,

\[
b_0(y) := M(\mu^\star) - \mu^\star = \mathbf 0
\]

under exact symmetry. The penalty centroid satisfies \(\mu_\delta = \mu^\star\).

### 5.2 Parity of the kernel mean about \(\mu^\star\)

When (1)–(3) hold for **all** \(\gamma\) near \(\mu^\star\),

\[
M(\mu^\star + u) - \mu^\star
=
-\,\bigl(M(\mu^\star - u) - \mu^\star\bigr),
\]

i.e.\ \(M(\gamma) - \mu^\star\) is **odd** in \(u = \gamma - \mu^\star\). Taylor
expansion at \(\mu^\star\):

\[
M(\gamma)
=
\mu^\star + \mathcal{L}\,( \gamma - \mu^\star) + O\bigl(\|\gamma - \mu^\star\|^3\bigr),
\]

with **no constant offset** \(b_0\). The Gaussian affine map is the **first-order**
symmetric approximation; GLMs generally deviate at order \(\|\gamma-\mu^\star\|^3\).

### 5.3 Parity of the log gap at \(\gamma = \mu^\star\)

If \(q(\cdot, \mu^\star)\) is even about \(\mu^\star\) and \(q_Q\) is Normal
\(N(\mu^\star, \Sigma^\star)\), then \(g(\cdot, \mu^\star)\) is **even** in
\(\gamma' - \mu^\star\). If additionally \(g(\cdot, \mu^\star)\) is strictly
convex (§6), the inner minimizer is

\[
\gamma'_\star(\mu^\star) = \mu^\star.
\]

At the center fiber, \(\delta(\mu^\star) = g(\mu^\star, \mu^\star)\) — typically
**not** the closed Normal Schur value unless \(q(\cdot, \mu^\star)\) is Gaussian.

---

## 6. Inner problem without Gaussian \(q\)

### 6.1 Precision dominance (certification condition)

From §2.2.4,

\[
\Lambda_q(\gamma', \gamma)
=
\Sigma^{\star-1}
-
\Sigma^{\star-1}\,
\widetilde{\mathrm{Cov}}_{\beta}\bigl(m(\beta)\,\big|\,\gamma', \gamma\bigr)\,
\Sigma^{\star-1},
\qquad
\Lambda_Q := \Sigma^{\star-1}.
\]

**Require on \(C_d \times \mathbb R^q\):**

\[
\boxed{
\Lambda_q(\gamma', \gamma) \;\prec\; \Lambda_Q
\quad\Longleftrightarrow\quad
\widetilde{\mathrm{Cov}}_{\beta}\bigl(m(\beta)\,\big|\,\gamma', \gamma\bigr) \;\succ\; 0.
}
\]

Then \(g(\cdot, \gamma)\) is **strictly convex** in \(\gamma'\) and the inner
minimizer \(\gamma'_\star(\gamma)\) is unique. At the minimizer,
\(\widetilde E[m(\beta)\mid \gamma'_\star,\gamma] = \mu^\star\) (§2.2.6).

**Normal case (Gaussian closure only).** \(\Lambda_q = \Sigma^{-1}\) **constant**;
dominance \(\Lambda_q \prec \Lambda_Q\) is equivalent to **variance** dominance
\(\Sigma \succ \Sigma^\star\). The tilted covariance in §2.2.4 then agrees with
§2.2.7; this equivalence is **not** available for GLMs.

**GLM case:** evaluate \(\widetilde{\mathrm{Cov}}_{\beta}(m)\) by importance sampling
(§2.2.8) or bound it analytically on \(C_d\).

### 6.2 Log-concavity route (sufficient)

If \((\gamma', \beta) \mapsto \log \kappa(\gamma', \beta) + \log \pi(\beta \mid \gamma, y)\)
is **jointly log-concave** for each fixed \(\gamma \in C_d\), then
\(q(\gamma', \gamma)\) is log-concave in \(\gamma'\) (marginalization). With
\(\Lambda_q \prec \Lambda_Q\), the inner problem is convex. Log-concavity of
\(\pi(\beta \mid \gamma, y)\) holds for many GLM full conditionals; verify per
model.

### 6.3 Computing \(\gamma'_\star(\gamma)\)

No closed form in general. Work in **precisions**, not second-moment inverses.

**Preferred:** solve \(\widetilde E[m(\beta)\mid\gamma',\gamma] = \mu^\star\) (§2.2.6), or
equivalently \(\nabla_{\gamma'} g = 0\), using importance weights from §2.2.8.

Alternatively enforce the log-density FOC

\[
\nabla_{\gamma'} \log q(\gamma', \gamma)
=
\Lambda_Q\,(\gamma' - \mu^\star),
\qquad
\Lambda_Q = \Sigma^{\star-1},
\]

which is the same as \(\gamma' - \widetilde m(\gamma') = \gamma' - \mu^\star\), i.e.\
\(\widetilde m(\gamma') = \mu^\star\).

**Optional initializer (Gaussian surrogate only).** If you **locally** approximate
\(\log q(\cdot, \gamma)\) by a Normal with precision \(\widehat\Lambda_q(\gamma)\)
(e.g.\ \(-\nabla_{\gamma'}^2 \log q\) at a provisional mode), then

\[
\gamma'_\star(\gamma)
\approx
\bigl(\Lambda_Q - \widehat\Lambda_q(\gamma)\bigr)^{-1}
\Bigl(\Lambda_Q \mu^\star - \widehat\Lambda_q(\gamma)\, \widehat m_q(\gamma)\Bigr),
\]

where \(\widehat m_q\) is the provisional kernel mode — **not** \(M(\gamma)\) from
moments unless verified. Refine by Newton on the **true** FOC. **Do not** set
\(\widehat\Lambda_q = \Sigma(\gamma)^{-1}\) from §2 without Gaussian closure.

---

## 7. Outer problem and penalty ellipsoid \(C_d\)

### 7.1 What breaks from the Gaussian note

`MINORIZATION_GAUSSIAN_REFRESH.md` §5.1 gives a closed \(\delta(\gamma)\) using
**both** a log-determinant ratio **and** a Mahalanobis term with
\((\Sigma + \Sigma^\star)^{-1}\). That entire display is valid **only** when
\(q = \phi_q(\cdot; M(\gamma), \Sigma)\). For GLMs:

- compute \(\delta(\gamma) = g(\gamma'_\star(\gamma), \gamma)\) after the **inner**
  precision-based solve (§6);
- do **not** replace \(\Lambda_q\) by \(\Sigma^{-1}\) or port
  \((\Sigma + \Sigma^\star)^{-1}\) using moment matrices from §2.

### 7.2 Symmetric center simplifies the penalty geometry

Under S1–S4, \(b_0 = 0\) and \(\mu_\delta = \mu^\star\). The **refresh center** and
**penalty centroid** coincide.

**Penalty ellipsoid (matched definition, §5.1.2 of Gaussian note).** Fix \(d>0\).
Define \(\Gamma_\gamma \succ 0\) (see §7.3) and

\[
\boxed{
C_d
=
\bigl\{
\gamma :
(\gamma - \mu^\star)^\top \Gamma_\gamma (\gamma - \mu^\star) \le d
\bigr\}.
}
\]

On \(C_d\), if one has a **quadratic lower bound**

\[
\delta(\gamma)
\;\ge\;
c_0 - \tfrac{1}{2}(\gamma - \mu^\star)^\top \Gamma_\gamma (\gamma - \mu^\star),
\qquad
c_0 = \text{const},
\]

then \(\min_{C_d} \delta \ge c_0 - d/2\) and
\(\varepsilon_d \ge \exp(c_0 - d/2)\). Equality holds when \(\delta\) **is** that
quadratic (Gaussian case) or when the worst \(\gamma\) lies on \(\partial C_d\) and
achieves the bound.

### 7.3 Choosing \(\Gamma_\gamma\) for \(C_d\) (precision only)

All options below are **precisions** (Hessians / curvature). Inverses are local
variance **surrogates**, not objects to substitute into inner formulas.

**Option A — Posterior precision at the mode (recommended default).**

\[
\Gamma_\gamma
:=
-\nabla_\gamma^2 \log \pi_\gamma(\gamma)\big|_{\gamma = \mu^\star}.
\]

Under log-concavity of \(\pi_\gamma\), this is the curvature of \(-\log \pi\) at
the symmetry center. **Gaussian special case:** if \(\pi_\gamma =
N(\mu^\star, \Sigma_{\mathrm{post}})\), then \(\Gamma_\gamma = \Sigma_{\mathrm{post}}^{-1}\).

**Option B — Outer Hessian of \(\delta\) at \(\mu^\star\).**

\[
\Gamma_\gamma
:=
-\nabla_\gamma^2 \delta(\gamma)\big|_{\gamma = \mu^\star},
\]

after \(\delta(\gamma) = g(\gamma'_\star(\gamma), \gamma)\) is defined by the
**true** inner problem (§6). This matches the penalty curvature actually used in
certification.

**Option C — Do not use (Gaussian-only trap).** The formula
\(\mathcal{L}^\top (\Sigma + \Sigma^\star)^{-1} \mathcal{L}\) from the Normal note
mixes a **Schur precision** valid only for Gaussian \(q\) with a second-moment
\(\Sigma\). **Avoid** for GLMs unless you have proved Gaussian closure.

**Option D — Empirical.** Fit a quadratic to numerically computed \(\delta(\gamma)\)
near \(\mu^\star\); the fitted second-order coefficient is a precision estimate.

For **exact** symmetric fixtures (§4), Options A and B should agree to first order
at \(\mu^\star\) when the inner problem is well approximated locally.

### 7.4 Level sets of \(\delta\) (canonical \(C_d\))

The outer certification plugs \(\gamma\) into the **inner log gap**

\[
\delta(\gamma)
\;=\;
g\bigl(\gamma'_\star(\gamma),\,\gamma\bigr),
\qquad
\varepsilon_d
\;=\;
\exp\Bigl(
\min_{\gamma \in C_d} \delta(\gamma)
\Bigr).
\]

The most direct definition of \(C_d\) is therefore a **level set of \(\delta\)** (or
of the **loss from the center fiber**), not a quadratic or deviance proxy unless
those coincide with \(\delta\).

Assume \(\delta(\mu^\star)\) is a **local maximum** of \(\delta\) in \(\gamma\)
(true under S4 + inner convexity when \(\gamma'_\star(\mu^\star)=\mu^\star\); in the
Gaussian note, \(\delta = c_0 - \text{quadratic}\)). Define the **gap loss**

\[
\boxed{
L_\delta(\gamma)
\;:=\;
\delta(\mu^\star) - \delta(\gamma)
\;\ge\; 0
\quad\text{near }\mu^\star.
}
\]

Fix a radius \(d > 0\). The **\(\delta\)-level small set** is

\[
\boxed{
C_d
\;:=\;
\bigl\{
\gamma \in \mathbb R^q :
L_\delta(\gamma) \le \tfrac{d}{2}
\bigr\}
\;=\;
\bigl\{
\gamma :
\delta(\gamma) \ge \delta(\mu^\star) - \tfrac{d}{2}
\bigr\}.
}
\]

Equivalently, \(C_d = \{\gamma : \delta(\gamma) \ge \tau\}\) with threshold
\(\tau := \delta(\mu^\star) - d/2\). The boundary is the **level contour**
\(\{\gamma : \delta(\gamma) = \tau\}\).

**Closed certified \(\varepsilon_d\) (by construction).** If \(C_d\) is exactly the
sublevel set above and \(\delta\) is continuous, then

\[
\boxed{
\min_{\gamma \in C_d} \delta(\gamma)
\;=\;
\delta(\mu^\star) - \tfrac{d}{2},
\qquad
\varepsilon_d
\;=\;
\exp\Bigl(
\delta(\mu^\star) - \tfrac{d}{2}
\Bigr),
}
\]

with the worst case on **any** point of \(\{\gamma : L_\delta(\gamma) = d/2\}\). No
quadratic geometry is required — only that \(\delta(\mu^\star)\) is the best value
on the neighborhood and worse values lie on the outer contour.

**Compactness (Rosenthal).** Require \(C_d\) compact, e.g.\ \(L_\delta(\gamma) \to +\infty\)
as \(\|\gamma\|\to\infty\) (typical if \(\delta \to -\infty\) far from \(\mu^\star\)).
If \(\{\delta \ge \tau\}\) is non-compact, replace by a **truncated** level set
\(C_d \cap \{\|\gamma-\mu^\star\| \le R\}\) or intersect with a drift sublevel set.

**Enlarging \(\widetilde C_d\) for small tail mass.** For Proposition R2
(`RESTRICTED_GIBBS_MINORIZATION_TV.md` §7.3), treat the certified level set as
\(\widetilde C_d\): increase the radius \(d\) (or deviance threshold) until
\(\pi(\widetilde C_d^c) < \delta_2\). Intersect with the **rate sublevel set**
\(E_{\lambda^\star} = \{\lambda_{\max}(A(\beta)) \le \lambda^\star\}\) (Gibbs rate
matrix \(A(\beta)\), `RESTRICTED_GIBBS_MINORIZATION_TV.md` §7.3 Step 1) to form
\(C_d = E_{\lambda^\star} \cap \widetilde C_d\); then \(\pi(C_d^c) \le \delta_1 +
\delta_2\). Minorization on \(\widetilde C_d\) restricts to \(C_d \subseteq
\widetilde C_d\).

**Gap deviance form.** With

\[
D_\delta(\gamma)
\;:=\;
2\,L_\delta(\gamma)
\;=\;
2\bigl[\delta(\mu^\star) - \delta(\gamma)\bigr],
\]

one has \(C_d = \{\gamma : D_\delta(\gamma) \le d\}\). This is the minorization analogue
of posterior deviance §7.5.1, but tied to the **function actually minimized** in
the outer step.

**Recovering the Gaussian ellipsoid (Gaussian note).** When
\(\delta(\gamma) = c_0 - \tfrac{1}{2}\|M(\gamma)-\mu^\star\|_{\Gamma_M}^2\),

\[
L_\delta(\gamma)
=
\tfrac{1}{2}
\Bigl(
\|M(\gamma)-\mu^\star\|_{\Gamma_M}^2
-
\|M(\mu^\star)-\mu^\star\|_{\Gamma_M}^2
\Bigr),
\qquad
C_d
=
\bigl\{
\|M(\gamma)-\mu^\star\|_{\Gamma_M}^2
\le
\|M(\mu^\star)-\mu^\star\|_{\Gamma_M}^2 + d
\bigr\}.
\]

The completed-square rewrite uses \(\mu_\delta\) from §5.1.1 but the anchor remains
\(\mu^\star\) via \(\delta(\mu^\star)\). When \(M(\mu^\star)=\mu^\star\),
\(\|M(\mu^\star)-\mu^\star\|=0\) and
\(C_d = \{(\gamma-\mu^\star)^\top \Gamma_\gamma (\gamma-\mu^\star) \le d\}\).
The quadratic is a **special case** of a \(\delta\) level set, not the other way
around.

**Non-Gaussian symmetric GLMs.** Compute \(\delta(\gamma)\) numerically (§6.3);
plot or search on \(\{\gamma : \delta(\gamma) = \tau\}\). Shape follows the **true**
outer objective; \(\varepsilon_d = \exp(\delta(\mu^\star) - d/2)\) remains closed
**once** \(C_d\) is defined as this level set.

### 7.5 Deviance and posterior level sets (alternative geometry)

The penalty ellipsoid §7.2 and posterior deviance §7.5.1 are **alternative**
geometries when you do **not** yet have \(\delta(\gamma)\) on a grid, or when you
want \(C_d\) to follow \(\pi(\gamma\mid y)\) instead of \(\delta\). The **canonical**
minorization choice is §7.4: level sets of \(\delta\) itself.

#### 7.5.1 Posterior relative deviance (recommended GLM shape)

Assume \(\pi_\gamma\) is unimodal and log-concave with unique mode at \(\mu^\star\)
(S1). Define the **relative posterior deviance**

\[
\boxed{
D_\pi(\gamma)
\;:=\;
2\Bigl[
\log \pi_\gamma(\mu^\star) - \log \pi_\gamma(\gamma)
\Bigr]
\;\ge\; 0.
}
\]

This is the usual “twice the log posterior ratio” to the mode (Wilks/LRT form at
the population block). The **deviance-based small set** is

\[
\boxed{
C_d^{(\pi)}
\;:=\;
\bigl\{
\gamma \in \mathbb R^q :
D_\pi(\gamma) \le d
\bigr\}.
}
\]

**Scalar signed deviance residual** (\(q=1\)):

\[
\boxed{
r_{\mathrm{dev}}(\gamma)
\;:=\;
\operatorname{sign}(\gamma - \mu^\star)\,
\sqrt{D_\pi(\gamma)},
\qquad
C_d^{(\pi)}
=
\bigl\{
\gamma :
\bigl|r_{\mathrm{dev}}(\gamma)\bigr| \le \sqrt{d}
\bigr\}.
}
\]

Under S1, \(D_\pi(\mu^\star + u) = D_\pi(\mu^\star - u)\): contours are **symmetric**
about \(\mu^\star\) but need **not** be ellipses when \(\pi_\gamma\) is non-Gaussian.

**Laplace link (when the quadratic is only local).** Taylor expansion at \(\mu^\star\):

\[
-\log \pi_\gamma(\gamma)
\;=\;
-\log \pi_\gamma(\mu^\star)
\;+\;
\tfrac{1}{2}(\gamma - \mu^\star)^\top \Gamma_{\mathrm{post}} (\gamma - \mu^\star)
\;+\;
O\bigl(\|\gamma - \mu^\star\|^3\bigr),
\]

with \(\Gamma_{\mathrm{post}} := -\nabla_\gamma^2 \log \pi_\gamma(\mu^\star)\) (§7.3
Option A). Hence

\[
D_\pi(\gamma)
\;=\;
(\gamma - \mu^\star)^\top \Gamma_{\mathrm{post}} (\gamma - \mu^\star)
\;+\;
O\bigl(\|\gamma - \mu^\star\|^3\bigr).
\]

So **small** \(d\): \(C_d^{(\pi)}\) agrees with the penalty ellipsoid
\(C_d = \{(\gamma-\mu^\star)^\top \Gamma_{\mathrm{post}}(\gamma-\mu^\star) \le d\}\).
**Moderate/large** \(d\): deviance contours follow the **true** posterior geometry;
the quadratic ellipsoid can **mis-size** \(C_d\) (too fat in tails or too thin in
curved directions).

**Calibrating \(d\).** Under a Laplace approximation,
\(D_\pi(\gamma) \approx (\gamma-\mu^\star)^\top \Gamma_{\mathrm{post}}(\gamma-\mu^\star)\)
and \(\Gamma_{\mathrm{post}}^{1/2}(\gamma-\mu^\star)\) is approximately standard
Normal at the mode, so \(D_\pi \sim \chi^2_q\) locally. One may take
\(d = \chi^2_{q,\,1-\alpha}\) (asymptotic \((1-\alpha)\)-contour). This calibrates
**posterior mass**, not the minorization constant directly.

#### 7.5.2 Gap deviance (= \(\delta\) level set)

With \(L_\delta(\gamma) = \delta(\mu^\star) - \delta(\gamma)\) as in §7.4,

\[
\boxed{
D_\delta(\gamma)
\;:=\;
2\,L_\delta(\gamma)
\;=\;
2\bigl[\delta(\mu^\star) - \delta(\gamma)\bigr]
\;\ge\; 0,
\qquad
C_d
\;=\;
\bigl\{
\gamma :
D_\delta(\gamma) \le d
\bigr\}.
}
\]

This is **identical** to the canonical §7.4 definition (not a separate choice).
When \(\delta\) is quadratic about \(\mu^\star\), \(D_\delta(\gamma) =
(\gamma-\mu^\star)^\top \Gamma_\gamma (\gamma-\mu^\star)\) and \(C_d\) is the
penalty ellipsoid §7.2.

#### 7.5.3 Which geometry for \(C_d\)?

| Small set | Definition | \(\varepsilon_d\) on \(C_d\) |
|---|---|---|
| **\(\delta\) level set §7.4** | \(\{\delta \ge \delta(\mu^\star) - d/2\}\) | **Closed:** \(\exp(\delta(\mu^\star)-d/2)\) |
| Penalty ellipsoid §7.2 | Quadratic **when** \(\delta\) is quadratic | Same as §7.4 in Gaussian note |
| \(C_d^{(\pi)}\) §7.5.1 | \(\{D_\pi \le d\}\) from \(\pi_\gamma\) | **Numeric** \(\min_{C_d}\delta\) in general |
| Intersection with drift | \(C_d \cap \{V_{\mathrm{drift}} \le \cdots\}\) | Numeric / intersect bounds |

**Inner problem unchanged.** Whichever \(C_d\) defines the **outer** set, the inner
step still uses **precisions** \(\Lambda_q\), \(\Lambda_Q\) and tilted expectations
(§2.2) — not deviance residuals.

#### 7.5.4 Deviance residuals vs precisions (do not mix)

| Object | Type | Role in \(C_d\) |
|---|---|---|
| \(D_\pi\), \(r_{\mathrm{dev}}\) | **Unitless** log-ratio / residual | Define **shape** of \(C_d^{(\pi)}\) |
| \(\Gamma_{\mathrm{post}}\), \(\Gamma_\gamma\) | **Precision** | Local quadratic **approximation** to \(D_\pi\) or \(\delta\) |
| \(\Lambda_q\), \(\Lambda_Q\) | **Precision** | Inner convexity / \(\gamma'_\star\) (§6) |

Deviance contours are **not** a substitute for \(\Lambda_q\). They answer “how far
from \(\mu^\star\) in posterior or gap units?” — not “how curved is \(\log q\) in
\(\gamma'\)?”.

#### 7.5.5 Outer minimum when \(C_d \neq \delta\) level set

On \(C_d^{(\pi)}\) (posterior deviance only), \(\delta(\gamma)\) is **not** constant
on \(\partial C_d^{(\pi)}\). Then

\[
\gamma_\star \in \arg\min_{\gamma \in C_d^{(\pi)}} \delta(\gamma),
\qquad
\varepsilon_d = \exp\bigl(\delta(\gamma_\star)\bigr),
\]

with a **numeric** search on the boundary \(\{\gamma : D_\pi(\gamma) = d\}\) (or over
a grid / MCMC hull). **Lower bound:** if \(\delta(\gamma) \ge \delta(\mu^\star) -
\tfrac{1}{2}(\gamma-\mu^\star)^\top \Gamma_\gamma(\gamma-\mu^\star)\) on
\(C_d^{(\pi)}\), then

\[
\min_{C_d^{(\pi)}} \delta
\;\ge\;
\delta(\mu^\star) - \frac{1}{2}
\max_{\gamma \in C_d^{(\pi)}}
(\gamma - \mu^\star)^\top \Gamma_\gamma (\gamma - \mu^\star),
\]

and the max Mahalanobis term on a **non-elliptic** deviance contour requires
optimization (unless \(d\) is small enough for Laplace).

**Symmetric balanced logit/probit (§4).** \(D_\pi\) is explicit in \(\gamma\) from
the marginal (integrate numerically); \(C_d^{(\pi)}\) is symmetric but typically
**narrower than** a Gaussian ellipsoid with the same \(d\) in the tails — a useful
conservative choice for coupling if the quadratic ellipsoid was too large.

---

## 8. End-to-end certification checklist (symmetric GLM)

| Step | Action |
|---|---|
| 1 | Verify S1–S3 for the fixture (§4) or numerically check \(\pi_\gamma(\mu^\star+u)=\pi_\gamma(\mu^\star-u)\). |
| 2 | Set \(\mu^\star = E[\gamma\mid y]\) (MCMC or quadrature); confirm \(\approx\) posterior mode. |
| 3 | Fix \(\Sigma^\star = \Sigma_{\mathrm{upd}}\); \(Q = N(\mu^\star, \Sigma^\star)\); record \(\Lambda_Q = \Sigma^{\star-1}\). |
| 4 | On \(C_d \times \mathbb R^q\), verify **precision dominance** \(\Lambda_q \prec \Lambda_Q\) (§6.1) — not \(\Sigma \succ \Sigma^\star\) unless Gaussian. |
| 5 | For each \(\gamma \in C_d\): solve inner FOC in **precision form** (§6.3); record \(\delta(\gamma)=g(\gamma'_\star,\gamma)\). |
| 6 | Define \(C_d\): **canonical** \(\delta\) level set §7.4; or proxy §7.5
  (\(\pi\) deviance) if \(\delta\) not yet on a grid. |
| 7 | Minimize \(\delta\) over \(C_d\) (numeric on deviance boundary if needed); set
  \(\varepsilon_d = \exp(\min \delta)\). |
| 8 | Plug \(\varepsilon_d\) into Rosenthal coupling (rate note §10.1). |
| 9 | For existence (R2): enlarge \(\widetilde C_d\); intersect with
  \(E_{\lambda^\star}=\{\lambda_{\max}(A(\beta))\le\lambda^\star\}\); define
  restricted Gibbs on \(C_d\)
  (`RESTRICTED_GIBBS_MINORIZATION_TV.md` §7.3). |

**Center fiber shortcut.** When S4 holds and inner convexity at \(\gamma=\mu^\star\),
check \(\gamma'_\star(\mu^\star)=\mu^\star\) and \(\delta(\mu^\star) = \log q(\mu^\star,
\mu^\star) - \log q_Q(\mu^\star)\) directly.

---

## 9. Comparison to the Gaussian note

| Feature | Normal conjugate (`MINORIZATION_GAUSSIAN_REFRESH.md`) | Symmetric non-Gaussian (this note) |
|---|---|---|
| \(q(\gamma', \gamma)\) | Gaussian | Non-Gaussian mixture / integral |
| Inner curvature | \(\Lambda_q = \Sigma^{-1}\) (constant) | \(\Lambda_q(\gamma', \gamma)\) from **Hessian of \(\log q\)** |
| Dominance check | \(\Sigma \succ \Sigma^\star\) **iff** \(\Lambda_q \prec \Lambda_Q\) | **Only** \(\Lambda_q \prec \Lambda_Q\) (moment inequality insufficient) |
| \(\mu^\star\) vs mode | Equal in general | Equal **by symmetry** (S1) |
| \(M(\mu^\star)\) | \(\neq \mu^\star\) in general | \(= \mu^\star\) under S4 |
| \(b_0\), \(\mu_\delta\) | Often \(\neq 0\) | \(b_0=0\), \(\mu_\delta=\mu^\star\) |
| Inner \(\gamma'_\star\) | Closed FOC via precisions | Numeric FOC; **precisions only** |
| \(\delta(\gamma)\) | Schur + \((\Sigma+\Sigma^\star)^{-1}\) term | Numeric \(g\) at inner min; no Schur port |
| \(C_d\) | Penalty ellipsoid = \(\delta\) level set when \(\delta\) quadratic | General \(\delta\) level set §7.4; \(\pi\) deviance proxy §7.5 |
| Closed \(\varepsilon_d\) | \(\exp(\delta(\mu^\star)-d/2)\) on §7.4 / matched ellipsoid | Same on §7.4; numeric if \(C_d\) uses only \(D_\pi\) |

---

## 10. Open items

1. **Sufficient conditions on balanced GLM designs** for S1 in multivariate
   \(\gamma\) with general \(\mathcal{W}_j\), \(D_j\) (`notation.md`).
2. **Uniform** \(\Lambda_q \prec \Lambda_Q\) on \(C_d\) for logit/probit without
   pointwise numeric search (curvature **precision** envelopes from score bounds §1–§6
   of rate note — not variance bounds on \(\Sigma\)).
3. **Poisson / cloglog:** symmetry requires different involutions (link asymmetry);
   balanced constructions are link-specific.
4. **R implementation:** `data-raw/` script on §4 fixtures computing \(\delta\),
   \(\varepsilon_d\), and dominance checks (not in `tests/testthat/` without approval).
5. **Drift intersection:** Foster–Lyapunov \(C_d\) at prior \(\mu_0\) vs penalty
   ellipsoid at \(\mu^\star\) when \(\mu_0 \neq \mu^\star\) even under approximate
   symmetry.
6. **Deviance \(C_d\):** closed-form \(D_\pi(\gamma)\) for balanced logit/probit
   fixtures (§4); compare \(\varepsilon_d\) on \(C_d^{(\pi)}\) vs penalty ellipsoid.

---

## 11. References

- `MINORIZATION_GAUSSIAN_REFRESH.md` — Hypotheses 1–2, §5.1.1–§5.1.2, §6.
- `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` — §8 drift, §10 coupling, §10.7 GLM precision.
- `notation.md` — \(\gamma\), \(\beta_j\), \(\mathcal{W}_j\), \(\Psi\).
- Lindley & Smith (1972); Rosenthal (1995) — cited in rate / coupling notes.
