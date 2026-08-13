# Gaussian refresh minorization (two-block Gibbs)

Standalone note on **Rosenthal minorization** with a **Gaussian refresh
density** \(Q\), in the **full Gaussian** (Normal-likelihood) case. Coupling
context: `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §10; notation aligned with
`inst/notation.md` and that draft.

**Population parameter** \(\gamma, \gamma' \in \mathbb R^q\). One sweep:
\(\beta \mid \gamma, y\), then \(\gamma' \mid \beta\).

---

## 1. Kernel and refresh densities

### 1.1 Block-2 conditional (always Gaussian in \(\gamma'\))

Given \(\beta\), the population draw in the two-block Gibbs sampler is

\[
\gamma' \mid \beta \;\sim\; N\bigl(m(\beta),\, \Sigma^\star\bigr),
\qquad
\kappa(\gamma', \beta) = \phi_q\bigl(\gamma';\, m(\beta),\, \Sigma^\star\bigr),
\]

with \(m(\beta)\) the Block~2 mean and \(\Sigma^\star \succ 0\) the **fixed**
update covariance from the hierarchy. **Scalar** §8.1 of the rate note:
\(m(\beta) = w\bar\beta + (1-w)\mu_0\),
\(\Sigma^\star = \sigma_\gamma^2 = (J\lambda_b + \lambda_\gamma)^{-1}\).
**Multivariate** (`notation.md`): with
\(P_\gamma := \Lambda_\gamma + \sum_j \mathcal{W}_j^\top \Psi^{-1}\mathcal{W}_j\),

\[
m(\beta)
=
P_\gamma^{-1}
\Bigl(
\Lambda_\gamma \mu_0
+
\sum_{j=1}^{J} \mathcal{W}_j^\top \Psi^{-1} \beta_j
\Bigr),
\qquad
\Sigma^\star = P_\gamma^{-1}.
\]

### 1.2 Marginal one-step kernel (Normal closure)

When \(\pi(\beta \mid \gamma, y)\) is Gaussian (Normal likelihood, conjugate
RE), integrating \(\beta \mid \gamma,y\) gives

\[
q(\gamma', \gamma)
=
\int \kappa(\gamma', \beta)\,\pi(\beta \mid \gamma, y)\,d\beta
=
\phi_q\bigl(\gamma';\, M(\gamma),\, \Sigma\bigr),
\]

\[
M(\gamma) := E[\gamma' \mid \gamma],
\qquad
\Sigma := \mathrm{Cov}(\gamma' \mid \gamma).
\]

**Total variance** (independent of the likelihood family for the second moment
identity):

\[
\boxed{
\Sigma
=
\Sigma^\star
+
\mathrm{Cov}\bigl(m(\beta) \mid \gamma, y\bigr).
}
\]

Scalar Normal closed form (rate note §10.5):

\[
\Sigma
=
\sigma_\gamma^2
+
w^2\,\mathrm{Var}(\bar\beta \mid \gamma, y).
\]

So \(\Sigma \succ \Sigma^\star\) whenever \(\mathrm{Cov}(m(\beta)\mid\gamma,y)\)
is nontrivial (equivalently \(\Lambda_q := \Sigma^{-1} \prec \Sigma^{\star-1}
=:\Lambda_Q\): the transition kernel has **smaller precision** than the block-2
conditional).

---

## 2. Hypotheses on the refresh density \(Q\)

Rosenthal’s minorization uses a **single** refresh measure \(Q\), valid for every
current state \(\gamma\) in the small set \(C_d\). We specify \(Q\) as multivariate
Normal with two hypotheses.

**Hypothesis 1 (covariance).** The refresh variance–covariance equals the
**block-2 conditional** variance–covariance (population draw given \(\beta\)):

\[
\boxed{
\mathrm{Cov}_Q \;:=\; \Sigma^\star.
}
\]

**Hypothesis 2 (mean).** The refresh mean equals the **posterior mean vector** of
the population parameter (marginal on \(\gamma\), given data):

\[
\boxed{
\mu^\star \;:=\; E[\gamma \mid y]
\;=\;
E_{\pi_\gamma}[\gamma]
\;=\;
\int \gamma\,\pi_\gamma(\gamma)\,d\gamma.
}
\]

This is **not** the prior hypermean \(\mu_0\) (used in drift / \(C_d\) in the rate
note) and **not** the state-dependent kernel mean \(M(\gamma) = E[\gamma'\mid\gamma]\).

**Refresh density.**

\[
Q(\gamma') = \phi_q(\gamma';\, \mu^\star,\, \Sigma^\star),
\qquad
q_Q(\gamma') = \phi_q(\gamma';\, \mu^\star,\, \Sigma^\star).
\]

**Minorization target.** Find \(\varepsilon_d \in (0,1]\) such that

\[
q(\gamma', \gamma) \;\ge\; \varepsilon_d\, q_Q(\gamma')
\qquad
\forall\,\gamma \in C_d,\;\gamma' \in \mathbb R^q.
\]

Subscript \(d\) marks dependence on the level-set radius defining \(C_d\); the
certified constant is \(\varepsilon_d = \exp(\min_{\gamma\in C_d}\delta(\gamma))\).

Equivalently, with the **log gap**

\[
g(\gamma', \gamma)
\;:=\;
\log q(\gamma', \gamma) - \log q_Q(\gamma'),
\]

require \(g(\gamma', \gamma) \ge \log\varepsilon_d\) for all \((\gamma,\gamma')\).
Certified constant:

\[
\log \varepsilon_d
=
\min_{\gamma \in C_d,\;\gamma' \in \mathbb R^q} g(\gamma', \gamma).
\]

---

## 3. Convexity of \(g(\cdot, \gamma)\)

Write

\[
\log q(\gamma', \gamma)
=
-\tfrac{q}{2}\log(2\pi)
-\tfrac{1}{2}\log\det\Sigma
-\tfrac{1}{2}(\gamma' - M(\gamma))^\top \Sigma^{-1}(\gamma' - M(\gamma)),
\]

\[
\log q_Q(\gamma')
=
-\tfrac{q}{2}\log(2\pi)
-\tfrac{1}{2}\log\det\Sigma^\star
-\tfrac{1}{2}(\gamma' - \mu^\star)^\top \Sigma^{\star-1}(\gamma' - \mu^\star).
\]

Then

\[
\nabla_{\gamma'} g(\gamma', \gamma)
=
-\,\Sigma^{-1}\bigl(\gamma' - M(\gamma)\bigr)
+
\Sigma^{\star-1}\bigl(\gamma' - \mu^\star\bigr),
\]

\[
\nabla_{\gamma'}^2 g
=
\Sigma^{\star-1} - \Sigma^{-1}
=
\Lambda_Q - \Lambda_q.
\]

By §1.2, \(\Sigma \succ \Sigma^\star\), hence \(\Lambda_q \prec \Lambda_Q\) and

\[
\nabla_{\gamma'}^2 g \succ 0.
\]

**Non-Gaussian \(q\) (expectation form).** When \(q\) is not Normal, \(\Lambda_q\) is
still defined by \(-\nabla_{\gamma'}^2 \log q\). For Block~2 Gaussian \(\kappa\),

\[
\Lambda_q
=
\Sigma^{\star-1}
-
\Sigma^{\star-1}\,
\widetilde{\mathrm{Cov}}_{\beta}\bigl(m(\beta)\mid\gamma',\gamma\bigr)\,
\Sigma^{\star-1}
\]

under the tilted measure \(\widetilde\pi \propto \kappa\,\pi(\beta\mid\gamma,y)\); see
`MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` §2.2. In the **Gaussian closure** above,
\(\widetilde{\mathrm{Cov}}_{\beta}(m)\) collapses so that \(\Lambda_q = \Sigma^{-1}\)
constant.

So for each fixed \(\gamma\), \(g(\cdot,\gamma)\) is **strictly convex** on
\(\mathbb R^q\); any critical point is the **global** minimizer.

---

## 4. Inner minimum: \(\gamma'_\star(\gamma)\)

For fixed \(\gamma\), solve \(\nabla_{\gamma'} g(\gamma', \gamma) = 0\):

\[
\Sigma^{\star-1}(\gamma' - \mu^\star)
=
\Sigma^{-1}(\gamma' - M(\gamma)).
\]

\[
\boxed{
\gamma'_\star(\gamma)
=
\bigl(\Sigma^{\star-1} - \Sigma^{-1}\bigr)^{-1}
\Bigl(
\Sigma^{\star-1}\,\mu^\star
-
\Sigma^{-1}\,M(\gamma)
\Bigr).
}
\]

**Equivalent forms** (same FOC):

\[
\gamma'_\star(\gamma)
=
\bigl(\Sigma - \Sigma^\star\bigr)^{-1}
\Bigl(
\Sigma\,\Sigma^{\star-1}\,\mu^\star
-
\Sigma^\star\,\Sigma^{-1}\,M(\gamma)
\Bigr),
\]

\[
\gamma'_\star(\gamma)
=
\mu^\star + \Sigma^\star\bigl(\Sigma - \Sigma^\star\bigr)^{-1}
\bigl(M(\gamma) - \mu^\star\bigr).
\]

**Scalar** (\(q=1\), write \(V=\Sigma\), \(\sigma_\star^2 = \Sigma^\star\)):

\[
\boxed{
\gamma'_\star(\gamma)
=
\frac{V\,\mu^\star - \sigma_\star^2\, M(\gamma)}{V - \sigma_\star^2}.
}
\]

---

## 5. Value of the log gap at the inner minimum

Define the **inner log deficit**

\[
\delta(\gamma)
\;:=\;
g\bigl(\gamma'_\star(\gamma),\, \gamma\bigr)
=
\log q\bigl(\gamma'_\star(\gamma) \mid \gamma\bigr)
-
\log q_Q\bigl(\gamma'_\star(\gamma)\bigr).
\]

Because \(g(\cdot,\gamma)\) is convex, \(g(\gamma',\gamma) \ge \delta(\gamma)\) for
**every** \(\gamma' \in \mathbb R^q\).

### 5.1 Closed form (multivariate)

Write \(\Gamma_M := (\Sigma + \Sigma^\star)^{-1}\) (Mahalanobis weight on the
\(M(\gamma)-\mu^\star\) shift). At the inner minimizer,

\[
\boxed{
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma^\star}{\det \Sigma}
\;-\;
\frac{1}{2}
\bigl(M(\gamma) - \mu^\star\bigr)^\top
\Gamma_M
\bigl(M(\gamma) - \mu^\star\bigr).
}
\]

**Derivation sketch.** Substitute
\(g(\gamma') = \text{const}
+ \tfrac{1}{2}(\gamma'-M)^\top\Sigma^{-1}(\gamma'-M)
- \tfrac{1}{2}(\gamma'-\mu^\star)^\top\Sigma^{\star-1}(\gamma'-\mu^\star)\)
and complete the square at \(\gamma' = \gamma'_\star(\gamma)\); the cross term
vanishes at the critical point, leaving the Schur complement
\((M-\mu^\star)^\top\Gamma_M(M-\mu^\star)\).

#### 5.1.1 Affine kernel mean — quadratic in \(\gamma\) alone

Assume the Normal-conjugate affine map from §9.1 (Block~1 + Block~2 GLS):

\[
M(\gamma) = \mathcal{L}\,\gamma + b(y),
\qquad
\mathcal{L} \in \mathbb R^{q\times q},
\qquad
b(y) \in \mathbb R^q,
\]

with refresh center \(\mu^\star = E[\gamma\mid y]\) (Hypothesis 2). Define the
**offset at \(\mu^\star\)** (data / prior only, no current \(\gamma\)):

\[
\boxed{
b_0(y)
:=
M(\mu^\star) - \mu^\star
=
\mathcal{L}\,\mu^\star + b(y) - \mu^\star,
}
\]

so \(M(\gamma)-\mu^\star = \mathcal{L}(\gamma-\mu^\star)+b_0(y)\).

Substitute into §5.1. Write \(\Gamma_M := (\Sigma + \Sigma^\star)^{-1}\) and expand
with \(u := \gamma - \mu^\star\):

\[
\bigl(M(\gamma)-\mu^\star\bigr)^\top \Gamma_M \bigl(M(\gamma)-\mu^\star\bigr)
=
\bigl(\mathcal{L} u + b_0\bigr)^\top \Gamma_M \bigl(\mathcal{L} u + b_0\bigr)
=
u^\top \mathcal{L}^\top \Gamma_M \mathcal{L}\, u
+
2\, b_0^\top \Gamma_M \mathcal{L}\, u
+
b_0^\top \Gamma_M b_0.
\]

**Re-center in \(\gamma\).** When \(\mathcal{L}^\top \Gamma_M \mathcal{L}\) is
invertible (e.g.\ \(\mathcal{L}\) full rank), define the **penalty centroid**
(data-only; not the refresh mean in general)

\[
\boxed{
\mu_\delta(y)
:=
\mu^\star
-
\bigl(\mathcal{L}^\top \Gamma_M \mathcal{L}\bigr)^{-1}
\mathcal{L}^\top \Gamma_M\, b_0(y),
}
\]

and the **\(\gamma\)-space precision**

\[
\boxed{
\Gamma_\gamma
:=
\mathcal{L}^\top \Gamma_M \mathcal{L}
=
\mathcal{L}^\top (\Sigma + \Sigma^\star)^{-1} \mathcal{L}.
}
\]

Then the shift satisfies \(\mathcal{L}(\gamma-\mu^\star)+b_0
= \mathcal{L}\bigl(\gamma-\mu_\delta\bigr)\) **iff** \(\mathcal{L}\) is
invertible; in general the Mahalanobis term **identically** re-centers as

\[
\boxed{
\bigl(M(\gamma)-\mu^\star\bigr)^\top \Gamma_M \bigl(M(\gamma)-\mu^\star\bigr)
=
\bigl(\gamma - \mu_\delta\bigr)^\top \Gamma_\gamma \bigl(\gamma - \mu_\delta\bigr).
}
\]

*Proof.* Expand the left-hand side as
\(u^\top \mathcal{L}^\top \Gamma_M \mathcal{L}\, u + 2 b_0^\top \Gamma_M \mathcal{L}\, u
+ b_0^\top \Gamma_M b_0\) with \(u=\gamma-\mu^\star\), and complete the square in
\(u\); the residual constant vanishes, leaving
\((\gamma-\mu_\delta)^\top \Gamma_\gamma(\gamma-\mu_\delta)\) with
\(\mu_\delta\) and \(\Gamma_\gamma\) as above. \(\square\)

**Full \(\delta(\gamma)\) in \(\gamma\)-centered form:**

\[
\boxed{
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma^\star}{\det \Sigma}
\;-\;
\frac{1}{2}
\bigl(\gamma - \mu_\delta\bigr)^\top \Gamma_\gamma \bigl(\gamma - \mu_\delta\bigr).
}
\]

| Symbol | Depends on \(\gamma\)? | Meaning |
|---|---|---|
| \(\mu_\delta(y)\) | no | centroid of inner penalty in \(\gamma\)-space |
| \(\Gamma_\gamma\) | no (Normal closure) | precision on \((\gamma-\mu_\delta)\) |
| \(\mu^\star\) | no | refresh mean \(E[\gamma\mid y]\); \(\mu_\delta=\mu^\star\) iff \(b_0=\mathbf 0\) |
| \(b_0(y)\) | no | \(M(\mu^\star)-\mu^\star\); drives \(\mu_\delta \neq \mu^\star\) |

**Zero offset \(b_0(y)=\mathbf 0\)** (equivalently \(M(\mu^\star)=\mu^\star\)):
\(\mu_\delta=\mu^\star\) and

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma^\star}{\det \Sigma}
\;-\;
\frac{1}{2}
(\gamma - \mu^\star)^\top
\mathcal{L}^\top (\Sigma + \Sigma^\star)^{-1} \mathcal{L}
(\gamma - \mu^\star).
\]

**Scalar intercept route** (§9.1.2, \(\mathcal{L}=\alpha_1\), \(V=\Sigma\),
\(\sigma_\star^2=\Sigma^\star\), \(\Gamma_M = 1/(V+\sigma_\star^2)\)):

\[
\Gamma_\gamma = \frac{\alpha_1^2}{V + \sigma_\star^2},
\qquad
\mu_\delta = \mu^\star - \frac{b_0(y)}{\alpha_1},
\]

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\sigma_\star^2}{V}
\;-\;
\frac{\alpha_1^2}{2\,(V + \sigma_\star^2)}
\bigl(\gamma - \mu_\delta\bigr)^2.
\]

**Convexity in \(\gamma\).** For fixed data, \(-\delta(\gamma)\) is concave quadratic
in \(\gamma\) when \(\Gamma_\gamma \succ 0\); the outer worst case on the
**penalty ellipsoid** \(C_d\) (§5.1.2) is on \(\partial C_d\) with a closed
\(\varepsilon_d\) (§6, §9.4).

#### 5.1.2 Penalty ellipsoid \(C_d\) (small set for coupling)

Fix a radius \(d > 0\). **Canonical definition:** \(C_d\) is a **level set of the
outer objective** \(\delta(\gamma)\) (the function plugged into the outer
minimum). With \(L_\delta(\gamma) := \delta(\mu^\star) - \delta(\gamma)\),

\[
C_d
\;:=\;
\bigl\{
\gamma :
\delta(\gamma) \ge \delta(\mu^\star) - \tfrac{d}{2}
\bigr\}
\;=\;
\bigl\{
\gamma :
L_\delta(\gamma) \le \tfrac{d}{2}
\bigr\},
\qquad
\varepsilon_d = \exp\bigl(\delta(\mu^\star) - \tfrac{d}{2}\bigr).
\]

See `MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` §7.4. No condition on \(b_0(y)\) enters
this definition: \(\mu^\star\) is always the refresh center and the **anchor** for
\(L_\delta\) (\(L_\delta(\mu^\star)=0\)).

**Gaussian closure (equivalent form).** From §5.1.1,
\(\delta(\gamma) = c_0 - \tfrac{1}{2}\|M(\gamma)-\mu^\star\|_{\Gamma_M}^2\) with
\(c_0 = \tfrac{1}{2}\log(\det\Sigma^\star/\det\Sigma)\) and
\(\|v\|_{\Gamma_M}^2 := v^\top\Gamma_M v\). Then

\[
L_\delta(\gamma)
=
\tfrac{1}{2}
\Bigl(
\|M(\gamma)-\mu^\star\|_{\Gamma_M}^2
-
\|M(\mu^\star)-\mu^\star\|_{\Gamma_M}^2
\Bigr),
\]

and \(C_d\) is the **shifted penalty ellipsoid**

\[
\boxed{
C_d
=
\bigl\{
\gamma :
\|M(\gamma)-\mu^\star\|_{\Gamma_M}^2
\le
\|M(\mu^\star)-\mu^\star\|_{\Gamma_M}^2 + d
\bigr\}.
}
\]

After completing the square (§5.1.1), the same set is
\((\gamma-\mu_\delta)^\top\Gamma_\gamma(\gamma-\mu_\delta)
\le d + (\mu^\star-\mu_\delta)^\top\Gamma_\gamma(\mu^\star-\mu_\delta)\);
\(\mu_\delta\) is **algebra only**, not a separate center for \(C_d\).

| Feature | Role |
|---|---|
| Anchor \(\mu^\star\) | Refresh mean; \(L_\delta(\mu^\star)=0\); threshold \(\delta(\mu^\star)-d/2\) |
| \(\|M(\gamma)-\mu^\star\|_{\Gamma_M}\) | Inner penalty in \(\gamma\) (§5.1.1) |
| \(\Gamma_\gamma\) | \(\mathcal{L}^\top(\Sigma+\Sigma^\star)^{-1}\mathcal{L}\); completed-square precision |
| Radius \(d\) | Certification knob (analogous to drift threshold in §8 of the rate note) |

**Scalar intercept** (\(\Gamma_\gamma = \alpha_1^2/(V+\sigma_\star^2)\),
\(M(\gamma)=\mu^\star+\alpha_1(\gamma-\mu^\star)+b_0(y)\)):

\[
C_d
=
\Bigl\{
\gamma :
\alpha_1^2(\gamma-\mu^\star)^2 + 2\alpha_1 b_0(y)(\gamma-\mu^\star)
\le d
\Bigr\}.
\]

Evaluate \(\mathcal{L}\), \(b_0(y)=M(\mu^\star)-\mu^\star\), \(\Sigma\), \(\Sigma^\star\)
(§5.1.1; \(\mathcal{L}\) from §9.1).

**On \(C_d\).** By construction,
\(\min_{\gamma\in C_d}\delta(\gamma)=\delta(\mu^\star)-d/2\) and
\(\varepsilon_d=\exp(\delta(\mu^\star)-d/2)\) (§6). On the boundary,
\(\delta(\gamma)=\delta(\mu^\star)-d/2\).

**Drift intersection (optional).** Foster–Lyapunov drift in
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §8 may use a sublevel set centered at the
**prior** \(\mu_0\). For coupling certified on states the chain actually visits,
take \(C_d \cap \{\gamma : V_{\mathrm{drift}}(\gamma) \le d_{\mathrm{drift}}\}\)
(or choose \(d\) so the penalty ellipsoid lies inside the drift set).

**Deviance alternative (GLM / non-quadratic \(\pi_\gamma\)).** When a global
quadratic mis-specifies the posterior shape, use a **posterior deviance** contour
\(C_d^{(\pi)} = \{\gamma : D_\pi(\gamma) \le d\}\) with
\(D_\pi(\gamma) = 2[\log\pi_\gamma(\mu^\star)-\log\pi_\gamma(\gamma)]\) (scalar:
\(|r_{\mathrm{dev}}(\gamma)| \le \sqrt{d}\)). Locally \(D_\pi \approx
(\gamma-\mu^\star)^\top \Gamma_{\mathrm{post}}(\gamma-\mu^\star)\); the penalty
ellipsoid is the **small-\(d\)** limit. See
`MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` §7.4–§7.5; outer \(\varepsilon_d\) is generally
**numeric** on deviance contours.

### 5.2 Scalar specialization

\[
\boxed{
\delta(\gamma)
=
\frac{1}{2}\log \frac{\sigma_\star^2}{V}
\;-\;
\frac{\bigl(M(\gamma) - \mu^\star\bigr)^2}{2\,(V + \sigma_\star^2)}.
}
\]

### 5.3 When \(M(\gamma) = \mu^\star\)

The Mahalanobis term vanishes:

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma^\star}{\det \Sigma},
\qquad
\exp\bigl(\delta(\gamma)\bigr)
=
\sqrt{\det\bigl(\Sigma^\star \Sigma^{-1}\bigr)}.
\]

Scalar: \(\exp(\delta) = \sigma_\star / \sqrt{V}\).

This occurs at \(\gamma\) with \(M(\gamma) = \mu^\star\), **not** necessarily at
\(\gamma = \mu^\star\) unless \(M(\mu^\star) = \mu^\star\).

### 5.4 Ratio form at the inner minimizer

\[
\exp\bigl(\delta(\gamma)\bigr)
=
\frac{
q\bigl(\gamma'_\star(\gamma) \mid \gamma\bigr)
}{
q_Q\bigl(\gamma'_\star(\gamma)\bigr)
}
=
\sqrt{\frac{\det \Sigma^\star}{\det \Sigma}}
\;
\exp\!\Bigl(
-\tfrac{1}{2}
\bigl\|M(\gamma) - \mu^\star\bigr\|_{\Gamma_M}^2
\Bigr)
=
\exp\!\Bigl(
-\tfrac{1}{2}
\bigl(\gamma - \mu_\delta\bigr)^\top \Gamma_\gamma \bigl(\gamma - \mu_\delta\bigr)
\Bigr)
\]

when \(M(\gamma)=\mathcal{L}\gamma+b(y)\) (§5.1.1).

---

## 6. Outer minimum (certified \(\varepsilon_d\))

The **certified minorization constant** for Hypotheses 1–2 is

\[
\boxed{
\gamma_\star \in \arg\min_{\gamma \in C_d} \delta(\gamma),
\qquad
\varepsilon_d
=
\exp\bigl(\delta(\gamma_\star)\bigr)
=
\frac{
q\bigl(\gamma'_\star(\gamma_\star) \mid \gamma_\star\bigr)
}{
q_Q\bigl(\gamma'_\star(\gamma_\star)\bigr)
}.
}
\]

Binding triple: \((\gamma_\star,\, \gamma'_\star(\gamma_\star),\, \delta(\gamma_\star))\).

**\(\delta\)-level set \(C_d\) (§5.1.2).** By definition of \(C_d\),

\[
\boxed{
\min_{\gamma \in C_d} \delta(\gamma)
=
\delta(\mu^\star) - \tfrac{d}{2},
\qquad
\gamma_\star \in \partial C_d
\ \text{(any boundary point)}.
}
\]

\[
\boxed{
\varepsilon_d
=
\exp\Bigl(\delta(\mu^\star) - \tfrac{d}{2}\Bigr).
}
\]

**Gaussian closed form.** With
\(\delta(\mu^\star)
= \tfrac{1}{2}\log(\det\Sigma^\star/\det\Sigma)
- \tfrac{1}{2}\|M(\mu^\star)-\mu^\star\|_{\Gamma_M}^2\),

\[
\boxed{
\varepsilon_d
=
\sqrt{\det\bigl(\Sigma^\star \Sigma^{-1}\bigr)}
\;
\exp\!\Bigl(
-\tfrac{1}{2}\|M(\mu^\star)-\mu^\star\|_{\Gamma_M}^2
-\tfrac{d}{2}
\Bigr).
}
\]

When \(M(\mu^\star)=\mu^\star\) this reduces to
\(\sqrt{\det(\Sigma^\star\Sigma^{-1})}\,e^{-d/2}\).

*Proof.* \(L_\delta(\gamma)\le d/2\) on \(C_d\) with equality on \(\partial C_d\);
\(\delta=\delta(\mu^\star)-L_\delta\). Substitute §5.1.1. \(\square\)

For all \(\gamma \in C_d\), \(\gamma' \in \mathbb R^q\):

\[
q(\gamma', \gamma) \;\ge\;
\varepsilon_d\, q_Q(\gamma')
\qquad\text{with}\qquad
q_Q = \phi_q(\cdot;\,\mu^\star,\,\Sigma^\star).
\]

---

## 7. Plug-in summary (Hypotheses 1–2)

| Quantity | Formula |
|---|---|
| \(Q\) | \(\phi_q(\gamma';\,\mu^\star,\,\Sigma^\star)\) |
| \(q(\cdot\mid\gamma)\) | \(\phi_q(\cdot;\,M(\gamma),\,\Sigma)\) |
| Convexity | \(\Sigma \succ \Sigma^\star\) |
| Inner argmin | \(\gamma'_\star(\gamma)\) from §4 |
| Inner log gap | \(\delta(\gamma)\) from §5 (§5.1.1: \((\gamma-\mu_\delta)^\top\Gamma_\gamma(\gamma-\mu_\delta)\) term) |
| Small set \(C_d\) | Penalty ellipsoid §5.1.2: \((\gamma-\mu_\delta)^\top\Gamma_\gamma(\gamma-\mu_\delta)\le d\) |
| Certified \(\varepsilon_d\) | \(\sqrt{\det(\Sigma^\star\Sigma^{-1})}\,e^{-d/2}\) on \(C_d\) (§6) |

**Scalar checklist** with
\(\sigma_\star^2 = \Sigma^\star\), \(V = \Sigma\):

\[
\gamma'_\star(\gamma)
=
\frac{V\mu^\star - \sigma_\star^2 M(\gamma)}{V - \sigma_\star^2},
\]

\[
\delta(\gamma)
=
\tfrac{1}{2}\log(\sigma_\star^2/V)
-
\tfrac{\Gamma_\gamma}{2}\bigl(\gamma-\mu_\delta\bigr)^2,
\qquad
\varepsilon_d = (\sigma_\star/\sqrt{V})\,e^{-d/2}
\ \text{on}\ C_d.
\]

---

## 8. Remarks

1. **Why \(\mu^\star\) and not \(M(\gamma)\)?** At fixed \(\gamma\), the best
   center for \(Q\) would be \(M(\gamma)\) (zeros the penalty at that fiber), but
   \(Q\) must be **fixed** for Rosenthal coupling. \(\mu^\star = E[\gamma\mid y]\)
   centers the refresh on the **posterior** location of the population parameter.

2. **Why \(\Sigma^\star\) and not \(\Sigma\)?** \(\Sigma^\star\) is the conditional
   spread of \(\gamma'\mid\beta\); \(\Sigma\) adds \(\mathrm{Cov}(m(\beta)\mid\gamma,y)\).
   Using \(\Sigma^\star\) for \(Q\) ensures \(\Lambda_q \prec \Lambda_Q\) and an
   **interior** inner minimizer \(\gamma'_\star(\gamma)\).

3. **Non-Gaussian likelihoods.** Block 2 remains
   \(\phi_q(\gamma';m(\beta),\Sigma^\star)\); \(q\) is not Gaussian. Keep
   Hypotheses 1–2 for \(Q\); verify \(\Lambda_q \prec \Sigma^{\star-1}\)
   numerically and solve \(\nabla_{\gamma'} g = 0\) (rate note §10.7).

4. **Coupling.** Plug \(\varepsilon_d\) from §6 into Rosenthal’s bound
   (`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §10.1; `LOGIT_SINGLE_GROUP_SAFE_REGION.md`
   §4.6).

6. **Penalty ellipsoid.** \(C_d\) uses the **same** \((\mu_\delta,\Gamma_\gamma)\)
   as the \(\gamma\)-quadratic in \(\delta\), so the outer min is closed and does
   not require a separate shape matrix \(G\) or generalized eigenvalues (§9.4.1
   for misaligned alternatives).

---

## 9. Small set \(C_d\): ellipsoid and outer minimum on the boundary

**Notation (package models).** Symbols follow `inst/notation.md` (`lmerb`
hierarchy): \(\beta_j \in \mathbb R^P\), level-2 design \(\mathcal{W}_j \in
\mathbb R^{P\times q}\) (`Wlist` in code), stage-2 covariance \(\Psi \succ 0\),
population \(\gamma \in \mathbb R^q\), prior \(\gamma \sim N(\mu_0,
\Lambda_\gamma^{-1})\) with precision \(\Lambda_\gamma\) (scalar route:
\(\Lambda_\gamma = \lambda_\gamma I\), \(\Sigma_\gamma = \lambda_\gamma^{-1}\)).
Block~2 uses \(\Sigma^\star = \Sigma_{\mathrm{upd}}\) and mean \(m(\beta)\) as in
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §8.1 / §10.1 (scalar:
\(m(\beta)=w\bar\beta+(1-w)\mu_0\), \(w=J\lambda_b/(J\lambda_b+\lambda_\gamma)\),
\(\Sigma^\star=\sigma_\gamma^2\)).
**Do not confuse** \(\mathcal{W}_j\) here (level-2 design, `notation.md`) with
\(\mathcal{W}_j = Z_j^\top W_j Z_j\) in `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md`
(data-weight matrix on RE space; there \(W_j=\mathrm{diag}(w_{j,i})\)).
**Disambiguation:** the matrix \(\mathcal{L}\in\mathbb R^{q\times q}\) below in
\(M(\gamma)=\mathcal{L}\gamma+b(y)\) is **not** the rate safe set \(A\subseteq(\gamma,\beta)\)-space nor the
contraction matrix \(A(\gamma,\beta)\) of `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md`.

The outer problem §6 minimizes \(\delta(\gamma)\) over \(\gamma \in C_d\). With
affine \(M(\gamma)\) (§9.1),

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma^\star}{\det \Sigma}
\;-\;
\frac{1}{2}
\bigl(\gamma - \mu_\delta\bigr)^\top \Gamma_\gamma \bigl(\gamma - \mu_\delta\bigr)
\]

(§5.1.1). The **default** small set \(C_d\) uses **this same quadratic**
(§5.1.2), so the outer minimum is closed (§6).

### 9.1 Kernel mean \(M(\gamma)\) — multivariate Normal conjugate

Throughout §9, \(\gamma, \gamma', \mu_0, \mu^\star, M(\gamma) \in \mathbb R^q\),
\(\Sigma, \Sigma^\star \in \mathbb R^{q\times q}\). Block~2 is
\(\gamma' \mid \beta \sim N(m(\beta), \Sigma^\star)\) with \(m(\beta) \in \mathbb R^q\)
(§1.1).

#### 9.1.1 `lmerb` hierarchy (`notation.md`)

**Stage 2 (prior on group coefficients).**

\[
\beta_j \mid \gamma \;\sim\; N\bigl(\mathcal{W}_j \gamma,\, \Psi\bigr),
\qquad
j = 1,\ldots,J,
\]

with \(\mathcal{W}_j \in \mathbb R^{P\times q}\) the level-2 design (`Wlist`;
block-diagonal across coefficient dimensions) and \(\Psi \in \mathbb R^{P\times P}\)
the stage-2 covariance (scalar RE route: \(P=1\), \(\mathcal{W}_j=1\),
\(\Psi=\lambda_b^{-1}\)).

**Stage 1 (Normal likelihood).**

\[
y_j \mid \beta_j \;\sim\; N\bigl(D_j \beta_j,\, \sigma_j^2 I_{n_j}\bigr),
\]

\(D_j \in \mathbb R^{n_j\times P}\) the likelihood design (`Dmat` restricted to
group \(j\)). Conjugate Block~1 gives \(\beta_j \mid \gamma, y_j \sim
N(m_j(\gamma), V_j)\) with

\[
V_j
=
\Bigl(\tfrac{1}{\sigma_j^2} D_j^\top D_j + \Psi^{-1}\Bigr)^{-1},
\qquad
m_j(\gamma)
=
V_j\Bigl(
\tfrac{1}{\sigma_j^2} D_j^\top y_j + \Psi^{-1}\mathcal{W}_j \gamma
\Bigr),
\]

hence **affine** group means

\[
\boxed{
m_j(\gamma)
=
\mathcal{L}_j\,\gamma + b_j(y),
\qquad
\mathcal{L}_j
:=
V_j\,\Psi^{-1}\mathcal{W}_j
\;\in\; \mathbb R^{P\times q},
\qquad
b_j(y)
:=
V_j\,\tfrac{1}{\sigma_j^2} D_j^\top y_j
\;\in\; \mathbb R^P.
}
\]

**Block~2 mean given \(\beta\)** (multivariate extension of §8.1; GLS across groups
with prior \(\gamma \sim N(\mu_0, \Lambda_\gamma^{-1})\)):

\[
\boxed{
P_\gamma
:=
\Lambda_\gamma
+
\sum_{j=1}^{J} \mathcal{W}_j^\top \Psi^{-1} \mathcal{W}_j,
\qquad
m(\beta)
=
P_\gamma^{-1}
\Bigl(
\Lambda_\gamma \mu_0
+
\sum_{j=1}^{J} \mathcal{W}_j^\top \Psi^{-1} \beta_j
\Bigr),
\qquad
\Sigma^\star
=
P_\gamma^{-1}
=
\Sigma_{\mathrm{upd}}.
}
\]

(Scalar intercept route §0.1 / §8.1: \(P=1\), \(\mathcal{W}_j=1\),
\(\Psi=\lambda_b^{-1}\), \(\Lambda_\gamma=\lambda_\gamma\), \(\sigma_j^2=\phi\),
\(D_j=\mathbf 1_{n_j}\) recovers \(P_\gamma = J\lambda_b+\lambda_\gamma\),
\(m(\beta)=w\bar\beta+(1-w)\mu_0\), \(\Sigma^\star=\sigma_\gamma^2\).)

**Marginal kernel mean** \(M(\gamma) := E[\gamma' \mid \gamma] = E[m(\beta)\mid\gamma,y]\):

\[
\boxed{
M(\gamma)
=
\mathcal{L}\,\gamma + b(y),
\qquad
\mathcal{L}
:=
P_\gamma^{-1}
\sum_{j=1}^{J} \mathcal{W}_j^\top \Psi^{-1} \mathcal{L}_j
\;\in\; \mathbb R^{q\times q},
\qquad
b(y)
:=
P_\gamma^{-1}
\Bigl(
\sum_{j=1}^{J} \mathcal{W}_j^\top \Psi^{-1} b_j(y)
+
\Lambda_\gamma \mu_0
\Bigr)
\;\in\; \mathbb R^q.
}
\]

**Shift from refresh center \(\mu^\star\).**

\[
\boxed{
M(\gamma) - \mu^\star
=
\mathcal{L}\bigl(\gamma - \mu^\star\bigr) + b_0(y),
\qquad
b_0(y)
:=
\mathcal{L}\,\mu^\star + b(y) - \mu^\star.
}
\]

| Object | Role of \(\gamma\) | Role of \(y\) | Role of \(\mu_0\) |
|---|---|---|---|
| \(\mathcal{L}_j, b_j(y)\) | via \(\mathcal{L}_j\gamma\) | via \(D_j^\top y_j\) | — |
| \(\mathcal{L}, b(y)\) | pooled \(\mathcal{L}\gamma\) | pooled \(b_j\) | \(\Lambda_\gamma\mu_0\) in \(b(y)\) |
| \(\mu^\star\) | center of \(C_d\) | \(E[\gamma\mid y]\) | not prior mean |
| \(b_0(y)\) | — | data + \(\mu_0,\mu^\star\) | offset at \(\gamma=\mu^\star\) |

\(b_0(y)=\mathbf 0\) iff \(M(\mu^\star)=\mu^\star\). In general \(b_0 \neq \mathbf 0\).

#### 9.1.2 Scalar intercept RE (`RATE_Ac_BOUNDS` §0.1)

Same as §9.1.1 with \(P=1\), \(\mathcal{W}_j=1\), \(\Psi=\lambda_b^{-1}\),
shared \(\phi\), and group **likelihood precision**
\(\mathscr{W}_j := n_j/\phi\) (rate note §1; **not** the level-2 design
\(\mathcal{W}_j\)). Then

\[
\mathcal{L}_j = \kappa_j := \frac{\lambda_b}{\lambda_b + \mathscr{W}_j},
\qquad
b_j(y) = \frac{\mathscr{W}_j\,\bar y_j}{\lambda_b + \mathscr{W}_j},
\]

\[
\mathcal{L} = \alpha_1 := w\,\bar\kappa,
\quad
\bar\kappa := \frac{1}{J}\sum_{j=1}^{J}\kappa_j,
\quad
w := \frac{J\lambda_b}{J\lambda_b + \lambda_\gamma},
\]

\[
b(y) = w\,\kappa_0(y) + (1-w)\mu_0,
\qquad
\kappa_0(y) := \frac{1}{J}\sum_{j=1}^{J}\frac{\mathscr{W}_j\,\bar y_j}{\lambda_b + \mathscr{W}_j},
\]

recovering \(M(\gamma) = \alpha_0(y) + \alpha_1 \gamma\) with \(\alpha_0=b(y)\),
\(\Sigma^\star=\sigma_\gamma^2=(J\lambda_b+\lambda_\gamma)^{-1}\), and
\(\Sigma = \sigma_\gamma^2 + w^2\,\mathrm{Var}(\bar\beta\mid\gamma,y)\) (§1.2).

### 9.2 Penalty ellipsoid \(C_d\) (§5.1.2)

After §9.1 and §5.1.1, \(\mu_\delta(y)\) and \(\Gamma_\gamma\) are fixed by data.
The **coupling small set** is

\[
\boxed{
C_d
=
\bigl\{
\gamma :
\bigl(\gamma - \mu_\delta\bigr)^\top \Gamma_\gamma \bigl(\gamma - \mu_\delta\bigr)
\le d
\bigr\}.
}
\]

This is the sublevel set of the **penalty term** in \(\delta(\gamma)\); it is
**not** required to be centered at the refresh mean \(\mu^\star\) unless
\(b_0(y)=\mathbf 0\). Choosing \(d\) sets how much log-gap slack is certified on
\(C_d\) (§6: \(\varepsilon_d = \sqrt{\det(\Sigma^\star\Sigma^{-1})}\,e^{-d/2}\)).

**Scalar intercept:** \(C_d = [\mu_\delta - \sqrt{d/\Gamma_\gamma},\,
\mu_\delta + \sqrt{d/\Gamma_\gamma}]\).

### 9.3 \(\delta(\gamma)\) in \(\gamma\)-centered form

With affine \(M(\gamma)=\mathcal{L}\gamma+b(y)\) from §9.1, the substitution and
re-centering are in **§5.1.1**:

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma^\star}{\det \Sigma}
\;-\;
\frac{1}{2}
\bigl(\gamma - \mu_\delta\bigr)^\top \Gamma_\gamma \bigl(\gamma - \mu_\delta\bigr),
\]

\[
\Gamma_\gamma
=
\mathcal{L}^\top (\Sigma + \Sigma^\star)^{-1} \mathcal{L},
\qquad
\mu_\delta
=
\mu^\star
-
\bigl(\mathcal{L}^\top (\Sigma + \Sigma^\star)^{-1} \mathcal{L}\bigr)^{-1}
\mathcal{L}^\top (\Sigma + \Sigma^\star)^{-1}\, b_0(y).
\]

**Equivalent shift form** (before re-centering): with
\(\Gamma_M := (\Sigma + \Sigma^\star)^{-1}\) and \(u=\gamma-\mu^\star\),

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma^\star}{\det \Sigma}
\;-\;
\frac{1}{2}
\bigl\|\mathcal{L} u + b_0(y)\bigr\|_{\Gamma_M}^{2}.
\]

**Scalar intercept** (\(q=1\), \(\mathcal{L}=\alpha_1\), \(V=\Sigma\),
\(\sigma_\star^2=\Sigma^\star\), \(\mu_\delta = \mu^\star - b_0/\alpha_1\)):

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\sigma_\star^2}{V}
\;-\;
\frac{\alpha_1^2}{2\,(V + \sigma_\star^2)}
\bigl(\gamma - \mu_\delta\bigr)^2.
\]

### 9.4 Outer minimum and certified \(\varepsilon_d\)

When \(C_d\) is the **penalty ellipsoid** (§9.2 = §5.1.2), the outer problem
from §6 is **closed**:

\[
\boxed{
\gamma_\star \in \partial C_d
\ \text{(any boundary point)},
\qquad
\delta(\gamma_\star)
=
\frac{1}{2}\log \frac{\det \Sigma^\star}{\det \Sigma}
\;-\;
\frac{d}{2},
\qquad
\varepsilon_d
=
\sqrt{\det\bigl(\Sigma^\star \Sigma^{-1}\bigr)}
\;
\exp\!\Bigl(-\frac{d}{2}\Bigr).
}
\]

*Proof.* On \(C_d\), \((\gamma-\mu_\delta)^\top\Gamma_\gamma(\gamma-\mu_\delta)\le d\);
\(\delta\) decreases in that quadratic, so the minimum over \(C_d\) is at
\(\partial C_d\) with value \(c_0-d/2\). \(\square\)

Report the binding triple
\((\gamma_\star,\, \gamma'_\star(\gamma_\star),\, \varepsilon_d)\) with
\(\gamma'_\star(\gamma)\) from §4.

#### 9.4.1 Misaligned outer set (optional)

If \(C_d\) is instead taken as \(\{\gamma : (\gamma-\mu^\star)^\top G(\gamma-\mu^\star)
\le d\}\) with **independent** \(G \succ 0\) (e.g.\ \(\Lambda_\gamma\) or a drift
Lyapunov shape), the penalty quadratic is **not** constant on \(C_d\). Then
\(\min_{C_d}\delta(\gamma)\) requires maximizing
\((\gamma-\mu_\delta)^\top\Gamma_\gamma(\gamma-\mu_\delta)\) over that ellipsoid
(generalized eigenvalue / Lagrange methods as in earlier drafts). The **default**
coupling route uses the **matched** penalty ellipsoid §9.2 so §6 applies without
extra optimization.

---

## References

- `inst/notation.md` — \(\beta_j\), \(\mathcal{W}_j\), \(\Psi\), \(\gamma\), Block~1
  conditionals.
- `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` — §0.1 scalar RE, §8.1 Block~2, §10.1 coupling,
  §10.3 two-stage min, §10.5 Normal closure, §10.7 multivariate precision form.
- `MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` — symmetric GLM population
  posteriors (logit/probit fixtures); extension when \(q\) is not Gaussian.
- `LOGIT_SINGLE_GROUP_SAFE_REGION.md` — merge construction and Theorem 12.
