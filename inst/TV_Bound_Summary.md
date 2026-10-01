# TV Distance Bounds for Two-Block Gibbs Samplers: A Summary of Key Results

## Overview

This document summarizes a progression of results bounding the total
variation (TV) distance between the $k$-step distribution of a two-block
Gibbs sampler and its target, in the context of hierarchical generalized
linear models (GLMMs). The results progress from exact Gaussian results
(Nygren 2020) through eigenvalue-based bounds for hierarchical Gaussian
models (Foundational Lemmas) to a general theorem for a broad class of
symmetric hierarchical GLMMs (Main Theorem document).

---

## 1. Nygren (2020): TV bounds for Gaussian two-block Gibbs samplers

**Reference.** Nygren, K. (2020). *On the total variation distance between
multivariate normal densities with applications to two-block Gibbs samplers.*

**Setting.** A two-block Gibbs sampler applied to a multivariate normal
target $N(\cdot|\mu,\Sigma)$ of dimension $n$, with blocks $X_1$
(dimension $q$) and $X_2$ (dimension $n-q$). The joint transition kernel
alternates between drawing $X_2|X_1$ and $X_1|X_2$, both Gaussian. The
precision matrix is partitioned as:
$$P = \Sigma^{-1} = \begin{pmatrix}P_{11}&P_{12}\\P_{21}&P_{22}\end{pmatrix}$$

**Key matrix.** Define $A := P_{11}^{-1/2}P_{12}P_{22}^{-1}P_{21}P_{11}^{-1/2}$,
the standardized between-block precision product. Its eigenvalues
$0\le a_1^2\le\cdots\le a_q^2 < 1$ are the squared canonical correlations
between the two blocks. The maximal eigenvalue $\lambda^* := \max_i a_i^2$
governs the geometric convergence rate.

**Theorem 3 (Nygren 2020).** *After $l$ iterations of the two-block Gibbs
sampler started at $x_1^{(0)}$:*
$$\bigl\|Q_{X_1}^{(l)}(x_1^{(0)},\cdot) - N(\cdot|\mu_1,\Sigma_{11})\bigr\|_{TV}
\;\le\; \sum_{i=1}^q d_i^{(l)} \;+\;
\mathrm{erf}_1\!\left(\frac{\sqrt{0.5\,(\lambda^*)^l\,(x_1^{(0)}-\mu_1)^T\Sigma_{11}^{-1}(x_1^{(0)}-\mu_1)}}{\sqrt{2}}\right)$$

*where $d_i^{(l)}$ are terms involving the $n$-dimensional error function
$\mathrm{erf}_n$ applied to the eigenvalue ratios:*
$$r_i^{(l)} = \frac{1-a_{i-1}^{2l}}{1-a_i^{2l}}, \qquad
d_i^{(l)} = \mathrm{erf}_{q+1-i}\!\left(\sqrt{\frac{r_i^{(l)}(q+1-i)\ln r_i^{(l)}}{2(r_i^{(l)}-1)}}\right)
- \mathrm{erf}_{q+1-i}\!\left(\sqrt{\frac{(q+1-i)\ln r_i^{(l)}}{2(r_i^{(l)}-1)}}\right)$$

**Corollary 1 (Geometric ergodicity).** *The bound decays geometrically in
$l$ at rate $(\lambda^*)^l$, establishing geometric ergodicity:*
$$\bigl\|Q_{X_1}^{(l)}(x_1^{(0)},\cdot)-N(\cdot|\mu_1,\Sigma_{11})\bigr\|_{TV}
\;\le\; C(\delta_2)\cdot(\lambda^*)^{l/2}$$

*for an explicit constant $C(\delta_2)$ depending on the starting point.*

**Key structural results** underlying Theorem 3:

*Lemma 1 (same covariance).* For two normal densities sharing covariance $\Sigma$:
$$\|N(\cdot|\mu_1,\Sigma) - N(\cdot|\mu_2,\Sigma)\|_{TV}
= \mathrm{erf}_1\!\left(\frac{\sqrt{0.5(\mu_2-\mu_1)^T\Sigma^{-1}(\mu_2-\mu_1)}}{\sqrt{2}}\right)$$

*Lemma 2 (same mean, ordered precision).* For $\Sigma_1^{-1}\succeq\Sigma_2^{-1}$,
with eigenvalues $1\le k_1\le\cdots\le k_n$ of $\Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2}$
and ratios $r_i=k_i/k_{i-1}$:
$$\|N(\cdot|\mu,\Sigma_1)-N(\cdot|\mu,\Sigma_2)\|_{TV} \;\le\; \sum_{i=1}^n d_i$$

**Implication for hierarchical models.** At fixed starting point
$x_1^{(0)}=\mu_1$ (i.e. starting at stationarity mean):
$$\bigl\|Q_{X_1}^{(l)}(\mu_1,\cdot)-N(\cdot|\mu_1,\Sigma_{11})\bigr\|_{TV}
\;\le\; \sum_{i=1}^q d_i^{(l)}$$

This is the **floor bound** — the TV distance between the $l$-step chain
and its Gaussian target, expressed entirely in terms of the eigenvalues
$a_i^2$ of $A$. The bound is tight and achieves the geometric rate
$(\lambda^*)^l$ as $l\to\infty$.

---

## 2. Foundational Lemmas: eigenvalue monotonicity of TV distances between Gaussians

**Document.** *Foundational Lemmas and a General Theorem on Total Variation
Distance Between Multivariate Normal Densities.*

**Setting.** The document establishes that the TV distance between two
multivariate normal densities depends only on the sorted eigenvalues of their
precision-ratio matrix, and that TV distance is monotone in those eigenvalues.
This provides the key structural result linking the floor model's eigenvalues
to a bound on the real model's TV distance.

**Lemma 1 (Reduction to diagonal form; TV depends only on eigenvalues).**
*Let $\Sigma_1,\Sigma_2$ be positive definite $n\times n$ matrices, and let
$M:=\Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2}=C\Lambda C^\top$ be an
eigendecomposition. Then:*

*(a)* The matrices $\Sigma_3:=\Lambda^{-1}$, $\Sigma_4:=I_n$ satisfy
$\Sigma_4^{1/2}\Sigma_3^{-1}\Sigma_4^{1/2}=\Lambda$ — diagonal with the
same eigenvalues as $M$.

*(b)*
$$\bigl\|N(0,\Sigma_1)-N(0,\Sigma_2)\bigr\|_{TV}
\;=\;\bigl\|N(0,\Sigma_3)-N(0,\Sigma_4)\bigr\|_{TV}
\;=\;\bigl\|N(0,\Lambda^{-1})-N(0,I)\bigr\|_{TV}$$

*Proof.* The invertible linear map $T(x):=C^\top\Sigma_2^{-1/2}x$ transforms
$N(0,\Sigma_2)\mapsto N(0,I)$ and $N(0,\Sigma_1)\mapsto N(0,\Lambda^{-1})$.
TV distance is invariant under any bijective measurable transformation.
$\blacksquare$

*Remark.* The TV distance between two mean-zero normals depends on $(\Sigma_1,\Sigma_2)$
**only through the eigenvalues** of $\Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2}$ —
the eigenbasis $C$ and the overall scale of $\Sigma_2$ are TV-irrelevant.

**Lemma 2 (TV is non-decreasing in a single eigenvalue).**
*Let $M=\mathrm{diag}(k_1,\ldots,k_n)$ and
$M'=\mathrm{diag}(k_1,\ldots,k_{j-1},k_j',k_{j+1},\ldots,k_n)$ be diagonal,
agreeing in every coordinate except $j$, where $k_j'\ge k_j$. Then:*
$$\bigl\|N(0,\Sigma_1)-N(0,\Sigma_2)\bigr\|_{TV}
\;\le\;\bigl\|N(0,\Sigma_3)-N(0,\Sigma_4)\bigr\|_{TV}$$
*strictly except in a measure-zero degenerate case. No sign or magnitude
restriction on the shared coordinates $k_i$, $i\ne j$, is required.*

*Proof.* By Lemma 1, it suffices to show $\phi(k_j)\le\phi(k_j')$ where
$\phi(\kappa):=\|N(0,\Lambda_\kappa^{-1})-N(0,I)\|_{TV}$ with $\Lambda_\kappa$
equal to $M$ except coordinate $j$ set to $\kappa$. The TV-optimal crossing
region at parameter $\kappa>1$ is $A_\kappa=\{x:x_j^2\le\tau(x_{-j},\kappa)\}$
for an explicit $\tau$. By the envelope theorem (differentiating while holding
$A_\kappa$ fixed) and the truncated-normal second-moment identity:
$$\phi'(\kappa) = \int h(x_{-j})\,\sqrt{\tau(x_{-j},\kappa)/\kappa}\;
\varphi\bigl(\sqrt{\kappa\,\tau(x_{-j},\kappa)}\bigr)\,
\mathbf{1}(\tau\ge0)\,dx_{-j}\;\ge\;0$$
with equality only when $\tau=0$ a.e. — the measure-zero degenerate case.
$\blacksquare$

*Remark.* This is a statement about the **true** TV distance, not the
triangle-inequality upper bound of Nygren (2020) Lemma 2. The two can move
in opposite directions under a single-eigenvalue increase.

**Lemma 3 (TV is non-decreasing under full coordinatewise eigenvalue
domination).**
*Let $M=\mathrm{diag}(k_1,\ldots,k_n)$ and
$M'=\mathrm{diag}(k_1',\ldots,k_n')$ with $k_i\le k_i'$ for every $i$.
Then:*
$$\bigl\|N(0,\Sigma_1)-N(0,\Sigma_2)\bigr\|_{TV}
\;\le\;\bigl\|N(0,\Sigma_3)-N(0,\Sigma_4)\bigr\|_{TV}$$

*Proof.* Define a chain $\Lambda^{(0)}=M,\Lambda^{(1)},\ldots,\Lambda^{(n)}=M'$
updating one coordinate at a time, each step replacing $k_m$ with $k_m'\ge k_m$.
Apply Lemma 2 at each step and chain by transitivity. $\blacksquare$

**Theorem 1 (TV monotonicity for arbitrary, non-diagonal covariance
matrices).**
*Let $N(\cdot|\mu_1,\Sigma_1)$, $N(\cdot|\mu_1,\Sigma_2)$,
$N(\cdot|\mu_2,\Sigma_3)$, $N(\cdot|\mu_2,\Sigma_4)$ be multivariate normals
(pairs may have different means). Let $k_1\le\cdots\le k_n$ be the sorted
eigenvalues of $\Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2}$ and
$k_1'\le\cdots\le k_n'$ be the sorted eigenvalues of
$\Sigma_4^{1/2}\Sigma_3^{-1}\Sigma_4^{1/2}$, with $k_i\le k_i'$ for every $i$.
Then:*
$$\bigl\|N(\cdot|\mu_1,\Sigma_1)-N(\cdot|\mu_1,\Sigma_2)\bigr\|_{TV}
\;\le\;\bigl\|N(\cdot|\mu_2,\Sigma_3)-N(\cdot|\mu_2,\Sigma_4)\bigr\|_{TV}$$

*Proof.* Lemma 1 reduces each pair to its diagonal form (against the identity
reference), and Lemma 3 then gives the result from the sorted eigenvalue
domination. The mean vectors do not affect the TV distance within each pair
(same-mean comparison). $\blacksquare$

**Key implication for hierarchical models with Gaussian data.** In the
two-block Gibbs sampler for a Gaussian hierarchical model, the $k$-step chain
$Q^{(k)}$ and its target $\pi$ are both Gaussian (Nygren 2020, Remark 7).
The precision-ratio matrix for the pair $(Q^{(k)},\pi)$ has eigenvalues
$1/(1-a_i^{2k})$ where $a_i^2$ are the eigenvalues of
$A=P_{11}^{-1/2}P_{12}P_{22}^{-1}P_{21}P_{11}^{-1/2}$.

Now consider two models: a **real** model (ridge-penalized, with likelihood
precision) and a **floor** model (ridge-only, Gaussian, no likelihood). The
floor model has **larger** $P_{22}$ (lower conditional precision for $\beta$)
and hence **smaller** $a_i^{0,2}$ — the floor eigenvalues dominate the real
eigenvalues: $a_i^2\ge a_i^{0,2}$ for every $i$. Consequently the precision
ratios $1/(1-a_i^{2k})\ge 1/(1-a_i^{0,2k})$ — the real model's pair has
**larger** precision-ratio eigenvalues than the floor model's pair. By
Theorem 1:

$$\boxed{\bigl\|Q_0^{(k)}(\cdot;\gamma_0)-\pi_0\bigr\|_{TV}
\;\le\;
\bigl\|Q_1^{(k)}(\cdot;\gamma_0)-\pi_1\bigr\|_{TV}}$$

That is, the **floor sampler converges more slowly** than the real sampler —
the floor model's TV distance bounds the real model's from above. The floor
model's TV distance is computable analytically from Nygren (2020) Theorem 3.

**Gradient and Hessian structure.** The document also derives the exact
curvature structure of the marginal posterior and sampler kernel for GLMM
posteriors, showing both share the form
$c(\gamma)=P_{11}-P_{12}V(\gamma)P_{21}$ — differing only in whether
$V(\gamma)=\mathrm{Cov}(\beta|\gamma,y)$ (posterior, single-point conditioning)
or $U(\gamma,\gamma')=\mathrm{Cov}(\beta|\gamma,\gamma',y)$ (sampler kernel,
bridge conditioning). This structure identifies the ratio comparison needed
to extend Lemmas 1–3 beyond the Gaussian case as a comparison of
$V-U$ (bridge vs single-point gap) between the real and floor models — the
key open quantity whose control the non-Gaussian section shows cannot be
achieved by pointwise conditions on the curvature functions alone.

**Limits of the non-Gaussian extension.** The document proves (Lemma 2$'$)
that the analog of Lemma 2 holds when the reference density is **exactly
shared** ($f_2=f_4$), but establishes by explicit numerical counterexamples
that the following natural conditions are each **insufficient** when
$f_2\ne f_4$:
- Pointwise curvature-ratio domination $\rho_2(x)\ge\rho_1(x)$.
- Pointwise curvature-gap domination $c_3(x)-c_4(x)\ge c_1(x)-c_2(x)$.
- Both conditions simultaneously.
- Both conditions with pair 2 forced to be exactly Gaussian.

This establishes that the Gaussian TV monotonicity of Theorem 1 does not
extend to general log-concave densities via any combination of pointwise
curvature conditions, motivating the alternative approach of the Main Theorem
document.

---

## 3. Main Theorem: extension to symmetric hierarchical GLMMs

**Document.** Ridge vs Floor TV Bound (current document).

**Setting.** A two-block Gibbs sampler for a one-dimensional ($q=1$)
hierarchical GLMM with $J\ge2$ groups, symmetric around $\gamma^*$ (A3).
The $\gamma$-marginal is non-Gaussian. The **floor model** replaces each
log-likelihood $\ell_i(\beta_i)$ with the ridge penalty alone, giving a
Gaussian floor $\pi_0(\gamma)$. All chains are started at $\gamma^*$.

**Assumptions.** (A1) proper group prior; (A4) log-concavity; (A2)
identifiability (marginal posterior mode $\hat\beta_i^{\rm marginal}$
finite); (A3) symmetry around $\gamma^*$; (A$_k$) finite-$k$ convergence
rate $A_k^{(1)}(\gamma)$ non-decreasing in $|\gamma-\gamma^*|$; (A$_{VS}$)
aggregate posterior variance $VS(\gamma)=\mathrm{Var}(S\mid\gamma)$
non-decreasing in $|\gamma-\gamma^*|$; (A5$'$) scalar convergence rate
monotonicity $A_1^{(1)}(\gamma)\le\cdots\le A_k^{(1)}(\gamma)\le\cdots
\le A_\infty^{(1)}(\gamma)$.

**Standardization.** Define scaling factors
$\lambda_t := 1/(\sqrt{2\pi}\pi_t(\gamma^*))$ and standardized densities
$\tilde\pi_t(z):=\lambda_t\pi_t(\lambda_tz)$, $\tilde Q_t^{(k)}(z):=
\lambda_tQ_t^{(k)}(\lambda_tz;\gamma^*)$, so that
$\tilde\pi_0(\gamma^*)=\tilde\pi_1(\gamma^*)=1/\sqrt{2\pi}$.

**Main Theorem.** *Under (A1)–(A4), (A$_k$), and (A5$'$), for every
$k\ge1$ and every $k'\in\{1,\ldots,k\}$:*
$$\boxed{\|Q_1^{(k')}(\cdot;\gamma^*)-\pi_1\|_{TV}
\;<\;
\|\tilde Q_0^{(k')}(\cdot;\gamma^*)-\tilde\pi_0\|_{TV}
\;+\;
\|\tilde\pi_0-\tilde\pi_1\|_{TV}}$$

*The bound holds for all fixed model parameters — no degeneration or
smallness condition on any parameter is required.*

**Interpretation.** The bound has two terms:

- **Geometric convergence term** $\|\tilde Q_0^{(k')}-\tilde\pi_0\|_{TV}$:
  the floor sampler's convergence to its own Gaussian target. Computable
  from Nygren (2020) / Foundational Lemmas. Decays geometrically in $k'$
  at rate $A_0^{2k'}$ where $A_0 = P_{12}P_{22,0}^{-1}P_{21}/P_{11}$
  is the floor model's scalar convergence factor.

- **Non-normality penalty** $\|\tilde\pi_0-\tilde\pi_1\|_{TV}$: the TV
  distance between the standardized Gaussian floor target and the
  standardized non-Gaussian real target. A **fixed constant** independent
  of $k'$, measuring the inherent non-Gaussianity of the posterior for
  $\gamma$. Present even in the symmetric case.

**Proof structure.** The proof proceeds via three single-crossing conditions:

*(SC1) Floor vs real sampler:* $\mathcal D_{k'}(z):=\log\tilde Q_0^{(k')}(z)
-\log\tilde Q_1^{(k')}(z)$ crosses zero exactly once on $(0,\infty)$.
Established from Lemma A.1 (curvature dichotomy, using A1–A4 and A$_k$)
and Lemma A.2 (center inequality $\mathcal D_{k'}(\gamma^*)>0$,
using A5$'$).

*(SC2) Real sampler vs real target:* $h_{k'}(z):=\log\tilde Q_1^{(k')}(z)
-\log\tilde\pi_1(z)$ crosses zero exactly once on $(0,\infty)$.
Established from (A5$'$) giving $h_{k'}''\le0$ (concavity), (A3) giving
$h_{k'}'(\gamma^*)=0$, and Lemma A.2 giving $h_{k'}(\gamma^*)>0$.

*(SC3) Nesting:* The crossing points satisfy $r_0(k')<r_1(k')$.
Established from the shape and center inequality without any parameter
conditions.

**The single-crossing lemma** then gives
$\|\tilde Q_1^{(k')}-\tilde\pi_1\|_{TV} < \|\tilde Q_0^{(k')}-\tilde\pi_1\|_{TV}$,
and the triangle inequality bounds the right-hand side by the two terms
in the Main Theorem.

**Key lemmas.**

*Lemma A.1 (curvature dichotomy).* Under (A1)–(A4) and (A$_k$):
$\mathcal D_{k'}$ is even with $\mathcal D_{k'}(z)\to-\infty$ in the
tails, and changes second-derivative sign at most once on $(0,\infty)$
from positive to negative.

*Lemma A.2 (center inequality).* Under (A1)–(A4) and (A5$'$):
$\mathcal D_{k'}(\gamma^*)>0$ for every finite $k'\ge1$.

*Lemma B.2 (variance monotonicity for symmetric logit/probit).* For
the symmetric two-group logit or probit ($y_2=1-y_1$,
$\hat\beta_2=-\hat\beta_1$): (A$_{VS}$) holds, i.e. $VS(\gamma)$ is
non-decreasing in $|\gamma-\gamma^*|$. Proved via the score identity
$\frac{d}{d\gamma}VS(\gamma)=(H/\Psi)^2[\mu_3^{(1)}(\gamma)-\mu_3^{(1)}(-\gamma)]/\Psi$
and the mirror symmetry $\mu_3^{(2)}(\gamma)=-\mu_3^{(1)}(-\gamma)$.

**"Essentially carries over" — what this means.** The Gaussian floor bound
does not transfer exactly — the non-normality penalty
$\|\tilde\pi_0-\tilde\pi_1\|_{TV}$ is an additional fixed cost. But the
bound retains the key qualitative structure of the Gaussian result:
the **geometric convergence rate** is the same as for the floor model, and
the penalty is a fixed additive constant independent of $k'$. As $k'\to\infty$,
$\|\tilde Q_0^{(k')}-\tilde\pi_0\|_{TV}\to0$ and the bound tightens to
the non-normality penalty alone.

**Symmetry requirement.** The theorem requires the model to be symmetric
around $\gamma^*$ (A3) — i.e. $\pi_0$, $\pi_1$, $Q_0^{(k)}$, $Q_1^{(k)}$
are all even functions of $\gamma-\gamma^*$. The symmetric two-group logit
and probit with $y_i=1-y_j$, $\hat\beta_j=-\hat\beta_i$ is one important
class satisfying this condition, but (A3) may hold more broadly. Extension
to the asymmetric case is an open question; see §6.

---

## 4. Sufficiency of the conditions and extensions

**The assumptions are sufficient but not necessary.** The single-crossing
conditions (SC1)–(SC3) are the necessary and sufficient conditions for the
TV inequality. The assumptions (A1)–(A4), (A$_k$), (A5$'$) are sufficient
for (SC1)–(SC3) but need not be necessary. Models that fail some sufficient
conditions may still satisfy (SC1)–(SC3):

- **Weakly asymmetric models** where $|y_i-1/2|$ is small may satisfy
  (SC1)–(SC3) even without exact symmetry (A3).
- **Models with non-log-concave likelihoods** may satisfy (A$_k$)
  empirically even if (A4) fails.
- **Multi-group models** with heterogeneous group sizes or loadings $H_i$
  may satisfy (A$_{VS}$) through cancellation effects even if individual
  (A$_V$) conditions fail.

**On (A3) and (A$_{VS}$) as sufficient but not necessary conditions for
(A5$'$).** The proof of (A5$'$) proceeds through three conditions on
$D_k(\gamma):=A_k^{(1)}(\gamma)-A_{k-1}^{(1)}(\gamma)$:

- *(i)* $D_k(\gamma^*)>0$ — follows from (A1)–(A4) alone, no (A3) needed.
- *(ii)* $D_k'(\gamma^*)=0$ — in the proof, this follows from (A3) making
  $D_k$ even. But all that is needed is that $\gamma^*$ is a critical point
  of $D_k$, which could hold without full symmetry.
- *(iii)* $D_k''(\gamma)\ge0$ — in the proof, this uses (A3) to establish
  a U-shaped likelihood ratio and (A$_{VS}$) to give $VS(\gamma)$
  non-decreasing in $|\gamma-\gamma^*|$. But the underlying requirement is
  just that $\int VS(\gamma)\,\delta(\gamma)\,d\gamma\ge0$ — which could
  hold under weaker conditions.

Together (i)–(iii) are sufficient for (A5$'$), but (A5$'$) itself is the
key assumption for the Main Theorem. A model could satisfy (A5$'$) directly
— or satisfy (i)–(iii) by some other route — without requiring (A3) or
(A$_{VS}$). In particular, a model that is not exactly symmetric around
$\gamma^*$ may still have $A_k^{(1)}(\gamma)$ non-decreasing in $k$ at
every $\gamma$ — (A5$'$) is a property of the chain's convergence behavior,
not of the model's symmetry.

**On the most technical condition.** (A$_k$) — monotonicity of the
finite-$k$ convergence rate in $|\gamma-\gamma^*|$ — is the most
substantive assumption. It is natural but not trivial to verify. Under
(A$_k$), (A5$'$) follows from (A$_{VS}$) — which is proved analytically
for symmetric logit/probit in Lemma B.2.

---

## 5. Relationship between the three results

The three results form a logical chain:
$$\underbrace{\text{Nygren (2020)}}_{\substack{\text{Gaussian, exact}\\\text{TV bound via erf}_n}}
\;\longrightarrow\;
\underbrace{\text{Foundational Lemmas}}_{\substack{\text{Gaussian hierarchical}\\\text{eigenvalue monotonicity}}}
\;\longrightarrow\;
\underbrace{\text{Main Theorem}}_{\substack{\text{Non-Gaussian, symmetric}\\\text{floor bound + non-normality penalty}}}$$

Each step generalizes the previous at a specific cost:
- Nygren (2020): exact Gaussian bound in terms of $\mathrm{erf}_n$ and eigenvalues.
- Foundational Lemmas: establishes that TV distance is monotone in eigenvalues of the precision-ratio matrix, giving the floor bound $\|Q_0^{(k)}-\pi_0\|_{TV}\le\|Q_1^{(k)}-\pi_1\|_{TV}$ for Gaussian hierarchical models — the floor sampler bounds the real sampler from above.
- Main Theorem: extends the floor-bound structure to non-Gaussian symmetric models, at the cost of adding the non-normality penalty $\|\tilde\pi_0-\tilde\pi_1\|_{TV}$.

The **floor model** — ridge-only Gaussian sampler — is the common
reference across all three results. Its TV distance is analytically
computable from Nygren (2020), monotonically bounds the Gaussian
hierarchical sampler via the Foundational Lemmas, and bounds the
non-Gaussian symmetric sampler (up to the non-normality penalty)
via the Main Theorem.

---

## 6. Open questions

*(i) Asymmetric models.* The theorem requires symmetry (A3). Extending
to asymmetric likelihoods — where $y_i\ne1-y_j$ — requires establishing
the single-crossing properties (SC1)–(SC3) without symmetry, which remains
an open problem.

*(ii) Matrix monotonicity (A5).* The stronger condition
$V_{k-1}(\gamma)\preceq V_k(\gamma)$ (matrix order) would imply (A5$'$)
but remains unproved for non-Gaussian likelihoods.

*(iii) Quantitative non-normality penalty.* The penalty
$\|\tilde\pi_0-\tilde\pi_1\|_{TV}$ is shown to be finite but not
explicitly bounded in terms of $n_i$, $\Psi_i$, $H_i$. A quantitative
bound would make the Main Theorem fully explicit.

*(iv) Multi-group asymmetric extensions.* Extending Lemma B.2 to
$J>2$ groups with asymmetric likelihoods requires generalizing the
mirror symmetry argument.

---

## 7. Asymmetric extension and prior calibration

The symmetric assumption (A3) is replaced by explicit consistency
conditions, and the main theorem extended to asymmetric models via
prior calibration. The full development is in the asymmetric document;
the key results are summarized here.

### 7.1 Replacing (A3): consistency conditions

Three conditions replace the symmetry assumption:

**(C1) Prior calibration.** Set $\mu_\gamma=\gamma^{\dagger}:=\arg\min_\gamma VS(\gamma)$.
Since $VS(\gamma)$ is independent of $\mu_\gamma$ (a property of the
posterior density alone), $\gamma^{\dagger}$ is computed first and
$\mu_\gamma$ set to match. This aligns the posterior mode with the
minimizer of $VS$ — the point of fastest mixing.

*Statistical interpretation (logit/probit).* At $\mu_\gamma=\gamma^{\dagger}$
the weighted sum of third central moments vanishes:
$\sum_i(H_i/\Psi_i)^2\mu_3^{(i)}(\gamma^{\dagger})=0$. For the
logit this means success probabilities are evenly distributed on
both sides of $0.5$ in a precision-weighted sense — the point of
maximum aggregate likelihood information.

**(C2) Mode-matched starting values.** For each $k\ge1$ and model
$t\in\{0,1\}$, $\gamma^{t,*}_k$ is the unique starting value such
that the mode of $Q_t^{(k)}(\cdot;\gamma^{t,*}_k|\mu_\gamma)$ equals
$\gamma^*_t$. Existence from the IVT; uniqueness from strict
monotonicity of the mode map $F_k^t$ under (A4).

**(C3) Integral condition.**
$\int VS(\gamma)\,\delta^{(k)}(\gamma;\gamma^{1,*}_k,\gamma^{1,*}_{k-1})\,d\gamma\ge0$
where $\delta^{(k)}:=p_{\rm sm}^k(\cdot;\gamma^{1,*}_k)-
p_{\rm sm}^{k-1}(\cdot;\gamma^{1,*}_{k-1})$. Automatic under (A3);
primitive assumption in the asymmetric case. The mode-matching (C2)
makes it plausible — both smoothing distributions are centered at
$\gamma^*$ — but does not prove it.

### 7.2 Main theorem — asymmetric case

With separate standardization factors $\lambda_t(\mu)$ shifting each
model's mode to $z=0$, the same TV inequality holds:
$$\bigl\|Q_1^{(k')}(\cdot;\gamma^{1,*}_{k'}|\mu_\gamma=\gamma^{\dagger})
-\pi_1(\cdot|\mu_\gamma=\gamma^{\dagger})\bigr\|_{TV}$$
$$<\;
\bigl\|\tilde Q_0^{(k')}(\cdot;\gamma^{0,*}_{k'}|\mu_\gamma=\gamma^{\dagger})
-\tilde\pi_0(\cdot|\mu_\gamma=\gamma^{\dagger})\bigr\|_{TV}
\;+\;
\bigl\|\tilde\pi_0(\cdot|\mu_\gamma=\gamma^{\dagger})
-\tilde\pi_1(\cdot|\mu_\gamma=\gamma^{\dagger})\bigr\|_{TV}$$
The floor convergence term (first) decays geometrically via Nygren
(2020). The non-normality penalty (second) is a fixed constant,
now incorporating both non-Gaussianity and the mode shift between
$\pi_0$ and $\pi_1$.

### 7.3 Bound for the real prior $\tilde\mu(\gamma^*)$

In practice the prior mean is $\tilde\mu(\gamma^*)$ — the value
making $\gamma^*$ the posterior mode — rather than the calibrated
$\gamma^{\dagger}$. Denote:
$$L_k := \bigl\|\tilde Q^{1,z}_{k'}(\cdot|\gamma^{\dagger})
-\tilde\pi_{1,z}(\cdot|\gamma^{\dagger})\bigr\|_{TV}, \qquad
R_k := \bigl\|\tilde Q^{1,z}_{k'}(\cdot|\tilde\mu(\gamma^*))
-\tilde\pi_{1,z}(\cdot|\tilde\mu(\gamma^*))\bigr\|_{TV}$$
By the triangle inequality with the difference-of-errors as pivot:
$$R_k \;\le\; L_k \;+\; \|\Delta_k\|_{TV}$$
where $\Delta_k:=[\tilde Q^{1,z}_{k'}(\cdot|\tilde\mu)-\tilde\pi_{1,z}(\cdot|\tilde\mu)]
-[\tilde Q^{1,z}_{k'}(\cdot|\gamma^{\dagger})-\tilde\pi_{1,z}(\cdot|\gamma^{\dagger})]$
— the change in the error (sampler minus target) as the prior mean
moves. Both brackets vanish as $k'\to\infty$, so both $L_k\to0$
and $\|\Delta_k\|_{TV}\to0$.

**Gaussian data:** $\Delta_k\equiv0$ exactly — the prior mean has
no effect on convergence in $z$-space when the likelihood is
Gaussian, since the $z$-coordinate absorbs the location-scale
change exactly.

**Gaussian approximation of targets only.** Adding and subtracting
Gaussian target approximations $\hat\pi_{1,z}(\cdot|\mu):=
N(0,1/(P_{11}\lambda_1^2(\mu)))$:
$$\|\Delta_k\|_{TV} \;\le\;
\underbrace{\|\hat\Delta_k\|_{TV}}_{\text{(C)}}
\;+\;
\underbrace{\|\tilde\pi_{1,z}(\cdot|\tilde\mu)-\hat\pi_{1,z}(\cdot|\tilde\mu)\|_{TV}}_{\text{(D)}}
\;+\;
\underbrace{\|\tilde\pi_{1,z}(\cdot|\gamma^{\dagger})-\hat\pi_{1,z}(\cdot|\gamma^{\dagger})\|_{TV}}_{\text{(E)}}$$
where $\hat\Delta_k:=[\tilde Q^{1,z}_{k'}(\cdot|\tilde\mu)-\hat\pi_{1,z}(\cdot|\tilde\mu)]
-[\tilde Q^{1,z}_{k'}(\cdot|\gamma^{\dagger})-\hat\pi_{1,z}(\cdot|\gamma^{\dagger})]$
uses exact samplers and Gaussian targets.

### 7.4 Combined bound

$$\boxed{R_k \;\le\;
\underbrace{\|\tilde Q_0^{(k')}(\cdot;\gamma^{0,*}_{k'}|\gamma^{\dagger})
-\tilde\pi_0(\cdot|\gamma^{\dagger})\|_{TV}}_{\text{(A) floor convergence, geometric}}
\;+\;
\underbrace{\|\tilde\pi_0(\cdot|\gamma^{\dagger})-\tilde\pi_1(\cdot|\gamma^{\dagger})\|_{TV}}_{\text{(B) non-normality penalty, fixed}}
\;+\;
\underbrace{\|\hat\Delta_k\|_{TV}}_{\text{(C) difference in non-normality}}
\;+\;
\underbrace{\|\tilde\pi_{1,z}(\cdot|\tilde\mu)-\hat\pi_{1,z}(\cdot|\tilde\mu)\|_{TV}}_{\text{(D) real target non-normality, fixed}}
\;+\;
\underbrace{\|\tilde\pi_{1,z}(\cdot|\gamma^{\dagger})-\hat\pi_{1,z}(\cdot|\gamma^{\dagger})\|_{TV}}_{\text{(E) calibrated target non-normality, fixed}}}$$

**Interpretation of terms:**
- **(A):** Analytic from Nygren (2020). Decays geometrically in $k'$.
- **(B):** Fixed. Approximated for large $k'$ by
  $\|\tilde\pi_0(\cdot|\gamma^{\dagger})-\tilde Q_1^{(k')}(\cdot|\gamma^{\dagger})\|_{TV}$
  — floor model vs real chain at calibrated prior, error $L_k\to0$.
- **(C):** Difference in non-normality across prior means. For small
  $k'$ driven by eigenvalue difference
  $(VS(\gamma^*)-VS(\gamma^{\dagger}))/P_{11}>0$; for large $k'$
  converges to (D)$-$(E). Exact since $A_\infty^{(1)}(\gamma^*)>A_\infty^{(1)}(\gamma^{\dagger})$
  — the real sampler has a larger convergence matrix.
- **(D):** Non-normality of real posterior. Fixed; larger than (E)
  since $\tilde\pi_{1,z}(\cdot|\tilde\mu)$ is more asymmetric.
  Approximated by $\|\tilde Q^{1,z}_{k'}(\cdot|\tilde\mu)-\hat\pi_{1,z}(\cdot|\tilde\mu)\|_{TV}$,
  error $R_k\to0$.
- **(E):** Non-normality of calibrated posterior. Fixed; minimized
  over all $\mu_\gamma$ by (C1). Approximated by
  $\|\tilde Q^{1,z}_{k'}(\cdot|\gamma^{\dagger})-\hat\pi_{1,z}(\cdot|\gamma^{\dagger})\|_{TV}$,
  error $L_k\to0$.

**Estimation.** All five terms are estimable from **two chain runs**
— one at $\tilde\mu(\gamma^*)$ and one at $\gamma^{\dagger}$ — plus
the analytical Gaussian floor model. Term (C) is the difference of
the same chain-vs-Gaussian quantities used for (D) and (E). Total
approximation error $L_k+R_k\to0$ as $k'\to\infty$.

**Behavior as $k'\to\infty$:** (A)$\to0$, (C)$\to$(D)$-$(E)$>0$;
the bound converges to (B)$+$(D) — the non-normality penalty at
$\gamma^{\dagger}$ plus the non-normality of the real target. In
the Gaussian case: (C)$=$0, (D)$=$0, (E)$=$0, recovering the main
theorem at the calibrated prior exactly.

---

## 7.5 Extension to the multivariate case under coordinate independence

Under coordinate independence each $\gamma_j$ runs its own
independent univariate chain, and $\mathbf{A}_\infty^{(1)}(\gamma)$
is diagonal with entries $\rho_j(\mu):=VS_j(\gamma^*_j(\mu))/P_{11,j}$.
The eigenvalues of the precision-ratio matrix of the $k'$-step
sampler vs target are exactly these diagonal entries — one per
coordinate — and the Foundational Lemmas of §2 apply coordinate
by coordinate via Lemma 3's chaining argument.

**Generalization of the one-eigenvalue-at-a-time formula.** The
combined bound on $R_k$ (§7.4) involves TV distances between
distributions whose precision-ratio eigenvalues differ by prior
mean. Under coordinate independence, moving the prior mean from
$\gamma^{\dagger}$ to $\tilde\mu(\gamma^*)$ changes each diagonal
eigenvalue from $\rho_j(\gamma^{\dagger})$ to $\rho_j(\tilde\mu(\gamma^*))>\rho_j(\gamma^{\dagger})$
(since $VS_j(\gamma^*_j(\tilde\mu))>VS_j(\gamma^{\dagger}_j)$ by
the univariate argument). By Lemma 3 (coordinatewise domination by
chaining), the TV distance increases monotonically as the eigenvalues
are moved one coordinate at a time.

Define the chain of diagonal precision-ratio matrices:
$$M^{(0)} \;=\; \mathrm{diag}\bigl(\rho_1(\gamma^{\dagger}),\ldots,\rho_p(\gamma^{\dagger})\bigr)$$
$$M^{(m)} \;=\; \mathrm{diag}\bigl(\rho_1(\tilde\mu),\ldots,\rho_m(\tilde\mu),
\rho_{m+1}(\gamma^{\dagger}),\ldots,\rho_p(\gamma^{\dagger})\bigr),
\quad m=1,\ldots,p$$
$$M^{(p)} \;=\; \mathrm{diag}\bigl(\rho_1(\tilde\mu),\ldots,\rho_p(\tilde\mu)\bigr)$$
At each step only coordinate $m$ changes:
$\rho_m(\gamma^{\dagger})\to\rho_m(\tilde\mu)>\rho_m(\gamma^{\dagger})$.
Lemma 2 applied at each step gives a non-negative increment, and
chaining gives:
$$\bigl\|N(0,\Sigma_1^{(0)})-N(0,\Sigma_2^{(0)})\bigr\|_{TV}
\;\le\;
\bigl\|N(0,\Sigma_1^{(p)})-N(0,\Sigma_2^{(p)})\bigr\|_{TV}$$
— the TV distance at the calibrated eigenvalues $M^{(0)}$ is
bounded above by the TV distance at the real eigenvalues $M^{(p)}$.

**Term (A) under coordinate independence.** The floor chain is
a multivariate Gaussian under coordinate independence — diagonal
covariance throughout. Its $k'$-step precision-ratio matrix against
the floor target is $\mathrm{diag}(1/(1-ho_{j,0}^{2k'}))_{j=1}^p$
with $ho_{j,0}:=VS_{0,j}/P_{11,j}$ the $j$-th floor eigenvalue.
By Nygren (2020) Theorem 3 applied to the diagonal case:
$$\text{(A)} \;=\;
\left\|N\!\left(0,\mathrm{diag}\!\left(\frac{1}{1-\rho_{j,0}^{2k'}}\right)_{j=1}^p\right)
- N(0,I_p)\right\|_{TV}$$
exactly, computable from the $\mathrm{erf}_n$ formula with the $p$
diagonal eigenvalue ratios. This is directly the multivariate
generalization of the univariate floor bound — the $p$ coordinates
each contribute their own Nygren term, assembled into the single
$p$-dimensional TV distance via the diagonal structure.

**Term (B) under coordinate independence.** The non-normality
penalty measures the TV distance between the standardized floor
and real stationary distributions at the calibrated prior
$\gamma^{\dagger}$. Under coordinate independence both factorize
and the TV distance satisfies:
$$\text{(B)} \;=\;
\bigl\|\tilde\pi_0(\cdot|\gamma^{\dagger})-\tilde\pi_1(\cdot|\gamma^{\dagger})\bigr\|_{TV}
\;\le\; \sum_{j=1}^p
\bigl\|\tilde\pi_{0,j}(\cdot|\gamma^{\dagger}_j)
-\tilde\pi_{1,j}(\cdot|\gamma^{\dagger}_j)\bigr\|_{TV}$$
via the tensorization inequality. Each coordinate's contribution
is a fixed univariate non-normality penalty of the same type as
the univariate term (B) — the TV distance between the $j$-th
coordinate's standardized floor and real stationary densities at
$\gamma^{\dagger}_j$, bounded by the foundational lemmas
applied coordinate-wise.

**Term (C) under coordinate independence.** Term (C)
$=\|\hat\Delta_k\|_{TV}$ involves the **exact non-Gaussian samplers**
against Gaussian target approximations — it is not a Gaussian TV
distance and is not computable in closed form. Under coordinate
independence it is bounded by a coordinate-wise sum:
$$\|\hat\Delta_k\|_{TV} \;\le\; \sum_{j=1}^p
\bigl\|[\tilde Q^{1,z}_{k',j}(\cdot|\tilde\mu_j)
-\hat\pi_{1,z,j}(\cdot|\tilde\mu_j)]
-[\tilde Q^{1,z}_{k',j}(\cdot|\gamma^{\dagger}_j)
-\hat\pi_{1,z,j}(\cdot|\gamma^{\dagger}_j)]\bigr\|_{TV}$$
— $p$ coordinate-wise terms of the same non-Gaussian character as
the univariate term (C). Each coordinate's contribution measures
how much the non-normality of the exact $k'$-step chain in that
coordinate differs across the two prior means.

**Gaussian approximation to (C) for small $k'$.** When both chains
are approximately Gaussian (small $k'$, both close to their
mode-matched starts), each coordinate's contribution to (C) is
approximately the difference of two univariate Gaussian TV distances,
giving $p$ increments:
$$\Delta_{k,m}^{\rm approx} \;:=\;
\bigl\|N(0,\sigma^2_{k',m}(\tilde\mu_m))
-N(0,1/(P_{11,m}\lambda_{1,m}^2(\tilde\mu_m)))\bigr\|_{TV}$$
$$-\;
\bigl\|N(0,\sigma^2_{k',m}(\gamma^{\dagger}_m))
-N(0,1/(P_{11,m}\lambda_{1,m}^2(\gamma^{\dagger}_m)))\bigr\|_{TV}
\;\ge\;0$$
where $\sigma^2_{k',m}(\mu):=(1-\rho_m^{2k'}(\mu))/(P_{11,m}\lambda_{1,m}^2(\mu))$
is the residual variance fraction in coordinate $m$. Each increment
is analytically computable from the Nygren (2020) $\mathrm{erf}_n$
formula applied to the single changed eigenvalue $\rho_m$, and the
one-eigenvalue-at-a-time chain of the Foundational Lemmas confirms
each increment is non-negative: $\rho_m(\tilde\mu)>\rho_m(\gamma^{\dagger})$
gives larger residual variance and hence larger TV from the Gaussian
approximation target.

**Terms (D) and (E) under coordinate independence.** Both terms
measure the non-normality of a stationary target against its
Gaussian approximation. Under coordinate independence the
standardized target factorizes as a product of $p$ independent
univariate densities, and the TV distance satisfies:
$$\text{(D)} \;=\;
\bigl\|\tilde\pi_{1,z}(\cdot|\tilde\mu)-\hat\pi_{1,z}(\cdot|\tilde\mu)\bigr\|_{TV}
\;\le\; \sum_{j=1}^p
\bigl\|\tilde\pi_{1,z,j}(\cdot|\tilde\mu_j)
-N\!\left(0,\frac{1}{P_{11,j}\lambda_{1,j}^2(\tilde\mu_j)}\right)\bigr\|_{TV}$$
$$\text{(E)} \;=\;
\bigl\|\tilde\pi_{1,z}(\cdot|\gamma^{\dagger})-\hat\pi_{1,z}(\cdot|\gamma^{\dagger})\bigr\|_{TV}
\;\le\; \sum_{j=1}^p
\bigl\|\tilde\pi_{1,z,j}(\cdot|\gamma^{\dagger}_j)
-N\!\left(0,\frac{1}{P_{11,j}\lambda_{1,j}^2(\gamma^{\dagger}_j)}\right)\bigr\|_{TV}$$
via the tensorization inequality for TV distances under independence.
Each coordinate's contribution is a univariate non-normality term
of the same type as the univariate (D) and (E) — the TV distance
between a standardized univariate posterior and its Gaussian
approximation, fixed in $k'$, and bounded by the foundational
lemmas applied coordinate-wise. The ordering (D)$_j>$(E)$_j$ holds
for each $j$ by the same argument as in the univariate case: the
posterior at $\tilde\mu_j$ is more asymmetric than at $\gamma^{\dagger}_j$
since the coordinate-wise balance condition holds at $\gamma^{\dagger}_j$
but not at $\gamma^*_j(\tilde\mu_j)\ne\gamma^{\dagger}_j$.

**$\gamma^{\dagger}$ as a vector.** Under coordinate independence
$\gamma^{\dagger}_j=\arg\min_{\gamma_j}VS_j(\gamma_j)$ for each
$j$ independently, satisfying the coordinate-wise balance condition:
$$\sum_i\left(\frac{H_{i,j}}{\Psi_{i,j}}\right)^2\mu_{3,j}^{(i)}(\gamma^{\dagger}_j)
\;=\; 0 \qquad j=1,\ldots,p$$
— one univariate balance condition per coordinate, each determining
$\mu_{\gamma,j}=\gamma^{\dagger}_j$ independently.

**Full multivariate bound under coordinate independence.**
Writing $\rho_{j,0}:=VS_{0,j}/P_{11,j}$ for the floor chain's
$j$-th eigenvalue and $\rho_j(\mu):=VS_j(\gamma^*_j(\mu))/P_{11,j}$
for the real chain's $j$-th eigenvalue at prior mean $\mu$:

$$\boxed{R_k \;\le\;
\underbrace{\left\|N\!\left(0,\mathrm{diag}\!
\left(\frac{1}{1-\rho_{j,0}^{2k'}}\right)_{j=1}^p\right)
- N(0,I_p)\right\|_{TV}}_{\text{(A) Nygren (2020) floor, exact}}
\;+\;
\underbrace{\sum_{j=1}^p
\bigl\|\tilde\pi_{0,j}(\cdot|\gamma^{\dagger}_j)
-\tilde\pi_{1,j}(\cdot|\gamma^{\dagger}_j)\bigr\|_{TV}}_{\text{(B) coordinate non-normality penalty}}
\;+\;
\underbrace{\sum_{j=1}^p\bigl\|[\tilde Q^{1,z}_{k',j}(\cdot|\tilde\mu_j)-\hat\pi_{1,z,j}(\cdot|\tilde\mu_j)]-[\tilde Q^{1,z}_{k',j}(\cdot|\gamma^{\dagger}_j)-\hat\pi_{1,z,j}(\cdot|\gamma^{\dagger}_j)]\bigr\|_{TV}}_{\text{(C) coordinate-wise non-Gaussianity difference}}
\;+\;
\underbrace{\sum_{j=1}^p
\left\|\tilde\pi_{1,z,j}(\cdot|\tilde\mu_j)
-N\!\left(0,\frac{1}{P_{11,j}\lambda_{1,j}^2(\tilde\mu_j)}\right)\right\|_{TV}}_{\text{(D) real target non-normality}}
\;+\;
\underbrace{\sum_{j=1}^p
\left\|\tilde\pi_{1,z,j}(\cdot|\gamma^{\dagger}_j)
-N\!\left(0,\frac{1}{P_{11,j}\lambda_{1,j}^2(\gamma^{\dagger}_j)}\right)\right\|_{TV}}_{\text{(E) calibrated target non-normality}}}$$

Each term has the same character as its univariate counterpart
in §7.4: (A) is analytic from Nygren (2020); (B), (D), (E) are
fixed non-normality terms bounded coordinate-wise; (C) is the
coordinate-wise sum of exact non-Gaussian differences — not
computable in closed form, but approximated for small $k'$ by
$p$ Nygren (2020) increments $\Delta_{k,m}^{\rm approx}\ge0$,
one per coordinate, each driven by the eigenvalue increase
$\rho_m(\tilde\mu)>\rho_m(\gamma^{\dagger})$; and
(D)$_j>$(E)$_j$ for each $j$ since the real posterior at
$\tilde\mu_j$ is more asymmetric than the calibrated posterior
at $\gamma^{\dagger}_j$.
