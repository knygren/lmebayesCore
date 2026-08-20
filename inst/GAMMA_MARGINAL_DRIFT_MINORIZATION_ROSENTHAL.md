# \(\gamma\)-marginal drift, minorization, and Rosenthal TV bound

Companion to `restricted_gibbs_minorization _v4.md` (Theorem 1),
`GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md`, `JOINT_GAMMA_BETA_TV_CERTIFICATE.md`,
and `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md`.

**Status:** draft design note — constants are structurally explicit; some
inequalities still require a formal lemma on \(\widetilde B(\delta_2)\).

---

## 0. Scope

This note certifies the **\(\beta\)-restricted** two-block Gibbs chain on the
\(\beta\)-safe set \(\widetilde B(\delta_2)\) (block 2 updated second, matching
Chapter C03). One sweep of the **\(\gamma\)-marginal kernel**
\(P_{\gamma\mid\widetilde B(\delta_2)}\):

\[
\gamma_{n-1}
\;\xrightarrow{\;\beta_n\sim\pi(\beta\mid\gamma_{n-1},y)\ \text{on }\widetilde B(\delta_2)\;}
\;\xrightarrow{\;\gamma_n\sim\pi(\gamma\mid\beta_n,y)\;}
\gamma_n.
\]

Each block is a genuine Gibbs update for the **\(\beta\)-truncated target**

\[
\pi(\gamma,\beta\mid y,\,\beta\in\widetilde B(\delta_2))
\;=\;
\frac{\pi(\gamma,\beta\mid y)\,\mathbf 1_{\beta\in\widetilde B(\delta_2)}}
{\pi(\beta\in\widetilde B(\delta_2)\mid y)}
\]

(`restricted_gibbs_minorization _v4.md` Definition 1; \(\beta\)-slice only — not
“run unrestricted and repair excursions”). Write the induced **\(\gamma\)-target**

\[
\pi_{\gamma\mid\widetilde B(\delta_2)}(\cdot)
\;:=\;
\pi_\gamma\bigl(\cdot\mid\widetilde B(\delta_2),y\bigr)
\;=\;
\frac{\displaystyle\int_{\widetilde B(\delta_2)}\pi(\gamma,\beta\mid y)\,d\beta}
{\displaystyle\int_{\mathbb R^q}\!\int_{\widetilde B(\delta_2)}\pi(\gamma,\beta\mid y)\,d\beta\,d\gamma}.
\]

**Rosenthal TV object (this note).**

\[
\bigl\|P_{\gamma\mid\widetilde B(\delta_2)}^{\,k}(\gamma_0,\cdot)
-\pi_{\gamma\mid\widetilde B(\delta_2)}\bigr\|_{TV}.
\]

**Full \(\gamma\)-posterior (optional outer step).** Triangle inequality against
\(\pi_\gamma(\cdot\mid y)\) adds a **static \(\beta\)-truncation** term, e.g.\
\(\|\pi_{\gamma\mid\widetilde B(\delta_2)}-\pi_\gamma\|_{TV}\le \pi_\beta(\widetilde B(\delta_2)^{\,c})\le\delta_2\)
(Lemma 1 pattern on the \(\beta\) marginal). That step is separate from §3.1.

We state **minorization** and **geometric drift** for
\(P_{\gamma\mid\widetilde B(\delta_2)}\) on
\(\widetilde{\mathcal R}=\widetilde C_d\times\widetilde B(\delta_2)\), with constants
from the floor spectrum \(\kappa_i^{\mathrm{LB}}(\delta_2)\), and apply
**Rosenthal (1995)** to \(\pi_{\gamma\mid\widetilde B(\delta_2)}\).

Notation follows C05: \(\gamma^\star\) from `population_mode()`, \(P_{11}\),
\(H_j\), \(P_b\), \(V_j(\gamma)\), \(\tilde J\), \(\kappa_i\), \(\Psi\),
\(\widetilde C_d=\{\Psi\le d\}\). Shorthand:
\(\pi_{\gamma\mid\widetilde B}:=\pi_\gamma(\cdot\mid\widetilde B(\delta_2),y)\),
\(P_\gamma:=P_{\gamma\mid\widetilde B(\delta_2)}\) when \(\delta_2\) is fixed.

### Terminology (read this first)

| Phrase | Meaning in this note |
|--------|----------------------|
| **One sweep** | One application of \(P_{\gamma\mid\widetilde B(\delta_2)}\): draw \(\beta_n\in\widetilde B(\delta_2)\), then \(\gamma_n\) (§0). |
| **Foster / drift condition** | A **single-sweep** inequality: \((P_\gamma V)(\gamma)=\mathbb E[V(\gamma_n)\mid\gamma_{n-1}=\gamma]\le \lambda V(\gamma)+b\,\mathbb I_{C^c}(\gamma)\). |
| **§2** | Derives \((\lambda,b)\) by expanding that **one** expectation (with \(\beta_n\) integrated out). |
| **§3 (Rosenthal)** | Uses \((\lambda,b,\varepsilon)\) to bound **\(k\) sweeps** — \(P_\gamma^k\), not a separate “\(k\)-step drift” formula. |

### 0.2 \(\widetilde B(\delta_2)\): marginal mode + \(r_{\mathrm{Gauss}}\) on true \(\Xi\)

Full design note: **`inst/BETA_MARGINAL_MODE_LEVELSET.md`**.

**Anchor.** \(\beta^\dagger=\arg\max_\beta \widetilde\pi(\beta\mid y)\) after integrating
\(\gamma\) (`LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §5). **Not**
\(b(\gamma^\star)=E[\beta\mid\gamma^\star,y]\) (wrong origin for \(\Xi\)) and not
the joint ICM mode unless it equals \(\beta^\dagger\).

**Set.**

\[
\Xi(\beta)=f(\beta)-f(\beta^\dagger),\quad
\widetilde B(\delta_2)=\{\Xi\le r_{\mathrm{Gauss}}(n,\delta_2)\},\quad
r_{\mathrm{Gauss}}(n,\delta_2)=\tfrac12\chi^2_{n,\,1-\delta_2},\ n=Jp_{\mathrm{re}}.
\]

**Mass.**

- **Asymptotic** (\(n_j\to\infty\)): \(\pi_\beta(\widetilde B^c)=\delta_2+O(N^{-1/2})\)
  under BvM / Laplace at \(\beta^\dagger\).
- **Finite \(n\):** \(\delta_2\) is the Laplace **design** target; validate true tail
  by integration if a proof is required (not via Prop (P2), which is loose at fixed \(n\)).

**Floors.** Same \(\widetilde B(\delta_2)\): support functions on \(\{\Xi\le r\}\) →
\(\omega_{j,i}(\delta_2)\), \(P_{22,j}^{\mathrm{LB}}\), \(\kappa_i^{\mathrm{LB}}(\delta_2)\)
(§2). Implementation: `R/c05_beta_marginal_set.R`, `group_precision_floor()` migration
(BETA note §7).

**Not what we mean by “one-step”.**

- Not the \(\gamma\)-draw **alone** with \(\beta\) held fixed (that is only an internal tower step in the proof).
- Not an \(n\)-step drift bound like \(\lambda^n V + b/(1-\lambda)\) (that comes **after** Foster, in escape/coupling lemmas — Rosenthal §3.3).

**Logic flow.**

```text
β-restricted sweep  P_{γ|B̃}  ──►  Foster (λ(δ₂), b(δ₂))     [§1.2, §2]
       │
       k sweeps ──►  Rosenthal TV to π_γ(·|B̃(δ₂),y)          [§3]
       │
       (optional) + π_β(B̃^c) ≈ δ₂  ──►  full π_γ              [§0.2, §8]
       │
       (optional) C05 γ-only ──►  (1−ε(d))^n + δ_γ             [§4]
```

---

## 1. Fundamental conditions (Meyn–Tweedie / Rosenthal)

Let \(P_\gamma:=P_{\gamma\mid\widetilde B(\delta_2)}\) be the **\(\beta\)-restricted**
marginal kernel (§0) and
\(\pi_{\gamma\mid\widetilde B}:=\pi_\gamma(\cdot\mid\widetilde B(\delta_2),y)\) its
Gibbs-invariant target.

### 1.1 Minorization (Doeblin on a small set)

There exist a compact **small set** \(C\subseteq\mathbb R^q\), a constant
\(\varepsilon\in(0,1]\), and a probability measure \(\nu\) on \(\mathbb R^q\) such
that

\[
\boxed{
P_\gamma(\gamma,\cdot)\ \ge\ \varepsilon\,\nu(\cdot)
\qquad\text{for all }\gamma\in C.
}
\]

**C05 instantiation on \(\widetilde{\mathcal R}\).** Lemma 17(a) of
`restricted_gibbs_minorization _v4.md` applies to the **\(\beta\)-restricted**
one-sweep kernel \(P_\gamma\) on \(\widetilde B(\delta_2)\) (minorization verified
on \(\gamma\in\widetilde C_d\); Foster constants from §2 use
\(\kappa_i^{\mathrm{LB}}(\delta_2)\)):

\[
C=\widetilde C_d=\{\gamma:\Psi(\gamma)\le d\},
\qquad
\nu=Q,
\qquad
\varepsilon=\varepsilon(d),
\]

with \(\varepsilon(d)=e^{-d}\varepsilon(\gamma^\star)\) and **untruncated** refresh
\(Q=N(\gamma^\star,P_{11}^{-1})\). On \(\gamma\in\widetilde C_d\),

\[
P_\gamma(\gamma,A)\ \ge\ \varepsilon(d)\,Q(A)
\qquad\text{for all Borel }A.
\]

(Lemma 17(b) / `certificate()` use the **doubly** restricted chain on
\(\widetilde C_d\times\widetilde B(\delta_2)\) with
\(\varepsilon_{\mathrm{rest}}(d)=\varepsilon(d)\,Q(\widetilde C_d)\) against
\(Q_{\widetilde C_d}\) — a different normalization.)

**Design.** \(\widetilde C_d\) is the **small set where minorization is verified**;
\(\widetilde B(\delta_2)\) is where the **sampler and target are conditioned**.
Shrinking \(d\) shrinks \(\widetilde C_d\) and **raises** \(\varepsilon(d)\).
Rosenthal §3 bounds drift off \(\widetilde C_d\) via \(b(\delta_2)\) and
\(V(\gamma_0)\). The gap to the **untruncated** \(\pi_\gamma(\cdot\mid y)\) is
\(\pi_\beta(\widetilde B(\delta_2)^{\,c})\le\delta_2\) (§0), not a term inside §3.1.

### 1.2 Geometric drift (Foster–Lyapunov)

There exist \(V:\mathbb R^q\to[1,\infty)\), \(\lambda\in(0,1)\), and \(b\in[0,\infty)\)
such that

\[
\boxed{
(P_\gamma V)(\gamma)
:= \mathbb E\bigl[V(\gamma_n)\mid\gamma_{n-1}=\gamma\bigr]
\ \le\
\lambda\,V(\gamma)+b\,\mathbb I_{C^c}(\gamma)
\qquad\text{for all }\gamma.
}
\]

(Standard convention: the constant \(b\) is added on the **complement** \(C^c\); on
\(C\) the inequality is \(PV\le\lambda V\) without extra \(b\).)

**Lyapunov function (MT shift).**

\[
\boxed{
V(\gamma)=1+V_0(\gamma),
\qquad
V_0(\gamma)=\tfrac12\|\gamma-\gamma^\star\|_{P_{11}}^2.
}
\]

Complementarity gives \(\Psi(\gamma)\le V_0(\gamma)\), so
\(\{\|\gamma-\gamma^\star\|_{P_{11}}^2\le 2d\}\subseteq\widetilde C_d\).

---

## 2. Verifying Foster: constants \((\lambda,b)\) for one sweep

Section 1.2 asks for constants such that **one application of \(P_\gamma\)** satisfies
\((P_\gamma V)(\gamma)\le \lambda V(\gamma)+b\,\mathbb I_{C^c}(\gamma)\). This section
computes \((\lambda,b)\). Rosenthal (§3) then uses the same constants for **\(k\)
sweeps**.

**Floor coupling spectrum (depends on \(\delta_2\)).** Fix the compact
\(\beta\)-safe set \(\widetilde B(\delta_2)\) from Step 1 (§6). Weight floors
\(\omega_{j,i}(\delta_2):=\inf_{\beta\in\widetilde B(\delta_2)} w_{j,i}(\beta)\) determine
\(P_{22,j}^{\mathrm{LB}}(\delta_2)=P_b+\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)\),
hence \(S^{\mathrm{LB}}(\delta_2)\), \(A^{\mathrm{LB}}(\delta_2)\), and eigenvalues

\[
\kappa_1^{\mathrm{LB}}(\delta_2),\ldots,\kappa_q^{\mathrm{LB}}(\delta_2)
=\mathrm{eig}\bigl(A^{\mathrm{LB}}(\delta_2)\bigr),
\qquad
\kappa_{\max}^{\mathrm{LB}}(\delta_2):=\max_i\kappa_i^{\mathrm{LB}}(\delta_2).
\]

All drift constants below are **functions of \(\delta_2\)** through this spectrum
(unless the prior-only fallback is used). **Monotonicity (typical):** smaller
\(\delta_2\Rightarrow\) larger \(\widetilde B(\delta_2)\Rightarrow\) lower
\(\omega_{j,i}(\delta_2)\Rightarrow\) larger \(\kappa_i^{\mathrm{LB}}(\delta_2)\) and
worse Foster/Rosenthal drift — not better. C05 deficiency B-weights
\(w_i=\kappa_i/(1-\kappa_i)\) at \((\gamma^\star,\beta)\) are a separate operating
spectrum; do not confuse them with \(\kappa_i^{\mathrm{LB}}(\delta_2)\).

Fix \(\gamma_{n-1}=\gamma\). Conditional on a drawn \(\beta_n\), the \(\gamma\)-refresh is

\[
\gamma_n\mid\beta_n,\gamma
\sim
N\!\bigl(m(\beta_n),\,P_{11}^{-1}\bigr),
\qquad
m(\beta)=P_{11}^{-1}\Bigl(\Lambda_\gamma\mu_0+\sum_{j=1}^J H_j^\top P_b\,\beta_j\Bigr).
\]

Define \(M(\gamma)\) as the posterior mean map (same formula with
\(b_j(\gamma):=E[\beta_j\mid\gamma,y]\) replacing \(\beta_j\)).

### 2.1 The Foster expectation \((P_\gamma V)(\gamma)\)

By definition of the marginal chain, \((P_\gamma V)(\gamma)=\mathbb E[V(\gamma_n)\mid\gamma_{n-1}=\gamma]\) integrates **both** random draws in one sweep. Using the tower property (first condition on \(\beta_n\), then on \(\gamma_n\)):

\[
\boxed{
\mathbb E[V_0(\gamma_n)\mid\gamma]
=
\mathbb E_{\beta_n}\Bigl[\tfrac12\|m(\beta_n)-\gamma^\star\|_{P_{11}}^2\Bigr]
+\tfrac q2,
}
\]

\[
\boxed{
\mathbb E[V(\gamma_n)\mid\gamma]
=
1+\mathbb E_{\beta_n}\Bigl[\tfrac12\|m(\beta_n)-\gamma^\star\|_{P_{11}}^2\Bigr]
+\tfrac q2.
}
\]

### 2.2 Contraction rate \(\lambda\)

On \(\widetilde B(\delta_2)\), uniform weight floors give a coupling ceiling
\(\kappa_{\max}^+(\delta_2)<1\) from

\[
A=P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2},
\qquad
\kappa_{\max}^+(\delta_2):=\sup_{\beta\in\widetilde B(\delta_2)}\lambda_{\max}(A).
\]

Chapter C03 Claim 2 / EM linearization bound the **mean map**. With
\(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\) from the floor blocks (and
\(\kappa_{\max}^+(\delta_2)\ge\kappa_{\max}^{\mathrm{LB}}(\delta_2)\) from a
worst-case scan on \(\widetilde B(\delta_2)\)):

\[
\|M(\gamma)-\gamma^\star\|_{P_{11}}^2
\ \le\
2\lambda(\delta_2)\,V_0(\gamma),
\qquad
\boxed{
\lambda(\delta_2)=\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2
\;=\;
\bigl(\max_i \kappa_i^{\mathrm{LB}}(\delta_2)\bigr)^2.
}
\]

### 2.3 Split of \(C_\beta\)

Because \(m(\beta)\) is linear in \(\beta\) and \(\mathbb E_{\beta_n}[m(\beta_n)\mid\gamma]=M(\gamma)\),

\[
\mathbb E_{\beta_n}\|m(\beta_n)-\gamma^\star\|_{P_{11}}^2
=
\|M(\gamma)-\gamma^\star\|_{P_{11}}^2
+
\mathbb E_{\beta_n}\|m(\beta_n)-M(\gamma)\|_{P_{11}}^2.
\]

The second term is the \(\beta\)-fluctuation of the refresh center:

\[
m(\beta_n)-M(\gamma)
=
P_{11}^{-1}\sum_{j=1}^J H_j^\top P_b\bigl(\beta_{j,n}-b_j(\gamma)\bigr),
\]

\[
\boxed{
C_\beta^{(\mathrm{var})}(\gamma)
:=
\tfrac12\,\mathbb E_{\beta_n}\|m(\beta_n)-M(\gamma)\|_{P_{11}}^2
=
\tfrac12\,\mathrm{tr}\!\Bigl(
\sum_{j=1}^J H_j^\top P_b\,V_j(\gamma)\,P_b H_j\,P_{11}^{-1}
\Bigr).
}
\]

On \(\widetilde B(\delta_2)\), Brascamp–Lieb under (H2) gives \(V_j(\gamma)\preceq
\bigl(\text{precision of }\pi(\beta_j\mid\gamma,y)\bigr)^{-1}\). The **prior-only**
ceiling \(V_j\preceq P_b^{-1}\) is always valid but **loose** when data contribute
precision. On \(\widetilde B\), Lemma 4.3 of `JOINT_GAMMA_BETA_TV_CERTIFICATE.md`
lower-bounds data precision, so block-2 precision satisfies

\[
P_{22,j}(\beta_j)\ \succeq\ P_{22,j}^{\mathrm{LB}}(\delta_2)
:= P_b+\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2),
\qquad
\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)
:= D_j^\top\operatorname{diag}(\omega_{j,i}(\delta_2))\,D_j,
\]

with \(\omega_{j,i}(\delta_2)=\inf_{\widetilde B(\delta_2)} w_{j,i}(\beta)\). Hence the **tight uniform
upper bound**

\[
\boxed{
V_j(\gamma,\beta)\ \preceq\ \bigl(P_{22,j}^{\mathrm{LB}}(\delta_2)\bigr)^{-1}
\;=\;\bigl(P_b+\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)\bigr)^{-1}
\ \prec\ P_b^{-1}
\qquad
\forall\beta\in\widetilde B(\delta_2).
}
\]

(Larger data precision \(\Rightarrow\) smaller conditional covariance.) Define
\(S^{\mathrm{LB}}(\delta_2)\), \(A^{\mathrm{LB}}(\delta_2)\), and \(C_\beta^{+}(\delta_2)\)
from the floor blocks:

\[
S^{\mathrm{LB}}(\delta_2)
:= \sum_{j=1}^J
H_j^\top P_b\,\bigl(P_{22,j}^{\mathrm{LB}}(\delta_2)\bigr)^{-1}P_b H_j,
\qquad
A^{\mathrm{LB}}(\delta_2)
:= P_{11}^{-1/2}\,S^{\mathrm{LB}}(\delta_2)\,P_{11}^{-1/2},
\qquad
\tilde J^{\mathrm{LB}}(\delta_2) := P_{11}^{-1}S^{\mathrm{LB}}(\delta_2).
\]

Then \(\mathrm{tr}(A^{\mathrm{LB}}(\delta_2))=\sum_i\kappa_i^{\mathrm{LB}}(\delta_2)\) and

\[
\boxed{
C_\beta^{+}(\delta_2)
=
\tfrac12\,\mathrm{tr}(A^{\mathrm{LB}}(\delta_2))
=
\tfrac12\,\mathrm{tr}(\tilde J^{\mathrm{LB}}(\delta_2))
=
\tfrac12\sum_{i=1}^q \kappa_i^{\mathrm{LB}}(\delta_2).
}
\]

(Matrix form, equivalent:)
\[
C_\beta^{+}(\delta_2)
=
\tfrac12\,\mathrm{tr}\!\Bigl(
\sum_{j=1}^J
H_j^\top P_b\,\bigl(P_{22,j}^{\mathrm{LB}}(\delta_2)\bigr)^{-1}P_b H_j
\,P_{11}^{-1}
\Bigr).
\]

**Prior-only fallback** (no \(\widetilde B\) / no weight floor): replace
\((P_b+\underline{\mathcal P}_{j,\mathrm{data}})^{-1}\) by \(P_b^{-1}\), giving the
looser \(C_\beta^{+,\mathrm{prior}}=\tfrac12\mathrm{tr}(\sum_j H_j^\top P_b H_j P_{11}^{-1})\).

\[
C_\beta^{(\mathrm{con})}(\gamma)
:=
\tfrac12\|M(\gamma)-\gamma^\star\|_{P_{11}}^2-\lambda(\delta_2) V_0(\gamma)
\ \le\ 0
\]

when the C03 bound is tight. Any positive slack from a cruder contraction constant
is absorbed into \(C_\beta\).

**Total \(\beta\)-contribution:**

\[
\boxed{
C_\beta(\gamma)
:=
\mathbb E_{\beta_n}\Bigl[\tfrac12\|m(\beta_n)-\gamma^\star\|_{P_{11}}^2\Bigr]
-\lambda(\delta_2) V_0(\gamma)
=
C_\beta^{(\mathrm{con})}(\gamma)+C_\beta^{(\mathrm{var})}(\gamma).
}
\]

Uniform on \(\widetilde B(\delta_2)\): \(C_\beta(\gamma)\le C_\beta^{+}(\delta_2)\).

### 2.4 Closing the Foster inequality

§2.1 gives \((P_\gamma V_0)(\gamma)\). §2.2–2.3 bound the \(\beta\)-averaged refresh-center term so that

\[
\mathbb E[V(\gamma_n)\mid\gamma]
\le
\lambda(\delta_2) V(\gamma)+\bigl(1-\lambda(\delta_2)\bigr)+\tfrac q2+C_\beta(\gamma).
\]

On \(\widetilde{\mathcal R}\) (or on \(C^c\) in the standard Foster form):

\[
\boxed{
b(\delta_2)=1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+C_\beta^{+}(\delta_2),
\qquad
\lambda(\delta_2)=\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2,
\qquad
C_\beta^{+}(\delta_2)=\tfrac12\sum_{i=1}^q \kappa_i^{\mathrm{LB}}(\delta_2).
}
\]

| Term | Meaning |
|------|---------|
| \(1-\lambda(\delta_2)\) | \(1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\) — MT shift |
| \(q/2\) | Gaussian refresh variance (\(q=\dim\gamma\)) |
| \(C_\beta^{+}(\delta_2)\) | \(\tfrac12\sum_i \kappa_i^{\mathrm{LB}}(\delta_2)=\tfrac12\,\mathrm{tr}(A^{\mathrm{LB}}(\delta_2))\) |

This is exactly the Foster condition §1.2 with
\(b=b(\delta_2)=1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+q/2+C_\beta^{+}(\delta_2)\)
on the
region where the bounds hold (typically \(\gamma\in C^c\) or uniformly on \(\widetilde{\mathcal R}\)).
Iterating Foster **\(k\) times** (standard MT) yields escape bounds
\(P_\gamma^k(\gamma,C^c)\lesssim \lambda^k V(\gamma)+b/(1-\lambda)\); Rosenthal packages
that together with minorization into §3.1.

**Requires:** compact \(\widetilde B(\delta_2)\) with weight floor \(\omega_{j,i}(\delta_2)>0\) so
\(\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)\succ 0\),
\(\kappa_{\max}^+(\delta_2)<1\), and \(C_\beta^{+}(\delta_2)<\infty\) uniformly.

---

## 3. Rosenthal (1995) total-variation bound

**Reference.** Rosenthal, J.S. (1995). *Minorization conditions and convergence
rates for Markov chain Monte Carlo.* JASA 90(430), 558–566, Theorem 12.

Assume §1.1–1.2 for the **\(\beta\)-restricted** kernel \(P_\gamma\) (§0) with small set
\(C=\widetilde C_d\). Write \(\lambda=\lambda(\delta_2)\), \(b=b(\delta_2)\) from §2.
For any start \(\gamma_0\), any \(k\in\mathbb N\), and any
\(\alpha\in\bigl(\lambda(\delta_2),1\bigr)\),

\[
\boxed{
\bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\bigr\|_{TV}
\ \le\
(1-\varepsilon)^{\lfloor\alpha k\rfloor}
\;+\;
\frac{U}{\alpha}
\left(
\frac{1+2b+\lambda V(\gamma_0)}{1+2b/(1-\lambda)}
\right)
\alpha^{k}.
}
\]

### 3.1 Complete Rosenthal formula (all constants explicit)

**Total-variation bound** (Rosenthal 1995, Theorem 12;
\(\beta\)-restricted \(\gamma\)-marginal chain, \(k\ge 1\) sweeps, start \(\gamma_0\)):

\[
\boxed{
\bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\bigr\|_{TV}
\ \le\
\underbrace{(1-\varepsilon)^{\lfloor\alpha k\rfloor}}_{\text{minorization term}}
\;+\;
\underbrace{
\frac{U}{\alpha}\,
\frac{1+2b+\lambda V(\gamma_0)}{1+2b/(1-\lambda)}\,
\alpha^{k}
}_{\text{drift term}}.
}
\]

Here \(\pi_{\gamma\mid\widetilde B}=\pi_\gamma(\cdot\mid\widetilde B(\delta_2),y)\) and
\(P_\gamma=P_{\gamma\mid\widetilde B(\delta_2)}\) (§0).

**Free tuning parameter:** \(\alpha\in\bigl(\lambda(\delta_2),1\bigr)\) (minimize RHS at
target \(k\)).

---

#### Minorization constants

| Symbol | Definition |
|--------|------------|
| \(P_\gamma\) | \(P_{\gamma\mid\widetilde B(\delta_2)}\) — \(\gamma\)-marginal of Gibbs for \(\pi(\cdot\mid y,\beta\in\widetilde B(\delta_2))\) (§0) |
| \(\pi_{\gamma\mid\widetilde B}\) | \(\pi_\gamma(\cdot\mid\widetilde B(\delta_2),y)\) — invariant \(\gamma\)-target of \(P_\gamma\) |
| \(d\) | Deficiency level \(d=d(\delta)\) from `deficiency_calibrate(\delta, ...)` |
| \(\widetilde C_d\) | \(\{\gamma:\Psi(\gamma)\le d\}\), \(\Psi=\log(\varepsilon(\gamma^\star)/\varepsilon(\gamma))\) |
| \(\varepsilon(\gamma^\star)\) | C05 closure value \(\det(I+\tilde J)^{-1/2}\) at `population_mode()` |
| \(\varepsilon(d)\) | \(e^{-d}\,\varepsilon(\gamma^\star)\) — kernel floor on \(\widetilde C_d\) (Lemma 17(a)) |
| \(Q\) | Untruncated refresh \(N(\gamma^\star,P_{11}^{-1})\) |
| \(\varepsilon\) | \(\varepsilon(d)\) — Doeblin constant for \(P_\gamma\) on \(\widetilde C_d\) |
| \(Q(\widetilde C_d)\) | \(Q\)-mass on the certified set (chi-sq lower bound); used for **restricted**-chain \(\varepsilon_{\mathrm{rest}}(d)=\varepsilon(d)\,Q(\widetilde C_d)\) in `certificate()`, not in §3.1 |

---

#### Drift / Lyapunov constants

| Symbol | Definition |
|--------|------------|
| \(q\) | \(\dim(\gamma)\) |
| \(\gamma^\star\) | `population_mode()` fixed point |
| \(P_{11}\) | \(\Lambda_\gamma+\sum_j H_j^\top P_b H_j\) (anchor metric; see §8) |
| \(V(\gamma)\) | \(1+\tfrac12\|\gamma-\gamma^\star\|_{P_{11}}^2\) |
| \(V(\gamma_0)\) | \(1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\) |
| \(\widetilde B(\delta_2)\) | Compact \(\beta\)-safe set; \(\pi_\beta(\widetilde B^{\,c})\le\delta_2\) |
| \(\omega_{j,i}(\delta_2)\) | \(\inf_{\beta\in\widetilde B(\delta_2)} w_{j,i}(\beta)\) |
| \(\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)\) | \(D_j^\top\operatorname{diag}(\omega_{j,i}(\delta_2))\,D_j\) |
| \(P_{22,j}^{\mathrm{LB}}(\delta_2)\) | \(P_b+\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)\) |
| \(A^{\mathrm{LB}}(\delta_2)\) | \(P_{11}^{-1/2}\,S^{\mathrm{LB}}(\delta_2)\,P_{11}^{-1/2}\), \(S^{\mathrm{LB}}(\delta_2)=\sum_j H_j^\top P_b\,(P_{22,j}^{\mathrm{LB}}(\delta_2))^{-1}P_b H_j\) |
| \(\kappa_i^{\mathrm{LB}}(\delta_2)\) | \(\mathrm{eig}(A^{\mathrm{LB}}(\delta_2))\); floor coupling spectrum (§2) |
| \(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\) | \(\max_i\kappa_i^{\mathrm{LB}}(\delta_2)\) |
| \(\lambda(\delta_2)\) | \(\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\) |
| \(C_\beta^{+}(\delta_2)\) | \(\tfrac12\sum_i \kappa_i^{\mathrm{LB}}(\delta_2)=\tfrac12\,\mathrm{tr}(A^{\mathrm{LB}}(\delta_2))\) |
| \(b(\delta_2)\) | \(1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+C_\beta^{+}(\delta_2)\) |

---

#### Rosenthal auxiliary constants

| Symbol | Definition |
|--------|------------|
| \(U(\delta_2,d)\) | \(1+2b(\delta_2)+\lambda(\delta_2)\,V_{\sup}(d)\) |
| \(V_{\sup}(d)\) | \(\sup_{\gamma\in\widetilde C_d}V(\gamma)\); finite since \(\widetilde C_d\) compact. Bound \(V_{\sup}(d)\le 1+d+\sup_{\widetilde C_d}\bar\Phi\) |

**Drift small-set convention:** Foster (§1.2) uses \(C=\widetilde C_d\); the indicator
\(b\,\mathbb I_{C^c}\) adds \(b\) only when \(\gamma\notin\widetilde C_d\). Rosenthal
uses the same \(C\) for regeneration.

---

#### Expanded drift term (single line)

Define the **drift prefactor**

\[
\boxed{
D(\gamma_0)
:=
\frac{U}{\alpha}\,
\frac{1+2b+\lambda V(\gamma_0)}{1+2b/(1-\lambda)}
=
\frac{1+2b+\lambda V_{\sup}}{\alpha}\,
\frac{1+2b+\lambda V(\gamma_0)}{1+2b/(1-\lambda)}.
}
\]

Then

\[
\bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\bigr\|_{TV}
\le
(1-\varepsilon)^{\lfloor\alpha k\rfloor}+D(\gamma_0)\,\alpha^{k}.
\]

---

#### One-line substitution (all certificate constants)

**Minorization (depends on \(d\)).** Certified set \(\widetilde C_d=\{\Psi\le d\}\) with
\(d=d(\delta)\) from `deficiency_calibrate(\delta,\ldots)` at \(\beta\in\widetilde B(\delta_2)\):

\[
\varepsilon(d):=e^{-d}\,\varepsilon(\gamma^\star),
\qquad
V_{\sup}(d):=\sup_{\gamma\in\widetilde C_d}V(\gamma)
\ \le\ 1+d+\sup_{\gamma\in\widetilde C_d}\bar\Phi(\gamma)
\]

(\(\bar\Phi=V_0-\Psi\ge0\); on \(\partial\widetilde C_d=\{\Psi=d\}\), evaluate
\(V_{\sup}(d)\) numerically or via the envelope).

**Drift (depends on \(\delta_2\)).** Floor spectrum on \(\widetilde B(\delta_2)\):
\(\kappa_i^{\mathrm{LB}}(\delta_2)=\mathrm{eig}(A^{\mathrm{LB}}(\delta_2))\),
\(\kappa_{\max}^{\mathrm{LB}}(\delta_2):=\max_i\kappa_i^{\mathrm{LB}}(\delta_2)\).
Start: \(V(\gamma_0)=1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\).

\[
\boxed{
\begin{aligned}
\bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\bigr\|_{TV}
\ \le\
&\bigl(1-\varepsilon(d)\bigr)^{\lfloor\alpha k\rfloor} \\[4pt]
&\quad+
\frac{
1+2\Bigl[
1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
\Bigr]
+\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2 V_{\sup}(d)
}{\alpha} \\[4pt]
&\quad\quad\times
\frac{
1+2\Bigl[
1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
\Bigr]
+\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\Bigl(1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\Bigr)
}{
1+2\Bigl[
1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
\Bigr]\Big/\Bigl(1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\Bigr)
}\,
\alpha^{k},
\end{aligned}
}
\]

with \(\alpha\in\bigl(\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2,\,1\bigr)\).

**Monotonicity in \(d\).** Larger \(d\) enlarges \(\widetilde C_d=\{\Psi\le d\}\) (more
\(\gamma\) satisfy minorization). Accordingly:

\[
d_1<d_2
\quad\Longrightarrow\quad
\varepsilon(d_1)>\varepsilon(d_2),
\qquad
V_{\sup}(d_1)\le V_{\sup}(d_2),
\]

with \(\varepsilon(d)=e^{-d}\varepsilon(\gamma^\star)\) **decreasing** and
\(V_{\sup}(d)\) **increasing** in \(d\). Within §3.1 to
\(\pi_{\gamma\mid\widetilde B}\), smaller certified \(d\) is uniformly better
(subject to calibration feasibility).

**Limits in \(d\).** As \(d\downarrow 0\), \(\widetilde C_d\downarrow\{\gamma:\Psi(\gamma)\le0\}\),
which collapses to \(\{\gamma^\star\}\) (since \(\Psi\ge0\) with \(\Psi(\gamma^\star)=0\)).
Hence

\[
\varepsilon(d)\ \longrightarrow\ \varepsilon(\gamma^\star),
\qquad
V_{\sup}(d)\ \longrightarrow\ V(\gamma^\star)=1
\qquad(d\downarrow 0),
\]

i.e.\ the profile floor \(\varepsilon(\gamma^\star)\) and the Lyapunov baseline at the
mode. As \(d\uparrow\infty\), \(\widetilde C_d\uparrow\mathbb R^q\),
\(\varepsilon(d)\downarrow 0\), and typically \(V_{\sup}(d)\uparrow\infty\) (unless
\(\bar\Phi\) is uniformly bounded).

**Optimal start.** For fixed \((d,\delta_2,\alpha)\), the drift factor is minimized at
\(\gamma_0=\gamma^\star\), since \(V(\gamma_0)=1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\ge1\)
with equality only at \(\gamma^\star\).

**Sharpest displayed certificate** (combine \(\gamma_0=\gamma^\star\) with the
\(d\downarrow0\) limits \(\varepsilon(d)\to\varepsilon(\gamma^\star)\),
\(V_{\sup}(d)\to1\); same expanded \(\kappa_i^{\mathrm{LB}}(\delta_2)\) drift block):

\[
\boxed{
\begin{aligned}
\bigl\|P_\gamma^{k}(\gamma^\star,\cdot)-\pi_{\gamma\mid\widetilde B}\bigr\|_{TV}
\ \le\
&\bigl(1-\varepsilon(\gamma^\star)\bigr)^{\lfloor\alpha k\rfloor} \\[4pt]
&\quad+
\frac{
1+2\Bigl[
1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
\Bigr]
+\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2
}{\alpha} \\[4pt]
&\quad\quad\times
\frac{
1+2\Bigl[
1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
\Bigr]
+\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2
}{
1+2\Bigl[
1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
\Bigr]\Big/\Bigl(1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\Bigr)
}\,
\alpha^{k}.
\end{aligned}
}
\]

(At \(\gamma_0=\gamma^\star\) the two drift numerators coincide because
\(V(\gamma^\star)=V_{\sup}(d)=1\) in the \(d\downarrow0\) limit; for finite \(d>0\) use
the general box above with \(V_{\sup}(d)\) and
\(V(\gamma_0)=1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\).)

**Impact and limits of \(d\) and \(\delta_2\).**

§3.1 bounds \(\|P_\gamma^k-\pi_{\gamma\mid\widetilde B}\|_{TV}\) for the
**\(\beta\)-restricted** sampler (§0). The **full** marginal
\(\pi_\gamma(\cdot\mid y)\) needs the extra step
\(\|\pi_{\gamma\mid\widetilde B}-\pi_\gamma\|_{TV}\le\pi_\beta(\widetilde B(\delta_2)^{\,c})\le\delta_2\).

| Knob | What it controls | Typical effect |
|------|------------------|----------------|
| \(\delta\) (γ tail) | Calibration input for \(d=d(\delta)\) at \(\beta\in\widetilde B(\delta_2)\) | Larger \(\delta\) ⇒ larger \(d\) ⇒ **smaller** \(\varepsilon(d)\), **larger** \(V_{\sup}(d)\) |
| \(d\) | Size of \(\widetilde C_d\) and regeneration set \(C\) | **Smaller \(d\)** ⇒ **larger** \(\varepsilon(d)\), **lower** \(V_{\sup}(d)\) in §3.1 |
| \(\delta_2\) (β tail) | \(\pi_\beta(\widetilde B^{\,c})\le\delta_2\); \(\widetilde B(\delta_2)\) and \(\omega_{j,i}(\delta_2)\) | **Smaller \(\delta_2\)** ⇒ larger \(\widetilde B\) ⇒ **higher** \(\kappa_i^{\mathrm{LB}}(\delta_2)\) (**worse** drift); **tighter** gap to full \(\pi_\gamma\). **Larger \(\delta_2\)** ⇒ **better** \(\kappa^{\mathrm{LB}}\) but \(\pi_{\gamma\mid\widetilde B}\) farther from \(\pi_\gamma\) |

**Note on \(V_{\sup}(d)\).** This is a sup **on** \(\widetilde C_d\) (where minorization
holds), not on \(\widetilde C_d^{\,c}\). Off-set mass enters via \(V(\gamma_0)\) and
\(b(\delta_2)\) in the drift term, not via \(V_{\sup}(d)\).

**Calibration order (recommended).** Fix \(\delta_2\) and build \(\widetilde B(\delta_2)\);
compute \(\kappa_i^{\mathrm{LB}}(\delta_2)\); define \(P_\gamma\) and
\(\pi_{\gamma\mid\widetilde B}\); calibrate \(d(\delta)\) (hence \(\varepsilon(d)\),
\(V_{\sup}(d)\)) at \(\beta\in\widetilde B(\delta_2)\). To report against the **full**
\(\pi_\gamma\), add \(\pi_\beta(\widetilde B^{\,c})\le\delta_2\) (§0). C05-only
\(\gamma\)-truncation \(\pi_\gamma(\widetilde C_d^{\,c})\) remains in §4.

**Spectrum link (same matrix family as C05/C03).**

| Drift constant | Eigenvalue form (all depend on \(\delta_2\)) |
|----------------|---------------------------------------------|
| \(1-\lambda(\delta_2)\) | \(1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\) — MT shift |
| \(\lambda(\delta_2)\) | \(\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\) — mean contraction |
| \(C_\beta^{+}(\delta_2)\) | \(\tfrac12\sum_i \kappa_i^{\mathrm{LB}}(\delta_2)=\tfrac12\,\mathrm{tr}(A^{\mathrm{LB}}(\delta_2))\) |
| C05 B-weights | \(w_i=\kappa_i/(1-\kappa_i)\) at \((\gamma^\star,\beta)\) — deficiency sizing (**not** \(\kappa_i^{\mathrm{LB}}(\delta_2)\)) |

Compute \(\kappa_i^{\mathrm{LB}}(\delta_2)\) from `deficiency_spectrum()` with
\(P_{22,j}^{\mathrm{LB}}(\delta_2)\) built from \(\omega_{j,i}(\delta_2)\) on
\(\widetilde B(\delta_2)\).

Write \(b(\delta_2)\) numerically, then evaluate the bound (or optimize
\(\alpha\in\bigl(\lambda(\delta_2),1\bigr)\) at target \(k\)).

---

For \(d_V>2b/(1-\lambda)\), take drift sublevel \(C_V=\{\gamma:V(\gamma)\le d_V\}\) with
\(\widetilde C_d\subseteq C_V\). Define

\[
\alpha^{-1}:=\frac{1+2b+\lambda d_V}{1+d_V}<1,
\qquad
\Lambda:=1+2(\lambda d_V+b).
\]

Then for any \(r\in(0,1)\),

\[
P(T>k)\le
(1-\varepsilon)^{rk}+\alpha^{-k}(\alpha\Lambda)^{rk}\Bigl[1+\frac{b}{1-\lambda}+V(\gamma_0)\Bigr],
\]

and \(\bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\bigr\|_{TV}\le P(T>k)\)
coupling inequality (`LOGIT_SINGLE_GROUP_SAFE_REGION.md` §4.6).

---

#### C05-only special case (no drift hypothesis)

If only minorization on \(\widetilde C_d\) is used (Theorem 1, **restricted** chain
started in \(\widetilde C_d\); rate uses \(\varepsilon_d\) directly):

\[
\boxed{
\|P_\gamma^n(\gamma,\cdot\mid\widetilde C_d)-\pi_\gamma\|_{TV}
\le
(1-\varepsilon_d)^n+\delta,
\qquad
\delta=\pi_\gamma(\widetilde C_d^{\,c}).
}
\]

This is **not** the Rosenthal bound above; it has no \((\lambda,b)\) and no drift term.

**Default:** \(\alpha=(1+\lambda)/2\); minimize the full RHS over \(\alpha\in(\lambda,1)\) at
target \(k\).

### 3.2 Reading the two terms

| Term | Mechanism | Certificate input |
|------|-----------|-------------------|
| \((1-\varepsilon)^{\lfloor\alpha k\rfloor}\) | **Minorization** on \(\widetilde C_d\) | \(\varepsilon(d)=e^{-d}\varepsilon(\gamma^\star)\) for \(P_\gamma\) on \(\widetilde B(\delta_2)\) |
| \(\frac{U}{\alpha}(\cdots)\alpha^k\) | **Drift** toward \(\widetilde C_d\) | \(\lambda(\delta_2)\), \(C_\beta^{+}(\delta_2)\), \(b(\delta_2)\), \(V_{\sup}(d)\) from \(\kappa_i^{\mathrm{LB}}(\delta_2)\) |

Both \(\varepsilon\) and \(\lambda\) enter. Weakening either slows convergence.

### 3.3 Coupling-time form (equivalent)

Rosenthal's proof also gives, for coupling time \(T\) to regeneration on
\(\widetilde C_d\), and any \(r\in(0,1)\),

\[
P(T>k)\ \le\
(1-\varepsilon)^{rk}
\;+\;
\alpha^{-k}(\alpha\Lambda)^{rk}
\Bigl[1+\tfrac{b}{1-\lambda}+V(\gamma_0)\Bigr],
\]

with explicit \(\alpha^{-1}<1\) when \(C=\{V\le d_{\mathrm{lev}}\}\) and
\(d_{\mathrm{lev}}>2b/(1-\lambda)\). TV follows from \(P(T>k)\) via the standard
coupling inequality. See `LOGIT_SINGLE_GROUP_SAFE_REGION.md` §4.6.

---

## 4. C05 Theorem 1 (minorization + truncation only)

When only **minorization on \(\widetilde C_d\)** is invoked (no drift hypothesis),
`restricted_gibbs_minorization _v4.md` Theorem 1 gives the simpler **two-term**
bound for the restricted chain started in \(\widetilde C_d\):

\[
\boxed{
\|P_\gamma^n(\gamma,\cdot\mid\widetilde C_d)-\pi_\gamma\|_{TV}
\ \le\
(1-\varepsilon_d)^n+\delta,
\qquad
\delta=\pi_\gamma(\widetilde C_d^{\,c}).
}
\]

| Term | Role |
|------|------|
| \((1-\varepsilon_d)^n\) | Doeblin refresh on the **restricted** chain in \(\widetilde C_d\) |
| \(\delta\) | **Static** truncation: \(\|\pi(\cdot\mid\widetilde C_d)-\pi\|_{TV}\) |

This is **not** the Foster plateau \(b/(1-\lambda)\). Rosenthal §3 adds the drift
term when \(\gamma_0\notin\widetilde C_d\) or when dynamic escape must be controlled.

**Relationship.** C05 supplies \((\varepsilon(d),\delta_\gamma)\) for \(\gamma\)-only
truncation (§4). §2 supplies \((\lambda(\delta_2),b(\delta_2))\) on
\(\widetilde B(\delta_2)\). Rosenthal §3.1 bounds
\(\|P_\gamma^k-\pi_{\gamma\mid\widetilde B}\|_{TV}\). The gap
\(\|\pi_{\gamma\mid\widetilde B}-\pi_\gamma\|_{TV}\le\delta_2\) is the **\(\beta\)-conditioning**
step (§0).

---

## 5. Meyn–Tweedie simplified bound (optional)

For qualitative geometric ergodicity without optimizing Rosenthal constants:

\[
\bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\bigr\|_{TV}
\ \le\
M\,V(\gamma_0)\,\rho^k,
\qquad
\rho=\max\Bigl((1-\varepsilon)^{1/m_0},\,\alpha\Bigr),
\]

with block length \(m_0\) and \(\alpha\in(\lambda,1)\). The prefactor \(M\) is
typically looser than §3.1 at a fixed \(k\). See
`GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md` §4.

---

## 6. Application pipeline (draft)

Work on \(\widetilde{\mathcal R}=\widetilde C_d\times\widetilde B(\delta_2)\).

### Step 1 — Build \(\widetilde B(\delta_2)\)

See **`inst/BETA_MARGINAL_MODE_LEVELSET.md`**.

1. Integrate \(\gamma\) → \(\Lambda_\beta\), \(\mu_\beta\); find **marginal mode**
   \(\beta^\dagger\) (Newton on \(f(\beta)\)).
2. Set \(r=r_{\mathrm{Gauss}}(Jp_{\mathrm{re}},\delta_2)\);
   \(\widetilde B(\delta_2)=\{\Xi(\beta)\le r\}\) with \(\Xi=f-f(\beta^\dagger)\).
3. Support functions on \(\widetilde B(\delta_2)\) → \(\omega_{j,i}(\delta_2)>0\);
   build \(P_{22,j}^{\mathrm{LB}}(\delta_2)\);
   \(\kappa_i^{\mathrm{LB}}(\delta_2)=\mathrm{eig}(A^{\mathrm{LB}}(\delta_2))\).
4. **Tail for full \(\pi_\gamma\):** report \(\delta_2\) as asymptotic Laplace mass;
   optional MC/bracket for true \(\pi_\beta(\widetilde B^c)\) (BETA note §7 Phase 4).

### Step 2 — Drift constants (functions of \(\delta_2\))

\[
\lambda(\delta_2)=\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2,
\qquad
C_\beta^{+}(\delta_2)=\tfrac12\sum_{i=1}^q \kappa_i^{\mathrm{LB}}(\delta_2),
\qquad
b(\delta_2)=1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+C_\beta^{+}(\delta_2).
\]

### Step 3 — Minorization constants (C05)

At \(\beta\in\widetilde B(\delta_2)\), `deficiency_calibrate(\delta,\ldots)` yields
\(d=d(\delta)\), \(\varepsilon(d)=e^{-d}\varepsilon(\gamma^\star)\), \(V_{\sup}(d)\), and

\[
\varepsilon=\varepsilon(d)
\quad\text{(Rosenthal §3.1)},
\qquad
\varepsilon_{\mathrm{rest}}(d)=\varepsilon(d)\,Q(\widetilde C_d)
\quad\text{(`certificate()`, restricted chain)},
\qquad
\delta=\pi_\gamma(\widetilde C_d^{\,c}).
\]

Smaller \(d\) shrinks \(\widetilde C_d\), **increases** \(\varepsilon(d)\), and
**decreases** \(V_{\sup}(d)\) in the bound to \(\pi_{\gamma\mid\widetilde B}\).

### Step 4 — Rosenthal TV bound

Plug \(\bigl(\lambda(\delta_2),b(\delta_2),\varepsilon(d),V(\gamma_0)\bigr)\) into §3.1
for \(\|P_\gamma^k(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\|_{TV}\). Optional:
add \(\delta_2\) for \(\pi_\gamma(\cdot\mid y)\) (§0).

### Step 5 — Report (template)

| Quantity | Formula | Source |
|----------|---------|--------|
| \(\kappa_i^{\mathrm{LB}}(\delta_2)\) | \(\mathrm{eig}(A^{\mathrm{LB}}(\delta_2))\) | Step 1 on \(\widetilde B(\delta_2)\) |
| \(\lambda(\delta_2)\) | \(\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\) | Floor spectrum §2 |
| \(C_\beta^{+}(\delta_2)\) | \(\tfrac12\sum_i \kappa_i^{\mathrm{LB}}(\delta_2)\) | Same spectrum |
| \(b(\delta_2)\) | \(1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+q/2+C_\beta^{+}(\delta_2)\) | §2.4 |
| \(\varepsilon(d)\) | \(e^{-d}\varepsilon(\gamma^\star)\) | Lemma 17(a); Rosenthal §3.1 |
| \(V_{\sup}(d)\) | \(\sup_{\gamma\in\widetilde C_d}V(\gamma)\) | §3.1; **increases** with \(d\) (sup on a larger set) |
| \(\varepsilon_{\mathrm{rest}}(d)\) | \(\varepsilon(d)\,Q(\widetilde C_d)\) | `certificate()` (restricted chain) |
| \(\delta_2\) | \(\pi_\beta(\widetilde B^{\,c})\approx\delta_2\) (asymptotic; BETA note) | Gap \(\pi_{\gamma\mid\widetilde B}\to\pi_\gamma\) (§0) |
| \(\delta_\gamma\) | \(\pi_\gamma(\widetilde C_d^{\,c})\) | C05-only §4 (doubly restricted) |
| TV bound (restricted) | §3.1 | \(\|P_\gamma^k-\pi_{\gamma\mid\widetilde B}\|_{TV}\) |
| TV to full \(\pi_\gamma\) | §0 + §3.1 | §3.1 bound \(+\,\delta_2\) |

---

## 7. Gaussian exact case (Chapter C03)

When the target is multivariate normal and \(\lambda^\star=\kappa_{\max}(A)<1\)
globally, Claim 2 and Theorem 3 give explicit decay without Foster–Lyapunov; the
matrix \(A\) supplies \(\lambda^\star\) directly. The present note is the
**non-Gaussian / restricted** route where \(\widetilde B(\delta_2)\) replaces global
\(\lambda^\star<1\).

---

## 8. Open items

1. **Formal lemma** for §2.2–2.4 on \(\widetilde B(\delta_2)\) (Foster on \(C^c\)).
2. **Asymptotic lemma:** BvM + \(r_{\mathrm{Gauss}}\) on \(\{\Xi\le r\}\) ⇒
   \(\pi_\beta(\widetilde B^c)=\delta_2+O(N^{-1/2})\) (`BETA_MARGINAL_MODE_LEVELSET.md` §4).
3. **Metric:** \(P_{11}\) fixed at \((\gamma^\star,\beta^\star)\) vs worst-case on
   \(\widetilde B\).
4. **Joint TV:** add \(\delta_2\), \(\delta_{12}\) from
   `JOINT_GAMMA_BETA_TV_CERTIFICATE.md` §6.
5. **Implementation** (`BETA_MARGINAL_MODE_LEVELSET.md` §7): **`beta_marginal_safe_set()`**,
   **`floor_coupling_spectrum()`**, **`rosenthal_tv_bound()`**, **`gamma_beta_tv_certificate()`**
   (sharpest displayed box §691–722); `group_precision_floor()` gains `mode_method` /
   `level_method`; tail MC check in `data-raw/` remains optional.

---

## 9. References

- `inst/restricted_gibbs_minorization _v4.md` — C05 minorization, Theorem 1
- `inst/GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md` — general MT/Rosenthal catalog
- `inst/JOINT_GAMMA_BETA_TV_CERTIFICATE.md` — \(\widetilde B(\delta_2)\), product region
- `inst/BETA_MARGINAL_MODE_LEVELSET.md` — marginal \(\beta^\dagger\), \(r_{\mathrm{Gauss}}\), asymptotics
- `inst/LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` — scalar worked Foster proof
- `inst/CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` — matrix \(A\), \(\kappa_{\max}\)
- Rosenthal (1995); Meyn & Tweedie (1993)
