# \(\gamma\)-marginal drift, minorization, and Rosenthal TV bound

Companion to `restricted_gibbs_minorization _v4.md` (Theorem 1),
`GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md`, `JOINT_GAMMA_BETA_TV_CERTIFICATE.md`,
and `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md`.

**Status:** draft design note — constants are structurally explicit. The formal
existence result for \(\widetilde B(\delta_2)\) is **Lemma B-Cert** (§0.3); what
remains asymptotic is the identification of the design level
\(r_{\mathrm{Gauss}}\) with the exact level \(r^\star(\delta_2)\) (Remark B.3).

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
- **Existence:** an *exact* level \(r^\star(\delta_2)\) with
  \(\pi_\beta(\widetilde B^{\,c})=\delta_2\) always exists and is unique — **Lemma B-Cert**
  (§0.3, part (iii)). \(r_{\mathrm{Gauss}}\) is the computable surrogate for it
  (Remark B.3).

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

### 0.3 Lemma B-Cert: existence of a certified \(\beta\)-safe set

This is the \(\beta\)-side counterpart of **Theorem 2** of
`restricted_gibbs_minorization _v4.md` (§3.4, proved in §6 / A.4), which asserts the
existence of a certified \(\gamma\)-set \(\widetilde C_d\) carrying an attained
minorization constant. Here the object is \(\widetilde B(\delta_2)\), and the certified
quantity is the **weight floor**, which is what makes the Foster constants
\((\lambda(\delta_2),b(\delta_2))\) of §2 well defined and strictly sub-unit.

**Notation for this subsection.** Stack the hyper-design and the random-effect prior
precision as \(\mathcal H:=[H_1^\top,\ldots,H_J^\top]^\top\),
\(\mathcal P_b:=I_J\otimes P_b\), so \(P_{11}^{\mathrm{RE}}=\mathcal H^\top\mathcal P_b\mathcal H\)
and \(P_{11}=\Lambda_\gamma+P_{11}^{\mathrm{RE}}\). Integrating \(\gamma\) out of the
prior gives the \(\gamma\)-marginal \(\beta\)-target of §0.2,

\[
\widetilde\pi(\beta\mid y)\;\propto\;
\exp\Bigl(\sum_j\ell_j(\beta_j)\Bigr)\,
\exp\Bigl(-\tfrac12\|\beta-\mu_\beta\|_{\Lambda_\beta}^2\Bigr),
\qquad
\mu_\beta=\mathcal H\mu_0,
\]

with \(f(\beta)=-\log\widetilde\pi(\beta\mid y)+\mathrm{const}\),
\(\Xi=f-f(\beta^\dagger)\), \(\widetilde B(r)=\{\Xi\le r\}\), and \(n=Jp_{\mathrm{re}}\).
Write \(\eta_{j,i}(\beta_j)=d_{j,i}^\top\beta_j+o_{j,i}\) for the group linear
predictors (\(d_{j,i}^\top\) the rows of \(D_j\)), \(w_{j,i}(\eta)\) for the IRLS
weights, and

\[
\mathcal I(\beta):=\operatorname{blockdiag}_j\bigl(D_j^\top W_j(\beta_j)D_j\bigr),
\qquad
W_j(\beta_j)=\operatorname{diag}\bigl(w_{j,i}(\eta_{j,i}(\beta_j))\bigr).
\]

---

> **Lemma B-Cert (existence of a certified \(\beta\)-safe set).**
>
> **Assumptions (model only).** These are the \(\beta\)-side reading of (H1)–(H3) in
> **Proposition R-Cert** (§3.1):
>
> 1. **(B1 — proper random-effect prior)** \(P_b=\Psi^{-1}\succ0\), and either
>    \(\Lambda_\gamma\succ0\), or \(\Lambda_\gamma=0\) together with (B3b).
> 2. **(B2 — log-concave group likelihoods)** Each \(\ell_j\) is \(C^2\) and concave in
>    \(\beta_j\), with IRLS weights \(w_{j,i}(\cdot)\) continuous and **strictly positive
>    at every finite \(\eta\)** (true for every canonical GLM family in the package:
>    \(\mu(1-\mu)>0\) for logit/probit, \(\mu=e^\eta>0\) for Poisson-log).
> 3. **(B3 — rank / estimability)** **(B3a)** each \(D_j\) has full column rank
>    \(p_{\mathrm{re}}\) and the stacked hyper-design has rank \(q\), i.e.
>    \(P_{11}^{\mathrm{RE}}\succ0\); **(B3b)** a finite group-wise MLE exists for each
>    \(j\) (no complete separation) — needed only when \(\Lambda_\gamma=0\).
>
> **Statement.** Under (B1)–(B3):
>
> **(i) Anchor and profile.** \(\Lambda_\beta\succ0\) when \(\Lambda_\gamma\succ0\), and
> \(f\) is finite, \(C^2\), **strictly convex and coercive** on \(\mathbb R^n\) in both
> prior regimes. Hence \(\beta^\dagger=\arg\min f\) exists and is **unique**, and
> \[
> \Xi(\beta)=f(\beta)-f(\beta^\dagger)\ \ge\ 0,
> \qquad
> \Xi(\beta)=0\iff\beta=\beta^\dagger,
> \]
> with \(\Xi\) convex, finite everywhere, and \(\nabla\Xi(\beta^\dagger)=0\).
>
> **(ii) Geometry.** For every \(r>0\), \(\widetilde B(r)=\{\Xi\le r\}\) is **convex,
> compact**, with \(\beta^\dagger\) in its interior — under (B1)–(B3) alone, in both
> prior regimes. If in addition \(\Lambda_\gamma\succ0\), the **outer** ellipsoid below
> holds; if further the weights are bounded above, \(w_{j,i}\le\bar w_{j,i}\) (automatic
> for binomial: \(\bar w_{j,i}=n_{j,i}/4\)), then with
> \(\Lambda^{\mathrm{UB}}:=\Lambda_\beta+\operatorname{blockdiag}_j(D_j^\top\bar W_jD_j)\)
> the **inner** ellipsoid holds too:
> \[
> \boxed{\;
> \bigl\{\|\beta-\beta^\dagger\|^2_{\Lambda^{\mathrm{UB}}}\le2r\bigr\}
> \ \subseteq\
> \widetilde B(r)
> \ \subseteq\
> \bigl\{\|\beta-\beta^\dagger\|^2_{\Lambda_\beta}\le2r\bigr\}.\;}
> \]
>
> **(iii) Exhaustion and exact mass calibration.** \(r\mapsto\widetilde B(r)\) is
> nondecreasing with \(\bigcup_{r>0}\widetilde B(r)=\mathbb R^n\); the escape function
> \(G(r):=\pi_\beta\bigl(\widetilde B(r)^{\,c}\bigr)\) is **continuous and strictly
> decreasing** on \([0,\infty)\), with \(G(0)=1\), \(G>0\) everywhere, and
> \(G(r)\downarrow0\). Hence for every
> \(\delta_2\in(0,1)\) there is a **unique** level \(r^\star(\delta_2)<\infty\) with
> \[
> \boxed{\;
> \pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr)=\delta_2
> \quad\text{exactly, at }\ \widetilde B(\delta_2):=\widetilde B\bigl(r^\star(\delta_2)\bigr),\;}
> \]
> and \(r^\star(\delta_2)\le r_{\mathrm{P2}}(n,\delta_2)<\infty\), the finite level
> certified by the distribution-free radial budget (P2)
> (`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` §3A.4.1).
>
> **(iv) Attained floors; non-degenerate drift constants.** For every \((j,i)\) the
> weight floor \(\omega_{j,i}(\delta_2)=\inf_{\widetilde B(\delta_2)}w_{j,i}\) is
> **attained and strictly positive**, hence
> \(\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)\succ0\),
> \(P_{22,j}^{\mathrm{LB}}(\delta_2)\succ P_b\), and
> \[
> \boxed{\;
> S^{\mathrm{LB}}(\delta_2)\ \prec\ P_{11}^{\mathrm{RE}}\ \preceq\ P_{11}
> \quad\Longrightarrow\quad
> \kappa_{\max}^{\mathrm{LB}}(\delta_2)<1 .\;}
> \]
> Consequently the Foster pair of §2 is well defined and non-degenerate:
> \(\lambda(\delta_2)=\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2<1\),
> \(b(\delta_2)<\infty\), and the deficiency weight
> \(w_{\max}^{\mathrm{LB}}=\max_i\kappa_i^{\mathrm{LB}}/(1-\kappa_i^{\mathrm{LB}})<\infty\).
>
> **(v) Monotonicity (the \(\delta_2\) trade-off).** \(\delta_2\downarrow\) implies
> \(r^\star(\delta_2)\uparrow\), \(\widetilde B(\delta_2)\uparrow\),
> \(\omega_{j,i}(\delta_2)\downarrow\), \(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\uparrow\)
> and \(\lambda(\delta_2)\uparrow\): a smaller \(\beta\)-truncation term is paid for with
> a slower certified drift rate, never the reverse.
>
> In particular (i)–(iii) **discharge (H4)** of **Proposition R-Cert**: the marginal
> \(\beta\)-mode, the convex profile with coercive level sets, and the mass-calibrated
> set are *consequences* of (B1)–(B3), not additional hypotheses.

---

#### Proof of Lemma B-Cert

**Step 1 (the marginal precision \(\Lambda_\beta\)).** Under the prior
\(\beta\mid\gamma\sim N(\mathcal H\gamma,\mathcal P_b^{-1})\),
\(\gamma\sim N(\mu_0,\Lambda_\gamma^{-1})\), integrating \(\gamma\) gives
\(\beta\sim N\bigl(\mathcal H\mu_0,\ \Sigma_\beta\bigr)\) with
\(\Sigma_\beta=\mathcal P_b^{-1}+\mathcal H\Lambda_\gamma^{-1}\mathcal H^\top\). Since
\(\mathcal P_b^{-1}\succ0\) and the second term is positive semidefinite,
\(\Sigma_\beta\succ0\), so \(\Lambda_\beta=\Sigma_\beta^{-1}\succ0\). Woodbury gives the
form used throughout,

\[
\Lambda_\beta=\mathcal P_b-\mathcal P_b\mathcal H\,P_{11}^{-1}\mathcal H^\top\mathcal P_b ,
\qquad P_{11}=\Lambda_\gamma+\mathcal H^\top\mathcal P_b\mathcal H .
\]

At \(\Lambda_\gamma=0\) the same expression reads
\(\Lambda_\beta=\mathcal P_b^{1/2}(I-\Pi)\mathcal P_b^{1/2}\succeq0\), where \(\Pi\) is the
orthogonal projector onto \(\operatorname{range}(\mathcal P_b^{1/2}\mathcal H)\); hence
\(\Lambda_\beta\) is singular with \(\ker\Lambda_\beta=\operatorname{range}(\mathcal H)\)
— the population-shift directions carry no prior information in the flat limit. (This
is the \(\beta\)-side analogue of v4 Lemma 7: the prior degenerates exactly along the
hyper-design range, and (B3) must replace it there.) \(\square\)

**Step 2 (strict convexity of \(f\)).** \(f(\beta)=\sum_j(-\ell_j(\beta_j))
+\tfrac12\|\beta-\mu_\beta\|^2_{\Lambda_\beta}\) is \(C^2\) with

\[
\nabla^2f(\beta)=\Lambda_\beta+\mathcal I(\beta),
\qquad
\mathcal I(\beta)\succeq0 \ \text{ by (B2)} .
\]

*Proper population prior.* \(\nabla^2f\succeq\Lambda_\beta\succ0\) (Step 1) — done.

*Flat limit.* Let \(0\ne v\in\mathbb R^n\). If \(v\notin\ker\Lambda_\beta\) then
\(v^\top\nabla^2f\,v\ge v^\top\Lambda_\beta v>0\). If \(v\in\ker\Lambda_\beta\), Step 1
gives \(v=\mathcal Hu\) with \(u\ne0\) (\(\mathcal H\) is injective because
\(\mathcal H^\top\mathcal P_b\mathcal H=P_{11}^{\mathrm{RE}}\succ0\) by (B3a)). By (B3a)
there is \(j_0\) with \(H_{j_0}u\ne0\), and \(D_{j_0}\) has full column rank, so
\(D_{j_0}H_{j_0}u\ne0\); with \(W_{j_0}\succ0\) by (B2),

\[
v^\top\nabla^2f(\beta)\,v
\ \ge\ v^\top\mathcal I(\beta)v
\ \ge\ (D_{j_0}H_{j_0}u)^\top W_{j_0}(\beta_{j_0})\,(D_{j_0}H_{j_0}u)\ >\ 0 .
\]

So \(\nabla^2f(\beta)\succ0\) at **every** finite \(\beta\), in both regimes — strictly,
not almost everywhere. (Same mechanism as v4 Lemma 8: (B3a) supplies a group that the
direction actually moves, (B2) supplies positive curvature there.) \(\square\)

**Step 3 (coercivity; existence and uniqueness of \(\beta^\dagger\)).**

*Proper population prior.* From \(\nabla^2f\succeq\Lambda_\beta\succ0\),
\(f(\beta)\ge f(\beta_0)+\nabla f(\beta_0)^\top(\beta-\beta_0)
+\tfrac12\|\beta-\beta_0\|^2_{\Lambda_\beta}\), which \(\to\infty\) as
\(\|\beta\|\to\infty\).

*Flat limit.* \(f\) is convex, so it is coercive iff its recession function is strictly
positive in every nonzero direction, i.e. iff \(t\mapsto f(\beta_0+tv)\to\infty\) for
every \(v\ne0\). For \(v\notin\ker\Lambda_\beta\) the quadratic term grows
quadratically while \(\sum_j(-\ell_j)\) is bounded below (each \(\ell_j\) is concave and,
by (B3b), attains a finite maximum). For \(v=\mathcal Hu\in\ker\Lambda_\beta\), pick
\(j_0\) with \(H_{j_0}u\ne0\) as in Step 2; (B3b) plus full column rank of \(D_{j_0}\)
makes \(-\ell_{j_0}\) coercive on \(\mathbb R^{p_{\mathrm{re}}}\) (a concave GLM
log-likelihood with a finite maximiser has no direction of recession), so the \(j_0\)
term \(\to\infty\) while the remaining terms are bounded below and the quadratic term is
constant along \(v\). Either way \(f(\beta_0+tv)\to\infty\).

A continuous, strictly convex (Step 2), coercive function attains its minimum at a
unique point \(\beta^\dagger\), where \(\nabla f(\beta^\dagger)=0\). Setting
\(\Xi:=f-f(\beta^\dagger)\) gives \(\Xi\ge0\) with equality only at \(\beta^\dagger\),
convex and finite everywhere. This proves **(i)**. \(\square\)

**Step 4 (geometry of \(\widetilde B(r)\)).** Convexity: \(\widetilde B(r)\) is a
sublevel set of the convex \(\Xi\). Closedness: \(\Xi\) is continuous. Boundedness:
\(\Xi\) is coercive (Step 3), so every sublevel set is bounded; hence
\(\widetilde B(r)\) is compact. Interior: \(\Xi\) is continuous with
\(\Xi(\beta^\dagger)=0<r\), so \(\widetilde B(r)\) contains a ball around
\(\beta^\dagger\) and has nonempty interior.

*Ellipsoid sandwich.* Since \(\Xi(\beta^\dagger)=0\) and \(\nabla\Xi(\beta^\dagger)=0\),
Taylor with integral remainder along the segment \(\beta_t=\beta^\dagger+t(\beta-\beta^\dagger)\)
gives, with \(h:=\beta-\beta^\dagger\),

\[
\Xi(\beta)=\int_0^1(1-t)\,h^\top\nabla^2\Xi(\beta_t)\,h\,dt .
\]

When \(\Lambda_\gamma\succ0\), \(\nabla^2\Xi\succeq\Lambda_\beta\) (Step 2), and
\(\int_0^1(1-t)\,dt=\tfrac12\), so
\(\Xi(\beta)\ge\tfrac12\|h\|^2_{\Lambda_\beta}\), i.e.
\(\widetilde B(r)\subseteq\{\|\beta-\beta^\dagger\|^2_{\Lambda_\beta}\le2r\}\) — the outer
ellipsoid. If moreover \(w_{j,i}\le\bar w_{j,i}\), then
\(\mathcal I(\beta)\preceq\operatorname{blockdiag}_j(D_j^\top\bar W_jD_j)\) uniformly in
\(\beta\), so \(\nabla^2\Xi\preceq\Lambda^{\mathrm{UB}}\) and the same identity gives
\(\Xi(\beta)\le\tfrac12\|h\|^2_{\Lambda^{\mathrm{UB}}}\), i.e.
\(\{\|\beta-\beta^\dagger\|^2_{\Lambda^{\mathrm{UB}}}\le2r\}\subseteq\widetilde B(r)\) —
the inner ellipsoid. Together these are the boxed sandwich, proving **(ii)**. \(\square\)

**Step 5 (exhaustion; continuity and strict monotonicity of \(G\)).** \(\Xi(\beta)<\infty\)
for every \(\beta\) (Step 3), so \(\beta\in\widetilde B(r)\) as soon as
\(r\ge\Xi(\beta)\): the family is nondecreasing with union \(\mathbb R^n\). Under
(B1)–(B3) the posterior is proper — automatic when \(\Lambda_\gamma\succ0\), and v4
Lemma 2 at \(\Lambda_\gamma=0\) — so \(\pi_\beta\) is a probability measure with a
density with respect to Lebesgue measure, and continuity from above on the nonincreasing
sets \(\widetilde B(r)^{\,c}\) with empty intersection gives \(G(r)\downarrow0\).

*No atoms.* \(G(r^-)-G(r)=\pi_\beta(\Xi=r)\), and \(\{\Xi=r\}\subseteq\partial\widetilde B(r)\)
is the boundary of a compact convex body with nonempty interior (Step 4), hence Lebesgue-null;
since \(\pi_\beta\ll\) Lebesgue, \(\pi_\beta(\Xi=r)=0\). So \(G\) is continuous on
\((0,\infty)\), and \(G(0)=1-\pi_\beta(\{\beta^\dagger\})=1\) for the same reason.

*Strict decrease.* Suppose \(0\le r_1<r_2\) with \(G(r_1)=G(r_2)\), i.e.
\(\pi_\beta(r_1<\Xi\le r_2)=0\). Fix any unit direction \(u\); \(s\mapsto\Xi(\beta^\dagger+su)\)
is continuous, equals \(0\) at \(s=0\) and \(\to\infty\) (coercivity, Step 3), so by the
intermediate value theorem it takes every value in \((r_1,r_2)\). Hence
\(\{r_1<\Xi<r_2\}\) is nonempty, and it is open because \(\Xi\) is continuous. Since
\(f<\infty\) everywhere, \(\widetilde\pi(\cdot\mid y)>0\) everywhere, so that open set
carries strictly positive \(\pi_\beta\)-mass — a contradiction. The same argument with
\(r_2=\infty\) shows \(G(r)>0\) for every \(r\).

Hence \(G:[0,\infty)\to(0,1]\) is continuous, strictly decreasing, equal to \(1\) at
\(0\), and tends to \(0\). By the intermediate value theorem, for every
\(\delta_2\in(0,1)\) there is a unique \(r^\star(\delta_2)\) with
\(G(r^\star(\delta_2))=\delta_2\). \(\square\)

**Step 6 (a finite certified level, without asymptotics).** \(\Xi\) is convex on
\(\mathbb R^n\) with a unique minimum \(\Xi(\beta^\dagger)=0\) (Step 3), which is exactly
the hypothesis of the radial budget (P2)
(`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` §3A.4.1, Proposition 2):

\[
\pi_\beta\bigl(\widetilde B(r)^{\,c}\bigr)
\ \le\
\frac{\Gamma(n,r)}{\gamma(n,r)}
\ \xrightarrow[r\to\infty]{}\ 0 .
\]

Let \(r_{\mathrm{P2}}(n,\delta_2)\) solve \(\Gamma(n,r)/\gamma(n,r)=\delta_2\). Then
\(G(r_{\mathrm{P2}})\le\delta_2\), so by monotonicity of \(G\),
\(r^\star(\delta_2)\le r_{\mathrm{P2}}(n,\delta_2)<\infty\). This proves **(iii)**, and
gives a *proved* finite level with no appeal to Laplace or BvM — the Gaussian design
level \(r_{\mathrm{Gauss}}(n,\delta_2)=\tfrac12\chi^2_{n,1-\delta_2}\) of §0.2 is the
asymptotically tight surrogate for \(r^\star(\delta_2)\), exact under the Laplace
reference and \(\delta_2+O(N^{-1/2})\) under BvM, but is not needed for existence. \(\square\)

**Step 7 (attained floors; \(\kappa_{\max}^{\mathrm{LB}}(\delta_2)<1\)).** Fix
\(\delta_2\) and write \(\widetilde B:=\widetilde B(\delta_2)\), compact by Step 4. Each
\(\eta_{j,i}\) is a continuous (affine) function of \(\beta\), so
\(\eta_{j,i}(\widetilde B)\) is a compact interval; \(w_{j,i}\) is continuous and
strictly positive on it by (B2), so by Weierstrass the infimum

\[
\omega_{j,i}(\delta_2)=\inf_{\beta\in\widetilde B}w_{j,i}\bigl(\eta_{j,i}(\beta)\bigr)
\]

is **attained** and \(>0\). Hence
\(\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)=D_j^\top\operatorname{diag}(\omega_{j,i}(\delta_2))D_j\succ0\)
by (B3a), and \(P_{22,j}^{\mathrm{LB}}(\delta_2)=P_b+\underline{\mathcal P}_{j,\mathrm{data}}(\delta_2)\succ P_b\).
Inverting reverses the order, \(\bigl(P_{22,j}^{\mathrm{LB}}\bigr)^{-1}\prec P_b^{-1}\).
For \(0\ne u\in\mathbb R^q\),

\[
u^\top\bigl(P_{11}^{\mathrm{RE}}-S^{\mathrm{LB}}(\delta_2)\bigr)u
=\sum_{j}(P_bH_ju)^\top
\Bigl(P_b^{-1}-\bigl(P_{22,j}^{\mathrm{LB}}(\delta_2)\bigr)^{-1}\Bigr)
(P_bH_ju)\ >\ 0,
\]

because every summand is nonnegative and, by (B3a), some \(H_{j_0}u\ne0\) makes the
\(j_0\) summand strictly positive. So
\(S^{\mathrm{LB}}(\delta_2)\prec P_{11}^{\mathrm{RE}}\preceq P_{11}\), hence
\(A^{\mathrm{LB}}(\delta_2)=P_{11}^{-1/2}S^{\mathrm{LB}}(\delta_2)P_{11}^{-1/2}\prec I\)
and \(\kappa_{\max}^{\mathrm{LB}}(\delta_2)<1\). Therefore
\(\lambda(\delta_2)=\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2<1\),
\(b(\delta_2)=1-\lambda(\delta_2)+\tfrac q2+\tfrac12\sum_i\kappa_i^{\mathrm{LB}}(\delta_2)<\infty\),
and \(w_{\max}^{\mathrm{LB}}=\max_i\kappa_i^{\mathrm{LB}}/(1-\kappa_i^{\mathrm{LB}})<\infty\).
This proves **(iv)**. \(\square\)

**Step 8 (monotonicity).** \(G\) is continuous and strictly decreasing (Step 5), so its
inverse \(\delta_2\mapsto r^\star(\delta_2)\) is strictly decreasing:
\(\delta_2\downarrow\) forces \(r^\star(\delta_2)\uparrow\) and hence
\(\widetilde B(\delta_2)\uparrow\) by nesting. An infimum over a larger set is no
larger, so \(\omega_{j,i}(\delta_2)\downarrow\); then
\(\underline{\mathcal P}_{j,\mathrm{data}}\downarrow\),
\(P_{22,j}^{\mathrm{LB}}\downarrow\), \(\bigl(P_{22,j}^{\mathrm{LB}}\bigr)^{-1}\uparrow\),
\(S^{\mathrm{LB}}\uparrow\), and \(A^{\mathrm{LB}}\uparrow\) in the positive semidefinite
order, so \(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\uparrow\) by Weyl monotonicity and
\(\lambda(\delta_2)\uparrow\). All the inequalities of (iv) remain strict at every
\(\delta_2\in(0,1)\), because Step 7 used only compactness of \(\widetilde B(\delta_2)\)
together with (B2) and (B3a), which hold at every level. This proves **(v)**, and
completes the proof. \(\blacksquare\)

---

**Remark B.1 (what each assumption buys).** (B1) supplies \(\Lambda_\beta\succ0\) and
with it strict convexity and the explicit outer ellipsoid; (B2) supplies both the
convexity of \(-\ell_j\) and the *strictly positive* weights that make the floors
non-vacuous; (B3a) is used three times — injectivity of \(\mathcal H\), positive
curvature along \(\ker\Lambda_\beta\) in the flat limit, and strictness of
\(S^{\mathrm{LB}}\prec P_{11}^{\mathrm{RE}}\); (B3b) enters exactly once, to supply
coercivity along \(\operatorname{range}(\mathcal H)\) when \(\Lambda_\gamma=0\). This is
the same division of labour as v4 Remark 5.2 on the \(\gamma\) side.

**Remark B.2 (the sandwich is oriented opposite to \(\widetilde C_d\)).** On the
\(\gamma\) side, v4 Lemma 14 gives
\(\{\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\le2d\}\subseteq\widetilde C_d
\subseteq\{\|\gamma-\gamma^\star\|^2_{S_\flat}\le2d\}\): the *prior-side* metric bounds
\(\widetilde C_d\) from **inside**. Here the prior-side metric \(\Lambda_\beta\) bounds
\(\widetilde B(r)\) from **outside**. The reason is structural, not a sign slip:
\(\Psi\) is the deficiency gap of the *concave* profile \(\mathcal D\), and satisfies
\(\Psi\le\tfrac12\|\cdot-\gamma^\star\|^2_{P_{11}}\) through \(\bar\Phi\ge0\) (v4 Lemma
13(3)), whereas \(\Xi\) is the negative-log-density gap itself, bounded **below** by its
own strong-convexity modulus. Compactness therefore comes from opposite halves of the
two sandwiches.

**Remark B.3 (design level versus certified level).** Step 5 gives an *exact* level
\(r^\star(\delta_2)\), and Step 6 a *proved finite* upper bound \(r_{\mathrm{P2}}\), but
neither is what the pipeline evaluates: §0.2 uses
\(r_{\mathrm{Gauss}}(n,\delta_2)=\tfrac12\chi^2_{n,1-\delta_2}\), which reproduces
\(\delta_2\) exactly under the Laplace reference and to \(O(N^{-1/2})\) under BvM. Since
\(r_{\mathrm{Gauss}}\ll r_{\mathrm{P2}}\) at fixed \(n\) (roughly
\(r_{\mathrm{P2}}(n,\delta)\approx r_{\mathrm{Gauss}}(2n,\delta)\)), the design route is
the sharper of the two but carries an asymptotic remainder; the certified route is
conservative but unconditional. Only the design route is currently wired into the drift
constants, which is why the \(\pi_\beta(\widetilde B(\delta_2)^{\,c})\le\delta_2\) step in
**Proposition R-Cert** is stated as a design target with an asymptotic justification.

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

That these floors are **attained and strictly positive**, and that the resulting
spectrum satisfies \(\kappa_{\max}^{\mathrm{LB}}(\delta_2)<1\) — so that
\(\lambda(\delta_2)<1\) and \(b(\delta_2)<\infty\) below are non-degenerate — is
**Lemma B-Cert (iv)** (§0.3).

All drift constants below are **functions of \(\delta_2\)** through this spectrum
(unless the prior-only fallback is used). **Monotonicity** (Lemma B-Cert (v), proved
rather than typical): smaller
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

#### Certificate notation (for Proposition R-Cert)

Fix a GLMM with population block \(\gamma\in\mathbb R^q\), group effects
\(\beta=(\beta_1,\ldots,\beta_J)\), proper priors, and data \(y\). Write
\(\pi(\gamma,\beta\mid y)\) for the joint posterior. The symbols in the minorization,
drift, and auxiliary tables above are used; additionally:

| Object | Definition |
|--------|------------|
| \(\beta^\dagger\) | Marginal \(\beta\)-mode after integrating \(\gamma\) (`beta_marginal_mode()`; `BETA_MARGINAL_MODE_LEVELSET.md`) |
| \(\Xi(\beta)\) | \(f(\beta)-f(\beta^\dagger)\) with \(f(\beta)=-\log\widetilde\pi(\beta\mid y)\) |
| \(\widetilde B(\delta_2)\) | \(\{\beta:\Xi(\beta)\le r_{\mathrm{Gauss}}(Jp_{\mathrm{re}},\delta_2)\}\); design target \(\pi_\beta(\widetilde B^{\,c})\approx\delta_2\) |
| \(\pi_\gamma\) | Full \(\gamma\)-marginal posterior \(\pi_\gamma(\cdot\mid y)\) (outer sample space) |
| \(P_\gamma\) | One-sweep \(\gamma\)-marginal kernel of **\(\beta\)-restricted** Gibbs (§0): draw \(\beta_n\in\widetilde B(\delta_2)\), then \(\gamma_n\sim\pi(\gamma\mid\beta_n,y)\) |
| \(S^{\mathrm{LB}}(\delta_2)\) | \(\sum_j H_j^\top P_b\,(P_{22,j}^{\mathrm{LB}}(\delta_2))^{-1}P_b H_j\) |

All other symbols (\(\gamma^\star\), \(\Psi\), \(\widetilde C_d\), \(\varepsilon(d)\),
\(\kappa_i^{\mathrm{LB}}(\delta_2)\), \(V_{\sup}(d)\), etc.) are as in the tables above.

---

> **Proposition R-Cert (Rosenthal \(\gamma\)-marginal TV certificate).**
>
> **Assumptions (model only).** Assume **(H1)–(H4)** for the GLMM of §0 / v4 §3.1:
>
> 1. **(H1 — proper posterior)** \(\pi(\gamma,\beta\mid y)\) is a probability measure
>    (automatic when \(\Lambda_\gamma\succ0\) and the likelihood is bounded; at
>    \(\Lambda_\gamma=0\), derived from (H2)+(H3), v4 §3.2).
> 2. **(H2 — log-concave conditionals)** For each \(j\), \(\pi(\beta_j\mid\gamma,y_j)\) is
>    log-concave in \(\beta_j\) (canonical GLM families supported by the package).
> 3. **(H3 — rank / estimability)** **(H3a)** \(P_{11}^{\mathrm{RE}}=\sum_j H_j^\top P_b H_j\succ0\);
>    **(H3b)** a finite group-wise MLE exists for each \(j\) (no complete separation; needed
>    when \(\Lambda_\gamma=0\)).
> 4. **(Marginal \(\beta\) profile)** After integrating \(\gamma\), the marginal
>    log-density \(f(\beta)=-\log\widetilde\pi(\beta\mid y)\) has a mode \(\beta^\dagger\) and
>    convex profile \(\Xi(\beta)=f(\beta)-f(\beta^\dagger)\) with coercive level sets.
>    This is **not an extra hypothesis**: it is *derived* from (H1)–(H3) in
>    **Lemma B-Cert** (§0.3, parts (i)–(iii)), and is listed here only to name the
>    objects used below (`BETA_MARGINAL_MODE_LEVELSET.md`; Prop 3.1 / Prop 4.1 in
>    `JOINT_GAMMA_BETA_TV_CERTIFICATE.md`).
>
> **Statement (existence of a certified restricted sampler).**
>
> Fix tail budgets \(\delta_2\in(0,1)\) and \(\delta>0\). Under **(H1)–(H4)**, there exist:
>
> - a compact convex **\(\beta\)-restricted set** \(\widetilde B(\delta_2)\) anchored at
>   \(\beta^\dagger\) with design target \(\pi_\beta(\widetilde B(\delta_2)^{\,c})\approx\delta_2\);
> - a deficiency level \(d=d(\delta)\) and certified **small set**
>   \(\widetilde C_d=\{\gamma:\Psi(\gamma)\le d\}\) with \(\gamma^\star\in\widetilde C_d\);
> - a **restricted two-block Gibbs sampler** whose one-sweep \(\gamma\)-marginal kernel is
>   \(P_\gamma=P_{\gamma\mid\widetilde B(\delta_2)}\) (§0: draw \(\beta_n\in\widetilde B(\delta_2)\),
>   then \(\gamma_n\sim\pi(\gamma\mid\beta_n,y)\)),
>
> such that **both** Rosenthal hypotheses hold for \(P_\gamma\) with small set
> \(C=\widetilde C_d\) and Lyapunov function \(V(\gamma)=1+\tfrac12\|\gamma-\gamma^\star\|_{P_{11}}^2\):
>
> - **Minorization (proved):** \(P_\gamma(\gamma,\cdot)\ge\varepsilon(d)\,Q(\cdot)\) for all
>   \(\gamma\in\widetilde C_d\), with \(\varepsilon(d)=e^{-d}\varepsilon(\gamma^\star)\) and
>   \(Q=N(\gamma^\star,P_{11}^{-1})\) (C05 Lemma 17(a); v4 §5).
> - **Foster drift (proved):** \((P_\gamma V)(\gamma)\le\lambda(\delta_2)V(\gamma)+b(\delta_2)\,
>   \mathbb I_{\widetilde C_d^c}(\gamma)\) with \(\lambda(\delta_2)=\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\)
>   and \(b(\delta_2)=1-\lambda(\delta_2)+q/2+C_\beta^{+}(\delta_2)\) from the floor coupling
>   spectrum on \(\widetilde B(\delta_2)\) (§2).
>
> Consequently, for any start \(\gamma_0\in\mathbb R^q\), any integer \(k\ge1\), and any
> \(\alpha\in\bigl(\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2,\,1\bigr)\), the
> **one-line substitution** bound for the **full** marginal total variation
> \(\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_\gamma\|_{TV}\) holds:
>
> \[
> \boxed{\small
> \begin{aligned}
> \bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_\gamma\bigr\|_{TV}
> \ \le\
> &\bigl(1-\varepsilon(d)\bigr)^{\lfloor\alpha k\rfloor} \\[4pt]
> &\quad+
> \frac{
> 1+2\Bigl[
> 1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
> \Bigr]
> +\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2 V_{\sup}(d)
> }{\alpha} \\[4pt]
> &\quad\quad\times
> \frac{
> 1+2\Bigl[
> 1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
> \Bigr]
> +\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\Bigl(1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\Bigr)
> }{
> 1+2\Bigl[
> 1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
> \Bigr]\Big/\Bigl(1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\Bigr)
> }\,
> \alpha^{k} \\[4pt]
> &\quad\quad\quad+\;\pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr),
> \end{aligned}
> }
> \]
>
> with \(\varepsilon(d)=e^{-d}\varepsilon(\gamma^\star)\),
> \(V_{\sup}(d)=\sup_{\gamma\in\widetilde C_d}\bigl(1+\tfrac12\|\gamma-\gamma^\star\|_{P_{11}}^2\bigr)\), and
> \(\pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr)\le\delta_2\) by construction of
> \(\widetilde B(\delta_2)\) (design target; asymptotically under Laplace / BvM at
> \(\beta^\dagger\), `BETA_MARGINAL_MODE_LEVELSET.md` §4).
>
> *Proof sketch.* (H1)–(H3) supply C05 mode \(\gamma^\star\), closure \(\varepsilon(\gamma^\star)\),
> and Lemma 17(a). (H4) and `deficiency_calibrate(\delta,\ldots)` give \(\widetilde C_d\) with
> \(\gamma^\star\in\widetilde C_d\) (v4 Lemma 23). The \(\beta\)-safe set \(\widetilde B(\delta_2)\)
> and restricted sampler are constructed in `BETA_MARGINAL_MODE_LEVELSET.md` /
> `beta_marginal_safe_set()`. §2 verifies Foster drift from weight floors on
> \(\widetilde B(\delta_2)\). Apply Rosenthal (1995), Theorem 12, to bound
> \(\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\|_{TV}\), substitute the floor
> spectrum for \((\lambda,b)\), then apply the triangle inequality
> \(\|P_\gamma^{k}-\pi_\gamma\|_{TV}\le
> \|P_\gamma^{k}-\pi_{\gamma\mid\widetilde B}\|_{TV}+\|\pi_{\gamma\mid\widetilde B}-\pi_\gamma\|_{TV}\)
> with \(\|\pi_{\gamma\mid\widetilde B}-\pi_\gamma\|_{TV}\le\pi_\beta(\widetilde B(\delta_2)^{\,c})\le\delta_2\)
> (§0). \(\square\)
>
> **Optimality (within the certified family).**
>
> For fixed \((\delta_2,\delta,\alpha,k)\), the **dynamic** part of the bound (all terms
> before \(+\,\pi_\beta(\widetilde B(\delta_2)^{\,c})\)) is **minimized** (among feasible
> calibrations meeting the \(\gamma\)-tail budget \(\delta\)) by taking \(\widetilde C_d\) as
> **small as possible** — equivalently, the smallest \(d=d(\delta)\) with
> \(\pi_\gamma(\widetilde C_d^{\,c})\le\delta\) — because \(\varepsilon(d)=e^{-d}\varepsilon(\gamma^\star)\)
> is increasing and \(V_{\sup}(d)\) is decreasing in \(d\) on the certified range (§3.1 below).
> For fixed \((d,\delta_2,\alpha,k)\), the drift factor is minimized by initializing at
> \(\gamma_0=\gamma^\star\), since \(V(\gamma_0)=1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\ge1\)
> with equality only at the mode. The \(\beta\)-truncation term
> \(\pi_\beta(\widetilde B(\delta_2)^{\,c})\le\delta_2\) is fixed by the choice of
> \(\widetilde B(\delta_2)\).
>
> **Sharpest certified certificate (one-line substitution).**
>
> In the simultaneous limit **smallest feasible** \(\widetilde C_d\) (profile limit \(d\downarrow0\),
> so \(\varepsilon(d)\to\varepsilon(\gamma^\star)\) and \(V_{\sup}(d)\to1\)) and start
> \(\gamma_0=\gamma^\star\), the full-marginal bound reduces to
>
> \[
> \boxed{\small
> \begin{aligned}
> \bigl\|P_\gamma^{k}(\gamma^\star,\cdot)-\pi_\gamma\bigr\|_{TV}
> \ \le\
> &\bigl(1-\varepsilon(\gamma^\star)\bigr)^{\lfloor\alpha k\rfloor} \\[4pt]
> &\quad+
> \frac{
> 1+2\Bigl[
> 1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
> \Bigr]
> +\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2
> }{\alpha} \\[4pt]
> &\quad\quad\times
> \frac{
> 1+2\Bigl[
> 1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
> \Bigr]
> +\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2
> }{
> 1+2\Bigl[
> 1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2+\tfrac q2+\tfrac12\sum_{i=1}^q\kappa_i^{\mathrm{LB}}(\delta_2)
> \Bigr]\Big/\Bigl(1-\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\Bigr)
> }\,
> \alpha^{k} \\[4pt]
> &\quad\quad\quad+\;\pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr),
> \end{aligned}
> }
> \]
>
> with \(\pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr)\le\delta_2\). This is the
> **display form** used when \(\gamma_0=\gamma^\star\) and the certified set is taken as
> tight as the \(\delta\)-budget allows. For finite \(d>0\) or \(\gamma_0\neq\gamma^\star\),
> use the general bound in the statement above (with \(\varepsilon(d)\), \(V_{\sup}(d)\), and
> \(V(\gamma_0)\)); when \(\gamma_0\) is only approximate, bound \(V(\gamma_0)\) by the
> matrix mode-distance bound on \(\|\gamma_0-\gamma^\star\|_{\Pi^{\mathrm{LB}}}\)
> (**Remarks 2** and **4** below; derivation in §3.1).

**Remark 1 (the marginal mode and \(\varepsilon(\gamma^\star)\) are not directly
computable).** Outside the Gaussian **closure** regime, the marginal posterior mode
\(\gamma^\star\) is the unique fixed point of the mean map (v4 §3.3, Lemma 3),

\[
\gamma^\star=M(\gamma^\star),
\qquad
M(\gamma)=P_{11}^{-1}\Bigl(\Lambda_\gamma\mu_0
+\sum_j H_j^\top P_b\,\mathbb E[\beta_j\mid\gamma,y]\Bigr),
\]

so locating \(\gamma^\star\) is equivalent to solving \(\nabla\Phi(\gamma)=0\) for the
population score \(\nabla\Phi(\gamma)=P_{11}\bigl(\gamma-M(\gamma)\bigr)\). Every
evaluation of either object requires the **conditional expectations**
\(\mathbb E[\beta_j\mid\gamma,y]\), which for canonical GLM group blocks have no closed
form and are approximated by **Monte Carlo** from \(\pi(\beta_j\mid\gamma,y_j)\). The
iterate \(\tilde\gamma\) returned by a simulation-based search is therefore an
**approximate mode**, not \(\gamma^\star\).

The same obstruction reaches the minorization constant. For non-Gaussian models the
one-step kernel \(q(\gamma'\mid\gamma)=\int\phi_q(\gamma';m(\beta),P_{11}^{-1})\,
\pi(\beta\mid\gamma,y)\,d\beta\) is a **Gaussian mixture**, so the profile
\[
\varepsilon(\gamma)
=\exp\Bigl\{\inf_{\gamma'\in\mathbb R^q} g(\gamma'\mid\gamma)\Bigr\},
\qquad
g(\gamma'\mid\gamma):=\log q(\gamma'\mid\gamma)-\log q_Q(\gamma'),
\]
has **no** closed form at the mode (Definition 5, v4 §5). The scalar
\(\varepsilon^\star=\varepsilon(\gamma^\star)\) is obtained only by **numerical convex
optimization** (Lemma 10: strict convexity of \(\gamma'\mapsto g(\gamma'\mid\gamma)\);
Lemma 11: attainment), and it must be evaluated at a centre that is itself known only
approximately. In practice one **substitutes** \(\tilde\gamma\) and MC-estimated mixing
laws:
\[
\widehat\varepsilon^\star
:= \varepsilon(\tilde\gamma)
\approx
\exp\Bigl\{\min_{\gamma'} g(\gamma'\mid\tilde\gamma)\Bigr\},
\]
with \(g(\cdot\mid\tilde\gamma)\) built from sample means \(\hat b_j\) in place of
\(\mathbb E[\beta_j\mid\tilde\gamma,y]\). The deficiency level \(d=d(\delta)\) and
certified set \(\widetilde C_d=\{\gamma:\Psi(\gamma)\le d\}\),
\(\Psi(\gamma)=\mathcal D(\gamma^\star)-\mathcal D(\gamma)\),
\(\mathcal D(\gamma)=\log\varepsilon(\gamma)\), are calibrated from the **same**
plug-in centre (in practice \(\gamma^\star\) is replaced by \(\tilde\gamma\) in
\(\mathcal D(\gamma^\star)\), in \(J(\gamma)\), and in \(Q\); coupling spectrum of
\(J(\tilde\gamma)=P_{11}^{-1}\sum_j H_j^\top P_b V_j(\tilde\gamma)P_b H_j\)), so the
bound is evaluated with
\[
\varepsilon(d)=e^{-d}\,\widehat\varepsilon^\star,
\qquad
Q=N(\tilde\gamma,P_{11}^{-1})
\]
instead of the true refresh \(Q=N(\gamma^\star,P_{11}^{-1})\). **Remarks 2–5** quantify
what this substitution costs on each side of the bound: Remarks 2 and 4 treat the drift
side, Remarks 3 and 5 the minorization side.

**Remark 2 (certified \(V(\gamma_0)\) when the gradient is known).** The Rosenthal drift
factor depends on \(\gamma^\star\) only through
\(V(\gamma_0)=1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\), and this can be bounded
from the gradient **at \(\gamma_0\) alone**, without locating \(\gamma^\star\). The
marginal Hessian obeys the data-precision floor

\[
\nabla^2\Phi\ \succeq\ \Pi^{\mathrm{LB}}:=P_{11}-S^{\mathrm{LB}}(\delta_2)\ \succ\ 0,
\]

positive definite because \(\kappa_{\max}^{\mathrm{LB}}(\delta_2)<1\)
(**Lemma B-Cert** (iv), §0.3). Since \(\Phi\) is \(\Pi^{\mathrm{LB}}\)-strongly convex
and \(\nabla\Phi(\gamma^\star)=0\), strong convexity gives, writing
\(\Delta:=\gamma_0-\gamma^\star\) (written \(d\) in §3.1; renamed here to avoid collision
with the deficiency level \(d=d(\delta)\)) and \(g_0:=\nabla\Phi(\gamma_0)\),

\[
\Delta^\top\Pi^{\mathrm{LB}}\Delta
\ \le\
g_0^\top(\Pi^{\mathrm{LB}})^{-1}g_0
=:B_{\mathrm{obs}}
\]

(§3.1, “Correct matrix bound (no eigenvalue collapse)”). Converting to the \(P_{11}\)
metric used by \(V\) through the exact split
\(P_{11}=\Pi^{\mathrm{LB}}+S^{\mathrm{LB}}\), with worst-case alignment
\(\lambda_{\max}\bigl((\Pi^{\mathrm{LB}})^{-1}S^{\mathrm{LB}}\bigr)=w_{\max}^{\mathrm{LB}}\),

\[
\boxed{\small
V(\gamma_0)
= 1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2
\ \le\
1+\tfrac12\,B_{\mathrm{obs}}\,\bigl(1+w_{\max}^{\mathrm{LB}}\bigr),
\qquad
B_{\mathrm{obs}}=g_0^\top(\Pi^{\mathrm{LB}})^{-1}g_0 .
}
\]

This bound is **deterministic**: it holds exactly whenever \(g_0\) is known exactly, as in
the Gaussian closure regime, and it requires no knowledge of \(\gamma^\star\) beyond the
fact that it is the stationary point. **Computation:** one Cholesky of
\(\Pi^{\mathrm{LB}}\) and one solve \(v=(\Pi^{\mathrm{LB}})^{-1}g_0\), then \(g_0^\top v\);
no eigendecomposition of \(\Phi\) is needed, and \(w_{\max}^{\mathrm{LB}}\) comes from
`floor_coupling_eigenvalues()$kappa_lb`. Pair the gradient with
\((\Pi^{\mathrm{LB}})^{-1}\), **not** \(P_{11}^{-1}\): these are different norms.

**Remark 3 (any anchor yields a valid \(\varepsilon\)).** Unlike \(V(\gamma_0)\), the
minorization constant needs no correction for anchor error. The profile
\(\varepsilon(\cdot)\) is **maximal** at the marginal mode (Lemma 13(2)), so for every
\(\gamma\in\mathbb R^q\)

\[
\varepsilon(\gamma)\ \le\ \varepsilon(\gamma^\star).
\]

Hence \(\widehat\varepsilon^\star=\varepsilon(\tilde\gamma)\) evaluated at any convenient
anchor is a **valid** minorization constant, not merely an approximation to one: the
Rosenthal hypothesis \(P_\gamma(\gamma,\cdot)\ge\varepsilon\,Q(\cdot)\) still holds on the
small set, and since the bound is monotone in \(\varepsilon\) through
\((1-\varepsilon(d))^{\lfloor\alpha k\rfloor}\), the certificate remains **true** with only
the **rate** degraded. The anchor error is therefore in the safe direction, and its size is
observable rather than something that must itself be certified — which is why no analogue
of \(B_{\mathrm{cert}}\) is needed here.

Two conditions attach. First, the reference measure \(Q=N(\tilde\gamma,P_{11}^{-1})\) and
the small set \(\widetilde C_d\) must be anchored at the **same** \(\tilde\gamma\), so that
the minorization and drift hypotheses refer to one set; this is exactly what the common
plug-in centre of **Remark 1** secures. Second, the statement concerns the *exact* profile
\(\varepsilon(\tilde\gamma)\). When that profile is itself estimated by simulation, the
estimate can exceed \(\varepsilon(\tilde\gamma)\) and validity is no longer automatic; see
**Remark 5**.

**Remark 4 (certified \(V(\gamma_0)\) when the gradient is available only by simulation).**
Initializing at \(\gamma_0=\tilde\gamma\) inflates the Lyapunov start value
\(V(\gamma_0)=1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\), so the Rosenthal drift
factor is understated if \(V(\gamma_0)\) is set to the sharp value \(1\) that is valid only
at \(\gamma_0=\gamma^\star\). When the score \(g=\nabla\Phi(\tilde\gamma)\) is itself a
Monte Carlo estimate, \(B_{\mathrm{obs}}\) of **Remark 2** inherits that noise and must be
replaced by an upper confidence limit. Propagating the MC uncertainty of the mean map
through the quadratic form (§3.1, “Monte Carlo: plug-in bound and high-probability
version”) gives

\[
B_{\mathrm{cert}}
:= B_{\mathrm{obs}}
+z_{1-\alpha}\,\widehat{\mathrm{se}}\bigl(g^\top(\Pi^{\mathrm{LB}})^{-1}g\bigr),
\]

\[
\boxed{\small
V(\gamma_0)
= 1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2
\ \le\
1+\tfrac12\,B_{\mathrm{cert}}\,\bigl(1+w_{\max}^{\mathrm{LB}}\bigr),
}
\]

with \(\Pi^{\mathrm{LB}}=P_{11}-S^{\mathrm{LB}}(\delta_2)\) and
\(w_{\max}^{\mathrm{LB}}=\max_i\kappa_i^{\mathrm{LB}}/(1-\kappa_i^{\mathrm{LB}})\)
for \(\kappa_i^{\mathrm{LB}}=\mathrm{eig}_i\bigl(P_{11}^{-1/2}S^{\mathrm{LB}}P_{11}^{-1/2}\bigr)\).
This is the statement of **Remark 2** with \(B_{\mathrm{obs}}\) replaced by its upper
confidence limit, and it holds at level \(1-\alpha\) **over the simulation**, not
deterministically. (Here and in **Remark 5**, \(\alpha\) in \(z_{1-\alpha}\) is the
**confidence level**, unrelated to the Rosenthal tuning parameter
\(\alpha\in(\lambda(\delta_2),1)\) of the bound itself; the document reuses the symbol in
both roles, as in §3.1.) Use this certified \(V(\gamma_0)\) — equivalently the acceptance test
\(B_{\mathrm{cert}}\le\tau_B^2\) — in the general bound above whenever \(\gamma_0\) is
simulation-based; see §3.1 (“Bound on \(V(\gamma_0)\) for Rosenthal drift”, recommended
stopping rule). Because this cost enters the drift term **additively** through
\(V(\gamma_0)\), a loose \(B_{\mathrm{cert}}\) is comparatively cheap.

**Remark 5 (confidence bound for \(\varepsilon(\gamma)\) under a Monte Carlo E-step).**
**Remark 3** makes the exact profile \(\varepsilon(\tilde\gamma)\) valid at any anchor, but
with a Monte Carlo E-step one observes only an estimate
\(\widehat\varepsilon(\tilde\gamma)\), and simulation noise can push that estimate
**above** \(\varepsilon(\tilde\gamma)\) — the unsafe direction. A lower confidence bound is
therefore required. With \(\beta_1,\ldots,\beta_n\) drawn i.i.d. from
\(\pi(\beta\mid\tilde\gamma,y)\), write

\[
\varphi_m(\gamma')=\phi_q\bigl(\gamma';m(\beta_m),P_{11}^{-1}\bigr),
\qquad
\hat q(\gamma')=\frac1n\sum_{m=1}^n\varphi_m(\gamma'),
\qquad
\hat\gamma'=\arg\min_{\gamma'}\bigl[\log\hat q(\gamma')-\log q_Q(\gamma')\bigr],
\]

and let \(\hat V=\frac1{n-1}\sum_{m}\bigl(\varphi_m(\hat\gamma')-\hat q(\hat\gamma')\bigr)^2\)
be the sample variance of the mixture summands at \(\hat\gamma'\). Then, at level
\(1-\alpha\),

\[
\boxed{\small
\begin{aligned}
\log\varepsilon^{\mathrm{LCB}}(\tilde\gamma)
&=\Bigl[\log\hat q(\hat\gamma')-\log q_Q(\hat\gamma')\Bigr]
-z_{1-\alpha}\,\frac{\hat V^{1/2}}{\hat q(\hat\gamma')\sqrt n},
\\[4pt]
\varepsilon(d)^{\mathrm{LCB}}
&=e^{-d}\,\varepsilon^{\mathrm{LCB}}(\tilde\gamma).
\end{aligned}
}
\]

Both **systematic** errors already run in the safe direction: \(\log\) is concave, so
\(\mathbb E[\log\hat q]\le\log q\), and the infimum of a noisy function lies below the
infimum of its mean. The interval is therefore correcting **realized** noise rather than
bias. Two limitations. The interval is **pointwise** at \(\hat\gamma'\), so on its own it
does not certify the infimum over all \(\gamma'\) (§8, Open items). And because
\(\varepsilon\) enters **multiplicatively** through
\((1-\varepsilon(d))^{\lfloor\alpha k\rfloor}\) while \(V(\gamma_0)\) enters additively, a
loose lower bound here costs far more sweeps than a loose \(B_{\mathrm{cert}}\) in
**Remark 4**; simulation effort should be allocated accordingly.

**Error budget across Remarks 4 and 5.** The two are separate probabilistic statements. To
report a certificate valid at overall level \(1-\alpha\), split the budget between them
(for instance \(\alpha/2\) each, Bonferroni) rather than using \(\alpha\) in both.

The subsections below restate the same one-line substitutions for implementation
(`gamma_beta_tv_certificate()`, `optimal_rosenthal_tv_bound()`). Rosenthal compact
forms \((U/\alpha)(\cdots)\alpha^k\) appear there for cross-reference only;
**Proposition R-Cert** is the certificate statement.

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

*Restates the **Proposition R-Cert** conclusion (full \(\pi_\gamma\) target with
\(+\,\pi_\beta(\widetilde B(\delta_2)^{\,c})\)); notation in §3.1 above.*

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
\boxed{\small
\begin{aligned}
\bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_\gamma\bigr\|_{TV}
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
\alpha^{k} \\[4pt]
&\quad\quad\quad+\;\pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr),
\end{aligned}
}
\]

with \(\alpha\in\bigl(\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2,\,1\bigr)\) and
\(\pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr)\le\delta_2\).

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

**Mode distance when \(\gamma_0\neq\gamma^\star\).** Sharp display sets \(V(\gamma_0)=1\);
general display needs \(\|\gamma_0-\gamma^\star\|_{P_{11}}\). The gradient and the
data-precision floor supply a **matrix** bound without collapsing to a single
eigenvalue.

From §2.3, \(S^{\mathrm{LB}}(\delta_2)\preceq S(\gamma)\) on \(\widetilde B(\delta_2)\), so
marginal observed information is bounded below by

\[
\Pi^{\mathrm{LB}}(\delta_2)
:= P_{11}-S^{\mathrm{LB}}(\delta_2)
\;\preceq\;
\nabla^2\Phi(\gamma)
\]

(Lemma 5 / C05: \(\nabla^2\Phi=P_{11}-S\).) The EM/Fisher identity is
\(\nabla\Phi(\gamma)=P_{11}(\gamma-M(\gamma))\).

**Correct matrix bound (no eigenvalue collapse).** Let
\(\tilde\gamma\) be an approximate mode (e.g.\ `population_mode()$fixef`) and

\[
g:=\nabla\Phi(\tilde\gamma)=P_{11}\bigl(\tilde\gamma-M(\tilde\gamma)\bigr).
\]

If \(\nabla^2\Phi\succeq\Pi^{\mathrm{LB}}\) (from the data-precision floor), then for
\(d:=\tilde\gamma-\gamma^\star\),

\[
\boxed{
\|d\|_{\Pi^{\mathrm{LB}}}^2
=
d^\top \Pi^{\mathrm{LB}} d
\;\le\;
g^\top (\Pi^{\mathrm{LB}})^{-1} g
=
\|\nabla\Phi(\tilde\gamma)\|_{(\Pi^{\mathrm{LB}})^{-1}}^2.
}
\]

This is the standard strong-convexity inequality: at the minimizer
\(\nabla\Phi(\gamma^\star)=0\), with \(\nabla^2\Phi\succeq\Pi^{\mathrm{LB}}\),

\[
\|\tilde\gamma-\gamma^\star\|_{\Pi^{\mathrm{LB}}}
\;\le\;
\|\nabla\Phi(\tilde\gamma)\|_{(\Pi^{\mathrm{LB}})^{-1}}.
\]

No eigenvalues appear. **Computation:** Cholesky \(\Pi^{\mathrm{LB}}\), one solve
\(v=(\Pi^{\mathrm{LB}})^{-1}g\), then \(g^\top v\). In the package,
`floor_coupling_eigenvalues()` forms \(S^{\mathrm{LB}}\); \(\Pi^{\mathrm{LB}}=P_{11}-S^{\mathrm{LB}}\)
uses the same `population_mode()$p11$P11`. The EM screen residual
\(\|\tilde\gamma-M(\tilde\gamma)\|_{P_{11}}=\|\nabla\Phi(\tilde\gamma)\|_{P_{11}^{-1}}\) is a
**different** norm — pair the gradient with \((\Pi^{\mathrm{LB}})^{-1}\), not \(P_{11}^{-1}\).
Conversion to Rosenthal \(V(\gamma_0)\) is below (“Bound on \(V(\gamma_0)\) for Rosenthal drift”).

**Monte Carlo: plug-in bound and high-probability version.** With
`estep = "mc"`, \(M(\tilde\gamma)\) uses sample means \(\hat b_j\) from \(n\) draws
per group, so the bound above is **random**: the deterministic inequality applies to
the **true** gradient, while the computed object uses the plug-in
\(g_{\mathrm{obs}}:=P_{11}(\tilde\gamma-M(\hat b))\).

The mean map is affine in \(b\):
\(M(b)=P_{11}^{-1}(\Lambda_\gamma\mu_0+\sum_j H_j^\top P_b b_j)\), hence
\(g(b)=g_0-Tb\) for fixed \(\tilde\gamma\). With \(\hat b\) the sample mean of \(n\)
iid draws,

\[
\mathrm{Cov}(\hat b)\ \approx\ \frac{1}{n}\,\Sigma_\beta(\tilde\gamma),
\qquad
\mathrm{Cov}(g)\ =\ T\,\mathrm{Cov}(\hat b)\,T^\top
\ \propto\ \frac{1}{n}.
\]

The package records elementwise MC standard errors
\(\widehat{\mathrm{se}}(\hat b_{jk})=\mathrm{sd}(\text{draws})/\sqrt{n}\) as
`b_mc_se`, and propagates them through the mean map (`.c05_mc_delta_floor()` in
`R/c05_em.R`):

\[
\mathrm{Cov}_\gamma
:= \sum_j A_j\,\mathrm{diag}(\widehat{\mathrm{se}}^2)\,A_j^\top,
\qquad
A_j := P_{11}^{-1}H_j^\top P_b,
\]

\[
\boxed{
\sigma_\delta
:= \texttt{mc\_delta\_floor}
= \sqrt{\mathrm{tr}(P_{11}\,\mathrm{Cov}_\gamma)}
\ \propto\ n^{-1/2}.
}
\]

This is the MC standard deviation of the EM residual
\(\delta=\|\tilde\gamma-M(\hat b)\|_{P_{11}}=\|g\|_{P_{11}^{-1}}\). Current stopping
uses \(\texttt{tol\_eff}=\max(\texttt{tol},\, z_{1-\alpha/2}\,\sigma_\delta)\) with
`mc_alpha` (see `population_mode()`).

Write \(\Pi:=\Pi^{\mathrm{LB}}\) and \(\Sigma_g(n):=\mathrm{Cov}(g)\approx P_{11}\,\mathrm{Cov}_\gamma\,P_{11}\).

**Plug-in (random) bound.**

\[
B_{\mathrm{obs}}
:= g_{\mathrm{obs}}^\top \Pi^{-1} g_{\mathrm{obs}},
\qquad
\|d\|_{\Pi}^2 \ \text{(estimated)}\ \le\ B_{\mathrm{obs}}.
\]

**Conservative high-probability version.** For a **proved** statement at level
\(1-\alpha\), inflate the quadratic form by its MC sampling uncertainty. With
\(\mu_g:=\mathbb E[g]\) and \(\Sigma_g(n)\) from the propagation above,

\[
\mathrm{Var}(g^\top\Pi^{-1}g)
\ \approx\
4\,\mu_g^\top\Pi^{-1}\Sigma_g(n)\Pi^{-1}\mu_g
+ 2\,\mathrm{tr}\bigl((\Pi^{-1}\Sigma_g(n))^2\bigr),
\]

so the standard error of the plug-in bound is \(O(n^{-1/2})\) when \(\mu_g\neq0\)
and \(O(n^{-1})\) near a fixed point (\(\mu_g\approx0\)). A one-sided conservative
certificate is

\[
\boxed{
\|d\|_{\Pi}^2
\ \le\
B_{\mathrm{obs}}
+ z_{1-\alpha}\,
\widehat{\mathrm{se}}(g^\top\Pi^{-1}g),
}
\]

with \(\widehat{\mathrm{se}}\) from the delta method and `b_mc_se`, or from paired
MC replicates (`mc_stable` \(>1\)). Since \(\Pi\preceq P_{11}\), hence
\((\Pi^{\mathrm{LB}})^{-1}\succeq P_{11}^{-1}\),

\[
\|g\|_{(\Pi^{\mathrm{LB}})^{-1}}^2
\ \ge\
\|g\|_{P_{11}^{-1}}^2
=\delta^2,
\]

so MC noise in \(B_{\mathrm{obs}}\) is **not smaller** than that in \(\delta\).

**Bound on \(V(\gamma_0)\) for Rosenthal drift.** Rosenthal enters through
\(V(\gamma_0)=1+V_0(\gamma_0)\) with
\(V_0(\gamma_0)=\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\) (§2.1, MT shift).
The matrix bound controls \(\|d\|_{\Pi^{\mathrm{LB}}}\) for
\(d:=\gamma_0-\gamma^\star\); convert to the \(P_{11}\) metric used by \(V\).

Write \(\Pi:=\Pi^{\mathrm{LB}}\) and the certified quadratic bound

\[
B_{\mathrm{cert}}
:= B_{\mathrm{obs}}+z_{1-\alpha}\,\widehat{\mathrm{se}}(g^\top\Pi^{-1}g),
\qquad
d^\top\Pi d\le B_{\mathrm{cert}}.
\]

Use the exact split \(P_{11}=\Pi+S^{\mathrm{LB}}\):

\[
\|d\|_{P_{11}}^2
= d^\top\Pi d + d^\top S^{\mathrm{LB}} d.
\]

Given only \(d^\top\Pi d\le B_{\mathrm{cert}}\), the worst-case alignment for the
\(S^{\mathrm{LB}}\) term is

\[
\max_{d^\top\Pi d\le B_{\mathrm{cert}}} d^\top S^{\mathrm{LB}} d
= B_{\mathrm{cert}}\,
\lambda_{\max}\!\bigl(\Pi^{-1}S^{\mathrm{LB}}\bigr).
\]

With \(A^{\mathrm{LB}}=P_{11}^{-1/2}S^{\mathrm{LB}}P_{11}^{-1/2}\) and
\(\kappa_i^{\mathrm{LB}}=\mathrm{eig}_i(A^{\mathrm{LB}})\) from
`floor_coupling_eigenvalues()`,

\[
\lambda_{\max}(\Pi^{-1}S^{\mathrm{LB}})
= \max_i \frac{\kappa_i^{\mathrm{LB}}}{1-\kappa_i^{\mathrm{LB}}}
=: w_{\max}^{\mathrm{LB}}
\]

(the C05 deficiency weights \(w_i=\kappa_i/(1-\kappa_i)\), not \(1/(1-\kappa_{\max})\)).
Hence

\[
\boxed{
\|\gamma_0-\gamma^\star\|_{P_{11}}^2
\ \le\
B_{\mathrm{cert}}\,\bigl(1+w_{\max}^{\mathrm{LB}}\bigr),
}
\]

and

\[
\boxed{
V(\gamma_0)
= 1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2
\ \le\
1+\frac{B_{\mathrm{cert}}}{2}\,\bigl(1+w_{\max}^{\mathrm{LB}}\bigr).
}
\]

Plug-in only: replace \(B_{\mathrm{cert}}\) by \(B_{\mathrm{obs}}\). **Computation:**
Cholesky \(\Pi^{\mathrm{LB}}\), solve for \(B_{\mathrm{obs}}\); form
\(w_{\max}^{\mathrm{LB}}\) from `floor_coupling_eigenvalues()$kappa_lb`.

**Start \(\gamma_0\neq\tilde\gamma\).** Split
\(d=(\gamma_0-\tilde\gamma)+(\tilde\gamma-\gamma^\star)=d_0+d_m\). Then

\[
\|d\|_{P_{11}}\le\|d_0\|_{P_{11}}+\|d_m\|_{P_{11}},
\qquad
V_0=\tfrac12\|d\|_{P_{11}}^2
\le \tfrac12\bigl(\|d_0\|_{P_{11}}+\|d_m\|_{P_{11}}\bigr)^2,
\]

with \(\|d_0\|_{P_{11}}\) from the declared start and
\(\|d_m\|_{P_{11}}^2\le B_{\mathrm{cert}}(1+w_{\max}^{\mathrm{LB}})\) at
\(\tilde\gamma=\) `population_mode()$fixef`. Conservative:
\((a+b)^2\le 2(a^2+b^2)\).

**Versus `display_mode = "general"`.** `rosenthal_drift_constants()` in general mode
centres \(V_0\) on `mode$fixef` (\(\hat\gamma\)), so \(\gamma_0=\hat\gamma\) gives
\(V_0=1\) even when \(\hat\gamma\neq\gamma^\star\). The bound above is w.r.t.\
the **true** mode and is the conservative input when certifying drift; pass the
result as `V_sup` / inflated `V_gamma_0` once wired (not implemented in v1).

**Recommended stopping rule for `population_mode()` (MC E-step).** The present
implementation stops when the iterate change
\(\delta_k=\|\gamma^{(k)}-\gamma^{(k+1)}\|_{P_{11}}\) (or the ICM-screen fixed-point
residual \(\|\tilde\gamma-M(\hat b)\|_{P_{11}}\)) satisfies
\(\delta\le\texttt{tol\_eff}=\max(\texttt{tol},\, z_{1-\alpha/2}\,\sigma_\delta)\)
with \(\sigma_\delta=\texttt{mc\_delta\_floor}\propto n^{-1/2}\) and optional
`mc_stable` consecutive passes. That resolves **fixed-point noise in the
\(P_{11}\) metric** only. The bounds above motivate a **two-layer** rule that also
certifies **mode distance** and **\(V(\gamma_0)\)** for Rosenthal.

**Quantities at candidate \(\tilde\gamma\)** (current `fixef`, after an E-step).
Compute \(g_{\mathrm{obs}}\), \(\Pi=\Pi^{\mathrm{LB}}\), \(B_{\mathrm{obs}}\),
\(B_{\mathrm{cert}}\), \(w_{\max}^{\mathrm{LB}}\), and \(V_{\mathrm{cert}}\) as in
the preceding subsections; \(\delta=\|g_{\mathrm{obs}}\|_{P_{11}^{-1}}\).

**Acceptance (certificate mode).** Require **all** of the following for
`mc_stable` consecutive MC replicates (independent seeds at the same
\(\tilde\gamma\)):

\[
\boxed{\small
\begin{aligned}
\text{(A) Fixed point:}\quad &
\delta \le \texttt{tol\_eff}
= \max(\texttt{tol},\, z_\delta\,\sigma_\delta), \\[4pt]
\text{(B) Mode distance:}\quad &
B_{\mathrm{cert}} \le \tau_B^2, \\[4pt]
\text{(C) Rosenthal drift:}\quad &
V_{\mathrm{cert}}
= 1+\tfrac12 B_{\mathrm{cert}}(1+w_{\max}^{\mathrm{LB}})
\le V_{\max}=1+V_{0,\max}.
\end{aligned}
}
\]

(C) is equivalent to (B) when \(\tau_B\) is derived from \(V_{0,\max}\) below.
Use the same `mc_alpha` for \(z_\delta\) and \(z_B\) (e.g.\
\(z_\delta=z_B=z_{1-\alpha/2}\)).

**Linking thresholds to MC (`n`, `mc_delta_floor`).**

1. **MC-native deterministic tolerance.** Under MC, a tiny fixed `tol` (e.g.\
`1e-10`) is not meaningful. Prefer
\(\texttt{tol}=\rho_\delta\,\sigma_\delta\) with \(\rho_\delta\in[0,1)\) (e.g.\
`0.1`), so \(\texttt{tol\_eff}=z_\delta\,\sigma_\delta\) whenever
\(z_\delta\ge\rho_\delta\). On the Gaussian exact path, \(\sigma_\delta=0\) and
`tol` remains purely deterministic.

2. **Mode bound linked to fixed-point MC floor.** Since
\(\|g\|_{(\Pi^{\mathrm{LB}})^{-1}}^2\ge\delta^2\),
\[
\boxed{
\tau_B = \eta_B\,\texttt{tol\_eff},
\qquad \eta_B\ge 1,
}
\]
(e.g.\ \(\eta_B=1\): mode certificate at the same scale as the MC fixed-point
resolution; \(\eta_B<1\) is stricter in the \(\Pi^{\mathrm{LB}}\) metric).

3. **Certificate-linked thresholds from \(V_{0,\max}\).** If the Rosenthal budget
is a target Lyapunov inflation \(V_{0,\max}\) above the sharp baseline \(V=1\),
invert the bound on \(V(\gamma_0)\):
\[
\boxed{
\tau_B^2 = \frac{2\,V_{0,\max}}{1+w_{\max}^{\mathrm{LB}}}.
}
\]
Optionally cap the fixed-point tolerance by the same budget:
\[
\texttt{tol\_eff}\le
\frac{\sqrt{2\,V_{0,\max}}}{1+w_{\max}^{\mathrm{LB}}}.
\]

**Use cases.**

| Mode | Conditions | Typical settings |
|------|------------|------------------|
| Exploratory MC EM | (A) only | Current default: `tol_eff`, `mc_stable` |
| Certificate pipeline | (A) + (B) [+(C)] | Set `V0_max` / `v0_cert_tol`; derive \(\tau_B^2\) |
| Strict certificate | (A)+(B), `mc_stable`\(\ge2\), adaptive `n` | Increase `n` if (A) passes but (B) fails |

**EM loop vs ICM screen.**

- **ICM screen:** already evaluates \(\delta\) at fixed `fixef`; add
\(B_{\mathrm{obs}}\), \(B_{\mathrm{cert}}\) at the same point.
- **EM loop:** \(\delta_k\) is the iterate step; at convergence,
\(\delta_k\approx\|g\|_{P_{11}^{-1}}\). On acceptance, run one final E-step at
\(\tilde\gamma\) and evaluate (B)/(C).

**Adaptive `n`.** If \(\delta\le\texttt{tol\_eff}\) but
\(B_{\mathrm{cert}}>\tau_B^2\), the fixed point is MC-resolved but mode distance
is not certificate-small: increase `n`, re-evaluate (B) at the same
\(\tilde\gamma\) without further EM unless the mean map moves.

**Planned API (not wired in v1).**

```r
population_mode(
  ...,
  tol = NULL,              # NULL => MC-only: tol_eff = z_delta * mc_delta_floor
  mc_alpha = 0.05,
  mc_stable = 2L,
  v0_cert_tol = NULL,      # if set: tau_B^2 = 2*v0/(1 + w_max^LB)
  mode_stop = c("fixed_point", "certificate")  # or both
)
```

Return diagnostics (proposed): `B_obs`, `B_cert`, `V_cert`, `tau_B`, `w_max_lb`,
`mode_stop_pass`. See `R/c05_em.R` for current `tol_eff`, `mc_delta_floor`,
`icm_screen`, and `em_route`.

**Sharpest certified certificate** (*optimality display in **Proposition R-Cert***; smallest feasible
\(\widetilde C_d\), \(\gamma_0=\gamma^\star\), profile limit \(d\downarrow0\)):

\[
\boxed{\small
\begin{aligned}
\bigl\|P_\gamma^{k}(\gamma^\star,\cdot)-\pi_\gamma\bigr\|_{TV}
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
\alpha^{k} \\[4pt]
&\quad\quad\quad+\;\pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr),
\end{aligned}
}
\]

with \(\pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr)\le\delta_2\).

(At \(\gamma_0=\gamma^\star\) the two drift numerators coincide because
\(V(\gamma^\star)=V_{\sup}(d)=1\) in the \(d\downarrow0\) limit; for finite \(d>0\) use
the general box above with \(V_{\sup}(d)\) and
\(V(\gamma_0)=1+\tfrac12\|\gamma_0-\gamma^\star\|_{P_{11}}^2\).)

**Impact and limits of \(d\) and \(\delta_2\).**

§3.1 bounds \(\|P_\gamma^k-\pi_\gamma\|_{TV}\) for the **\(\beta\)-restricted**
sampler (§0), with \(\pi_\beta(\widetilde B(\delta_2)^{\,c})\) on the RHS and
\(\pi_\beta(\widetilde B(\delta_2)^{\,c})\le\delta_2\). The **dynamic** Rosenthal part
targets \(\pi_{\gamma\mid\widetilde B}\); the additive term is the **\(\beta\)-truncation**
gap \(\|\pi_{\gamma\mid\widetilde B}-\pi_\gamma\|_{TV}\le\pi_\beta(\widetilde B(\delta_2)^{\,c})\).

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
5. **Implementation** (`BETA_MARGINAL_MODE_LEVELSET.md` §7): eight-file API —
   **`beta_marginal_mode()`**, **`beta_marginal_safe_set()`**, **`group_precision_floor()`**,
   **`floor_coupling_eigenvalues()`**, **`rosenthal_drift_constants()`**, **`rosenthal_tv_bound()`**,
   **`optimal_rosenthal_tv_bound()`**, **`gamma_beta_tv_certificate()`**.
6. **Uniformity in \(\gamma'\) for \(\varepsilon^{\mathrm{LCB}}\)** (**Remark 5**): the
   stated interval is pointwise at the fitted \(\hat\gamma'\). Certifying
   \(\inf_{\gamma'}\) requires a deviation bound holding uniformly over \(\gamma'\),
   which in turn needs a certified search radius; the summands
   \(\varphi_m(\gamma')/q_Q(\gamma')\) are bounded, so a Bernstein-type argument over a
   compact \(\gamma'\)-region is the natural route.

---

## 9. References

- `inst/restricted_gibbs_minorization _v4.md` — C05 minorization, Theorem 1
- `inst/GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md` — general MT/Rosenthal catalog
- `inst/JOINT_GAMMA_BETA_TV_CERTIFICATE.md` — \(\widetilde B(\delta_2)\), product region
- `inst/BETA_MARGINAL_MODE_LEVELSET.md` — marginal \(\beta^\dagger\), \(r_{\mathrm{Gauss}}\), asymptotics
- `inst/LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` — scalar worked Foster proof
- `inst/CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` — matrix \(A\), \(\kappa_{\max}\)
- Rosenthal (1995); Meyn & Tweedie (1993)
