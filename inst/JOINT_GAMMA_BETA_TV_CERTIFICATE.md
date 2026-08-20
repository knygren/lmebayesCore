# Joint \(\gamma\)–\(\beta\) restriction for usable TV bounds

Companion to `vignette("Chapter-C05")`, `inst/CHAPTER_C05_IMPLEMENTATION.md`,
`inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md`, `inst/GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md`,
and `inst/SAFE_UNSAFE_TV_DECOMPOSITION.md`.

**Status:** draft design note — not implemented, not cited as proved.

---

## 0. Motivation

Chapter C05 produces a **restricted** two-block Gibbs certificate on a
\(\gamma\)-set \(\widetilde C_d\):

\[
q_n(\gamma,\cdot\mid \widetilde C_d)\ \text{ minorizes toward }\ \pi_\gamma
\qquad\text{with multiplier }\ \varepsilon = \varepsilon_d\,Q(\widetilde C_d),
\]

and Theorem 2 gives

\[
\bigl\|q_n(\cdot\mid \widetilde C_d)-\pi_\gamma\bigr\|_{TV}
\ \le\ (1-\varepsilon)^n + \delta,
\]

where \(\delta\) budgets escape from \(\widetilde C_d\) under \(\pi_\gamma\).

**This does not, by itself, control total variation for the full joint chain
\((\gamma,\beta)\)** on the original state space, for three independent reasons:

1. **The certified set is \(\gamma\)-only.** Block 1 (\(\beta\)) is never
   restricted in the shipped C05 pipeline. Curvature objects \(V_j(\gamma,\beta)\),
   \(S(\gamma,\beta)\), and \(\tilde J\) depend on the **current** \(\beta\)
   through the E-step, not only on \(\gamma^\star\).

2. **\(\varepsilon\) can be numerically tiny while \(\widetilde C_d\) is
   "reasonable" in \(\gamma\).** When coupling eigenvalues \(\kappa_i\) are
   near 1, the \(B\)-weights \(w_i=\kappa_i/(1-\kappa_i)\) are large, the
   calibrated \(d\) is large, and \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\)
   collapses — even if \(\varepsilon(\gamma^\star)\) itself is moderate. Minorization
   on \(\widetilde C_d\) alone then yields no practical sweep count.

3. **Off \(\widetilde C_d\), the chain is uncertified.** Production samplers
   do not confine \(\gamma\) to \(\widetilde C_d\). Any TV claim for the
   **unrestricted** kernel requires a separate escape term for \(\widetilde C_d^{\,c}\)
   *and* control of pathological \(\beta\) regions where precision floors fail.

**Conclusion.** A usable general TV bound likely requires a **joint**
restriction \((\gamma,\beta)\in \widetilde C_d \times \widetilde B(\delta_2)\),
with \(\delta_2\) budgeting \(\beta\)-mass outside a compact convex set on which
**data precision is uniformly bounded below**. Minorization on \(\gamma\) alone
is one ingredient, not the full certificate.

---

## 1. Notation (aligned with existing inst notes)

Joint state \(x=(\gamma,\beta)\) with \(\beta=(\beta_1,\ldots,\beta_J)\).

| Symbol | Meaning |
|--------|---------|
| \(\pi(\gamma,\beta\mid y)\) | Joint posterior |
| \(\pi_\gamma\), \(\pi_\beta\) | Marginals |
| \(\gamma^\star\) | C05 population mode (EM fixed point of conditional means) |
| \(\beta^\star\) | **Marginal posterior mode** of \(\beta\) (§3) |
| \(P_{11}(\beta)\), \(P_{22}(\beta)\) | Block refresh / conditional precisions at fixed \(\beta\) |
| \(H_j\), \(P_b\), \(V_j(\gamma,\beta)\) | C05 hyper-design and conditional covariances |
| \(S(\gamma,\beta)=P_{11}(\beta)\tilde J(\gamma,\beta)\) | Coupling matrix |
| \(\kappa_i(\gamma,\beta)\) | Symmetrized coupling eigenvalues (§4A) |
| \(\widetilde C_d\) | Certified \(\gamma\)-set, \(\{\Psi(\gamma)\le d\}\) |
| \(\widetilde B(\delta_2)\) | Certified \(\beta\)-set (§4) |
| \(\delta\) | \(\gamma\)-escape budget, \(\pi_\gamma(\widetilde C_d^{\,c})\) |
| \(\delta_2\) | \(\beta\)-escape budget, \(\pi_\beta(\widetilde B(\delta_2)^{\,c})\) |

Write observation-level weights \(w_{j,i}(\beta)\) for the GLM curvature
(e.g. \(p(1-p)\) for logit), so the group data contribution to the
\(\beta\)-block precision is

\[
\mathcal P_{j,\mathrm{data}}(\beta_j)
= D_j^\top \operatorname{diag}\!\big(w_{j,i}(\beta)\big)\, D_j .
\]

The **marginal** \(\gamma\)-precision at fixed \(\beta\) is

\[
\Pi(\gamma,\beta) = P_{11}(\beta) - S(\gamma,\beta),
\qquad
P_{11}(\beta) = \Lambda_\gamma + \sum_j H_j^\top P_b H_j .
\]

---

## 2. Identifiability and estimability (preflight layer)

> **Assumption H-id (rank and estimability).**
> For every group \(j\):
> - \(D_j\) has full column rank \(p_j\);
> - \(W_j\) (population design for \(\gamma\)) has full column rank \(q_j\);
> - the stacked population design has full rank \(q\);
> - for non-Gaussian families, the group-wise MLE \(\hat\beta_j\) exists and is
>   finite (no complete/quasi-complete separation; no degenerate sufficient
>   statistic).
>
> These are the same screens already enforced by `check_identifiability()` /
> `PREFLIGHT_model_setup.md` before sampling.

Under **H-id**:

- The joint posterior is proper on \(\mathbb R^q \times \mathbb R^{\sum_j p_j}\)
  (modulo the usual flat-prior caveats for \(\gamma\)).
- Group and population effects are **estimable** in the Fisher / IRLS sense at
  interior points.
- Local Laplace approximations at \(\hat\beta_j\) are valid starting points.

**What H-id does *not* give:** uniform curvature bounds on all of
\(\mathbb R^{\sum_j p_j}\). That requires a **restricted \(\beta\)-set** (§4).

---

## 3. Existence of a marginal \(\beta\)-mode

> **Proposition 3.1 (mode existence, conditional on H-id).**
> Suppose the log-likelihood is upper semicontinuous, coercive on each
> \(\beta_j\)-fiber after integrating out \(\gamma\) is not required — work on the
> **joint** log-posterior \(\log\pi(\gamma,\beta\mid y)\).
> If \(\log\pi(\gamma,\beta\mid y)\to -\infty\) as \(\|(\gamma,\beta)\|\to\infty\)
> along every ray (guaranteed for canonical GLMMs with proper priors on \(\Psi\)
> under H-id), then a **joint** maximizer \((\gamma^\dagger,\beta^\dagger)\) exists.
> If additionally the \(\beta\)-marginal is unimodal (strong log-concavity in
> \(\beta\) given \(\gamma\), or verified numerically), define
> \[
> \beta^\star \in \operatorname{argmax}_{\beta}\ \log\pi_\beta(\beta\mid y).
> \]

**Operational anchor.** In practice:

- \(\gamma^\star\) from `population_mode()` (C05 EM at conditional means);
- \(\beta^\dagger\) from Newton on the \(\gamma\)-integrated marginal
  \(\widetilde\pi(\beta\mid y)\) (`BETA_MARGINAL_MODE_LEVELSET.md`; contrast
  joint ICM / `glmerb_posterior_mode()` for a working point only).

The certificate is anchored at \((\gamma^\star,\beta^\dagger)\) for
\(\widetilde B(\delta_2)\) and \(\Xi\); matrix blocks may still use
\((\gamma^\star,\beta^{\mathrm{ICM}}_\star)\) until unified.

---

## 4. The \(\beta\)-safe set \(\widetilde B(\delta_2)\)

### 4.1 Existence of a mass-calibrated compact convex set

> **Proposition 4.1 (tail calibration).**
> Let \(\pi_\beta\) be a probability measure on \(\mathbb R^d\).
> For every \(\delta_2\in(0,1)\), there exists a compact convex set
> \(K_{\delta_2}\) with \(\pi_\beta(K_{\delta_2}^{\,c})\le \delta_2\)
> (take a large enough ball).
>
> **Goal:** construct \(K_{\delta_2}\) *algorithmically* so that a **uniform
> data-precision floor** holds on \(K_{\delta_2}\), not merely mass control.

Two constructive routes (compatible with existing inst notes):

**(C) Marginal mode + \(r_{\mathrm{Gauss}}\) on true \(\Xi\) (recommended).**
Integrate \(\gamma\) from the prior; find the **marginal mode**
\(\beta^\dagger=\arg\max_\beta \widetilde\pi(\beta\mid y)\); define
\(\Xi(\beta)=f(\beta)-f(\beta^\dagger)\) and
\[
\widetilde B(\delta_2)
=\bigl\{\beta:\ \Xi(\beta)\le r_{\mathrm{Gauss}}(n,\delta_2)\bigr\},
\qquad
r_{\mathrm{Gauss}}(n,\delta_2)=\tfrac12\chi^2_{n,\,1-\delta_2},
\ n=Jp_{\mathrm{re}}.
\]
Under Laplace / BvM (\(n_j\to\infty\)),
\(\pi_\beta(\widetilde B(\delta_2)^{\,c})\to\delta_2\). Full design:
`inst/BETA_MARGINAL_MODE_LEVELSET.md`; helpers in `R/c05_beta_marginal_set.R`.

**(A) Sublevel set of a convex profile (Gaussian-majorization route).**
When \(\Xi(\beta)\) is convex with \(\Xi(\beta^\dagger)=0\) and
\(\Xi(\beta)\to\infty\) as \(\|\beta\|\to\infty\) (see
`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` §7), set

\[
\widetilde B_r := \{\beta:\ \Xi(\beta)\le r\}.
\]

Then \(\widetilde B_r\) is compact and convex, \(\widetilde B_r\uparrow\mathbb R^d\)
as \(r\uparrow\infty\). **Do not** calibrate \(r\) by inverting Prop (P2) if a
tight tail is needed; prefer route **(C)** or MC validation on a fixed \(r\).

**(A′) Legacy (P2) calibration.** Choose \(r=r(\delta_2)\) so
\(\pi_\beta(\widetilde B_{r(\delta_2)}^{\,c})\le\delta_2\) via
\(\Gamma(n,r)/\gamma(n,r)\) — valid but loose at fixed \(n\).

**(B) Precision-floor safe region (safe/unsafe route).**
Fix \(\kappa_j>0\). Define

\[
\widetilde B_{\mathrm{safe}}
:= \bigcap_j \Bigl\{\beta_j:\ \lambda_{\min}\bigl(H_j(\beta_j)\bigr)\ge \kappa_j\Bigr\}.
\]

This is convex (sublevel of convex functions) but **not** automatically
mass-calibrated; intersect with a ball of radius \(R(\delta_2)\) chosen so
\(\pi_\beta(\widetilde B_{\mathrm{safe}}^c)\le \delta_2\).

**Definition 4.2 (certified \(\beta\)-set).**

\[
\boxed{\ \widetilde B(\delta_2)\ :=\ \text{compact convex set anchored at }\beta^\star
\ \text{with}\ \pi_\beta\bigl(\widetilde B(\delta_2)^{\,c}\bigr)\le \delta_2\ }
\]

and, crucially,

\[
\inf_{\beta\in \widetilde B(\delta_2)} w_{j,i}(\beta)\ \ge\ \omega_{j,i}>0
\qquad\text{for all observations } (j,i).
\]

### 4.2 Uniform data-precision lower bound

> **Lemma 4.3 (PD data precision on \(\widetilde B(\delta_2)\)).**
> Under H-id and the weight floor
> \[
> w_{j,i}(\beta)\ \ge\ \omega_{j,i}>0
> \qquad\forall\beta\in \widetilde B(\delta_2),
> \]
> define
> \[
> \underline{\mathcal P}_{j,\mathrm{data}}
> := D_j^\top \operatorname{diag}(\omega_{j,i})\, D_j .
> \]
> Then \(\underline{\mathcal P}_{j,\mathrm{data}}\succ 0\) (full rank) and
> \[
> \mathcal P_{j,\mathrm{data}}(\beta_j)\ \succeq\ \underline{\mathcal P}_{j,\mathrm{data}}
> \qquad\forall\beta\in \widetilde B(\delta_2).
> \]
> Consequently the **block-2 precision**
> \[
> P_{22,j}(\beta_j) = P_b + \mathcal P_{j,\mathrm{data}}(\beta_j)
> \ \succeq\ P_b + \underline{\mathcal P}_{j,\mathrm{data}}
> \ =:\ P_{22,j}^{\mathrm{LB}} .
> \]

This is the input required by `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` and by
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` for logit/probit/Poisson floors.

> **Lemma 4.4 (marginal \(\gamma\)-precision floor on \(\widetilde B(\delta_2)\)).**
> With \(P_{22,j}^{\mathrm{LB}}\) as above, define
> \[
> B^{\mathrm{UB}} := \sum_j H_j^\top P_b\bigl(P_{22,j}^{\mathrm{LB}}\bigr)^{-1} P_b H_j,
> \qquad
> \Pi^{\mathrm{LB}} := P_{11} - B^{\mathrm{UB}} .
> \]
> If \(\Pi^{\mathrm{LB}}\succ 0\), then for every \(\beta\in \widetilde B(\delta_2)\)
> and every \(\gamma\),
> \[
> \Pi(\gamma,\beta)\ \succeq\ \Pi^{\mathrm{LB}} .
> \]

Thus **on the product set** \(\widetilde C_d \times \widetilde B(\delta_2)\),
the \(\gamma\)-curvature used in C05 cannot collapse merely because \(\beta\)
wandered into a separation tail — provided \(\beta\) stays in
\(\widetilde B(\delta_2)\).

### 4.3 Uniform coupling ceiling on \(\widetilde B(\delta_2)\)

Write \(\kappa_{\max}^+(\delta_2)\) for the **worst-case** coupling ceiling over the
certified \(\beta\)-set:

\[
\kappa_{\max}^+(\delta_2)
\ :=\
\sup_{\beta\in \widetilde B(\delta_2)}\ \max_i \kappa_i(\gamma^\star,\beta),
\qquad
w_{\max}^+(\delta_2)\ :=\ \frac{\kappa_{\max}^+(\delta_2)}{1-\kappa_{\max}^+(\delta_2)}.
\]

On \(\widetilde B(\delta_2)\), Lemma 4.3–4.4 bound \(V_j(\gamma,\beta)\preceq
(P_{22,j}^{\mathrm{LB}})^{-1}\) uniformly, so \(S(\gamma,\beta)\) and
\(\tilde J(\gamma,\beta)\) are **uniformly bounded** in operator norm. Under H-id
and a positive weight floor, \(\kappa_{\max}^+(\delta_2)\) is typically **strictly
smaller** than the global supremum over all \(\beta\in\mathbb R^{\sum p_j}\)
(where separation drives \(\kappa\to 1\)).

**Consequence for sizing.** Calibrate \(d(\delta)\) from the \(B\)-spectrum at
\((\gamma^\star,\beta)\) with \(\beta\in \widetilde B(\delta_2)\) (worst case or
fixed anchor). For the **same** \(\delta\),

\[
d(\delta,\,\widetilde B(\delta_2))
\ \le\
d(\delta,\,\text{unrestricted }\beta),
\qquad\text{with strict inequality when }\kappa_{\max}^+\ll 1.
\]

Since \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\), a **smaller** certified
\(d\) yields a **larger** kernel floor \(\varepsilon_d\) on \(\widetilde C_d\).

---

## 5. Drift toward \(\widetilde C_d\) and geometric ergodicity

This section records why \(\widetilde B(\delta_2)\) is not only a TV bookkeeping
device but also the missing hypothesis for **Foster–Lyapunov drift** of the
\(\gamma\)-block toward the C05 certified set.

### 5.1 What fails without \(\beta\) restriction

`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §8–§11 (retained as a cautionary draft) attempted
Foster–Lyapunov drift on the \(\gamma\)-chain using a fixed refresh center and
IRLS weights evaluated at \(\hat\beta_j\). That argument is **not valid** for
logit/GLM when the unsafe set \(A^c\neq\varnothing\): as \(\beta\) enters a
separation tail, observation weights vanish, local curvature collapses, and no
uniform drift constant exists on all of \(\mathbb R^{\sum p_j}\).

The C05 deficiency gap

\[
\Psi(\gamma)
=\log\frac{\varepsilon(\gamma^\star)}{\varepsilon(\gamma)}
=\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}-\bar\Phi(\gamma)
\]

depends on \(\beta\) through \(\bar\Phi(\gamma)=\Phi(\gamma)-\Phi(\gamma^\star)\)
and through \(V_j(\gamma,\beta)\) in \(S(\gamma,\beta)\). Without confining
\(\beta\), the vector field \(\nabla\Psi(\gamma)=P_{11}\bigl(M(\gamma)-\gamma^\star\bigr)\)
(Chapter C05 §7) does not see a uniform log-concavity or precision floor.

### 5.2 Rosenthal and Meyn–Tweedie conditions (reference)

Full statements, notation conventions, and the relation to C05 Theorem 1 are in
**`inst/GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md`**. Summary:

**Minorization.** Small set \(C\), \(\varepsilon\in(0,1]\), probability \(\nu\):
\(P(x,\cdot)\ge \varepsilon\,\nu(\cdot)\) for all \(x\in C\).

**Geometric drift.** Lyapunov \(V:\mathcal X\to[1,\infty)\), \(\lambda\in(0,1)\),
\(b<\infty\):
\[
(PV)(x)\ \le\ \lambda V(x) + b\,\mathbb I_{C^c}(x).
\]
(Standard convention: constant \(b\) on the **complement** of \(C\), not on \(C\).)

**Rosenthal (1995) TV bound** (both hypotheses): for the **\(\beta\)-restricted**
\(\gamma\)-marginal chain \(P_{\gamma\mid\widetilde B(\delta_2)}\) to
\(\pi_\gamma(\cdot\mid\widetilde B(\delta_2),y)\),
\[
\bigl\|P_\gamma^{k}(\gamma_0,\cdot)-\pi_{\gamma\mid\widetilde B}\bigr\|_{TV}
\ \le\
(1-\varepsilon)^{\lfloor \alpha k\rfloor}
\;+\;
\frac{U}{\alpha}
\left(\frac{1+2b+\lambda V(\gamma_0)}{1+2b/(1-\lambda)}\right)\alpha^{k},
\qquad \alpha\in(\lambda,1),
\]
with \(U=1+2b+\lambda\sup_{x\in C} V(x)\). Add \(\pi_\beta(\widetilde B^{\,c})\le\delta_2\)
for the full \(\pi_\gamma\) (see `GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` §0).

**Meyn–Tweedie alternative:** \(\|P^k(x_0,\cdot)-\pi\|_{TV}\le M V(x_0)\rho^k\),
with \(\rho=\max\bigl((1-\varepsilon)^{1/m_0},\,\alpha\bigr)\) combining both
mechanisms in one rate.

**C05 Theorem 1** (`restricted_gibbs_minorization _v4.md` §2.4) is the
**minorization-only** special case
\(\|P_n(x,\cdot\mid C_d)-\pi\|_{TV}\le (1-\varepsilon)^n+\pi(C_d^c)\): first term
= Doeblin on \(C_d\); second term = **truncation** \(\delta=\pi(C_d^c)\), not Foster
drift.

### 5.3 Certificate instantiation

On \(\widetilde{\mathcal R}=\widetilde C_d\times \widetilde B(\delta_2)\):

| Object | Choice |
|--------|--------|
| Small set \(C\) | \(\widetilde C_d=\{\Psi(\gamma)\le d\}\) (× \(\widetilde B\) in product kernel) |
| \(V\ge 1\) | \(V(\gamma)=1+\tfrac12\|\gamma-\gamma^\star\|_{P_{11}}^2=\Psi+\bar\Phi\) |
| \(\lambda(\delta_2)\) | \(\bigl(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\bigr)^2\) from floor spectrum on \(\widetilde B(\delta_2)\) |
| \(b(\delta_2)\) | \(1-\lambda(\delta_2)+q/2+C_\beta^{+}(\delta_2)\), \(C_\beta^{+}(\delta_2)=\tfrac12\sum_i\kappa_i^{\mathrm{LB}}(\delta_2)\) |
| \(\varepsilon\) (Rosenthal / unrestricted \(P_\gamma\)) | \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\) (Lemma 17(a)) |
| \(\varepsilon_{\mathrm{rest}}\) (`certificate()`) | \(\varepsilon_d\,Q(\widetilde C_d)\) — restricted chain vs \(Q_{\widetilde C_d}\) |
| \(\pi(C^c)\) | \(\delta=\pi_\gamma(\widetilde C_d^{\,c})\) (+ \(\delta_2\) for joint TV) |

**Lyapunov vs certified set.** Minorization is tied to \(\widetilde C_d=\{\Psi\le d\}\).
Drift uses the **refresh quadratic** \(V=1+\tfrac12\|\cdot-\gamma^\star\|_{P_{11}}^2\)
(Chapter C03 Claim 2 metric), not \(\Psi\) alone. Since \(\bar\Phi\ge 0\),
\(\Psi\le V_0\) and \(\{\|\gamma-\gamma^\star\|_{P_{11}}^2\le 2d\}\subseteq\widetilde C_d\)
(Corollary 15).

> **Remark 5.1 (compact \(\widetilde B\)).** Uniform \(\lambda<1\) and finite \(b\)
> require compact \(\widetilde B(\delta_2)\) with weight floor
> \(\inf_{\widetilde B} w_{j,i}\ge\omega_{j,i}>0\). Without it,
> \(\kappa_{\max}\to 1\) and Foster drift fails (`RATE_Ac` §8).

> **Remark 5.2 (which bound when).** Use \((1-\varepsilon)^n+\delta\) for the shipped
> C05 certificate; add Rosenthal §3.1 or MT §4.1 when controlling **dynamic** escape
> to \(\widetilde C_d^{\,c}\) or starts outside \(C\). See
> `GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md` §0.

### 5.4 Shrinking \(\widetilde C_d\) and raising \(\varepsilon_d\)

Two distinct "small set" benefits appear once \(\beta\in \widetilde B(\delta_2)\):

**(A) Smaller calibration radius for fixed \(\delta\).** As in §4.3,
\(\kappa_{\max}^+(\delta_2)\) and \(w_{\max}^+(\delta_2)\) shrink relative to the
unrestricted worst case, so the \(B\)-spectrum inversion gives

\[
r(\delta,\,\widetilde B)^2\ \le\ r(\delta,\,\text{global})^2,
\qquad
d(\delta,\,\widetilde B)\ =\ \tfrac12 r(\delta,\,\widetilde B)^2
\]

**decreases** at the same tail budget \(\delta\). The certified \(\gamma\)-set
\(\widetilde C_d=\{\Psi\le d\}\) is therefore **smaller** (tighter minorization
region) while still meeting \(\pi_\gamma(\widetilde C_d^{\,c})\lesssim\delta\) under
the Gaussian-reference calibration.

**(B) Larger kernel floor on that set.** Since \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\),

\[
d\ \downarrow
\quad\Longrightarrow\quad
\varepsilon_d\ \uparrow
\]

on the **unrestricted** marginal kernel (Rosenthal minorization rate). The shipped
`certificate()` constant
\(\varepsilon_{\mathrm{rest}}=\varepsilon_d\,Q(\widetilde C_d)\) also rises when
\(d\downarrow\) (modulo the mild \(Q\)-mass increase — usually secondary).

**(C) Drift pulls mass into the smaller set.** Foster drift on
\(V=1+\tfrac12\|\cdot-\gamma^\star\|_{P_{11}}^2\) (§5.3) shows the **chain** also
concentrates toward \(\widetilde C_d\), not merely the target \(\pi_\gamma\).

**Design implication.** One should choose \(\delta_2\) first (fix \(\widetilde B\),
compute \(\kappa_{\max}^+(\delta_2)\), **then** calibrate \(d(\delta)\) from the
\(B\)-spectrum at that \(\beta\)-safe curvature — not from a worst-case \(\beta\)
over all of \(\mathbb R^{\sum p_j}\). This ordering aligns minorization, drift,
and geometric ergodicity on the same product set \(\widetilde{\mathcal R}\).

---

## 6. Joint restricted certificate (sketch)

Define the **joint certified region**

\[
\widetilde{\mathcal R}
:= \widetilde C_d \times \widetilde B(\delta_2).
\]

Consider the **doubly restricted** Gibbs kernel: update \(\beta\) only when
\(\gamma\in \widetilde C_d\), update \(\gamma\) only when \(\beta\in \widetilde B(\delta_2)\),
(or restrict proposals to \(\widetilde{\mathcal R}\) after refresh).

**Stage 1 (C05, on \(\gamma\)).** At \(\beta\in \widetilde B(\delta_2)\), compute
\(\kappa_i(\gamma^\star,\beta)\) with uniform ceiling \(\kappa_{\max}^+(\delta_2)\)
(§4.3), calibrate \(d(\delta)\) from the \(B\)-spectrum, and form \(\varepsilon_d\),
\(\varepsilon\) as in `CHAPTER_C05_IMPLEMENTATION.md`. Prefer **worst-case**
\(\beta\in \widetilde B(\delta_2)\) for \(d\) only if it remains certifiable;
otherwise anchor at \(\beta^\star\).

**Stage 2 (new, on \(\beta\)).** At fixed \(\gamma\in \widetilde C_d\), use
Lemma 4.3–4.4 to certify a **Block-1 minorization / majorization** constant on
\(\widetilde B(\delta_2)\) (exact for Gaussian blocks; envelope for GLM — see
`MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md`, `LOGIT_SINGLE_GROUP_SAFE_REGION.md`).

**Stage 3 (TV assembly).** Decompose joint TV as in
`SAFE_UNSAFE_TV_DECOMPOSITION.md`:

\[
\|K^n(x,\cdot)-\pi\|_{TV}
\ \le\
\underbrace{d_{\widetilde{\mathcal R}}(x)}_{\text{restricted convergence}}
\ +\
\underbrace{K^n(x,\widetilde{\mathcal R}^{\,c})}_{\text{chain escape}}
\ +\
\underbrace{\pi(\widetilde{\mathcal R}^{\,c})}_{\text{target escape}} .
\]

Split the target escape:

\[
\pi(\widetilde{\mathcal R}^{\,c})
\ \le\
\pi_\gamma(\widetilde C_d^{\,c})
\ +\
\pi_\beta(\widetilde B(\delta_2)^{\,c})
\ +\
\pi\bigl(\widetilde C_d \times \widetilde B(\delta_2)^{\,c}\bigr)
\ \le\
\delta + \delta_2 + \delta_{12},
\]

where \(\delta_{12}\) is the **cross-term** ( \(\gamma\) inside, \(\beta\) outside).
Bounding \(\delta_{12}\) requires either a product-structure assumption or a
Gaussian majorization on the joint tail — open (§7).

The **restricted** term \(d_{\widetilde{\mathcal R}}\) uses the **product**
minorization constant

\[
\varepsilon_{\mathrm{joint}}
\ \approx\
\varepsilon_\gamma(\widetilde C_d)\ \times\
\varepsilon_\beta(\widetilde B(\delta_2))
\]

(under independence of refresh blocks — to be made precise).

---

## 7. Why \(\gamma\)-only minorization fails in practice (summary table)

| Quantity | Depends on | Fixed at certificate time? | Without \(\widetilde B\)? |
|----------|------------|----------------------------|---------------------------|
| \(\varepsilon(\gamma^\star)\) | \(\tilde J\) at \((\gamma^\star,\beta_{\mathrm{E\text{-}step}})\) | Yes, one E-step | \(\beta\) in E-step must be specified |
| \(\kappa_i\), \(d\), \(\varepsilon_d\) | \(S(\gamma^\star,\beta)\) | Uses terminal E-step | \(\kappa\to 1\) if \(\beta\) in tail |
| \(Q(\widetilde C_d)\) | \(\kappa\)-ellipsoid | Yes | OK on \(\gamma\) margin |
| \(\varepsilon\) | \(\varepsilon_d Q(\widetilde C_d)\) | Yes | Can be \(\ll 1\) |
| Block-1 conditioning | \(w_{j,i}(\beta)\) | **No** | Weights \(\to 0\) in tails |
| \(\Pi^{\mathrm{LB}}\) | \(P_{22,j}^{\mathrm{LB}}\) | Requires \(\widetilde B\) | Fails off \(\widetilde B\) |

**Takeaway:** \(\widetilde C_d\) controls the **\(\gamma\)** refresh geometry;
\(\widetilde B(\delta_2)\) controls **likelihood curvature** feeding back into
\(\gamma\) through \(V_j(\beta)\) and \(P_{22}(\beta)\). Both are needed for a
joint TV bound that survives GLM tails.

---

## 8. Open problems (implementation-facing)

1. **Cross-term \(\delta_{12}\).** Bound
   \(\pi(\widetilde C_d \times \widetilde B^{\,c})\) without assuming
   \(\gamma\)-\(\beta\) independence.

2. **Anchor consistency.** Should the E-step at \(\gamma^\star\) use
   \(\beta^\star\), conditional means, or modes? C05 currently uses conditional
   means; \(\beta^\star\) may differ.

3. **Constructing \(\Xi(\beta)\)** for automatic \(\widetilde B(\delta_2)\) in
   every family — logit is done (`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` §7);
   Poisson-log is not.

4. **Computability of \(\omega_{j,i}\)** on \(\widetilde B(\delta_2)\): support
   function of \(\widetilde B_r\) vs grid over groups (`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md`).

5. **Production vs certified chain.** The doubly restricted kernel is a proof
   device; mapping bounds to `rlmerb` / `rglmerb` requires the safe/unsafe
   split in `SAFE_UNSAFE_TV_DECOMPOSITION.md`.

6. **Choosing \((\delta,\delta_2)\).** Trade-off on the \(\beta\) side: **larger**
   \(\delta_2\) allows a tighter \(\widetilde B(\delta_2)\) with higher weight floors
   \(\omega_{j,i}(\delta_2)\), hence **lower** \(\kappa_i^{\mathrm{LB}}(\delta_2)\) and
   **better** Rosenthal drift; **smaller** \(\delta_2\) inflates \(\widetilde B\) to meet
   \(\pi_\beta(\widetilde B^{\,c})\le\delta_2\), typically **raising** \(\kappa_{\max}^{\mathrm{LB}}(\delta_2)\)
   and hurting drift. On the \(\gamma\) side, smaller certified \(d\) (given \(\delta_2\))
   still **raises** \(\varepsilon(d)\) in Rosenthal §3.1. Calibrate \(\delta_2\) before
   \(\delta\) (§5.4), accepting the \(\beta\)-mass vs drift-constant trade-off.

7. **Formal drift proof on \(\widetilde{\mathcal R}\).** Turn §5.3 into a lemma with
   explicit \((\lambda,b)\) in terms of \(\Pi^{\mathrm{LB}}\),
   \(\kappa_{\max}^+(\delta_2)\), and refresh constants; cite Rosenthal (1995) or
   `GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md` for the assembled TV bound.

---

## 9. Proposed pipeline (not shipped)

1. **Preflight:** H-id via `check_identifiability()`.
2. **Modes:** \(\gamma^\star\) (`population_mode`), \(\beta^\star\) (TBD API).
3. **\(\beta\)-set:** build \(\widetilde B(\delta_2)\), verify weight floors
   \(\omega_{j,i}\), certify \(\pi_\beta(\widetilde B^{\,c})\le \delta_2\); record
   \(\kappa_{\max}^+(\delta_2)\).
4. **\(\gamma\)-set:** `deficiency_spectrum` / `deficiency_calibrate` at
   \((\gamma^\star,\beta)\) on \(\widetilde B(\delta_2)\) — yields smaller
   \(d(\delta)\) and larger \(\varepsilon_d\) than unrestricted \(\beta\) (§5.4).
5. **Drift check (optional):** verify Foster–Lyapunov for
   \(V=1+\tfrac12\|\gamma-\gamma^\star\|_{P_{11}}^2\) on \(\widetilde{\mathcal R}\)
   (§5.3; `GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md` §2.2).
6. **Joint constants:** \(\varepsilon_\gamma\), \(\varepsilon_\beta\),
   \(\varepsilon_{\mathrm{joint}}\).
7. **TV / ergodicity bound:** C05 \((1-\varepsilon)^n+\delta\) plus Rosenthal/MT
   when dynamic escape matters (§5.2); \(\delta + \delta_2 + \delta_{12}\) for joint
   escape.

---

## 10. References

- `inst/GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` — explicit \((\lambda,b,C_\beta,\varepsilon)\)
  and Rosenthal TV assembly for the \(\gamma\)-marginal chain.
- `inst/GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md` — general MT/Rosenthal catalog
- `inst/restricted_gibbs_minorization _v4.md` — Theorem 1 \((1-\varepsilon)^n+\pi(C^c)\).
- `inst/MVN_CALIBRATED_MINORIZATION_SET.md` — \(\varepsilon_d\) vs \(\varepsilon\).
- `inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` — \(\widetilde B_r\), \(P_{22}^{\mathrm{LB}}\).
- `inst/SAFE_UNSAFE_TV_DECOMPOSITION.md` — partition / three-term TV.
- `inst/RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` — per-family weight floors; §8 drift
  **invalid** without \(\widetilde B\) (contrast §5.1 here).
- `inst/BLOCK_GIBBS_ERGODICITY.md` — geometric contraction in the Gaussian case.
- `inst/HIERARCHICAL_GENERALIZED_LINEAR_MODEL_NOTATION.md` — \(D_j\), \(W_j\), weights.
- `inst/PREFLIGHT_model_setup.md` — estimability gates.
