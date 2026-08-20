# Laplace–profile group marginal bound for certified precision floors

Draft companion to `inst/LOGIT_MARGINAL_INTEGRATE_GAMMA.md` (integrate
\(\gamma\)), `inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` (Proposition 2,
level sets, support functions), and `R/group_precision_floor.R` (implementation
sketch).

**Goal.** Certify per-group data-precision lower bounds \(\Gamma_j^{\mathrm{LB}}\)
from a total escape budget \(\epsilon\), using an **approximate marginal in
\(\beta_j\)** that is **not** multivariate normal. The Gaussian structure lives
in the **other block** \(\beta_{-j}\) (all coordinates except group \(j\)),
conditional on \(\beta_j\). The marginal in \(\beta_j\) acquires a
\(\tfrac12\log\det\) **volume factor** from integrating that conditional.

This note states **Assumption (L)**, derives the approximate marginal
\(\widetilde m_j\), bounds its outer tail with \(\epsilon/J\), and connects the
certified \(\beta_j\)-set to scalar weight floors \( \underline w_{ij}\).

---

## 1. Setup after integrating \(\gamma\)

Stack \(\beta = (\beta_1^\top,\ldots,\beta_J^\top)^\top \in \mathbb R^n\) with
\(n = J p_{\mathrm{re}}\). After integrating population parameters \(\gamma\)
(exact Gaussian step; see `LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §4),

\[
\widetilde\pi(\beta \mid y)
\;\propto\;
\exp\Bigl(\sum_{j=1}^J \ell_j(\beta_j)\Bigr)\,
\exp\Bigl(-\tfrac12(\beta - \mu_\beta)^\top \Lambda_\beta (\beta - \mu_\beta)\Bigr),
\]

with \(\Lambda_\beta\) the Schur-complement prior precision on \(\beta\).

Define the **joint deficiency** at the posterior mode \(\beta^\dagger\):

\[
\Xi(\beta) := f(\beta) - f(\beta^\dagger),
\qquad
f(\beta) := -\log \widetilde\pi(\beta \mid y) + \text{const},
\]

so \(\Xi\) is convex, \(\Xi(\beta^\dagger)=0\), and \(\Xi>0\) off the mode.
Write \(\beta = (\beta_j, \beta_{-j})\) with \(\beta_j \in \mathbb R^{p_{\mathrm{re}}}\)
and \(\beta_{-j} \in \mathbb R^m\), \(m = n - p_{\mathrm{re}}\).

For group \(j\), observations \(i=1,\ldots,n_j\), design rows \(z_{ij}^\top\),
weights \(w_{ij}(\eta) = -d^2\ell_{j,i}/d\eta^2\) at \(\eta = z_{ij}^\top\beta_j\)
(logit: \(w_{ij} = n_{ij} p(1-p)\); Poisson-log: monotone; etc.).

---

## 2. Profile and the conditional Gaussian in \(\beta_{-j}\)

### 2.1 Profile deficiency

\[
\boxed{
\Xi_j^{\mathrm{prof}}(\beta_j)
\;:=\;
\min_{\beta_{-j}} \Xi(\beta_j, \beta_{-j}),
}
\qquad
\hat\beta_{-j}(\beta_j) \in \arg\min_{\beta_{-j}} \Xi(\beta_j,\cdot).
\]

This is computed implicitly by the KKT system for support functions on
\(\{\Xi \le R\}\): when maximizing a functional supported on block \(j\), the
other blocks sit at \(\hat\beta_{-j}(\beta_j)\) (their conditional mode given
\(\beta_j\)). **Do not** freeze \(\beta_{-j}\) at the joint mode \(\beta_{-j}^\dagger\)
for certification — that is anti-conservative.

### 2.2 Assumption (L): Laplace–Gaussian in the other block

**Assumption (L).** For each fixed \(\beta_j\), the conditional law of
\(\beta_{-j}\) induced by \(\exp(-\Xi(\beta_j,\cdot))\) is approximated by

\[
\boxed{
\beta_{-j} \,\big|\, \beta_j
\;\approx\;
N\!\Bigl(\hat\beta_{-j}(\beta_j),\; H_{-j}(\beta_j)^{-1}\Bigr),
}
\]

where

\[
H_{-j}(\beta_j)
\;:=\;
\nabla^2_{\beta_{-j}\beta_{-j}}\,
\Xi\bigl(\beta_j, \hat\beta_{-j}(\beta_j)\bigr).
\]

This is the standard **Laplace** approximation: Gaussian in \(\beta_{-j}\), **not**
in \(\beta_j\). It is exact when \(\Xi(\beta_j,\cdot)\) is quadratic in
\(\beta_{-j}\); for logit it is local at the profile point.

### 2.3 Hessian sandwich (certified)

Write \(\bar G_{-j}\) for the block-diagonal matrix of maximal observation
Fisher curvatures on groups \(\neq j\) (logit: \((n_{k,i}/4)\, z_{ki} z_{ki}^\top\)
assembled into the \(-j\) block). Then on the profile manifold,

\[
\boxed{
\Lambda_\beta^{-j,-j}
\;\preceq\;
H_{-j}(\beta_j)
\;\preceq\;
\Lambda_\beta^{-j,-j} + \bar G_{-j}.
}
\]

This is the same uniform bound used for the crude \(\kappa\) constant in
`data-raw/_chk_group_marginal_bound.R`.

---

## 3. Approximate group marginal (not MVN in \(\beta_j\))

Integrate the conditional Gaussian over \(\beta_{-j}\):

\[
V(\beta_j)
\;:=\;
\int_{\mathbb R^m}
\exp\Bigl(
-\bigl[\Xi(\beta_j,\beta_{-j}) - \Xi_j^{\mathrm{prof}}(\beta_j)\bigr]
\Bigr)\, d\beta_{-j}.
\]

Under (L), \(V(\beta_j) \approx (2\pi)^{m/2} \det(H_{-j}(\beta_j))^{-1/2}\). The
**approximate marginal** in \(\beta_j\) alone is therefore

\[
\boxed{
\widetilde m_j(\beta_j)
\;\propto\;
\exp\!\Bigl(
-\Psi_j(\beta_j)
\Bigr),
\qquad
\Psi_j(\beta_j)
\;:=\;
\Xi_j^{\mathrm{prof}}(\beta_j)
\;+\;
\tfrac12 \log\det H_{-j}(\beta_j),
}
\]

(up to \(\beta_j\)-independent constants). **Key point:** \(\widetilde m_j\) is
**not** \(N(\beta_j^\dagger, S_j^{-1})\). The only MVN is in \(\beta_{-j}\);
\(\Psi_j\) adds a **log-volume** term from integrating that block.

The **true** marginal \(\widetilde\pi_j\) shares the profile factor
\(e^{-\Xi_j^{\mathrm{prof}}}\) but has a **different** volume factor \(V(\beta_j)\).
Certification compares sets via bounds on \(\log V\), not via \(r_{\mathrm{Gauss}}(p_{\mathrm{re}})\)
as if \(\beta_j\) were normal.

---

## 4. Bounding variation of the log-volume term

Only the **variation** of \(\tfrac12\log\det H_{-j}\) matters for level-set
inclusions (the absolute normalization cancels between ratio statements).

### 4.1 Certified upper bound (crude)

For any \(\beta_j^{(a)}, \beta_j^{(b)}\),

\[
\bigl|\log V(\beta_j^{(a)}) - \log V(\beta_j^{(b)})\bigr|
\;\le\;
\kappa_j^{\mathrm{cert}},
\]

with

\[
\boxed{
\kappa_j^{\mathrm{cert}}
\;:=\;
\tfrac12 \log\det\!\Bigl(
I + (\Lambda_\beta^{-j,-j})^{-1}\bar G_{-j}
\Bigr).
}
\]

This prices integrating a **uniformly** larger Gaussian over \(\beta_{-j}\)
(Hessian sandwich). It is **safe** but can be very loose when \(m = n - p_{\mathrm{re}}\)
is large (see `data-raw/_chk_group_marginal_bound.R`: mean crude \(\approx 18\)
vs Laplace \(\approx 0.4\) on the \(J=10\), \(p_{\mathrm{re}}=3\) example).

### 4.2 Diagnostic (Laplace along the extremal path)

\[
\kappa_j^{\mathrm{Lap}}
\;:=\;
\tfrac12
\Bigl|
\log\det H_{-j}(b_{\mathrm{ext}}) - \log\det H_{-j}(\beta^\dagger)
\Bigr|,
\]

where \(b_{\mathrm{ext}}\) is the \(\beta\) vector at a support-function solve on
\(\{\Xi \le R\}\). This matches the **actual** determinant ratio under (L) but
is **not** a proved uniform bound unless supplemented by a Lipschitz argument.

### 4.3 Profile–volume identity

At any \(\beta_j\),

\[
\Psi_j(\beta_j) - \Xi_j^{\mathrm{prof}}(\beta_j)
\;=\;
\tfrac12 \log\det H_{-j}(\beta_j).
\]

Hence for any two points on a path,

\[
\bigl|\Psi_j(\beta_j^{(a)}) - \Psi_j(\beta_j^{(b)})\bigr|
\;\le\;
\bigl|\Xi_j^{\mathrm{prof}}(\beta_j^{(a)}) - \Xi_j^{\mathrm{prof}}(\beta_j^{(b)})\bigr|
\;+\;
2\kappa_j,
\]

when \(|\tfrac12\Delta\log\det| \le \kappa_j\) on that path.

---

## 5. Tail budget: Proposition 2 on the profile part only

Proposition 2 (`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` §3A.4.1) applies to **any**
convex \(\Xi\) in dimension \(d\):

\[
\pi\bigl(\{\Xi > r\}\bigr)
\;\le\;
\frac{\Gamma(d,r)}{\gamma(d,r)}
\;=\;
\frac{\Pr(\mathrm{Gamma}(d,1)>r)}{\Pr(\mathrm{Gamma}(d,1)\le r)}.
\]

Apply it to **\(\Xi_j^{\mathrm{prof}}\)** in dimension \(d = p_{\mathrm{re}}\), with
**union bound** over \(j\):

\[
\boxed{
r_j
\;:=\;
r_{\mathrm{P2}}\!\Bigl(p_{\mathrm{re}},\,\frac{\epsilon}{J}\Bigr),
}
\]

the unique solution of

\[
\frac{\Gamma(p_{\mathrm{re}}, r_j)}{\gamma(p_{\mathrm{re}}, r_j)}
\;=\;
\frac{\epsilon}{J}.
\]

This controls outer mass for \(\exp(-\Xi_j^{\mathrm{prof}})\) only. It is a
**certified upper bound** on escape probability (convex worst case), **not** the
exact tail of \(\widetilde m_j\) or \(\widetilde\pi_j\).

**Contrast with \(r_{\mathrm{Gauss}}\).** For a **true** \(N(0,I_d)\) law in
\(\mathbb R^d\),

\[
r_{\mathrm{Gauss}}(d,\delta)
= q_{1-\delta}\bigl(\mathrm{Gamma}(d/2,1)\bigr)
= \tfrac12 \chi^2_{d,\,1-\delta}.
\]

At fixed \(\delta\), \(r_{\mathrm{P2}}(d,\delta) \approx r_{\mathrm{Gauss}}(2d,\delta)\):
Prop 2 charges **twice the effective dimension**. Neither \(r_{\mathrm{P2}}\) nor
\(r_{\mathrm{Gauss}}\) is the exact threshold for \(\widetilde m_j\) under (L);
§6 links them to \(\Psi_j\).

---

## 6. Inclusion (L′): from \(\Psi_j\) to a profile level set

Define the **certified profile superlevel set**

\[
\boxed{
\mathcal B_j(R)
\;:=\;
\bigl\{\beta_j : \Xi_j^{\mathrm{prof}}(\beta_j) \le R\bigr\}.
}
\]

**Volume inflation.** If \(|\tfrac12\log\det H_{-j}(\beta_j) - \tfrac12\log\det H_{-j}(\beta_j^\dagger)| \le \kappa_j\)
along the relevant path, then

\[
\Psi_j(\beta_j) \le R
\;\Longrightarrow\;
\Xi_j^{\mathrm{prof}}(\beta_j) \le R + \kappa_j,
\]

and conversely \(\Xi_j^{\mathrm{prof}}(\beta_j) \le R - \kappa_j\) implies
\(\Psi_j(\beta_j) \le R\) when the log-det is \(\ge -\kappa_j\) above its value at
the mode. For an **outer-tail** certificate on \(\widetilde m_j\), use the **larger**
profile set

\[
\boxed{
R_j := r_j + 2\kappa_j.
}
\]

Then on \(\mathcal B_j(R_j)^c\),

\[
\widetilde m_j(\beta_j)
\;\propto\;
e^{-\Psi_j(\beta_j)}
\;\le\;
e^{\kappa_j}\, e^{-\Xi_j^{\mathrm{prof}}(\beta_j)}
\quad\text{when } \tfrac12\log\det H \ge \tfrac12\log\det H^\dagger - \kappa_j.
\]

**Mass bound (approximate marginal under (L)).** If \(\kappa_j\) certifies
\(|\tfrac12\Delta\log\det| \le \kappa_j\) on \(\mathcal B_j(R_j)\) and
\(R_j \ge r_j + \kappa_j\), then

\[
\widetilde m_j\bigl(\mathcal B_j(R_j)^c\bigr)
\;\le\;
e^{\kappa_j}\,
\Pr\bigl[\Xi_j^{\mathrm{prof}} > r_j\bigr]_{\exp(-\Xi^{\mathrm{prof}})}
\;\le\;
e^{\kappa_j}\,\frac{\epsilon}{J},
\]

where the last step is Proposition 2 on \(\Xi_j^{\mathrm{prof}}\) in dimension
\(p_{\mathrm{re}}\). Equivalently, absorb the volume inflation into the budget:

\[
r_j
\;=\;
r_{\mathrm{P2}}\!\Bigl(
p_{\mathrm{re}},\,\frac{\epsilon}{J\,e^{\kappa_j}}
\Bigr),
\qquad
R_j = r_j + 2\kappa_j.
\]

When \(\kappa_j^{\mathrm{cert}}\) is large (crude sandwich), the per-group budget
\(\epsilon/(J e^{\kappa_j})\) is tiny and \(R_j\) may exceed the **joint**
\(r_{\mathrm{P2}}(n,\epsilon)\) — the per-group scheme can be worse than certifying
jointly (`data-raw/_chk_group_marginal_bound.R`). Tighter \(\kappa_j\) is the
main open improvement (§8).

---

## 7. Proposition L (statement)

> **Proposition L (Laplace–profile group marginal bound).**
> Assume (L) and the Hessian sandwich §2.3. Let \(\widetilde m_j \propto
> \exp(-\Psi_j)\) with \(\Psi_j = \Xi_j^{\mathrm{prof}} + \tfrac12\log\det H_{-j}\).
> For escape budget \(\epsilon \in (0,1)\), define \(\kappa_j\) as in §4.1
> (certified) or §4.2 (diagnostic), and
> \[
> r_j = r_{\mathrm{P2}}\!\bigl(p_{\mathrm{re}},\, \epsilon/(J e^{\kappa_j})\bigr),
> \qquad
> R_j = r_j + 2\kappa_j.
> \]
> Then under the approximate marginal \(\widetilde m_j\),
> \[
> \widetilde m_j\bigl(\mathcal B_j(R_j)^c\bigr) \le \epsilon/J.
> \]
> On \(\mathcal B_j(R_j)\), the logit (or GLM) weights satisfy
> \(w_{ij}(z_{ij}^\top\beta_j) \ge \underline w_{ij}\) with \(\underline w_{ij}\)
> the minimum of \(w_{ij}\) over the **image interval** of \(\eta_{ij} = z_{ij}^\top\beta_j\)
> on that set (§9).

*Proof sketch.* §6: volume factor bounded by \(e^{\kappa_j}\) on the relevant tail;
profile tail bounded by Prop 2 with deflated budget \(\epsilon/(J e^{\kappa_j})\);
set \(R_j = r_j + 2\kappa_j\) so \(\{\Psi_j > R_j\} \subseteq \{\Xi_j^{\mathrm{prof}} > r_j\}\)
up to the certified log-det variation. \(\blacksquare\)

**What is *not* claimed.** (i) \(\widetilde m_j = \widetilde\pi_j\); (ii) \(\beta_j\)
is multivariate normal; (iii) \(r_{\mathrm{Gauss}}(p_{\mathrm{re}})\) is the
correct level for \(\widetilde m_j\). The Gaussian piece is **only** in
\(\beta_{-j}\).

---

## 8. Weight floor on \(\mathcal B_j(R_j)\)

For each observation \((j,i)\), let \(\eta_{ij} = z_{ij}^\top \beta_j\). On
\(\mathcal B_j(R_j)\),

\[
[\eta_{ij}^-, \eta_{ij}^+]
\;=\;
\Bigl\{
z_{ij}^\top \beta_j : \beta_j \in \mathcal B_j(R_j)
\Bigr\},
\]

computed by **support functions** on the convex set \(\{\Xi \le R_j\}\) with
\(c^\top\beta\) supported on block \(j\) (KKT / Newton; same as `supp()` in
`data-raw/_ex_logit_floor_J10_p3.R`).

For logit, \(w_{ij}(\eta)\) is log-concave in \(\eta\), so

\[
\underline w_{ij}
\;=\;
\min_{\eta \in [\eta_{ij}^-, \eta_{ij}^+]} w_{ij}(\eta)
\;=\;
\min\bigl\{w_{ij}(\eta_{ij}^-),\, w_{ij}(\eta_{ij}^+)\bigr\}.
\]

Assemble the **certified data-precision lower bound**

\[
\boxed{
\Gamma_j^{\mathrm{LB}}
\;=\;
\sum_{i=1}^{n_j} \underline w_{ij}\, z_{ij} z_{ij}^\top.
}
\]

Optional contraction diagnostic (same as `group_precision_floor()`):

\[
\omega_j
\;=\;
\lambda_{\max}\!\Bigl(
(P_b + \Gamma_j^{\mathrm{LB}})^{-1} P_b
\Bigr),
\]

with \(P_b\) the within-group prior precision.

---

## 9. End-to-end pipeline (per group \(j\))

```text
  ε  (total escape budget)
       │
       ▼
  κ_j   ←  §4.1 certified (or §4.2 diagnostic)
       │
       ▼
  r_j = r_P2(p_re, ε / (J e^κ_j))     ←  Prop 2 on Ξ_j^prof only
       │
       ▼
  R_j = r_j + 2 κ_j                   ←  profile set for Ψ_j / weights
       │
       ▼
  support functions on {Ξ ≤ R_j}      ←  others at β̂_{-j}(β_j)
       │
       ▼
  Γ_j^LB = Σ_i w̲_ij z_ij z_ij'      ←  endpoint rule on η intervals
```

Implementation: `group_precision_floor()` in `R/group_precision_floor.R` follows
this pipeline with `kappa_method = "crude"|"laplace"|"none"` and
`.gamma_level_for_budget()` for \(r_{\mathrm{P2}}\). Refinements in §10 align
documentation with Proposition L (budget deflation \(e^{-\kappa_j}\), explicit
Assumption (L)).

---

## 10. Open improvements (same approximation class)

1. **Tighter certified \(\kappa_j\).** Replace §4.1 with a determinant-**ratio**
   bound using concavity of \(\log\det\) and \(\mathrm{tr}(H^{-1}\Delta G)\),
   combined with \(|w'| \le w\) for logit — target \(\kappa_j \approx 0.4\)–\(0.9\)
   instead of \(\approx 18\) on the \(J=10\) example.

2. **Strong convexity on \(\Xi_j^{\mathrm{prof}}\).** When data are informative,
   replace \(r_{\mathrm{P2}}\) by a \(\mathrm{Gamma}(p_{\mathrm{re}}/2)\)-type
   level (Gaussian growth along rays) on the profile part only.

3. **True vs approximate marginal.** Proposition L certifies \(\widetilde m_j\)
   under (L). Linking \(\widetilde\pi_j\) to \(\widetilde m_j\) requires a
   separate Laplace error bound (not automatic from Prop 2).

---

## 11. Numerical checks

| Script | Role |
|---|---|
| `data-raw/_chk_group_marginal_bound.R` | Per-group vs joint \(r\); crude vs Laplace \(\kappa\); floor at various \(R\) |
| `data-raw/_chk_group_precision_floor_pkg.R` | Package internals vs scratch fixture |
| `data-raw/_chk_hpd_density_level.R` | \(r_{\mathrm{Gauss}}(n,\epsilon)\) reference (not used in Prop L) |
| `data-raw/_ex_logit_floor_J10_p3.R` | Full logit example; profile coupling |

---

## References

- `inst/LOGIT_MARGINAL_INTEGRATE_GAMMA.md` — integrate \(\gamma\), \(\Lambda_\beta\)
- `inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` — Proposition 2, support functions, \(\kappa\) (M)
- `R/group_precision_floor.R` — production entry point (v1)
