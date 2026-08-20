# Marginal-mode \(\beta\) level sets and \(r_{\mathrm{Gauss}}\) calibration

Companion to `inst/GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` (Rosenthal +
\(\delta_2\) outer step), `inst/JOINT_GAMMA_BETA_TV_CERTIFICATE.md` §4,
`inst/LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §4–§5,
`inst/LAPLACE_PROFILE_GROUP_MARGINAL_BOUND.md`, and
`inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` §3A.4.

**Status:** design note + implementation plan. Documents the **recommended**
construction of \(\widetilde B(\delta_2)\) for the step
\(\|\pi_{\gamma\mid\widetilde B}-\pi_\gamma\|_{TV}\le \pi_\beta(\widetilde B^c)\).

---

## 1. Problem

The full-\(\pi_\gamma\) certificate adds a **static \(\beta\)-truncation** term
(`GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` §0):

\[
\|P_\gamma^k-\pi_\gamma\|_{TV}
\;\le\;
\underbrace{\|P_\gamma^k-\pi_{\gamma\mid\widetilde B}\|_{TV}}_{\text{Rosenthal §3.1}}
\;+\;
\underbrace{\pi_\beta(\widetilde B(\delta_2)^{\,c})}_{\text{must be usable}}.
\]

Proposition (P2) (`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` §3A.4.1) gives a
**valid** upper bound on \(\pi_\beta(\widetilde B_r^c)\) but is **not tight** at
fixed state dimension \(n=Jp_{\mathrm{re}}\) (slack often \(10^2\)–\(10^3\times\);
see `data-raw/_chk_logit_levelset_bounds.R` Check 7). Inverting (P2) to hit
\(\delta_2\) forces an unnecessarily large level \(r\) and destroys weight floors.

**This note:** calibrate the level with a **Gaussian / Laplace reference** at the
**marginal \(\beta\) mode**, but define the set with the **true** convex
deficiency \(\Xi\).

---

## 2. Target marginal and anchor

After integrating population parameters \(\gamma\) out of the Gaussian prior
(`LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §4),

\[
\widetilde\pi(\beta\mid y)
\;\propto\;
\exp\Bigl(\sum_j \ell_j(\beta_j)\Bigr)\,
\exp\Bigl(-\tfrac12(\beta-\mu_\beta)^\top \Lambda_\beta (\beta-\mu_\beta)\Bigr),
\qquad
\beta\in\mathbb R^n,\ n=Jp_{\mathrm{re}}.
\]

**Marginal mode (required anchor).**

\[
\beta^\dagger
\;\in\;
\operatorname{argmax}_\beta \log \widetilde\pi(\beta\mid y),
\qquad
f(\beta):=-\log\widetilde\pi(\beta\mid y)+\text{const}.
\]

**True deficiency** (unique minimum \(0\) at \(\beta^\dagger\)):

\[
\boxed{
\Xi(\beta) := f(\beta)-f(\beta^\dagger)\ \ge\ 0.
}
\]

**Certified \(\beta\)-safe set.**

\[
\boxed{
\widetilde B(\delta_2)
\;:=\;
\bigl\{\beta:\ \Xi(\beta)\ \le\ r_{\mathrm{Gauss}}(n,\delta_2)\bigr\}.
}
\]

**Do not** center \(\Xi\) at \(b(\gamma^\star)=E[\beta\mid\gamma^\star,y]\) or at a
joint ICM mode unless it coincides with \(\beta^\dagger\): then \(\Xi(\beta^\dagger)<0\)
and level sets are not HPD superlevel sets of \(\widetilde\pi(\beta\mid y)\).

| Anchor | Role |
|--------|------|
| \(\gamma^\star\) | C05 mean-map fixed point (`population_mode()`) — centre of \(Q\), \(\varepsilon(\gamma^\star)\) |
| \(\beta^\dagger\) | **This note** — origin of \(\Xi\), Laplace Hessian, \(\widetilde B(\delta_2)\) |
| \(\beta^{\mathrm{ICM}}_\star\) | Optional joint working point for \(P_{22}\) at \((\gamma^\star,\beta)\); \(\to\beta^\dagger\) as data grow |

`CHAPTER_C05_IMPLEMENTATION.md` §2.3: never substitute conditional **modes** for
**means** in the \(\gamma\) E-step. That warning is separate from \(\beta^\dagger\).

---

## 3. Level calibration: \(r_{\mathrm{Gauss}}\)

Laplace at \(\beta^\dagger\): \(H:=\nabla^2 f(\beta^\dagger)\),

\[
\widetilde\pi(\beta\mid y)
\approx
c\,\exp\Bigl(-\tfrac12(\beta-\beta^\dagger)^\top H(\beta-\beta^\dagger)\Bigr).
\]

Under \(N(\beta^\dagger,H^{-1})\), with \(\Xi_{\mathrm{Lap}}(\beta)=\tfrac12(\beta-\beta^\dagger)^\top H(\beta-\beta^\dagger)\),

\[
r_{\mathrm{Gauss}}(n,\delta_2)
\;:=\;
\tfrac12\,\chi^2_{n,\,1-\delta_2}
\;=\;
q_{1-\delta_2}\bigl(\mathrm{Gamma}(n/2,1)\bigr).
\]

Then **under the Gaussian reference**,

\[
P_{\mathrm{Lap}}\bigl(\Xi_{\mathrm{Lap}} > r_{\mathrm{Gauss}}\bigr)=\delta_2
\quad\text{exactly}.
\]

Same convention as `mvn_calibrate()` / `.c05_mvn_d()` on the \(\gamma\) side
(`R/c05_mvn.R`). **Not** the same as inverting Proposition (P2):
\(r_{\mathrm{P2}}(n,\delta)\approx r_{\mathrm{Gauss}}(2n,\delta)\) at fixed \(\delta\)
(`LAPLACE_PROFILE_GROUP_MARGINAL_BOUND.md` §5; `GAUSSIAN_MAJORIZATION` §1658).

**Radius reporting.** \(s=\sqrt{2r}\) is the deficiency radius in the same units
as `group_precision_floor()$level$s` (`GAUSSIAN_MAJORIZATION` §1682 ff.).

---

## 4. Asymptotics (why this is the right large-sample route)

Fix \(J,p_{\mathrm{re}}\), priors. Let trial counts \(n_j\to\infty\) (equivalently
total \(N\to\infty\)).

**Bernstein–von Mises / Laplace.** Under standard GLMM regularity,

\[
\widetilde\pi(\beta\mid y)
\;\Rightarrow\;
N\bigl(\beta^\dagger,\,H_N^{-1}\bigr)
\quad\text{in total variation (locally)}.
\]

**Proposition (asymptotic mass calibration — template).** With
\(r_N=r_{\mathrm{Gauss}}(n,\delta_2)\) and \(\widetilde B_N=\{\Xi_N\le r_N\}\),

\[
\pi_{\beta,N}\bigl(\widetilde B_N(\delta_2)^{\,c}\bigr)
\;=\;
\delta_2 + O(N^{-1/2})
\]

(typically \(O(N^{-1})\) with uniform third-derivative bounds on a
\(O(\sqrt{\log N})\) Laplace neighbourhood).

**Interpretation.**

- **Finite \(n\):** treat \(\delta_2\) as the **design target** under the Laplace
  law; validate true \(\pi_\beta(\widetilde B^c)\) by integration when a proof is
  needed (`LOGIT_STATIC_TAIL_CERTIFICATION.md` §4.5 pattern).
- **Large \(N\):** \(r_{\mathrm{Gauss}}\) on **true** \(\Xi\) level sets is
  asymptotically **tight** for the tail mass; (P2) slack does **not** vanish with
  \(N\) at fixed \(n\).

**Anchors converge.** As \(N\to\infty\),
\(\beta^\dagger\), \(E[\beta\mid y]\), \(b(\gamma^\star)\), and
\(\beta^{\mathrm{ICM}}_\star\) agree at rate \(O(N^{-1/2})\); using \(\beta^\dagger\)
minimizes the Laplace remainder in the \(\Xi\) coordinate.

**Weak data (\(n_j\) fixed, small).** Asymptotics do not apply; profile/marginal
geometry may still be handled per group (`LAPLACE_PROFILE`); tail mass may differ
from \(\delta_2\) by \(O(1)\) — integration or conservative inflation required.

---

## 5. Floors on \(\widetilde B(\delta_2)\)

Mass calibration and **weight floors** use the **same** set \(\{\Xi\le r\}\).

On \(\widetilde B(\delta_2)\), for each \((j,i)\),

\[
\omega_{j,i}(\delta_2)
\;=\;
\inf_{\beta\in\widetilde B(\delta_2)} w_{j,i}(\eta_{j,i}(\beta)),
\qquad
\eta_{j,i}=z_{ij}^\top\beta_j + \text{offset}.
\]

For logit, \(w(\eta)=np(1-p)\) is quasi-concave in \(\eta\); the minimum is at
**support-function endpoints** on \(\partial\widetilde B\) (KKT / `supp()` —
`GAUSSIAN_MAJORIZATION` §3A.4.2; `data-raw/_ex_logit_floor_J10_p3.R`).

**Joint level (this note).** One \(r=r_{\mathrm{Gauss}}(n,\delta_2)\) for all groups;
support solves on \(\{\Xi\le r\}\) with other blocks at the **profile** point
(not frozen at the joint mode).

**Per-group union (optional).** Budget \(\delta_2/J\), dimension \(p_{\mathrm{re}}\),
profile sets \(\{\Xi^{\mathrm{prof}}_j\le r_j\}\) — see `LAPLACE_PROFILE` and
`group_precision_floor()`. Differs from the **joint** \(\widetilde B\) above but
shares \(r_{\mathrm{Gauss}}\) machinery.

Assemble \(P_{22,j}^{\mathrm{LB}}=P_b+\underline{\mathcal P}_{j,\mathrm{data}}\),
then \(\kappa_i^{\mathrm{LB}}(\delta_2)\) and Rosenthal drift (`GAMMA` §2).

---

## 6. What we do **not** use for \(\delta_2\)

| Route | Issue |
|-------|--------|
| Invert (P2) \(\Gamma(n,r)/\gamma(n,r)=\delta_2\) | Valid but loose; \(r\approx 2\times r_{\mathrm{Gauss}}\); true mass \(\ll\delta_2\) |
| Laplace ellipsoid \(\{\Xi_{\mathrm{Lap}}\le r\}\) only | Wrong shape for logit floors |
| MC quantile of \(\Xi\) without certificate | Good estimate, not a proof |
| Center \(\Xi\) at \(b(\gamma^\star)\) | Breaks \(\Xi\ge0\) with min at origin |

---

## 7. Implementation path

### Phase 0 — Internal numerics (done / in progress)

| Item | Location | Notes |
|------|----------|--------|
| \(r_{\mathrm{Gauss}}(n,\delta_2)\) | `R/c05_beta_marginal_set.R` `.c05_beta_r_gauss_level()` | \(r=\tfrac12 q_{1-\delta_2}(\chi^2_n)\); mirrors `.c05_mvn_d()` |
| Laplace tail reference | `.c05_beta_laplace_tail_mass()` | \(P(\Xi_{\mathrm{Lap}}>r)\) under \(N(0,I_n)\) reference |
| Marginal engine | `R/group_precision_floor_engine.R` `.group_floor_engine()` | \(f\), \(\Xi\), `He_f`, Newton — **already marginal** target |
| Prop (P2) level (legacy) | `.gamma_level_for_budget()` | Keep for `group_precision_floor()` until migrated |

### Phase 1 — Marginal mode (replace ICM start)

| Item | Action |
|------|--------|
| Mode finder | Newton from \(\mu_\beta\): `engine$newton(mu_b)` on \(f(\beta)\) |
| Validation | Compare \(\|\beta^\dagger_{\mathrm{Newton}}-\beta^{\mathrm{ICM}}_\star\|\); warn if large (small \(n_j\)) |
| `group_precision_floor()` | Add `mode_method = c("marginal_newton", "icm")` (default **`marginal_newton`**) |
| Return object | `mode$beta_dagger`, `mode$hessian`, `mode$method` |

**No new public arguments on `population_mode()` or `certificate()`** without a
separate API review.

### Phase 2 — Level scheme switch

| Item | Action |
|------|--------|
| `level_method` | `"r_gauss_joint"` (default target), `"prop2_groups"`, `"prop2_joint"` (legacy) |
| Joint \(\widetilde B\) | Single \(R_j\leftarrow r_{\mathrm{Gauss}}(n,\delta_2)\) for **all** support solves on \(\{\Xi\le R\}\) |
| Per-group legacy | Current `group_precision_floor()` behaviour under `"prop2_groups"` |
| Report | `level$r_gauss`, `level$laplace_tail_mass`, `level$scheme` |

### Phase 3 — Drift / Rosenthal wiring (**implemented**)

| Item | Location | Notes |
|------|----------|--------|
| Marginal safe set | **`beta_marginal_safe_set()`** | `R/beta_marginal_safe_set.R`; Newton mode + joint `r_gauss_joint` |
| Floor spectrum | **`floor_coupling_spectrum()`** | `R/c05_floor_spectrum.R`; \(\kappa_i^{\mathrm{LB}}(\delta_2)\) from \(\Gamma_j^{\mathrm{LB}}\) |
| Rosenthal bound | **`rosenthal_tv_bound()`** | `R/c05_rosenthal_tv.R`; sharp vs general `display_mode` |
| Full certificate | **`gamma_beta_tv_certificate()`** | `R/gamma_beta_tv_certificate.R`; inner + \(\delta_2\) full-\(\pi_\gamma\) row |
| Legacy floors | `group_precision_floor()` | `mode_method`, `level_method` for parity; default unchanged (`icm` + Prop 2) |
| Scratch check | `data-raw/_chk_sharpest_tv_certificate.R` | Gaussian smoke + spectrum parity |

`certificate()` is unchanged (restricted \(\gamma\)-only Theorem 2 route).

### Phase 4 — Tail mass (asymptotic + finite-\(n\) validation)

| Item | Action |
|------|--------|
| **Asymptotic report** | `laplace_tail_mass = delta_2` by construction at \(r_{\mathrm{Gauss}}\) |
| **True mass (optional)** | `data-raw/_chk_beta_marginal_tail_mass.R`: grid \(\gamma\mid y\) + conditional \(\beta_j\mid\gamma\) MC; bracket if needed |
| **Full TV row** | Document \(\|P^k-\pi_\gamma\|_{TV}\le\) Rosenthal \(+\,\widehat\pi_\beta(\widetilde B^c)\) with explicit asymptotic disclaimer |

### Dependency graph

```text
model_setup + pfamily_list
        │
        ▼
.group_floor_integrate_gamma()  ──► Lambda_beta, mu_beta
        │
        ▼
.group_floor_engine() + Newton    ──► beta_dagger, H = He_f(beta_dagger)
        │
        ▼
.c05_beta_r_gauss_level(n, delta_2) ──► r
        │
        ▼
support on {Xi <= r}              ──► omega_ji, Gamma_j^LB
        │
        ├─► kappa_i^LB(delta_2) ──► GAMMA Rosenthal (lambda, b)
        │
        └─► delta_2 asymptotic / MC tail ──► full pi_gamma TV
```

### Tests and checks (not `tests/testthat/` until approved)

- Extend `data-raw/_ex_logit_floor_J10_p3.R`: floors at \(r_{\mathrm{Gauss}}\) vs (P2).
- New `data-raw/_chk_beta_marginal_tail_mass.R`: true vs Laplace mass at fixed \(r\).
- Reuse `_chk_group_marginal_bound.R` coupling diagnostics.

---

## 8. References

- `inst/GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` — Rosenthal + \(\delta_2\)
- `inst/JOINT_GAMMA_BETA_TV_CERTIFICATE.md` §4 — certified \(\beta\)-set definition
- `inst/CHAPTER_C05_IMPLEMENTATION.md` §2.3 — mean vs mode trap (\(\gamma\))
- `R/group_precision_floor.R` — current floors (ICM + Prop 2); migration target
- `R/c05_mvn.R` — \(\chi^2\) calibration precedent on \(\gamma\)
