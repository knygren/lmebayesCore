# Rosenthal TV-bounds for Restricted Two-Block Gibbs Samplers

Companion to:

- `inst/BETA_MARGINAL_MODE_LEVELSET.md` — marginal β mode and B̃(δ₂) construction
- `inst/GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` — Rosenthal drift and minorization
- `inst/JOINT_GAMMA_BETA_TV_CERTIFICATE.md` — joint restriction context

**Status:** implementation reference for the eight-layer exported API.
`certificate()` is unchanged (restricted γ-only Theorem 2 route).

**Formal statement.** The certified TV bound evaluated by
`rosenthal_tv_bound()` / `optimal_rosenthal_tv_bound()` is **Proposition R-Cert**
in `inst/GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` (one-line substitution
boxes only — general and sharpest displayed).

---

## Pipeline (categories 1–8)

```text
model_setup + pfamily_list + family + delta_2
        │
        ▼
(3) beta_marginal_mode()       ──► β† (Newton on γ-integrated marginal)
        │
        ▼
(4) beta_marginal_safe_set()   ──► B̃(δ₂) = {Ξ(β) ≤ r_Gauss(n, δ₂)}
        │
        ▼
(1) population_mode(beta_set)  ──► γ*, P₁₁, tilde J (restricted mean map;
        │                           ICM init + MC screen on non-Gaussian)
        ▼
(2) epsilon_star() / optimize  ──► ε* minorization at restricted γ*
        │
        ▼
(5) group_precision_floor()    ──► Γ_j^LB on certified level set
        │
        ▼
(6) floor_coupling_eigenvalues() ──► κ_i^LB from floor spectrum
        │
        ▼
(7) rosenthal_drift_constants()  ──► λ^LB, b, U
    rosenthal_tv_bound()         ──► bound at fixed α, k
        │
        ▼
(8) optimal_rosenthal_tv_bound() ──► optimize α and/or invert for k
        │
        ▼
    gamma_beta_tv_certificate()  ──► orchestrator (all steps + print method)
```

---

## Functions

| Cat | Function | File | Role |
|-----|----------|------|------|
| 1 | `population_mode()` | `R/c05_em.R` | EM / ICM mode γ* and block-1 precision |
| 1′ | `rNormal_reg_group_safe()` | `R/simfunction_group_safe.R` | Gaussian Block-1 + widetilde-B rejection |
| 1″ | `rNormalGLM_reg_group_safe()` | same | GLM Block-1 + widetilde-B rejection (primary) |
| 1‴ | `group_effects_conditional_mean()` | `R/group_effects_conditional_mean.R` | E[β_j \| γ, y] (aggregates safe draws) |
| 1⁴ | `beta_in_marginal_safe_set()` | same | Membership test Ξ(β) ≤ r_Gauss |
| 2 | `epsilon_star()` | `R/c05_epsilon_star.R` | Minorization ε* at γ* |
| 3 | `beta_marginal_mode()` | `R/beta_marginal_mode.R` | Marginal β† only (γ integrated out) |
| 4 | `beta_marginal_safe_set()` | `R/beta_marginal_safe_set.R` | Certified set B̃(δ₂); geometry only |
| 5 | `group_precision_floor()` | `R/group_precision_floor.R` | Per-group data precision lower bounds Γ_j^LB |
| 6 | `floor_coupling_eigenvalues()` | `R/floor_coupling_eigenvalues.R` | Coupling eigenvalues κ_i^LB |
| 6′ | `floor_coupling_spectrum()` | same | Back-compat wrapper (eigenvalues + drift) |
| 7 | `rosenthal_drift_constants()` | `R/rosenthal_tv_bound.R` | Foster drift constants (λ, b, U) |
| 7 | `rosenthal_tv_bound()` | `R/rosenthal_tv_bound.R` | Rosenthal TV bound at fixed `alpha`, `k` |
| 8 | `optimal_rosenthal_tv_bound()` | `R/optimal_rosenthal_tv_bound.R` | Optimize α; invert k for `inner_tol` / `total_tol` |
| — | `gamma_beta_tv_certificate()` | `R/gamma_beta_tv_certificate.R` | End-to-end orchestrator |

---

## Typical usage

### One-shot certificate

```r
cert <- gamma_beta_tv_certificate(
  design = setup,
  pfamily_list = plist,
  family = gaussian(),
  delta_2 = 0.01,
  k = 50L
)
print(cert)
```

Returns class `"gamma_beta_tv_certificate"` with components:
`mode`, `eps`, `beta_mode`, `beta_set`, `floor`, `eigenvalues`, `optimal`,
`delta_2`, `display_mode`, `certified`.

### Step-by-step (inspect intermediates)

```r
bmode <- beta_marginal_mode(design, pfamily_list, family = gaussian())
bset  <- beta_marginal_safe_set(beta_mode = bmode, delta_2 = 0.01)
mode  <- population_mode(design, pfamily_list, family = gaussian(),
                         beta_set = bset, estep = "exact")
eps   <- epsilon_star(mode, method = "optimize")
bstep <- group_effects_conditional_mean(mode = mode)
floor <- group_precision_floor(beta_mode = bmode, beta_set = bset,
                               pfamily_list = plist, family = gaussian())
ev    <- floor_coupling_eigenvalues(mode, floor)
opt   <- optimal_rosenthal_tv_bound(mode, ev, eps = eps$eps_star,
                                    inner_tol = 0.01)
```

For non-Gaussian families, pass `n` (draws per group) to
`group_effects_conditional_mean()` and `population_mode()`.

### β-restriction (widetilde B)

`rNormal_reg_group_safe()` (Gaussian) and **`rNormalGLM_reg_group_safe()`**
(non-Gaussian) draw via `rNormal_reg_group` / `rNormalGLM_reg_group`, then reject
when outside B̃(δ₂) using `beta_in_marginal_safe_set()`. C05 callers route by
`family` through `.c05_safe_group_draw()`.

```r
# one restricted Gibbs sweep (non-Gaussian)
b_draw <- rNormalGLM_reg_group_safe(
  n = 1L, y = design$y, x = design$D, group = design$group,
  prior_list = plist, family = binomial(), beta_set = bset
)$coefficients
```

### Full π_γ bound

When `include_full_pi_gamma = TRUE` (default in `gamma_beta_tv_certificate`
and `optimal_rosenthal_tv_bound`), the reported full bound is:

**inner Rosenthal bound + δ₂**

(triangle inequality against the β-truncation term; see
`GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` §0).

### Mode distance and `display_mode`

`display_mode = "sharp"` sets \(V(\gamma_0)=1\), which assumes \(\gamma_0=\gamma^\star\).
When the chain starts elsewhere, use `"general"` or inflate \(V(\gamma_0)\) via the
**matrix** strong-convexity bound in
`GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` (§3.1, “Mode distance when
\(\gamma_0\neq\gamma^\star\)”): \(\|d\|_{\Pi^{\mathrm{LB}}}^2\le g^\top(\Pi^{\mathrm{LB}})^{-1}g\)
with \(g=\nabla\Phi(\tilde\gamma)\) and \(\Pi^{\mathrm{LB}}=P_{11}-S^{\mathrm{LB}}\) from
`floor_coupling_eigenvalues()`.

With MC E-steps, the plug-in bound is random; see the same doc for
\(n^{-1/2}\) propagation via `mc_delta_floor`, the high-probability inflation,
certified \(V(\gamma_0)\le 1+\tfrac12 B_{\mathrm{cert}}(1+w_{\max}^{\mathrm{LB}})\),
and the **recommended stopping rule** (fixed-point layer (A) + mode certificate
(B)–(C), MC-linked `tol_eff` and \(\tau_B\), optional `v0_cert_tol`).

---

## Key arguments

| Argument | Used in | Meaning |
|----------|---------|---------|
| `delta_2` | safe set, optimal, certificate | Tail budget for B̃(δ₂) via r_Gauss(n, δ₂) |
| `kappa_method` | `group_precision_floor`, certificate | `"crude"` (certified), `"laplace"`, `"none"` |
| `inner_tol` | optimal, certificate | Target for β-restricted Rosenthal bound |
| `total_tol` | optimal, certificate | Target for full π_γ bound (uses total_tol − δ₂) |
| `alpha` | rosenthal, optimal | Rosenthal tuning; optimized if NULL |
| `display_mode` | drift, rosenthal, optimal | `"sharp"` (default) or `"general"` |

---

## Return object classes

| Class | Main fields |
|-------|-------------|
| `rNormal_reg_group_safe` / `rNormalGLM_reg_group_safe` | `coefficients`, `coefficients_all`, `accept_rate`, `restricted` |
| `group_effects_conditional_mean` | `b_mean`, `V_list`, `b_mc_se`, `beta_set`, `restricted` |
| `beta_marginal_mode` | `beta`, `hessian`, `engine`, `P_b` |
| `beta_marginal_safe_set` | `level`, `beta_mode`, `delta_2`, `P_b` |
| `group_precision_floor` | `Gamma_lb`, `per_group`, `level` |
| `floor_coupling_eigenvalues` | `kappa_lb`, `kappa_max_lb`, `S_lb`, `q` |
| `optimal_rosenthal_tv_bound` | `bound`, `k`, `alpha`, `rosenthal`, `drift`, `full_bound` |
| `gamma_beta_tv_certificate` | all of the above, plus `certified` metadata |

Each has a `print` method except `floor_coupling_eigenvalues` and
`rosenthal_drift_constants`.

---

## Validation

Scratch check (not run by `testthat`):

```r
source("data-raw/_chk_sharpest_tv_certificate.R")
```

Uses the Gaussian `bayesrules::big_word_club` fixture; verifies pipeline parity
and expected κ_max^LB ≈ 0.93 after the Gaussian Fisher / dispersion fix.

---

## Related docs

- `inst/BETA_MARGINAL_MODE_LEVELSET.md` §7 — phase-3 implementation table
- `inst/GAMMA_MARGINAL_DRIFT_MINORIZATION_ROSENTHAL.md` §8 — Rosenthal wiring
- `inst/JOINT_GAMMA_BETA_TV_CERTIFICATE.md` — joint restriction theory
