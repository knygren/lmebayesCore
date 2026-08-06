# The Joint Hessian for (Fixed Effects, Random Effects, Precision) and the State-Dependent Extension of Claims 1/3

This document is the concrete derivation promised as "future work" at the end of
`inst/BLOCK_GIBBS_ERGODICITY.md`. It extends the multivariate-normal two-block
theory (Nygren 2020, as implemented in `R/two_block_ergodicity.R` and used by
`two_block_rate()`) to a model where the random-effects precision
\(\lambda_p := 1/\tau^2_p\) is **also** sampled, under the
`dIndependent_Normal_Gamma()` ("ING") Block~2 prior.

**Scope.** This is a theory note only — no code changes are implied or made
here. It derives the Hessian of the negative log-posterior over
\((\gamma, \beta, \lambda)\) jointly, shows exactly which blocks of the
existing \(P_{11}\)/\(P_{12}\)/\(P_{22}\) construction stay the same and which
are new, and explains precisely why the resulting analogue of the Claim 1 /
Claim 3 contraction matrix \(A\) is **state-dependent** — i.e. why it is not a
constant, unlike the fixed-\(\Psi\) case `two_block_rate()` currently handles.

---

## 1. Notation

This note follows `inst/notation.md`'s "centered" Lindley & Smith (1972)
parameterization (\(\beta_j\), \(\gamma\), \(\mathcal{W}_j\), \(\Psi\)), since
that is the notation used in `rLMM_reg.R`'s current documentation. A
translation table to the code symbols used in `R/two_block_ergodicity.R`
(\(Z\), \(b\), \(H_j\), \(P_b\)) is given at the end of §2.

| Symbol | Role | Dimension |
|---|---|---|
| \(y_j\) | response, group \(j\) | \(n_j \times 1\) |
| \(D_j\) | likelihood design, group \(j\) | \(n_j \times P\) |
| \(\Omega_j\) | per-observation likelihood precision (diag; \(1/\sigma^2\) for Gaussian, Fisher weights otherwise) | \(n_j \times n_j\) |
| \(\beta_j\) | group-level coefficients (Block 1) | \(P \times 1\) |
| \(\beta_{jp}\) | \(p\)-th entry of \(\beta_j\) | scalar |
| \(\gamma\) | population / hyper coefficients (Block 2), stacked over \(p = 1,\dots,P\) | \(q \times 1,\ q = \sum_p q_p\) |
| \(\gamma_p\) | sub-vector of \(\gamma\) for RE dimension \(p\) | \(q_p \times 1\) |
| \(\mathcal{W}_j\) | level-2 design, group \(j\); block-diag\((W_{1j},\dots,W_{Pj})\) | \(P \times q\) |
| \(W_p\) | stacked level-2 design for dimension \(p\) (\(J\times q_p\); row \(j\) is \(W_{pj}\)) | \(J \times q_p\) |
| \(\Psi\) | covariance of \(\beta_j\) about \(\mathcal{W}_j\gamma\) | \(P \times P\) |
| \(\lambda_p := (\Psi^{-1})_{pp}\) | precision of RE dimension \(p\) | scalar |
| \(u_{jp} := \beta_{jp} - W_{pj}\gamma_p\) | Block~1/Block~2 residual | scalar |
| \(\mu_p^0, V_p\) | Normal prior mean/covariance for \(\gamma_p\) | \(q_p\times 1,\ q_p\times q_p\) |
| \(a_p^0, r_p^0\) | Gamma **shape/rate** prior for \(\lambda_p\) (`dIndependent_Normal_Gamma(shape=, rate=)`) | scalar |
| \(\mathrm{ING} \subseteq \{1,\dots,P\}\) | RE dimensions with a sampled (ING) precision; the complement has fixed, known \(\lambda_p\) (`dNormal()`) | — |

Package convention: **\(\Psi^{-1} = \mathrm{diag}(\lambda_1,\dots,\lambda_P)\)
is always diagonal** — there is no cross-dimension RE covariance in the
current model family. This is what makes the block structure below tractable
(no dense \(P\times P\) precision to track — only the \(P\) scalars
\(\lambda_p\)).

---

## 2. The model, and why \(\gamma\) and \(\lambda\) form one Gibbs block

**Stage 1 (likelihood):** \(y_j \mid \beta_j \sim N(D_j\beta_j, \Omega_j^{-1})\)
(or a GLM likelihood, with \(\Omega_j\) the working/Fisher precision at the
current \(\beta_j\); this note keeps the Gaussian case exact and notes the GLM
extension only informally, exactly as `two_block_rate()` already does via
`two_block_mode_weights()`).

**Stage 2 (RE prior / hyper-regression):**
\(\beta_{jp} \mid \gamma_p, \lambda_p \sim N(W_{pj}\gamma_p,\ 1/\lambda_p)\),
independent across \(p\) (diagonal \(\Psi^{-1}\)) and across \(j\) (independent
groups).

**Stage 3 (Block~2 priors):**

- \(\gamma_p \sim N(\mu_p^0, V_p)\) — **independent of \(\lambda_p\)** (this is
  the "Independent" in `dIndependent_Normal_Gamma`, as opposed to a classical
  conjugate Normal-Gamma prior in which \(\gamma_p\)'s prior precision would
  itself be scaled by \(\lambda_p\)).
- For \(p \in \mathrm{ING}\): \(\lambda_p \sim \mathrm{Gamma}(a_p^0, r_p^0)\)
  (shape/rate). For \(p \notin \mathrm{ING}\), \(\lambda_p\) is a fixed known
  constant (`dNormal()`), contributing no state dimension.

The sampler updates **Block 1** (\(\beta\)) first, then **Block 2**
(\((\gamma_p,\lambda_p)\) jointly, per ING component, via the exact Normal-Gamma
conjugate full conditional) second. This is exactly the point of this note:
**Block 2's state is \((\gamma,\lambda)\), not just \(\gamma\)**, once any
component is ING. The two-block Gibbs cycle is therefore already
\(x_1 = (\gamma,\lambda_{\mathrm{ING}})\), \(x_2 = \beta\) — the same
"\(x_1\) updated second / \(x_2\) updated first" convention
`R/two_block_ergodicity.R` uses, just with \(x_1\) enlarged.

**Translation to `R/two_block_ergodicity.R` code symbols:**

| This note | Code (`two_block_ergodicity.R`) |
|---|---|
| \(D_j\) (likelihood design) | `x[rows, ]` (a.k.a. \(Z_j\)) |
| \(\beta_j\) | `b[j, ]` |
| \(\mathcal{W}_j\), \(W_{pj}\) | `H_j`, `x_j[[k]]` |
| \(\gamma_p\) | `gamma[cols[[k]]]` |
| \(\Psi^{-1} = \mathrm{diag}(\lambda_p)\) | `P_b` |
| \(\Omega_j\) | `w[rows]` (diag) |
| \(V_p\) | `prior_list_block2[[k]]$Sigma` inverted |

---

## 3. The negative log-posterior

Dropping additive constants that do not involve \(\gamma\), \(\beta\), or
\(\lambda\):

$$
\ell(\gamma,\beta,\lambda) \;=\;
\underbrace{\tfrac12\sum_j (y_j - D_j\beta_j)'\Omega_j(y_j - D_j\beta_j)}_{\text{Block 1 likelihood}}
\;+\;
\underbrace{\tfrac12\sum_j\sum_p \lambda_p\, u_{jp}^2 \;-\; \tfrac{J}{2}\sum_{p \in \mathrm{ING}} \log\lambda_p}_{\text{RE prior } \beta_{jp}\mid\gamma_p,\lambda_p}
\;+\;
\underbrace{\tfrac12\sum_p (\gamma_p - \mu_p^0)'V_p^{-1}(\gamma_p - \mu_p^0)}_{\text{Block 2 prior on }\gamma}
\;+\;
\underbrace{\sum_{p \in \mathrm{ING}} \Big[(1-a_p^0)\log\lambda_p + r_p^0\lambda_p\Big]}_{\text{Block 2 prior on }\lambda}
$$

The \(-\tfrac{J}{2}\log\lambda_p\) term is the normalizing constant of the
Gaussian RE-prior density (it depends on \(\lambda_p\) only for
\(p \in \mathrm{ING}\); for fixed-\(\lambda_p\) dimensions it is a constant and
drops out of the Hessian entirely, so those dimensions reduce exactly to the
existing fixed-\(\Psi\) theory).

This is the **exact** joint negative log-density (up to an additive
constant); no Laplace approximation has been used yet. The only place a
local/quadratic approximation enters is in treating the resulting Hessian
as a *rate matrix* for a nonlinear one-step Gibbs map (§7) — the model itself
is unchanged.

---

## 4. First derivatives (needed for two things: the score equations that
   define the joint mode, and to compute the *second* derivatives)

$$
\frac{\partial \ell}{\partial \beta_{jp}}
= \big[D_j'\Omega_j(D_j\beta_j - y_j)\big]_p \;+\; \lambda_p u_{jp}
$$

$$
\frac{\partial \ell}{\partial \gamma_p}
= -\lambda_p \sum_j u_{jp} W_{pj}' \;+\; V_p^{-1}(\gamma_p - \mu_p^0)
= -\lambda_p W_p' u_p + V_p^{-1}(\gamma_p-\mu_p^0), \quad u_p := (u_{1p},\dots,u_{Jp})'
$$

$$
\frac{\partial \ell}{\partial \lambda_p}
= \tfrac12\sum_j u_{jp}^2 - \tfrac{J}{2\lambda_p} + \frac{1-a_p^0}{\lambda_p} + r_p^0
\qquad (p \in \mathrm{ING})
$$

Setting the last equation's negative to zero and solving reproduces exactly
the Block~2 Gibbs step's exact conjugate full conditional
\(\lambda_p \mid \beta,\gamma \sim \mathrm{Gamma}\big(a_p^0 + J/2,\ r_p^0 + \tfrac12\sum_j u_{jp}^2\big)\)
— i.e. the same `shape2 = shape + n_w/2.0` update implemented in
`src/rIndepNormalGammaReg.cpp` / `src/block_rIndepNormalGammaReg.cpp`. This is
a useful correctness check before trusting the second derivatives below.

---

## 5. Second derivatives (the Hessian blocks)

Write \(r_{jk\text{-type}}\) terms out fully; every cross term not listed
below is **exactly zero** (see the "why zero" column) — this is what keeps the
extended Hessian nearly as sparse as the original.

| Block | Formula | New vs. original? | Why (non)zero |
|---|---|---|---|
| \(\partial^2\ell/\partial\beta_{jp}^2\) (same \(j,p\)) | \([D_j'\Omega_jD_j]_{pp} + \lambda_p\) | same *formula*; \(\lambda_p\) now a variable | direct |
| \(\partial^2\ell/\partial\beta_{jp}\partial\beta_{jp'}\), \(p\ne p'\) | \([D_j'\Omega_jD_j]_{pp'}\) | unchanged (likelihood-only) | RE-prior term is diagonal in \(p\) |
| \(\partial^2\ell/\partial\beta_{jp}\partial\beta_{j'p}\), \(j\ne j'\) | \(0\) | unchanged | groups conditionally independent |
| \(\partial^2\ell/\partial\beta_{jp}\partial\gamma_p'\) | \(-\lambda_p W_{pj}'\) | same *formula*, \(\lambda_p\) now a variable | direct |
| \(\partial^2\ell/\partial\beta_{jp}\partial\gamma_{p'}'\), \(p\ne p'\) | \(0\) | unchanged | \(u_{jp}\) does not involve \(\gamma_{p'}\) |
| \(\partial^2\ell/\partial\beta_{jp}\partial\lambda_p\) | \(u_{jp}\) | **new** | \(\partial(\lambda_p u_{jp})/\partial\lambda_p = u_{jp}\) |
| \(\partial^2\ell/\partial\beta_{jp}\partial\lambda_{p'}\), \(p\ne p'\) | \(0\) | **new (trivially zero)** | \(\beta_{jp}\) only enters the \(p\)-th RE-prior term |
| \(\partial^2\ell/\partial\gamma_p\partial\gamma_p'^\top\) | \(\lambda_p W_p'W_p + V_p^{-1}\) | same *formula*, \(\lambda_p\) now a variable | direct |
| \(\partial^2\ell/\partial\gamma_p\partial\gamma_{p'}'^\top\), \(p\ne p'\) | \(0\) | unchanged | diagonal \(\Psi^{-1}\) + independent priors |
| \(\partial^2\ell/\partial\gamma_p\partial\lambda_p\) | \(-W_p'u_p\) | **new** | \(\partial(-\lambda_pW_p'u_p)/\partial\lambda_p\), holding \(u_p\) fixed |
| \(\partial^2\ell/\partial\gamma_p\partial\lambda_{p'}\), \(p\ne p'\) | \(0\) | **new (trivially zero)** | independence across \(p\) |
| \(\partial^2\ell/\partial\lambda_p^2\) | \(\dfrac{J/2 + a_p^0 - 1}{\lambda_p^2}\) | **new** | see derivation and sanity check below |
| \(\partial^2\ell/\partial\lambda_p\partial\lambda_{p'}\), \(p\ne p'\) | \(0\) | **new (trivially zero)** | independent Gamma priors, no shared term |

**Deriving \(\partial^2\ell/\partial\lambda_p^2\):** the RE-prior quadratic term
\(\tfrac12\lambda_p\sum_j u_{jp}^2\) is *linear* in \(\lambda_p\), so it
contributes nothing to \(\partial^2/\partial\lambda_p^2\) (only to the
\(\beta\)/\(\gamma\)-\(\lambda\) cross terms above). What remains is purely the
two \(\log\lambda_p\) terms: \(-\tfrac{J}{2}\log\lambda_p\) contributes
\((J/2)/\lambda_p^2\), and \((1-a_p^0)\log\lambda_p\) contributes
\((a_p^0-1)/\lambda_p^2\). Summing gives \((J/2 + a_p^0 - 1)/\lambda_p^2\).

**Sanity check.** This must equal minus the second derivative of the log of
the exact Gamma **full conditional** density
\(\mathrm{Gamma}(a_p^0+J/2,\ \cdot)\) with respect to its own argument, since
partial derivatives don't care what else is held fixed. For
\(X \sim \mathrm{Gamma}(a,\text{rate})\), \(-\partial^2\log f(x)/\partial x^2 = (a-1)/x^2\).
With \(a = a_p^0 + J/2\) this is exactly \((a_p^0 + J/2 - 1)/\lambda_p^2\) —
matches. \(\blacksquare\)

---

## 6. Assembling the extended \(P_{11}\), \(P_{12}\), \(P_{22}\)

Let \(x_1 = (\gamma,\ \lambda_{\mathrm{ING}})\), dimension
\(q + |\mathrm{ING}|\), and \(x_2 = \beta\), dimension \(JP\) (matching
`two_block_ergodicity.R`'s \(x_1\)/\(x_2\) convention, just with \(x_1\)
enlarged by the ING precisions).

**Extended \(P_{22}\)** (\(JP \times JP\), block-diagonal in \(j\), unchanged
in *form*):

$$
P_{22} = \mathrm{blockdiag}_j\big(B_j\big), \qquad B_j = D_j'\Omega_jD_j + \Psi^{-1} = D_j'\Omega_jD_j + \mathrm{diag}(\lambda_1,\dots,\lambda_P)
$$

Identical to the original `B_j <- crossprod(Z_j, Z_j * w_j) + P_b`. The
*only* difference is that \(\mathrm{diag}(\lambda)\) is now the **current
state**, not a fixed constant plugged in once.

**Extended \(P_{12}\)** (size \((q+|\mathrm{ING}|) \times JP\)): for each
group \(j\),

$$
P_{12}[\,\cdot,\ \beta_j\,] =
\begin{pmatrix} -\mathcal{W}_j'\,\Psi^{-1} \\[4pt] u_{j,\mathrm{ING}}' \end{pmatrix}
$$

where the top \(q \times P\) block is exactly the original
\(-H_j'P_b\), and the **new** bottom \(|\mathrm{ING}| \times P\) block has row
\(p \in \mathrm{ING}\) equal to \(u_{jp}\) in column \(p\) and \(0\) elsewhere
(from \(\partial^2\ell/\partial\beta_{jp}\partial\lambda_p = u_{jp}\),
\(\partial^2\ell/\partial\beta_{jp'}\partial\lambda_p = 0\) for \(p'\ne p\)).

**Extended \(P_{11}\)** (size \((q+|\mathrm{ING}|)\times(q+|\mathrm{ING}|)\)),
in \((\gamma,\lambda)\) partitioned form:

```
        ┌ P11^(γγ)     P11^(γλ) ┐
P11  =  │                        │
        └ P11^(γλ)'    P11^(λλ) ┘
```

- \(P_{11}^{(\gamma\gamma)}\): same formula as the original \(P_{11}\)
  (block-diagonal in \(p\): \(\lambda_p W_p'W_p + V_p^{-1}\)), \(\lambda\) now
  state-dependent.
- \(P_{11}^{(\lambda\lambda)}\): **new**, diagonal,
  \(\mathrm{diag}_{p\in\mathrm{ING}}\big((J/2+a_p^0-1)/\lambda_p^2\big)\).
- \(P_{11}^{(\gamma\lambda)}\): **new**, block-sparse — column \(p\in\mathrm{ING}\)
  has \(-W_p'u_p\) in the \(\gamma_p\) rows and \(0\) in every \(\gamma_{p'}\)
  row for \(p'\ne p\).

Every zero pattern above is a *structural* consequence of (i) diagonal
\(\Psi^{-1}\), (ii) independent Block~2 priors across \(p\), and (iii)
independent groups — exactly the same three assumptions the original theory
already relies on for its own block-diagonal-in-\(p\) structure. Nothing
extra had to be assumed to get this sparsity.

**Comparison table** (what changes vs. `.two_block_S_P11()`'s current
assembly):

| Quantity | Fixed-\(\Psi\) theory (current) | ING extension (this note) |
|---|---|---|
| \(x_1\) dimension | \(q\) | \(q + \lvert\mathrm{ING}\rvert\) |
| \(P_{22}\) | \(B_j = D_j'\Omega_jD_j + \Psi^{-1}\), \(\Psi\) fixed | same formula, \(\Psi^{-1}\) = current state |
| \(P_{12}\) | \(-\mathcal{W}_j'\Psi^{-1}\) only | \(+\) new \(u_{j,\mathrm{ING}}'\) row block |
| \(P_{11}\) | \(\sum_j \mathcal{W}_j'\Psi^{-1}\mathcal{W}_j + \mathrm{blockdiag}(V_p^{-1})\) | \(+\) new \(\lambda\lambda\) diagonal block, \(+\) new \(\gamma\lambda\) cross block |
| Constant in state? | **Yes** (given fixed \(\Psi\)) | **No** — every block above depends on \((\gamma,\beta,\lambda)\) through \(\lambda\) itself and/or the residuals \(u_{jp}\) |

---

## 7. The extended rate matrix \(A\), and why Claims 1/3 no longer apply as stated

Formally define, exactly as `.two_block_gen_eigen()` does,

$$
A(\gamma,\beta,\lambda) \;=\; P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2}
$$

using the *extended* blocks of §6, evaluated at the current state
\((\gamma,\beta,\lambda)\). The original Claims 1/3 (`two_block_rate()`'s
\(\lambda^*\), and the geometric TV bounds of `two_block_tv_bound()`) rely
critically on \(A\) being a **single, state-independent matrix**, because that
is exactly what makes the two-block Gibbs kernel an *exactly linear* Gaussian
map: each sweep contracts the (centered) state by the same matrix \(A\)
regardless of where the chain currently is, which is what licenses writing
the \(l\)-step contraction as \(A^l\)/\(A^{2l}\) and citing Nygren's exact
finite-sample TV bounds.

Once \(\lambda\) is part of the sampled state, \(A(\gamma,\beta,\lambda)\) is a
genuine **function of the current draw** — through \(\lambda\) itself (in
\(P_{22}\) and \(P_{11}^{(\gamma\gamma)}\)) and through the residuals
\(u_{jp}\) (in the new \(P_{12}\) and \(P_{11}^{(\gamma\lambda)}\) blocks).
The one-step Gibbs map is therefore **nonlinear** in \((\gamma,\beta,\lambda)\)
jointly (consistent with the earlier observation, in this same conversation,
that the mean-map for \(\lambda\) is nonlinear and its fixed point need not
coincide with any single plug-in "starting value"). \(A(\cdot)\) above is
precisely the **local linearization (Jacobian) of that nonlinear map** around
a given state — the same role a Jacobian plays in classical local stability
analysis of a nonlinear dynamical system, rather than the *global* transition
operator Claims 1/3 need.

Concretely, this means:

- **A single global \(\lambda^*\)** (in the sense of `two_block_rate()`) does
  not exist for the joint \((\gamma,\beta,\lambda)\) chain in general; at best
  one can define a **local** rate \(A(\gamma,\beta,\lambda)\) at a reference
  point (e.g. the joint mode, or the current sweep's draw) and interpret its
  spectral radius as an *instantaneous* contraction diagnostic, not a
  certified global TV bound.
- Nygren's Theorem 3 / Corollary 1 (`two_block_tv_bound()`) are proved for an
  *exactly* Gaussian target with a state-independent kernel; they do not
  apply to \(A(\cdot)\) as defined above without re-derivation for a
  state-dependent (non-Gaussian-kernel) Gibbs chain. No such re-derivation is
  attempted here.
- `plot_var_convergence()`/`plot_mean_convergence()`'s existing
  `component = "precision"` mode (which only tracks each \(1/\tau_p^2\)'s own
  marginal cross-chain variance against its empirical last-sweep value) is
  therefore, in this light, best understood as tracking one **diagonal
  entry** of \(P_{11}^{(\lambda\lambda)-1}\) computed *locally* — a
  reasonable and honest diagnostic, but not (yet) backed by a global
  contraction theorem the way `two_block_rate()`'s \(\lambda^*\) is for the
  fixed-\(\Psi\) case.

---

## 8. Special cases and sanity checks

**Reduction to the original theory.** If \(\mathrm{ING} = \varnothing\) (all
RE dimensions use `dNormal()`, i.e. fixed \(\Psi\)), every "new" row/column
above vanishes identically — \(x_1\) collapses back to \(\gamma\) alone, and
\(P_{11}\), \(P_{12}\), \(P_{22}\) reduce exactly to
`.two_block_S_P11()`'s current formulas, with \(A(\gamma,\beta,\lambda)\)
collapsing to the fixed, state-independent \(A\) of `two_block_rate()`. The
extension is a strict generalization, not a different theory.

**When does the \(\gamma_p\)-\(\lambda_p\) cross term vanish?** At a
stationary point of \(\ell\) with respect to \(\gamma_p\)
(\(\partial\ell/\partial\gamma_p = 0\)), \(\lambda_p W_p'u_p = V_p^{-1}(\gamma_p-\mu_p^0)\).
If the Block~2 prior on \(\gamma_p\) is flat/improper (\(V_p^{-1}=0\)) — or
negligible relative to the data — this forces \(W_p'u_p = 0\) exactly, so the
new \(P_{11}^{(\gamma\lambda)}\) cross-block vanishes *at the mode* even
though it is nonzero away from it. With a genuinely informative
\(\gamma_p\)-prior, the cross term persists at the mode too.

**When does the \(\lambda_p\)-\(\beta\) cross term vanish?** Row \(j\) of the
new \(P_{12}\) block for component \(p\) is \(u_{jp}\), the *residual for
group \(j\)*. At the joint mode, \(\partial\ell/\partial\beta_{jp}=0\) gives
\(\lambda_p u_{jp} = -[D_j'\Omega_j(D_j\beta_j-y_j)]_p\): the coupling
survives at the mode too, in general, and vanishes only when group \(j\)
carries **no likelihood information** in dimension \(p\) (the
\([D_j'\Omega_j(\cdot)]_p = 0\) case) — the same "data-free block" failure
mode already discussed for \(C_b=0\) in §2 of `BLOCK_GIBBS_ERGODICITY.md`,
now recurring one level up in the \(\lambda\)-coupling.

---

## 9. What this does *not* yet give (still future work)

This note supplies the algebra; it deliberately stops short of:

- **Cross-chain covariance capture.** §6-7 need the *cross-chain covariance*
  between precision components, and between precision and \(\gamma_{\mathrm{full}}\),
  to build a whitened/eigen-decomposed diagnostic analogous to
  `whitened = TRUE` for `component = "fixef"`. Only the marginal variance of
  each \(1/\tau_p^2\) is currently captured (`disp_table`).
- **A certified, re-derived TV bound.** As discussed in §7, Nygren's exact
  finite-\(l\) bounds do not carry over to a state-dependent \(A(\cdot)\)
  without new theory (or a different proof strategy, e.g. a
  drift/minorization argument for the nonlinear kernel, which is a
  meaningfully different undertaking from re-scaling Nygren's Gaussian
  result).
- **An implementation.** No `.two_block_S_P11_ing()`-style helper has been
  written; §6 is written so that, if implemented later, it would extend
  `.two_block_rate_inputs()`/`.two_block_S_P11()` by (a) appending the ING
  \(\lambda\)'s to `gamma_cols`/`gamma_names`, (b) extending each group's
  \(H_j\) with a residual row per ING component, and (c) plugging the current
  \(\lambda\) into \(B_j\) and \(P_{11}^{(\gamma\gamma)}\) in place of the
  fixed \(P_b\) used today. The needed \(u_{jp}\) residuals are already
  computed as a byproduct of the existing Block~2 ING update.
- **A marginal (Omega-integrated-out) alternative to the joint extension.**
  §16 derives, for the \(\Omega_j\)-only case, an exact PSD/log-concavity
  criterion on the residuals by integrating \(\Omega_j\) out analytically
  instead of extending the joint Hessian — a genuinely different diagnostic
  from §5-8/§14 above, and (like the rest of this list) not implemented.

---

## 10. A certified worst-case bound for the \((\gamma,\beta)\)-given-\((\lambda,\Omega)\) sub-chain

Sections 7-9 concluded that the *fully joint* \((\gamma,\beta,\lambda)\)
chain has a genuinely state-dependent rate matrix \(A(\gamma,\beta,\lambda)\),
with no certified global TV bound. This section proves a **weaker but
rigorous and immediately useful** claim: because the package's own prior
specification already *bounds* \(\lambda\) from above and (for an
ING-estimated measurement dispersion) \(\Omega\) from below, plugging in
those bounds to the **original, fixed-\((\lambda,\Omega)\)** theory gives a
**valid worst-case rate for the \((\gamma,\beta)\) sub-chain, uniformly over
every value \((\lambda,\Omega)\) the sampler can ever actually visit.**

### 10.1 This is already implemented — this section supplies the missing proof

- **RE-precision side.** `dIndependent_Normal_Gamma()` requires a finite
  `disp_lower` (§2 `disp_lower`/`disp_upper` truncation window on \(\tau_p^2\)
  for every ING component), and `two_block_rate_from_pfamily_list()` already
  plugs in `disp_lower` (not the prior mean) for every ING component:
  *"for `dIndependent_Normal_Gamma` components the conservative `disp_lower`
  plug-in is used (the rate is then an upper bound over the truncated
  `tau^2` range)"* (`R/two_block_ergodicity.R`).
- **Measurement-dispersion side.** When the observation-level dispersion is
  itself ING-estimated, `.rLMM_measurement_disp_upper_for_rate()`
  (`R/rLMM_reg.R`) plugs in `disp_upper` (not the prior mean `shape/rate`
  used for ICM) specifically "for rate calibration" — i.e. exactly the
  \(\Omega \succeq \Omega_{\min} = 1/\text{disp\_upper}\) bound this section
  needs.
- **Combined, plus an empirical safeguard.** `.two_block_pilot_ub_from_coefficients()`
  (`R/rLMM_reg.R`, `R/two_block_glmm_pilot_helpers.R`) calls
  `two_block_rate_from_pfamily_list()` with the conservative `dispersion`
  plug-in *at every pilot draw's `b_mode`* and keeps the **maximum**
  observed `lambda_star` across pilot chains (`pmax(max_eigenvalues, ...)`,
  `lambda_star_best`) — i.e. it already combines the analytic bound with an
  empirical max-over-realized-states safeguard, then certifies
  `m_min_upper` from that maximum via `two_block_l_for_tv()`.

What was missing — and is the point of this section — is a proof that this
practice is *sound*: that no value of \((\lambda,\Omega)\) the sampler could
visit gives a **worse** rate than the one computed at the bounds.

### 10.2 The claim

Fix the design (\(D,\mathcal{W},V\)) and consider the **original** (Nygren
2020 / `two_block_rate()`) rate matrix
\(A(\lambda,\Omega) = P_{11}^{-1/2}(\lambda)\,S(\lambda,\Omega)\,P_{11}^{-1/2}(\lambda)\)
for *fixed* \(\Psi^{-1}=\mathrm{diag}(\lambda_1,\dots,\lambda_P)=:\Lambda\) and
*fixed* likelihood precision \(\Omega\) (block-diagonal over groups) — i.e.
exactly what `two_block_rate()` computes for a plugged-in \(\Lambda,\Omega\),
**not** the extended \((\gamma,\beta,\lambda)\)-Hessian of §5-7.

**Claim.** \(\lambda^*(\Lambda,\Omega) := \lambda_{\max}(A(\Lambda,\Omega))\)
is, for every \(p\), non-decreasing in \(\lambda_p\) (holding \(\Omega\) and
every \(\lambda_{p'}\), \(p'\ne p\), fixed) and Loewner-non-increasing in
\(\Omega\) (holding \(\Lambda\) fixed). Consequently, for
\(\lambda_p \in [\lambda_p^{\min},\lambda_p^{\max}]\ \forall p\) and
\(\Omega \succeq \Omega_{\min}\):

$$
\lambda^*(\Lambda,\Omega) \;\le\; \lambda^*(\Lambda^{\max},\Omega_{\min}), \qquad
\Lambda^{\max} := \mathrm{diag}(\lambda_1^{\max},\dots,\lambda_P^{\max})
$$

i.e. evaluating at the corner of the box (every precision at its own upper
truncation bound, the likelihood precision at its lower bound) gives a valid
upper bound on the rate for **every** \((\Lambda,\Omega)\) consistent with
the priors — exactly what `two_block_rate_from_pfamily_list()` +
`.rLMM_measurement_disp_upper_for_rate()` compute.

### 10.3 Proof: monotonicity in \(\Omega\) (general, complete)

Recall from §6 (dropping \(\lambda\) as a variable here — this is the
*original*, fixed-\(\Lambda\) theory) that \(P_{11}\) does **not** depend on
\(\Omega\) at all: \(P_{11} = \mathcal{W}'\Lambda_{\text{full}}\mathcal{W} + V^{-1}\)
(\(\Lambda_{\text{full}} = I_J\otimes\Lambda\)). Marginalizing \(\beta\) out
of the joint Gaussian gives, exactly (not approximately — this is a linear
Gaussian marginalization identity),

$$
M(\Omega) \;=\; V^{-1} + \underbrace{\mathcal{W}'D'\big(D\Lambda_{\text{full}}^{-1}D' + \Omega^{-1}\big)^{-1}D\mathcal{W}}_{K(\Omega)}
$$

the Schur-complement / marginal precision of \(\gamma\), with
\(\lambda^*(\Omega) = 1 - \lambda_{\min}\big(M(\Omega), P_{11}\big)\) (generalized
eigenvalue of the pencil \((M,P_{11})\); recall \(A = I - P_{11}^{-1/2}MP_{11}^{-1/2}\)
since \(M = P_{11}-S\)).

If \(\Omega_1 \preceq \Omega_2\) (Loewner), then \(\Omega_1^{-1}\succeq\Omega_2^{-1}\)
(matrix inversion reverses the Loewner order on the PD cone), so
\(D\Lambda_{\text{full}}^{-1}D'+\Omega_1^{-1} \succeq D\Lambda_{\text{full}}^{-1}D'+\Omega_2^{-1}\),
so their inverses reverse again, and congruence by \(\mathcal{W}'D'(\cdot)D\mathcal{W}\)
preserves the order: \(K(\Omega_1)\preceq K(\Omega_2)\), hence
\(M(\Omega_1)\preceq M(\Omega_2)\). For any fixed unit \(v\),
\(v'M(\Omega_1)v \le v'M(\Omega_2)v\), so
\(v'M(\Omega_1)v/v'P_{11}v \le v'M(\Omega_2)v/v'P_{11}v\) for every \(v\)
(only the numerator moves; \(P_{11}\) is fixed), hence the minimum over
\(v\) satisfies \(\lambda_{\min}(M(\Omega_1),P_{11}) \le \lambda_{\min}(M(\Omega_2),P_{11})\)
(the min of a pointwise-smaller function is itself smaller — no envelope
theorem or extra regularity needed). Therefore
\(\lambda^*(\Omega_1) = 1-\lambda_{\min}(M(\Omega_1),P_{11}) \ge 1-\lambda_{\min}(M(\Omega_2),P_{11}) = \lambda^*(\Omega_2)\):
**\(\lambda^*\) is Loewner-non-increasing in \(\Omega\).** \(\blacksquare\)

This holds for **any** \(D\), any correlation structure in \(\Omega\)
(heterogeneous per-group/per-observation precisions included), and any
\(\Lambda\) — it is not restricted to diagonal or scalar \(\Omega\). It is
exactly why "more/better data can only help mixing" (§3 of
`BLOCK_GIBBS_ERGODICITY.md`) generalizes cleanly, and it fully justifies
`.rLMM_measurement_disp_upper_for_rate()`'s use of `disp_upper` (the least
favorable, most diffuse admissible measurement dispersion).

### 10.4 Proof: monotonicity in \(\lambda_p\) (complete for \(p_{re}=1\); cited for \(p_{re}>1\))

Write \(\lambda^* = \max_v v'Sv/v'P_{11}v\), \(v'P_{11}v = \sum_{p'} \lambda_{p'} c_{p'}(v) + \kappa(v)\)
(linear, non-decreasing in each \(\lambda_{p'}\); \(c_{p'}(v)=v_{p'}'W_{p'}'W_{p'}v_{p'}\ge0\)),
and \(v'Sv = \sum_j (\mathcal{W}_jv)'\,g_j(\Lambda)\,(\mathcal{W}_jv)\) with
\(g_j(\Lambda) = \Lambda(\Lambda + C_j)^{-1}\Lambda\), \(C_j := D_j'\Omega_jD_j\)
(the per-group likelihood precision restricted to \(\beta\)-space; this is
exactly `.two_block_S_P11()`'s `C_j <- P_b %*% solve(B_j, P_b)`, generalized
to a state-dependent \(\Lambda\)).

**\(p_{re}=1\) case (single RE dimension per group, e.g. random-intercept-only
models).** Here \(\Lambda=\lambda\) and \(C_j\) are both **scalars**, so
\(g_j(\lambda) = \lambda^2/(\lambda+C_j)\), and directly
\(g_j'(\lambda) = \lambda(2C_j+\lambda)/(\lambda+C_j)^2 > 0\) for \(\lambda,C_j>0\).
So \(v'Sv\) is strictly increasing in \(\lambda\) for every \(v\), *and*
\(v'P_{11}v\) is non-decreasing in \(\lambda\) with the same sign — but we only
need the **ratio**'s numerator to move favorably relative to a
*non-decreasing* denominator applied *pointwise in \(v\)*, which is exactly
Section 3 of `BLOCK_GIBBS_ERGODICITY.md`'s scalar formula
\(r = A/(A+C_b)\) (there \(A\leftrightarrow\lambda\)): a direct calculation
(\(dr/dA = C_b/(A+C_b)^2 > 0\)) already proves this case, and the multivariate
Rayleigh-quotient machinery above just re-derives the same scalar fact
generically per group before summing. \(\lambda^*\) is therefore rigorously
non-decreasing in \(\lambda\) whenever \(p_{re}=1\). \(\blacksquare\)

**\(p_{re}>1\) case (correlated RE dimensions, e.g. random intercept + slope).**
\(g(\Lambda)=\Lambda(\Lambda+C)^{-1}\Lambda\) for a general (non-diagonal)
\(C\succ0\) is exactly the classical **information-filter / Riccati precision
update** studied in Kalman filtering and Riccati-equation theory: it is the
precision contributed to a linear-Gaussian estimate after combining a "prior"
precision \(\Lambda\) with a "measurement" precision \(C\), then marginalizing
one side out — the same object whose *operator monotonicity* (non-decreasing
in \(\Lambda\) along any Loewner/diagonal-coordinate direction, for fixed
\(C\succ0\), and correspondingly non-increasing in \(C\) for fixed \(\Lambda\))
is a standard, well-established result in that literature (see e.g. Anderson
& Moore, *Optimal Filtering*, on monotonicity of the Riccati recursion in the
noise covariances; equivalently, \(g(\Lambda) = \Lambda - (\Lambda^{-1}+C^{-1})^{-1}\)
in terms of the operator-monotone matrix **parallel sum** \((\Lambda^{-1}+C^{-1})^{-1}\)
of Anderson & Duffin (1969)). This note does not re-derive that operator-monotonicity
result from scratch for general non-commuting \((\Lambda,C)\); it is cited
here as the fact needed to extend the \(p_{re}=1\) proof above to the general
diagonal-\(\Lambda\), correlated-\(C_j\) case. Given it, the same "pointwise-in-\(v\),
coordinatewise-in-\(\lambda_p\)" argument as §10.3 applies verbatim (each term
\((\mathcal{W}_jv)'g_j(\Lambda)(\mathcal{W}_jv)\) is non-decreasing in each
\(\lambda_p\) for fixed \(v\), fixed other \(\lambda_{p'}\), fixed \(\Omega\)),
giving the same "min/max of pointwise-monotone functions is itself monotone"
conclusion.

### 10.5 Assembling the box bound

§10.3-10.4 show \(\lambda^*(\Lambda,\Omega)\) is monotone in **each**
coordinate of \((\lambda_1,\dots,\lambda_P,\Omega)\) separately, for **every**
fixed value of the others in the box — not just at one reference point. A
function monotone coordinatewise everywhere on a box attains its supremum at
a corner, so

$$
\sup_{\lambda_p \in [\lambda_p^{\min},\lambda_p^{\max}],\ \Omega\succeq\Omega_{\min}} \lambda^*(\Lambda,\Omega)
\;=\; \lambda^*(\Lambda^{\max}, \Omega_{\min})
$$

exactly the corner `two_block_rate_from_pfamily_list()` (`disp_lower` per
ING component \(\Rightarrow \lambda_p^{\max}\)) and
`.rLMM_measurement_disp_upper_for_rate()` (`disp_upper` \(\Rightarrow \Omega_{\min}\))
already evaluate. `.two_block_pilot_ub_from_coefficients()`'s
max-over-pilot-draws is then a (harmless, conservative) *additional* margin
on top of this already-valid analytic bound — it does not by itself
certify anything beyond what §10.2-10.4 already prove, but it does guard
against the bound derivation above having missed some other source of
state-dependence (e.g. non-Gaussian `weights` computed at a specific
`b_mode`, which is a local-Gaussian heuristic, not itself covered by Nygren's
exact theorem).

### 10.6 Scope: what this does, and does not, give

- **What it gives:** a *certified*, non-asymptotic statement that the
  \((\gamma,\beta)\) two-block sub-chain's rate is at most
  \(\lambda^*(\Lambda^{\max},\Omega_{\min})\), for **any** value of
  \((\lambda,\Omega)\) actually realized during sampling (since the ING
  truncation windows force every sampled \(\lambda_p,\Omega\) to lie in the
  box by construction). This is exactly the quantity needed to pick a safe
  `m_convergence` for the \(\beta\)-given-\((\gamma,\lambda)\) /
  \(\gamma\)-given-\((\beta,\lambda)\) alternation *conditional on whatever
  \(\lambda\) is currently plugged in* — precisely the pilot-calibration use
  case `.two_block_pilot_ub_from_coefficients()` serves.
- **What it does not give:** a bound on how fast \(\lambda\) *itself* mixes,
  nor on the full **joint** \((\gamma,\beta,\lambda)\) chain's TV distance to
  stationarity (the open problem of §7-9). Freezing \(\lambda\) at a
  worst-case corner and bounding the resulting *linear* sub-chain says
  nothing about the *nonlinear* coupling terms (the \(\gamma\)-\(\lambda\) and
  \(\beta\)-\(\lambda\) Hessian blocks of §5-6) that make the full chain's
  actual rate matrix state-dependent in the first place. The two results are
  complementary, not substitutes: §10 certifies the sub-chain conditional on
  \(\lambda\); §7-9 explain why no analogous global certificate yet exists
  for \(\lambda\)'s own dynamics or the fully joint chain.

---

## 11. Freezing the residuals at zero gives the *smallest*, not largest, eigenvalue

A natural question about §5-7's genuinely new blocks — \(P_{11}^{(\gamma\lambda)} = -W_p'u_p\)
and \(P_{12}\)'s \(\lambda\)-rows \(= u_{jp}\), both linear in the residuals
\(u_{jp} = \beta_{jp}-W_{pj}\gamma_p\) — is whether setting them to zero
(the *only* natural "no correlation" reference point, since there is no
prior-truncation bound on \(u\) the way §10 has for \(\lambda,\Omega\)) gives
a worst case, analogous to evaluating \(\Lambda,\Omega\) at their bounds.
**It does not — it gives the exact opposite extreme.**

**Claim.** For fixed \(\Lambda,\Omega\), \(\lambda_{\max}(A(u))\) — the top
eigenvalue of the *extended* rate matrix as a function of the residual
vector \(u = (u_{jp})_{j,p\in\mathrm{ING}}\) — attains its **minimum** at
\(u=0\), where it equals **exactly** \(\lambda^*_{\text{original}}(\Lambda,\Omega)\)
(§10's already-bounded quantity, with the \(\lambda\)-directions contributing
nothing). For every other \(u\), \(\lambda_{\max}(A(u)) \ge \lambda^*_{\text{original}}(\Lambda,\Omega)\).

**Proof.** Write \(x_1=(\gamma,\lambda)\), \(x_2=\beta\), and use
\(\lambda_{\max}(A(u)) = \max_{z=(z_\gamma,z_\lambda)\ne0}\, z'S(u)z / z'P_{11}(u)z\),
\(S(u):=P_{12}(u)P_{22}^{-1}P_{21}(u)\) (\(P_{22}\) does not depend on \(u\)
at all — only on \(\Lambda,\Omega\)). Restrict the max to the subspace
\(z_\lambda = 0\):

- **Numerator.** \(P_{12}(u)\)'s \(\gamma\)-rows are \(-\Lambda\mathcal{W}\)
  (fixed) and its \(\lambda\)-rows are \(u\); with \(z_\lambda=0\) the
  \(\lambda\)-rows are never touched, so
  \(z'S(u)z = z_\gamma'\big(P_{12,\gamma}P_{22}^{-1}P_{12,\gamma}'\big)z_\gamma\)
  — **exactly** the original (non-extended) \(S_{\gamma\gamma}\), independent of \(u\).
- **Denominator.** Likewise \(z'P_{11}(u)z = z_\gamma'P_{11}^{(\gamma\gamma)}z_\gamma\)
  (the \(-W_p'u_p\) cross-block is never touched since \(z_\lambda=0\)) —
  independent of \(u\).

So on \(\{z_\lambda=0\}\) the ratio is **identically**
\(\lambda^*_{\text{original}}(\Lambda,\Omega)\) for *every* \(u\). Two
consequences: (i) at \(u=0\), any \(z_\lambda\ne0\) strictly *hurts* — the
numerator gains nothing from \(z_\lambda\) (the \(\lambda\)-rows of
\(P_{12}(0)\) are literally zero) while the denominator strictly gains
\(z_\lambda'P_{11}^{(\lambda\lambda)}z_\lambda>0\) — so the unrestricted max at
\(u=0\) is attained at \(z_\lambda=0\), giving
\(\lambda_{\max}(A(0)) = \lambda^*_{\text{original}}(\Lambda,\Omega)\) exactly.
(ii) For *any* \(u\), the unrestricted max (over all \(z\)) is a max over a
superset of \(\{z_\lambda=0\}\), hence \(\ge\) the restricted value
\(\lambda^*_{\text{original}}(\Lambda,\Omega)\), which is the same for every
\(u\). \(\blacksquare\)

**Interpretation.** The residual terms are exactly what §7 already said they
are — the *only* channel through which \(\lambda\) is correlated with
\((\gamma,\beta)\) in this local picture. At \(u=0\) that channel is silent
and \(\lambda\) decouples completely (its own directions contribute zero to
the top eigenvalue); switching it on can only *add* correlation where there
was none, consistent with the theory's central theme (more cross-block
correlation \(\Rightarrow\) slower mixing, §3 of `BLOCK_GIBBS_ERGODICITY.md`).
So \(u=0\) is the *best* case, not the worst one.

**Consequence for §10.** Freezing \(u=0\) is **not** a valid way to extend
the §10 box bound to cover the new residual-driven blocks — doing so would
*understate* the true worst-case rate for the extended system. A genuine
worst-case bound accounting for \(u\) would need an upper bound on the size
of \(u\) (or the relevant quadratic forms it enters), and unlike \(\Lambda,\Omega\)
there is no prior-truncation window that bounds \(u\) directly. §12 outlines
an empirical safeguard for this, mirroring the package's existing treatment
of the analogous problem for non-Gaussian likelihoods.

---

## 12. Proposed empirical safeguard for \(u\), mirroring the non-Gaussian pilot machinery

The package already faces a structurally identical problem for non-Gaussian
families: the likelihood precision \(\Omega\) is not fixed either (Fisher
information depends on the current fit), and there is no way to bound it
from prior truncation the way \(\Lambda\) is bounded. It is handled with two
mechanisms, which this section proposes reusing verbatim for \(u\) (not
implemented here — an outline only, per §9's scope).

| | Non-Gaussian \(\Omega\) (existing) | ING residual \(u\) (proposed) |
|---|---|---|
| Unknown quantity | \(\Omega\) (state-dependent Fisher info) | \(u_{jp}\) (no prior bound exists) |
| Mechanism 1: plug-in at a state | `two_block_mode_weights()` evaluates IRLS/Fisher weights at a specific `b_mode` | Compute \(u_{i,jp} = b_{i,jp} - W_{pj}\gamma_{i,p}\) at a specific pilot draw \((\gamma_i,b_i)\) |
| Mechanism 2: worst-of-pilot search | `.two_block_pilot_eigenvalue_ub()` loops over `n_pilot` draws, `pmax`s the eigenvalue spectrum | Same loop, calling the *extended* rate (§6-7) at each \(u_i\) instead of the plain rate |
| Certification | `two_block_l_for_tv()` on the worst observed spectrum \(\to\) `m_min_upper` | Same, on the extended worst spectrum \(\to\) `m_min_upper_ing` |
| Caveat already on record | "local-Gaussian heuristic... no theorem applies" (`two_block_rate()` docs) | Same caveat, **plus** an additional one: unlike \(\Omega\), \(u\) has *no* prior-enforced ceiling, so even this heuristic is only as good as the pilot sample, with no fallback certified bound the way `disp_lower`/`disp_upper` still back the \(\Omega\)-heuristic |

**Step 1 (plug-in, cheaper than the \(\Omega\) case — exact arithmetic, no
family-specific formula).** `pilot_raw$fixef_draws[[p]]` (an
\(n_{\text{pilot}}\times q_p\) matrix; row \(i\) is draw \(i\)'s \(\gamma_{i,p}\))
and `pilot_raw$coefficients` (row-blocked by draw, giving \(b_i\)) are
**already both collected** by the same pilot chains feeding
`.two_block_pilot_ub_from_coefficients()` today — no new sampling is needed,
just \(u_{i,jp} = b_i[j,p] - W_p[j,]\cdot\gamma_i[[p]]\) per draw, per ING
component.

**Step 2 (extended local rate at that draw).** Build \(P_{11}(u_i)\),
\(P_{12}(u_i)\) (§6, with \(\Lambda\) pinned at \(\Lambda^{\max}\) and
\(\Omega\) at \(\Omega_{\min}\)/local Fisher weights — §10's already-certified
corner, now perturbed by the residual blocks) and compute
\(\lambda_{\max}(A(u_i))\) plus its full spectrum via the
`.two_block_S_P11_ing()`-style helper sketched in §9.

**Step 3 (search over pilot draws).** Loop \(i=1,\dots,n_{\text{pilot}}\);
track `max_eigenvalues_ing <- pmax(max_eigenvalues_ing, eigenvalues(A(u_i)))`
and the argmax draw index — the same accumulation pattern as
`.two_block_pilot_eigenvalue_ub()`'s existing
`pmax(max_eigenvalues, rate_i$eigenvalues)`.

**Step 4 (certify `m_convergence`).** Feed the worst observed spectrum into
`two_block_l_for_tv()` \(\to\) `m_min_upper_ing`; use
`m_convergence_used <- max(m_min_upper, m_min_upper_ing)` for the main run,
combining with (not replacing) §10's certified box bound.

This closes the gap left by §11 with a concrete, implementable design that
reuses existing pilot infrastructure end to end, but it should be documented
with a caveat *stronger* than the existing \(\Omega\)-heuristic's: because
\(u\) has no prior-enforced ceiling, this mechanism is a practical safeguard
against an *observed* worst case, not a certificate against every possible
one, however small the resulting spectrum may empirically look.

---

## 13. Possible future implementation: enforcing a *hard* bound on \(u\) via truncated sampling

§12's safeguard is empirical because nothing in the model actually *prevents*
\(u\) from exceeding whatever was seen in the pilot. A more ambitious (and
not yet attempted) alternative — noted here as a design idea, not
implemented — would make the bound on \(u\) genuinely **enforced by the
sampler**, the same way `disp_lower`/`disp_upper` already enforce a hard
ceiling on \(\lambda_p\) rather than just observing one empirically.

### 13.1 The idea

Conditional on \((\gamma_p,\lambda_p)\), and independent across groups
\(j\) and (given the "Independent" Normal-Gamma structure) across ING
components \(p\), \(u_{jp} = \beta_{jp} - W_{pj}\gamma_p \sim N(0, 1/\lambda_p)\)
exactly — this is just the RE prior itself, re-centered. So:

1. **Bound \(u_{jp}\) directly**, e.g. to a symmetric window
   \(u_{jp} \in [-c_p, c_p]\), with \(c_p\) derived from a quantile of
   \(N(0,\tau_p^{2,\max})\) where \(\tau_p^{2,\max} = \text{disp\_upper}_p\) is
   the RE component's **already-existing** truncation upper bound (the
   *widest* admissible conditional spread of \(u_{jp}\)) — reusing exactly
   the quantile-window construction `.lmebayes_ing_prior_quantile_window()`
   already applies to \(\tau_p^2\) itself, just applied to the *residual*
   distribution instead. **No new user-facing parameter is needed.**
2. **Remap to \(\beta_{jp}\).** Since \(u_{jp}=\beta_{jp}-W_{pj}\gamma_p\), the
   bound on \(u_{jp}\) is exactly a bound on \(\beta_{jp}\):
   \(\beta_{jp} \in [W_{pj}\gamma_p - c_p,\ W_{pj}\gamma_p + c_p]\) — a window
   of **fixed width** \(2c_p\) whose *center* tracks the current
   \(\gamma_p\) (and the specific group's hyper-design row \(W_{pj}\)) every
   sweep. This is precisely "restrictions on the random effects that are
   dependent on the fixed-effect predictors."
3. **Enforce it by truncating the \(\beta_j\) full-conditional draw** (which
   is already an exact Gaussian, §4 of `inst/notation.md`) to this window,
   componentwise per \(p\in\mathrm{ING}\). The package already implements
   inverse-CDF truncated draws for exactly this kind of problem —
   `.rLMM_ing_sample_sigma2()`'s truncated-Gamma draw via
   `qgamma(F_lo + runif()*(F_hi-F_lo), ...)` is the same pattern applied to
   \(\tau_p^2\); a truncated-normal analogue (`qnorm(pnorm(lo) + runif()*(pnorm(hi)-pnorm(lo)))`)
   is simpler, not harder.

### 13.2 What this would buy, and what it would still need

If implemented, this converts \(u\) from an *unbounded* nuisance quantity
into a **genuinely truncation-enforced** one, exactly like \(\lambda,\Omega\)
already are — removing §12's "no ceiling" caveat and making a `disp_lower`/
`disp_upper`-style hard box for \(u\) a real property of the chain, not an
empirical observation about one pilot run.

It does **not**, by itself, hand back a certified worst-case bound the way
§10 does for \((\Lambda,\Omega)\). §11 only proves \(u=0\) is the *global
minimizer* of \(\lambda_{\max}(A(u))\) — it does not establish that
\(\lambda_{\max}(A(u))\) is *monotone in the size* of \(u\) within the
window, which is the extra fact §10.5's "coordinatewise-monotone-on-a-box
attains its extremum at a corner" argument would need to conclude that the
window's boundary (\(u_{jp}=\pm c_p\) everywhere) is the worst case within
the box, rather than some interior point. That monotonicity is plausible
(§11's mechanism — more coupling is worse — suggests \(|u|\) larger should
generally be worse) but is not proven here, and would be the natural next
piece of theory to attempt before treating a hard \(u\)-bound as a full
substitute for §12's pilot search rather than a complement to it.

### 13.3 Implementation risk to flag

A truncation window that is too tight would **bias** the posterior (the
whole point of the current `disp_lower`/`disp_upper` windows is to be wide
enough — a 98-99% prior-mass window — that they essentially never bind in a
well-behaved run; the same discipline would need to apply here). The
truncation is on \(u_{jp}\), not on \(\beta_{jp}\) directly, so as
\(\gamma_p\) moves across sweeps the effective window for \(\beta_{jp}\)
itself moves with it — the constraint should be understood as "how far a
random effect may deviate from its current hyper-mean," not as a fixed
absolute range for \(\beta_{jp}\).

---

## 14. The \(\Omega\) (measurement-precision) extension

Sections 5-13 extend the Hessian for a sampled **random-effects** precision
\(\lambda_p\) (`dIndependent_Normal_Gamma()` on \(\tau_p^2\)). This section
derives the exactly analogous extension for a sampled **measurement**
(observation-level) precision \(\Omega := 1/\sigma^2\) — pooled (`dGamma()`,
one value shared by every group) or per-group \(\Omega_j := 1/\sigma_j^2\)
(`dGamma_list()`, one independent value per group, as used by
`rLMMindepNormalGamma_reg_*()`). Everywhere below, \(\Omega\) means "pooled
\(\Omega\)" unless a subscript \(j\) is present.

### 14.1 Setup

Let \(e_j := y_j - D_j\beta_j\) (the \(n_j\times1\) **data** residual for group
\(j\) — not to be confused with \(u_{jp}\), the *RE-level* residual of §1-2).
Write \(n_j = \dim(y_j)\), \(n = \sum_j n_j\). The likelihood term of §3
becomes, once \(\Omega\) (or \(\Omega_j\)) is written out explicitly instead
of held fixed:

**Pooled case** (\(\Omega_j = \Omega I_{n_j}\) for every \(j\), one shared
scalar \(\Omega\), prior \(\Omega\sim\mathrm{Gamma}(a^0,r^0)\)):

$$
\ell \;\supset\; \frac{\Omega}{2}\sum_j e_j'e_j \;-\; \frac{n}{2}\log\Omega
\;+\; (1-a^0)\log\Omega + r^0\Omega
$$

**Per-group case** (\(\Omega_j = \Omega_j I_{n_j}\), independent
\(\Omega_j\sim\mathrm{Gamma}(a_j^0,r_j^0)\) per group):

$$
\ell \;\supset\; \sum_j\Big[\frac{\Omega_j}{2} e_j'e_j - \frac{n_j}{2}\log\Omega_j
+ (1-a_j^0)\log\Omega_j + r_j^0\Omega_j\Big]
$$

Both reduce to the fixed-\(\Omega\) likelihood term of §3 when \(\Omega\)
(resp. every \(\Omega_j\)) is held constant — this section only makes its
dependence explicit, exactly as §3-5 did for \(\lambda_p\).

### 14.2 First derivatives

$$
\frac{\partial\ell}{\partial\beta_{jp}} = \big[D_j'\Omega_j(D_j\beta_j-y_j)\big]_p + \lambda_p u_{jp}
\;=\; -\big[D_j'\Omega_je_j\big]_p + \lambda_p u_{jp}
$$

(same formula as §4, just substituting \(D_j\beta_j - y_j = -e_j\); \(\Omega_j\)
is \(\Omega\) in the pooled case). For the precision itself,

$$
\frac{\partial\ell}{\partial\Omega} = \frac12\sum_j e_j'e_j - \frac{n}{2\Omega} + \frac{1-a^0}{\Omega} + r^0
\qquad\text{(pooled)}
$$

$$
\frac{\partial\ell}{\partial\Omega_j} = \frac12 e_j'e_j - \frac{n_j}{2\Omega_j} + \frac{1-a_j^0}{\Omega_j} + r_j^0
\qquad\text{(per-group)}
$$

Setting either to zero reproduces the exact conjugate Gamma full conditional
\(\Omega \mid \beta \sim \mathrm{Gamma}(a^0+n/2,\ r^0+\tfrac12\sum_je_j'e_j)\)
(pooled) or \(\Omega_j\mid\beta_j\sim\mathrm{Gamma}(a_j^0+n_j/2,\ r_j^0+\tfrac12e_j'e_j)\)
(per-group) — the same `shape2 = shape + n_w/2.0`-style update used by
`dGamma()`/`dGamma_list()`'s Block~1 Gibbs step, and the same sanity-check
pattern used for \(\lambda_p\) in §4.

### 14.3 Second derivatives, and why the cross-terms with \((\gamma,\lambda)\) vanish exactly

| Block | Formula | Why (non)zero |
|---|---|---|
| \(\partial^2\ell/\partial\beta_{jp}\partial\Omega\) (pooled) | \(-[D_j'e_j]_p\) | \(\partial(-[D_j'\Omega e_j]_p)/\partial\Omega\), holding \(e_j\) fixed; nonzero for **every** group \(j\) (pooled \(\Omega\) enters every group's likelihood) |
| \(\partial^2\ell/\partial\beta_{jp}\partial\Omega_{j'}\) (per-group), \(j = j'\) | \(-[D_j'e_j]_p\) | same derivation, restricted to group \(j\)'s own precision |
| \(\partial^2\ell/\partial\beta_{jp}\partial\Omega_{j'}\) (per-group), \(j\ne j'\) | \(0\) | group \(j\)'s likelihood term does not involve \(\Omega_{j'}\) |
| \(\partial^2\ell/\partial\gamma_p\partial\Omega\) (or \(\Omega_j\)) | \(0\) | \(\gamma_p\) only appears in the RE-prior/Block~2 terms (§3), which do not involve \(\Omega\) at all |
| \(\partial^2\ell/\partial\lambda_p\partial\Omega\) (or \(\Omega_j\)) | \(0\) | \(\lambda_p\) only appears in the RE-prior term (§3); \(\Omega\) only in the likelihood term — the two never multiply the same additive piece of \(\ell\) |
| \(\partial^2\ell/\partial\Omega^2\) (pooled) | \((n/2+a^0-1)/\Omega^2\) | same derivation as §5's \(\lambda_p^2\) case, with \(n\) (not \(J\)) as the log-Jacobian count, since the Gaussian normalizing constant here scales with the number of *observations*, not groups |
| \(\partial^2\ell/\partial\Omega_j^2\) (per-group) | \((n_j/2+a_j^0-1)/\Omega_j^2\) | same, per group |
| \(\partial^2\ell/\partial\Omega_j\partial\Omega_{j'}\), \(j\ne j'\) (per-group) | \(0\) | independent Gamma priors, no shared likelihood term |

The two zero rows are the key structural fact used below: **\(\Omega\)
(pooled or per-group) is Hessian-orthogonal to both \(\gamma\) and
\(\lambda\)** — it interacts with the joint system *only* through its cross-
terms with \(\beta\). This is because \(\Omega\) multiplies only the
likelihood piece of \(\ell\), which contains \(\beta\) but never \(\gamma\)
or \(\lambda\) directly (those enter only via the RE-prior/Block~2 terms).

### 14.4 Assembling the extended blocks

Let \(x_1 = (\gamma,\ \lambda_{\mathrm{ING}},\ \Omega_{\mathrm{ING}})\) —
appending the \(\Omega\)-extension of this section onto the \(\lambda\)-
extension of §6 (§14.5 below covers stacking both explicitly; here consider
\(\Omega\) added to the *plain*, non-ING-\(\lambda\), \(x_1=\gamma\) system
for clarity first). \(P_{22}\) is **unchanged in form**
(\(B_j = D_j'\Omega_jD_j + \Psi^{-1}\); \(\Omega_j\) is simply now the current
state rather than a fixed plug-in, exactly as \(\lambda\) was in §6).

**Pooled \(\Omega\):** one new row/column. New \(P_{12}\) row, for every group
\(j\): \(-[D_j'e_j]'\) (a \(1\times p_{re}\) block, dense — nonzero in every
column of that group's \(p_{re}\)-wide strip, since pooled \(\Omega\) touches
every \(\beta_{jp}\)). New \(P_{11}\) diagonal entry:
\((n/2+a^0-1)/\Omega^2\). New \(P_{11}^{(\gamma,\Omega)}\) and
\(P_{11}^{(\lambda,\Omega)}\) cross blocks: **exactly zero** (§14.3).

**Per-group \(\Omega_j\):** \(J\) new rows/columns (one per group), each
**block-local**: row \(j\)'s new \(P_{12}\) entries are \(-[D_j'e_j]'\) within
group \(j\)'s own \(p_{re}\)-wide strip and **zero** in every other group's
strip (unlike the pooled case, which is dense across all \(J\) strips in its
single new row). New \(P_{11}\) diagonal entries:
\(\mathrm{diag}_j\big((n_j/2+a_j^0-1)/\Omega_j^2\big)\). Cross blocks with
\(\gamma\) and \(\lambda\): again exactly zero, and additionally
\(P_{11}^{(\Omega_j,\Omega_{j'})}=0\) for \(j\ne j'\) (§14.3) — the per-group
\(\Omega\)-block is itself diagonal.

### 14.5 Stacking both extensions

Because \(P_{11}^{(\gamma,\Omega)}=P_{11}^{(\lambda,\Omega)}=0\) **exactly**
(§14.3) — not approximately, and not only at a stationary point (contrast
with the \(\gamma\)-\(\lambda\) cross term of §8, which vanishes only at the
mode under a flat \(\gamma\)-prior) — the \(\lambda\)-extension (§5-6) and the
\(\Omega\)-extension (§14.1-14.4) are **Hessian-orthogonal** and can be
stacked with no further cross-derivation needed. For
\(x_1 = (\gamma,\ \lambda_{\mathrm{ING}},\ \Omega_{\mathrm{ING}})\) the
extended \(P_{11}\) is, in partitioned form,

```
        ┌ P11^(γγ)     P11^(γλ)      0        ┐
P11  =  │ P11^(γλ)'    P11^(λλ)      0        │
        └ 0            0             P11^(ΩΩ) ┘
```

i.e. the new \(\Omega\)-block is **block-diagonal** with the
\((\gamma,\lambda)\)-block of §6 — appended on the diagonal, contributing no
new off-diagonal structure beyond its own \(P_{12}\) rows (§14.4) into
\(\beta\)-space. This is exactly the structure `.two_block_S_P11_ing()`
(§15) exploits: the \(\lambda\)- and \(\Omega\)-extensions are assembled
independently and simply concatenated, never requiring a genuinely new
cross-derivation between them.

All of §7-13's conclusions (state-dependence, the §10 corner-bound
construction for the \((\gamma,\beta)\) sub-chain conditional on
\((\lambda,\Omega)\), and §11-13's residual-freezing/safeguard discussion) go
through for \(\Omega\) with \(u_{jp}\to e_j\) and \(\lambda_p\to\Omega\)
(or \(\Omega_j\)) throughout — no new argument is needed, since §14.3's exact
zero cross-terms are the only structural fact §7-13's arguments actually use
beyond what already held for \(\lambda\).

---

## 15. Implementation note: \code{.two_block_S_P11_ing()}

Implemented in `R/two_block_ergodicity_ing.R` as a **diagnostic-only**
extension of `.two_block_S_P11()`, evaluated once at a single reference
state \((\hat\gamma,\hat\beta,\hat\lambda,\hat\Omega)\) (not a pilot-draw
scan; see §12 for that still-unimplemented direction) — reusing
`.two_block_rate_inputs()`/`.two_block_S_P11()` unchanged for the base
\((\gamma,\beta)\) blocks (evaluated at whatever \((\Lambda,\Omega)\)
reference the caller supplies, typically the same one the corresponding
`rLMM*_reg_*()` engine's own pilot-rate calibration already uses) and
appending the §6/§14 blocks for whichever of \(\lambda_{\mathrm{ING}}\),
\(\Omega_{\mathrm{ING}}\) apply to the model at hand. `two_block_rate_ing()`
is the corresponding user-facing export (mirroring `two_block_rate()`'s
shape), documented as a **local, uncertified diagnostic** (§7, §9): it
reports how much the base rate would move if the (currently ignored)
\(\beta\)-\((\lambda,\Omega)\) coupling terms of §6/§14 were included, holding
everything else fixed — not a new certified bound. §16 below derives a
different (unimplemented) diagnostic for the \(\Omega_j\) case that
integrates \(\Omega_j\) out analytically instead of extending the joint
Hessian this section's helper computes.

---

## 16. Integrating \(\Omega_j\) out exactly: a marginal Hessian and a log-concavity criterion on the residuals

Sections 14-15 extend the *joint* \((\gamma,\lambda,\Omega,\beta)\) Hessian by
appending \(\Omega_j\) as its own coordinate and evaluating the cross terms at
a supplied reference \((\hat\beta_j,\hat\Omega_j)\) pair. This section takes a
different route: integrate \(\Omega_j\) out of the model **analytically**
instead, leaving a Hessian in \(\beta_j\) alone with no \(\Omega_j\)
coordinate at all, and ask when *that* Hessian is guaranteed positive
semi-definite (PSD) — i.e. when the marginal negative log-density is locally
log-concave in \(\beta_j\). Diagnostic/theoretical only, like §12-13: nothing
here is implemented in `R/two_block_ergodicity_ing.R`.

### 16.1 Setup: marginalizing instead of extending

From §14.1-14.2, \(\Omega_j\mid\beta_j\sim\mathrm{Gamma}\big(a_j^0+n_j/2,\;
r_j^0+\tfrac12e_j'e_j\big)\) exactly, with \(e_j=y_j-D_j\beta_j\). Integrating
\(\Omega_j\) out of the joint density \(p(y_j,\Omega_j\mid\beta_j)\) (the
standard Normal-Gamma-mixture identity) leaves a Student-t kernel in
\(\beta_j\):

$$
p(y_j\mid\beta_j)\;\propto\;\Big(r_j^0+\tfrac12e_j'e_j\Big)^{-(a_j^0+n_j/2)}
$$

i.e. \(y_j\mid\beta_j\) is (up to scale) multivariate-\(t\) with
\(2(a_j^0+n_j/2)\) degrees of freedom, location \(D_j\beta_j\), in place of
the Gaussian likelihood the fixed-\(\Omega_j\) system uses. This is a
genuinely different sampler/diagnostic variant from §14-15's joint extension
(that one keeps \(\Omega_j\) as a sampled state variable; this one removes it
from the model entirely by exact marginalization) — not a further extension
of it.

### 16.2 The marginal Hessian

Write

$$
\Omega_j^{\mathrm{eff}}(\beta_j)\;:=\;\frac{a_j^0+n_j/2}{r_j^0+\tfrac12e_j'e_j}
\;=\;E[\Omega_j\mid\beta_j]
$$

— literally the mean of the same Gamma full conditional §14.2 already
derives (the one `dGamma()`/`dGamma_list()`'s Block~1 Gibbs step samples
from), evaluated at \(\beta_j\) rather than at a fresh draw. Differentiating
\(-\log p(y_j\mid\beta_j)\) twice with respect to \(\beta_j\) gives the closed
form

$$
H_j(\beta_j)\;=\;\Omega_j^{\mathrm{eff}}D_j'D_j\;-\;
\frac{\big(\Omega_j^{\mathrm{eff}}\big)^2}{a_j^0+n_j/2}(D_j'e_j)(D_j'e_j)'
$$

(verified numerically against a brute-force finite-difference Hessian of
\(-\log p(y_j\mid\beta_j)\) on a synthetic design). The first term is exactly the current
fixed-\(\Omega_j\) plug-in \(\Omega_jD_j'D_j\) with \(\Omega_j\) replaced by
its conditional mean; the second is a rank-one, negative-semidefinite
correction that vanishes exactly at \(e_j=0\), reproducing §11's finding that
\(e_j=0\) (like \(u_{jp}=0\)) is the best-behaved reference state and that
moving away from it can only add negative curvature, never positive.

### 16.3 Exact log-concavity criterion

\(H_j(\beta_j)\) has the form \(M-cvv'\) with \(M=\Omega_j^{\mathrm{eff}}D_j'D_j\)
(PD, assuming \(D_j\) has full column rank \(p_{re}\)), \(v=D_j'e_j\), and
\(c=\big(\Omega_j^{\mathrm{eff}}\big)^2/(a_j^0+n_j/2)\). The standard
rank-one-perturbation fact (\(M-cvv'\) is PSD iff \(cv'M^{-1}v\le1\)) reduces,
after substituting \(\Omega_j^{\mathrm{eff}}\), to

$$
q_j\;:=\;(D_j'e_j)'(D_j'D_j)^{-1}(D_j'e_j)\;\le\;r_j^0+\tfrac12e_j'e_j
$$

Writing \(\hat\beta_j^{\mathrm{ols}}=(D_j'D_j)^{-1}D_j'y_j\) and using the
Pythagorean OLS decomposition \(e_j=\hat{e}_j^{\mathrm{ols}}+
D_j(\hat\beta_j^{\mathrm{ols}}-\beta_j)\) (with \(D_j'\hat{e}_j^{\mathrm{ols}}=0\)
by the OLS normal equations), \(q_j\) is exactly the \(D_j'D_j\)-weighted
squared distance from \(\beta_j\) to the group's own OLS fit, and
\(e_j'e_j=\mathrm{RSS}_j^{\mathrm{ols}}+q_j\). Substituting gives the
equivalent, more interpretable form:

$$
\big\|\beta_j-\hat\beta_j^{\mathrm{ols}}\big\|^2_{D_j'D_j}\;\le\;
2r_j^0+\mathrm{RSS}_j^{\mathrm{ols}}
$$

i.e. \(H_j(\beta_j)\) stays PSD (the marginal density stays locally
log-concave) exactly as long as \(\beta_j\) is not too far — in the
\(D_j'D_j\)-weighted sense — from the group's own OLS fit, with the allowed
distance growing with both the prior's informativeness (\(r_j^0\)) and how
poorly that group's OLS fit already explains its own data
(\(\mathrm{RSS}_j^{\mathrm{ols}}\)). Both criteria were checked numerically
against a direct eigenvalue computation of \(H_j\), scanning \(\beta_j\) away
from \(\hat\beta_j^{\mathrm{ols}}\): the sign of the smallest eigenvalue of
\(H_j\) matched the sign of \(r_j^0+\tfrac12e_j'e_j-q_j\) at every point
tested.

### 16.4 Two practical "buckets"

Since \(q_j\le e_j'e_j\) always (the residual's projection onto
\(\mathrm{col}(D_j)\) cannot exceed its own squared norm — \(H_{D_j}:=
D_j(D_j'D_j)^{-1}D_j'\) is an orthogonal projection), two simple corollaries
follow, requiring no knowledge of \(D_j\)'s leverage structure or of
\(\hat\beta_j^{\mathrm{ols}}\):

- **Design-independent bucket.** If \(e_j'e_j\le2r_j^0\) (the current
  residual sum of squares is at most twice the prior Gamma rate), then
  \(H_j(\beta_j)\) is PSD *regardless* of \(D_j\) or how \(e_j\) is oriented
  relative to \(\mathrm{col}(D_j)\).
- **Leverage bucket.** Writing \(\rho_j:=q_j/e_j'e_j\in[0,1]\) (the fraction
  of the residual's squared norm that lies in \(\mathrm{col}(D_j)\), i.e. \(0\)
  exactly at the OLS fit and growing as \(\beta_j\) moves away from it along
  directions well-explained by the design), \(\rho_j\le\tfrac12\) guarantees
  \(H_j(\beta_j)\) is PSD for *any* residual magnitude; only \(\rho_j>\tfrac12\)
  admits a finite counterexample, at \(e_j'e_j>r_j^0/(\rho_j-\tfrac12)\).

### 16.5 Relationship to the current implementation and the Ex_13 empirical finding

This is a *different* diagnostic angle from §14-15's joint extension, not a
refinement of it. §16.1-16.4's untruncated formulas are, since §16.6, present
in code only as a conservative comparison branch alongside the exact
truncated one (never as the primary criterion — see §16.6's "Practical
implication"). It does, however, explain a finding from the empirical "scan
every main-stage draw" diagnostic
(the `.two_block_rate_ing_over_draws()` helper used by
`demo("Ex_13_...")`/`demo("Ex_13b_...")`/`demo("Ex_14_...")`): that diagnostic
plugs each draw's independently-sampled \((\beta_j^{(i)},\Omega_j^{(i)})\)
pair into the *joint* extended Hessian of §14, rather than
\(\Omega_j^{\mathrm{eff}}(\beta_j^{(i)})\) — so \(\Omega_j^{(i)}\) is not
constrained to satisfy the §16.3 bound relative to that same draw's
\(e_j^{(i)}\) the way the marginal system's \(\Omega_j^{\mathrm{eff}}\) always
is by construction. A large fraction of draws showing \(\lambda^*\ge1\) in
that scan is consistent with — and partially explained by — evaluating an
uncontrolled joint reference state rather than the marginal one this section
derives.

### 16.6 Accounting for the actual *truncated* prior: exact marginal and Hessian in terms of \(\texttt{disp\_lower}\)/\(\texttt{disp\_upper}\)

**Status: implemented.** §16.1-16.5 integrate \(\Omega_j\) out assuming the
*untruncated* \(\mathrm{Gamma}(a_j^0,r_j^0)\) prior. But the prior that
`dGamma_list()` actually hands the sampler is truncated to a bounded window
\(\sigma_j^2\in[\texttt{disp\_lower}_j,\texttt{disp\_upper}_j]\) (an
accept/reject envelope over that range, not the full \((0,\infty)\) support),
which the untruncated derivation never accounts for.
This subsection re-derives §16.2's marginal and Hessian accounting for that
truncation exactly, in closed form, using only quantities already computed by
`Prior_Setup_GLMM()`/`dGamma_list()`. The exact formulas below are now
wired into the eigenvalue/`lambda_star_marginal` machinery at three call
sites: the internal helper `.two_block_truncated_omega_moments()`
(`R/two_block_ergodicity_ing_marginal.R`), which
`.two_block_lambda_star_marginal_over_draws()` (the real safeguard inside
`rLMMindepNormalGamma_reg_known_vcov()`, `known_vcov` engine only) now calls
in place of the old fixed-\(\Omega_j^{\mathrm{eff}}\) formula; its scratch
mirror `.tmp_lambda_star_marginal_over_draws()`
(`data-raw/_scratch_lambda_star_marginal_over_draws.R`); and the per-draw
ellipsoid test `.tmp_rss_ellipsoid_test()`
(`data-raw/_scratch_rss_ellipsoid_test.R`), which now reports both the
UNTRUNCATED (legacy, `pct_draws_outside`) and EXACT (`pct_draws_outside_exact`)
violation rates side by side. The untruncated formulas are kept everywhere as
a conservative comparison column, per the "Practical implication" paragraph
below.

**Setup.** Since \(\Omega_j=1/\sigma_j^2\), the truncation on precision is the
window with endpoints swapped:

$$
\omega_L:=1/\texttt{disp\_upper}_j,\qquad\omega_U:=1/\texttt{disp\_lower}_j.
$$

Reuse \(s_j:=a_j^0+n_j/2\) (shape, fixed) and \(t_j(\beta_j):=r_j^0+\tfrac12
e_j'e_j\) (rate, the only \(\beta_j\)-dependence).

**Exact marginal.** The joint density in \(\Omega_j\) given \(\beta_j\),
restricted to \([\omega_L,\omega_U]\), is \(\omega^{s_j-1}
e^{-t_j(\beta_j)\omega}\). Using the incomplete-gamma identity
\(\int_0^x u^{k-1}e^{-tu}\,du=\Gamma(k)t^{-k}\,\texttt{pgamma}(x;k,t)\)
(`pgamma(x, shape = k, rate = t)`, R's parameterization),

$$
p(y_j\mid\beta_j)\;\propto\;I_j(t_j)\;:=\;\Gamma(s_j)\,t_j^{-s_j}\,
\Big[\texttt{pgamma}(\omega_U;s_j,t_j)-\texttt{pgamma}(\omega_L;s_j,t_j)\Big],
$$

i.e. the familiar Student-t kernel \(t_j^{-s_j}\) from §16.1, divided by a
*truncation-fraction correction* \(\Delta P(s_j,t_j):=\texttt{pgamma}
(\omega_U;s_j,t_j)-\texttt{pgamma}(\omega_L;s_j,t_j)\) — the probability the
untruncated \(\mathrm{Gamma}(s_j,t_j)\) puts on \([\omega_L,\omega_U]\) — that
carries all the truncation information and depends on \(\beta_j\) (through
\(t_j\)) exactly like the kernel itself does.

**Exact Hessian.** The same exponential-tilting argument that reproduces
§16.2's formula in the untruncated case gives, for *any* mixing measure on
\(\Omega_j\) (truncated or not),

$$
H_j(\beta_j)\;=\;E_t[\Omega_j]\,D_j'D_j\;-\;\mathrm{Var}_t[\Omega_j]\,
(D_j'e_j)(D_j'e_j)',
$$

where \(E_t[\cdot],\mathrm{Var}_t[\cdot]\) are the mean/variance of
\(\Omega_j\) under the tilted-and-truncated density \(\propto
\omega^{s_j-1}e^{-t_j\omega}\) on \([\omega_L,\omega_U]\) — i.e. the truncated
\(\mathrm{Gamma}(s_j,t_j(\beta_j))\) posterior for \(\Omega_j\) given
\(\beta_j\), the one `dGamma()`'s Block~1 step actually draws from. The
standard "raise the shape by 1 (or 2)" identity for incomplete-gamma moments
gives these in closed form:

$$
E_t[\Omega_j]=\frac{s_j}{t_j}\cdot\frac{\Delta P(s_j+1,t_j)}{\Delta P(s_j,t_j)},
\qquad
E_t[\Omega_j^2]=\frac{s_j(s_j+1)}{t_j^2}\cdot\frac{\Delta P(s_j+2,t_j)}{\Delta P(s_j,t_j)},
\qquad
\mathrm{Var}_t[\Omega_j]=E_t[\Omega_j^2]-\big(E_t[\Omega_j]\big)^2.
$$

Three `pgamma()` calls (shapes \(s_j\), \(s_j+1\), \(s_j+2\), all at rate
\(t_j(\beta_j)\)) give the exact truncated Hessian — no numerical integration,
no approximation.

**Consistency check.** As \(\texttt{disp\_lower}_j\to0\) and
\(\texttt{disp\_upper}_j\to\infty\) (\(\omega_L\to0\), \(\omega_U\to\infty\)),
\(\texttt{pgamma}(\omega_U;k,t)\to1\) and \(\texttt{pgamma}(\omega_L;k,t)\to0\)
for every \(k\), so \(\Delta P(k,t)\to1\) regardless of \(k\), and the formulas
collapse to \(E_t[\Omega_j]\to s_j/t_j\), \(\mathrm{Var}_t[\Omega_j]\to
s_j/t_j^2\) — exactly §16.2's untruncated moments. This is a strict
generalization, not a different model.

**The degenerate limit: recovering the fixed-\(\Omega_j\) plug-in exactly.**
At the opposite extreme, as the window collapses to a point
(\(\texttt{disp\_lower}_j\uparrow\texttt{disp\_upper}_j\to\omega_0^{-1}\), i.e.
\(\omega_L,\omega_U\to\omega_0\)), \(\Omega_j\) is forced to equal \(\omega_0\)
with certainty regardless of \(t_j(\beta_j)\): \(E_t[\Omega_j]\to\omega_0\)
and \(\mathrm{Var}_t[\Omega_j]\to0\). The rank-one correction term vanishes
entirely and

$$
H_j(\beta_j)\;\to\;\omega_0\,D_j'D_j,
$$

exactly the ordinary *fixed*-\(\Omega_j\) Gaussian Hessian §14's joint
extension (and the plain, non-ING fixed-vcov engine) already uses — globally
PD whenever \(D_j\) has full column rank, with no finite-radius ellipsoid at
all. So §16.6's family of Hessians interpolates continuously between two
known, previously-derived extremes as the truncation window narrows or widens:
the plug-in fixed-precision case at one end, §16.2's full (untruncated)
Student-t at the other.

**Which moment does the window width actually move?** \(E_t[\Omega_j]\) and
\(\mathrm{Var}_t[\Omega_j]\) do not respond to the window the same way:

- \(E_t[\Omega_j]\) is a *mean*. For the quantile-based windows
  `Prior_Setup_GLMM()` actually constructs (symmetric-ish quantiles of the
  same untruncated \(\mathrm{Gamma}(s_j,t_j)\) around its own center),
  widening or narrowing the window symmetrically moves probability mass off
  *both* tails at once, and the two changes largely cancel in their effect on
  the mean — so \(E_t[\Omega_j]\) should move comparatively little as
  \(\texttt{max\_disp\_perc}\)/\(\texttt{pwt\_measurement}\) are varied,
  provided the window stays reasonably centered on the untruncated
  distribution's own mass.
- \(\mathrm{Var}_t[\Omega_j]\), by contrast, is exactly what truncation acts
  on directly: for nested truncation intervals sharing (approximately) the
  same center, a *narrower* window can only remove probability mass further
  from that center, which — for the unimodal Gamma densities in play here —
  can only shrink the variance monotonically; a *wider* window restores mass
  further out and grows it back, up to the untruncated \(s_j/t_j^2\) ceiling.
  (This monotonicity is the standard "truncating a unimodal density to a
  smaller interval around its mode cannot increase its variance" fact; it is
  stated here as a working remark, not reproven in full generality for every
  possible pair of nested intervals.)

Combined with the §16.3-style criterion \(q_j\le E_t[\Omega_j]/
\mathrm{Var}_t[\Omega_j]\) (dividing the PSD condition
\(\mathrm{Var}_t[\Omega_j]\,q_j\le E_t[\Omega_j]\) through), this gives a
direct, Hessian-level explanation for the tail-probability-level finding
already in `inst/omega-ing-marginal-multivariate-t.md` (a sharper prior
shrinking the expected fraction of draws outside the ellipsoid): tightening
`disp_lower`/`disp_upper` mainly *lowers* \(\mathrm{Var}_t[\Omega_j]\) (with
\(E_t[\Omega_j]\) roughly unchanged), which directly *raises* the effective
threshold \(E_t[\Omega_j]/\mathrm{Var}_t[\Omega_j]\) on the right-hand side —
making the log-concavity violation harder to trigger — while widening the
window pushes the threshold back down toward §16.3's untruncated value. The
two limits above are exactly the two ends of that same continuum: a point
window (\(\mathrm{Var}_t\to0\)) sends the threshold to \(+\infty\) (no finite
\(\beta_j\) can violate it), and the fully untruncated window
(\(\mathrm{Var}_t\to s_j/t_j^2\)) recovers §16.3's finite-radius ellipsoid
exactly.

**Practical implication.** The `lambda_star_marginal`/ellipsoid diagnostics
(`R/two_block_ergodicity_ing_marginal.R`,
`data-raw/_scratch_lambda_star_marginal_over_draws.R`,
`data-raw/_scratch_rss_ellipsoid_test.R`) now compute the *exact*, truncation-
aware §16.6 Hessian/threshold by default (three extra `pgamma()` calls per
group per draw via `.two_block_truncated_omega_moments()`), and additionally
report the old *untruncated* §16.2/16.3 numbers alongside for comparison
(`h_violates_count` vs. `h_violates_count_exact` in the `lambda_star_marginal`
diagnostics; `pct_draws_outside` vs. `pct_draws_outside_exact` in the
ellipsoid test). Per the consistency check above, the untruncated numbers are
systematically *conservative* relative to the exact ones (they overstate how
far into the tail \(\mathrm{Var}_t[\Omega_j]\) actually reaches, hence
understate the true threshold and overstate the true violation probability).
Empirically (see `demo("Ex_13b_...")`/`demo("Ex_13c_...")`'s Section 7d/2b),
the effect is often large: under the current `Prior_Setup_GLMM()` default
window (\(\texttt{max\_disp\_perc}=0.8\), i.e. the central 60% prior-mass
window), even previously-flagged outlier groups can show an exact violation
rate of essentially \(0\%\) although their *untruncated* rate was several
percent — i.e. once the truncation is accounted for exactly, the per-group
`pwt_measurement` tailoring of §9 (`inst/omega-ing-marginal-multivariate-t.md`)
may no longer be necessary at the package's default truncation window; see
`data-raw/_scratch_group_pwt_measurement_noncentral.R`'s `w_star_exact` column
for the corresponding per-group re-derivation.

### 16.7 The bad region is a single (untruncated) or bounded double (truncated) crossing — never more

**Status: theory-only, no code changes are made here.** Fix any single
direction \(u\) and parameterize \(\beta_j=\hat\beta_j^{\mathrm{ols}}+su\) by
the signed distance \(s\). Because \(D_j'\hat e_j^{\mathrm{ols}}=0\), both
\(q_j(s)\) and \(t_j(s)=r_j^0+\tfrac12e_j'e_j\) depend on \(s\) only through
\(s^2\) (Pythagoras, §16.3), so the whole problem along this ray is an *even*
function of \(s\) — it suffices to track \(s\ge0\).

**Untruncated case: exactly one crossing, then convex forever, decaying to
zero.** §16.2's criterion is \(q_j(s)\le E[\Omega_j\mid\beta_j]/
\mathrm{Var}[\Omega_j\mid\beta_j]=t_j(s)/(a_j^0+n_j/2)\) — the right side is
*linear* in \(t_j(s)\), and \(t_j(s)\), \(q_j(s)\) are both \(\Theta(s^2)\) for
large \(s\), so the gap between the two sides, once it opens (at
\(s=\sqrt{2r_j^0+\mathrm{RSS}_j^{\mathrm{ols}}}\) per §16.3), never closes
again — there is exactly one crossing per side. The magnitude of the
resulting negative eigenvalue does not diverge, though: it is
\(\propto\mathrm{Var}[\Omega_j\mid\beta_j]\|D_j'e_j\|^2\propto
s^2/t_j(s)^2\propto1/s^2\to0\). (This is the exact multivariate analogue of
the univariate Student-t fact that \(f''(x)=-(\nu+1)(\nu-x^2)/(\nu+x^2)^2\)
crosses zero once at \(x=\pm\sqrt\nu\) and then decays like \((\nu+1)/x^2\)
without ever re-crossing; \(f'(x)=-(\nu+1)x/(\nu+x^2)\) stays strictly
negative throughout, i.e. the log-density never "turns upward" — it only
flattens, consistent with the Student-t's polynomial, strictly monotone
tail.)

**Truncated case: a second crossing restores PD far out.** With a *fixed*,
nontrivial window \([\omega_L,\omega_U]\), \(\omega_L>0\), §16.6's threshold
is \(E_t[\Omega_j]/\mathrm{Var}_t[\Omega_j]\), and this no longer stays linear
in \(t_j\) as \(t_j\to\infty\). As \(t_j\to\infty\) the tilted density
\(\propto\omega^{s_j-1}e^{-t_j\omega}\) concentrates at the left edge
\(\omega_L\); substituting \(\omega=\omega_L+v/t_j\), the factor
\((\omega_L+v/t_j)^{s_j-1}\to\omega_L^{s_j-1}\) pointwise (for *any* fixed
\(s_j\), this needs no large-shape approximation) while \(e^{-t_j\omega}=
e^{-t_j\omega_L}e^{-v}\) sharpens into a plain \(\mathrm{Exponential}(1)\)
shape in \(v\) — a standard boundary-Laplace argument. To leading order in
\(1/t_j\),

$$
E_t[\Omega_j]\to\omega_L+\frac1{t_j},\qquad
\mathrm{Var}_t[\Omega_j]\to\frac1{t_j^2},\qquad\text{so}\qquad
\frac{E_t[\Omega_j]}{\mathrm{Var}_t[\Omega_j]}\;\sim\;\omega_L\,t_j^2.
$$

The threshold now grows **quadratically in \(t_j\)** (quartically in \(s\)),
while \(q_j(s)\) still only grows like \(t_j\) (quadratically in \(s\)). A
quartic eventually overtakes a quadratic, so the criterion \(q_j\le
E_t[\Omega_j]/\mathrm{Var}_t[\Omega_j]\) — violated on some intermediate
range — is satisfied again for \(s\) large enough. Along each ray, the sign
pattern of \(H_j\)'s bad eigendirection is therefore

$$
\underbrace{\text{PSD}}_{|s|<x_1}\ \Big|\ \underbrace{\text{not PSD}}_{x_1\le|s|\le x_2}\ \Big|\ \underbrace{\text{PSD}}_{|s|>x_2},
$$

a **bounded** annulus \([x_1,x_2]\) per side, not a half-line — the direct
consequence of \(\omega_L>0\) strictly (an untruncated or zero-lower-bound
window degenerates back to the single-crossing case above, since then the
threshold's growth rate in \(t_j\) drops from quadratic back to linear).

### 16.8 A chord majorizes the bounded annulus, exactly licensed by Nygren and Nygren (2006)

**Status: theory-only.** Because \(f\) (equivalently \(\log p(y_j\mid\beta_j)\)
as a function of \(s\) along the fixed ray of §16.7) is convex on the bounded
interval \([x_1,x_2]\), the chord joining \((x_1,f(x_1))\) and
\((x_2,f(x_2))\) is a valid upper bound there:

$$
f(s)\;\le\;\ell(s):=f(x_1)+\frac{f(x_2)-f(x_1)}{x_2-x_1}(s-x_1),\qquad s\in[x_1,x_2],
$$

and since \(f\) is decreasing across the whole annulus (§16.7), \(\ell\) has
negative slope, so \(p(y_j\mid\beta_j)\propto e^{f(s)}\le e^{\ell(s)}\) decays
exponentially over the annulus with an explicit rate — no numerical
integration needed, since \(f(x_1)\), \(f(x_2)\) come straight from §16.6's
closed-form \(I_j(t_j)\).

This device is not an ad hoc patch on top of Nygren and Nygren (2006,
*Likelihood Subgradient Densities*, JASA 101(475):1144-1156, the paper behind
`glmbayes`'s `EnvelopeBuild_c()`/`Gridtype` machinery documented in
`vignettes/Chapter-A05.Rmd`) — it is a direct instance of that paper's own
most general construction. Definition 2 there only requires a bounding
function \(g\ge f\) for which \(-\log g\) has a subgradient (i.e. is convex)
at the chosen point; \(g\) need not equal \(f\). Fact A.1/Remark A.2 make the
existence criterion explicit: a valid (generalized) likelihood-subgradient
density exists on a set if and only if \(f\) is bounded above there by a
function of the form \(\exp(a+b^Ts)\) — exactly the form \(e^{\ell(s)}\)
takes. Using \(g=f\) directly (as the paper's plain "likelihood-subgradient
density" special case does, and as §16's tangent pieces below implicitly do)
fails specifically on \([x_1,x_2]\), because \(-\log f\) is not convex there
by construction — but \(g=e^{\ell(s)}\) succeeds, satisfying Fact A.1 exactly,
with \(-\log g=-\ell(s)\) already affine (trivially convex, with constant
"subgradient" equal to \(-\ell\)'s own slope everywhere on the piece).

Two further points push this well past mere formal admissibility:

- **Claim 1's algebra does not care where \(c\) came from.** For a
  multivariate normal prior, any fixed subgradient vector \(c(\bar\theta)\)
  produces a restricted density that is again multivariate normal with mean
  \(\tilde\mu=\mu-\Sigma c(\bar\theta)\), regardless of whether \(c\) is a
  genuine gradient \(-\nabla\log f(\bar\theta)\) (a tangent, valid on the
  concave pieces) or the chord's constant slope (valid, via Fact A.1, on the
  bounded convex annulus). So Remark 5's restriction/renormalization,
  Remark 6's explicit mixture-weight formula, and Example 2's inverse-
  transform sampling of a box-truncated normal (Claim 2's whole apparatus)
  all apply *unchanged* to a chord-anchored piece — no new sampling machinery
  is required, only a different (constant, not state-dependent) \(c\).
- **The resulting partition is a direct five-piece generalization of §3.2's
  three-piece construction**, alternating tangent/chord/tangent/chord/tangent
  along each ray: two far-tail tangent pieces (PD restored, §16.7), two
  bounded chord pieces (the annuli), and one center tangent piece — mirroring
  \(\{\theta^\ast-\omega,\theta^\ast,\theta^\ast+\omega\}\)'s three anchors
  with the five anchors \(\{-x_2,-x_1,\theta^\ast,x_1,x_2\}\).

**What does *not* transfer for free.** Theorem 2/3's specific asymptotic
result (\(\tilde a\to2/\sqrt\pi\) as \(N\to\infty\)) is proved (Appendix,
Claims A.1-A.4) via a Mills'-ratio argument tied to the tangent-only,
Gaussian-data construction; the analogous tightness limit for a chord piece
(as \(n_j\to\infty\), or as \([\omega_L,\omega_U]\) is varied) would need its
own derivation from §16.6's exact marginal, not inherited automatically.
Likewise, §3.3's multivariate "standard form" diagonalizes the Hessian
*once*, at the mode — a fixed matrix — whereas here the bad eigendirection
\(v(\beta_j)=D_j'D_j(\hat\beta_j^{\mathrm{ols}}-\beta_j)\) is state-dependent
(rotates with \(\beta_j\) in general, though it stays fixed in direction along
any single ray from \(\hat\beta_j^{\mathrm{ols}}\), per §16.7's setup); a
literal multivariate five-piece grid would need either a direction-uniform
worst-case pair of annulus radii or a genuinely per-direction construction,
neither of which is derived here. This section remains, like §16.1-16.7, a
diagnostic/theoretical note only — nothing here is implemented in
`R/two_block_ergodicity_ing_marginal.R` or elsewhere.

### 16.9 Combining the marginal safeguard with the disp_upper plug-in envelope, rather than replacing it

**Status: implemented** (`.two_block_combine_rate_envelopes()`,
`R/two_block_ergodicity_ing_marginal.R`; called from
`.rLMMIngNormal_reg_run_with_pilot()` in `R/rLMM_reg.R`).

Earlier versions of the pilot safeguard treated
\(\lambda^\star_{\mathrm{marginal}}\) (§16, computed at each pilot draw's own
\(\beta_j\), maxed over draws) and \(\lambda^\star_{\mathrm{upper}}\) (§14,
the fixed \(\omega_L=1/\texttt{disp\_upper}_j\) plug-in envelope, D0 = 0) as
mutually exclusive: if the marginal envelope was computable and valid
(`marginal_rate_valid == TRUE`), it *replaced* the plug-in one outright.
Neither envelope is an unconditionally valid upper bound on its own,
however, each has a distinct blind spot:

- \(\lambda^\star_{\mathrm{upper}}\) fixes the per-group precision at
  \(\omega_L D_j'D_j\) — the degenerate, zero-variance-at-\(\omega_L\) corner
  of the exact §16.6 family (§16.6's "degenerate limit" specialized to
  \(\omega_0=\omega_L\)) — regardless of what the pilot draws' own residuals
  actually look like. It is only guaranteed to lower-bound the true
  \(H_j(\beta_j)\) (and hence upper-bound \(\lambda^\star\)) on the *tighter*
  sub-ellipsoid \(q_j\le K_j':=(E_t[\Omega_j]-\omega_L)/\mathrm{Var}_t[\Omega_j]\),
  strictly smaller than the ordinary log-concavity ellipsoid
  \(q_j\le K_j:=E_t[\Omega_j]/\mathrm{Var}_t[\Omega_j]\) that §16.3/§16.6
  already check — so a draw can sit inside the checked ellipsoid (no PD
  failure) yet still, in principle, have a true local rate above what
  \(\lambda^\star_{\mathrm{upper}}\) certifies.
- \(\lambda^\star_{\mathrm{marginal}}\) uses the exact, state-dependent
  \(H_j(\beta_j)\) at every observed pilot draw, but is only a \(\max\) over
  the *finite* pilot sample (\(n_{\mathrm{pilot}}\) draws); a worse state the
  pilot happened not to visit is not certified against.

Rather than picking one, `.two_block_combine_rate_envelopes()` builds a
single conservative envelope: sort each spectrum ascending (the order
`two_block_tv_bound()`/Lemma 2 already require internally), take the
componentwise (rank-matched) maximum, and run Theorem 3's TV bound once on
the combined spectrum at the pilot-mean-start \(D_0\)
(`two_block_d0_pilot_start()`) both envelopes are already certified against.
Because the componentwise max of two rank-\(i\) non-decreasing sequences is
itself non-decreasing, the combined object's own \(\lambda^\star\) is simply
its largest (last) entry, and the resulting \(m_{\min,\mathrm{combined}}\)
is \(\ge\) both \(m_{\min,\mathrm{upper}}\) and \(m_{\min,\mathrm{marginal}}\)
individually — in fact typically strictly larger than either, since Theorem
3's bound sums over the *whole* spectrum, not just the top eigenvalue, so a
rank-\(i\) crossover elsewhere in the spectrum (one envelope's rank-\(i\)
eigenvalue exceeding the other's, even where the *top* eigenvalues favor the
other envelope) still adds to the combined bound. `m_convergence_used` is
then the max of this combined \(m_{\min}\) and the plain-rate floor computed
at pilot-start (mirroring the existing "only ever increase" pattern already
used to fold `pilot_ub$m_min_upper` into `m_convergence_used`).
`fit$convergence_info$lambda_star_combined`/`eigenvalues_combined`/
`m_min_combined` record the combined spectrum's own summary; when the
marginal safeguard is invalid, behavior is unchanged (falls back to the
plain-rate envelope alone, with the same warning as before).

---

## 17. A split-support heuristic for a revised end-of-simulation TV bound

**Status: implemented (diagnostic, `known_vcov` engine only).** §§14-16 all
make the same point from different angles: once \(\Omega_j\) is state
dependent, the contraction matrix \(A\) that Theorem 3 needs is not a
constant, so a single \(\lambda^*\) calibrated once (at the pilot mean, or as
a pilot-draws \(\max\) via the marginal safeguard of §16/`R/
two_block_ergodicity_ing_marginal.R`) is only a *local* rate. This section
gives a post-hoc, end-of-run correction that turns the empirical per-draw scan
already used in §7d of `demo("Ex_13b_...")`/`demo("Ex_13c_...")` into a
revised total-variation estimate, rather than leaving the calibrated
\(\texttt{tv\_tol}\) (e.g. 0.01) as the reported number unconditionally.

### 17.1 Why the calibrated bound can be optimistic here

`two_block_l_for_tv()` picks `m_convergence` so that, **if** the chain were
the fixed-\(\Omega\) Gaussian two-block Gibbs sampler with contraction matrix
\(A(\Lambda,\Omega^\star)\) (spectral radius \(\lambda^\star\)) *everywhere*,
Theorem 3 guarantees
\(d_{\mathrm{TV}}(P^m(x_0,\cdot),\pi)\le\texttt{tv\_tol}\). In the ING
`known_vcov` engine, \(\Omega_j\) is instead resampled every sweep, so the
*local* rate at the state actually visited varies draw to draw (§14), and can
in principle exceed the \(\lambda^\star\) used to size \(m\) — or the local
quadratic approximation can fail to be log-concave at all (§16.3's \(H_j\)
PSD criterion; equivalently \(\Lambda+H_j(\beta_j)\) failing to be PD, the
inversion the marginal scan needs). Theorem 3's proof gives no guarantee for
sweeps that land in such a state. The calibrated `tv_tol` is therefore only
demonstrably valid *conditional on staying in the region where the local rate
behaves*; this section estimates how much slack to add for the rest.

### 17.2 Splitting the support

Let \(\lambda^\star_{\mathrm{used}}\) be whichever rate actually calibrated
`m_convergence` for the run being diagnosed —
\(\texttt{fit\$convergence\_info\$lambda\_star\_combined}\) if the marginal
safeguard was valid for that pilot (`marginal_rate_valid == TRUE`; §16.9
below), else the base/extended \(\lambda^\star_{\text{upper}}\) it fell back
to (§14, `two_block_rate_ing()`/`two_block_rate()`). Partition the state
space into

$$
A \;=\; \Big\{\,\beta:\ \Lambda+H_j(\beta_j)\succ0\ \ \forall j,\ \ \text{and}\ \ \lambda^\star_{\mathrm{marginal}}(\beta)\le\lambda^\star_{\mathrm{used}}\,\Big\},
\qquad A^c \;=\; \text{its complement,}
$$

using the §16.2 per-draw marginal Hessian \(H_j(\beta_j)\) — exactly the
quantity `.two_block_lambda_star_marginal_over_draws()`
(`R/two_block_ergodicity_ing_marginal.R`, mirrored for demo use by
`.tmp_lambda_star_marginal_over_draws()` in
`data-raw/_scratch_lambda_star_marginal_over_draws.R`) already computes for
every main-stage draw in §7d. A draw is in \(A^c\) iff *either* some group's
block failed the PD check (`skipped[i]`) *or* its computed
`lambda_star_vec[i] > lambda_star_used`.

### 17.3 A triangle-inequality decomposition

For any probability measures \(\mu,\nu\) and measurable \(A\), and any
\(B\subseteq\mathcal X\),

$$
|\mu(B)-\nu(B)| \;=\; \big|\mu(B\cap A)+\mu(B\cap A^c)-\nu(B\cap A)-\nu(B\cap A^c)\big|
\;\le\; \big|\mu(B\cap A)-\nu(B\cap A)\big| \;+\; \mu(A^c) \;+\; \nu(A^c).
$$

Taking the supremum over \(B\) (the definition of \(d_{\mathrm{TV}}\)) with
\(\mu=P^m(x_0,\cdot)\), \(\nu=\pi\):

$$
d_{\mathrm{TV}}\big(P^m(x_0,\cdot),\pi\big) \;\le\;
\underbrace{\sup_B\big|P^m(x_0,B\cap A)-\pi(B\cap A)\big|}_{\text{"on-}A\text{" gap}}
\;+\; P^m(x_0,A^c) \;+\; \pi(A^c).
$$

**Heuristic assumption.** The "on-\(A\)" gap is bounded by `tv_tol`: i.e. that
restricted to trajectories that stay inside the region where the local rate
is \(\le\lambda^\star_{\mathrm{used}}\) and every block is PD, the chain
behaves enough like the homogeneous-rate chain Theorem 3 was calibrated
against that its guarantee still applies. This is *not* a proven fact for the
state-dependent system (that would require its own homogeneous-in-\(A\)
contraction argument); it is the natural reading of what `tv_tol` was
supposed to buy, restricted to where the argument is on firm ground.

### 17.4 Estimating the \(A^c\) terms and the resulting formula

\(P^m(x_0,A^c)\) and \(\pi(A^c)\) are two conceptually different quantities —
the chain's actual excursion probability starting from \(x_0\), and the
*stationary* excursion probability — that this diagnostic cannot tell apart
from a single stream of stored draws. Since the whole point of calibrating
`m_convergence` is that the stored main-stage draws are already close to
\(\pi\), both are approximated by the same empirical fraction

$$
\hat p \;=\; \frac{1}{n}\sum_{i=1}^n \mathbf{1}\{\text{draw }i\in A^c\},
$$

over the \(n\) main-stage draws, giving

$$
\mathrm{TV}_{\mathrm{revised}} \;=\; \texttt{tv\_tol} + 2\hat p.
$$

A looser, one-sided variant \(\mathrm{TV}_{\mathrm{revised,\,loose}} =
\texttt{tv\_tol}+\hat p\) — dropping either the \(P^m(x_0,A^c)\) or the
\(\pi(A^c)\) term — is also reported as a sensitivity check; it is not the
quantity the triangle-inequality argument above actually derives, but bounds
how much of \(\mathrm{TV}_{\mathrm{revised}}\) comes from doubling a single
empirical estimate versus from the \(\hat p\) itself.

### 17.5 Caveats

- This is a **diagnostic augmentation of the reported number, not a new
  certified theorem** — precisely because §17.3's "on-\(A\)" step is a
  heuristic, not something proven here or elsewhere in this document. It
  should be read the same way as the existing `known_vcov` marginal
  safeguard: informative, not a substitute for keeping \(\hat p\) small by
  construction (e.g. removing outlier groups, as in `demo("Ex_13c_...")`).
- \(\hat p\) is computed from the **main-stage** draws — the same draws whose
  distribution `TV_revised` is trying to characterize — so this is
  intrinsically self-referential/post-hoc, not a pre-run guarantee the way
  `two_block_l_for_tv()`'s own `tv_tol` calibration is. A large \(\hat p\)
  is itself informative (it says the run spent a non-trivial fraction of its
  time outside the region the calibration trusted) independent of whether
  \(2\hat p\) is the "right" multiplier.
- Scope matches the marginal safeguard of §16/`R/
  two_block_ergodicity_ing_marginal.R`: `known_vcov` engine only, where
  \(\Lambda\) (RE precision) is fixed/known and only \(\Omega_j\) is
  state-dependent.
- See `data-raw/_scratch_tv_bound_revised.R` for the scratch implementation
  (temporary, investigation only) and `demo("Ex_13b_...")`/
  `demo("Ex_13c_...")` §7e for its use immediately after §7d's per-draw
  marginal scan, reusing that scan's output rather than recomputing it.

---

## References

- Anderson, W.N. and Duffin, R.J. (1969). "Series and Parallel Addition of
  Matrices." *Journal of Mathematical Analysis and Applications*, 26(3),
  576-594. (operator monotonicity of the matrix parallel sum, used in §10.4.)
- Anderson, B.D.O. and Moore, J.B. (1979). *Optimal Filtering*. Prentice-Hall.
  (monotonicity of the Riccati/information-filter precision update in the
  noise covariances, used in §10.4 for \(p_{re}>1\).)
- Lindley, D.V. and Smith, A.F.M. (1972). "Bayes Estimates for the Linear
  Model." *Journal of the Royal Statistical Society, Series B*, 34(1), 1-41.
  (source of the \(\beta_j,\gamma,\mathcal{W}_j,\Psi\) parameterization; see
  `inst/notation.md`.)
- Nygren, K. (2020). *On the total variation distance between multivariate
  normal densities with applications to two-block Gibbs samplers.*
  Unpublished manuscript. (source of Claims 1-3, Remark 8, Theorem 3,
  Corollary 1; see `R/two_block_ergodicity.R`.)
- Nygren, K.N. and Nygren, L.M. (2006). "Likelihood Subgradient Densities."
  *Journal of the American Statistical Association*, 101(475), 1144-1156.
  (source of Definition 1/2, Fact A.1, Claim 1/2, Theorem 1/2/3 cited in
  §16.8's chord-majorization argument; implemented for `glmbayes`'s
  log-concave GLM samplers in `EnvelopeBuild_c()`, documented in
  `vignettes/Chapter-A05.Rmd`.)
- See also `inst/BLOCK_GIBBS_ERGODICITY.md` (rank-deficiency / identifiability
  theory this note's §8 data-free-group remark connects back to) and
  `inst/notation.md` (the \(D,\mathcal{W},\beta,\gamma,\Psi\) notation used
  throughout).
