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
- See also `inst/BLOCK_GIBBS_ERGODICITY.md` (rank-deficiency / identifiability
  theory this note's §8 data-free-group remark connects back to) and
  `inst/notation.md` (the \(D,\mathcal{W},\beta,\gamma,\Psi\) notation used
  throughout).
