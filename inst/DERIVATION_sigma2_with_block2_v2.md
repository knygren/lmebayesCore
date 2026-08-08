# v2 Gibbs partition: observation dispersion with population block

**Status.** Theory + API stubs only. Production sampling for
`rLMMindepNormalGamma_reg_*_v2()` is **not implemented** yet.

**Related code / docs**

| Item | Role today (v1) | Role in v2 |
|---|---|---|
| `rLMMindepNormalGamma_reg_{known,estimated}_vcov` | Joint ING Block~1 \((\beta,\sigma^2)\) | Unchanged |
| `rLMMindepNormalGamma_reg_{known,estimated}_vcov_v2` | — | New partition (stubs) |
| `.lmebayes_ing_measurement_prior_list[_group]` | Builds ING Block~1 `prior_list` (`mu`/`Sigma`/`shape`/`rate`/bounds) | Still used by v1 |
| `.lmebayes_ing_measurement_prior_list[_group]_v2` | — | Measurement Gamma hyperparameters only (no `mu`/`Sigma`) |
| `inst/BLOCK_GIBBS_ERGODICITY_ING.md` §10, §14, §16 | Corner bound + Ω Hessian extension | Mathematics reused below |
| `.rLMM_measurement_rate_inputs` | \(w_i = 1/\texttt{disp\_upper}\) | Same corner inputs initially; interpretation changes |

---

## 1. Why move \(\sigma^2\) out of the group block

### v1 partition (current)

\[
\begin{aligned}
\text{Block G (group):}&\quad
(\beta,\sigma^2)\mid(\gamma,\Psi)
&&\text{(joint ING envelope / }rindepNormalGamma\text{)}\\
\text{Block P (population):}&\quad
(\gamma,\Psi)\mid\beta.
\end{aligned}
\]

Group-level accept–reject for the joint \((\beta_j,\sigma_j^2)\) step often has
very low acceptance when the truncation window
\([\texttt{disp\_lower},\texttt{disp\_upper}]\) is wide or \(n_j\) is large.

### v2 partition (proposed)

\[
\begin{aligned}
\text{Block G:}&
\quad\beta\mid(\gamma,\Psi,\Omega)
&&\text{(Gaussian / fixed-dispersion full conditional)}\\
\text{Block P:}&
\quad(\gamma,\Psi,\Omega)\mid\beta
&&\text{(population + observation precision)}.
\end{aligned}
\]

Here \(\Omega_j = \sigma_j^{-2} I_{n_j}\) (pooled: one shared \(\Omega\);
per-group: one \(\Omega_j\) per \(j\)). Given current \(\Omega\), Block G is
exactly the existing fixed-dispersion LMM Block~1 update
(`rNormal` / `rGLMM_sweep` without ING measurement). Observation dispersion
is drawn in Block P from its Gamma full conditional given residuals
\(e_j = y_j - D_j\beta_j\) (pooled or per-group), jointly or sequentially
with \((\gamma,\tau^2)\).

This matches the **Hessian block naming** already used in
`BLOCK_GIBBS_ERGODICITY_ING.md` §14:

\[
x_1 = (\gamma,\ \lambda_{\mathrm{ING}},\ \Omega),\qquad
x_2 = \beta,
\]

whereas v1 *samples* \(\Omega\) inside \(x_2\)'s joint ING step and only uses
the §10 corner bound as if \(\Omega\) were fixed.

---

## 2. Model and full conditionals (Gaussian LMM)

Likelihood and hierarchical prior (package notation):

\[
y_j\mid\beta_j,\Omega_j\sim N(D_j\beta_j,\Omega_j^{-1}),\qquad
\beta_j\mid\gamma,\Psi\sim N(\mathcal{W}_j\gamma,\Psi).
\]

Gamma prior on precision (shape/rate), pooled or per group:

\[
\Omega\sim\mathrm{Gamma}(a^0,r^0)
\quad\text{or}\quad
\Omega_j\sim\mathrm{Gamma}(a_j^0,r_j^0),
\]

truncated to the same windows that v1 stores as `disp_lower`/`disp_upper`
(on the variance scale \(\sigma^2=1/\Omega\)).

### Block G (v2)

Standard Gaussian conjugate update for \(\beta\) at **fixed** current
\((\gamma,\Psi,\Omega)\). No envelope for \(\sigma^2\).

### Block P (v2) — observation precision

Given \(\beta\),

\[
\Omega_j\mid\beta_j
\sim\mathrm{Gamma}\Bigl(a_j^0+\tfrac{n_j}{2},\;
r_j^0+\tfrac12\|y_j-D_j\beta_j\|^2\Bigr)
\]

(truncated to the prior window). Pooled case uses
\(a^0+\sum_j n_j/2\) and \(\sum_j\|e_j\|^2/2\).

Independently (given \(\beta\)), population updates for
\((\gamma_p,\lambda_p)\) are unchanged from today's Block~2
(Normal or ING on \(\tau^2_p=1/\lambda_p\)).

**Hessian orthogonality (§14.3).** Cross terms
\(P_{11}^{(\gamma\Omega)}=P_{11}^{(\lambda\Omega)}=0\) exactly, so
\((\gamma,\lambda)\) and \(\Omega\) may be drawn in either order inside
Block P without changing the two-block partition \(\{x_1,x_2\}\).

---

## 3. Rate matrix and \(\lambda^*\): what stays, what changes

Recall the fixed-\((\Lambda,\Omega)\) contraction (`two_block_rate()`):

\[
A = P_{11}^{-1/2}\,S\,P_{11}^{-1/2},\qquad
S = P_{12}\,P_{22}^{-1}\,P_{21},\qquad
\lambda^* = \lambda_{\max}(A),
\]

with

\[
P_{22}=\mathrm{blockdiag}_j(B_j),\quad
B_j = D_j'\Omega_j D_j + \Psi^{-1},
\]

and \(P_{11}\) the population precision in \(\gamma\) (and \(\lambda\) when
ING). Production v1 certifies Theorem~3 for the
**\((\gamma,\beta)\mid(\Lambda,\Omega)\) sub-chain** by the §10 corner

\[
\Omega_{\min} = 1/\texttt{disp\_upper}
\quad\text{(least favorable / most diffuse admissible }\Omega\text{)},
\]

plus \(\Lambda^{\max}\) from population `disp_lower` when \(\tau^2\) is ING
(`.rLMM_measurement_rate_inputs`, `two_block_rate_from_pfamily_list`).

### 3.1 Conditional-on-\(\Omega\) sub-chain (still valid inside v2)

Conditional on the current \(\Omega\) (and \(\Lambda\)), Block G is still
the same Gaussian \(\beta\)-update, and the \((\gamma,\beta)\) rate formulas
are unchanged. So:

- The **algebra** of \(B_j\), \(S\), \(A\), \(\lambda^*\) is identical.
- Code path `two_block_rate()` / `.two_block_S_P11()` needs **no formula
  rewrite** for that conditional kernel.

What changes is the **meaning of the certified bound for the full kernel**:
v1's Theorem~3 `m_convergence` treats \(\Omega\) as fixed at the corner;
v2's full kernel *moves* \(\Omega\) in Block P, so that certificate alone is
not a full proof for the joint chain (same caveat already recorded for
population \(\lambda\) in §7–§10).

### 3.2 Extended \(P_{11}\) when \(\Omega\in x_1\) (v2 production target)

Following §14, set \(x_1=(\gamma,\lambda_{\mathrm{ING}},\Omega)\) and
\(x_2=\beta\). Then

\[
P_{11}
=
\begin{bmatrix}
P_{11}^{(\gamma\gamma)} & P_{11}^{(\gamma\lambda)} & 0 \\
P_{11}^{(\gamma\lambda)\top} & P_{11}^{(\lambda\lambda)} & 0 \\
0 & 0 & P_{11}^{(\Omega\Omega)}
\end{bmatrix},
\qquad
P_{12}
=
\begin{bmatrix}
P_{12}^{(\gamma\beta)} \\
P_{12}^{(\lambda\beta)} \\
P_{12}^{(\Omega\beta)}
\end{bmatrix}.
\]

**Unchanged blocks** (same formulas as today / §5–§6):

- \(P_{22}\): \(B_j = D_j'\Omega_j D_j + \mathrm{diag}(\lambda)\).
- \(P_{11}^{(\gamma\gamma)}\), \(P_{11}^{(\gamma\lambda)}\), \(P_{11}^{(\lambda\lambda)}\).
- \(P_{12}^{(\gamma\beta)}\), \(P_{12}^{(\lambda\beta)}\) (state-dependent through
  \(\lambda\) and residuals \(u_{jp}\)).

**New / now-on-the-production-path blocks** (§14):

| Block | Formula (pooled \(\Omega\); per-group analogous) |
|---|---|
| \(P_{11}^{(\Omega\Omega)}\) | \(\dfrac{a^0 + N/2 - 1}{\Omega^2}\) (shape/rate Gamma Hessian) |
| \(P_{12}^{(\Omega\beta)}\) rows | from \(\partial^2\ell/\partial\Omega\,\partial\beta = D'(y-D\beta)\) (residual score) |
| \(P_{11}^{(\gamma\Omega)}\), \(P_{11}^{(\lambda\Omega)}\) | **exactly 0** |

Extended rate:

\[
A_{\mathrm{v2}}(\gamma,\beta,\lambda,\Omega)
=
P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2}.
\]

This is **state-dependent**. Diagnostic code already exists as
`two_block_rate_ing()` / `.two_block_S_P11_ing()` (§15) — not yet used for
production `m_convergence`.

### 3.3 Practical bound strategies for v2 (candidates)

Until a full state-dependent Theorem~3 analogue is certified, v2 can use
(layered, increasingly tight):

1. **Reuse §10 corner on the \((\gamma,\beta)\) sub-chain**  
   Keep `.rLMM_measurement_rate_inputs` → \(w=1/\texttt{disp\_upper}\) and
   (if estimated vcov) RE `disp_lower`.  
   **Interpretation change:** this bounds the \(\beta\leftrightarrow\gamma\)
   coupling at the worst admissible \(\Omega\), but **does not** account for
   extra mixing from moving \(\Omega\) in Block P. Conservative for that
   sub-chain; incomplete for the full kernel (document as such).

2. **Pilot envelope of \(\lambda_{\max}(A_{\mathrm{v2}})\)**  
   Evaluate the §14/§15 extended \(A\) at pilot draws
   \((\hat\gamma,\hat\beta,\hat\lambda,\hat\Omega)\) and take componentwise
   `pmax` spectra (mirror of today's
   `.two_block_pilot_ub_from_coefficients`).  
   Needs production wiring of `.two_block_S_P11_ing` (or equivalent) with
   measurement \(\Omega\) as part of \(x_1\).

3. **§16 marginal-\(\Omega\) Hessian (alternative diagnostic)**  
   Integrate \(\Omega\) out → Student-\(t\) likelihood in \(\beta\); use
   \(H_j(\beta)\) with \(\Omega^{\mathrm{eff}}\) and the rank-one correction.  
   Different model (no sampled \(\Omega\)); useful cross-check, not the v2
   sampler itself.

4. **Known-vcov + fixed population \(\tau^2\)**  
   \(x_1=(\gamma,\Omega)\) only. Extended \(P_{11}\) is block-diagonal
   \(\mathrm{diag}(P_{11}^{(\gamma\gamma)},P_{11}^{(\Omega\Omega)})\); rate
   formulas simplify further (§14.5 stacking).

---

## 4. Eigenvalue monotonicity reminders (unchanged algebra)

For fixed \((\Lambda,\Omega)\), \(\lambda^*\) is Loewner-nonincreasing in
\(\Omega\) and nondecreasing in RE precisions \(\lambda_p\)
(`BLOCK_GIBBS_ERGODICITY_ING.md` §10). Hence:

- Worst measurement corner remains \(\Omega_{\min}=1/\texttt{disp\_upper}\).
- Worst RE corner (ING \(\tau^2\)) remains \(\lambda_p^{\max}\) from
  `disp_lower` on \(\tau^2\).

v2 does **not** change those inequalities for the conditional
\((\gamma,\beta)\) rate; it changes which state variables are free when
forming the **full** \(A_{\mathrm{v2}}\).

---

## 5. Implementation checklist (later)

1. **Block G:** fixed-dispersion sweep (`rGLMM_sweep` / `rNormal` path) with
   current \(\Omega\) from the chain state (not ING joint).
2. **Block P:** existing `.two_block_block2_all_chains` **plus** Gamma draw(s)
   for \(\Omega\) / \(\Omega_j\) given \(\beta\) (pooled vs
   `shape_group`/`rate_group`).
3. **Priors:** `.lmebayes_ing_measurement_prior_list*_v2` — `shape`/`rate`/
   windows only (no `mu`/`Sigma` for joint ING).
4. **Rate:** start with §10 corner + document incompleteness; then wire
   pilot \(A_{\mathrm{v2}}\) via `.two_block_S_P11_ing`.
5. **Returns:** keep `group.dispersion` for observation \(\sigma^2\) draws;
   do not conflate with `popef.dispersion` (\(\tau^2\)).
6. **Routing:** do not replace v1 in `REG_ROUTE_TABLE` until v2 is validated;
   opt-in via `*_v2` exports / future `sim_method`.

---

## 6. Summary

| Question | Answer |
|---|---|
| Do \(B_j\), \(S\), conditional \(\lambda^*\) formulas change? | **No** for fixed \((\Lambda,\Omega)\). |
| Does the production certificate for the full kernel change? | **Yes** — \(\Omega\) joins \(x_1\); need §14 extended \(A\) (or a documented interim corner). |
| What is the main code delta for \(\lambda^*\)? | Treat measurement \(\Omega\) like ING \(\lambda\): append \(P_{11}^{(\Omega\Omega)}\) and \(P_{12}^{(\Omega\beta)}\); keep §10 corner as a **sub-chain** bound only. |
| Why should acceptance improve? | Block G is Gaussian (exact); no joint \((\beta,\sigma^2)\) envelope. |
