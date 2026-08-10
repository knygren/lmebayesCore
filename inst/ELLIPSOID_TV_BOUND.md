# A Three-Term Total Variation Bound via Safe/Unsafe Region Decomposition

Companion note to `BLOCK_GIBBS_ERGODICITY.md`, `P_MATRIX_MARGINAL_PRECISIONS.md`,
and the Nygren total-variation paper (`Multivariate_Normal_Distances_07_20_20.pdf`).

**Revision note (this version).** Sections 1–4 are now stated at the level
of generality needed to cover *any* block whose conditional negative
log-density has a Hessian that can be bounded via a precision floor and a
full-rank/estimability condition on its design matrix — not just the
truncated-Gamma/multivariate-t case. Sections 5–8 specialize this common
machinery to four likelihood families: multivariate-t (§5), logit (§6),
probit (§7), and Poisson (§8). Each of §5–§8 mirrors the structure of the
old §5 (drift bound on the sampler escape probability), §6 (assembled
bound), and §7 (open items) as subsections X.2–X.4, preceded by a new
subsection X.1 identifying that likelihood's safe region and curvature
floor.

Throughout, `x = (gamma, beta)` denotes the full block-2/block-1 state,
`P^n(x, .)` the `n`-step transition kernel starting from `x`, and `pi(.)`
the target density (for §5, the marginal with precisions integrated out;
for §6–§8, the target conditional on fixed dispersion/precision
hyperparameters, with a fixed Gaussian prior on the regression block).

---

## 1. General setup and notation

For a generic block `beta_j` (dimension `p_j`, design matrix `D_j`, group
size `n_j`), write the block's contribution to the negative log-target as
```
-log pi(beta_j | rest) = -sum_i log L(y_i ; eta_i)  +  0.5 * tau_j^{-2} ||beta_j||^2  + const,   eta_i = D_j[i,] beta_j
```
with local Hessian
```
H_j(beta_j) = tau_j^{-2} I_{p_j}  +  D_j' diag( c(eta_i) ) D_j,     c(eta) := -d^2/deta^2 log L(y;eta)
```
where `c(.)` is the likelihood's curvature function (`= Omega_j` for the
Gaussian/t case with precision `Omega_j`; `= p(1-p)` for logit; the probit
analogue in §7; `= exp(eta)` for Poisson).

**Assumption 1 (full-rank / estimability).** `D_j` has full column rank
`p_j`, and (for the non-Gaussian likelihoods) the group-level MLE
`beta_j_hat` exists and is finite — i.e. the group `j` data do not exhibit
complete or quasi-complete separation (logit/probit) and do not have a
degenerate sufficient statistic (Poisson: not all `y_i=0` in directions
spanned by `D_j`). This is the same condition your code already checks via
a full-rank / estimability screen on `D_j` before attempting to fit a
group-level model, and it is what guarantees the anchor point (`beta_j_hat`
for GLM likelihoods, `beta_j_ols` for the Gaussian/t case) used throughout
below is well-defined.

**Definition 1 (safe region via a precision floor).** Fix a threshold
`kappa_j > 0`. Block `j` is *safe* at `beta_j` if
```
lambda_min( H_j(beta_j) )  >=  kappa_j
```
The **safe region** is `A := intersect_j { beta_j : lambda_min(H_j(beta_j)) >= kappa_j }`
(intersected, where relevant, with analogous per-block conditions for
population/RE blocks `gamma_p`); `A^c` is its complement. This definition
specializes to the ellipsoid condition of `P_MATRIX_MARGINAL_PRECISIONS.md`
§1.3 when `c(eta) = Omega_j` is itself random (t-likelihood case), and to
concrete eigenvalue conditions on `D_j'diag(c(eta_i))D_j` in the GLM cases
below.

---

## 2. The safe/unsafe decomposition

**Lemma 1 (TV splits across a partition).**
*For any measurable partition `Omega = A cup A^c`,*
```
|| P^n(x,.) - pi(.) ||_TV  <=  d_A(x)  +  P^n(x, A^c)  +  pi(A^c)
```
*where `d_A(x) := sup_{B subseteq A} | P^n(x,B) - pi(B) |`.*

*Proof.* For any measurable `B`, write `B = (B cap A) cup (B cap A^c)`, so
```
P^n(x,B) - pi(B) = [P^n(x,B cap A) - pi(B cap A)] + [P^n(x,B cap A^c) - pi(B cap A^c)].
```
Take absolute values and `sup_B`; the first term's supremum over subsets of
`A` is `d_A(x)` by definition. For the second, since `B cap A^c subseteq
A^c`, `P^n(x,B cap A^c) <= P^n(x,A^c)` and `pi(B cap A^c) <= pi(A^c)`, so
its contribution is bounded by `P^n(x,A^c) + pi(A^c)`. **QED**

This decomposition and its proof are entirely likelihood-agnostic — they
only use that `A, A^c` partition the state space.

---

### 2.1 A sharper split: converting the escape-mass term into an estimable quantity

The Corollary to Lemma 6 (§3, derived below) produces a term
`2 P^n(x,A^c) + 2 \pi(A^c)` in the final bound. Rewrite it algebraically
(trivially, for any two numbers):
```
2 P^n(x,A^c) + 2 \pi(A^c)  =  4 P^n(x,A^c) + 2\big[\pi(A^c) - P^n(x,A^c)\big]
```
This is useful only if the bracketed gap can be shown to vanish quickly —
otherwise it's a no-op. It can, but not for free: bounding
`|\pi(A^c)-P^n(x,A^c)|` by the raw TV distance `\|P^n(x,.)-\pi(.)\|_{TV}`
would be circular (that's the quantity being bounded in the first place).
Instead, the bound below uses the *stronger* Meyn–Tweedie `V`-norm (a.k.a.
`f`-norm) geometric ergodicity theorem — which is precisely what
drift-plus-minorization, already assembled in every likelihood section
below, is built to deliver, so no new hypotheses are needed beyond what
each section already establishes.

**Claim 5 (geometric decay of the escape-probability gap).** *Suppose the
Lyapunov function `V` (as in Lemma 4, §5.2/§7.2/§8.2, or the whole-space
argument of §6.2) satisfies a geometric drift condition off a small set `S`
on which minorization holds — i.e. the standard Foster–Lyapunov pair
`(V,\lambda,b,S,\epsilon)` used throughout this note. Then (Meyn & Tweedie
1993, Thm 16.0.1) there exist `R<\infty`, `\rho\in(0,1)` such that*
```
\sup_{|f|\le V} \big| P^n f(x) - \pi f \big|  \le  R \rho^n V(x)
```
*for every `n`. Since `\mathbf 1_{A^c}(y) \le V(y)/K^\star` pointwise (by
the definition of `K^\star` in Lemma 4), taking `f = K^\star \mathbf
1_{A^c}` (a valid test function, `|f|\le V`) gives*
```
\big| \pi(A^c) - P^n(x,A^c) \big|  \le  \frac{R \rho^n V(x)}{K^\star}
```

**Corollary (the refined escape term).**
```
2 P^n(x,A^c) + 2 \pi(A^c)  \;\le\;  4 P^n(x,A^c) \;+\; \frac{2 R \rho^n V(x)}{K^\star}
```

**Practical consequence.** `P^n(x,A^c)` is, by definition, the probability
that the chain — started at `x`, run for `n` sweeps — lands in the unsafe
region. That is exactly what a running MCMC chain lets you estimate
directly: the empirical proportion of post-burn-in draws falling in `A^c`
(with a Monte Carlo standard error, or a concentration-based confidence
band if a rigorous interval is wanted) is a direct estimate of
`P^n(x,A^c)`. So rather than relying solely on the closed-form but
possibly loose analytic bound on `\pi(A^c)` (Claim 2, §4 — built from a
union bound and, for logit (§6), not even in closed form), you can
substitute the empirical estimate for the dominant `4P^n(x,A^c)` term and
rely on analysis only for the correction `2R\rho^n V(x)/K^\star`, which is
provably geometric and — for `n` past a modest burn-in — typically
negligible relative to the empirical term itself.

**Caveat: `(R,\rho)` are rarely explicit outside the Gaussian case.** For
§5 (multivariate-t), the *exact* Nygren machinery (Theorem 3 / Corollary 1)
gives genuinely explicit constants in place of the abstract `(R,\rho)` from
the general Meyn–Tweedie theorem, so Claim 5 can be sharpened there without
appeal to the abstract theorem at all. For §7 (probit) and §8 (Poisson),
`(R,\rho)` exist by the cited theorem but are not computed explicitly in
this note (they depend on the specific drift/minorization constants
`\lambda,b,\epsilon`, none of which have been reduced to closed form yet —
see §9). For §6 (logit), the *global* strong log-concavity argument gives
something strictly better than Claim 5: because the whole space is a
single small set, the chain is *uniformly* (not merely geometrically)
ergodic, so `|\pi(A^c)-P^n(x,A^c)| \le \|P^n(x,.)-\pi(.)\|_{TV} \le M
\rho_{\text{glob}}^n` directly, with `\rho_{\text{glob}}` computable from
the global sandwich `[\tau_j^{-2}, M_j^{\text{theory}}]` — no `V`-norm
detour needed.

---

## 3. Bounding `d_A`: a normalization-corrected version

**3.0 Why naively applying a Gaussian-comparison theorem to `d_A` is
invalid.** `P^n(x,\cdot)` and `pi(\cdot)` restricted to `B subseteq A` are
sub-probability measures integrating to `p := P^n(x,A) <= 1` and
`q := pi(A) <= 1`, generally `p != q`. A theorem stated for two properly
normalized densities (e.g. the Nygren MVN comparison theorem, or any
comparable device for a GLM likelihood) cannot be applied directly to these
unnormalized restrictions.

**Definitions.** For `B subseteq A`, `P^n_A(x,B) := P^n(x,B)/p`,
`pi_A(B) := pi(B)/q` — the conditional (properly normalized) laws on `A`.

**Lemma 5 (three-term decomposition).** *For `B subseteq A`,*
```
P^n(x,B) - pi(B) = [P^n_A(x,B) - pi_A(B)] + [P^n(x,B) - P^n_A(x,B)] - [pi(B) - pi_A(B)]
```
*Proof.* Add and subtract `P^n_A(x,B)` and `pi_A(B)`; identical algebra
regardless of the likelihood generating `P^n`, `pi`. **QED**

**Lemma 6 (bound on `d_A`).**
```
d_A(x)  <=  || P^n_A(x,.) - pi_A(.) ||_TV  +  P^n(x,A^c)  +  pi(A^c)
```
*Proof.* `sup_{B subseteq A}` of Lemma 5, triangle inequality across the
three terms; `sup_B |P^n(x,B)-P^n_A(x,B)| = (1-p)` and
`sup_B|pi(B)-pi_A(B)| = (1-q)` as in the Gaussian-case derivation, since the
algebra `P^n(x,B) = p \cdot P^n_A(x,B)` (and similarly for `pi`) holds for
any likelihood. **QED**

**Corollary (cost of the fix).** Combining Lemma 1 and Lemma 6,
```
|| P^n(x,.) - pi(.) ||_TV  <=  || P^n_A(x,.) - pi_A(.) ||_TV  +  2 P^n(x,A^c)  +  2 pi(A^c)
```

**3.1 Bounding the conditional TV term.** `P^n_A(x,\cdot)` and
`pi_A(\cdot)` are the sampler/target laws truncated to `A` and renormalized
— not, in general, literally Gaussian (or any other tractable closed form)
even when the untruncated density is. Two routes, likelihood-agnostic in
statement though not in execution:

*Route 1 — Pinsker + truncated KL:*
`||P^n_A-pi_A||_TV <= sqrt(0.5 D_KL(P^n_A || pi_A))`, computing the KL
divergence between the two densities truncated to the common region `A`.

*Route 2 — domination by an untruncated comparison density:* bound the true
truncated target on `A` above/below by an untruncated comparison density
built from the curvature floor `kappa_j` (a Gaussian, in every case
considered here, since `H_j(beta_j) \succeq kappa_j I` on `A` by
definition), then apply a Gaussian-comparison theorem to that comparison
pair, absorbing the domination gap as additional slack.

Neither route is carried to a closed form in this note; §5–§8 each note the
concrete comparison object Route 2 would use for that likelihood.

---

## 4. Bounding `pi(A^c)`: the target escape probability

**Lemma 2 (Gaussian tail domination via the curvature floor).**
*Suppose the block's negative log-density satisfies*
```
-log pi(beta_j | rest)  >=  0.5 kappa_j ||beta_j - beta_j_hat||^2_{D_j'D_j}  +  const
```
*for `beta_j` outside a neighborhood of the anchor `beta_j_hat` (this holds
whenever `H_j(beta_j) \succeq kappa_j D_j'D_j / \lambda_max(D_j'D_j)`-type
bounds are available, or more simply whenever `c(eta) \geq 0` and the prior
alone supplies `kappa_j = \tau_j^{-2}`, as in every likelihood below).
Then*
```
pi(beta_j | rest)  <=  C_j exp( -0.5 kappa_j ||beta_j-beta_j_hat||^2_{D_j'D_j} )
```
*for a constant `C_j` not depending on `beta_j`.*

*Proof.* Integrate the bound on `-log pi` and exponentiate; the constant
absorbs the normalizer and any additive terms independent of `beta_j`.
**QED**

**Lemma 3 (per-block escape probability).**
*With `q_j(beta_j) := \|beta_j-beta_j\_hat\|^2_{D_j'D_j}` and any fixed
`K_j^\star > 0`,*
```
pi( q_j(beta_j) > K_j^\star )  <=  1 - erf_{p_j}( sqrt( kappa_j K_j^\star / 2 ) )
```
*Proof.* As in the t-likelihood derivation: change of variables
`z=(D_j'D_j)^{1/2}(beta_j-beta_j_hat)` reduces `q_j` to `\|z\|^2`, Lemma 2
dominates the density of `z` by `N(0,\kappa_j^{-1}I_{p_j})`, and Remark 3 of
the Nygren paper converts the resulting Gaussian tail probability into
`erf_{p_j}`. **QED**

**Claim 2 (`pi(A^c)` via union bound).**
```
pi(A^c)  <=  sum_j [ 1 - erf_{p_j}( sqrt(kappa_j K_j^\star/2) ) ]  +  (analogous terms for population/RE blocks)
```
*Proof.* Union (Boole's) bound over the per-block "outside" events, each
bounded by Lemma 3. **QED**

This machinery — Lemmas 2–3 and Claim 2 — is stated generically in terms of
a curvature floor `kappa_j`; §5–§8 supply the concrete value of `kappa_j`
and the resulting escape bound for each likelihood.

---

## 5. Multivariate-t likelihood (truncated-Gamma scale mixture)

### 5.1 Safe region and curvature floor

`c(eta_i) = Omega_j` (a random, tilted-truncated-Gamma precision, integrated
out per `P_MATRIX_MARGINAL_PRECISIONS.md` §1). The Hessian is
```
H_j^{(\Omega)}(\beta_j) = E_t[\Omega_j]\,D_j'D_j - \mathrm{Var}_t[\Omega_j]\,(D_j'e_j)(D_j'e_j)'
```
— a rank-one **subtraction**, not merely a shrinking positive term. This is
the key structural difference from every other likelihood in this note
(§6–§8): **`H_j^{(\Omega)}(\beta_j)` is *not* PSD on the whole space.** It
is PSD if and only if `\beta_j` lies inside the ellipsoid
```
A_j = \{\beta_j : q_j(\beta_j) \le K_j^{(\Omega)}(\beta_j)\},
\qquad K_j^{(\Omega)}(\beta_j) := E_t[\Omega_j]/\mathrm{Var}_t[\Omega_j]
```
(§1.3, `P_MATRIX_MARGINAL_PRECISIONS.md`) — the "matched" condition of that
document's §5.2 is precisely the statement that `H_j^{(\Omega)}\succeq 0`
*on `A_j`*, not a global property. Outside `A_j` the rank-one term can
dominate `E_t[\Omega_j]D_j'D_j`, and the Hessian is genuinely indefinite —
this is qualitatively worse than a curvature *floor* that merely shrinks
toward zero (as in Poisson/probit, §7–§8), because an indefinite Hessian
means the local map can be **expansive**, not just slow-mixing (see the
`S=2>1` counterexample, §5.2 of `P_MATRIX_MARGINAL_PRECISIONS.md`, and item
4 of §9 below).

So the safe region `A` for this section is *defined by* the PSD boundary
itself — `A = \bigcap_j A_j` (intersected with the analogous population/RE
ellipsoids) — rather than by an arbitrarily chosen sublevel set of a
curvature that is PD everywhere but merely small outside the threshold.
`A^c` is genuinely "the region with no PD guarantee at all," not "the
region with a weaker PD guarantee." Everything downstream (§5.2's drift
argument, in particular) has to supply its own reason for the chain not to
diverge on `A^c` — it cannot lean on a residual positive-definite floor the
way §7–§8 do.

The **fixed** (state-independent) radius `K_j^\star` used in §4's generic
machinery is the conservative floor
`K_j^\star := \inf_{\beta_j} K_j^{(\Omega)}(\beta_j) > 0` from Remark 1/2 of
the previous revision of this note; the associated curvature floor plugged
into Lemma 2 (§4) is `\kappa_j = \omega_{L,j}`, the left endpoint of the
truncation window on `\Omega_j` — but note this `\omega_{L,j}` bound is
what controls the **target's** tail (§4's Gaussian-domination argument,
which only ever uses `E_t[\Omega_j]D_j'D_j`, the always-PSD part), and is
logically separate from the PSD/indefiniteness question above, which
governs the **sampler's** local contraction, not the target density's tail
mass.

### 5.2 Bounding `P^n(x,A^c)`: drift/Lyapunov argument

**Lemma 4 (Markov on the ellipsoid-score Lyapunov function).**
Let `V(\gamma,\beta) := 1+\sum_j q_j(\beta_j)+\sum_p q_p^{(\lambda)}`,
`K^\star := \min(\min_j K_j^\star, \min_p K_p^{(\lambda)\star})`. Then
`P^n(x,A^c) \le E[V(X_n)\mid x]/K^\star` (Markov's inequality, since
`A^c \subseteq \{V>K^\star\}`).

**Claim 3 (corner-domination drift bound).** If the truncation-box-corner
spectral radius `\lambda^\star_{corner}` (§10, `BLOCK_GIBBS_ERGODICITY_ING.md`)
dominates the true state-dependent contraction throughout the box,
```
E[V(X_n)\mid x] \le (\lambda^\star_{corner})^{2n} V(x) + V_\infty,
\quad V_\infty := 1+\sum_j tr(D_j'D_j\,\Sigma_{11,j}^{\infty,corner})+\sum_p tr(\Sigma_{11,p}^{\infty,corner})
```
via Claims 2–3 of the Nygren paper applied at the corner covariance.

**Theorem 1 (t-likelihood sampler escape bound).**
```
P^n(x,A^c)  <=  [ (\lambda^\star_{corner})^{2n} V(x) + V_\infty ] / K^\star
```
*Status:* the marginal-vs-conditional composition step is open (§8 of the
prior revision; carried forward as item 2 in §9 below).

### 5.3 Assembled bound

Using the Corollary to Lemma 6 (§3) as originally stated:
```
|| P^n(x,.) - pi(.) ||_TV
   <=  || P^n_A(x,.) - pi_A(.) ||_TV
     + 2 [ (\lambda^\star_{corner})^{2n} V(x) + V_\infty ] / K^\star
     + 2 \sum_j [ 1 - erf_{p_j}(\sqrt{\omega_{L,j} K_j^\star/2}) ]
     + 2 \sum_p [ 1 - erf_J(\sqrt{\lambda_{L,p} K_p^{(\lambda)\star}/2}) ]
```
Or, via the refined split of §2.1 (Claim 5), replacing the last three terms
(the `2P^n(x,A^c)+2\pi(A^c)` piece) with `4P^n(x,A^c)+2R\rho^nV(x)/K^\star`:
```
|| P^n(x,.) - pi(.) ||_TV
   <=  || P^n_A(x,.) - pi_A(.) ||_TV
     + 4 [ (\lambda^\star_{corner})^{2n} V(x) + V_\infty ] / K^\star
     + 2 R \rho^n V(x) / K^\star
```
where the first bracketed term (`4\times` the Theorem 1 bound) can now be
**estimated empirically** as `4\times`(proportion of post-burn-in draws
outside the safe region), rather than relying on the closed-form
`\mathrm{erf}`-based bound on `\pi(A^c)`; the remaining `2R\rho^nV(x)/K^\star`
is the provably-geometric correction, and for this section `(R,\rho)` can
in principle be sharpened past the abstract Meyn–Tweedie constants using
the exact corner-covariance rate `\lambda^\star_{corner}` already computed
for Theorem 1, since the underlying chain is exactly Gaussian conditional
on the precision path.

### 5.4 Open items (t-likelihood specific)

1. Composition/domination step, §5.2 (marginal vs. fixed-point sub-chain).
2. The conditional TV term of §3.1 (Route 1/2) not yet carried out.
3. Union bound (Claim 2, §4) discards inter-block correlation via `\gamma`.
4. Unmatched-block regime: if `H_j^{(\Omega)}` can be indefinite inside the
   chosen ellipsoid, `A` must be shrunk or the matched condition verified
   analytically at the `\omega_{L,j}` floor.

---

## 6. Logit (logistic) likelihood

### 6.1 Safe region via a pilot-calibrated eigenvalue ceiling

As before, `c(\eta) = p(1-p) \in [0,1/4]` for every finite `\eta`, so the
**floor** is global and unconditional:
```
H_j(\beta_j) \;\succeq\; \tau_j^{-2} I_{p_j} \qquad \text{for every } \beta_j
```
— this part of the earlier argument still holds everywhere and is never in
question for logit (contrast §5, §7, §8, where the floor itself can fail or
require a separate tail argument). What changes here is the **ceiling**.
The theoretical worst case,
```
M_j^{\text{theory}} := \tau_j^{-2} + \tfrac14\lambda_{\max}(D_j'D_j),
```
is attained only if `c(\eta_i)=1/4` for every `i` simultaneously — i.e.
every observation in group `j` sits exactly at `\eta_i=0` — which is a
worst case essentially never realized in practice. Rather than certifying
mixing against this loose global ceiling, use a **pilot-phase-estimated**
ceiling vector
```
M_j^{\text{pilot}} := \max_{1\le l\le L} \lambda_{\max}\!\Big(D_j'\,\mathrm{diag}\big(c(\eta_i^{(l)})\big)\,D_j\Big)
```
i.e. the largest data-curvature eigenvalue actually observed for block `j`
across `L` pilot sweeps, `l=1,\dots,L` (one scalar per block `j`, so
`M^{\text{pilot}} = (M_1^{\text{pilot}},\dots,M_m^{\text{pilot}})` is a
vector across blocks). Define the **safe region**
```
A_j := \big\{\beta_j : \lambda_{\max}\big(D_j'\,\mathrm{diag}(c(\eta_i))\,D_j\big) \le M_j^{\text{pilot}}\big\}, \qquad A := \bigcap_j A_j
```
and `A^c` its complement — exactly the safe/unsafe split of Definition 1
(§1), now with the split governed by an **upper** eigenvalue threshold
rather than a lower one, since the lower bound is never at risk here. On
`A`, the sandwich tightens to
```
\tau_j^{-2}I \;\preceq\; H_j(\beta_j) \;\preceq\; \big(\tau_j^{-2}+M_j^{\text{pilot}}\big) I
```
giving a sharper comparison-Gaussian ceiling than `M_j^{\text{theory}}`,
hence (via the §6.2 sandwich argument below) a sharper mixing-rate
certificate on `A`. On `A^c`, only the looser global bound
`H_j(\beta_j)\preceq(\tau_j^{-2}+\tfrac14\lambda_{\max}(D_j'D_j))I` is
available, but — critically — this bound is **still valid** there; unlike
§5, §7, §8, `A^c` for logit never threatens loss of positive-definiteness
or expansion, only a looser (but still finite, still `<1`-rate) contraction
certificate.

**Caveat on `M_j^{\text{pilot}}` as a threshold.** Because
`M_j^{\text{pilot}}` is itself estimated from a finite pilot run, it is a
random quantity, and treating it as a *fixed* threshold when bounding
`P^n(x,A^c)` for a subsequent (or continued) run is only rigorous if either
(a) the pilot run is treated as separate from, and independent of, the run
being certified (a sample-splitting argument), or (b) `M_j^{\text{pilot}}`
is inflated by a concentration-based safety margin over its pilot-sample
estimate. Neither is carried out here; §6.4 records this as an open item
rather than treating the pilot ceiling as a proof-grade constant.

### 6.2 Bounding `P^n(x,A^c)` and `\pi(A^c)`

Because the *global* sandwich `[\tau_j^{-2}, M_j^{\text{theory}}]` holds
**unconditionally** — on both `A` and `A^c` — the chain is uniformly
ergodic at the (loose) global rate regardless of whether the pilot ceiling
is ever exceeded (the whole-space strongly-log-concave argument of the
previous revision, Mengersen & Tweedie 1996, still applies globally using
`M_j^{\text{theory}}`). This has a useful consequence: unlike §5, §7, §8,
you do **not** need a separate Foster–Lyapunov drift argument to rule out
divergence on `A^c` — global strong log-concavity already does that. What
you need instead is simply a bound on *how often* the chain visits `A^c`,
since that determines how often you're stuck using the loose
`M_j^{\text{theory}}` rate instead of the sharp `M_j^{\text{pilot}}` rate.

**Claim 4 (escape bound via the eigenvalue functional's tail).** Write
`g_j(\beta_j) := \lambda_{\max}(D_j'\,\mathrm{diag}(c(\eta_i))\,D_j)`, a
continuous, bounded (`\le\tfrac14\lambda_{\max}(D_j'D_j)`) function of
`\beta_j`. Then
```
\pi(A_j^c) = \pi\big(g_j(\beta_j) > M_j^{\text{pilot}}\big)
```
is, by construction of `M_j^{\text{pilot}}` as an empirical maximum over
`L` pilot draws from (approximately) `\pi` itself, estimable directly as an
**empirical tail frequency** — e.g. via a binomial/Clopper–Pearson interval
on the proportion of a held-out validation sample exceeding
`M_j^{\text{pilot}}` — rather than via the closed-form `\mathrm{erf}_{p_j}`
device of Lemma 3 (§4), which does not directly apply here since `g_j` is a
nonlinear (eigenvalue) functional of `\beta_j` rather than the quadratic
form `q_j(\beta_j)` that Lemma 2's Gaussian-domination argument was built
for. `P^n(x,A^c)` is bounded analogously as an empirical trace statistic
from the pilot chain itself (proportion of pilot sweeps with
`g_j(\beta_j^{(l)}) > M_j^{\text{pilot}}` — by construction this is exactly
`0` for the pilot sample defining the max, so a genuine bound requires the
independent/held-out sample of the Caveat in §6.1).

*Status:* this is the one departure from the closed-form, non-simulation
philosophy of §4–§5 in this note — the eigenvalue functional `g_j` does not
admit the same Gaussian-tail closed form as the quadratic ellipsoid scores
`q_j` used elsewhere, so §4's machinery is not directly reusable here
without further work (a possible closed-form route, not pursued: bound
`g_j(\beta_j) \le \tfrac14\lambda_{\max}\big(D_j'D_j\big)\cdot
\max_i\{4c(\eta_i)\}` and separately tail-bound `\max_i c(\eta_i)` via a
union bound over the `n_j` observations in group `j`, each `c(\eta_i)`
individually controllable via a Gaussian tail on `\eta_i` — flagged for
future work rather than carried out).

### 6.3 Assembled bound

```
|| P^n(x,.) - \pi(.) ||_TV
   \le \|P^n_A(x,.)-\pi_A(.)\|_{TV}
     + P^n(x,A^c) + \pi(A^c)
```
via Lemma 1 (§2) directly (the normalization-correction machinery of §3 can
be layered on top identically to §5, §7, §8, at the same cost of a factor
of 2 on the escape terms), where `\|P^n_A(x,.)-\pi_A(.)\|_{TV}` is bounded
via the *sharp* sandwich `[\tau_j^{-2}, \tau_j^{-2}+M_j^{\text{pilot}}]` on
`A` (§6.1) using the Route-2 comparison-Gaussian argument of §3.1.

Applying the §2.1 refinement here is unusually clean, because logit gets
the *uniform*-ergodicity shortcut noted at the end of §2.1's Caveat rather
than needing the abstract `V`-norm theorem:
```
2 P^n(x,A^c) + 2 \pi(A^c)
   = 4 P^n(x,A^c) + 2\big[\pi(A^c)-P^n(x,A^c)\big]
   \le 4 P^n(x,A^c) + 2 M \rho_{\text{glob}}^n
```
using `|\pi(A^c)-P^n(x,A^c)| \le \|P^n(x,.)-\pi(.)\|_{TV} \le M
\rho_{\text{glob}}^n` directly from global strong log-concavity (no `V` or
minorization-on-a-small-set bookkeeping required, since the whole space
already plays that role). So `4P^n(x,A^c)` — estimable empirically as `4
\times` the post-burn-in proportion of draws with
`g_j(\beta_j)>M_j^{\text{pilot}}` for some block `j` — replaces the
currently-empirical `\pi(A^c)` term of §6.2 as the dominant piece, with
`2M\rho_{\text{glob}}^n` a fully closed-form, explicitly computable
correction (unlike §5/§7/§8, where `(R,\rho)` are not yet made explicit).

### 6.4 Open items (logit-specific)

1. **Pilot-threshold rigor (§6.1 Caveat).** Either a sample-splitting
   argument or a concentration-based inflation of `M_j^{\text{pilot}}` is
   needed before this ceiling can be used in a proof-grade (rather than
   diagnostic) bound.
2. **Closed-form tail bound for `g_j` (§6.2).** The eigenvalue functional
   does not reduce to the quadratic-form machinery of §4; the union-bound
   sketch at the end of §6.2 is the natural next step but is not carried
   out.
3. **Sandwiching-Gaussian comparison (carried over from the previous
   revision).** Explicit coordinate-by-coordinate application of Nygren's
   Theorem 2 to the `[\tau_j^{-2},\tau_j^{-2}+M_j^{\text{pilot}}]` sandwich
   on `A` has not been written out.
4. **Full-rank / estimability (Assumption 1).** Doing real work here: if
   violated, `\hat\beta_j` (the MLE) need not exist (complete separation);
   the floor is still `\tau_j^{-2}` from the prior, but the anchor point
   for the sandwich argument must be replaced by the posterior mode
   instead of the MLE.

---

## 7. Probit likelihood

### 7.1 Safe region and curvature floor — two-region split needed

For binary `y_i` with `\eta_i = D_j[i,]\beta_j`, the probit log-likelihood
contribution has curvature
```
c(\eta) = \frac{\phi(\eta)^2}{\Phi(\eta)\Phi(-\eta)}\quad(y=1\text{ case; symmetric for }y=0)
```
(`\phi,\Phi` the standard normal pdf/cdf). This is bounded **above** by a
universal constant (`\sup_\eta c(\eta) \approx 0.6366 = 2/\pi`, attained at
`\eta=0`), but — unlike logit — is **not** bounded away from `0` uniformly:
as `\eta\to-\infty` (the "wrong-direction" tail for `y=1`), `c(\eta)\to0`,
via the Mills-ratio asymptotic `\phi(\eta)/\Phi(\eta) \sim -\eta` as
`\eta\to-\infty`, giving `c(\eta)\sim\eta^2\phi(\eta)/\Phi(\eta) \to 0`
faster than any polynomial rate but with a heavier per-observation tail
footprint than logit's `p(1-p)\sim e^{-|\eta|}` (probit's tail is Gaussian,
`\phi(\eta)`, times a polynomial correction, so it decays even faster in
`\eta` but the *event* `\eta\to-\infty` is reached at smaller `\|\beta_j\|`
for a Gaussian-tailed link than for logit's exponential-tailed link — the
qualitative structure, not the precise rate, is what matters here). As with
Poisson (§8), the prior alone supplies the floor
`H_j(\beta_j)\succeq\tau_j^{-2}I` for **every** finite `\beta_j` (since
`c(\eta)\ge0` always), but the *data* contribution to curvature can vanish,
so a genuine safe/unsafe split — exactly as in §8 — is needed if you want a
tighter floor than the bare prior precision on the "good" region.

Define, as in Definition 1, `A = \{\beta_j : \lambda_{\min}(D_j'\,
\mathrm{diag}(c(\eta_i))\,D_j) \ge \kappa\}` for a chosen threshold
`\kappa>0`; outside `A`, fall back to the universal prior floor
`\tau_j^{-2}`, which — critically, exactly as in the Poisson case — is
never actually violated (curvature is always `\ge0`), so `H_j(\beta_j)
\succ \tau_j^{-2}I` strictly for every finite `\beta_j`, and the "unsafe"
region is where this floor is the *only* thing you can rely on, not a
region where the chain can expand.

### 7.2 Bounding `P^n(x,A^c)`: drift via the density tail, not the curvature comparison

Since the data-curvature contribution is always nonnegative and the target
density's tail is controlled by the same argument as Poisson (§8.2): the
prior's quadratic term `-\|\beta_j\|^2/(2\tau_j^2)` dominates the
likelihood contribution (bounded above by `0`, since `\log\Phi(\cdot)\le0`)
in every direction as `\|\beta_j\|\to\infty`, giving genuinely Gaussian
tails for the *marginal* `\pi(\beta_j)`, exactly as derived generically in
Lemma 2 (§4) with `\kappa_j=\tau_j^{-2}`. The Lyapunov function
`V(\beta_j)=\exp(a\|\beta_j\|^2)` for `a<1/(2\tau_j^2)` is legitimate here
by the same tail argument as §8.2, and the resulting drift bound has the
same form as Theorem 1 of §5, with `\kappa_j=\tau_j^{-2}` in place of
`\omega_{L,j}` and a probit-specific one-step contraction constant (bounded
by the same corner-type argument, now over the compact "safe" curvature
region only, since the unsafe region's contraction cannot be certified
below `1` — only the density-tail argument, not a curvature-comparison
argument, supplies drift there).

### 7.3 Assembled bound

Structurally identical to §5.3 and §8.3, with `\kappa_j=\tau_j^{-2}` and
`p_j` in place of the t-likelihood's `\omega_{L,j}`, `K_j^\star`. As
originally stated:
```
|| P^n(x,.) - \pi(.) ||_TV
  <=  \|P^n_A(x,.)-\pi_A(.)\|_{TV}
    + 2\,P^n(x,A^c)
    + 2\sum_j\Big[1-\mathrm{erf}_{p_j}\!\big(\sqrt{\tau_j^{-2}K_j^\star/2}\big)\Big]
```
Via the §2.1 refinement, the last two terms (`2P^n(x,A^c)+2\pi(A^c)`)
become `4P^n(x,A^c)+2R\rho^nV(x)/K^\star`, with `4P^n(x,A^c)` estimable
from the empirical post-burn-in proportion of draws in `A^c` and
`(R,\rho)` the (not yet explicit — see §9) Meyn–Tweedie constants attached
to the drift/minorization pair of §7.2.

### 7.4 Open items (probit-specific)

1. A precise (rather than asymptotic/qualitative) uniform-in-`\beta_j`
   lower bound on `c(\eta)` restricted to the safe region `A` has not been
   derived in closed form (unlike logit's exact `p(1-p)` bound); this is
   needed to make `\kappa` (the safe-region threshold) and the
   corner-contraction constant of §7.2 concrete numbers rather than
   qualitative statements.
2. Whether the Mills-ratio tail decay is fast enough to give a *better*
   (not just qualitatively similar) escape-probability bound than
   Poisson's is not investigated here.

---

## 8. Poisson likelihood

### 8.1 Safe region and curvature floor

`c(\eta)=e^{\eta}`, unbounded above (helps drift) but `\to0` as
`\eta\to-\infty` (curvature vanishes toward, but never below, the prior
floor). `H_j(\beta_j) = \tau_j^{-2}I + D_j'\mathrm{diag}(e^{\eta_i})D_j
\succ \tau_j^{-2}I` strictly for every finite `\beta_j`. Define
`A = \{\beta_j : \lambda_{\min}(D_j'\,\mathrm{diag}(e^{\eta_i})\,D_j)\le
\kappa\}^c \cap \{\|\beta_j\|\le R\}^c` — i.e. the *unsafe* set is the
intersection of "low data curvature" and "not yet rescued by the density
tail", per the three-regime discussion earlier in this conversation:
(i) `\eta` large positive — strong drift; (ii) `\eta` moderate/negative but
`\|\beta_j\|` bounded — the genuine problem region, handled by minorization
on a compact set; (iii) `\|\beta_j\|\to\infty` in any direction — prior
quadratic term wins regardless of curvature, drift resumes.

### 8.2 Bounding `P^n(x,A^c)`: drift via the density tail

`V(\beta_j)=\exp(a\|\beta_j\|^2)`, `a<\tau_j^{-2}/2`, is a valid Lyapunov
function on the *density* tail (regime (iii)) since
`\pi(\beta_j)\propto\exp(-\|\beta_j\|^2/(2\tau_j^2)+o(\|\beta_j\|^2))` —
the `-e^{\eta}` and `y\eta` terms are respectively bounded above by `0` and
growing at most linearly, so the Gaussian tail from the prior dominates.
Combined with the compact-set minorization on regime (ii) (as in §5.2 of
the earlier `P_MATRIX...` discussion — continuity and positivity of the
truncated-Gamma-free Poisson-Gaussian conditional density on any compact
set), this gives a genuine Foster–Lyapunov pair, and Markov's inequality
(Lemma 4, §5.2, restated for this `V`) gives
```
P^n(x,A^c)  <=  E[V(X_n)\mid x] / K^\star
```
with `E[V(X_n)\mid x]` bounded geometrically via the drift inequality
`PV(x)\le\lambda V(x)+b\,\mathbb 1_C(x)` (standard Meyn–Tweedie iteration,
not the exact Gaussian-AR machinery of §5, since the Poisson conditional is
not Gaussian) rather than the corner-covariance closed form of §5.2.

### 8.3 Assembled bound

As originally stated:
```
|| P^n(x,.) - \pi(.) ||_TV
  <=  \|P^n_A(x,.)-\pi_A(.)\|_{TV}
    + 2\Big[\frac{\lambda^n V(x)+b/(1-\lambda)}{K^\star}\Big]
    + 2\sum_j\Big[1-\mathrm{erf}_{p_j}\!\big(\sqrt{\tau_j^{-2}K_j^\star/2}\big)\Big]
```
(using the standard geometric-drift consequence `E[V(X_n)]\le\lambda^n V(x)
+b/(1-\lambda)` in place of §5's exact corner formula). Via the §2.1
refinement:
```
|| P^n(x,.) - \pi(.) ||_TV
  <=  \|P^n_A(x,.)-\pi_A(.)\|_{TV}
    + 4\Big[\frac{\lambda^n V(x)+b/(1-\lambda)}{K^\star}\Big]
    + 2 R \rho^n V(x)/K^\star
```
with the first bracketed term now interpretable as `4\times`(empirical
post-burn-in proportion of draws in `A^c`) in practice, and the correction
`2R\rho^nV(x)/K^\star` sharing the same not-yet-explicit-constants caveat
as §7.3.

### 8.4 Open items (Poisson-specific)

1. Making regime (ii)'s minorization constant and regime (i)/(iii)'s drift
   constants `(\lambda,b)` explicit, rather than asserted to exist by
   continuity/compactness — this is exactly the "make $R$ concrete" task
   raised earlier in this conversation and not yet carried out.
2. Verifying regimes (i)–(iii) actually cover `A^c` with no gap (also
   raised earlier, not yet completed).

---

## 9. Cross-cutting open items (all likelihoods)

1. **Conditional TV term (§3.1).** Neither Route 1 nor Route 2 is carried
   to a closed form for any likelihood; §6.2 sketches Route 2 furthest (for
   logit) but does not finish the coordinate-by-coordinate application of
   Nygren's Theorem 2 to a sandwiching pair.
2. **t-likelihood composition step (§5.2).** Marginal-vs-fixed-point
   domination not yet verified.
3. **Union bound looseness (Claim 2, §4).** Discards inter-block
   correlation via `\gamma` in every likelihood section.
4. **Unmatched-block / indefinite-Hessian regime (t-likelihood only, §5.4
   item 4).** Does not arise for logit (globally PSD), probit, or Poisson
   (both always `\succeq \tau_j^{-2}I \succ 0`), since only the
   rank-one-*subtraction* structure of the t-likelihood's marginalized
   Hessian can produce indefiniteness.
5. **Explicit `(R,\rho)` constants for the §2.1 refinement.** Claim 5's
   Meyn–Tweedie `V`-norm theorem guarantees `R,\rho` exist given
   drift+minorization, but they are not reduced to closed form here for
   §7 (probit) or §8 (Poisson) — both currently depend on drift constants
   `(\lambda,b)` and a minorization constant `\epsilon` on a small set that
   are themselves flagged as open (§7.4, §8.4). §5 (t-likelihood) is the
   one case where this could plausibly be sharpened past the abstract
   theorem, using the exact corner-covariance rate already in hand.

---

## References

- Nygren, K. *On the total variation distance between multivariate normal
  densities with applications to two-block Gibbs samplers.*
  (`Multivariate_Normal_Distances_07_20_20.pdf`) — Definitions 1–7, Lemma 1,
  Lemma 2, Theorem 2, Claim 2, Claim 3, Theorem 3, Corollary 1.
- `BLOCK_GIBBS_ERGODICITY.md` — scalar/vector two-block spectral radius
  derivation, full-rank identifiability conditions, Future Work §.
- `P_MATRIX_MARGINAL_PRECISIONS.md` §1, §5 — ellipsoid scores, tilted
  truncated-Gamma moments, matched vs. unmatched Hessian conditions.
- `BLOCK_GIBBS_ERGODICITY_ING.md` §10 (cited, not reproduced) — worst-case-
  over-truncation-box monotonicity/Riccati argument.
- Mengersen, K.L. and Tweedie, R.L. (1996). *Rates of convergence of the
  Hastings and Metropolis algorithms.* Annals of Statistics — strongly
  log-concave / uniform ergodicity result invoked in §6.2.
