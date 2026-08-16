# Hierarchical Generalized Linear Model Notation

This document defines the notation used internally by `glmerb` and its
relationship to standard mixed-model notation (`lme4::glmer`) and to the
classical Bayesian hierarchical-model literature (Lindley & Smith, 1972).

It is the generalized-linear counterpart of
`inst/HIERARCHICAL_LINEAR_MODEL_NOTATION.md`. **Stages 2 and 3 are
unchanged**; only the first-stage likelihood is replaced. That single
substitution is what forces every downstream difference, and §7 collects
them in one table. Read the two documents side by side: sections are
numbered to correspond.

## 1. Model

For group $j = 1,\dots,J$, with $n_j$ observations per group and $P$
group-varying coefficients:

**Likelihood (stage 1):** the observations are conditionally independent
draws from an exponential-dispersion family with linear predictor
$\eta_j = D_j\beta_j + o_j$ and link $g$:

$$p(y_{ij} \mid \beta_j, \phi_j) = \exp\left\{\frac{a_{ij}\big(y_{ij}\theta_{ij} - b(\theta_{ij})\big)}{\phi_j} + c\big(y_{ij}, \phi_j/a_{ij}\big)\right\}$$

$$\mu_{ij} = b'(\theta_{ij}) = E[y_{ij}\mid\beta_j], \qquad g(\mu_{ij}) = \eta_{ij} = d_{ij}'\beta_j + o_{ij}, \qquad \operatorname{Var}(y_{ij}\mid\beta_j) = \frac{\phi_j}{a_{ij}}V(\mu_{ij})$$

Setting $g = \text{identity}$, $V \equiv 1$, $a_{ij}\equiv1$,
$\phi_j = \sigma_j^2$ recovers the linear model
$y_j\mid\beta_j \sim N(D_j\beta_j, \sigma_j^2 I_{n_j})$ exactly.

**Hierarchical prior (stage 2)** — *identical to the linear case*:

$$\beta_j \mid \gamma, \Psi \ \sim\ N\big(\mathcal{W}_j\gamma,\ \Psi\big)$$

**Population prior (stage 3)** — *also unaffected by the family*:

$$\gamma \ \sim\ N\big(\mu_0,\ \Lambda_\gamma^{-1}\big)$$

(The linear document leaves this stage implicit, quoting only the
flat-prior GLS form for $\gamma$; it is stated explicitly here because
§4 and §8 both need $\Lambda_\gamma$, and because $\Lambda_\gamma \to 0$
is exactly where the two documents' $\gamma$ blocks coincide.)

The full, non-centered coefficient vector $\beta_j$ — not a mean-zero
deviation — enters the **linear predictor** directly. This is a
**centered parameterization** in the sense of Lindley & Smith (1972):
$\beta_j$ is the first-stage parameter, $\gamma$ is the second-stage
(population/hyper-mean) parameter, and $\Psi$ is the second-stage
covariance.

> **The stage-2 prior is Gaussian even though the stage-1 likelihood is
> not.** This is the structural fact the whole sampler is built on. The
> non-Gaussianity is confined to the $\beta$ block; the $\gamma$ and
> $\Psi$ blocks see the data *only through* $\beta$, so they are
> conjugate exactly as in the linear model. See §4 and §7.

## 2. Symbol table

Rows marked **new** have no counterpart in the linear document; rows
marked **changed** replace a linear-model symbol.

| Symbol | Name | Dimension | Role |
|---|---|---|---|
| $y_j$ | response, group $j$ | $n_j \times 1$ | data |
| $D_j$ | likelihood design matrix, group $j$ | $n_j \times P$ | multiplies $\beta_j$ **in the linear predictor**, not in the mean |
| $\beta_j$ | group-level coefficients | $P \times 1$ | **not mean zero**; sampled directly against $y_j$ |
| $\eta_j = D_j\beta_j + o_j$ | linear predictor, group $j$ | $n_j\times1$ | **new** — the layer at which the design acts |
| $o_j$ | offset, group $j$ | $n_j\times1$ | **new** — known additive term on the link scale |
| $g,\ g^{-1}$ | link and inverse link | — | **new** — $\mu_{ij}=g^{-1}(\eta_{ij})$ |
| $\mu_j$ | conditional mean, group $j$ | $n_j\times1$ | **new** — nonlinear in $\beta_j$ unless $g$ is the identity |
| $V(\cdot)$ | variance function | — | **new** — family-determined |
| $a_{ij}$ | prior weight (binomial trials, `weights=`) | scalar | **new** |
| $\phi_j$ | dispersion, group $j$ | scalar | **changed** — replaces $\sigma_j^2$; fixed at $1$ for binomial/Poisson, free for Gaussian/Gamma |
| $\ell_j(\beta_j)$ | group $j$ log-likelihood | scalar | **new** — $\sum_i \log p(y_{ij}\mid\beta_j,\phi_j)$; concave in $\beta_j$ for canonical links |
| $\Omega_j(\beta_j)$ | GLM weight matrix | $n_j\times n_j$ diagonal | **new** — $\omega_{ij} = \dfrac{a_{ij}}{\phi_j}\dfrac{(d\mu/d\eta)_{ij}^2}{V(\mu_{ij})}$ |
| $\mathcal{G}_j(\beta_j) = D_j'\Omega_j(\beta_j)D_j$ | group information | $P\times P$ | **new** — replaces the *constant* $\sigma_j^{-2}D_j'D_j$ |
| $\mathcal{W}_j$ | level-2 design matrix, group $j$ | $P \times q$ | block-diagonal; links $\gamma \to E[\beta_j]$ |
| $\gamma$ | population (hyper-mean) coefficients | $q \times 1$ | shared across all groups |
| $\mu_0,\ \Lambda_\gamma$ | population prior mean and precision | $q\times1$, $q\times q$ | stage 3; $\Lambda_\gamma = 0$ is the flat limit |
| $\Psi$ | covariance of $\beta_j$ about $\mathcal{W}_j\gamma$ | $P \times P$ | shrinkage/hierarchical covariance |
| $u_j := \beta_j - \mathcal{W}_j\gamma$ | derived deviation | $P \times 1$ | mean zero; **not sampled directly** — byproduct used only in the $\Psi$ update |

> **Notation clash to watch.** The GLM literature writes $W$ for the
> diagonal IRLS weight matrix, and this document (following
> `notation.md`) writes $\mathcal{W}_j$ for the **level-2 design**. They
> are unrelated. To keep them apart the GLM weights are written
> $\Omega_j$ here. Other notes in `inst/` use the conventional $W_j$ for
> the weights — and `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` goes
> further and reuses $\mathcal{W}_j$ for $Z_j^\top W_jZ_j$, the group
> information. §8 gives the full dictionary; the collision is real and
> has bitten before.

$\mathcal{W}_j$ is block-diagonal across the $P$ coefficient dimensions,
exactly as in the linear case:

$$\mathcal{W}_j = \text{blockdiag}\big(W_{1j}, W_{2j}, \dots, W_{Pj}\big), \qquad W_{pj} \in \mathbb{R}^{1\times q_p},\quad q = \sum_{p=1}^P q_p$$

Each coefficient dimension $p$ may have its own set (and number) of
level-2 predictors — $q_p$ need not be equal across $p$, and $q_p = 1$
(intercept only) is the natural default when a coefficient has no
level-2 predictors.

### Score and information

Written out for reference, since these replace the sufficient statistics
$D_j'y_j$ and $D_j'D_j$ of the linear model:

$$\nabla\ell_j(\beta_j) = D_j'\,s_j(\beta_j), \qquad s_{ij} = \frac{a_{ij}\,(y_{ij}-\mu_{ij})}{\phi_j\,V(\mu_{ij})}\left(\frac{d\mu}{d\eta}\right)_{ij}$$

$$-E\big[\nabla^2\ell_j(\beta_j)\big] = \mathcal{G}_j(\beta_j) = D_j'\Omega_j(\beta_j)D_j$$

For a **canonical** link $(d\mu/d\eta) = V(\mu)$, so
$\omega_{ij} = a_{ij}V(\mu_{ij})/\phi_j$, the observed and expected
information coincide, and $\ell_j$ is concave in $\beta_j$. Under the
identity link with $V\equiv1$ this collapses to
$\mathcal{G}_j = \sigma_j^{-2}D_j'D_j$ — constant in $\beta_j$, which is
precisely why the linear model is conjugate.

## 3. Internal storage

Unchanged from the linear case except for the family/dispersion objects:

| Quantity | Storage | Shape |
|---|---|---|
| Stacked likelihood design | `Dmat` | $n \times P$ |
| Group index | `group` | length-$n$ vector, $g_i = j$ |
| Level-2 designs, one per coefficient | `Wlist = list(W_1, ..., W_P)` | `Wlist[[p]]` is $J \times q_p$ |
| Dimension bookkeeping | `qvec = sapply(Wlist, ncol)`, `q = sum(qvec)` | — |
| Family (link, variance function, `aic`) | `family` | `stats::family` object |
| Prior weights $a_{ij}$ / offsets $o_{ij}$ | `weights`, `offset` | length-$n$ |
| Dispersion mode for $\phi_j$ | `dispersion_ranef` | `"fixed"`, `"fixed_vector"`, `"gamma"`, `"gamma_list"` |
| Block-1 prior objects (one per group) | `pfamily_list` / `prior_list` | list of `dNormal()` etc. |

Recovery rules (identical):

- $D_j$ = `Dmat[group == j, ]`
- $\mathcal{W}_j$ = `blockdiag(W_1[j, ], W_2[j, ], ..., W_P[j, ])`, assembled on demand (row $j$ from each list element)

`Dmat` is stored row-stacked, **not** as the block-diagonal
$n \times PJ$ object `lme4` implicitly uses for $Z$; the block-diagonal
structure is recovered implicitly via `group` rather than materialized.

`Wlist`'s heterogeneous column widths ($q_p$ varying by $p$) are the
natural fix for the raggedness of level-2 predictor sets across
coefficient dimensions — no padding needed, unlike a 3-index array.

## 4. Full conditionals (Gibbs)

**$\beta_j$** (**not** conjugate — this is the one block that changes):

$$\log p(\beta_j \mid \cdot) = \ell_j(\beta_j) - \tfrac12(\beta_j - \mathcal{W}_j\gamma)'\Psi^{-1}(\beta_j - \mathcal{W}_j\gamma) + \text{const}$$

There is no closed-form normal draw. What survives instead:

- **Log-concavity.** For a canonical link $\ell_j$ is concave, so the
  conditional is a log-concave density on $\mathbb{R}^P$ — strongly
  log-concave, since the stage-2 term contributes $-\Psi^{-1}$. This is
  what makes accept-reject sampling against a Gaussian envelope viable,
  and it is the route `glmbayesCore::rglmb()` takes.
- **Curvature.** $-\nabla^2\log p(\beta_j\mid\cdot) = \mathcal{G}^{\text{obs}}_j(\beta_j) + \Psi^{-1} \succeq \Psi^{-1}$,
  giving the Brascamp–Lieb bound
  $\operatorname{Cov}(\beta_j\mid\cdot) \preceq \Psi$: the conditional is
  never more dispersed than the stage-2 prior.
- **Gaussian analogue.** The mode-based approximation
  $N\big(\hat\beta_j, (\mathcal{G}_j(\hat\beta_j) + \Psi^{-1})^{-1}\big)$
  is the exact conditional in the linear case and a proposal/envelope
  otherwise; note $\mathcal{G}_j(\hat\beta_j)$ has replaced the constant
  $\sigma_j^{-2}D_j'D_j$.

**$\gamma$** (**still conjugate normal** — unchanged from the linear
case): given $\{\beta_j\}$ and $\Psi$, the likelihood contributes nothing,
because $y \perp \gamma \mid \beta$. So with the stage-3 prior
$\gamma\sim N(\mu_0,\Lambda_\gamma^{-1})$,

$$\gamma \mid \cdot \ \sim\ N\big(m,\ P_{11}^{-1}\big), \qquad P_{11} = \Lambda_\gamma + \sum_j \mathcal{W}_j'\Psi^{-1}\mathcal{W}_j, \qquad m = P_{11}^{-1}\Big(\Lambda_\gamma\mu_0 + \sum_j \mathcal{W}_j'\Psi^{-1}\beta_j\Big)$$

and in the flat limit $\Lambda_\gamma \to 0$ this is the GLS estimator of
the linear document,

$$\hat\gamma = \Big(\sum_j \mathcal{W}_j'\Psi^{-1}\mathcal{W}_j\Big)^{-1}\sum_j \mathcal{W}_j'\Psi^{-1}\beta_j$$

with covariance $\big(\sum_j \mathcal{W}_j'\Psi^{-1}\mathcal{W}_j\big)^{-1}$.
**No family, link, dispersion or weight appears in this block.**

Because $\mathcal{W}_j$ is block-diagonal, $\mathcal{W}_j'\Psi^{-1}\mathcal{W}_j$
and $\mathcal{W}_j'\Psi^{-1}\beta_j$ decompose into $P\times P$ /
$P\times1$ sub-blocks built directly from `Wlist[[p]][j, ]` and the
$(p,p')$ entries of $\Psi^{-1}$ — implement this way rather than
assembling the dense $\mathcal{W}_j$ each iteration.

**$\Psi$** (unchanged): updated from $u_j = \beta_j - \mathcal{W}_j\gamma$
(derived, not sampled) via the usual inverse-Wishart / truncated-gamma
conjugate update. As with $\gamma$, the family never enters.

**$\phi_j$** (family-dependent): fixed at $1$ for binomial and Poisson,
so there is no block at all. For Gaussian and Gamma it is a free
dispersion and behaves as $\sigma_j^2$ did, under whichever
`dispersion_ranef` mode is in force.

## 5. Correspondence to `lme4::glmer`

The design decomposition still holds **exactly**, but one layer down — at
the linear predictor rather than at the response:

$$\eta_j = D_j\beta_j = D_j(\mathcal{W}_j\gamma + u_j) = \underbrace{(D_j\mathcal{W}_j)}_{X_j}\gamma + \underbrace{D_j}_{Z_j}u_j$$

which is precisely the `glmer` linear predictor $X_j\gamma + Z_j u_j$.
The mapping of objects is therefore identical to the linear case:

| `glmerb` object | `lme4::glmer` object | Relationship |
|---|---|---|
| $D_j$ | — | consumed into **both** $X_j$ and $Z_j$ below |
| $\mathcal{W}_j$ | — | consumed into $X_j$ |
| $D_j\mathcal{W}_j$ | $X_j$ | fixed-effects design block, group $j$ |
| $D_j$ | $Z_j$ | random-effects design block, group $j$ |
| $\gamma$ | $\beta$ | fixed effects |
| $u_j = \beta_j - \mathcal{W}_j\gamma$ | $b_j$ | random effects (mean zero) |
| $\Psi$ | per-group block of $G$ | $G = I_J \otimes \Psi$ |
| $\phi_j$ | `sigma(fit)^2` (Gaussian/Gamma) or fixed at $1$ | dispersion |
| $\beta_j = u_j + \mathcal{W}_j\gamma$ | $\hat\beta + \hat b_j$ | coincides with `lme4::coef()` when $\mathcal{W}_j$ is intercept-only ($q_p = 1\ \forall p$), and generalizes it (via the $\mathcal{W}_j$ weighting) otherwise |

**Note the asymmetry** (unchanged): $D_j$ has no single `lme4`
counterpart — it is the common ancestor of both $X_j$ (after
right-multiplication by $\mathcal{W}_j$) and $Z_j$ (used as-is).
`glmerb`'s parameterization retains the full, non-centered $\beta_j$ in
the linear predictor, whereas `glmer` works with $b_j$ and $\beta$.

**What does *not* carry over: marginalization.** In the linear document
§5 substitutes stage 2 into stage 1 and lands on a model that is still
Gaussian, so the marginal likelihood of $\gamma$ is available in closed
form and `lmer` maximizes it exactly. Here the substitution happens
inside $g^{-1}$, so

$$p(y_j \mid \gamma, \Psi) = \int p(y_j \mid \beta_j)\, N(\beta_j; \mathcal{W}_j\gamma, \Psi)\, d\beta_j$$

has no closed form. `glmer` approximates this integral (Laplace, or
adaptive Gauss–Hermite with `nAGQ > 1`); the sampler here does not
approximate it at all, but instead avoids it by conditioning on
$\beta_j$. This is the single most important practical difference
between the two documents, and it drives §6.

## 6. Equivalence test (for package test suite)

When $\mathcal{W}_j$ is intercept-only for every coefficient dimension
$p$ (i.e., no level-2 predictors, $q_p = 1\ \forall p$), the model
reduces to a standard random-intercept/random-slope `glmer` model. In
this special case:

- $\hat\gamma$ (posterior mean) should track `glmer`'s fixed-effect
  estimates
- $\hat\beta_j - \hat\gamma$ should track `lme4::ranef()`
- $\hat\beta_j$ should track `lme4::coef()`

**But the tolerance must be looser than in the linear case, and the
disagreement is not all sampler error.** Three distinct gaps:

1. **Integration error in the reference.** `glmer` reports the maximizer
   of a *Laplace-approximated* marginal likelihood by default. Tighten
   with `nAGQ` (scalar random effects only) before blaming the sampler.
2. **Posterior mean vs. mode.** `glmer` returns a (penalized) maximizer;
   the sampler returns a posterior mean. These differ at $O(1/J)$ for
   skewed conditionals — small counts, near-boundary binomial cells.
3. **The prior is not flat.** Stage 3 contributes $\Lambda_\gamma$ and
   the $\Psi$ update contributes its own hyperprior, neither of which
   `glmer` has. Set $\Lambda_\gamma \to 0$ and a weak $\Psi$ prior for a
   like-for-like comparison.

For a **tight** check, use `family = gaussian(link = "identity")`, where
the model reduces exactly to the linear document and all three gaps
vanish; `glmerb(family = gaussian())` is routed through the identical LMM
engine (`README_LMERB_GLMERB_FRONT_DOOR.md` §1), so this doubles as a
route-consistency test.

## 7. What changes, and what does not

The single substitution in stage 1 propagates as follows.

| Component | Linear model | Generalized linear model |
|---|---|---|
| Stage-1 mean | $E[y_j] = D_j\beta_j$ | $E[y_j] = g^{-1}(D_j\beta_j + o_j)$ |
| Stage-1 sufficient statistics | $D_j'y_j$, $D_j'D_j$ | none; replaced by $\nabla\ell_j$, $\mathcal{G}_j(\beta_j)$ — **state-dependent** |
| Stage 2, stage 3 | $\beta_j\mid\gamma\sim N(\mathcal{W}_j\gamma,\Psi)$; $\gamma\sim N(\mu_0,\Lambda_\gamma^{-1})$ | **identical — the family never appears** |
| $\beta_j$ block | conjugate normal, closed form | log-concave, no closed form; accept-reject |
| $\gamma$ block | $N(m, P_{11}^{-1})$ | **identical** — $y\perp\gamma\mid\beta$ |
| $\Psi$ block | conjugate | **identical** |
| Dispersion | $\sigma_j^2$ | $\phi_j$; fixed at $1$ for binomial/Poisson |
| Design decomposition $X_j = D_j\mathcal{W}_j$, $Z_j = D_j$ | at the response | **at the linear predictor** — same algebra |
| Marginal likelihood $p(y\mid\gamma,\Psi)$ | closed form (Gaussian) | intractable integral |
| `lme4` reference | `lmer`, exact | `glmer`, Laplace/AGQ-approximate |
| $\operatorname{Cov}(\beta_j\mid\cdot)$ | $(\sigma_j^{-2}D_j'D_j + \Psi^{-1})^{-1}$, constant | $\preceq \Psi$ (Brascamp–Lieb), varies with $\gamma$ |

The two-line summary: **the hierarchy is untouched, the likelihood block
is not.** Everything that made $\gamma$ and $\Psi$ easy in the linear
model was a property of stages 2 and 3, which are Gaussian in both
documents; everything that was closed-form about $\beta_j$ was a property
of stage 1, which is not.

## 8. Cross-reference to the theory notes

The ergodicity and minorization notes in `inst/` use a compressed
notation for the same model, and they do **not** all agree with each
other. Read the warning below before the dictionary.

> **$\mathcal{W}_j$ means two different things across `inst/`.** In this
> document and in `RESTRICTED_GIBBS_MINORIZATION_TV.md` it is the
> **level-2 (hyper) design**, $P\times q$, there written interchangeably
> as $H_j$. In `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §I.2 it is
> the **group data-precision Gram matrix**
> $\mathcal{W}_j := Z_j^\top W_j Z_j$, which is $p_{\mathrm{re}}\times
> p_{\mathrm{re}}$ and is what this document calls $\mathcal{G}_j$. The
> two are never the same object and rarely even the same shape. When
> reading C03, take $H_j$ as the hyper-design and treat its
> $\mathcal{W}_j$ as $\mathcal{G}_j$.

| Here | Theory notes | Notes |
|---|---|---|
| $\mathcal{W}_j$ | $H_j$ (both notes); also $\mathcal{W}_j$ in `RESTRICTED_GIBBS_...` only | hyper-design. **Prefer $H_j$ when writing across notes** |
| $D_j$ | $Z_j$ | within-group GLM design; agrees with the `lme4` $Z_j$ of §5 |
| $P$ | $p_{\mathrm{re}}$ | number of group-varying coefficients |
| $\Psi^{-1}$ | $P_b$ | stage-2 precision |
| $\Omega_j(\beta_j)$ | $W_j = \mathrm{diag}(w_i)$ | GLM weights, $w_i = \mathrm{wt}_i[\mu'(\eta_i)]^2/(V(\mu_i)\phi)$ — the package convention, computed by `two_block_mode_weights()` |
| $\mathcal{G}_j(\beta_j) = D_j'\Omega_jD_j$ | $\mathcal{G}_j$ in `RESTRICTED_GIBBS_...`; $\mathcal{W}_j = Z_j^\top W_jZ_j$ in `CHAPTER_C03` | group information — see the warning above |
| $\mathcal{G}_j + \Psi^{-1}$ | $B_j$ | one-step RE conditional precision; block $j$ of $P_{22}$ |
| $\sum_j \mathcal{W}_j'\Psi^{-1}\mathcal{W}_j$ | $P_{11}^{\mathrm{RE}} = \sum_j H_j^\top P_bH_j$ | random-effect part of the $\gamma$ precision |
| $\Lambda_\gamma + P_{11}^{\mathrm{RE}}$ | $P_{11}$ | full $\gamma$ precision; $\Sigma^\star = P_{11}^{-1}$ |
| $m$ (conditional mean of $\gamma$ given $\beta$) | $m(\beta)$ | the **mean map**; $M(\gamma) = E[m(\beta)\mid\gamma,y]$ |
| $\ell_j(\beta_j)$ | $\ell_j(\beta_j)$ | group log-likelihood |
| $\operatorname{Cov}(\beta_j\mid\gamma,y)$ | $V_j(\gamma)$ | drives the curvature results of `RESTRICTED_GIBBS_MINORIZATION_TV.md` §7.3.2A |

Relevant hypotheses in those notes, stated in the symbols above:

- **(H2)** each $\ell_j$ is concave in $\beta_j$ — automatic for
  canonical links.
- **(H3a)** $P_{11}^{\mathrm{RE}} = \sum_j\mathcal{W}_j'\Psi^{-1}\mathcal{W}_j \succ 0$,
  i.e. $\bigcap_j\ker\mathcal{W}_j = \{0\}$ — the hyper-design has full
  rank across groups. Required for the flat limit $\Lambda_\gamma = 0$
  to have a proper $\gamma$ block at all.
- **(H3b)** group-wise estimability — each $\ell_j$ has a finite
  maximizer, i.e. a finite `glm()` fit with positive weights.

---

*Reference: Lindley, D.V. and Smith, A.F.M. (1972). "Bayes Estimates
for the Linear Model." Journal of the Royal Statistical Society,
Series B, 34(1), 1-41.*

*See also: `inst/HIERARCHICAL_LINEAR_MODEL_NOTATION.md` (the Gaussian
case), `inst/README_LMERB_GLMERB_FRONT_DOOR.md` (the shared call chain),
`inst/RESTRICTED_GIBBS_MINORIZATION_TV.md` (ergodicity of the two-block
sampler in this notation).*
