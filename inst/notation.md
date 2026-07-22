# Notation Reference: `lmerb` Model Specification

This document defines the notation used internally by `lmerb` and its
relationship to standard mixed-model notation (`lme4`/Laird-Ware) and to
the classical Bayesian hierarchical-model literature (Lindley & Smith,
1972).

## 1. Model

For group $j = 1,\dots,J$, with $n_j$ observations per group and $P$
group-varying coefficients:

**Likelihood (stage 1):**

$$y_j \mid \beta_j, \sigma_j^2 \ \sim\ N\big(D_j\beta_j,\ \sigma_j^2 I_{n_j}\big)$$

**Hierarchical prior (stage 2):**

$$\beta_j \mid \gamma, \Psi \ \sim\ N\big(\mathcal{W}_j\gamma,\ \Psi\big)$$

The full, non-centered coefficient vector $\beta_j$ — not a mean-zero
deviation — appears directly in the likelihood. This is a **centered
parameterization** in the sense of Lindley & Smith (1972): $\beta_j$
is the first-stage parameter, $\gamma$ is the second-stage
(population/hyper-mean) parameter, and $\Psi$ is the second-stage
covariance.

## 2. Symbol table

| Symbol | Name | Dimension | Role |
|---|---|---|---|
| $y_j$ | response, group $j$ | $n_j \times 1$ | data |
| $D_j$ | likelihood design matrix, group $j$ | $n_j \times P$ | multiplies $\beta_j$ in the likelihood |
| $\beta_j$ | group-level coefficients | $P \times 1$ | **not mean zero**; sampled directly against $y_j$ |
| $\sigma_j^2$ | residual variance, group $j$ | scalar | per-group (or shared) dispersion |
| $\mathcal{W}_j$ | level-2 design matrix, group $j$ | $P \times q$ | block-diagonal; links $\gamma \to E[\beta_j]$ |
| $\gamma$ | population (hyper-mean) coefficients | $q \times 1$ | shared across all groups |
| $\Psi$ | covariance of $\beta_j$ about $\mathcal{W}_j\gamma$ | $P \times P$ | shrinkage/hierarchical covariance |
| $u_j := \beta_j - \mathcal{W}_j\gamma$ | derived deviation | $P \times 1$ | mean zero; **not sampled directly** — byproduct used only in the $\Psi$ update |

$\mathcal{W}_j$ is block-diagonal across the $P$ coefficient dimensions:

$$\mathcal{W}_j = \text{blockdiag}\big(W_{1j}, W_{2j}, \dots, W_{Pj}\big), \qquad W_{pj} \in \mathbb{R}^{1\times q_p},\quad q = \sum_{p=1}^P q_p$$

Each coefficient dimension $p$ may have its own set (and number) of
level-2 predictors — $q_p$ need not be equal across $p$, and $q_p = 1$
(intercept only) is the natural default when a coefficient has no
level-2 predictors.

## 3. Internal storage

| Quantity | Storage | Shape |
|---|---|---|
| Stacked likelihood design | `Dmat` | $n \times P$ |
| Group index | `group` | length-$n$ vector, $g_i = j$ |
| Level-2 designs, one per coefficient | `Wlist = list(W_1, ..., W_P)` | `Wlist[[p]]` is $J \times q_p$ |
| Dimension bookkeeping | `qvec = sapply(Wlist, ncol)`, `q = sum(qvec)` | — |

Recovery rules:

- $D_j$ = `Dmat[group == j, ]`
- $\mathcal{W}_j$ = `blockdiag(W_1[j, ], W_2[j, ], ..., W_P[j, ])`, assembled on demand (row $j$ from each list element)

`Dmat` is stored row-stacked, **not** as the block-diagonal
$n \times PJ$ object `lme4` implicitly uses for $Z$; the block-diagonal
structure is recovered implicitly via `group` rather than materialized.

`Wlist`'s heterogeneous column widths ($q_p$ varying by $p$) are the
natural fix for the raggedness of level-2 predictor sets across
coefficient dimensions — no padding needed, unlike a 3-index array.

## 4. Full conditionals (Gibbs)

**$\beta_j$** (conjugate normal, per group):

$$\beta_j \mid \cdot \ \sim\ N(m_j, V_j)$$
$$V_j = \Big(\tfrac{1}{\sigma_j^2}D_j'D_j + \Psi^{-1}\Big)^{-1}, \qquad m_j = V_j\Big(\tfrac{1}{\sigma_j^2}D_j'y_j + \Psi^{-1}\mathcal{W}_j\gamma\Big)$$

**$\gamma$** (GLS across groups):

$$\hat\gamma = \Big(\sum_j \mathcal{W}_j'\Psi^{-1}\mathcal{W}_j\Big)^{-1}\sum_j \mathcal{W}_j'\Psi^{-1}\beta_j$$

Because $\mathcal{W}_j$ is block-diagonal, $\mathcal{W}_j'\Psi^{-1}\mathcal{W}_j$
and $\mathcal{W}_j'\Psi^{-1}\beta_j$ decompose into $P\times P$ /
$P\times1$ sub-blocks built directly from `Wlist[[p]][j, ]` and the
$(p,p')$ entries of $\Psi^{-1}$ — implement this way rather than
assembling the dense $\mathcal{W}_j$ each iteration.

**$\Psi$**: updated from $u_j = \beta_j - \mathcal{W}_j\gamma$ (derived,
not sampled) via the usual inverse-Wishart / truncated-gamma conjugate
update.

## 5. Correspondence to `lme4` / Laird-Ware

Marginalizing (substituting the stage-2 equation into stage 1) recovers
the standard `lme4` form:

$$y_j = D_j(\mathcal{W}_j\gamma + u_j) + \varepsilon_j = \underbrace{(D_j\mathcal{W}_j)}_{X_j}\gamma + \underbrace{D_j}_{Z_j}u_j + \varepsilon_j$$

| `lmerb` object | `lme4` / Laird-Ware object | Relationship |
|---|---|---|
| $D_j$ | — | consumed into **both** $X_j$ and $Z_j$ below |
| $\mathcal{W}_j$ | — | consumed into $X_j$ |
| $D_j\mathcal{W}_j$ | $X_j$ | fixed-effects design block, group $j$ |
| $D_j$ | $Z_j$ | random-effects design block, group $j$ |
| $\gamma$ | $\beta$ | fixed effects |
| $u_j = \beta_j - \mathcal{W}_j\gamma$ | $b_j$ | random effects (mean zero) |
| $\Psi$ | per-group block of $G$ | $G = I_J \otimes \Psi$ |
| $\beta_j$ | *(no direct analog)* | closest is `lme4::coef()` output, $\hat\beta + \hat b_j$, computed post hoc |

**Note the asymmetry**: $D_j$ has no single `lme4` counterpart — it is
the common ancestor of both $X_j$ (after right-multiplication by
$\mathcal{W}_j$) and $Z_j$ (used as-is). This reflects a structural
difference, not just relabeling: `lmerb`'s parameterization retains the
full, non-centered $\beta_j$ in the likelihood, whereas `lme4` works
with the marginal model in which only the mean-zero deviation $b_j$
and the pooled fixed effect $\beta$ survive.

## 6. Equivalence test (for package test suite)

When $\mathcal{W}_j$ is intercept-only for every coefficient dimension
$p$ (i.e., no level-2 predictors, $q_p = 1\ \forall p$), the model
reduces to a standard random-intercept/random-slope `lmer` model. In
this special case:

- $\hat\gamma$ (posterior mean) should recover `lme4`'s fixed-effect
  estimates
- $\hat\beta_j - \hat\gamma$ should track `lme4::ranef()`
- $\hat\beta_j$ should track `lme4::coef()`

This provides a natural automated equivalence check against `lmer` fits.

---

*Reference: Lindley, D.V. and Smith, A.F.M. (1972). "Bayes Estimates
for the Linear Model." Journal of the Royal Statistical Society,
Series B, 34(1), 1-41.*
