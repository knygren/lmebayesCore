# Building likelihood subgradient densities for the marginal population density

Draft design note for **lmebayesCore** / **lmebayes** two-block Gibbs and iid
sampling. It connects the JASA likelihood-subgradient-density (LSD) theory
(Nygren & Nygren, 2006) to a **marginal target on population parameters**
\(\gamma\), and records cost trade-offs the user flagged: expensive build,
cheap proposal, expensive accept–reject.

**Primary references (on this machine).**

| Source | Location |
|---|---|
| Nygren & Nygren (2006), *Likelihood Subgradient Densities*, JASA 101(475):1144–1156 | User reference folder for **lmebayes**; also `inst/REFERENCES.bib` (`Nygren2006`) in **glmbayes** / **lmebayes** |
| Implementation walkthrough | **glmbayes** `vignette("Chapter-A05")` — standard form, `EnvelopeBuild`, grid types |
| Envelope map / PLSD details | **glmbayes** `vignette("Chapter-A08")` |
| Block 1 accept–reject in GLMMs | **lmebayes** `vignette("Chapter-09")` (Poisson GLMM example) |
| Code | **glmbayesCore** `EnvelopeBuild()`, `EnvelopeSize()`; **glmbayes** `.rNormalGLM_cpp` pipeline |

**Related notes in this repo:** `LOGIT_MARGINAL_INTEGRATE_GAMMA.md` (marginal on
\(\beta\) after \(\int\gamma\)), `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md`
(Prékopa + per-group marginals; majorization on \(\gamma\) for tail budgets),
`RESTRICTED_GIBBS_MINORIZATION_TV_v2.md` (Prékopa on bridge / integrated terms),
`BLOCK_GIBBS_ERGODICITY_ING.md` §16 (generalized LSD when \(-\log f\) is not
convex on a **piece** — different issue from marginal log-concavity).

---

## 1. What problem we are trying to solve

In the hierarchical GLMM notation (`notation.md`, `LOGIT_MARGINAL_INTEGRATE_GAMMA.md`),

\[
\beta_j \mid \gamma \sim N(\mathcal W_j \gamma,\,\Psi),
\qquad
\gamma \sim N(\mu_0,\,\Lambda_\gamma^{-1}),
\qquad
y_j \mid \beta_j \sim \text{GLM family (log-concave)}.
\]

**Population block (Block 2)** targets the marginal

\[
\pi_\gamma(\gamma \mid y)
\;\propto\;
p(\gamma)\,
\underbrace{
\int \exp\!\Bigl(\sum_j \ell_j(\beta_j)\Bigr)\,
\exp\!\Bigl(-\tfrac12(\beta-\mathcal W\gamma)^\top P_b^{\mathrm{stack}}(\beta-\mathcal W\gamma)\Bigr)\,d\beta
}_{M(\gamma)\ \text{(integrated likelihood)}}.
\]

**Likelihood-subgradient theory** (JASA paper) applies when the target log-density
decomposes as

\[
\boxed{
-\log \pi(\theta \mid y)
= \underbrace{\tfrac12(\theta-\mu)^\top P(\theta-\mu)}_{\text{multivariate normal prior term}}
+ \underbrace{(-\log L(\theta))}_{\text{log-concave ``likelihood'' term}},
}
\]

with **log-concave** \(L\). Then one can build an **envelope mixture** of
**restricted multivariate normal** densities (the LSD components) that
**majorizes** \(\pi\) and supports exact accept–reject sampling.

**Today’s production path** uses this decomposition for **Block 1** draws of
\(\beta \mid \gamma, y\) (conditional on fixed \(\gamma\)), not for
\(\pi_\gamma(\gamma\mid y)\). This note asks whether the **same machinery** can
be centered on the **marginal population density** \(\pi_\gamma\), and what that
costs.

---

## 2. LSD recap in the notation the code uses

After centering the prior mean to zero and **standardizing** (`glmb_Standardize_Model()`;
Chapter-A05 §5), the working target is

\[
-\log \pi(\theta^\star \mid y)
= \tfrac12 \|\theta\|_2^2 + \mathrm{NLL}(\theta),
\qquad
\theta := \theta^\star - \mu,
\]

with **identity** prior precision and **diagonal** Hessian of \(-\log\pi\) at the
mode (Definition 3 / Remarks 11–15, Nygren & Nygren 2006).

Let \(\theta^\dagger\) be the posterior mode. At a **tangency point**
\(\bar\theta_j\) (one of \(\theta^\dagger_i \pm \omega_i\) per coordinate in the
\(3^p\) grid), let \(c_j := \nabla \mathrm{NLL}(\bar\theta_j)\) (subgradient of
the negative log-likelihood piece; in the interior of a concave piece this is the
ordinary gradient).

A **single LSD component** (Chapter-A05 / `EnvelopeBuild` roxygen) is

\[
f_j(\theta)
\;=\;
\frac{
\exp\!\Bigl(
-\tfrac12 \theta^\top A \theta + c_j^\top(\theta-\bar\theta_j)
\Bigr)
}{
(2\pi)^{p/2}\,|A|^{-1/2}\,\mathrm{MGF}_A(c_j)
},
\qquad
\mathrm{MGF}_A(c_j) = \exp\!\bigl(\tfrac12 c_j^\top A^{-1} c_j\bigr),
\]

with diagonal **data precision** \(A\) in standard form. Each \(f_j\) is a
**tilted, renormalized** multivariate normal that **dominates** the posterior on
its restriction region.

The **envelope** is a mixture over grid faces \(j=1,\ldots,K\) (\(K\le 3^p\)):

\[
\boxed{
g(\theta) = \sum_{j=1}^K p_j\, f_j(\theta),
\qquad
\pi(\theta \mid y) \le g(\theta)\ \ \text{(majorization)}.
}
\]

Face probabilities **PLSD** (\(p_j\)) come from log-CDF differences on each
face’s box (`loglt`, `logrt`, `cbars` in `EnvelopeBuild`). Sampling: pick
\(J\sim \mathrm{PLSD}\), draw \(\theta\sim f_J\), accept with probability
\(\pi(\theta)/\bigl(M_J g(\theta)\bigr)\) where \(M_J\) is the face-specific
slack (legacy `rNormalGLM_std` / `lg_prob_factor` algebra).

**Optimally centered** means: choose \(\theta^\dagger = \arg\max \pi(\theta\mid y)\)
and build the grid around that mode. All tangency widths \(\omega_i\) and
subgradients are computed relative to that center.

---

## 3. Log-concavity of \(\pi_\gamma\) (Prékopa) and the LSD template

Write

\[
\log \pi_\gamma(\gamma \mid y)
= \log p(\gamma) + \log M(\gamma) + \text{const},
\qquad
\log p(\gamma) = -\tfrac12(\gamma-\mu_0)^\top \Lambda_\gamma(\gamma-\mu_0),
\]

\[
M(\gamma)
=
\int
\exp\!\Bigl(\sum_j \ell_j(\beta_j)\Bigr)\,
\exp\!\Bigl(-\tfrac12(\beta-\mathcal W\gamma)^\top P_b^{\mathrm{stack}}(\beta-\mathcal W\gamma)\Bigr)\,d\beta .
\]

### 3.1 Prékopa: \(\pi_\gamma\) is log-concave under standard GLMM assumptions

**Prékopa’s theorem** (marginalization preserves log-concavity): if
\(F(\beta,\gamma)\) is log-concave in \((\beta,\gamma)\), then
\(\int F(\beta,\gamma)\,d\beta\) is log-concave in \(\gamma\).

For canonical GLM stage-1 likelihoods (logit, Poisson, …), each \(\ell_j(\beta_j)\)
is **concave**. Define the \(\beta\)-block integrand (before the explicit prior on
\(\gamma\)):

\[
L(\beta,\gamma)
:=
\sum_j \ell_j(\beta_j)
-\tfrac12\sum_j(\beta_j-\mathcal W_j\gamma)^\top P_b(\beta_j-\mathcal W_j\gamma).
\]

Then \(L\) is **jointly concave** in \((\beta,\gamma)\): a sum of concave functions
in \(\beta_j\), plus a concave quadratic in the affine map \(\beta_j-\mathcal W_j\gamma\).
The same argument appears in `RESTRICTED_GIBBS_MINORIZATION_TV_v2.md` for the bridge
law \(q(\gamma'\mid\gamma)\) and in `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` for
per-group marginals \(\widetilde\pi_j\) after \(\int\beta_{-j}\).

Therefore:

\[
\boxed{
M(\gamma) = \int e^{L(\beta,\gamma)}\,d\beta
\ \text{is log-concave in }\gamma,
}
\]

and since \(\log p(\gamma)\) is concave,

\[
\boxed{
\pi_\gamma(\gamma\mid y)\ \text{is log-concave in }\gamma
\quad\text{(under log-concave GLM likelihoods + Gaussian stage-2 prior).}
}
\]

**Corollary (JASA decomposition).** Define the integrated “likelihood” piece
\(\mathrm{NLL}_M(\gamma) := -\log M(\gamma)\). Then

\[
\boxed{
-\log \pi_\gamma(\gamma\mid y)
=
\underbrace{\tfrac12(\gamma-\mu_0)^\top \Lambda_\gamma(\gamma-\mu_0)}_{\text{population prior}}
+
\underbrace{\mathrm{NLL}_M(\gamma)}_{\text{convex (integrated likelihood)}}
+ \text{const},
}
\]

so the **plain** Nygren & Nygren (2006) setup applies in principle: multivariate
normal prior term + log-concave remainder. No generalized chord envelope is
required for log-concavity itself.

**What Prékopa does *not* give:** a closed form for \(M(\gamma)\), \(\nabla M\), or
\(\nabla^2 M\). Log-concavity is an **existence** result for LSD theory, not a
**computational** shortcut.

### 3.2 Gaussian first stage (closure)

If every \(\ell_j\) is quadratic, \(M(\gamma)\) is proportional to a multivariate
normal density in \(\gamma\) and \(\mathrm{NLL}_M\) is quadratic. Envelopes are
tight; accept rates approach the \(2/\sqrt\pi\) limits in Theorem 2–3 of the paper.
Gradients of \(\mathrm{NLL}_M\) are explicit — this is the easy end of the same
Prékopa story.

### 3.3 Same fact on \(\widetilde\pi(\beta\mid y)\) after \(\int\gamma\)

Marginalizing \(\gamma\) instead of \(\beta\) gives
\(\widetilde\pi(\beta\mid y)\) in `LOGIT_MARGINAL_INTEGRATE_GAMMA.md` — again
log-concave by Prékopa (also visible directly as Schur quadratic prior plus
\(\sum_j \ell_j\)). Per-group marginals \(\widetilde\pi_j\) after \(\int\beta_{-j}\)
are log-concave for the same reason (`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` §10).

### 3.4 Where the obstruction actually is

| Question | Answer |
|---|---|
| Is \(\pi_\gamma\) log-concave? | **Yes** (Prékopa), under standard GLMM assumptions |
| Does plain LSD theory apply? | **Yes**, as MN prior + \(\mathrm{NLL}_M\) |
| Is \(\mathrm{NLL}_M\) explicit? | **No** (except Gaussian closure) |
| Is \(\nabla \mathrm{NLL}_M\) cheap? | **No** — needs \(\mathbb E[\cdot\mid\gamma,y]\) over \(\beta\) |
| Is building the envelope cheap? | **No** — many face evaluations of \(f_2\)/\(f_3\) at grid points |

So the open problem is **not** log-concavity; it is **evaluating** (and
differentiating) the integrated object that Prékopa guarantees is well-behaved.

### 3.5 Functional form: integral, then conditional expectations

Split the population posterior into an **explicit prior** and an **integrated
likelihood** (same factorization as §1). With \(P_b=\Psi^{-1}\) and stacked
\(\mathcal W=\mathrm{blockdiag}(\mathcal W_1,\ldots,\mathcal W_J)\),

\[
\boxed{
M(\gamma)
=
\int_{\mathbb R^{JP}}
\exp\!\Bigl(
L(\beta,\gamma)
\Bigr)\,d\beta,
}
\]

\[
\boxed{
L(\beta,\gamma)
=
\sum_{j=1}^{J}\ell_j(\beta_j)
-\tfrac12(\beta-\mathcal W\gamma)^\top P_b(\beta-\mathcal W\gamma).
}
\]

Stage-1 log-likelihood \(\ell_j\) depends **only on \(\beta_j\)**; all
\(\gamma\)-dependence in \(L\) is through the **Gaussian coupling**
\((\beta-\mathcal W\gamma)\). The population prior enters separately:

\[
\log \pi_\gamma(\gamma\mid y)
=
-\tfrac12(\gamma-\mu_0)^\top \Lambda_\gamma(\gamma-\mu_0)
+
\underbrace{\log M(\gamma)}_{\text{no closed form for logit/Poisson}}
+ \text{const}.
\]

For LSD, the **“log-likelihood” piece** in the JASA decomposition is

\[
\boxed{
\mathrm{NLL}_M(\gamma) := -\log M(\gamma).
}
\]

There is **no** simpler closed-form \(\mathrm{NLL}_M(\gamma)\) for general GLMs
— only the **integral definition** above (or equivalents below).

#### Conditional posterior that defines the expectations

For each fixed \(\gamma\), define the **\(\beta\)-block conditional**

\[
\boxed{
\pi(\beta\mid\gamma,y)
=
\frac{
\exp\bigl(L(\beta,\gamma)\bigr)
}{
Z(\gamma)
},
\qquad
Z(\gamma) := M(\gamma).
}
\]

This is exactly the target Block 1 would sample from in `glmerb()` (same
\(L(\beta,\gamma)\), different outer loop). All population derivatives of
\(\log M\) are **expectations under this conditional**.

#### First derivative (gradient)

Since \(\nabla_\gamma \ell_j(\beta_j)=0\),

\[
\nabla_\gamma L(\beta,\gamma)
=
\mathcal W^\top P_b(\beta-\mathcal W\gamma).
\]

Therefore

\[
\boxed{
\nabla_\gamma \log M(\gamma)
=
\mathbb E_{\beta\mid\gamma,y}\bigl[
\mathcal W^\top P_b(\beta-\mathcal W\gamma)
\bigr]
=
\mathcal W^\top P_b\Bigl(
\mathbb E_{\beta\mid\gamma,y}[\beta] - \mathcal W\gamma
\Bigr).
}
\]

Per group \(j\), with \(P_b\) block-diagonal across groups,

\[
\bigl[\nabla_\gamma \log M(\gamma)\bigr]_k
=
\sum_{j=1}^{J}
\mathcal W_{jk}^\top P_b
\Bigl(
\mathbb E_{\beta\mid\gamma,y}[\beta_j] - \mathcal W_j\gamma
\Bigr).
\]

So **yes**: at each \(\gamma\), the gradient of the integrated piece is a
**linear functional of the conditional mean** \(\mathbb E[\beta_j\mid\gamma,y]\)
(and hence of conditional means of predictors \(\eta_{j,i}=d_{j,i}^\top\beta_j\)
if you expand through \(\mathcal W_j\)).

Equivalently, with the **prior score in \(\beta\)** at fixed \(\gamma\),

\[
\nabla_\gamma \log M(\gamma)
=
\mathbb E_{\beta\mid\gamma,y}\bigl[
\nabla_\gamma \log p(\beta\mid\gamma)
\bigr],
\]

because \(\log p(\beta\mid\gamma)= -\tfrac12(\beta-\mathcal W\gamma)^\top P_b(\beta-\mathcal W\gamma)+\text{const}\).

The **full** population gradient for \(-\log\pi_\gamma\) is

\[
\nabla_\gamma\bigl[-\log\pi_\gamma(\gamma\mid y)\bigr]
=
\Lambda_\gamma(\gamma-\mu_0)
-
\nabla_\gamma \log M(\gamma).
\]

#### Second derivative (Hessian piece)

Differentiate again (log-derivative formula). Since
\(\nabla^2_\gamma L(\beta,\gamma) = -\mathcal W^\top P_b \mathcal W\) is
**constant in \(\beta\)**,

\[
\boxed{
\nabla^2_\gamma \log M(\gamma)
=
-\mathcal W^\top P_b \mathcal W
+
\mathrm{Cov}_{\beta\mid\gamma,y}\!\Bigl(
\mathcal W^\top P_b(\beta-\mathcal W\gamma)
\Bigr).
}
\]

So the Hessian of the integrated piece needs the **conditional covariance** of
\(\beta\) (or of the prior score \(\mathcal W^\top P_b(\beta-\mathcal W\gamma)\)),
not just the mean. Envelope construction (`glmb_Standardize_Model`, face
subgradients) uses **Hessian information at the mode**, so each outer iteration
needs second-order conditional summaries (or a numerical Hessian built from many
gradient calls).

#### What LSD calls \(\mathrm{NLL}_M\) at a grid point

At tangency \(\bar\gamma_j\), the subgradient vector for the **integrated**
piece is

\[
c_j := \nabla \mathrm{NLL}_M(\bar\gamma_j)
=
-\mathcal W^\top P_b\Bigl(
\mathbb E_{\beta\mid\bar\gamma_j,y}[\beta] - \mathcal W\bar\gamma_j
\Bigr).
\]

Each of the \(K\) envelope faces needs its **own** conditional expectation under
\(\pi(\beta\mid\bar\gamma_j,y)\) — potentially a **fresh** Block-1 accuracy
problem per face.

#### Laplace functional form (approximate, not exact)

A common **approximate** functional form (Laplace at the conditional mode
\(\hat\beta(\gamma)=\arg\min_\beta L(\beta,\gamma)\)) is

\[
\log M(\gamma)
\approx
L\bigl(\hat\beta(\gamma),\gamma\bigr)
+
\tfrac{JP}{2}\log(2\pi)
-\tfrac12\log\det H_\beta(\gamma),
\]

\[
H_\beta(\gamma)
=
\nabla^2_\beta L(\beta,\gamma)\Big|_{\hat\beta(\gamma)}
=
\mathrm{blockdiag}\bigl(-\nabla^2\ell_j(\hat\beta_j)\bigr) + P_b.
\]

This is **explicit in \(\gamma\)** only through \(\hat\beta(\gamma)\) and
\(\det H_\beta(\gamma)\) — still an implicit function, but no integral. It is
**not** what Prékopa or exact LSD require; it is a separate approximation route
(cf. `LAPLACE_PROFILE_GROUP_MARGINAL_BOUND.md` on the \(\beta\)-marginal side).

**Summary.** The exact “log-likelihood” for the population block is
\(\mathrm{NLL}_M(\gamma)=-\log\int e^{L(\beta,\gamma)}d\beta\). Its derivatives
are **conditional expectations and covariances** of \(\beta\mid\gamma,y\) — the
same objects Block 1 sampling is designed to approximate.

### 3.6 Nygren majorization (correct normalization) and estimating \(M(\gamma)\)

Fix \(\gamma\). Write the **unnormalized** conditional kernel and normalized
posterior as in §3.5:

\[
f(\beta) := e^{L(\beta,\gamma)} ,
\qquad
M(\gamma) = \int f(\beta)\,d\beta ,
\qquad
\pi(\beta\mid\gamma,y) = \frac{f(\beta)}{M(\gamma)}.
\]

**Important:** The JASA envelope does **not** majorize \(f\) directly. It
majorizes the **normalized** density (Nygren & Nygren 2006, Theorem 1;
`glmbayes` `vignette("Chapter-A08")` §2.3–§2.6).

Let \(g_{\bar\theta}(\beta)\) be a **likelihood-subgradient density** with
\(\int g_{\bar\theta}\,d\beta = 1\) (restricted, renormalized on face \(A\), then
mixed with PLSD weights). Let \(\tilde a(\bar\theta)\) denote the **restricted
constant** from Remark 5 — a theoretical slack in the majorization below, **not**
a quantity the shipped sampler evaluates explicitly.

#### Correct dominance relation

On each face, for \(\beta\in A\),

\[
\boxed{
\pi(\beta\mid\gamma,y)
\;=\;
\tilde a(\bar\theta)\,
h_{\bar\theta}(\beta)\,
\tilde q_{\bar\theta}(\beta),
\qquad
0\le h_{\bar\theta}(\beta)\le 1,
}
\]

and (Theorem 1, before restriction)

\[
\boxed{
\pi(\beta\mid\gamma,y)
\;\le\;
a(\bar\theta)\, q_{\bar\theta}(\beta).
}
\]

Equivalently in **unnormalized** form, with \(f(y)\) in the paper replaced by
\(M(\gamma)\) for this conditional block:

\[
\boxed{
\frac{f(\beta)}{M(\gamma)}
\;\le\;
\tilde a(\bar\theta)\, g(\beta),
\qquad
\int g(\beta)\,d\beta = 1.
}
\]

This matches \(f(\beta)/M(\gamma) \le \tilde a\, g(\beta)\) with \(\int g=1\).

#### Sampler implementation: \(\tilde a\) need not be known

**Key point:** Correct accept–reject for \(\pi(\beta\mid\gamma,y)\) does **not**
require knowing \(\tilde a\) (nor \(M(\gamma)\)) as explicit scalars at run time.
`EnvelopeBuild` / `setlogP_C2` fold the normalization into face tables
(`logP`, `LLconst`, PLSD); `rNormalGLM.cpp` implements

\[
\text{accept iff}\quad
\mathrm{LLconst}[J] + c_J^\top\beta - \mathrm{NLL}(\beta) \ge \log U,
\]

which is algebraically equivalent to accepting with probability
\(\pi(\beta\mid\gamma,y)/(\tilde a(\bar\theta_J)\, g(\beta))\) **without ever
computing** \(\tilde a\) or \(M(\gamma)\). This is why Block‑1 sampling ships as
a self-contained envelope draw: build once, propose from \(g\), apply the stored
accept test.

In practice \(\tilde a\) is **not evaluated at run time** as a user-visible scalar;
after trials, \(\widehat p_{\mathrm{acc}}\approx 1/\tilde a_{\mathrm{eff}}\) — envelope
tightness learned **after the fact**, not required as sampler input. For
**certification**, the JASA paper gives **explicit** face constants (Theorem 1 /
Remark 6) and \(\tilde a\le (2/\sqrt\pi)^k\) in multivariate Normal standard form
(Theorems 2–3; §7.3). In **practice** the same constants apply to log-concave
non-Gaussian GLMs (logit, Poisson, …) with good empirical behaviour; see §7.3.

#### What acceptance actually estimates

The shipped sampler draws \(\beta\sim g\) (mixture of restricted LSD components)
and accepts with probability (Chapter-A08 §2.3; `rNormalGLM.cpp`)

\[
\frac{\pi(\beta\mid\gamma,y)}{\tilde a(\bar\theta_J)\, g(\beta)}
\;=\;
\frac{f(\beta)}{M(\gamma)\,\tilde a(\bar\theta_J)\, g(\beta)},
\]

using the accept test \(\mathrm{LLconst}[J] + c_J^\top\beta - \mathrm{NLL}(\beta)
\ge \log U\). Therefore

\[
\boxed{
\mathbb E_{\beta\sim g}[\text{accept}]
\;=\;
\int \frac{f(\beta)}{M(\gamma)\,\tilde a\, g(\beta)}\, g(\beta)\,d\beta
\;=\;
\frac{1}{\tilde a_{\mathrm{eff}}},
}
\]

where \(\tilde a_{\mathrm{eff}}\) is the appropriate face/mixture slack (Theorem 1:
**acceptance probability \(= 1/a(\bar\theta)\)**). The expected rate estimates
**\(1/a\)** (or \(1/\tilde a\) on the restricted face), **not** \(M(\gamma)\) by
itself — but see below.

An earlier draft incorrectly used \(f/g\) and claimed \(\mathbb E[\alpha]=M(\gamma)\);
that omitted the **\(1/M(\gamma)\)** factor in the correct majorization.

#### From estimated \(a(\bar\theta)\) to \(M(\gamma)\) (Theorem 1)

Acceptance rate alone is **not** \(M(\gamma)\). Once we **also** have an estimate
\(\widehat a(\bar\theta)\) of the face constant, Theorem 1 (Chapter-A08 §2.3) gives
\(M(\gamma)\) immediately. At tangency \(\bar\theta\) on a face,

\[
a(\bar\theta)
\;=\;
\frac{
g(\bar\theta)\,\mathrm{MGF}\bigl(-c(\bar\theta)\bigr)
}{
M(\gamma)\,
\exp\!\bigl(-c(\bar\theta)^\top\bar\theta\bigr)
},
\]

with \(g(\bar\theta)\) the LSD mixture density at \(\bar\theta\), \(c(\bar\theta)\) the
face subgradient from the build, and \(\mathrm{MGF}(-c)=\exp(\tfrac12\|c\|^2)\) on the
standardized grid. Rearranging,

\[
\boxed{
\widehat M(\gamma)
\;=\;
\frac{
g(\bar\theta)\,\mathrm{MGF}\bigl(-c(\bar\theta)\bigr)
}{
\widehat a(\bar\theta)\,
\exp\!\bigl(-c(\bar\theta)^\top\bar\theta\bigr)
}.
}
\]

**Estimating \(\widehat a\).** The sampler does not need \(a\) at run time, but after
\(N_{\mathrm{in}}\) inner trials,

\[
\widehat a(\bar\theta)
\;\approx\;
\frac{1}{\widehat p_{\mathrm{acc}}},
\qquad
\widehat p_{\mathrm{acc}}
=
\frac{1}{N_{\mathrm{in}}}\sum_{k=1}^{N_{\mathrm{in}}} \mathbf 1\{\text{accept}_k\},
\]

for the **mixture** slack \(\tilde a_{\mathrm{eff}}\) (overall acceptance rate). For a
**specific face** \(\bar\theta_j\), restrict to trials with \(J_k=j\) (or run face‑conditional
acceptance) to estimate \(\widehat a(\bar\theta_j)\). Use the **mode face**
\(\bar\theta=\beta^\dagger\) when \(h_{\bar\theta}(\bar\theta)=1\).

**What the build supplies vs what MC supplies.**

| Ingredient | Source |
|---|---|
| \(g(\bar\theta)\), \(c(\bar\theta)\), \(\bar\theta\), PLSD | Step 1 `EnvelopeBuild` (known from tables) |
| \(a(\bar\theta)\) | **Estimated** post hoc: \(\widehat a \approx 1/\widehat p_{\mathrm{acc}}\) |
| \(M(\gamma)\) | **Theorem 1 inversion** above with \(\widehat a\) |

So the nested pipeline can recover \(\widehat M(\bar\gamma_j)\) at each outer tangency
**without** ever evaluating \(a\) or \(M\) during inner sampling: build once, run trials,
estimate \(\widehat a\), invert Theorem 1. The inner sampler itself still never needs
\(\tilde a\) or \(M\) explicitly — only the **post-processing** step does.

**Payoff (population LSD).** This is not an end in itself. Theorem 1 + \(\widehat a\)
supplies an **approximate** value of \(\mathrm{NLL}_M(\bar\gamma_j)=-\log M(\bar\gamma_j)\)
at each outer tangency. Together with **approximate subgradients**
\(\widehat c_j \approx \nabla\mathrm{NLL}_M(\bar\gamma_j)\) from accepted inner
\(\beta\) (§3.5), these are exactly the **face ingredients** Chapter-A05 needs to
build a **population likelihood-subgradient density** \(g_\gamma(\gamma)\) on
\(\pi_\gamma(\gamma\mid y)\). The outer envelope is therefore **derivable** from
nested Block‑1 machinery plus MC — close to the true population LSD when inner
draw budgets are large enough that \(\widehat M\) and \(\widehat{\nabla\log M}\) are
accurate. The outer sampler, like Block 1, still need not know its own slack constant
explicitly once `LLconst` / PLSD tables are built.

#### Monte Carlo role for \(\tilde a\) and \(M\)

Run \(N_{\mathrm{in}}\) inner envelope trials to estimate **acceptance rate**
\(\widehat p_{\mathrm{acc}}\). Under a correct build,

\[
\widehat p_{\mathrm{acc}} \approx \frac{1}{\tilde a_{\mathrm{eff}}},
\qquad
\widehat a \approx \frac{1}{\widehat p_{\mathrm{acc}}},
\]

then Theorem 1 yields \(\widehat M(\gamma)\) as above. Use accepted draws (§3.5, §4.4)
for \(\nabla\log M(\gamma)\).

| Quantity | Role |
|---|---|
| Inner sampler | `LLconst` + PLSD — **no explicit \(\tilde a\) or \(M\)** needed |
| Acceptance rate | Estimates **\(1/a\)**; with build tables, feeds **\(\widehat M(\gamma)\)** via Theorem 1 |
| \(\nabla\log M(\gamma)\) | MC expectations from accepted \(\beta\) (§3.5) |
| **Population LSD** \(g_\gamma\) | Outer `EnvelopeBuild` from approximate \(\mathrm{NLL}_M\), \(\nabla\mathrm{NLL}_M\) at \(\gamma\) tangencies (§4) |

**Population gradient vs integral.** \(\widehat M\) comes from \(\widehat a\) + Theorem 1;
\(\nabla\log M\) comes from accepted \(\beta\) (§3.5). Both feed the **same inner run**
and together approximate the outer face data needed to **derive** \(g_\gamma\).

**Generalized LSD** (`BLOCK_GIBBS_ERGODICITY_ING.md` §16) remains relevant when
\(-\log f\) fails to be **convex on a restriction piece** (chord majorants on
annuli), not because the **global** marginal lacks log-concavity.

### 3.7 Targets compared

| Target | Log-concave? | Plain LSD? | Main difficulty |
|---|---|---|---|
| \(\pi(\beta\mid\gamma,y)\) | Yes | Yes | **None** (shipped Block 1) |
| \(\pi_\gamma(\gamma\mid y)\) | Yes (Prékopa) | Yes (in principle) | **Approximate** face \(\mathrm{NLL},\nabla\mathrm{NLL}\) via nested inner envelopes |
| \(\widetilde\pi(\beta\mid y)\) after \(\int\gamma\) | Yes (Prékopa) | Yes | Dimension \(Jp_{\mathrm{re}}\); explicit \(\nabla\) |

The rest of this note assumes the Prékopa/JASA decomposition. **§4** records the
proposed **nested-envelope pipeline**; later sections discuss cost and alternatives.

---

## 4. Proposed pipeline: nested \(\beta\)-envelopes → approximate population LSD

This section documents the **intended construction path** for a population LSD
envelope \(g_\gamma(\gamma)\) on \(\pi_\gamma(\gamma\mid y)\), using existing
Block‑1 machinery at each inner \(\gamma\) tangency. The **goal** is not merely to
estimate \(M(\gamma)\) as a scalar, but to obtain the **face-level negative log-density
and subgradients** that Chapter-A05 requires — both **approximated** from inner
envelope draws — and then run the **outer** `EnvelopeBuild` to produce a mixture
density from which population parameters can be sampled (exactly up to MC error in
the face estimates, and up to envelope slack as in Block 1).

**What each inner tangency supplies to the outer LSD.**

\[
-\log \pi_\gamma(\bar\gamma_j\mid y)
=
\underbrace{\tfrac12(\bar\gamma_j-\mu_0)^\top \Lambda_\gamma(\bar\gamma_j-\mu_0)}_{\text{exact prior}}
+
\underbrace{\widehat{\mathrm{NLL}}_M(\bar\gamma_j)}_{\text{from }\widehat a + \text{Theorem 1}}
+ \text{const},
\]

\[
\widehat c_j
\approx
\nabla\mathrm{NLL}_\gamma(\bar\gamma_j)
=
\underbrace{\Lambda_\gamma(\bar\gamma_j-\mu_0)}_{\text{exact prior gradient}}
+
\underbrace{(-\widehat{\nabla\log M(\bar\gamma_j)})}_{\text{from accepted inner }\beta}.
\]

With \((\bar\gamma_j,\,\widehat c_j,\,-\log\pi_\gamma(\bar\gamma_j))\) on the outer
grid, the shipped envelope algebra (standard form, restricted tilted MVNs, PLSD)
yields an **approximate population LSD** \(g_\gamma\). Refining inner draw counts
makes \(g_\gamma\) approach the LSD one would obtain from exact \(\mathrm{NLL}_M\)
and \(\nabla\mathrm{NLL}_M\).

**Outer target (population).** By §3.1–§3.5,

\[
-\log \pi_\gamma(\gamma\mid y)
=
\tfrac12(\gamma-\mu_0)^\top \Lambda_\gamma(\gamma-\mu_0)
+ \mathrm{NLL}_M(\gamma) + \text{const},
\qquad
\mathrm{NLL}_M(\gamma) = -\log M(\gamma).
\]

**Inner target (coefficients, \(\gamma\) fixed).** At each candidate \(\gamma\),

\[
\pi(\beta\mid\gamma,y) \propto e^{L(\beta,\gamma)},
\qquad
M(\gamma) = \int e^{L(\beta,\gamma)}\,d\beta.
\]

The outer envelope is built in the same pattern as Chapter-A05: mode
\(\gamma^\dagger\), standard form, tangency grid \(\bar\gamma_1,\ldots,\bar\gamma_K\),
face subgradients \(c_j\), then `EnvelopeBuild` / PLSD on \(\gamma\)-space. Steps
1–3 below **populate** those face tables with MC approximations; step 4 **assembles**
\(g_\gamma\).

### 4.1 Overview (four steps per outer tangency point)

For each outer tangency point \(\bar\gamma_j\) on the population grid:

| Step | Action | Output |
|---|---|---|
| **1** | Build Block‑1 LSD envelope(s) for \(\pi(\beta\mid\bar\gamma_j,y)\) | Inner `Env_out`: `cbars`, `PLSD`, `LLconst`, `logP`, … |
| **2** | Inner trials → \(\widehat a\approx 1/\widehat p_{\mathrm{acc}}\) → Theorem 1 → \(\widehat M(\bar\gamma_j)\) | \(\log M\) at tangency (§3.6); sampler needs neither \(a\) nor \(M\) |
| **3** | Use inner accepted \(\beta\) draws for conditional expectations | \(\widehat{\nabla\log M(\bar\gamma_j)}\), optional \(\widehat{\nabla^2\log M}\) |
| **4** | Outer `EnvelopeBuild` on \(\gamma\) from stored face data | **Approximate population LSD** \(g_\gamma\); sample \(\gamma\) as in Block 1 |

After **all** outer tangencies are processed, run the **outer** envelope build on
\(\gamma\) (standard form at \(\gamma^\dagger\), PLSD on population faces) —
same algebra as `EnvelopeBuild`, but `f2`/`f3` call the integrated
\(\mathrm{NLL}_M\) oracle defined by steps 1–3.

```text
  γ†  = mode of π_γ(γ | y)
        │
        ▼
  For each outer tangency γ̄_j  (j = 1 … K_γ)
        │
        ├─► (1) EnvelopeBuild( β | γ̄_j, y )     ← shipped Block 1
        │
        ├─► (2) N_in inner trials:  p̂_acc → â ≈ 1/p̂_acc
        │         Theorem 1 → M̂(γ̄_j)  (build supplies g(θ̄), c, MGF)
        │
        ├─► (3) Accepted β (or importance-weighted trials):
        │         ∇̂ log M(γ̄_j) = (1/n_acc) Σ W'P(β - W γ̄_j)
        │         c_j = -∇̂ log M(γ̄_j)   [outer subgradient]
        │
        └─► store (γ̄_j, c_j, log M̂_j) for outer face j
        │
        ▼
  EnvelopeBuild( γ | y )  →  approximate population LSD g_γ
        │
        ▼
  Sample γ ~ g_γ, accept w.r.t. π_γ  (Block 1 pattern on population block)
```

### 4.2 Step 1 — Inner envelope at \(\bar\gamma_j\)

Fix outer tangency \(\bar\gamma_j\). Construct the **conditional** LSD envelope
for \(\pi(\beta\mid\bar\gamma_j,y)\) using the shipped pipeline:

1. `glmb_Standardize_Model()` on the Block‑1 posterior at \(\gamma=\bar\gamma_j\).
2. `EnvelopeOpt` / `Gridtype` for inner grid size.
3. `EnvelopeBuild()` → inner envelope \(g_\beta(\beta)\) with `PLSD`, `LLconst`,
   `cbars_β`, `logP`, …

**Scope of “each \(\beta_j\)”.** Production `glmerb()` builds **one** envelope
for the **stacked** coefficient vector
\(\beta=(\beta_1^\top,\ldots,\beta_J^\top)^\top\) at fixed \(\gamma\), because
\(L(\beta,\gamma)\) couples groups only through the prior
\((\beta-\mathcal W\gamma)\) while \(\ell_j\) factorizes over \(j\). That single
inner build is the default in this pipeline.

Optionally (large \(J\), blockwise samplers in
`BLOCK_ING_RINDEPNORMALGAMMA_REG.md`), run **`BlockEnvelopeBuild`** per group:
one inner envelope per \(\beta_j\mid\bar\gamma_j,y\). Step 2–3 then run **per
block**; the outer population gradient still sums the \(\mathcal W_j\)-weighted
contributions in §3.5.

**Cost:** one full inner `EnvelopeBuild` **per outer tangency** \(\bar\gamma_j\)
(unless envelopes are interpolated/cached between nearby \(\gamma\) — not assumed here).

### 4.3 Step 2 — Estimate \(M(\bar\gamma_j)\) via \(\widehat a\) and Theorem 1

With inner unnormalized kernel \(f(\beta)=e^{L(\beta,\bar\gamma_j)}\), normalized
\(\pi(\beta\mid\bar\gamma_j,y)=f/M(\bar\gamma_j)\), and inner mixture envelope
\(g_\beta\) from step 1, the correct majorization is (§3.6)

\[
\frac{f(\beta)}{M(\bar\gamma_j)} \le \tilde a(\bar\theta)\, g(\beta),
\qquad \int g = 1.
\]

**Inner sampling does not need \(\tilde a\) or \(M\).** Step 1 suffices to run the
shipped accept test (`LLconst`, PLSD). Constants \(a\) and \(M\) enter only in
**post-processing** for the outer population envelope.

**Estimate \(a\), then invert Theorem 1.** Run \(N_{\mathrm{in}}\) inner envelope
trials at \(\bar\gamma_j\). Record acceptance indicators; optionally stratify by
drawn face \(J_k\). Then

\[
\widehat p_{\mathrm{acc}} = \frac{1}{N_{\mathrm{in}}}\sum_k \mathbf 1\{\text{accept}_k\},
\qquad
\widehat a(\bar\theta) \approx \frac{1}{\widehat p_{\mathrm{acc}}}
\]

(for mixture slack, or face‑conditional \(\widehat p_{\mathrm{acc},j}\) for
\(\widehat a(\bar\theta_j)\)). With build outputs \(g(\bar\theta)\), \(c(\bar\theta)\),
\(\bar\theta\) from step 1,

\[
\boxed{
\widehat M(\bar\gamma_j)
\;=\;
\frac{
g(\bar\theta)\,\mathrm{MGF}\bigl(-c(\bar\theta)\bigr)
}{
\widehat a(\bar\theta)\,
\exp\!\bigl(-c(\bar\theta)^\top\bar\theta\bigr)
}.
}
\]

Prefer the **mode face** \(\bar\theta=\beta^\dagger\) when \(h_{\bar\theta}(\bar\theta)=1\).
Store \(\log\widehat M(\bar\gamma_j)=-\mathrm{NLL}_M(\bar\gamma_j)\) for outer face
normalizers.

**Do not use \(f/g\).** An earlier draft used \(\widehat M=(1/N)\sum f/g\); that omitted
the \(1/M(\gamma)\) normalization. Acceptance rate estimates **\(1/a\)**, and Theorem 1
**links** \(\widehat a\) to \(\widehat M\) — not a direct identity \(\widehat p_{\mathrm{acc}}=M\).

**Sample size.** Propagate MC error in \(\widehat a\) into \(\widehat M\) (delta method
on \(\log \widehat M\)). Target \(\mathrm{SE}(\widehat p_{\mathrm{acc}})\le \tau_p\) so
\(\widehat a\) is stable; the same inner run can supply accepted \(\beta\) for step 3.

### 4.4 Step 3 — Gradients from inner accepted draws

Outer subgradient at tangency \(\bar\gamma_j\) (§3.5):

\[
c_j
:=
\nabla \mathrm{NLL}_M(\bar\gamma_j)
=
-\nabla_\gamma \log M(\gamma)\Big|_{\bar\gamma_j}
=
-\mathcal W^\top P_b\Bigl(
\mathbb E_{\beta\mid\bar\gamma_j,y}[\beta] - \mathcal W\bar\gamma_j
\Bigr).
\]

**Primary estimator (user path):** run inner accept–reject until \(n_{\mathrm{acc}}\)
**accepted** draws \(\beta^{(1)},\ldots,\beta^{(n_{\mathrm{acc}})}\sim
\pi(\beta\mid\bar\gamma_j,y)\), then

\[
\boxed{
\widehat{\nabla\log M(\bar\gamma_j)}
=
\frac{1}{n_{\mathrm{acc}}}\sum_{r=1}^{n_{\mathrm{acc}}
\mathcal W^\top P_b\bigl(\beta^{(r)} - \mathcal W\bar\gamma_j\bigr),
}
\]

\[
\boxed{
c_j = -\widehat{\nabla\log M(\bar\gamma_j)}.
}
\]

Per group \(j\): accumulate \(\mathcal W_{jk}^\top P_b(\beta_j^{(r)}-\mathcal W_j\bar\gamma_j)\)
as in §3.5.

**Optional (more efficient):** reuse inner A/R trials with self-normalized
importance sampling. With the **same accept statistic** as `rNormalGLM.cpp`
(\(\mathrm{LLconst}[J]+c_J^\top\beta-\mathrm{NLL}\ge\log U\), equivalent to
\(\pi/(\tilde a\, g)\) without evaluating \(\tilde a\)), with weights
\(w_k := \exp\bigl(\mathrm{LLconst}[J_k]+c_{J_k}^\top\beta^{(k)}-\mathrm{NLL}(\beta^{(k)})\bigr)\)
(the stored accept ratio; \(\le 1\) under majorization),

\[
\widehat{\mathbb E}_{\pi}[h(\beta)]
\approx
\frac{\sum_k h(\beta^{(k)})\, w_k}{\sum_k w_k},
\qquad
h(\beta) := \mathcal W^\top P_b(\beta-\mathcal W\bar\gamma_j),
\]

which uses proposals from \(g\) without discarding rejections. The user’s “accepted
draws only” path is **unbiased** and simpler; importance pooling reduces variance
when many inner trials are already run for step 3.

**Hessian (if outer standard form needs it):** use accepted draws for

\[
\widehat{\mathrm{Cov}}_{\beta\mid\bar\gamma_j,y}\bigl(
\mathcal W^\top P_b(\beta-\mathcal W\bar\gamma_j)
\bigr),
\]

plus the exact term \(-\mathcal W^\top P_b\mathcal W\) from §3.5, or build outer
\(H_\gamma\) by numerical differences of \(\widehat{\nabla\log M}\) on the
\(\gamma\) grid.

**Cost:** expect \(\sim n_{\mathrm{acc}}\times \mathbb E[\text{candidates per accept}]\)
inner `f2` calls per outer tangency — on top of the \(N_{\mathrm{in}}\) trials in
step 2 (which can **overlap** if one long inner run supplies both \(\widehat a\) and
accepts for step 3).

### 4.5 Step 4 — Assemble the approximate population LSD

Repeat steps 1–3 for every outer tangency
\(\bar\gamma_j\in\{\gamma^\dagger_i\pm\omega_i\}\) (or sparse subset from
`EnvelopeOpt` in \(\gamma\)-space). Collect **approximate** face data:

- tangency locations \(\bar\gamma_j\);
- outer subgradients \(\widehat c_j \approx \nabla\mathrm{NLL}_\gamma(\bar\gamma_j)\) (step 3);
- face log-weights from \(\widehat{\mathrm{NLL}}_\gamma(\bar\gamma_j)\), combining
  the **exact** population prior with \(\widehat{\mathrm{NLL}}_M(\bar\gamma_j)\)
  from step 2.

Then run the **outer** `EnvelopeBuild` / `setlogP_C2` / PLSD on \(\gamma\):

1. Standardize at \(\gamma^\dagger\) using outer \(H_\gamma\) (from Hessian
   estimates or pilot optimization).
2. Form restricted tilted MVN components on each outer face with \((\bar\gamma_j,\widehat c_j)\).
3. Normalize to outer PLSD; deliver `Env_out_γ`.

**Deliverable:** an **approximate** population LSD mixture \(g_\gamma(\gamma)\)
majorizing \(\pi_\gamma(\gamma\mid y)\) (up to MC error in the face \(\mathrm{NLL}\)
and subgradient estimates, plus the usual envelope slack). This is the density
needed to **propose** population parameters; accept–reject (or majorization bounds)
uses the same stored accept test as Block 1 — the outer slack constant need not be
known explicitly at run time. Quality improves as inner budgets refine
\(\widehat M\) and \(\widehat{\nabla\log M}\).

### 4.6 Cost summary of the nested path

| Layer | Dominant cost | Scales with |
|---|---|---|
| Step 1 (inner build) | `EnvelopeBuild` at each \(\bar\gamma_j\) | \(K_\gamma\times 3^{Jp_{\mathrm{re}}}\) (or `EnvelopeOpt` subset) |
| Step 2 (\(\widehat M\)) | \(N_{\mathrm{in}}\) inner trials → \(\widehat a\) → Theorem 1 | \(K_\gamma\times N_{\mathrm{in}}\) |
| Step 3 (gradient) | Accepted (or all IS) inner trials | \(K_\gamma\times n_{\mathrm{acc}}\times E[\text{candidates}]\) |
| Step 4 (outer build) | One \(\gamma\)-space `EnvelopeBuild` | \(3^{q}\) or `EnvelopeOpt` subset, \(q=\dim\gamma\) |

**Compared to §3.4:** Prékopa guarantees the outer LSD **exists**; this pipeline
**derives an approximation** to it by nested Block‑1 envelopes: inner trials
estimate \(\widehat a\) and Theorem 1 yields \(\widehat{\mathrm{NLL}}_M\); accepted
\(\beta\) yield \(\widehat{\nabla\mathrm{NLL}}_M\); outer `EnvelopeBuild` turns
those face approximations into \(g_\gamma\). Same technology as `glmerb()`, extended
to the population block.

### 4.7 Implementation sketch (not shipped)

| Component | Suggested hook |
|---|---|
| Inner envelope at \(\gamma\) | Existing `EnvelopeBuild` / `BlockEnvelopeBuild` with \(\gamma\) fixed in `f2`/`f3` |
| \(\widehat M(\gamma)\) | \(\widehat a\approx 1/\widehat p_{\mathrm{acc}}\) + Theorem 1 with build \(g,c\) |
| \(\widehat a\) / slack | Empirical inner acceptance rate (not needed at sampler run time) |
| Gradient accumulator | On accept, add \(\mathcal W^\top P_b(\beta-\mathcal W\gamma)\); or IS ratio from step 2 |
| Outer envelope / LSD | `PopulationEnvelopeBuild()` — same C++ as `EnvelopeBuild`, face tables from steps 1–3 |

Pilot tolerances \((\tau_M,\, n_{\mathrm{acc}})\) should be validated on a small
\(J,p,q\) instance before scaling.

---

## 5. Building an optimally centered LSD for \(\pi_\gamma\) (outer shell)

**§4** is the operational nested pipeline: inner envelopes at each \(\gamma\)
tangency approximate the face \(\mathrm{NLL}\) and subgradients; the outer build
produces \(g_\gamma\). This section restates the **outer** envelope steps abstractly
(same as Chapter-A05 on \(\gamma\)), assuming those face oracles — exact or MC —
are available.

\[
-\log \pi_\gamma(\gamma\mid y)
= \tfrac12(\gamma-\mu_0)^\top \Lambda_\gamma(\gamma-\mu_0)
+ \mathrm{NLL}_M(\gamma)
+ \text{const},
\qquad \mathrm{NLL}_M := -\log M,\ \ \mathrm{NLL}_M\ \text{convex}.
\]

The steps below are the same as Block 1 / Chapter-A05, but with population
\(\mathrm{NLL}_M\) in place of the conditional GLM negative log-likelihood.

### Step A — Mode and curvature

1. Find \(\gamma^\dagger = \arg\min \mathrm{NLL}_\text{total}(\gamma)\) (BFGS +
   user-supplied `f2`/`f3` analogues).
2. Approximate Hessian \(H_\gamma := \nabla^2\bigl(-\log\pi_\gamma(\gamma^\dagger\mid y)\bigr)\).

**Cost driver:** each `f3` evaluation may require the **nested integral** in §6
(unless the §4 oracle is used).

### Step B — Standard form (Chapter-A05 §5)

Apply the same pipeline as `glmb_Standardize_Model()`:

- Cholesky whitening so \(H_\gamma\) is diagonal at the mode.
- Diagonal shift \(P \mapsto D\) with \(P-D\succ 0\) (Remarks 12–13) so prior
  precision becomes identity in standard coordinates.

Output: diagonal data precision \(A\), mode \(\gamma^\dagger\), transforms
`L2Inv`, `L3Inv` for back-mapping.

### Step C — Grid size (`Gridtype`, `EnvelopeOpt`)

Per dimension \(i\), posterior precision in standard form is \(1+a_i\).
`Gridtype` 1–4 and `EnvelopeOpt(a,n)` implement the paper’s trade-off:

- **Build cost** \(\propto 3^{(\#\ \text{three-point dims})}\) face evaluations.
- **Sampling cost** \(\propto n \times \prod_i \text{scaleest}_i\) expected
  slope evaluations per accepted draw (Chapter-A05 §6).

For GPU builds, `n` is scaled by OpenCL core count when optimizing grid richness.

### Step D — Face construction (`EnvelopeBuild`)

For each face \(j\) in the \(3^p\) partition (or sparse subset):

1. Tangency point \(\bar\gamma_j\) from \(\gamma^\dagger \pm \omega\).
2. Evaluate \(\mathrm{NLL}_\gamma(\bar\gamma_j)\) and **subgradient**
   \(c_j = \nabla \mathrm{NLL}_\gamma(\bar\gamma_j)\).
3. Store box limits `Lint`, `cbars`, face log-weights `logP`.
4. Normalize to **PLSD** face probabilities.

**Deliverable:** envelope object (`Env_out`) usable by `rNormalGLM_std`-style
draws, but with **population** `f2`/`f3` implementations.

---

## 6. The nested expectation problem (why build is expensive without §4)

For the **true** marginal target,

\[
M(\gamma)
=
\int
\exp\!\Bigl(\sum_j \ell_j(\beta_j)\Bigr)\,
\exp\!\Bigl(-\tfrac12(\beta-\mathcal W\gamma)^\top P_b^{\mathrm{stack}}(\beta-\mathcal W\gamma)\Bigr)\,d\beta.
\]

Gradients w.r.t. \(\gamma\) involve **posterior expectations under
\(\pi(\beta\mid\gamma,y)\)**:

\[
\nabla_\gamma \log M(\gamma)
=
\mathbb E_{\beta\mid\gamma,y}\bigl[
\nabla_\gamma \log p(\beta\mid\gamma)
\bigr]
\]

(with explicit form from the \(\mathcal W_j\gamma\) mean structure). Each
**exact** gradient evaluation is an integration problem over \(\beta\).

**Recursive structure (resolved by §4):**

```
Build approximate population LSD g_γ for π_γ
  → for each outer tangency γ̄_j:
        EnvelopeBuild(β | γ̄_j)              step 1
        NLL̂_M(γ̄_j) via â + Theorem 1       step 2  (face log-density piece)
        ∇̂ NLL_M from accepted β             step 3  (face subgradient piece)
  → EnvelopeBuild(γ) → g_γ                   step 4  (population LSD mixture)
  → sample γ ~ g_γ, accept w.r.t. π_γ        (Block 1 pattern on γ)
```

Before §4, the recursion looked open-ended because each face evaluation seemed to
require a fresh inner integration. The pipeline above **closes the loop**: inner
envelopes plus MC approximations supply the face data; outer `EnvelopeBuild` **derives**
\(g_\gamma\), close to the true population LSD as inner draw budgets grow.

So the user’s intuition is right:

| Phase | Typical cost | Why |
|---|---|---|
| **(i) Derive / build LSD** | **High** | Mode + Hessian + \(3^p\) (or optimized subset) faces; each face needs **likelihood + subgradient**; subgradient needs **conditional expectations** that are themselves MC or inner envelopes |
| **(ii) Draw candidate from mixture** | **Low** | One PLSD index + one tilted MVN draw + box check; \(O(p)\) or \(O(\log K)\) after tables built |
| **(iii) Accept–reject** | **High** | Each rejection repeats **full** \(\log\pi_\gamma(\gamma\mid y)\) (and possibly inner accuracy for gradients if reused); expected repetitions \(\prod_i \sqrt{1+a_i}\) in favorable cases, worse with loose envelopes |

**Asymmetry:** build cost is **up-front**; A/R cost is **per draw** and scales with
sample size \(n\) and envelope tightness.

---

## 7. Two uses of the population LSD (and why the scaling factor matters)

The nested pipeline of §4 produces an **approximate** population LSD mixture
\(g_\gamma(\gamma)\) with \(\int g_\gamma=1\), built from MC face estimates of
\(-\log\pi_\gamma\) and \(\nabla(-\log\pi_\gamma)\). That object supports two
distinct uses — **sampling** and **certification** — with different requirements
on how well the **scaling factor** must be known.

Throughout, write the **correct** Nygren majorization (§3.6):

\[
\pi_\gamma(\gamma\mid y)
\;\le\;
a_\gamma(\bar\theta)\, g_\gamma(\gamma),
\qquad
\int g_\gamma\,d\gamma = 1,
\]

with face/mixture slack \(a_\gamma\) (or restricted \(\tilde a_\gamma\)). The
sampler folds \(a_\gamma\) into `LLconst` and need not evaluate it explicitly;
**tail certification** needs a **provable upper bound** \(A_\gamma\ge a_\gamma\).

### 7.1 Use A — iid sampling of population parameters

**Mechanism (same as Block 1 / `glmb()`).**

1. **Build** \(g_\gamma\) via §4: inner envelopes at each \(\gamma\) tangency,
   MC estimates of \(\widehat{\mathrm{NLL}}_M\) and \(\widehat{\nabla\mathrm{NLL}}_M\),
   outer `EnvelopeBuild`.
2. **Propose** \(\gamma\sim g_{\gamma,J}\): draw face \(J\sim\mathrm{PLSD}\), then
   the restricted tilted MVN on that face.
3. **Accept** with the stored test (equivalent to
   \(\pi_\gamma(\gamma\mid y)/\bigl(a_\gamma\, g_\gamma(\gamma)\bigr)\),
   implemented as `LLconst` + subgradient correction \(-\mathrm{NLL}_\gamma\)
   \(\ge \log U\)).

**Why this works without knowing \(a_\gamma\).** The accept probability is
**automatically** correct if the face tables are internally consistent — exactly
as in shipped Block 1. You never form \(a_\gamma\) as a number; you only observe
\(\widehat p_{\mathrm{acc}}\approx 1/a_{\mathrm{eff}}\) afterward.

**Role of MC approximations.** Face \(\mathrm{NLL}\) and subgradients from inner
draws introduce **build error** in \(g_\gamma\): the mixture may fail to majorize
\(\pi_\gamma\) pointwise, or may be loose. For **operational** iid sampling you
refine inner budgets until acceptance and face checks look stable. Exact
correctness of the majorization is not required for the **logic** of A/R — only
that the accept test matches the **target** \(\pi_\gamma\). (If the envelope is
wrong, acceptance rate drops or the test can fail; that is an engineering check,
not a certificate.)

**Output:** iid (or approximately iid, if face MC error is small) draws from
\(\pi_\gamma(\gamma\mid y)\) without running a Gibbs chain on \((\beta,\gamma)\).

### 7.2 Use B — upper bounds on tail probabilities (majorization)

**Mechanism (same pattern as `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md`).**

Fix a tail / escape set \(E\subseteq\Gamma\) (e.g. Chapter C05
\(\widetilde C_d^{\,c}=\{\Psi>d\}\), or a half-space, or \(\{\|\gamma-\gamma^\star\|>r\}\)).
Suppose on \(E\) you have a **pointwise majorization**

\[
\boxed{
\pi_\gamma(\gamma\mid y)
\;\le\;
A_\gamma\, g_\gamma(\gamma)
\qquad\text{for all }\gamma\in E,
}
\]

where \(g_\gamma\) is a **known** density (LSD mixture with explicit Gaussian
pieces) and \(A_\gamma\ge a_\gamma\) is a **certified upper bound** on the slack.
Then for any measurable \(A\subseteq E\),

\[
\boxed{
\pi_\gamma(A\mid y)
=
\int_A \pi_\gamma(\gamma\mid y)\,d\gamma
\;\le\;
A_\gamma \int_A g_\gamma(\gamma)\,d\gamma.
}
\]

In particular, **escape / tail mass**

\[
\pi_\gamma(E\mid y)
\;\le\;
A_\gamma\, G_\gamma(E),
\qquad
G_\gamma(E):=\int_E g_\gamma(\gamma)\,d\gamma.
\]

Because \(g_\gamma\) is a mixture of restricted multivariate normals on known
boxes, \(G_\gamma(E)\) is often computable (or further boundable) in closed form
— same technology as face CDFs in `EnvelopeBuild`.

**Link to drift-free TV (Chapter C05).** With a minorization on a central set \(C\)
and majorization on \(C^c=E\), Theorem 1 gives schematically

\[
\|q_n(\cdot\mid C)-\pi_\gamma\|_{TV}
\;\le\;
(1-\varepsilon)^n
+
\underbrace{A_\gamma\, G_\gamma(C^c)}_{\text{escape mass bound}}.
\]

The population LSD route replaces the **fixed Gaussian** \(q_2\) in
`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` with a **data-adaptive, mode-centered**
\(g_\gamma\), but the **logic is identical**: integrate a **dominating density**
times a **scaling factor**.

### 7.3 Why the scaling-factor bound is key

The scaling factor is the difference between “we have an envelope mixture” and
“we have a **certified** probability inequality.”

| Question | Needs \(g_\gamma\) only? | Needs bound \(A_\gamma\ge a_\gamma\)? |
|---|---|---|
| Propose candidates | Yes | No (`LLconst` suffices) |
| Accept w.r.t. \(\pi_\gamma\) | Yes (exact test) | No |
| Estimate operational slack | No | Empirical \(\widehat p_{\mathrm{acc}}\approx 1/a\) suffices |
| **Prove** \(\pi_\gamma(E)\le \delta\) | No | **Yes** |

**Without \(A_\gamma\).** Knowing \(g_\gamma\) and even knowing that some
\(a_\gamma\) exists with \(\pi_\gamma\le a_\gamma g_\gamma\) does **not** yield
a numeric tail bound: \(\int_E \pi_\gamma\) could be anywhere in
\([0,\int_E a_\gamma g_\gamma]\) unless you **control \(a_\gamma\) from above**.

**With \(A_\gamma\).** One integration (or Gaussian tail formula) gives a **finite,
auditable** upper bound. Chapter C05’s existential “pick \(d(\delta)\)” becomes
“pick \(d\) so \(A_\gamma G_\gamma(\widetilde C_d^{\,c})\le\delta\)” — provided
\(A_\gamma\) is **certified**, not merely estimated.

**Operational vs certified slack.**

- **Operational (sampling):** \(\widehat a\approx 1/\widehat p_{\mathrm{acc}}\) after
  inner trials — good for tuning and reporting when the closed-form build constants
  are not used explicitly.
- **Certified (tail bounds):** require \(A_\gamma\ge a_\gamma\) with
  \(\pi_\gamma(\gamma\mid y)\le A_\gamma g_\gamma(\gamma)\) **provably** on \(E\).

#### JASA: explicit scaling in the multivariate normal case

The earlier remark that \(a\) is “typically not known” refers to **run-time sampling**:
`LLconst` implements the accept test without forming \(a\) as a user-visible scalar.
For **certification**, the JASA paper (Nygren & Nygren 2006; `vignette("Chapter-A08")`)
does **not** leave the scaling factor as an unknown to be guessed from acceptance
rate alone — in **standard form** it is **computable and boundable**.

**Face constant (Theorem 1).** With prior \(\mathcal N(0,I)\), diagonal posterior
precision at the mode, and subgradient \(c(\bar\theta)\) at tangency \(\bar\theta\),

\[
a(\bar\theta)
=
\frac{
g(\bar\theta)\,\exp\!\bigl(\tfrac12\|c(\bar\theta)\|^2\bigr)
}{
f(y)\,\exp\!\bigl(-c(\bar\theta)^\top\bar\theta\bigr)
},
\]

with restricted \(\tilde a(\bar\theta)\) from Remark 5 obtained by multiplying by the
**factorized** box integrals (Remark 6) — products of normal CDF differences across
coordinates. This is exactly what `EnvelopeBuild` / `setlogP_C2` evaluate (`logU`,
`loglt`, `logrt`, `logP`, `LLconst`). So for Block 1 (and for an outer build once
face \(\mathrm{NLL}\)/subgradients are fixed), **\(a\) is known from the build**, not
only estimable post hoc.

**Universal MVN upper bound (Theorems 2–3).** For Normal regression in standard form
with the paper’s **three-interval** partition per coordinate:

\[
\tilde a
\;\le\;
\left(\frac{2}{\sqrt\pi}\right)^k
\qquad (k=\dim\theta),
\]

with univariate limit \(\tilde a^\ast(N)\to 2/\sqrt\pi\) as sample size \(N\to\infty\)
(Theorem 2). Chapter-A05 encodes this in `EnvelopeOpt`: single-point faces use
\(\sqrt{1+a_i}\) expected rejection burden; three-point faces use the constant
\(2/\sqrt\pi\) per promoted dimension. So the **size of the scaling factor** in the
multivariate normal case is a **paper theorem**, not an MC mystery.

**In practice: non-Gaussian log-concave likelihoods.** Theorem 3 is stated for
Normal regression, but the shipped pipeline (`glmb()`, Block 1 of `glmerb()`) applies
the **same** envelope algebra to logit, Poisson, and other **log-concave** GLM
likelihoods. Chapter-A08 (after Theorem 3) notes that posterior asymptotic normality
suggests **analogous practical tightness** as the posterior concentrates; Chapter-A05
(Example 1 discussion) expects similar acceptance behaviour when the data are
“normal-like.” Empirically:

- the **explicit** Remark 6 face constants from `EnvelopeBuild` (not merely
  \((2/\sqrt\pi)^k\)) are the operational majorization for logit/Poisson/etc.;
- the **universal** bound \(\tilde a\le (2/\sqrt\pi)^k\) often remains a **valid
  conservative** scaling factor in these models — frequently much looser than the
  build’s actual \(\tilde a\), but still usable for **Use B** when a certified
  upper bound is needed without re-deriving a likelihood-specific theorem;
- acceptance rates in production match the JASA prediction far more closely than
  naive \(f/g\) majorization would allow.

So for **Block 1** and (by analogy) the population block in §4, the MVN scaling-factor
**machinery** is the right default even when the “likelihood” term is integrated
\(\mathrm{NLL}_M\) rather than a Gaussian regression log-likelihood. What is
**proved** for non-Gaussian targets is Theorem 1 + Remark 6 pointwise dominance given
correct subgradients; what is **observed in practice** is that the MVN **bound**
\((2/\sqrt\pi)^k\) often remains safe. Small \(n\), strong curvature asymmetry, or
poor standardization can still loosen the envelope — then explicit build \(\tilde a\)
(or MC validation) matters more than the universal constant alone.

**Implication for the population block.**

| Layer | Scaling factor status |
|---|---|
| Block 1 \(\beta\mid\gamma,y\) | Explicit \(\tilde a\) from build; MVN bound \((2/\sqrt\pi)^{Jp_{\mathrm{re}}}\) (theorem for Normal data; **practically** valid for log-concave GLMs) |
| Outer \(\pi_\gamma(\gamma\mid y)\) | Same JASA template once \(\mathrm{NLL}_M\) plays the role of \(-\log L\); Prékopa + asymptotic normality suggest analogous tightness when \(\pi_\gamma\) is concentrated |
| §4 MC face approximations | **Additional** error on top of JASA algebra — must be bounded separately for strict certification |

For **Use B** (tail bounds), the JASA MVN results are the main reason a **certified**
\(A_\gamma\) is realistic: you can take \(A_\gamma = (2/\sqrt\pi)^q\) (or the explicit
\(\tilde a\) from the outer `EnvelopeBuild` tables) when the population envelope is
in standard form and the three-interval construction applies — **without** relying on
\(\widehat p_{\mathrm{acc}}\) alone. MC acceptance rate then **validates** the build
rather than **defining** the constant. Where the target is not close to MVN (heavy
tails, small \(n\)), Theorem 3’s bound may be loose but still **certified**; tighter
\(A_\gamma\) comes from the explicit Remark 6 face constants when subgradients are
certified (not just MC-estimated).

Other certified sources (when JASA MVN bounds are not tight enough):

- hybrid with `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` on thin grid / escape regions;
- generalized LSD chord pieces (§16) beyond the built grid.

**Approximate §4 build and certification.** MC error in \(\widehat{\mathrm{NLL}}_M\)
and \(\widehat{\nabla\mathrm{NLL}}_M\) means the constructed \(g_\gamma\) is
**not automatically** a proven majorant of \(\pi_\gamma\). For **Use A** (iid
sampling), you accept MC refinement. For **Use B** (tail bounds), you must **close
the gap**: either certify the face data, or bound the gap between approximate and
true envelopes — otherwise you only have a **plausible** majorant, not a bound.

**Why population parameters specifically.** Block 2 Gibbs conditions on \(\beta\);
its transition kernel is not the **marginal** \(\pi_\gamma\). Tail budgets for
**drift-free** certification concern the **population block’s** contribution to
TV distance — i.e. mass of \(\pi_\gamma\) in outer regions. An LSD built on the
**integrated** target \(\pi_\gamma(\gamma\mid y)\) aligns the majorant with that
block directly, whereas conditional Block‑2 updates can underestimate escape mass
because they ignore integration over \(\beta\).

### 7.4 Summary

```text
§4 nested pipeline  →  approximate population LSD  g_γ

Use A (iid sampling):     propose γ ~ g_γ ; accept via stored LLconst
                          →  no explicit a_γ needed at run time
                          →  MC refines face NLL / subgradients

Use B (tail bounds):      π_γ ≤ A_γ g_γ  on escape set E
                          →  π_γ(E) ≤ A_γ G_γ(E)
                          →  JASA: explicit ã or (2/√π)^q in MVN standard form
                          →  MC â validates; strict proof needs certified face NLL/∇
```

Both uses share the **same** \(g_\gamma\) machinery; they diverge on whether
**provable control of the scaling factor** is required.

---

## 8. Relation to current sampler architecture

```text
Two-block Gibbs (glmerb / lmerb)
────────────────────────────────
Block 1:  β | γ, y     →  LSD envelope (glmbayesCore EnvelopeBuild)  ✓ shipped
Block 2:  γ | β, y     →  conjugate / ING / specialized              ✓ shipped
Marginal: γ | y        →  LSD on π_γ via §4 pipeline              ○ proposed; not shipped
Marginal: β | y (∫γ)   →  LSD on π̃(β|y) (Prékopa)                ○ explicit ∇; not primary path
```

Block 1 works because **fixing \(\gamma\)** makes the conditional log-posterior
**exactly** MN prior on \(\beta\) + log-concave GLM likelihood — the JASA setup.

Block 2 in hierarchical models often uses **conjugate** structure for
\((\gamma,\tau^2)\) rather than LSD, because the conditional is already tractable.

A **marginal-population LSD** would replace or supplement Block 2 when:

- one wants **iid \(\gamma\)** draws without cycling \(\beta\), or
- one needs a **majorant** tied to the **integrated** target rather than the
  conditional.

---

## 9. Cost accounting summary (user hypothesis formalized)

Write \(C_{\mathrm{build}}\), \(C_{\mathrm{prop}}\), \(C_{\mathrm{AR}}\) for one
accepted draw’s amortized components.

**Build (once per envelope refresh):**

\[
C_{\mathrm{build}}
\approx
\sum_{j=1}^{K}
\Bigl(
C_{f}(\bar\gamma_j) + C_{\nabla f}(\bar\gamma_j)
\Bigr),
\]

where \(C_{\nabla f}\) hides the nested \(\mathbb E[\cdot\mid\gamma]\) chain.
With full \(3^p\) grid, \(K=3^p\); `EnvelopeOpt` reduces \(K\) at the price of
looser acceptance.

**Proposal (per candidate):**

\[
C_{\mathrm{prop}} = O(p) + O(\log K)
\]

after PLSD/CDF tables exist (tilted MVN draw + box indicator).

**Accept–reject (per accepted draw):**

\[
\mathbb E[\text{candidates per accept}] \approx \prod_{i=1}^p s_i,
\qquad
s_i \approx \sqrt{1+a_i}\ \text{(single tangent)}\
\text{or}\ \tfrac{2}{\sqrt\pi}\ \text{(three tangents)},
\]

times \(C_{\mathrm{eval}}(\log\pi_\gamma)\) per candidate. For marginal targets,
\(C_{\mathrm{eval}}\) is **not** separable across coordinates unless
\(\log M(\gamma)\) factorizes.

Hence:

1. **Building** the envelope is expensive when \(q=\dim\gamma\) is moderate–large
   or when inner expectations are needed for each subgradient.
2. **Generating** from the mixture is cheap.
3. **Accept–reject** is expensive when the envelope is loose **or** when each
   acceptance requires a costly marginal likelihood evaluation.

This matches the ordering (i)–(iii) in the user prompt.

---

## 10. Alternative: LSD on \(\widetilde\pi(\beta\mid y)\) instead of \(\pi_\gamma\)

If the scientific target is **group coefficients after integrating population
uncertainty**, use the **Schur marginal** (`LOGIT_MARGINAL_INTEGRATE_GAMMA.md`):

\[
-\log \widetilde\pi(\beta\mid y)
=
\tfrac12(\beta-\mu_\beta)^\top \Lambda_\beta(\beta-\mu_\beta)
+ \sum_j \bigl(-\ell_j(\beta_j)\bigr).
\]

**Advantages for LSD:**

- No nested \(\mathbb E[\cdot\mid\gamma]\) in the gradient; `f3` is explicit
  from GLM + \(\Lambda_\beta\).
- Same `EnvelopeBuild` / Block 1 machinery, different prior precision
  (\(\Lambda_\beta\) couples groups but remains quadratic).

**Disadvantages:**

- Dimension \(n = J p_{\mathrm{re}}\) may be large → \(3^n\) grid explosion unless
  `EnvelopeOpt` / blockwise envelopes (`BLOCK_ING_RINDEPNORMALGAMMA_REG.md`) apply.
- Does not directly sample \(\gamma\).

Tail certification for this target is developed in
`LAPLACE_PROFILE_GROUP_MARGINAL_BOUND.md` (Laplace + Prop 2), which is a
**different** majorization route than LSD.

---

## 11. Open research directions

1. **Implement §4** — `PopulationEnvelopeBuild()` with inner-oracle hooks; tune
   \((N_{\mathrm{in}}, n_{\mathrm{acc}}, \tau_M)\); validate on Gaussian closure
   where \(M(\gamma)\) is analytic.

2. **Efficient \(\nabla \mathrm{NLL}_M\)** — importance pooling across all inner
   trials (§4.4) vs accepted-only path; reuse one inner run for steps 2 and 3.

3. **Inner loop budget** — cache / interpolate inner envelopes between nearby
   \(\bar\gamma_j\); quantify bias vs full rebuild per tangency.

4. **Hybrid certification** — LSD majorant for **sampling**; Prop 2 / profile
   bounds for **drift-free TV** (`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md`); avoid
   requiring one envelope to serve both roles. For **minorization set sizing** on
   \(\pi_\gamma\), see also `MVN_CALIBRATED_MINORIZATION_SET.md` (tail budget
   \(\delta\), \(r(\delta)\), \(\widetilde C_d\), \(\varepsilon_d\), \(\varepsilon\)).

5. **Optimal centering** — mode is correct for standard theory; investigate
   **translation-invariant** re-centering when \(\Lambda_\gamma\) is flat in
   some directions (non-identifiable contrasts).

---

## 12. References

- Prékopa, A. (1973). Logarithmic concave measures with applications.
  *Acta Scientiarum Mathematicarum*, 32, 301–316. (Marginalization preserves
  log-concavity; see also `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md`,
  `RESTRICTED_GIBBS_MINORIZATION_TV_v2.md`.)
- Nygren, K. N., & Nygren, L. M. (2006). Likelihood subgradient densities.
  *Journal of the American Statistical Association*, 101(475), 1144–1156.
  [doi:10.1198/016214506000000357](https://doi.org/10.1198/016214506000000357)
- **glmbayes** `vignette("Chapter-A05")` — standardization, grid types, envelope build.
- **glmbayes** `vignette("Chapter-A08")` — envelope orchestration and PLSD.
- **lmebayes** `vignette("Chapter-09")` — Block 1 LSD in Poisson GLMM fits.
- **lmebayesCore** `inst/LOGIT_MARGINAL_INTEGRATE_GAMMA.md` — \(\int\gamma\) marginal on \(\beta\).
- **lmebayesCore** `inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` — alternative majorization on \(\gamma\).
- **lmebayesCore** `inst/BLOCK_GIBBS_ERGODICITY_ING.md` §16 — generalized LSD when \(-\log f\) is not convex.
