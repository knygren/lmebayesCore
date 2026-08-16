# Logit static tail certification (group-level, conditional on \(\gamma\))

Certify that **conditional on population parameters \(\gamma\)**, the group-level
posterior on \(\beta\) carries small mass in the rate-**unsafe** region — using
**per-group** tilted-Gaussian bounds and a **two-level partition** of
\(\gamma\)-space. Static, equilibrium **target** law only; no Gibbs chain, no
Foster–Lyapunov drift.

**Companion notes.**

| Note | Content |
|---|---|
| `CONDITIONAL_TAIL_BOUNDS_ONE_DIM_GLMM.md` | Prop 1 template all package likelihoods; \(G_j\), \(E_j(\omega_j)\) |
| `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0–§2 | Factorization given \(\gamma\); Lemmas T, G, Z; Prop L1 |
| `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §IV.2–IV.3 | Safe set \(A/A^c\), box certification |
| `LOGIT_MARGINAL_INTEGRATE_GAMMA.md` | Optional: \(\gamma\)-integrated \(\widetilde{\pi}(\beta\mid y)\) if a full \(\beta\)-marginal normalizer route is needed later |

**Scope.** **Group-level** results only: scalar intercept RE per group,
\(\beta_j \mid \gamma \sim N(x_j \gamma, \lambda_b^{-1})\) (often \(x_j \equiv 1\)).
All certification is built from **\(\pi(\beta_j \mid \gamma, y_j)\)** and
\(\pi(\beta \mid \gamma, y) = \prod_j \pi(\beta_j \mid \gamma, y_j)\).

---

## 1. Setup and definitions

### 1.1 Data and likelihood (group \(j\))

\[
y_{j,i} \mid \beta_j \sim \mathrm{Binomial}(n_{j,i}, p(\eta_{j,i})),\quad
p(\eta) = \mathrm{expit}(\eta),\quad
\eta_{j,i} = \mathrm{offset}_{j,i} + x_j \beta_j .
\]

Group log-likelihood (canonical logit):

\[
\ell_j(\beta_j)
=
\sum_{i=1}^{n_j}
\bigl[y_{j,i}\, x_j \beta_j - n_{j,i}\log(1+e^{x_j \beta_j})\bigr],
\qquad
\ell_j''(\beta_j) \le 0 .
\]

Write \(N_j := \sum_i n_{j,i}\) and

\[
L_j
:=
\sum_{i=1}^{n_j}
\max(y_{j,i}, n_{j,i}-y_{j,i})
\;\le\;
N_j .
\]

**Estimability (interpretation only).** Finite group MLE / no complete separation
is the usual reading of “unsafe-rate” tail stickiness. It is **not** required for
the score bound (valid for all \(y_{j,i} \in \{0,\ldots,n_{j,i}\}\)).

### 1.2 Conditional posterior on \(\beta\) given \(\gamma\)

The **working object** is the **target conditional** (same law as one Gibbs block
\(\beta \mid \gamma, y\), not the sampler path):

\[
\boxed{
\pi(\beta \mid \gamma, y)
=
\prod_{j=1}^{J}
\pi(\beta_j \mid \gamma, y_j),
\qquad
\pi(\beta_j \mid \gamma, y_j)
\propto
e^{\ell_j(\beta_j)}\,
\phi(\beta_j; x_j\gamma, \lambda_b^{-1}) .
}
\]

Groups **factorize given \(\gamma\)**. There is no need to integrate \(\gamma\)
out for the per-group tail machinery below.

### 1.3 Safe and unsafe sets on \(\beta\) (matrix form)

Fix certification threshold \(\tau \in (0,1)\). At a working point
\(\beta = (\beta_1,\ldots,\beta_J)\), use the **\(P\)-matrix** notation of
`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §I.2–§I.4 (likelihood enters
only through data weights).

**Per-group objects (group \(j\)).** Stack designs \(Z_j\) (data) and
\(H_j \in \mathbb{R}^{p_{\mathrm{re}}\times q}\) (level-2 / RE design linking
\(\gamma\) to \(\beta_j\); intercept RE: \(H_j = x_j\)). Write the diagonal
**observation-weight matrix**

\[
W_{j,\mathrm{obs}}(\beta)
:= \mathrm{diag}\bigl(w_{j,1}(\eta_{j,1}(\beta)),\,\ldots,\,w_{j,n_j}(\eta_{j,n_j}(\beta))\bigr),
\qquad
w_{j,i} = n_{j,i}\, p(\eta_{j,i})\bigl(1-p(\eta_{j,i})\bigr)
\ \text{(logit)}.
\]

The **data-precision Gram matrix** and **RE conditional precision block** are

\[
\mathcal{G}_j(\beta)
:= Z_j^\top W_{j,\mathrm{obs}}(\beta)\, Z_j
\;\in\; \mathbb{R}^{p_{\mathrm{re}}\times p_{\mathrm{re}}},
\qquad
B_j(\beta)
:= P_b + \mathcal{G}_j(\beta),
\]

with \(P_b = \Psi^{-1}\) the RE prior precision (scalar intercept:
\(P_b = \lambda_b\)). The **group shrinkage operator** is the middle block of
the coupling \(P_{12,j}P_{22,j}^{-1}P_{21,j}\):

\[
\boxed{
\Omega_j(\beta)
:= P_b\, B_j(\beta)^{-1}\, P_b
\;=\;
P_b\,\bigl(P_b + \mathcal{G}_j(\beta)\bigr)^{-1} P_b
\;\in\; \mathbb{R}^{p_{\mathrm{re}}\times p_{\mathrm{re}}}.
}
\]

**Gibbs blocks and rate matrix.** Stack the joint precision (Chapter-C03 §I.3):

\[
P_{22,j}(\beta) = B_j(\beta),
\qquad
P_{12,j} = -H_j^\top P_b,
\qquad
P_{21} = P_{12}^\top,
\]

\[
\boxed{
P_{11}(\lambda_\gamma)
\;=\;
\underbrace{\sum_{j=1}^{J} H_j^\top P_b H_j}_{P_{11}^{\mathrm{RE}}}
\;+\;
\underbrace{\mathrm{blockdiag}_{k}(V_k^{-1})}_{\text{population prior on }\gamma}.
}
\]

Likelihood enters only through \(P_{22}\) (via \(\mathcal{G}_j\)). The RE prior
\(P_b\) appears in \(P_{22}\), \(P_{12}\), and the first addend of \(P_{11}\).
The **population** prior \(V_k^{-1}\) (scalar: \(\lambda_\gamma\)) appears
**only** in the second addend of \(P_{11}\).

The **coupling matrix** is

\[
\boxed{
S(\beta)
\;:=\;
P_{12}(\beta)\, P_{22}(\beta)^{-1}\, P_{21}(\beta)
\;=\;
\sum_{j=1}^{J}
(-H_j^\top P_b)\,
\bigl(P_b + \mathcal{G}_j(\beta)\bigr)^{-1}\,
(-P_b H_j)
\;=\;
\sum_{j=1}^{J} H_j^\top \Omega_j(\beta)\, H_j,
}
\]

which is **independent of \(\lambda_\gamma\)**. The **Gibbs rate matrix** is

\[
\boxed{
A(\beta;\lambda_\gamma)
\;:=\;
P_{11}(\lambda_\gamma)^{-1/2}\,
S(\beta)\,
P_{11}(\lambda_\gamma)^{-1/2},
\qquad
\lambda^\star(\beta;\lambda_\gamma)
\;:=\;
\lambda_{\max}\bigl(A(\beta;\lambda_\gamma)\bigr).
}
\]

**Flat-\(\gamma\) limit (certification uses this).** Static certification
evaluates the rate at the **weak population-prior limit**
\(\lambda_\gamma \downarrow 0\): replace \(P_{11}(\lambda_\gamma)\) by

\[
\boxed{
P_{11}^{\mathrm{RE}}
\;:=\;
\lim_{\lambda_\gamma \downarrow 0} P_{11}(\lambda_\gamma)
\;=\;
\sum_{j=1}^{J} H_j^\top P_b H_j,
}
\]

and define the **limit rate matrix**

\[
\boxed{
A^{(0)}(\beta)
\;:=\;
{P_{11}^{\mathrm{RE}}}^{-1/2}\,
S(\beta)\,
{P_{11}^{\mathrm{RE}}}^{-1/2},
\qquad
\lambda^{\star,(0)}(\beta)
\;:=\;
\lambda_{\max}\bigl(A^{(0)}(\beta)\bigr).
}
\]

*Proof.* Fix \(\beta\). Since \(S(\beta)\) does not depend on \(\lambda_\gamma\)
and \(P_{11}(\lambda_\gamma)=P_{11}^{\mathrm{RE}}+\mathrm{blockdiag}(V_k^{-1})\)
with \(P_{11}^{\mathrm{RE}}\succ 0\),

\[
A(\beta;\lambda_\gamma)
=
\bigl(P_{11}^{\mathrm{RE}}+\mathrm{blockdiag}(V_k^{-1})\bigr)^{-1/2}
S(\beta)
\bigl(P_{11}^{\mathrm{RE}}+\mathrm{blockdiag}(V_k^{-1})\bigr)^{-1/2}
\;\longrightarrow\;
A^{(0)}(\beta)
\quad (\lambda_\gamma \downarrow 0).
\]

\(P_{22}\), \(P_{12}\), \(\Omega_j\), and \(S\) are **not** limited; only the
\(P_{11}^{-1/2}\) normalization drops the population term. With
\(\lambda_\gamma>0\), the full rate is **smaller** in the usual Loewner /
scalar order (tighter \(\gamma\) prior \(\Rightarrow\) faster certified
contraction).

*Scalar specialization* (\(p_{\mathrm{re}}=1\), \(q=1\), \(P_b=\lambda_b\)).
Keep \(H_j\), \(P_b\), \(B_j(\beta)\), and \(Z_j\) in the notation; write each
summand of \(S(\beta)\) from the blocks \(P_{12,j}=-H_j^\top P_b\),
\(P_{22,j}=B_j(\beta)\), \(P_{21,j}=-P_b H_j\) (left to right):

\[
S(\beta)
=
\sum_{j=1}^{J}
\underbrace{(-H_j^\top P_b)}_{P_{12,j}}
\;
\underbrace{\bigl(P_b + \mathcal{G}_j(\beta)\bigr)^{-1}}_{P_{22,j}^{-1}}
\;
\underbrace{(-P_b H_j)}_{P_{21,j}},
\qquad
B_j(\beta)=P_b+\mathcal{G}_j(\beta),
\qquad
\mathcal{G}_j(\beta)=Z_j^\top W_{j,\mathrm{obs}}(\beta)\, Z_j.
\]

When all rows of \(Z_j\) share one linear predictor
\(\eta_{j,i}=\mathrm{offset}_{j,i}+x_j\beta_j\) (intercept RE design
\(H_j=x_j\)), \(\mathcal{G}_j(\beta)=W_j^{\mathrm{sum}}(\beta)\),
\(B_j(\beta)=\lambda_b+W_j^{\mathrm{sum}}(\beta)\) (scalar; \(W_j^{\mathrm{sum}}\)
defined below), and the middle factor is
\(\bigl(\lambda_b + W_j^{\mathrm{sum}}(\beta)\bigr)^{-1}\). With
\(P_{11}^{\mathrm{RE}}=\sum_{\ell=1}^{J} H_\ell^\top P_b H_\ell\),
factor \(\bigl(P_{11}^{\mathrm{RE}}\bigr)^{-1/2}\) **outside** the sum (same
value as distributing it into each summand):

\[
\boxed{
A^{(0)}(\beta)
=
\Bigl(\sum_{\ell=1}^{J} H_\ell^\top P_b H_\ell\Bigr)^{-1/2}\,
\left\{
\sum_{j=1}^{J}
(-H_j^\top P_b^{1/2})\,
P_b^{1/2}\,
\bigl(P_b + \mathcal{G}_j(\beta)\bigr)^{-1}\,
P_b^{1/2}\,
\bigl(-P_b^{1/2} H_j\bigr)
\right\}
\Bigl(\sum_{\ell=1}^{J} H_\ell^\top P_b H_\ell\Bigr)^{-1/2},
}
\]
\[
\lambda^{\star,(0)}(\beta)
=
\lambda_{\max}\bigl(A^{(0)}(\beta)\bigr).
\]

The braces contain \(S(\beta)\), with each summand factored as
\((-H_j^\top P_b)=(-H_j^\top P_b^{1/2})P_b^{1/2}\) and
\((-P_b H_j)=P_b^{1/2}(-P_b^{1/2} H_j)\). Equivalently,
\(A^{(0)}=\sum_j
\bigl(P_{11}^{\mathrm{RE}}\bigr)^{-1/2}\,
P_{12,j}\,P_{22,j}^{-1}\,P_{21,j}\,
\bigl(P_{11}^{\mathrm{RE}}\bigr)^{-1/2}\).

Only now specialize \(H_j=x_j\), \(P_b=\lambda_b\),
\(\mathcal{G}_j=W_j^{\mathrm{sum}}\). When \(q=1\), the middle block is the
**group shrinkage weight**
\[
\omega_j(\beta)
\;:=\;
P_b^{1/2}\,
\bigl(P_b + \mathcal{G}_j(\beta)\bigr)^{-1}\,
P_b^{1/2}
\;=\;
\lambda_b\,\bigl(\lambda_b + W_j^{\mathrm{sum}}(\beta)\bigr)^{-1},
\]
and each \(j\)-th addend inside \(\{\,\cdot\,\}\) is
\[
(-H_j^\top P_b^{1/2})\,
\omega_j(\beta)\,
\bigl(-P_b^{1/2} H_j\bigr)
=
\omega_j(\beta)\,(H_j^\top P_b H_j)
=
\lambda_b\, x_j^2\,\omega_j(\beta).
\]

**Prior weight on group \(j\).** \(\Gamma_{j,\mathrm{prior}}:=H_j^\top P_b H_j\) is
the \(j\)-th summand of \(P_{11}^{\mathrm{RE}}\): the **prior precision with
which group \(j\)'s RE prior couples \(\beta_j\) to the population mean
\(\gamma\) (intercept RE: \(\Gamma_{j,\mathrm{prior}}=\lambda_b x_j^2\)). In
the factorized addend, \(\omega_j\) is the **shrinkage** factor (prior vs.\ data
on the \(b_j\) block) and \(H_j^\top P_b H_j\) is the **prior weight** on the
hyper side. Define relative weights
\(\rho_j:=\Gamma_{j,\mathrm{prior}}\big/\sum_{\ell=1}^{J}\Gamma_{\ell,\mathrm{prior}}\)
(\(\sum_j\rho_j=1\)). Distributing the outer sandwich gives per-group
contributions \(\rho_j\,\omega_j(\beta)\); when \(q=1\),
\[
\lambda^{\star,(0)}(\beta)
=
\sum_{j=1}^{J}\rho_j\,\omega_j(\beta).
\]
With \(x_j\equiv 1\), \(\rho_j\equiv 1/J\) and
\[
\lambda^{\star,(0)}(\beta)
=
J^{-1}\sum_{j=1}^{J}\omega_j(\beta)
=
\bar\omega(\beta).
\]
(The boxed \(\omega_j\) below agrees: \(\lambda_{\max}(\Theta_j)=\omega_j\) when
\(q=1\).)
With \(\lambda_\gamma>0\),
\(\lambda^\star(\beta;\lambda_\gamma)
=\lambda^{\star,(0)}(\beta)\,
\bigl(\sum_{\ell} H_\ell^\top P_b H_\ell\bigr)\,
\bigl(\sum_{\ell} H_\ell^\top P_b H_\ell + \lambda_\gamma\bigr)^{-1}
=\lambda^{\star,(0)}(\beta)\,
\bigl(1+\lambda_\gamma/\sum_{\ell} H_\ell^\top P_b H_\ell\bigr)^{-1}\)
(and \(\lambda^{\star,(0)}/(1+\lambda_\gamma/(J\lambda_b))\) when \(x_j\equiv 1\))
\(\le \lambda^{\star,(0)}(\beta)\).
The \(j\)-th addend inside \(\{\,\cdot\,\}\) is **not** \(\Theta_j(\beta)\):
both share the middle shrinkage block
\(P_b^{1/2}\,\bigl(P_b+\mathcal{G}_j\bigr)^{-1}\,P_b^{1/2}\)
(\(=\omega_j\) when \(q=1\)), but \(A^{(0)}\) uses **global**
\(\bigl(\sum_{\ell} H_\ell^\top P_b H_\ell\bigr)^{-1/2}\) outside the sum and
flanks \((-H_j^\top P_b^{1/2})\), \((-P_b^{1/2}H_j)\), whereas
\(\Theta_j\) uses **local** \(\Gamma_{j,\mathrm{prior}}^{-1/2}=
(H_j^\top P_b H_j)^{-1/2}\) (see below).

**Per-group coupling blocks.** Group \(j\) contributes the symmetric
\(q\times q\) block \(H_j^\top \Omega_j(\beta)\, H_j\) to \(S(\beta)\). Write
the **RE prior and total Gram matrices** on the hyper side,

\[
\Gamma_{j,\mathrm{prior}}
:= H_j^\top P_b H_j,
\qquad
\Gamma_{j,\mathrm{tot}}(\beta)
:= H_j^\top B_j(\beta)\, H_j
\;\in\; \mathbb{R}^{q\times q},
\]

both symmetric (assuming \(P_b\succ 0\), \(B_j\succ 0\)). When \(q=1\),
\(\Gamma_{j,\mathrm{prior}}\) and \(\Gamma_{j,\mathrm{tot}}\) are scalars.
Note \(\Gamma_{j,\mathrm{prior}}=H_j^\top P_b H_j\) is the \(j\)-th summand of
\(P_{11}^{\mathrm{RE}}\) — group \(j\)'s **prior coupling weight** on
\(\gamma\), not the population prior on \(\gamma\).

**Symmetric shrinkage matrix (multivariate case).** The limit rate
\(A^{(0)}= {P_{11}^{\mathrm{RE}}}^{-1/2} S {P_{11}^{\mathrm{RE}}}^{-1/2}\)
applies a **global** square-root normalization to the full coupling \(S\).
At group level, the analogous **local** congruence (same pattern, but using
\(\Gamma_{j,\mathrm{prior}}\) instead of the full \(P_{11}^{\mathrm{RE}}\))
is

\[
\boxed{
\Theta_j(\beta)
\;:=\;
\Gamma_{j,\mathrm{prior}}^{-1/2}\,
\bigl(H_j^\top \Omega_j(\beta)\, H_j\bigr)\,
\Gamma_{j,\mathrm{prior}}^{-1/2}
\;=\;
\Gamma_{j,\mathrm{prior}}^{-1/2}\,
(-H_j^\top P_b)\,
\bigl(P_b + \mathcal{G}_j(\beta)\bigr)^{-1}\,
(-P_b H_j)\,
\Gamma_{j,\mathrm{prior}}^{-1/2}
\;\in\; \mathbb{R}^{q\times q},
}
\]

with \(B_j(\beta)=P_b+\mathcal{G}_j(\beta)\) and
\(\mathcal{G}_j(\beta)=Z_j^\top W_{j,\mathrm{obs}}(\beta)\, Z_j\).

which is **symmetric** (hence has real eigenvalues). A plain ratio
\(\Gamma_{j,\mathrm{prior}} / \Gamma_{j,\mathrm{tot}}\) is **not** defined when
\(q>1\); the sandwich above is the multivariate replacement. Equivalently,

\[
\Theta_j(\beta)
\;=\;
\Gamma_{j,\mathrm{tot}}(\beta)^{-1/2}\,
\Gamma_{j,\mathrm{prior}}\,
\Gamma_{j,\mathrm{tot}}(\beta)^{-1/2},
\]

using \(\Gamma_{j,\mathrm{prior}} = H_j^\top P_b H_j\) and
\(\Gamma_{j,\mathrm{tot}} = H_j^\top B_j H_j\); the two displays agree when
\(q=1\). *Proof (equivalence).* With \(B_j=P_b+\mathcal{G}_j\),
\(\Omega_j=P_b B_j^{-1}P_b\), and \(\Gamma_{j,\mathrm{tot}}=H_j^\top B_j H_j\),
the Woodbury identity on the hyper-side projection gives the second display;
for \(q=1\) both reduce to \(\Gamma_{j,\mathrm{prior}}/\Gamma_{j,\mathrm{tot}}\).
\(\square\)

**Scalar summary per group.** Extract

\[
\boxed{
\omega_j(\beta)
\;:=\;
\lambda_{\max}\bigl(\Theta_j(\beta)\bigr),
}
\]

(or any fixed scalar functional of \(\Theta_j\) fixed once and for all when
\(q>1\); \(\lambda_{\max}\) matches the rate-matrix logic). When \(q=1\),
\(\omega_j(\beta)=P_b^{1/2}\,\bigl(P_b+\mathcal{G}_j(\beta)\bigr)^{-1}\,P_b^{1/2}\)
is the middle block in \(A^{(0)}\) above; equivalently, with
\(\Gamma_{j,\mathrm{prior}}=H_j^\top P_b H_j\),

\[
\omega_j(\beta)
=
\lambda_{\max}(\Theta_j)
=
P_b^{1/2}\,
\bigl(P_b + \mathcal{G}_j(\beta)\bigr)^{-1}\,
P_b^{1/2}
=
\lambda_b\,(\lambda_b + W_j^{\mathrm{sum}}(\beta))^{-1}
\quad\text{when }H_j=x_j,\;P_b=\lambda_b,\;
\mathcal{G}_j(\beta)=Z_j^\top W_{j,\mathrm{obs}}(\beta)\, Z_j=W_j^{\mathrm{sum}}(\beta).
\]

*Proof.* Factor
\((-H_j^\top P_b)\,\bigl(P_b+\mathcal{G}_j\bigr)^{-1}\,(-P_b H_j)
=(-H_j^\top P_b^{1/2})\,
P_b^{1/2}\,\bigl(P_b+\mathcal{G}_j\bigr)^{-1}\,P_b^{1/2}\,
(-P_b^{1/2} H_j)\).
With \(\Gamma_{j,\mathrm{prior}}^{-1/2}=(H_j^\top P_b H_j)^{-1/2}\),
\(\Theta_j=\Gamma_{j,\mathrm{prior}}^{-1/2}\,
(H_j^\top P_b B_j^{-1} P_b H_j)\,
\Gamma_{j,\mathrm{prior}}^{-1/2}\), i.e.\ the local congruence of the
coupling. When \(q=1\),
\(\lambda_{\max}(\Theta_j)=P_b^{1/2}\,\bigl(P_b+\mathcal{G}_j\bigr)^{-1}\,P_b^{1/2}\).
\(\square\)

**Do not** replace \(\Gamma_{j,\mathrm{prior}}\) by \(H_j^\top P_b^{-1} H_j\), or
use a non-symmetric sandwich — that breaks the multivariate congruence.

**Interpretation.** \(\Theta_j\) is the group-\(j\) **matrix shrinkage factor**
inside \(S\); it does **not** depend on \(\lambda_\gamma\). When
\(\mathcal{G}_j\) is small, \(B_j \approx P_b\), \(\Theta_j \approx I_q\), and
\(\omega_j \approx 1\); when \(\mathcal{G}_j\) is large, \(\Theta_j \prec I_q\)
and \(\omega_j\downarrow 0\).

Define the **group shrinkage vector**

\[
\boxed{
\boldsymbol{\omega}(\beta)
\;:=\;
\bigl(\omega_1(\beta_1),\,\ldots,\,\omega_J(\beta_J)\bigr)^\top
\;\in\;(0,1]^J .
}
\]

Each coordinate is \(\omega_j(\beta_j)=\lambda_{\max}(\Theta_j(\beta_j))\) from §1.3
(scalar \(q=1\):
\(\omega_j(\beta_j)=\lambda_b/(\lambda_b+W_j^{\mathrm{sum}}(\beta_j))\)).

Define the **average shrinkage weight** (equal-weight summary of
\(\boldsymbol{\omega}\); see prior weights \(\rho_j\) below for the rate-aligned
summary):

\[
\boxed{
\bar\omega(\beta)
:=
J^{-1}\sum_{j=1}^{J} \omega_j(\beta)
=
J^{-1}\,\mathbf{1}^\top \boldsymbol{\omega}(\beta).
}
\]

*Link to the limit rate.* For \(q=1\), \(A^{(0)}\) is scalar and
\(\lambda^{\star,(0)}=\lambda_{\max}(A^{(0)})\) is the **prior-weighted**
functional of \(\boldsymbol{\omega}\),

\[
\boxed{
\Psi(\boldsymbol{\omega})
\;:=\;
\sum_{j=1}^{J}\rho_j\,\omega_j
\;=\;
\boldsymbol{\rho}^\top \boldsymbol{\omega},
\qquad
\rho_j
:=
\frac{H_j^\top P_b H_j}{\sum_{\ell=1}^{J} H_\ell^\top P_b H_\ell},
\qquad
\lambda^{\star,(0)}(\beta)=\Psi\bigl(\boldsymbol{\omega}(\beta)\bigr).
}
\]

**Limit rate at the certification profile \(\boldsymbol{\omega}\).** Let
\(\boldsymbol{\omega}=(\omega_1,\ldots,\omega_J)\) denote a threshold vector
(as in §2, distinct from the state \(\boldsymbol{\omega}(\beta)\)). Say
\(\beta\) **matches** \(\boldsymbol{\omega}\) if each group sits on that
shrinkage level:

\[
\omega_j(\beta_j)=\omega_j
\qquad\forall j
\qquad\Longleftrightarrow\qquad
\boldsymbol{\omega}(\beta)=\boldsymbol{\omega}.
\]

When \(q=1\), \(\lambda^{\star,(0)}(\beta)\) depends on \(\beta\) only through
\(\boldsymbol{\omega}(\beta)\), so this profile rate is well-defined:

\[
\boxed{
\lambda^{\star,(0)}(\boldsymbol{\omega})
\;:=\;
\lambda^{\star,(0)}(\beta)\Big|_{\boldsymbol{\omega}(\beta)=\boldsymbol{\omega}}
\;=\;
\Psi(\boldsymbol{\omega})
\;=\;
\sum_{j=1}^{J}\rho_j\,\omega_j .
}
\]

(Write \(\lambda^\star(\boldsymbol{\omega})\) when the flat-\(\gamma\) limit is
clear.) On the **safe boundary**
\(\partial B(\boldsymbol{\omega})=\{\beta:\omega_j(\beta_j)=\omega_j\ \forall j\}\),
the limit rate equals \(\lambda^{\star,(0)}(\boldsymbol{\omega})\). Under
Lemma 3 (intercept RE, offsets \(0\)), on \(\partial B(\boldsymbol{\omega})\) one has
\(\beta_j=\pm b_j^\star(\omega_j)\) and \(\omega_j(\beta_j)=\omega_j\).

When \(x_j\equiv 1\), \(\rho_j=1/J\) and
\(\lambda^{\star,(0)}(\boldsymbol{\omega})=J^{-1}\mathbf{1}^\top\boldsymbol{\omega}\).
For \(q>1\), \(\lambda^{\star,(0)}(\beta)\) is not determined by
\(\boldsymbol{\omega}\) alone; use \(\lambda_{\max}(A^{(0)}(\beta))\) at a
matching \(\beta\).

*Rate reference (Chapter-C03).* The Gibbs **rate-unsafe** condition
\(\lambda^{\star,(0)}(\beta)\ge\tau\) (equivalently
\(\Psi(\boldsymbol{\omega}(\beta))\ge\tau\) when \(q=1\)) is a separate
aggregate diagnostic. The **certification unsafe set** \(A^c(\boldsymbol{\omega})\)
is defined from per-group tails \(E_j(\omega_j)\) in §2 — not from \(\tau\)
alone.

**Scalar intercept specialization (\(p_{\mathrm{re}}=1\), \(H_j=x_j\),
\(P_b=\lambda_b\)).** If all rows of \(Z_j\) share the same linear predictor
\(\eta_{j,i} = \mathrm{offset}_{j,i} + x_j \beta_j\), then
\(\mathcal{G}_j(\beta) = Z_j^\top W_{j,\mathrm{obs}}(\beta)\, Z_j
= W_j^{\mathrm{sum}}(\beta)\) with
\(W_j^{\mathrm{sum}}(\beta) := \sum_{i\in j} w_{j,i}(\eta_{j,i}(\beta))\)
(the scalar sum used in Lemma 3),
\(B_j(\beta)=\lambda_b+W_j^{\mathrm{sum}}(\beta)\), and

\[
H_j^\top \Omega_j H_j
=
(-H_j^\top P_b)\,
\bigl(\lambda_b + W_j^{\mathrm{sum}}(\beta)\bigr)^{-1}\,
(-P_b H_j),
\qquad
\omega_j(\beta)
=
\lambda_{\max}\bigl(\Theta_j(\beta)\bigr)
=
\lambda_b\,(\lambda_b + W_j^{\mathrm{sum}}(\beta))^{-1},
\qquad
\bar\omega(\beta) = J^{-1}\,\mathbf{1}^\top \boldsymbol{\omega}(\beta).
\]

*(Notation: \(W_j^{\mathrm{sum}}\) is the logit **curvature sum**; do not confuse
with the diagonal matrix \(W_{j,i}\) or the level-2 design \(\mathcal{W}_j\) in
`notation.md`.)*

### 1.4 Conditional unsafe mass (primary target)

Fix a threshold vector \(\boldsymbol{\omega}=(\omega_1,\ldots,\omega_J)\)
(§2; do **not** confuse with the state vector \(\boldsymbol{\omega}(\beta)\)).
The unsafe set \(A^c(\boldsymbol{\omega})\) is defined in §2; its
\(\beta\)-section is \(A^c(\boldsymbol{\omega})|_\beta=\bigcup_j \mathcal{E}_j(\omega_j)\).
Because this depends on \(\beta\) only given \(\boldsymbol{\omega}\),

\[
\boxed{
g(\gamma;\boldsymbol{\omega})
:= \pi\bigl(A^c(\boldsymbol{\omega})\big|_\beta \mid \gamma, y\bigr)
\in [0,1] .
}
\]

**Certification goal (given \(\gamma\)).** For tolerances \(\varepsilon_{1,j}\) and
aggregate \(\varepsilon_1^{\mathrm{eff}}\) (§4), show

\[
g(\gamma;\boldsymbol{\omega}) \le h(\gamma;\boldsymbol{\omega}) \le \varepsilon_1^{\mathrm{eff}}
\quad\text{on a large low-risk set of \(\gamma\),}
\]

with \(h(\gamma;\boldsymbol{\omega})\) the certified profile from §3. Control of
**high-risk** \(\gamma\) (§4) completes a bound on the **population average**
\(E_{\pi(\gamma \mid y)}[g(\gamma;\boldsymbol{\omega})]\) when that average is needed.

### 1.5 Optional: average over \(\pi(\gamma \mid y)\)

If the integrated-over-\(\gamma\) probability is required,

\[
\pi_{\beta \mid y}\bigl(A^c(\boldsymbol{\omega})\big|_\beta\bigr)
=
E_{\pi(\gamma \mid y)}[g(\gamma;\boldsymbol{\omega})] .
\]

This is a **single outer integral**; it is not the primary object of this note.
§4 gives a two-level bound on this average from per-\(\gamma\) certificates.

---

## 2. Per-group lemmas (drift-free)

### Lemma 1 (bounded score)

For every \(\beta_j\) and every admissible count \(y_{j,i}\),

\[
|\ell_j'(\beta_j)| \le L_j .
\]

*Proof.* \(y_{j,i}-n_{j,i}p(x_j\beta_j) \in (y_{j,i}-n_{j,i}, y_{j,i})\) since
\(p\in(0,1)\); bound each term and sum. \(\square\)

### Lemma 2 (tangent minorant)

For any \(b_0 \in \mathbb R\) and all \(\beta_j\),

\[
\ell_j(\beta_j)
\ge
\ell_j(b_0) - L_j|\beta_j - b_0| .
\]

*Proof.* Concavity plus Lemma 1 (rate note §0.3 Lemma T). \(\square\)

**Tail event \(E_j(\omega_j)\) (definition).** Fix a per-group threshold
\(\omega_j \in (0,1]\) (one component of the certification vector
\(\boldsymbol{\omega}\)). Recall \(\omega_j(\beta_j)=\lambda_{\max}(\Theta_j(\beta_j))\)
from §1.3 (scalar \(q=1\):
\(\omega_j(\beta_j)=\lambda_b/(\lambda_b+W_j^{\mathrm{sum}}(\beta_j))\)).
The **per-group tail event** is a subset of the coordinate space — scalar
\(q=1\): \(\beta_j\in\mathbb{R}\); multivariate RE: \(\beta_j\in\mathbb{R}^{p_{\mathrm{re}}}\):

\[
\boxed{
E_j(\omega_j)
\;:=\;
\{\,\beta_j : \omega_j(\beta_j) \ge \omega_j\,\}.
}
\]

(In the inequality, \(\omega_j(\beta_j)\) is shrinkage **at** \(\beta_j\); the
argument \(\omega_j\) of \(E_j(\cdot)\) is the certification **threshold** for
group \(j\). Do not confuse with the state vector \(\boldsymbol{\omega}(\beta)\).)
Under Lemma 3 (intercept RE, offsets \(0\)), the margin
\(b_j^\star(\omega_j)\) and floor \(\kappa_j(\omega_j)\) are **derived from**
\(\omega_j\) (not chosen independently).

**Multivariate RE (\(q>1\)).** The same display applies with
\(\omega_j(\beta_j)=\lambda_{\max}(\Theta_j(\beta_j))\); a fixed scalar
functional of \(\Theta_j\) may be used if chosen once for all certification
(§1.3). The safe slice on group \(j\) is
\(S_j(\omega_j) := \{\beta_j : \omega_j(\beta_j) < \omega_j\}\).

Write \(\mathcal{E}_j(\omega_j) := \{\,\beta=(\beta_1,\ldots,\beta_J) :
\beta_j\in E_j(\omega_j)\,\}\) for the same event in full \(\beta\)-space
(other coordinates free).

**Global unsafe set \(A^c(\boldsymbol{\omega})\).** Let
\(\boldsymbol{\omega}=(\omega_1,\ldots,\omega_J)\). In \((\gamma,\beta)\)-space,

\[
\boxed{
A^c(\boldsymbol{\omega})
\;:=\;
\bigl\{(\gamma,\beta) : \exists\, j\in\{1,\ldots,J\},\ \beta_j \in E_j(\omega_j)\bigr\}.
}
\]

(Companion notes sometimes write \((\lambda,\beta)\) for \((\gamma,\beta)\).)

Only \(\beta_j\) enters the condition; equivalently
\[
A^c(\boldsymbol{\omega})\big|_\beta
=
\bigcup_{j=1}^{J} \mathcal{E}_j(\omega_j)
=
\{\,\beta : \exists\, j,\ \beta_j \in E_j(\omega_j)\,\}.
\]
The **safe region** is the complement
\(A(\boldsymbol{\omega})|_\beta = B(\boldsymbol{\omega})
:= \{\beta : \beta_j \in S_j(\omega_j)\ \forall j\}\).

Define the **per-group conditional tail mass**

\[
\boxed{
g_j(\gamma;\omega_j)
\;:=\;
\pi\bigl(\beta_j \in E_j(\omega_j) \mid \gamma, y_j\bigr)
=
\pi\bigl(\mathcal{E}_j(\omega_j) \mid \gamma, y\bigr)
}
\]

(under \(\pi(\beta\mid\gamma,y)=\prod_j\pi(\beta_j\mid\gamma,y_j)\)).

**Modeller's choice.** \(\boldsymbol{\omega}\) fixes both \(E_j(\omega_j)\) and
\(A^c(\boldsymbol{\omega})\). Optional **rate cross-check** (Chapter-C03): require
\(\lambda^{\star,(0)}(\boldsymbol{\omega})<\tau\) for a chosen rate threshold \(\tau\).

### Lemma 3 (scalar intercept: margin certificate for \(E_j(\omega_j)\))

Write the **per-observation weight**
\(u(\eta) := p(\eta)\bigl(1-p(\eta)\bigr)\) and
\(\eta_{j,i}(\beta_j) := \mathrm{offset}_{j,i} + x_j \beta_j\), so

\[
W_j^{\mathrm{sum}}(\beta_j)
:=
\sum_{i=1}^{n_j} n_{j,i}\, u\bigl(\eta_{j,i}(\beta_j)\bigr)
\;=\;
\mathcal{G}_j(\beta_j)
\quad\text{(scalar intercept RE, §1.3)}.
\]

**Monotonicity.** \(u(\eta)\) is symmetric about \(0\) and **strictly decreasing**
in \(|\eta|\) for \(|\eta|>0\) (derivative \(u'(\eta)=p(\eta)(1-p(\eta))(1-2p(\eta))\);
see `LOGIT_SINGLE_GROUP_SAFE_REGION.md` §2).

**Intercept RE with common offset (group \(j\)).** Suppose
\(\mathrm{offset}_{j,i}\equiv 0\) and \(x_j \neq 0\), so
\(\eta_{j,i}(\beta_j)=x_j\beta_j\) and
\(W_j^{\mathrm{sum}}(\beta_j)=N_j\, u(x_j\beta_j)\).
Fix the certification threshold \(\omega_j \in (0,1]\) and define the **derived**
curvature floor and margin (all as functions of \(\omega_j\)):

\[
\boxed{
\kappa_j(\omega_j)
:=
\frac{\lambda_b\,(1-\omega_j)}{\omega_j},
\qquad
\eta_j^\star(\omega_j)
:=
\text{the unique }\eta \ge 0
\text{ with }
u(\eta) = \frac{\kappa_j(\omega_j)}{N_j},
}
\]

\[
\boxed{
b_j^\star(\omega_j)
:=
\frac{\eta_j^\star(\omega_j)}{|x_j|}
\qquad
\bigl(\text{equivalently }
\kappa_j(\omega_j) = N_j\, u\bigl(|x_j|\, b_j^\star(\omega_j)\bigr)
= N_j\, p\bigl(\eta_j^\star(\omega_j)\bigr)\bigl(1-p(\eta_j^\star(\omega_j))\bigr)\bigr).
}
\]

(Existence/uniqueness of \(\eta_j^\star(\omega_j)\): \(u(\eta)\) is symmetric,
maximal at \(0\), and strictly decreasing on \(\eta\ge 0\); require
\(\kappa_j(\omega_j)/N_j \le u(0)=1/4\), i.e.\ feasible \(\omega_j\) for group \(j\).)
Then **pointwise** on \(\beta_j\in\mathbb{R}\) (deterministic geometry),

\[
\boxed{
E_j(\omega_j)
=
\{\,\beta_j : \omega_j(\beta_j) \ge \omega_j\,\}
=
\{\,W_j^{\mathrm{sum}}(\beta_j) \le \kappa_j(\omega_j)\,\}
=
\{\,|\beta_j| > b_j^\star(\omega_j)\,\}.
}
\]

*Proof.* Strict monotonicity of \(u\) in \(|\eta|\) gives
\(|x_j\beta_j|> \eta_j^\star(\omega_j) \Leftrightarrow W_j^{\mathrm{sum}}(\beta_j)\le\kappa_j(\omega_j)\)
and \(|x_j\beta_j|< \eta_j^\star(\omega_j) \Leftrightarrow W_j^{\mathrm{sum}}(\beta_j)>\kappa_j(\omega_j)\).
At \(|x_j\beta_j|=\eta_j^\star(\omega_j)\), \(\omega_j(\beta_j)=\omega_j\); with
\(\omega_j(\cdot)=\lambda_b/(\lambda_b+W_j^{\mathrm{sum}})\) strictly decreasing in
\(W_j^{\mathrm{sum}}\), the three sets coincide.
\(\square\)

**Reading.** Lemma 3 certifies the intercept margin identity
\(E_j(\omega_j)=\{\,|\beta_j| > b_j^\star(\omega_j)\,\}\) from the primary
threshold definition §2. Prop L1 bounds \(g_j(\gamma;\omega_j)\) via
\(b_j^\star(\omega_j)\).

*(Offsets.* If \(\mathrm{offset}_{j,i}\) vary, \(E_j(\omega_j)\) need not match
\(\{|\beta_j|>b_j^\star(\omega_j)\}\); keep the primary definition
\(\{\omega_j(\beta_j)\ge\omega_j\}\) and certify \(\omega_j\) from margins on
\(\eta_{j,i}\).)

### Proposition L1 (per-group conditional tail)

**Scalar intercept route.** Recall from §2 that
\[
E_j(\omega_j)
=
\{\,\beta_j : \omega_j(\beta_j) \ge \omega_j\,\}.
\]
Under Lemma 3 (\(\mathrm{offset}_{j,i}\equiv 0\), \(x_j \neq 0\)), this is
\(\{\,|\beta_j| > b_j^\star(\omega_j)\,\}\). Write \(\gamma_j := x_j\gamma\). Then

\[
\boxed{
g_j(\gamma;\omega_j)
=
\pi\bigl(\beta_j \in E_j(\omega_j) \mid \gamma, y_j\bigr)
\le
G_j(\gamma;\omega_j)
}
\]

with closed-form \(G_j(\gamma;\omega_j)\) from rate note §2.2 Prop L1, using the
margin \(b_j^\star(\omega_j)\) from Lemma 3:

\[
G_j(\gamma;\omega_j)
=
\exp\!\Big(L_j\,\lambda_b^{-1/2}\sqrt{2/\pi} + \tfrac{L_j^2}{2\lambda_b}\Big)
\left[
\Phi\!\Big(-\tfrac{b_j^\star(\omega_j) - \gamma_j - L_j/\lambda_b}{\lambda_b^{-1/2}}\Big)
+
\Phi\!\Big(-\tfrac{b_j^\star(\omega_j) + \gamma_j - L_j/\lambda_b}{\lambda_b^{-1/2}}\Big)
\right].
\]

*Proof.* Lemma 2 + tilted Gaussian tail (Lemmas G, Z), integrating over
\(E_j(\omega_j)=\{\,|\beta_j| > b_j^\star(\omega_j)\,\}\) (Lemma 3). \(\square\)

**Link to \(\mathcal{G}_j\).** On \(E_j(\omega_j)\), Lemma 3 gives
\(W_j^{\mathrm{sum}}(\beta_j)\le \kappa_j(\omega_j)\) (intercept RE, offsets \(0\)).
Prop L1 certifies how much \(\gamma\)-conditional mass sits in the
**high-\(\omega_j\)** tail for group \(j\).

---

## 3. Union profile \(h(\gamma;\boldsymbol{\omega})\)

### 3.1 \(A^c(\boldsymbol{\omega})\) is the union of coordinate tails

By construction (§2),
\[
A^c(\boldsymbol{\omega})\big|_\beta
=
\bigcup_{j=1}^{J}\mathcal{E}_j(\omega_j)
=
B(\boldsymbol{\omega})^c,
\qquad
B(\boldsymbol{\omega})
=
\{\,\beta : \beta_j \in S_j(\omega_j)\ \forall j\,\}.
\]
Under \(\pi(\beta\mid\gamma,y)=\prod_j\pi(\beta_j\mid\gamma,y_j)\),

\[
\boxed{
g(\gamma;\boldsymbol{\omega})
=
\pi\bigl(A^c(\boldsymbol{\omega})\big|_\beta \mid \gamma, y\bigr)
=
\pi\bigl(B(\boldsymbol{\omega})^c \mid \gamma, y\bigr)
=
1 - \prod_{j=1}^{J}\bigl(1 - g_j(\gamma;\omega_j)\bigr) .
}
\]

**Scalar intercept certificate.** Under Lemma 3, with margins
\(b_j^\star(\omega_j)\) derived from \(\boldsymbol{\omega}\),

\[
B(\boldsymbol{\omega})
=
\prod_{j=1}^{J}\bigl[-b_j^\star(\omega_j),\, b_j^\star(\omega_j)\bigr].
\]
Optional **rate cross-check** (Chapter-C03 §IV.3): verify
\(B(\boldsymbol{\omega})\subseteq \{\beta:\lambda^{\star,(0)}(\beta)<\tau\}\),
equivalently \(\lambda^{\star,(0)}(\boldsymbol{\omega})<\tau\) when
\(\beta\in\partial B(\boldsymbol{\omega})\), for a chosen \(\tau\).

### 3.2 Certified upper bound

Define the **certified conditional profile**

\[
\boxed{
h(\gamma;\boldsymbol{\omega})
:=
1 - \prod_{j=1}^{J}\bigl(1 - G_j(\gamma;\omega_j)\bigr)
\;\ge\;
g(\gamma;\boldsymbol{\omega}) .
}
\]

*Proof.* Prop L1 gives \(g_j(\gamma;\omega_j)\le G_j(\gamma;\omega_j)\); substitute
into §3.1. \(\square\)

**Do not** use \(\sum_j G_j(\gamma;\omega_j)\) as the primary profile (valid union bound,
but inefficient).

---

## 4. Two-level \(\gamma\) partition

### 4.1 Low-risk \(\gamma\) set \(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\)

Fix per-group tolerances \(\varepsilon_{1,j} \in (0,1)\) and optional high-risk
mass budget \(\varepsilon_2 \in [0,1]\). For each group, the certified
\(\gamma\)-gate is

\[
\Gamma_{\mathrm{low},j}(\omega_j)
:= \{\gamma : G_j(\gamma;\omega_j) \le \varepsilon_{1,j}\}.
\]

The **safe population set** (all group gates pass) is

\[
\boxed{
\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})
:= \bigcap_{j=1}^{J}\Gamma_{\mathrm{low},j}(\omega_j)
= \{\gamma : G_j(\gamma;\omega_j) \le \varepsilon_{1,j}\ \forall j\}.
}
\]

(When \(\varepsilon_{1,j}\) vary by \(j\), write
\(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega},\{\varepsilon_{1,j}\})\); equal
caps \(\varepsilon_{1,j}\equiv\delta\) are common.)
Its complement is the **high-risk** union
\(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})^c
=\bigcup_{j=1}^{J}\Gamma_{\mathrm{low},j}(\omega_j)^{c}\).

On \(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\),

\[
\boxed{
g(\gamma;\boldsymbol{\omega}) \le h(\gamma;\boldsymbol{\omega}) \le \varepsilon_1^{\mathrm{eff}}
:= 1 - \prod_{j=1}^{J}(1-\varepsilon_{1,j}) .
}
\]

*Proof sketch.* On each \(\Gamma_{\mathrm{low},j}(\omega_j)\),
\(g_j(\gamma;\omega_j)\le G_j(\gamma;\omega_j)\le\varepsilon_{1,j}\) (Prop L1);
substitute into §3.1–§3.2. \(\square\)

Equal caps \(\varepsilon_{1,j}\equiv\delta\) give
\(\varepsilon_1^{\mathrm{eff}}=1-(1-\delta)^J\). To match a global budget
\(\varepsilon_1\), set \(\delta = 1-(1-\varepsilon_1)^{1/J}\) (not \(\varepsilon_1/J\)).

### 4.2 Main bound on the \(\gamma\)-average (optional outer step)

If \(\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})^c \mid y\bigr) \le \varepsilon_2\), then

\[
\boxed{
E_{\pi(\gamma \mid y)}[g(\gamma;\boldsymbol{\omega})]
\le
(1-\varepsilon_2)\,\varepsilon_1^{\mathrm{eff}}
+
\varepsilon_2 .
}
\]

*Proof.* Split the integral over \(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\) and
\(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})^c\); on the safe set use
\(g(\cdot;\boldsymbol{\omega}) \le h(\cdot;\boldsymbol{\omega}) \le \varepsilon_1^{\mathrm{eff}}\);
on the complement use \(g(\cdot;\boldsymbol{\omega}) \le 1\). \(\square\)

**Primary reading:** the pointwise certificate on \(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\)
is the group-level result; the boxed inequality in §4.2 is only if the
\(\beta\)-marginal probability \(\pi_{\beta \mid y}(A^c(\boldsymbol{\omega})|_\beta)\) is needed.

### 4.3 Tailoring by group

| Group profile | Knob |
|---|---|
| Large \(N_j\), central rate | Larger \(\omega_j\) (hence larger \(b_j^\star(\omega_j)\)), looser \(\varepsilon_{1,j}\) |
| Small \(N_j\), low curvature | Smaller \(\omega_j\) (hence smaller \(b_j^\star(\omega_j)\)), tighter \(\varepsilon_{1,j}\) |
| Extreme rate (\(L_j \approx N_j\)) | Tighter \(\varepsilon_{1,j}\) or larger certified safe margin |

Loosening a threshold component \(\omega_j\) (allowing larger \(\beta_j\)-tails) lowers
\(G_j(\gamma;\omega_j)\) and grows \(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\).
Optional rate cross-check: certify \(\lambda^{\star,(0)}(\boldsymbol{\omega})<\tau\).

### 4.4 Lower bound on \(\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega}) \mid y\bigr)\)

The pointwise certificate §4.1 does **not** need a \(\gamma\)-mass bound. The
optional average §4.2 does: it assumes a **high-risk budget**
\(\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})^c \mid y\bigr)\le\varepsilon_2\).
Equivalently, one needs a **lower bound**

\[
\boxed{
\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega}) \mid y\bigr)
\;\ge\;
1 - \varepsilon_2 .
}
\]

**Status.** §4.5 gives the preferred exact/bracketed evaluation; §4.4.1–§4.4.3 are
analytic shortcuts (union bound, inner interval, Laplace).

#### 4.4.1 Complement and union bound

\[
\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})^c
=
\bigcup_{j=1}^{J}
\Gamma_{\mathrm{low},j}(\omega_j)^c
=
\bigcup_{j=1}^{J}
\{\gamma : G_j(\gamma;\omega_j) > \varepsilon_{1,j}\}.
\]

Hence

\[
\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})^c \mid y\bigr)
\;\le\;
\sum_{j=1}^{J}
\pi_j^{\mathrm{high}}(\varepsilon_{1,j}),
\qquad
\pi_j^{\mathrm{high}}(\varepsilon)
:=
\int_{\{\gamma : G_j(\gamma;\omega_j) > \varepsilon\}}
\pi(\gamma \mid y)\,d\gamma,
\]

and therefore

\[
\boxed{
\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega}) \mid y\bigr)
\;\ge\;
1 - \sum_{j=1}^{J} \pi_j^{\mathrm{high}}(\varepsilon_{1,j}) .
}
\]

Each \(\pi_j^{\mathrm{high}}\) is a **one-dimensional tail integral** of the
\(\gamma\)-marginal (logit: \(Z_j(\gamma)\) not closed; quadrature or bracketing
via Lemma Z in the rate note). The union bound is **valid but often loose** when
several groups fail simultaneously with non-negligible \(\gamma\)-posterior correlation.

**Do not** use \(\prod_j \pi\bigl(\Gamma_{\mathrm{low},j}(\omega_j) \mid y\bigr)\) as a
lower bound on the intersection: \(\pi(\gamma \mid y)\) does not factorize across the
gates unless \(\gamma\)-components are independent under the marginal (generally false).

#### 4.4.2 Inner interval (scalar population \(\gamma\), \(x_j \equiv 1\))

Prop L1 gives \(G_j(\gamma;\omega_j)\) even in \(\gamma\) and typically **larger**
at large \(|\gamma|\) than near the origin (both tail \(\Phi\)-terms move toward
\(0\) or \(1\) as \(|\gamma|\to\infty\)). For \(\varepsilon_{1,j}\) in the
**feasible** range \(G_j(0;\omega_j) < \varepsilon_{1,j}\), define

\[
R_j(\varepsilon_{1,j})
:=
\sup\bigl\{\,|\gamma| : G_j(\gamma;\omega_j) \le \varepsilon_{1,j}\,\bigr\}
\;>\; 0 .
\]

Then \(\Gamma_{\mathrm{low},j}(\omega_j) \supseteq [-\,R_j,\, R_j]\) (equality when
\(G_j(\cdot;\omega_j)\) is unimodal about \(0\), the usual intercept case). With
\(R_\star := \min_j R_j(\varepsilon_{1,j})\),

\[
\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})
\;\supseteq\;
[-R_\star,\, R_\star],
\qquad
\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega}) \mid y\bigr)
\;\ge\;
\pi\bigl(|\gamma| \le R_\star \mid y\bigr)
=
1 - \pi\bigl(|\gamma| > R_\star \mid y\bigr).
\]

Certifying the right-hand side reduces to **one** symmetric tail of \(\pi(\gamma \mid y)\)
(often tighter than summing \(J\) one-sided gates). Tuning \(\boldsymbol{\omega}\) or
\(\varepsilon_{1,j}\) enlarges \(R_\star\) and the inner interval (§4.3 tradeoff).

#### 4.4.3 Posterior concentration (Laplace / deviance ellipsoid)

When \(\pi(\gamma \mid y)\) is concentrated, lower-bound mass by an **explicit subset**
\(\Gamma_{\mathrm{inner}} \subseteq \Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\).

**Mode inside the safe set.** If \(\hat\gamma := \arg\max_\gamma \pi(\gamma \mid y)\)
(or a certified mode \(\mu^\star\)) lies in \(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\),
and \(\Gamma_{\mathrm{safe}}\) contains a deviance sublevel set
\(\{\gamma : D_\pi(\gamma) \le d\}\) from `MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md`
§7.5.1, then

\[
\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega}) \mid y\bigr)
\;\ge\;
\pi\bigl(D_\pi(\gamma) \le d \mid y\bigr).
\]

Under a local Laplace approximation,
\(D_\pi(\gamma) \approx (\gamma-\mu^\star)^\top \Gamma_{\mathrm{post}}(\gamma-\mu^\star)\)
and the right-hand side is **asymptotically** \(\chi^2_q(d)\) (calibrate \(d\) to a
nominal \((1-\alpha)\) contour). This is a **lower bound on safe \(\gamma\)-mass** only
after verifying \(\{\,D_\pi \le d\,\} \subseteq \Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\)
(check \(G_j(\mu^\star;\omega_j)\le\varepsilon_{1,j}\) and curvature on the ellipsoid).

**Tight prior.** Large population precision \(\lambda_\gamma\) shrinks \(\pi(\gamma \mid y)\)
toward the prior mean; if that mean lies in \(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\)
and prior tails are Gaussian, \(\pi(\Gamma_{\mathrm{safe}} \mid y)\) is close to \(1\)
when \(\lambda_\gamma\) exceeds the likelihood's effective precision on
\(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})^c\) (see weak-prior discussion in
`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §IV.3).

#### 4.4.4 Workflow summary

| Goal | Certificate |
|---|---|
| Pointwise \(g(\gamma;\boldsymbol{\omega})\) on \(\Gamma_{\mathrm{safe}}\) | §4.1 only (no \(\gamma\)-mass bound) |
| \(\gamma\)-average \(E[g(\gamma;\boldsymbol{\omega})]\) | §4.1 + §4.5 mass on \(\Gamma_{\mathrm{safe}}\) |
| Loose analytic floor | §4.4.1 union: \(1 - \sum_j \pi_j^{\mathrm{high}}\) |
| Scalar \(\gamma\), symmetric gates | §4.4.2 + §4.5.3: \(\pi(|\gamma|\le R_\star \mid y)\) |
| Concentrated posterior | §4.4.3 inner ellipsoid + subset check |

**Preferred route:** §4.5 (exact or bracketed \(\gamma\)-marginal integral on the
gate geometry). §4.4.1–§4.4.3 remain useful **inner/outer** analytic shortcuts.

### 4.5 Exact and tight computation of \(\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega}) \mid y\bigr)\)

The object to evaluate is a **single probability** under the \(\gamma\)-marginal.
Once it is available (exactly or with certified brackets), §4.6 trades it against
\(\boldsymbol{\omega}\) and \(\varepsilon_{1,j}\).

#### 4.5.1 \(\gamma\)-marginal factorization (\(q=1\), scalar intercept RE)

Under §1.2,

\[
\pi(\gamma \mid y)
=
\frac{\widetilde\pi(\gamma \mid y)}{Z_\gamma(y)},
\qquad
\widetilde\pi(\gamma \mid y)
=
\pi(\gamma)\,
\prod_{j=1}^{J} Z_j(\gamma),
\]

\[
Z_j(\gamma)
:=
\int_{\mathbb R}
e^{\ell_j(\beta_j)}\,
\phi\bigl(\beta_j;\, x_j\gamma,\, \lambda_b^{-1}\bigr)\,d\beta_j,
\qquad
Z_\gamma(y)
:=
\int_{\mathbb R} \widetilde\pi(\gamma \mid y)\,d\gamma .
\]

All ingredients are **one-dimensional** when \(q=1\): each \(Z_j(\gamma)\) is an
integral over \(\beta_j\) only; \(Z_\gamma\) and any gate probability are integrals
over \(\gamma\) only. No closed form for \(Z_j\) under logit, but adaptive quadrature
with explicit error control is standard (`stats::integrate()` in R).

**Log domain (recommended).** Work with
\(\Lambda(\gamma) := \log \widetilde\pi(\gamma \mid y)
= \log\pi(\gamma) + \sum_j \log Z_j(\gamma)\) to avoid underflow; normalize once
via \(Z_\gamma = \int e^{\Lambda(\gamma)}d\gamma\).

#### 4.5.2 Gate geometry: invert \(G_j\) in closed form, integrate on intervals

Prop L1 gives **closed-form** \(G_j(\gamma;\omega_j)\). For scalar \(x_j\gamma\),
\(G_j\) is **even** in \(\gamma\). When \(G_j(0;\omega_j) < \varepsilon_{1,j}\),
define the **certified radius**

\[
R_j(\varepsilon_{1,j},\omega_j)
:=
\sup\bigl\{\,|\gamma| : G_j(\gamma;\omega_j) \le \varepsilon_{1,j}\,\bigr\}
\]

by **bisection on \(|\gamma|\)** (one root per group; no quadrature in this step).
Under the usual unimodality of \(G_j(\cdot;\omega_j)\) about \(0\),

\[
\Gamma_{\mathrm{low},j}(\omega_j)
=
[-R_j,\, R_j],
\qquad
\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})
=
\bigcap_{j=1}^{J} [-R_j,\, R_j]
=
[-R_\star,\, R_\star],
\qquad
R_\star := \min_{j} R_j(\varepsilon_{1,j},\omega_j).
\]

Then the target mass is **one symmetric interval integral** (not a union):

\[
\boxed{
\pi\bigl(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega}) \mid y\bigr)
=
\frac{
\int_{-R_\star}^{R_\star} \widetilde\pi(\gamma \mid y)\,d\gamma
}{
Z_\gamma(y)
}.
}
\]

**Efficiency.** Cost is \(O(J \cdot N_\beta)\) per \(\gamma\) grid point for
\(\{Z_j(\gamma)\}\) plus \(O(N_\gamma)\) for the outer integral, with
\(N_\beta, N_\gamma\) set by quadrature tolerances. Reuse \(\{Z_j(\gamma)\}\) when
sweeping \(\boldsymbol{\omega}\) or \(\varepsilon_{1,j}\) (only \(R_j\) and the
outer limits change; \(\widetilde\pi(\gamma \mid y)\) is unchanged).

**Tighter than §4.4.1.** The interval formula integrates on
\(\Gamma_{\mathrm{safe}}\) **exactly** (numerically), whereas the union bound
\(1 - \sum_j \pi_j^{\mathrm{high}}\) only uses that
\(\Gamma_{\mathrm{safe}} \supseteq [-R_\star,R_\star]\) in one direction and
\(\Gamma_{\mathrm{safe}}^c \subseteq \bigcup_j \{| \gamma| > R_j\}\) for the
complement — equality on the safe set when gates are exactly symmetric intervals.

#### 4.5.3 Exact numerical evaluation (default)

**Algorithm (scalar \(\gamma\)).**

1. For each \(j\), define `logZ_j(γ)` via `integrate()` on \(\beta_j\) (log-sum-exp
   if needed).
2. Form \(\Lambda(\gamma) = \log\pi(\gamma) + \sum_j \log Z_j(\gamma)\).
3. Compute \(Z_\gamma = \int e^{\Lambda(\gamma)}d\gamma\) on \(\mathbb R\) (same
   quadrature engine).
4. Invert each \(G_j(\cdot;\omega_j)=\varepsilon_{1,j}\) for \(R_j\); set
   \(R_\star = \min_j R_j\).
5. Return
   \(\pi(\Gamma_{\mathrm{safe}}\mid y) = Z_\gamma^{-1}\int_{-R_\star}^{R_\star}
   e^{\Lambda(\gamma)}d\gamma\).

Adaptive quadrature returns an **estimated error**; treat the evaluation as
**exact up to tolerance** \(\eta\) (e.g.\ \(\eta=10^{-8}\)) for certification
purposes. This is the recommended production path when \(q=1\).

#### 4.5.4 Certified analytic brackets (no floating-point trust)

When a fully analytic certificate is required, bracket each normalizer and
integrate the brackets.

**Per-group bounds** (rate note Lemma Z; `LOGIT_SINGLE_GROUP_SAFE_REGION.md`
Prop 1 numerator bound):

\[
Z_j^{\mathrm{lb}}(\gamma)
:=
\exp\!\Big(\ell_j(\gamma) - L_j\,\lambda_b^{-1/2}\sqrt{2/\pi}\Big),
\]

\[
Z_j^{\mathrm{ub}}(\gamma)
:=
\exp\!\Big(\ell_j(\gamma) + \tfrac{L_j^2}{2\lambda_b}\Big)
\int_{\mathbb R}
e^{L_j|\beta_j - x_j\gamma|}\,
\phi\bigl(\beta_j;\, x_j\gamma,\, \lambda_b^{-1}\bigr)\,d\beta_j
\]

(the integral is a sum of two closed-form \(\Phi\)-terms, as in Prop L1).

Pointwise on \(\gamma\),
\(Z_j^{\mathrm{lb}}(\gamma) \le Z_j(\gamma) \le Z_j^{\mathrm{ub}}(\gamma)\), hence

\[
\widetilde\pi^{\mathrm{lb}}(\gamma)
:=
\pi(\gamma)\prod_j Z_j^{\mathrm{lb}}(\gamma)
\;\le\;
\widetilde\pi(\gamma \mid y)
\;\le\;
\widetilde\pi^{\mathrm{ub}}(\gamma)
:=
\pi(\gamma)\prod_j Z_j^{\mathrm{ub}}(\gamma).
\]

For any measurable \(S \subseteq \mathbb R\),

\[
\boxed{
\frac{\int_S \widetilde\pi^{\mathrm{lb}}(\gamma)\,d\gamma}
     {\int_{\mathbb R} \widetilde\pi^{\mathrm{ub}}(\gamma)\,d\gamma}
\;\le\;
\pi(S \mid y)
\;\le\;
\frac{\int_S \widetilde\pi^{\mathrm{ub}}(\gamma)\,d\gamma}
     {\int_{\mathbb R} \widetilde\pi^{\mathrm{lb}}(\gamma)\,d\gamma}.
}
\]

Take \(S = [-R_\star, R_\star]\) for \(\Gamma_{\mathrm{safe}}\). Inner/outer
integrals are one-dimensional; outer integrals use the same \(\Phi\)-template as
Prop L1. Gap between lower and upper can be tightened by **Laplace** refinement
at \(\hat\gamma = \arg\max \Lambda(\gamma)\) on the bracket denominators only.

#### 4.5.5 Optional: \(\gamma\)-average of \(g\) without splitting \(\varepsilon_2\)

If §4.5.3 is implemented, the **optional** §4.2 split can be bypassed and

\[
E_{\pi(\gamma \mid y)}[g(\gamma;\boldsymbol{\omega})]
=
\int_{\mathbb R} g(\gamma;\boldsymbol{\omega})\,\pi(\gamma \mid y)\,d\gamma
\]

computed **directly** (outer integral), using the certified profile
\(g(\gamma;\boldsymbol{\omega}) \le h(\gamma;\boldsymbol{\omega})\) pointwise.
The two-level \((\varepsilon_1,\varepsilon_2)\) partition is then a **conservative
analytic shortcut**, not the only route.

### 4.6 Design tradeoffs: \(\boldsymbol{\omega}\), \(\varepsilon_{1,j}\), \(\varepsilon_2\)

With §4.5 in hand, \((\boldsymbol{\omega}, \{\varepsilon_{1,j}\})\) maps to
**computable** objects:

| Output | Depends on |
|---|---|
| \(R_j(\varepsilon_{1,j},\omega_j)\), \(R_\star\) | \(\boldsymbol{\omega}\), \(\varepsilon_{1,j}\) via \(G_j\) |
| \(\varepsilon_1^{\mathrm{eff}} = 1-\prod_j(1-\varepsilon_{1,j})\) | \(\varepsilon_{1,j}\) only |
| \(\pi(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\mid y)\) | §4.5 ( \(\widetilde\pi(\gamma\mid y)\) fixed; gates move) |
| \(\lambda^{\star,(0)}(\boldsymbol{\omega})=\Psi(\boldsymbol{\omega})\) | \(\boldsymbol{\omega}\) only (rate cross-check) |

#### 4.6.1 Monotonicity (how knobs move the safe set)

**Threshold \(\omega_j\)** (Lemma 3:
\(\kappa_j(\omega_j)=\lambda_b(1-\omega_j)/\omega_j\),
\(b_j^\star(\omega_j)\) from \(u(\eta_j^\star(\omega_j))=\kappa_j(\omega_j)/N_j\)):

- **Larger \(\omega_j\)** \(\Rightarrow\) smaller \(\kappa_j(\omega_j)\), larger
  \(b_j^\star(\omega_j)\), **smaller**
  tail event \(E_j(\omega_j)\), **smaller** \(G_j(\gamma;\omega_j)\) for all \(\gamma\).
- Hence **larger** \(\Gamma_{\mathrm{low},j}\), **larger** \(R_j\), **larger**
  \(\Gamma_{\mathrm{safe}}\), **smaller** \(\pi(\Gamma_{\mathrm{safe}}^c \mid y)\)
  (more \(\gamma\) mass certified safe).
- **Cost:** looser \(\beta\)-tail control on \(A^c(\boldsymbol{\omega})\) (more
  \(\beta\) allowed inside the certified box); typically **higher**
  \(\lambda^{\star,(0)}(\boldsymbol{\omega})\) (slower within-\(A\) mixing).

**Per-group cap \(\varepsilon_{1,j}\)**:

- **Larger \(\varepsilon_{1,j}\)** \(\Rightarrow\) larger \(R_j\), larger
  \(\Gamma_{\mathrm{safe}}\), larger \(\varepsilon_1^{\mathrm{eff}}\) on the safe set.
- Pure cap tightening with fixed \(\boldsymbol{\omega}\): shrinks \(\Gamma_{\mathrm{safe}}\)
  without changing the geometry of \(E_j(\omega_j)\).

**High-risk budget \(\varepsilon_2\)**:

- Not a free knob on geometry; \(\varepsilon_2\) should be set **equal to or above**
  the computed \(\pi(\Gamma_{\mathrm{safe}}^c \mid y)\) when using §4.2.
- With §4.5.5 direct integration, \(\varepsilon_2\) is **optional**.

#### 4.6.2 Constraint surface for the \(\gamma\)-average

§4.2 gives, for any certified \(\varepsilon_2 \ge \pi(\Gamma_{\mathrm{safe}}^c \mid y)\),

\[
E_{\pi(\gamma \mid y)}[g(\gamma;\boldsymbol{\omega})]
\le
(1-\varepsilon_2)\,\varepsilon_1^{\mathrm{eff}} + \varepsilon_2
=: B(\varepsilon_1^{\mathrm{eff}}, \varepsilon_2).
\]

Given a target \(B_\star\), feasible pairs satisfy

\[
\varepsilon_1^{\mathrm{eff}}
\le
\frac{B_\star - \varepsilon_2}{1-\varepsilon_2}
\qquad\text{whenever } \varepsilon_2 < B_\star .
\]

**Tight split (when §4.5.3 is exact).** Set
\(\varepsilon_2 = \pi(\Gamma_{\mathrm{safe}}^c \mid y)\) and
\(\varepsilon_1^{\mathrm{eff}}\) at its certified value on \(\Gamma_{\mathrm{safe}}\):
then \(B = \varepsilon_2 + (1-\varepsilon_2)\varepsilon_1^{\mathrm{eff}}\) is the
**least conservative** §4.2 bound for that \((\boldsymbol{\omega},\varepsilon_{1,j})\).

#### 4.6.3 Pareto tradeoff (mixing vs tail)

Modellers often want **small** \(E[g(\gamma;\boldsymbol{\omega})]\) and **small**
\(\lambda^{\star,(0)}(\boldsymbol{\omega})\) (fast mixing on \(A\)). These pull in
opposite directions on \(\boldsymbol{\omega}\):

| Pull | Direction on \(\omega_j\) |
|---|---|
| Small tail / small \(G_j\) | **Increase** \(\omega_j\) (wider certified \(\beta\) box) |
| Fast Gibbs on \(A\) | **Decrease** \(\omega_j\) (tighter \(\beta\) box, lower rate) |

**Practical sweep.** Fix design and data; grid over \(\boldsymbol{\omega}\) (or a
common \(\omega\) across groups); for each point compute \(R_\star\),
\(\pi(\Gamma_{\mathrm{safe}}\mid y)\), and \(\lambda^{\star,(0)}(\boldsymbol{\omega})\);
choose \(\varepsilon_{1,j}\) so \(G_j(\hat\gamma;\omega_j) \le \varepsilon_{1,j}\)
at the posterior mode \(\hat\gamma\) with margin. Plot the **Pareto frontier**
\(\bigl(\lambda^{\star,(0)}(\boldsymbol{\omega}),\ B(\varepsilon_1^{\mathrm{eff}},\varepsilon_2)\bigr)\)
and pick a knee point, or enforce a rate cap
\(\lambda^{\star,(0)}(\boldsymbol{\omega}) \le \tau\) and minimize
\(B\) subject to §4.5 mass.

#### 4.6.4 Equal-cap shortcut

Equal \(\varepsilon_{1,j} \equiv \delta\) gives
\(\varepsilon_1^{\mathrm{eff}} = 1-(1-\delta)^J\). To hit a global
\(\varepsilon_1^{\mathrm{eff}} = \eta\), set
\(\delta = 1-(1-\eta)^{1/J}\) (not \(\eta/J\)). Combine with §4.5 for
\(\pi(\Gamma_{\mathrm{safe}}\mid y)\) as a function of \((\boldsymbol{\omega},\delta)\).

---

## 5. Result map

| Result | Statement |
|---|---|
| \(\lambda^{\star,(0)}(\boldsymbol{\omega})\) | \(\Psi(\boldsymbol{\omega})\); rate when \(\boldsymbol{\omega}(\beta)=\boldsymbol{\omega}\), \(q=1\) |
| \(A^c(\boldsymbol{\omega})\) | \(\{(\gamma,\beta):\exists j,\ \beta_j\in E_j(\omega_j)\}\); §2 |
| \(A^c(\boldsymbol{\omega})\big|_\beta\) | \(\bigcup_j \mathcal{E}_j(\omega_j)=B(\boldsymbol{\omega})^c\) |
| \(g(\gamma;\boldsymbol{\omega})\) | \(\pi(A^c(\boldsymbol{\omega})|_\beta\mid\gamma,y)\); §1.4, §3.1 |
| \(\boldsymbol{\omega}(\beta)\) | state shrinkage vector; §1.3 |
| \(\Psi(\boldsymbol{\omega})\) | \(\boldsymbol{\rho}^\top\boldsymbol{\omega}\); rate reference when \(q=1\) |
| Lemma 1 | \(|\ell_j'|\le L_j\) |
| \(E_j(\omega_j)\) | \(\{\beta_j:\omega_j(\beta_j)\ge\omega_j\}\); §2 |
| \(g_j(\gamma;\omega_j)\) | \(\pi(E_j(\omega_j)\mid\gamma,y_j)\); §2 |
| Lemma 3 | \(E_j(\omega_j)=\{|\beta_j|>b_j^\star(\omega_j)\}\); \(\kappa_j(\omega_j)=\lambda_b(1-\omega_j)/\omega_j\) |
| Prop L1 | \(g_j(\gamma;\omega_j)\le G_j(\gamma;\omega_j)\) |
| Profile \(h\) | \(h(\gamma;\boldsymbol{\omega})\ge g(\gamma;\boldsymbol{\omega})\); §3.2 |
| Pointwise | \(g\le h\le \varepsilon_1^{\mathrm{eff}}\) on \(\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\) |
| \(\gamma\)-mass (optional) | \(\pi(\Gamma_{\mathrm{safe}} \mid y)\) via §4.5; §4.6 tradeoffs |
| Average (optional) | §4.2: \((1-\varepsilon_2)\varepsilon_1^{\mathrm{eff}}+\varepsilon_2\); or §4.5.5 direct |

**Logic.** Per-group Prop L1 \(\Rightarrow\) union identity §3.1 \(\Rightarrow\)
profile \(h(\gamma;\boldsymbol{\omega})\) \(\Rightarrow\)
per-group \(\gamma\) gates \(\Rightarrow\) §4.5 \(\gamma\)-mass \(\Rightarrow\)
optional §4.2 average or §4.5.5 direct integral; §4.6 trades \(\boldsymbol{\omega}\),
\(\varepsilon_{1,j}\), \(\varepsilon_2\).

---

## 6. What this note does **not** cover

- **\(\gamma\)-integrated** \(\widetilde{\pi}(\beta \mid y)\) and normalizer \(Z_{\mathrm{lb}}\)
  — see `LOGIT_MARGINAL_INTEGRATE_GAMMA.md` if that route is needed.
- **Full-target** pathwise Claims 2–3 and global TV without safe-slice restriction —
  §7.2 and `LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §11.
- Foster–Lyapunov §8–§11 of `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` — **invalid** for logit
  with nonempty \(A^c(\boldsymbol{\omega})\); not used in Approach A (§7.1).

---

## 7. Two program approaches

Static certification (§1–§5) supplies **geometry** and **equilibrium tail mass**
\(B(\boldsymbol{\omega})\). Two complementary ways to connect that to **sampler
convergence** differ in whether one restricts attention to the **safe slice** or
keeps the **full** Gibbs target.

### 7.1 Approach A — Large \(\boldsymbol{\omega}\), negligible tail, mixing on the safe slice

**Idea.** Choose certification thresholds \(\boldsymbol{\omega}\) large enough
(relative to design and priors) that the certified tail budget is negligible:

\[
B(\boldsymbol{\omega})
=
\varepsilon_2(\boldsymbol{\omega})
+
\bigl(1-\varepsilon_2(\boldsymbol{\omega})\bigr)\,\varepsilon_1^{\mathrm{eff}}
\;\ll\; 1,
\]

using §4.5–§4.6 (larger \(\omega_j\) \(\Rightarrow\) smaller \(E_j(\omega_j)\),
smaller \(G_j\), larger \(\Gamma_{\mathrm{safe}}\), smaller
\(\varepsilon_2(\boldsymbol{\omega})=\pi(\Gamma_{\mathrm{safe}}^c\mid y)\)). Then
analyze — or **run** — the sampler **as if the safe region were the state space**:

\[
A(\boldsymbol{\omega})\big|_\beta
=
B(\boldsymbol{\omega})
=
\{\,\beta : \beta_j \in S_j(\omega_j)\ \forall j\,\},
\qquad
S_j(\omega_j) := \{\beta_j : \omega_j(\beta_j) < \omega_j\},
\]

together with \(\gamma \in \Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\) when a
\(\gamma\)-gate is used. Equilibrium tail mass \(\pi(A^c)\) is certified small by
§4.2–§4.5; the remaining error is treated as negligible for inference or absorbed
into a user tolerance.

**Proposition 2 (Rosenthal, safe-space minorization).** Treat
\(A(\boldsymbol{\omega})\times\Gamma_{\mathrm{safe}}(\boldsymbol{\omega})\) as the
**state space** for the restricted sampler (or as the region where the chain is certified
to live after §4.5 tail control). On that set, the two-block Gibbs kernel is expected to
satisfy a **Doeblin / global minorization** on the safe space:

\[
q_\gamma(\gamma' \mid \gamma)
\;\ge\;
\varepsilon\, Q(\gamma')
\qquad
\forall\,\gamma,\gamma' \text{ in the safe slice},
\]

with fixed refresh measure \(Q\) (Gaussian refresh:
`MINORIZATION_GAUSSIAN_REFRESH.md`; construction in `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md`
§10). **Rosenthal (1995), Proposition 2** in `lmebayes/references/minor.pdf` (p. 5)
then gives, for the chain on that safe space,

\[
\boxed{
\bigl\|P^n(x,\cdot)-\pi(\cdot)\bigr\|_{TV}
\;\le\;
(1-\varepsilon)^{\lfloor n/k_0\rfloor}
}
\]

(with \(k_0=1\) when the minorization is one-step). **Applied to the safe region, this is
enough** for a certified mixing rate — no separate drift argument and **no**
rate-note §8.3 program.

**Cost (§4.6.3).** Larger \(\omega_j\) shrinks the tail but **raises**
\(\lambda^{\star,(0)}(\boldsymbol{\omega})=\Psi(\boldsymbol{\omega})\): slower
contraction **inside** \(A\). Approach A accepts that Pareto trade: pick
\(\boldsymbol{\omega}\) so \(B(\boldsymbol{\omega})\) is negligible **and** the safe
slice is large enough that minorization \(\varepsilon\) is usable.

**Implementation note.** Formally: (i) **restricted** kernel on
\(A\times\Gamma_{\mathrm{safe}}\) (Proposition 2 applies literally), or (ii) full Gibbs
with initialization in \(A\) and certified \(\pi(A^c)\approx 0\) from §4.5 (approximate
safe-space analysis). **Proposition R1** in `RESTRICTED_GIBBS_MINORIZATION_TV.md`
combines \((1-\varepsilon_d)^n\) with truncation \(\pi(C_d^c)\) in TV
(\(2\pi(C_d^c)\) in L1). **Proposition R2** (same note §7) states that under package
GLMM hypotheses, for every \(\delta>0\) there **exists** such a \(C_d\) with
\(\pi(C_d^c)<\delta\) and verified minorization — proof program, symmetric case in
`MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md`.

### 7.2 Approach B — Full target, static tail bounds, no safe-space restriction

**Idea.** Keep the **exact** Gibbs target on \((\beta,\gamma)\) and do **not**
restrict the state space. Use this note only to certify **equilibrium** tail mass
\(\pi(A^c)\) (and optional \(\gamma\)-averages §4.2) via Prop L1, profile \(h\),
and §4.5 quadrature — **without** minorization on all of \(\mathbb R\).

**TV decomposition.** Total variation is split as in `ELLIPSOID_TV_BOUND.md` /
`LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §7:

\[
\|\,P^n(x,\cdot) - \pi(\cdot)\,\|_{TV}
\;\le\;
\|P^n_A(x,\cdot) - \pi_A(\cdot)\|_{TV}
\;+\;
4\,P^n(x,A^c)
\;+\;
2\,\bigl|\pi(A^c) - P^n(x,A^c)\bigr|.
\]

The **on-\(A\)** term may still be handled heuristically
(`SAFE_UNSAFE_TV_DECOMPOSITION.md` §7.1); the
**escape** terms are controlled by the three **drift-free** claims in
`LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §11:

| Claim | Role |
|---|---|
| (1) \(\pi(A^c) \le B(\boldsymbol{\omega})\) | Equilibrium tail small (this note) |
| (2) \(P^n(x,A^c) \le \pi(A^c)\) | Path does not overshoot equilibrium tail |
| (3) \(\|\pi(A^c) - P^n(x,A^c)\|\) small | Finite-\(n\) gap closes |

Claims (2)–(3) remain **open** for logit; they must **not** use rate-note §8.3.
Approach B is the correct **full-target** program when one needs \(\pi(A^c)\) small
but cannot assume the chain stays in \(A\).

**What Approach B does not supply.** Without Doeblin minorization on the **full**
state space, Rosenthal's Proposition 2 does **not** apply globally — tail
stickiness in \(A^c\) breaks uniform minorization. Approach B certifies **how small
the tail can be made** (via \(\boldsymbol{\omega}\)), not **how fast** the
unrestricted chain mixes on all of \(\mathbb R\).

### 7.3 Comparison and when to use which

| | **Approach A** (safe space + Prop 2) | **Approach B** (full target + static tail) |
|---|---|---|
| **Primary output** | Rosenthal **Proposition 2** TV rate on safe space | Certified **\(\pi(A^c)\)** and escape terms in full TV |
| **Uses \(\boldsymbol{\omega}\)** | Large \(\omega_j\) so \(B(\boldsymbol{\omega})\approx 0\) | Same certification; tail may be small but not assumed zero |
| **Minorization** | Doeblin on entire safe slice (\(\varepsilon Q\)) | Not claimed on full \(\mathbb R\) |
| **`minor.pdf`** | **Proposition 2** (p. 5) | Not used for global rates |
| **Open work** | Certify \(\varepsilon\) on \(A(\boldsymbol{\omega})\) | Prove Claims 2–3 (no §8.3) |

**Practical recommendation.** If the goal is a **reportable MCMC burn-in / effective
sample size** guarantee on the region where the chain should live, use **Approach A**:
sweep \(\boldsymbol{\omega}\) (§4.6.3), enforce \(B(\boldsymbol{\omega}) \le \eta\), certify
minorization on the safe space, and apply **Rosenthal Proposition 2**
(`minor.pdf`). If the goal is **correctness of the exact posterior** with a
**theorem-level** bound on tail mass at equilibrium (and eventually on
\(P^n(x,A^c)\)), use **Approach B**.

Both approaches share the same **static** inputs from §1–§5 and
`CONDITIONAL_TAIL_BOUNDS_ONE_DIM_GLMM.md` (Prop 1 by likelihood); they differ only
in how those inputs connect to the **sampler**.

---

## References (in repo)

- `inst/CONDITIONAL_TAIL_BOUNDS_ONE_DIM_GLMM.md` — Prop 1 by likelihood
- `inst/RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0–§2 — Prop L1, Lemmas T/G/Z
- `inst/CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §IV.2–IV.3 — \(A/A^c\), box certification
- `inst/LOGIT_SINGLE_GROUP_SAFE_REGION.md` §4.6 — coupling template for Proposition 2
- `inst/MINORIZATION_GAUSSIAN_REFRESH.md` — refresh measure \(Q\), \(\varepsilon\)
- `inst/LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §11 — Claims 1–3 (Approach B)
- `inst/ELLIPSOID_TV_BOUND.md` — safe / unsafe TV split
- `inst/RESTRICTED_GIBBS_MINORIZATION_TV.md` — Proposition R1 (restricted kernel + \(\pi(C_d^c)\))
- `lmebayes/references/minor.pdf` — Rosenthal (1995), **Proposition 2** (Doeblin minorization)
