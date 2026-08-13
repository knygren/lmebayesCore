# Chapter-C03 §2.3: \(P\)- and \(A\)-matrix formulas by likelihood

Companion to `vignette("Chapter-C03")` §2.2–§2.3 (Claims 2–3, Remarks 8–10)
and to `.two_block_S_P11()` in `R/two_block_ergodicity.R`.

**Purpose.** Record the **full multivariate** joint precision blocks
\(P_{11},P_{12},P_{21},P_{22}\) and the rate matrix
\[
A \;=\; P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2},
\]
for the package's two-block hierarchical mixed model, **by likelihood family**.
Each likelihood section gives (i) the multivariate formulas, (ii) conditions
on group-level eigenvalues, and (iii) a scalar/diagonal specialization.

Likelihoods enter only through per-observation IRLS / Fisher weights \(w_i\)
(local Gaussian approximation at \(\eta_i\)).

**Update order (Chapter-C03 Definition 4).** Random effects \(\beta\) (block
updated **first**) \(\leftrightarrow\) paper's \(x_2\); population /
hyper parameters \(\gamma\) (updated **second**) \(\leftrightarrow\)
\(x_1\). The rate \(\lambda^\star = \lambda_{\max}(A)\) governs convergence
of the \(\gamma\)-chain.

---

## Part I — Multivariate framework (all likelihoods)

### I.1 State, dimensions, and designs

| Object | Dimension | Role |
|---|---|---|
| \(\gamma\) | \(q = \sum_{k=1}^{p_{\mathrm{re}}} q_k\) | Block 1 / hyper vector (updated second) |
| \(\beta\) | \(J p_{\mathrm{re}}\) | Stacked RE vectors \(\beta_j \in \mathbb{R}^{p_{\mathrm{re}}}\) |
| \(Z_j\) | \(n_j \times p_{\mathrm{re}}\) | RE design in group \(j\) |
| \(H_j\) | \(p_{\mathrm{re}} \times q\) | Maps \(\gamma\) to prior mean of \(\beta_j\) |
| \(P_b\) | \(p_{\mathrm{re}} \times p_{\mathrm{re}}\) | RE prior precision (Block 1) |
| \(V_k^{-1}\) | \(q_k \times q_k\) | Population prior precision, component \(k\) |

**Hyper design (`x_hyper`).** For RE component \(k\), row \(j\) of
\(X_k :=\) `x_hyper[[k]]` is a \(q_k\)-vector \(x_{j}^{(k)\top}\). Then
\((H_j)_{k,\,\mathrm{cols}(k)} = x_{j}^{(k)\top}\) and all other entries of
row \(k\) are zero. Schematically,
\[
H_j =
\begin{bmatrix}
x_j^{(1)\top} & 0 & \cdots & 0 \\
0 & x_j^{(2)\top} & \cdots & 0 \\
\vdots & & \ddots & \\
0 & \cdots & 0 & x_j^{(p_{\mathrm{re}})\top}
\end{bmatrix}
\in \mathbb{R}^{p_{\mathrm{re}} \times q}.
\]

**Linear predictors.** \(\eta_{j,i} = \mathrm{offset}_{j,i} + z_{j,i}^\top \beta_j\),
with \(z_{j,i}^\top\) row \(i\) of \(Z_j\).

### I.2 Observation weights (generic)

Package convention (`two_block_mode_weights()`):

\[
w_i \;=\; \frac{\mathrm{wt}_i\,\bigl[\mu'(\eta_i)\bigr]^2}{V(\mu_i)\,\phi},
\qquad
\eta_i = \mathrm{offset}_i + z_i^\top \beta_{j(i)}.
\]

For canonical links, \(w_i = -\partial^2 \ell_i / \partial \eta_i^2\)
(observed Hessian); otherwise Fisher information. Write
\(W_j := \mathrm{diag}(w_{j,1},\ldots,w_{j,n_j})\).

**Group data-precision Gram matrix** (the multivariate generalization of a
scalar weight sum):

\[
\mathcal{W}_j \;:=\; Z_j^\top W_j Z_j
\;=\;
\sum_{i=1}^{n_j} w_{j,i}\, z_{j,i} z_{j,i}^\top
\;\in\; \mathbb{R}^{p_{\mathrm{re}} \times p_{\mathrm{re}}}.
\]

**One-step RE conditional precision** (block \(j\) of \(P_{22}\)):

\[
B_j \;:=\; \mathcal{W}_j + P_b
\;\in\; \mathbb{R}^{p_{\mathrm{re}} \times p_{\mathrm{re}}}.
\]

### I.3 Full multivariate \(P\)-matrix (likelihood enters via \(\mathcal{W}_j\))

Stack \(\beta = (\beta_1^\top,\ldots,\beta_J^\top)^\top\). The joint
posterior precision (quadratic form of \(-\log\pi(\gamma,\beta\mid y)\) at a
working point) has blocks:

\[
\boxed{
\begin{aligned}
P_{11} &\in \mathbb{R}^{q \times q}, &
P_{11} &= \sum_{j=1}^{J} H_j^\top P_b H_j \;+\; \mathrm{blockdiag}_{k=1}^{p_{\mathrm{re}}}\!\bigl(V_k^{-1}\bigr), \\[6pt]
P_{22} &\in \mathbb{R}^{(Jp_{\mathrm{re}}) \times (Jp_{\mathrm{re}})}, &
P_{22} &= \mathrm{blockdiag}_{j=1}^{J}\!\bigl(B_j\bigr)
       \;=\; \mathrm{blockdiag}_{j=1}^{J}\!\bigl(Z_j^\top W_j Z_j + P_b\bigr), \\[6pt]
P_{12} &\in \mathbb{R}^{q \times (Jp_{\mathrm{re}})}, &
P_{12} &= \Big[\,-H_1^\top P_b \;\Big|\; \cdots \;\Big|\; -H_J^\top P_b \,\Big], \\[6pt]
P_{21} &\in \mathbb{R}^{(Jp_{\mathrm{re}}) \times q}, &
P_{21} &= P_{12}^\top.
\end{aligned}
}
\]

**Expanded \(P_{11}\) in hyper-design coordinates.** With
\(\gamma = (\gamma_1^\top,\ldots,\gamma_{p_{\mathrm{re}}}^\top)^\top\),

\[
\bigl(H_j^\top P_b H_j\bigr)_{k\ell}
\;=\;
(P_b)_{k\ell}\, x_j^{(k)} {x_j^{(\ell)}}^\top
\;\in\; \mathbb{R}^{q_k \times q_\ell},
\]
so
\[
P_{11} = \mathrm{blockdiag}(V_k^{-1})
+ \sum_{j=1}^{J}
\begin{bmatrix}
(P_b)_{1,1}\, x_j^{(1)}{x_j^{(1)}}^\top & \cdots & (P_b)_{1,p_{\mathrm{re}}}\, x_j^{(1)}{x_j^{(p_{\mathrm{re}})}}^\top \\
\vdots & \ddots & \vdots \\
(P_b)_{p_{\mathrm{re}},1}\, x_j^{(p_{\mathrm{re}})}{x_j^{(1)}}^\top & \cdots & (P_b)_{p_{\mathrm{re}},p_{\mathrm{re}}}\, x_j^{(p_{\mathrm{re}})}{x_j^{(p_{\mathrm{re}})}}^\top
\end{bmatrix}.
\]

**Coupling matrix** (Claim 2):

\[
P_{12} P_{22}^{-1} P_{21}
\;=\;
\sum_{j=1}^{J} H_j^\top \underbrace{\bigl(P_b\, B_j^{-1}\, P_b\bigr)}_{p_{\mathrm{re}} \times p_{\mathrm{re}}} H_j
\;\in\; \mathbb{R}^{q \times q}.
\]

The package also forms \(S := P_{12}P_{22}^{-1}P_{21}\) and passes \((S,P_{11})\)
to `.two_block_gen_eigen()` (equivalent to eigendecomposing \(A\)).

### I.4 Full multivariate \(A\)-matrix (Chapter-C03 §2.3)

\[
\boxed{
A \;=\; P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2}
\;\in\; \mathbb{R}^{q \times q},
\qquad
\lambda^\star \;:=\; \lambda_{\max}(A).
}
\]

Equivalent symmetric form used in code: with \(R = \mathrm{chol}(P_{11})\),

\[
A \;=\; R^{-\top}\, S\, R^{-1},
\qquad S = P_{12} P_{22}^{-1} P_{21}.
\]

**Mean contraction (Claim 2).**
\[
E[\gamma^{(l)} \mid \gamma^{(0)}] - \mu_\gamma
\;=\;
\bigl(P_{11}^{-1/2} A^l P_{11}^{1/2}\bigr)\bigl(\gamma^{(0)} - \mu_\gamma\bigr).
\]

**Variance iteration (Claim 3).**
\[
\Sigma^{(l)}_{11}
\;=\;
P_{11}^{-1/2}\left[\sum_{i=1}^{2l} A^{\,i-1}\right]P_{11}^{-1/2},
\qquad
\Sigma_{11}^{-1/2}\,\Sigma^{(l)}_{11}\,\Sigma_{11}^{-1/2} = I - A^{2l}.
\]

Geometric ergodicity of the \(\gamma\)-chain requires \(\rho(A) < 1\), i.e.
\(\lambda^\star < 1\).

### I.5 Invertibility of \(P_{11}\) and \(P_{22}\)

**\(P_{22}\).** Each block \(B_j = \mathcal{W}_j + P_b\) satisfies
\(\mathcal{W}_j \succeq 0\) and \(P_b \succ 0\), hence \(B_j \succ 0\).
Therefore \(P_{22} \succ 0\) and \(P_{22}^{-1}\) exists **unconditionally**
for every likelihood.

**\(P_{11}\).** With \(V_k^{-1} \succ 0\) and \(P_b \succ 0\),
\(P_{11} \succ 0\) provided the hyper design is estimable (the sum
\(\sum_j H_j^\top P_b H_j\) is not rank-deficient relative to the prior).

### I.6 Mixing rate vs.\ other curvature constraints

| Question | Matrix object | Bad direction | Typical cause |
|---|---|---|---|
| **Mixing rate** \(\lambda^\star < 1\) | Rate matrix \(A\) via \(P_{22}^{-1}\) | \(\lambda^\star \uparrow 1\) | \(\mathcal{W}_j\) **too small** (data precision \(\to 0\)) for many groups — \(B_j \approx P_b\), large \(P_{22}^{-1}\), strong coupling |

Part IV below develops the **rate partition** \(A / A^c\) from \(\lambda^\star\)
(average shrinkage too high). Pilot eigenvalue **ceilings** on \(\mathcal{W}_j\)
(too large) appear in `ELLIPSOID_TV_BOUND.md` and are **not** part of that
partition.

---

## Part II — Likelihood-specific formulas

Each section below repeats the **full multivariate** blocks from Part I with
that likelihood's \(w_i\), then states multivariate eigenvalue conditions,
then specializes.

---

## 1. Normal (Gaussian) regression

### 1.1 Weights

Homoscedastic Gaussian, dispersion \(\phi = \sigma^2\):

\[
w_i = \frac{1}{\phi}
\quad\text{(constant, independent of \(\eta\))}.
\]

Curvature: \(c(\eta) = 1/\phi\) (after scaling).

### 1.2 Multivariate \(P\)-blocks

\[
\mathcal{W}_j = \frac{1}{\phi}\, Z_j^\top Z_j,
\qquad
B_j = \frac{1}{\phi}\, Z_j^\top Z_j + P_b,
\]

\[
P_{11} = \sum_{j=1}^{J} H_j^\top P_b H_j + \mathrm{blockdiag}(V_k^{-1}),
\quad
P_{22} = \mathrm{blockdiag}_j(B_j),
\quad
P_{12} = \big[-H_1^\top P_b \mid \cdots \mid -H_J^\top P_b\big],
\quad
P_{21} = P_{12}^\top,
\]

\[
A = P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2}.
\]

### 1.3 Multivariate eigenvalue conditions

**Mixing (\(\lambda^\star < \tau\)).** No \(\eta\)-dependent collapse of
\(\mathcal{W}_j\); failure only if \(\phi \to \infty\). Require
\(\lambda_{\max}(A) < \tau\) at the working point (or via a comparison bound
with finite \(\phi\)).

**Safe region (ceiling).** \(\lambda_{\max}(Z_j^\top Z_j / \phi) \le
M_j^{\mathrm{pilot}}\) for all \(j\); vacuous for fixed \(\phi\) unless
dispersion is state-dependent.

### 1.4 Specialization: scalar \(P_{11}\), diagonal \(P_{22}\)

One RE column per group (\(p_{\mathrm{re}}=1\)), one hyper parameter
(\(q=1\)), intercept RE design \(z_{j,i}=1\), scalar hyper row \(H_j = x_j\),
scalar RE prior \(P_b = \lambda_b\), population precision \(V^{-1}=\lambda_\gamma\).

**Blocks.** Data precision is constant:
\(W_j = \sum_{i\in j} w_{j,i} = n_j/\phi\), independent of \(\beta\) and
\(\eta\). The joint precision specializes to
\[
P_{11} = \lambda_b \sum_{j=1}^{J} x_j^2 + \lambda_\gamma,
\quad
P_{22,jj} = \lambda_b + \frac{n_j}{\phi},
\quad
P_{12,j} = -\lambda_b x_j,
\quad
P_{21} = P_{12}^\top.
\]
Each \(P_{22,jj}\) is the total precision on \(b_j\) in the Gaussian
working model: RE prior \(\lambda_b\) plus \(n_j/\phi\) from the likelihood.

**Rate (substitute into §I.4).** Because \(P_{11}\) is scalar,
\(P_{11}^{-1/2} = P_{11}^{-1/2}\) and
\[
(P_{12} P_{22}^{-1} P_{21})_{11}
\;=\;
\sum_{j=1}^{J} \frac{(\lambda_b x_j)^2}{\lambda_b + n_j/\phi}
\;=\;
\lambda_b^2 \sum_{j=1}^{J} \frac{x_j^2}{\lambda_b + n_j/\phi},
\]
hence
\[
\boxed{
A \;=\;
\frac{\lambda_b^2}{P_{11}}
\sum_{j=1}^{J} \frac{x_j^2}{\lambda_b + n_j/\phi},
\qquad
\lambda^\star = A.
}
\]
Distributing \(\lambda_b^2 = \lambda_b\cdot(\lambda_b x_j^2)\) inside the sum
(as in §2.4) makes the numerator a **weighted sum of group coupling terms**
\(\lambda_b(\lambda_b x_j^2)/(\lambda_b + n_j/\phi)\), each divided by the
same \(P_{11}\).

**Group shrinkage weights (Normal template).** Write
\[
\omega_j
\;:=\;
\frac{\lambda_b\, x_j^2}{\lambda_b + n_j/\phi}
\;\in\;(0, x_j^2],
\]
the fraction of hyper–RE coupling that remains **prior-dominated** for group
\(j\) once data precision \(n_j/\phi\) is included. When \(n_j/\phi \ll
\lambda_b\), \(\omega_j \approx x_j^2\); when \(n_j/\phi \gg \lambda_b\),
\(\omega_j \approx 0\) and group \(j\) barely pulls on \(\gamma\). Then
\[
A \;=\; \frac{\lambda_b}{P_{11}} \sum_{j=1}^{J} \omega_j
\;=\;
\frac{\lambda_b\, J\,\bar\omega_x}{P_{11}},
\qquad
\bar\omega_x := \frac{1}{J}\sum_{j=1}^{J} \omega_j.
\]
The rate is **mean group shrinkage** times \(\lambda_b/P_{11}\): data pull
groups off \(\gamma\) (small \(\omega_j\)) and a tight \(\gamma\) prior (large
\(\lambda_\gamma\)) adds global contraction through \(P_{11}\). This is the
simplest instance of the \(\omega_j\) / \(\bar\omega\) picture in §IV.2; all
later likelihoods replace \(n_j/\phi\) by state-dependent \(W_j(\beta)\).

**Special case \(x_j \equiv 1\).** Then \(P_{11} = J\lambda_b + \lambda_\gamma\),
\(\omega_j = \lambda_b/(\lambda_b + n_j/\phi)\), and
\[
\lambda^\star
\;=\;
\frac{\lambda_b^2}{J\lambda_b + \lambda_\gamma}
\sum_{j=1}^{J} \frac{1}{\lambda_b + n_j/\phi}
\;=\;
\frac{\lambda_b}{J\lambda_b + \lambda_\gamma}
\sum_{j=1}^{J} \omega_j.
\]
If additionally \(n_j \equiv n\),
\[
\lambda^\star
\;=\;
\frac{J\lambda_b^2}{(J\lambda_b + \lambda_\gamma)\bigl(\lambda_b + n/\phi\bigr)}.
\]
Flat \(\gamma\) (\(\lambda_\gamma = 0\)): \(\lambda^\star =
\lambda_b/(\lambda_b + n/\phi) = \bar\omega\), the familiar **one-group
shrinkage factor** \(\lambda_b/(\lambda_b + \text{data precision})\), unchanged
when replicated across \(J\) identical groups — averaging does not dilute the
rate when every \(\omega_j\) is equal.

**Safe set \(A\) and unsafe \(A^c\) (§IV).** Because \(W_j = n_j/\phi\) does
**not** depend on \(\beta\), \(\lambda^\star\) is a **constant** for the whole
model at fixed \((\phi, \lambda_b, \lambda_\gamma, \{n_j, x_j\})\). The rate
partition is degenerate:
\[
\text{either all }\beta \in A
\quad\text{or}\quad
\text{all }\beta \in A^c,
\]
with no nontrivial geometry in \(\beta\)-space. Certification is a **design
and prior tuning** problem: increase \(n_j/\phi\), decrease \(\lambda_b\), or
tighten \(\lambda_\gamma\) until \(\lambda^\star < \tau\). There is no
separation tail to drive \(A^c\); the only mixing failure is \(\phi \to \infty\)
(data precision vanishes, \(\omega_j \uparrow\), \(\lambda^\star \to 1\)).

**Equal sample sizes, balanced hyper design (\(x_j \equiv 1\), \(n_j \equiv
n\)).** The certification inequality \(\lambda^\star < \tau\) becomes
\[
\frac{J\lambda_b}{J\lambda_b + \lambda_\gamma}\cdot
\frac{\lambda_b}{\lambda_b + n/\phi}
\;<\;
\tau,
\]
or, with flat \(\gamma\), simply \(\lambda_b/(\lambda_b + n/\phi) < \tau\).
Solving for sample size at fixed priors:
\[
\frac{n}{\phi}
\;>\;
\lambda_b\left(\frac{1}{\tau} - 1\right)
\quad\text{(flat }\gamma\text{, }\lambda_b < \tau\text{)}.
\]
Each group must carry enough Gaussian information (\(n_j/\phi\)) relative to
\(\lambda_b\) that the **average shrinkage** \(\bar\omega\) is below
\(\tau/\lambda_b\).

**Choosing \(\tau\) (certification threshold).** The target rate
\(\tau \in (0,1)\) is a **modeller's choice**, not a fixed constant (we use
\(\tau = 0.95\) only as a running example). It enters the safe set
\(A = \{\beta : \lambda^\star(\beta) < \tau\}\); the **achieved** rate
\(\lambda^\star\) is computed from \((\phi, \lambda_b, \lambda_\gamma,
\{n_j, x_j\})\). In the Normal scalar case these decouple:

1. **All \(\beta\) automatically safe or unsafe.** Because \(\lambda^\star\)
   does not depend on \(\beta\), either every working point lies in \(A\) or
   every one lies in \(A^c\). There is no partial geometry to navigate.

2. **Set \(\tau\) at the threshold.** After fitting priors and design, compute
   \(\lambda^\star\) once (e.g.\ `two_block_rate()`). A common Normal/Gamma
   practice is to choose \(\tau\) **just above** the achieved rate,
   \(\tau = \lambda^\star + \varepsilon\) with small \(\varepsilon > 0\), so
   the certification inequality is satisfied with **every** \(\beta \in A\)
   and Remark-8 applies globally on \(\beta\)-space. Equivalently: calibrate
   \((n_j/\phi, \lambda_b, \lambda_\gamma)\) so that the computed
   \(\lambda^\star\) sits below a pre-chosen \(\tau\) (e.g.\ \(0.95\)).

3. **Robust to a weak \(\gamma\) prior.** As \(\lambda_\gamma \downarrow\),
   \(P_{11}\) loses population anchoring and \(\lambda^\star \uparrow\)
   (§2.4). For \(x_j \equiv 1\), the **worst case** over data and
   \(\lambda_\gamma \ge 0\) at fixed \(\lambda_b\) is \(\lambda_\gamma = 0\)
   with minimal data precision \(n_j/\phi = 0\), giving \(\omega_j = 1\) and
   \(\lambda^\star = \lambda_b\). To certify **all** \(\beta \in A\) even when
   the \(\gamma\) prior is taken weak, choose \(\tau > \lambda_b\) (or, more
   sharply, \(\tau > \lambda^\star(\lambda_\gamma = 0)\) evaluated at the
   actual \(\{n_j\}\)). With adequate data,
   \(\lambda^\star(\lambda_\gamma = 0) = \lambda_b/( \lambda_b + n/\phi)\)
   (equal \(n_j\)) and a fixed \(\tau = 0.95\) already works once
   \(n/\phi > \lambda_b(1/\tau - 1)\) — no tight \(\gamma\) prior is needed.

4. **What cannot be done here but matters elsewhere.** Because \(\lambda^\star\)
   is constant, Normal offers the **only** likelihood in this note where one
   can make \(A^c = \varnothing\) by prior/design/threshold choice alone.
   GLMs have \(\lambda^\star(\beta)\) and \(A^c\) depending on \(\beta\); see
   §IV.5.

**Why this case still matters.** Normal is the **reference template** for
§2.4–§6 and Part IV: same \(P\)-block layout, but \(W_j\) replaced by
likelihood-specific sums. Logit/Poisson inherit the shrinkage-weight algebra
with \(W_j(\beta)\) in place of \(n_j/\phi\); only Normal makes \(A/A^c\) a
**global** (non-geometric) statement. Package evaluation:
`two_block_rate()` at the mode returns this \(\lambda^\star\) unchanged under
any \(\beta\)-update on a Gaussian fit because `two_block_mode_weights()` are
constant.

---

## 2. Logit (binomial, canonical logit link)

### 2.1 Weights

\(\mu_i = \mathrm{expit}(\eta_i)\), binomial count \(n_i\):

\[
w_i = n_i\, p_i(1-p_i),
\qquad
p_i = \mu_i,
\qquad
0 < w_i \le \frac{n_i}{4}.
\]

Curvature function: \(c(\eta) = n_i\, p(\eta)(1-p(\eta))\) per observation
(contributes to \(\mathcal{W}_j\) via \(w_{j,i} = c(\eta_{j,i})\) when
\(z_{j,i}\) is scalar; more generally \(w_{j,i}\) scales the rank-one piece
\(z_{j,i}z_{j,i}^\top\)).

### 2.2 Multivariate \(P\)-blocks

\[
\mathcal{W}_j = \sum_{i=1}^{n_j} n_{j,i}\, p_{j,i}(1-p_{j,i})\, z_{j,i} z_{j,i}^\top,
\qquad
B_j = \mathcal{W}_j + P_b,
\]

\[
P_{11} = \sum_{j} H_j^\top P_b H_j + \mathrm{blockdiag}(V_k^{-1}),
\quad
P_{22} = \mathrm{blockdiag}_j(B_j),
\quad
P_{12} = \big[-H_1^\top P_b \mid \cdots \mid -H_J^\top P_b\big],
\quad
P_{21} = P_{12}^\top,
\]

\[
P_{12} P_{22}^{-1} P_{21} = \sum_{j=1}^{J} H_j^\top P_b B_j^{-1} P_b H_j,
\qquad
A = P_{11}^{-1/2}\,(P_{12} P_{22}^{-1} P_{21})\,P_{11}^{-1/2}.
\]

### 2.3 Multivariate eigenvalue conditions

**Mixing (\(\lambda^\star < \tau\)).** Requires \(\lambda_{\max}(A) < \tau\).
This is threatened when **many groups** have \(\mathcal{W}_j \approx 0\)
(separation: \(| \eta_{j,i}|\to\infty\), \(p(1-p)\to 0\)), so
\(B_j \approx P_b\) and coupling through \(P_b B_j^{-1} P_b \approx P_b\)
is large. Sufficient guard (loose): \(\lambda_{\min}(\mathcal{W}_j) \ge \delta > 0\)
for all or most \(j\).

**Safe region (pilot ceiling on data curvature).** With
\(C_j(\beta_j) := Z_j^\top \mathrm{diag}\bigl(n_{j,i} p_{j,i}(1-p_{j,i})\bigr) Z_j\),

\[
\mathcal{C}_{M,j} = \bigl\{\beta_j : \lambda_{\max}(C_j(\beta_j)) \le M_j^{\mathrm{pilot}}\bigr\}.
\]

Because \(p(1-p)\) is symmetric and unimodal at \(\eta=0\), **exceeding a
ceiling** and **falling below a floor** are different; pilot-ceiling
\((\mathcal{C}_M)^c\) for **too-large** curvature is bounded (\(w_i \le n_i/4\)). The
**mixing** pathology (\(\mathcal{W}_j\) too small) is a **floor** problem on
**both tails**: \(| \eta_{j,i} |\) large \(\Rightarrow\) \(w_{j,i}\to 0\).

### 2.4 Specialization: scalar \(P_{11}\), diagonal \(P_{22}\)

\[
P_{11} = \sum_{j=1}^{J} \lambda_b x_j^2 + \lambda_\gamma,
\quad
P_{22,jj} = \lambda_b + W_j,
\quad
P_{12,j} = -\lambda_b x_j,
\qquad
W_j = \sum_{i\in j} n_{j,i}\, p_{j,i}(1-p_{j,i}).
\]

Substituting \(P_{11} = \sum_j \lambda_b x_j^2 + \lambda_\gamma\) into the scalar rate formula and distributing \(\lambda_b^2 = \lambda_b \cdot (\lambda_b x_j^2)\) inside the sum:

\[
A \;=\;
\frac{\displaystyle\sum_{j=1}^{J} \frac{\lambda_b\,\bigl(\lambda_b x_j^2\bigr)}{\lambda_b + W_j}}
     {\displaystyle\sum_{j=1}^{J} \lambda_b x_j^2 + \lambda_\gamma},
\qquad
\lambda^\star = A.
\]

The two sums are not interchangeable — the denominator is the explicit \(P_{11}\) with \(\lambda_b\) written inside the summation; the numerator applies the same \(\lambda_b\,(\lambda_b x_j^2)\) factor to each group’s coupling term \(x_j^2/(\lambda_b + W_j)\).

**Weak population prior (\(\lambda_\gamma \to 0\)).** The term \(\lambda_\gamma\) is the prior precision on \(\gamma\) (large \(\lambda_\gamma\) = tight prior, small = weak / diffuse). As \(\lambda_\gamma \downarrow\), the denominator \(P_{11} = \sum_j \lambda_b x_j^2 + \lambda_\gamma\) is dominated by the RE–hyper coupling sum and **shrinks**, while the numerator is unchanged at fixed \((\lambda_b, W_j, x_j)\). Hence \(A \uparrow\) and \(\lambda^\star \to 1\) from below: the \(\gamma\)-block update is less anchored by the prior, Gibbs coupling through \(P_{12}P_{22}^{-1}P_{21}\) carries more relative weight, and the certified contraction factor weakens (slow mixing). In the limit \(\lambda_\gamma = 0\) (improper flat prior on \(\gamma\), if defined),
\[
A \;\to\;
\frac{\displaystyle\sum_{j=1}^{J} \lambda_b x_j^2 / (\lambda_b + W_j)}
     {\displaystyle\sum_{j=1}^{J} x_j^2},
\]
provided \(\sum_j x_j^2 > 0\). If additionally \(x_j \equiv 1\),
\[
A \;\to\;
\frac{\lambda_b}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}
\;=\;
\frac{1}{J}\sum_{j=1}^{J} \frac{\lambda_b}{\lambda_b + W_j},
\]
which can still be \(< 1\) when the \(W_j\) are not all tiny (§IV.2), but is
**larger** than the same expression with any \(\lambda_\gamma > 0\) in the
denominator. A weak \(\gamma\) prior therefore never improves the rate bound;
it removes one source of contraction and pushes \(\lambda^\star\) toward the
boundary where Remark-8 certification fails (\(\lambda^\star = 1\)).

**Data precision requirement (\(W_j\) vs.\ \(\lambda_b\)).** Fast convergence requires \(\lambda^\star = A\) to be **meaningfully** below \(1\), not merely infinitesimally below it. Inspect each summand: \(\lambda_b(\lambda_b x_j^2)/(\lambda_b + W_j)\). When \(W_j \ll \lambda_b\), this is \(\approx \lambda_b x_j^2\) — the data contribution to \(P_{22,jj}\) is negligible and the group behaves as if only the RE prior \(\lambda_b\) were present. If **every** \(W_j\) is small relative to \(\lambda_b\), then
\[
\sum_{j=1}^{J} \frac{\lambda_b(\lambda_b x_j^2)}{\lambda_b + W_j}
\;\approx\;
\sum_{j=1}^{J} \lambda_b x_j^2,
\]
so with \(\lambda_\gamma = 0\) and \(x_j \equiv 1\) one obtains \(A \approx
\lambda_b\); with small \(\lambda_\gamma > 0\), \(A\) is only slightly below
that. **At least one** group must carry \(W_j\) that is **meaningfully large
compared with \(\lambda_b\)** — i.e.\ \(W_j \gtrsim \lambda_b\) or larger — so
that some factors \((\lambda_b + W_j)^{-1}\) materially shrink the numerator
and pull \(A\) away from \(\lambda_b\). For logit, that means at least one
group needs enough binomial information (large \(n_{j,i}\, p(1-p)\) at the
working \(\eta\)); groups in separation tails (\(W_j \approx 0\)) cannot
supply this by themselves.

---

## 3. Probit (binomial, probit link)

### 3.1 Weights

\(\mu_i = \Phi(\eta_i)\) (standard normal cdf):

\[
w_i = n_i \,\frac{\phi(\eta_i)^2}{\Phi(\eta_i)\,\Phi(-\eta_i)},
\qquad
w_i \le n_i \cdot \frac{2}{\pi} \approx 0.6366\, n_i.
\]

No uniform lower bound: \(w_i \to 0\) as \(\eta_i \to -\infty\).

### 3.2 Multivariate \(P\)-blocks

Same layout as §2.2 with

\[
\mathcal{W}_j = \sum_{i=1}^{n_j} w_{j,i}(\eta_{j,i})\, z_{j,i} z_{j,i}^\top,
\qquad
B_j = \mathcal{W}_j + P_b,
\]

and \(P_{11},P_{12},P_{21},P_{22},A\) as in §2.2.

### 3.3 Multivariate eigenvalue conditions

**Mixing.** \(\lambda^\star \uparrow 1\) when \(\lambda_{\min}(\mathcal{W}_j)\to 0\)
for many groups (Mills-ratio tail, typically \(\eta_{j,i}\) strongly
negative for \(y=1\) coding).

**Safe region.** Ceiling
\(\lambda_{\max}(\mathcal{W}_j) \le M_j^{\mathrm{pilot}}\): probit
\(w_i\) is bounded above, so **too-large** eigenvalues are less extreme than
Poisson; the **binding** safe-region issue is usually **too-small**
\(\mathcal{W}_j\) (left tail), which is **approximately one-sided** in
\(\eta\)-space.

### 3.4 Specialization: scalar \(P_{11}\), diagonal \(P_{22}\)

Same as §2.4 with \(W_j = \sum_{i\in j} w_{j,i}(\eta_{j,i})\).

---

## 4. Complementary log-log (binomial, cloglog link)

### 4.1 Weights

\(\mu_i = 1 - \exp(-\exp(\eta_i))\):

\[
w_i = n_i \,\frac{\exp(2\eta_i)\,(1-\mu_i)}{\mu_i}.
\]

As \(\eta_i \to -\infty\): \(w_i \to 0\). As \(\eta_i \to +\infty\):
\(w_i \to n_i\).

### 4.2 Multivariate \(P\)-blocks

Same multivariate structure as §2.2 with cloglog \(w_{j,i}\) in
\(\mathcal{W}_j = \sum_i w_{j,i}\, z_{j,i} z_{j,i}^\top\).

### 4.3 Multivariate eigenvalue conditions

**Mixing.** \(\mathcal{W}_j \to 0\) on the **lower** \(\eta\)-tail (one-sided).

**Safe region (ceiling — eigenvalues too large).**
\(\lambda_{\max}(\mathcal{W}_j) > M_j^{\mathrm{pilot}}\) occurs on the
**upper** \(\eta\)-tail (one-sided \((\mathcal{C}_M)^c\)). Lower tail is the mixing
problem, not the ceiling problem.

### 4.4 Specialization: scalar \(P_{11}\), diagonal \(P_{22}\)

Same as §2.4 with cloglog \(W_j\).

---

## 5. Poisson (log link, counts \(y_i > 0\))

### 5.1 Weights

\(\mu_i = e^{\eta_i}\), canonical log link, dispersion \(\phi=1\):

\[
w_i = e^{\eta_i} = \mu_i,
\qquad
c(\eta) = e^{\eta}.
\]

**Unbounded above** as \(\eta \to +\infty\); **vanishes** as \(\eta \to -\infty\).

### 5.2 Multivariate \(P\)-blocks

\[
\mathcal{W}_j = \sum_{i=1}^{n_j} e^{\eta_{j,i}}\, z_{j,i} z_{j,i}^\top,
\qquad
B_j = \mathcal{W}_j + P_b,
\]

\[
P_{11} = \sum_{j=1}^{J} H_j^\top P_b H_j + \mathrm{blockdiag}(V_k^{-1}),
\quad
P_{22} = \mathrm{blockdiag}_j(B_j),
\quad
P_{12} = \big[-H_1^\top P_b \mid \cdots \mid -H_J^\top P_b\big],
\quad
P_{21} = P_{12}^\top,
\]

\[
P_{12} P_{22}^{-1} P_{21} = \sum_{j=1}^{J} H_j^\top P_b B_j^{-1} P_b H_j,
\qquad
A = P_{11}^{-1/2}\,(P_{12} P_{22}^{-1} P_{21})\,P_{11}^{-1/2}.
\]

### 5.3 Multivariate eigenvalue conditions

**Safe region (ceiling — eigenvalues too large; your primary concern).**

\[
\mathcal{C}_{M,j} = \bigl\{\beta_j : \lambda_{\max}(\mathcal{W}_j(\beta_j)) \le M_j^{\mathrm{pilot}}\bigr\},
\qquad
\mathcal{C}_M = \bigcap_j \mathcal{C}_{M,j}.
\]

Because \(w_i = e^{\eta_i}\) grows without bound as \(\eta_{j,i} \to +\infty\),
\((\mathcal{C}_M)^c\) is **one-sided**: large positive linear predictors / coefficients
push \(\lambda_{\max}(\mathcal{W}_j)\) above the ceiling. For scalar
\(z_{j,i}=1\), \(\lambda_{\max}(\mathcal{W}_j) = \sum_i e^{\eta_{j,i}}\).

**Mixing (\(\lambda^\star < \tau\)) — opposite tail.** Here \(\mathcal{W}_j \to 0\)
when \(\eta_{j,i} \to -\infty\) (low rate). Then \(B_j \approx P_b\),
coupling grows, \(\lambda^\star \uparrow 1\). This is **not** the ceiling
problem; it is the low-precision pathology of §I.6.

**Ensuring \(\lambda_{\max}(A) < \tau\) (multivariate).** At a working point,
compute \(A\) explicitly; sufficient structural condition: keep
\(\lambda_{\min}(B_j)\) not too close to \(\lambda_{\min}(P_b)\) for too
many groups simultaneously — equivalently, avoid widespread \(\eta \ll 0\)
where \(e^{\eta}\) is negligible.

### 5.4 Specialization: scalar \(P_{11}\), diagonal \(P_{22}\)

\(p_{\mathrm{re}}=1\), \(q=1\), \(z_{j,i}=1\),

\[
P_{11} = \sum_{j=1}^{J} \lambda_b x_j^2 + \lambda_\gamma,
\quad
P_{22,jj} = \lambda_b + W_j,
\quad
P_{12,j} = -\lambda_b x_j,
\qquad
W_j = \sum_{i\in j} e^{\eta_{j,i}}.
\]

As in §2.4,

\[
A \;=\;
\frac{\displaystyle\sum_{j=1}^{J} \frac{\lambda_b\,\bigl(\lambda_b x_j^2\bigr)}{\lambda_b + W_j}}
     {\displaystyle\sum_{j=1}^{J} \lambda_b x_j^2 + \lambda_\gamma},
\qquad
\lambda^\star = A.
\]

**Weak population prior (\(\lambda_\gamma \to 0\)).** Identical mechanism to §2.4: a smaller \(\lambda_\gamma\) shrinks \(P_{11}\) without changing the numerator, so \(A \uparrow\) and \(\lambda^\star \to 1\). In the flat limit \(\lambda_\gamma = 0\),
\[
A \;\to\;
\frac{\displaystyle\sum_{j=1}^{J} \lambda_b x_j^2 / (\lambda_b + W_j)}
     {\displaystyle\sum_{j=1}^{J} x_j^2}.
\]

**Data precision requirement (\(W_j\) vs.\ \(\lambda_b\)).** Again each summand is \(\lambda_b(\lambda_b x_j^2)/(\lambda_b + W_j)\). For **meaningful** fast convergence (\(A\) well below \(1\)), at least one \(W_j\) must be **meaningfully large relative to \(\lambda_b\)**. For Poisson,
\[
W_j = \sum_{i\in j} e^{\eta_{j,i}} \;\to\; 0
\quad\text{as all } \eta_{j,i} \to -\infty
\]
(low expected counts / linear predictors in the **lower** tail). If every group sits in that tail, \(A \approx 1\) (exactly when \(\lambda_\gamma = 0\)). At least one group must have some observations with \(\eta_{j,i}\) not strongly negative — so that \(e^{\eta_{j,i}}\) is of order \(\lambda_b\) or larger — to pull \(A\) away from \(1\).

**Opposite tail (Poisson-specific tension).** Unlike logit, \(W_j\) is **unbounded above** as \(\eta_{j,i} \to +\infty\). Large \(W_j\) **helps** the mixing bound (\(\lambda_b + W_j\) in the denominator shrinks each summand), but simultaneously drives
\(\lambda_{\max}(\mathcal{W}_j) = W_j\) above any fixed pilot ceiling — the **one-sided** region \((\mathcal{C}_M)^c\) for eigenvalues **too large** (§5.3). So Poisson presents two distinct tail problems on **opposite** sides of the parameter line: **lower** \(\eta\) ( \(W_j \ll \lambda_b\) ) threatens \(\lambda^\star \approx 1\) (membership in \(A^c\)); **upper** \(\eta\) ( \(W_j \gg M_j^{\mathrm{pilot}}\) ) threatens \((\mathcal{C}_M)^c\), not the rate partition itself.

---

## 6. Gamma regression (log link)

### 6.1 Weights

\(\mu_i = e^{\eta_i}\), \(V(\mu_i) = \phi\,\mu_i^2\):

\[
w_i = \frac{(\mu_i')^2}{V(\mu_i)} = \frac{e^{2\eta_i}}{\phi\, e^{2\eta_i}} = \frac{1}{\phi}
\quad\text{(constant)}.
\]

### 6.2 Multivariate \(P\)-blocks

\[
\mathcal{W}_j = \frac{1}{\phi}\, Z_j^\top Z_j,
\qquad
B_j = \frac{1}{\phi}\, Z_j^\top Z_j + P_b,
\]

with \(P_{11},P_{12},P_{21},P_{22},A\) as in §1.2 (same structure as Normal).

### 6.3 Multivariate eigenvalue conditions

Identical to Normal §1.3: no \(\eta\)-dependent ceiling or floor in \(w_i\)
for log link. Failure modes only through \(\phi \to \infty\) or design rank.

### 6.4 Specialization: scalar \(P_{11}\), diagonal \(P_{22}\)

Same as §1.4 with \(n_j/\phi\) in place of \(n_j/\sigma^2\).

---

## Part III — Summary

### III.1 Multivariate \(P\)-blocks (all likelihoods)

| Block | Size | Formula |
|---|---|---|
| \(P_{11}\) | \(q \times q\) | \(\sum_j H_j^\top P_b H_j + \mathrm{blockdiag}(V_k^{-1})\) |
| \(P_{22}\) | \(Jp_{\mathrm{re}} \times Jp_{\mathrm{re}}\) | \(\mathrm{blockdiag}_j(Z_j^\top W_j Z_j + P_b)\) |
| \(P_{12}\) | \(q \times Jp_{\mathrm{re}}\) | \(\big[-H_1^\top P_b \mid \cdots \mid -H_J^\top P_b\big]\) |
| \(P_{21}\) | \(Jp_{\mathrm{re}} \times q\) | \(P_{12}^\top\) |
| \(A\) | \(q \times q\) | \(P_{11}^{-1/2} P_{12} P_{22}^{-1} P_{21} P_{11}^{-1/2}\) |

Likelihood enters only through \(W_j = \mathrm{diag}(w_{j,i})\).

### III.2 Per-likelihood weights and safe-region geometry

| Likelihood | \(w_i\) | Rate \(A^c\) (mixing; \(\mathcal{W}_j\) **too small**) |
|---|---|---|
| Normal | \(1/\phi\) | only if \(\phi \to \infty\) |
| Logit | \(n_i p(1-p)\) | both tails, separation |
| Probit | \(n_i \phi^2/(\Phi\Phi)\) | mostly lower tail |
| Cloglog | \(n_i e^{2\eta}(1-\mu)/\mu\) | lower tail |
| Poisson | \(e^{\eta_i}\) | lower tail |
| Gamma (log) | \(1/\phi\) | only if \(\phi \to \infty\) |

### III.3 Package implementation

| Object | Code |
|---|---|
| \(\mathcal{W}_j = Z_j^\top W_j Z_j\) | `crossprod(Z_j, Z_j * w[rows])` |
| \(P_{11},P_{12},P_{21},P_{22}\) | `.two_block_S_P11(inp)` |
| \(S = P_{12}P_{22}^{-1}P_{21}\), \(\lambda^\star\) | `.two_block_gen_eigen(S, P11)` → `two_block_rate()` |
| Weights at mode | `two_block_mode_weights()` |

---

## Part IV — Building \(A\) with \(\lambda^\star < \tau\) (example \(\tau = 0.95\))

**Safe / unsafe partition (\(A / A^c\)).** Following
`ELLIPSOID_TV_BOUND.md` and `LOGIT_SINGLE_GROUP_SAFE_REGION.md`, fix a
certification threshold \(\tau \in (0,1)\). This is a **choice** (any target
contraction factor strictly below \(1\); \(0.95\) is illustrative). Define the
**safe** and **unsafe** subsets of \((\gamma,\beta)\)-space by the Gibbs rate:
\[
A
\;:=\;
\bigl\{(\gamma,\beta) \text{ at a working point} :
\lambda_{\max}\!\bigl(A(\gamma,\beta)\bigr) < \tau \bigr\},
\qquad
A^c \;:=\; \{\text{complement of } A\},
\]
where \(A(\gamma,\beta)\) on the right is the rate **matrix** from §I.4
(context distinguishes matrix from set). Equivalently,
\(\lambda^\star(\beta) \ge \tau\) on \(A^c\). In the univariate scalar case
(§IV.2), this reduces to a **group average** above \(\tau\) — that is the
\(A^c\) used throughout the TV-bound notes.

This is where Remark-8 certifies contraction at rate strictly below \(1\) on
\(A\) (Corollary 1 of `Chapter-C03`).

### IV.1 General construction (multivariate)

1. **Fix** prior blocks \((P_b, V_k^{-1}, H_j)\) and designs \((Z_j)\).
2. **At each candidate** \((\gamma,\beta)\), compute \(w_i(\eta_i)\), hence
   \(\mathcal{W}_j = Z_j^\top W_j Z_j\), \(B_j = \mathcal{W}_j + P_b\), and
   \(P_{11},P_{12},P_{21},P_{22}\).
3. **Form** \(A = P_{11}^{-1/2} P_{12} P_{22}^{-1} P_{21} P_{11}^{-1/2}\).
4. **Define** \(A = \{(\gamma,\beta) : \lambda_{\max}(A_{\mathrm{mat}}) < \tau\}\)
   and \(A^c\) as its complement, writing \(A_{\mathrm{mat}}\) for the rate
   matrix when both objects appear in one display.

The boundary is implicit; check the mode, pilot draws, or a comparison bound.
The package evaluates \(\lambda^\star\) at `b_mode` via
`two_block_rate(weights = two_block_mode_weights(...))`.

### IV.2 Scalar specialization (\(q=1\), \(p_{\mathrm{re}}=1\), \(x_j \equiv 1\))

Then \(\lambda^\star = A\) and (§2.4), with \(x_j \equiv 1\),
\[
A
\;=\;
\frac{\lambda_b}{J\bigl(\lambda_b + \lambda_\gamma/J\bigr)}
\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}
\;=\;
\frac{\lambda_b^2}{J\lambda_b + \lambda_\gamma}
\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}.
\]
Dividing \(A < \tau\) by \(\lambda_b/J\) gives the **per-group average**
form (clearer when \(J\) is large):
\[
\boxed{
\frac{1}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j(\beta)}
\;<\;
c_\tau,
\qquad
c_\tau := \tau\left(1 + \frac{\lambda_\gamma}{J\lambda_b}\right),
}
\]
and \(A = \{\beta : \text{the inequality above holds}\}\), \(A^c = \{\beta :
\text{it fails}\}\).
For \(\tau = 0.95\), \(c_{0.95} = 0.95\bigl(1 + \lambda_\gamma/(J\lambda_b)\bigr)\).
Multiplying through by \(J\) recovers the unnormalized threshold
\(\sum_j 1/(\lambda_b + W_j) < \tau\bigl(J + \lambda_\gamma/\lambda_b\bigr)\).

**Shrinkage weights.** Each summand in the flat-\(\gamma\) limit is naturally
written as a **group shrinkage weight**
\[
\omega_j(\beta)
\;:=\;
\frac{\lambda_b}{\lambda_b + W_j(\beta)}
\;\in\;(0,1],
\]
the Gibbs coupling fraction for group \(j\): how strongly that group's random
effect is still **pulled toward the population mean** \(\gamma\) (through the
RE prior \(\lambda_b\)) relative to the data precision \(W_j\). When
\(W_j = 0\) (no data information), \(\omega_j = 1\) — full prior-dominated
shrinkage; as \(W_j \to \infty\), \(\omega_j \to 0\) — the group is
data-dominated and barely couples to \(\gamma\). Define the **average
shrinkage**
\[
\bar\omega(\beta)
\;:=\;
\frac{1}{J}\sum_{j=1}^{J} \omega_j(\beta).
\]
When \(\lambda_\gamma = 0\), \(\lambda^\star = \lambda_b\,\bar\omega\): the
rate is RE prior precision times mean shrinkage. The safe/unsafe split is:
**average shrinkage toward \(\gamma\) must not be so large that
\(\lambda_b\,\bar\omega\) exceeds \(\tau\).** On \(A\),
\(\bar\omega < \tau/\lambda_b\); on \(A^c\), \(\bar\omega \ge \tau/\lambda_b\)
— too much cross-group pulling on average, hence \(\lambda^\star\) too close
to \(1\) and slow mixing. With \(\lambda_\gamma > 0\),
\(\bar\alpha = \bar\omega/\lambda_b\) and
\(\lambda^\star = \dfrac{\lambda_b}{\lambda_b + \lambda_\gamma/J}\,\bar\alpha\);
the \(\gamma\)-prior factor is an extra global scale on top of the same
\(\omega_j\).

**Interpretation.** Each \(\omega_j \le 1\), with \(\omega_j = 1\) when
\(W_j = 0\). When \(\lambda_\gamma = 0\), uninformative groups contribute
\(\omega_j = 1\) and \(\lambda^\star = \lambda_b\,\bar\omega\); if every
group is uninformative (\(\bar\omega = 1\)), \(\lambda^\star = \lambda_b\).
For \(\tau = 0.95\) with flat \(\gamma\), one needs \(\lambda_b < 0.95\) for
\(\bar\omega = 1\) to remain safe; with larger \(\lambda_b\), at least one
(typically several) \(W_j\) must be **meaningfully \(\gtrsim \lambda_b\)** so
that some \(\omega_j\) fall materially below \(1\) and pull \(\bar\omega\)
below \(\tau/\lambda_b\). A weak \(\gamma\) prior (\(\lambda_\gamma
\downarrow\)) shrinks \(c_\tau \downarrow \tau\), **tightening** the cap on
\(\bar\omega\).

**Limit \(\lambda_\gamma \to 0\)** (flat / improper prior on \(\gamma\)). Then
\(c_\tau \to \tau\), \(\lambda^\star = \lambda_b\,\bar\omega\), and
\[
A = \{\beta : \bar\omega(\beta) < \tau/\lambda_b\},
\qquad
A^c = \{\beta : \bar\omega(\beta) \ge \tau/\lambda_b\}.
\]
For \(\tau = 0.95\): **average shrinkage must stay below \(\tau/\lambda_b\).**
Each \(\omega_j = 1\) when the data contribute nothing (\(W_j = 0\)) and
drops toward \(0\) as \(W_j \to \infty\). The worst case (\(W_j \equiv 0\))
gives \(\bar\omega = 1\); safe only if \(\lambda_b < \tau\). With \(\lambda_b\)
fixed above that, \(A\) is the set of \(\beta\) for which **enough groups
carry** \(W_j \gtrsim \lambda_b\) that \(\bar\omega\) stays below
\(\tau/\lambda_b\).

### IV.2.1 The unsafe region \(A^c\) (univariate)

Fix \(\tau \in (0,1)\). In the scalar setup above, \(\lambda^\star(\beta)\)
is a **single number**, so \(A^c\) is a **one-sided inequality on a group
average**, not a per-group eigenvalue check.

**General \(\lambda_\gamma\).** With
\[
\bar\alpha(\beta)
\;:=\;
\frac{1}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j(\beta)},
\qquad
\lambda^\star(\beta)
\;=\;
\frac{\lambda_b}{\lambda_b + \lambda_\gamma/J}\,\bar\alpha(\beta),
\]
\[
A = \{\beta : \bar\alpha(\beta) < c_\tau\}
\;=\;
\{\beta : \lambda^\star(\beta) < \tau\},
\]
\[
\boxed{
A^c
\;=\;
\{\beta : \bar\alpha(\beta) \ge c_\tau\}
\;=\;
\{\beta : \lambda^\star(\beta) \ge \tau\}.
}
\]
For any chosen \(\tau\), **\(A^c\) is exactly** the set of \(\beta\) at which
\(\bar\omega \ge \tau/\lambda_b\) (flat \(\gamma\); with prior,
\(\bar\alpha \ge c_\tau\)) — **average shrinkage toward \(\gamma\) is too
strong.** This is the object bounded in `LOGIT_SINGLE_GROUP_SAFE_REGION.md`
(single group: \(A^c = \{|b| > b^\star\}\), where separation drives
\(\omega \uparrow\); here: \(\bar\omega\) too high). Remark-8 gives no
contraction factor strictly below \(1\) on \(A^c\).

**Flat \(\gamma\) (\(\lambda_\gamma \to 0\)).** Then \(c_\tau \to \tau\),
\(\lambda^\star = \lambda_b\,\bar\omega\), and
\[
A = \{\beta : \bar\omega(\beta) < \tau/\lambda_b\},
\qquad
A^c = \{\beta : \bar\omega(\beta) \ge \tau/\lambda_b\}.
\]
Each \(\omega_j\) is that group's shrinkage weight; \(A^c\) is where their
**mean** is too large. Logit \(A^c\) is typically **two-sided** in \(\eta\)
(separation lowers \(W_j\), raising \(\omega_j\)).

### IV.3 Per likelihood: building \(A\) (example \(\tau = 0.95\))

Write \(W_j(\beta)\) from §II at \(\eta_{j,i} = \mathrm{offset}_{j,i} +
z_{j,i}^\top \beta_j\). Membership in \(A\) requires
\((1/J)\sum_j 1/(\lambda_b + W_j) < c_{0.95}\); \(\beta \in A^c\) when the
average meets or exceeds \(c_{0.95}\).

| Likelihood | \(W_j\) | How to build \(A\) (\(\tau=0.95\)) | Geometry of \(A^c\) |
|---|---|---|---|
| **Normal** | \(n_j/\phi\) (constant) | All \(\beta\) if \(\lambda^\star < 0.95\), else \(\varnothing\). Tune \(n_j/\phi\), \(\lambda_b\), \(\lambda_\gamma\). | No \(\eta\)-dependence |
| **Gamma (log)** | \(n_j/\phi\) (constant) | Same as Normal | No \(\eta\)-dependence |
| **Logit** | \(\sum_i n_{j,i} p_{j,i}(1-p_{j,i})\) | Inequality on \(\beta\); sufficient inner box: \(|\eta_{j,i}| \le b^\star(\kappa)\) with \(\kappa\) from \(p(1-p)\ge\kappa\), tuned so \((1/J)\sum_j 1/(\lambda_b + \kappa\sum_i n_{j,i}) < c_{0.95}\). | **Two-sided** (separation both tails) |
| **Probit** | \(\sum_i n_{j,i}\, \phi^2/(\Phi\Phi)\) | Same average inequality; sufficient: \(\eta_{j,i} \ge \eta^{\mathrm{floor}}_{j,i}\) from numerical \(w(\eta)\ge w_{\min}\). | **Mostly lower** tail |
| **Cloglog** | \(\sum_i n_{j,i}\, e^{2\eta}(1-\mu)/\mu\) | Same average for \(A\); lower floors on \(\eta\) (§IV.3 Poisson detail). | Lower tail |
| **Poisson** | \(\sum_i e^{\eta_{j,i}}\) | Same average; sufficient lower floors \(e^{\eta_{j,i}} \ge s_j > 0\) with \((1/J)\sum_j 1/(\lambda_b + n_j s_j) < c_{0.95}\) (tune \(s_j\)). | Lower tail |

**Logit detail.** With common \(\eta_j\) per group,
\(p(1-p) \ge \kappa \iff |\eta_j| \le b^\star(\kappa)\). Pick the largest
\(\kappa\) (hence smallest \(b^\star\)) such that
\[
\frac{1}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + \kappa \sum_{i\in j} n_{j,i}}
\;<\;
c_{0.95}.
\]
Then \(A \supseteq \{\beta : |\eta_{j,i}(\beta)| \le b^\star(\kappa)\ \forall j,i\}\)
(equivalently, that box is a **sufficient inner** subset of the safe region).

**Poisson detail.** Lower floors: choose \(s_j > 0\) for each group with
\[
\frac{1}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + n_j s_j} \;<\; c_{0.95}.
\]
Then \(\eta_{j,i} \ge \log s_j\) for all \(i \in j\) (when \(z_{j,i}=1\),
\(\beta_j \ge \log s_j\)) is a **sufficient** inner region for \(A\). Tune
\(s_j\) separately when \(n_j\) or offsets differ. **Special case:** equal
\(n_j\), exchangeable groups, and \(z=1\) allow a common \(s_j \equiv s\),
reducing the average to \(1/(\lambda_b + n s) < c_{0.95}\) and
\(\beta_j \ge \log s\) for all \(j\).

### IV.4 Multivariate (\(p_{\mathrm{re}} > 1\) or \(q > 1\))

Require \(\lambda_{\max}(A) < \tau\). Strategies:

1. Evaluate at `b_mode` and on pilot draws.
2. **Comparison bound** (Remark 10): upper-bound each \(\mathcal{W}_j\) on a
   box in \(\beta\)-space to get a conservative **superset** of \(A\).
3. **Sufficient:** \(\lambda_{\min}(\mathcal{W}_j) \ge \delta\lambda_b\) on
   enough groups that a Gershgorin bound on \(\lambda_{\max}(A)\) stays \(< \tau\).

### IV.5 Choosing \(\tau\): global certification (Normal) vs.\ small \(A^c\) (GLM)

The certification threshold \(\tau\) and the priors/design are chosen **together**.
The distinction is whether \(\lambda^\star\) depends on the working point
\(\beta\).

**Normal and Gamma (log) — empty \(A^c\) is achievable.** As in §1.4,
\(\lambda^\star\) is **constant** in \(\beta\). Workflow:

| Step | Action |
|---|---|
| 1 | Fix \((\lambda_b, \lambda_\gamma, \{n_j\}, \phi)\) from data and priors. |
| 2 | Compute \(\lambda^\star\) once (`two_block_rate()`). |
| 3 | Either pick \(\tau > \lambda^\star\) (often \(\tau = \lambda^\star + \varepsilon\)), **or** adjust priors/design until \(\lambda^\star < \tau\) for a pre-set \(\tau\) (e.g.\ \(0.95\)). |
| 4 | Optionally require step 3 to hold at \(\lambda_\gamma = 0\) (weak \(\gamma\)) by taking \(\tau\) above \(\lambda^\star(\lambda_\gamma = 0)\) or by enforcing enough \(n_j/\phi\) (§1.4). |

Then **every** \(\beta\) lies in \(A\) and \(A^c = \varnothing\): the TV-bound
machinery never needs to integrate over an unsafe tail in \(\beta\)-space for
rate reasons. This is the sense in which Normal is “simplest but still
interesting”: it is the only case where threshold, priors, and design can be
aligned to **eliminate** rate-\(A^c\) entirely.

**GLM (logit, probit, Poisson, cloglog) — \(A^c\) has geometry.** Here
\(W_j(\beta)\) varies with \(\eta_{j,i}(\beta)\), so \(\lambda^\star(\beta)\)
and \(A^c\) are **nonempty subsets** of \(\beta\)-space in general (§IV.3).
One **cannot** pick \(\tau\) and priors so that every \(\beta \in A\): separation
tails always enter \(A^c\) if \(\beta\) is allowed to range. What replaces
global certification:

1. **Pick \(\tau\)** (e.g.\ \(0.95\)) as the contraction factor to certify
   on \(A\), not as the achieved rate everywhere.

2. **Calibrate priors/design** so that a **sufficient inner region**
   (§IV.3 logit box, Poisson floor, etc.) lies inside \(A\), and the mode /
   typical draws fall there.

3. **Make \(A^c\) probabilistically small** under the quantity that matters for
   TV bounds: e.g.\ bound \(\pi(A^c \mid y)\) or \(\pi(A^c \mid \theta, y)\) in
   `LOGIT_SINGLE_GROUP_SAFE_REGION.md`, or require pilot-draw fraction
   \(\frac{1}{N}\sum_k \mathbf 1\{\beta^{(k)} \in A^c\}\) below a tolerance.
   Large \(b^\star\) (logit) or strong floors on \(\eta\) (Poisson) shrink
   \(A^c\); separation pushes mass into \(A^c\) but that mass can be shown small
   in \(\pi\).

4. **Weak \(\gamma\) prior** still raises \(\lambda^\star(\beta)\) for every
   \(\beta\) (§2.4), enlarging \(A^c\); there is no Normal-style “pick \(\tau\)
   above \(\lambda_b\) and ignore \(\beta\)” fix — one must either strengthen
   \(\lambda_\gamma\), increase information in \(W_j(\beta)\) on the support
   of \(\pi\), or accept a larger \(A^c\) and bound its probability.

**Summary.**

| Likelihood | Can \(A^c = \varnothing\)? | Role of \(\tau\) |
|---|---|---|
| Normal, Gamma (log) | **Yes** — tune \(\tau\), priors, \(n_j/\phi\) | Set at or above achieved \(\lambda^\star\); can robustify to weak \(\gamma\) |
| Logit, probit, Poisson, cloglog | **No** (nontrivial \(A^c\) in \(\beta\)) | Fixed certification target on \(A\); make \(\pi(A^c)\) small via design + inner safe box |

---

## References

- `vignettes/Chapter-C03.Rmd` — §2.2 (Claim 2, Remark 8), §2.3 (Claim 3, \(A\)).
- `R/two_block_ergodicity.R` — `.two_block_S_P11()`, `.two_block_gen_eigen()`.
- `inst/ELLIPSOID_TV_BOUND.md` — TV split and alternate safe-region definitions.
- `inst/P_MATRIX_MARGINAL_PRECISIONS.md` — extended \(P\)-blocks with dispersion priors.
