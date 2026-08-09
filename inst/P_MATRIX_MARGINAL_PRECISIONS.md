# Marginal \(P\)-blocks when precisions have truncated Gamma priors

A short companion to `BLOCK_GIBBS_ERGODICITY_ING.md` §14–§16. Keep only
the two-block coordinates \((\gamma,\beta)\) and integrate **both**
measurement precisions \(\Omega_j\) and population precisions \(\lambda_p\)
out under truncated-Gamma priors. Notation: `inst/notation.md`. Residuals
are **functions of the two-block coordinates**, never fixed plugs:

\[
e_j(\beta_j)\;:=\;y_j-D_j\beta_j,
\qquad
u_p(\beta_{\cdot p},\gamma_p)\;:=\;\beta_{\cdot p}-W_p\gamma_p.
\]

Moments \(E_t[\cdot]\), \(\mathrm{Var}_t[\cdot]\) are tilted-truncated
posterior moments (`BLOCK_GIBBS_ERGODICITY_ING.md` §16.6); they inherit
that dependence through the updated rates below.

---

## 1. The \(P\) matrix (both \(\Omega_j\) and \(\lambda_p\) truncated)

### 1.1 Residuals and truncated-Gamma moments

| Symbol | Definition |
|---|---|
| \(e_j(\beta_j)=y_j-D_j\beta_j\) | measurement residual, group \(j\) (length \(n_j\)) |
| \(u_p(\beta_{\cdot p},\gamma_p)=\beta_{\cdot p}-W_p\gamma_p\) | RE residual for component \(p\) (length \(J\)) |
| \(a_j^0,r_j^0\) | Gamma shape/rate for \(\Omega_j\) (before data update) |
| \(a_p^{0(\lambda)},r_p^{0(\lambda)}\) | Gamma shape/rate for \(\lambda_p\) |
| \([\omega_{L,j},\omega_{U,j}]\) | truncation window for \(\Omega_j\) (precision scale) |
| \([\lambda_{L,p},\lambda_{U,p}]\) | truncation window for \(\lambda_p\) |

Tilted truncated posterior moments \(E_t[\cdot]\), \(\mathrm{Var}_t[\cdot]\)
are as in `BLOCK_GIBBS_ERGODICITY_ING.md` §16.6 (three `pgamma` calls at
the updated rate). For \(\Omega_j\) the updated rate is

\[
t_j(\beta_j)
\;=\;
r_j^0+\tfrac12\,\|e_j(\beta_j)\|^2
\;=\;
r_j^0+\tfrac12\,\|y_j-D_j\beta_j\|^2;
\]

for \(\lambda_p\) it is

\[
t_p^{(\lambda)}(\beta_{\cdot p},\gamma_p)
\;=\;
r_p^{0(\lambda)}+\tfrac12\,\|u_p(\beta_{\cdot p},\gamma_p)\|^2
\;=\;
r_p^{0(\lambda)}+\tfrac12\,\|\beta_{\cdot p}-W_p\gamma_p\|^2.
\]

Thus \(E_t[\Omega_j]=E_t[\Omega_j;t_j(\beta_j)]\) and
\(E_t[\lambda_p]=E_t[\lambda_p;t_p^{(\lambda)}(\beta_{\cdot p},\gamma_p)]\)
(and likewise \(\mathrm{Var}_t\)).

### 1.2 Ellipsoid scores (as functions)

**Measurement (group \(j\)).** Assume \(D_j\) has full column rank. Let
\(\hat\beta_j^{\mathrm{ols}}=(D_j'D_j)^{-1}D_j'y_j\) (data-only; does not
depend on \(\beta_j\)). Define

\begin{align*}
q_j(\beta_j)
&\;:=\;
\bigl(D_j'e_j(\beta_j)\bigr)'(D_j'D_j)^{-1}\bigl(D_j'e_j(\beta_j)\bigr)
\\[4pt]
&\;=\;
\bigl\|\beta_j-\hat\beta_j^{\mathrm{ols}}\bigr\|^2_{D_j'D_j}
\;=\;
(\beta_j-\hat\beta_j^{\mathrm{ols}})'(D_j'D_j)\,(\beta_j-\hat\beta_j^{\mathrm{ols}}).
\end{align*}

With \(\mathrm{RSS}_j^{\mathrm{ols}}=\|y_j-D_j\hat\beta_j^{\mathrm{ols}}\|^2\)
fixed, the OLS split also gives
\(q_j(\beta_j)=\|e_j(\beta_j)\|^2-\mathrm{RSS}_j^{\mathrm{ols}}\).

**Population / RE (component \(p\)).** Define

\begin{align*}
q_p^{(\lambda)}(\beta_{\cdot p},\gamma_p)
&\;:=\;
\|u_p(\beta_{\cdot p},\gamma_p)\|^2
\;=\;
\|\beta_{\cdot p}-W_p\gamma_p\|^2
\\[4pt]
&\;=\;
\sum_{j=1}^{J}(\beta_{jp}-W_{pj}\gamma_p)^2.
\end{align*}

### 1.3 Inside the ellipsoid (both cases)

A rank-one Hessian of the form \(E_t[\theta]\,M-\mathrm{Var}_t[\theta]\,vv'\)
is positive (semi)definite iff the score that \(v\) induces on \(M\) does not
exceed the threshold \(E_t[\theta]/\mathrm{Var}_t[\theta]\)
(`BLOCK_GIBBS_ERGODICITY_ING.md` §16.3 / §16.6).

| Case | Inside the ellipsoid iff | Radius |
|---|---|---|
| \(\Omega_j\) | \(\displaystyle q_j(\beta_j) \;\le\; \frac{E_t[\Omega_j]}{\mathrm{Var}_t[\Omega_j]}\) | \(K_j^{(\Omega)}(\beta_j):=E_t[\Omega_j]/\mathrm{Var}_t[\Omega_j]\) (depends on \(\beta_j\) through \(e_j(\beta_j)\) and \(t_j(\beta_j)\)) |
| \(\lambda_p\) | \(\displaystyle q_p^{(\lambda)}(\beta_{\cdot p},\gamma_p) \;\le\; \frac{E_t[\lambda_p]}{\mathrm{Var}_t[\lambda_p]}\) | \(K_p^{(\lambda)}(\beta_{\cdot p},\gamma_p):=E_t[\lambda_p]/\mathrm{Var}_t[\lambda_p]\) (depends on \((\beta_{\cdot p},\gamma_p)\) through \(u_p(\beta_{\cdot p},\gamma_p)\) and \(t_p^{(\lambda)}(\beta_{\cdot p},\gamma_p)\)) |

**Limit \(\mathrm{Var}_t\to 0\)** (hold \(E_t\) fixed). Treat this as its
own limit (not as a statement about the truncation ends). With
\(E_t[\Omega_j]\) and \(E_t[\lambda_p]\) held fixed,

\[
K_j^{(\Omega)}(\beta_j)\;=\;\frac{E_t[\Omega_j]}{\mathrm{Var}_t[\Omega_j]}
\;\to\; +\infty,
\qquad
K_p^{(\lambda)}(\beta_{\cdot p},\gamma_p)\;=\;\frac{E_t[\lambda_p]}{\mathrm{Var}_t[\lambda_p]}
\;\to\; +\infty,
\]

so every finite \((\beta,\gamma)\) lies inside both ellipsoids: the PSD
constraints become vacuous and one recovers the Normal / fixed-precision
Hessians \(H_j^{(\Omega)}\to E_t[\Omega_j]\,D_j'D_j\),
\(H_u^{(p)}\to E_t[\lambda_p]\,I_J\) (rank-one corrections vanish).

**Prior-weight limit \(w\to 1\)** (hold \(E_t\) fixed; the interesting end).
Let \(w\in(0,1)\) be the Gamma prior weight on the precision
(`pwt_measurement` for \(\Omega_j\); analogously the population prior
weight for \(\lambda_p\)), entering through
\(n_{\mathrm{prior}}=w/(1-w)\,n\) into the tilted shape/rate
\((s_j,t_j)\). Truncated moments
(`BLOCK_GIBBS_ERGODICITY_ING.md` §16.6):

\[
E_t[\Omega_j]
=\frac{s_j}{t_j}\,\frac{\Delta P(s_j+1,t_j)}{\Delta P(s_j,t_j)},
\qquad
\mathrm{Var}_t[\Omega_j]
=E_t[\Omega_j^2]-\bigl(E_t[\Omega_j]\bigr)^2.
\]

Under the usual calibration the prior mean is anchored, so raising \(w\)
shrinks \(\mathrm{Var}_t\) at essentially fixed \(E_t\). Thus

\[
w\;\to\;1
\qquad\Rightarrow\qquad
n_{\mathrm{prior}}\;\to\;+\infty
\qquad\Rightarrow\qquad
\mathrm{Var}_t[\Omega_j]\;\to\;0
\qquad\Rightarrow\qquad
K_j^{(\Omega)}(\beta_j)=\frac{E_t[\Omega_j]}{\mathrm{Var}_t[\Omega_j]}
\;\to\;+\infty.
\]

The ellipsoid \(\{q_j(\beta_j)\le K_j^{(\Omega)}(\beta_j)\}\) therefore
expands without bound. In theory one can give the prior enough weight that
the posterior mass on the outside (where \(H_j^{(\Omega)}\) can fail to be
PSD) becomes arbitrarily small — a concrete route to the
\(\mathrm{Var}_t\to 0\) Normal end above. The same reading applies to
\(K_p^{(\lambda)}\) with \(E_t[\lambda_p]\) held fixed.

(The opposite end \(w\to 0\) only lowers \(K\) toward a finite floor; the
truncation window is a different \(\mathrm{Var}_t\)-knob, not a prior
weight.)

### 1.4 Hessians (PSD inside the ellipsoids)

Write \(e_j=e_j(\beta_j)\) and
\(u_p=u_p(\beta_{\cdot p},\gamma_p)\) for short. Moments and rank-one
directions both move with \((\gamma,\beta)\):

\begin{align*}
H_j^{(\Omega)}(\beta_j)
&= E_t[\Omega_j;t_j(\beta_j)]\,D_j'D_j
- \mathrm{Var}_t[\Omega_j;t_j(\beta_j)]\,
\bigl(D_j'e_j(\beta_j)\bigr)\bigl(D_j'e_j(\beta_j)\bigr)',
\\[4pt]
H_u^{(p)}(\beta_{\cdot p},\gamma_p)
&= E_t[\lambda_p;t_p^{(\lambda)}(\beta_{\cdot p},\gamma_p)]\,I_J
- \mathrm{Var}_t[\lambda_p;t_p^{(\lambda)}(\beta_{\cdot p},\gamma_p)]\,
u_p(\beta_{\cdot p},\gamma_p)\,u_p(\beta_{\cdot p},\gamma_p)'.
\end{align*}

| Hessian | PSD / PD precisely when |
|---|---|
| \(H_j^{(\Omega)}(\beta_j)\) | inside the \(\Omega_j\) ellipsoid of §1.3 (\(D_j\) full column rank \(\Rightarrow\) PD on the nose of the boundary) |
| \(H_u^{(p)}(\beta_{\cdot p},\gamma_p)\) | inside the \(\lambda_p\) ellipsoid of §1.3 |

Stack
\(H^{(\Omega)}(\beta):=\mathrm{blockdiag}_j\bigl(H_j^{(\Omega)}(\beta_j)\bigr)\)
and
\(H^{(u)}(\gamma,\beta):=\sum_p\mathrm{embed}_p\bigl(H_u^{(p)}(\beta_{\cdot p},\gamma_p)\bigr)\).

### 1.5 Joint precision blocks

All blocks below are functions of \((\gamma,\beta)\) through
\(e_j(\beta_j)\) and \(u_p(\beta_{\cdot p},\gamma_p)\).

| Block | Formula |
|---|---|
| \(P_{11}(\gamma,\beta)\) | \(\displaystyle\sum_p W_p'\,H_u^{(p)}(\beta_{\cdot p},\gamma_p)\,W_p \;+\; \mathrm{blockdiag}_p(V_p^{-1})\) |
| \(P_{12}(\gamma,\beta)\) | row-strips \(\bigl(-W_p'\,H_u^{(p)}(\beta_{\cdot p},\gamma_p)\bigr)_{p=1}^{P}\) over the \(\beta_{\cdot p}\) columns (zeros for \(p\neq p'\)) |
| \(P_{21}(\gamma,\beta)\) | \(P_{12}(\gamma,\beta)^{\top}\) |
| \(P_{22}(\gamma,\beta)\) | \(H^{(\Omega)}(\beta)+H^{(u)}(\gamma,\beta)\) |

\(\mathrm{embed}_p\bigl(H_u^{(p)}(\beta_{\cdot p},\gamma_p)\bigr)\) places
the \(J\times J\) matrix \(H_u^{(p)}(\beta_{\cdot p},\gamma_p)\) into the
\(\beta_{\cdot p}\) slots of the full \(JP\times JP\) array.

| Source | \(P_{11}\) | \(P_{12}\) | \(P_{21}\) | \(P_{22}\) |
|---|---|---|---|---|
| Truncated prior on \(\Omega_j\) \(\to H_j^{(\Omega)}\) | — | — | — | yes |
| Truncated prior on \(\lambda_p\) \(\to H_u^{(p)}\) | yes | yes | yes | yes |
| Fixed \(\gamma\) prior \(V_p^{-1}\) | yes | — | — | — |

**Normal / fixed-precision limit** (\(\mathrm{Var}_t\to 0\) of §1.3 with
\(E_t[\Omega_j]=\omega_j\), \(E_t[\lambda_p]=\lambda_p\) held fixed;
ellipsoids become all of \((\gamma,\beta)\)-space. Rank-one terms in
\(e_j(\beta_j)\) and \(u_p(\beta_{\cdot p},\gamma_p)\) drop out, so the
limiting \(P\) no longer depends on those residuals.)

| Block | Limit |
|---|---|
| \(P_{11}\) | \(\sum_p \lambda_p\,W_p'W_p + \mathrm{blockdiag}_p(V_p^{-1})\) |
| \(P_{12}\) | \(\bigl(-\lambda_p\,W_p'\bigr)_p\) |
| \(P_{21}\) | \(P_{12}^{\top}\) |
| \(P_{22}\) | \(\mathrm{blockdiag}_j(\omega_j D_j'D_j)+\mathrm{diag}(\lambda)\) on each \(\beta_j\) |

### 1.6 Rate matrix and sufficient condition

Whenever \(P_{11}(\gamma,\beta)\succ 0\) and \(P_{22}(\gamma,\beta)\succ 0\),
define

\[
A(\gamma,\beta)
\;=\;
P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2},
\qquad
\lambda^\star(\gamma,\beta):=\lambda_{\max}\bigl(A(\gamma,\beta)\bigr)
\]

(with each \(P_{ij}=P_{ij}(\gamma,\beta)\)).

**Sufficient condition (matched \(P_{12}=-W'H_u\)).** If every group is
inside its \(\Omega_j\) ellipsoid and every RE component is inside its
\(\lambda_p\) ellipsoid (§1.3), then
\(H^{(\Omega)}(\beta)\succ 0\) and \(H^{(u)}(\gamma,\beta)\succeq 0\).
With also \(P_{11}(\gamma,\beta)\succ 0\),

\[
\lambda_{\max}\bigl(A(\gamma,\beta)\bigr) \;<\; 1.
\]

Equivalently: \(P_{22}(\gamma,\beta)\succ 0\) alone is not enough — one
needs the measurement ellipsoids (so \(H^{(\Omega)}(\beta)\succ 0\)).
Matching through \(P_{12}/P_{21}\) does not remove the need for the
\(\lambda_p\) ellipsoids (so \(H^{(u)}(\gamma,\beta)\succeq 0\)); see
§5.2.

### 1.7 Special case: default `dGamma` priors (`Prior_Setup_GLMM`)

Focus: the default **Gamma priors on measurement and population
dispersion** (Block~1 `dGamma` / `group.ing_prior`; Block~2 Gamma factor
of `dIndependent_Normal_Gamma` / `pop.ing_prior`). Both are mean-matched
at A12 §3.3.4 marginal centers
(\(\hat\sigma^2\) pooled or per-group; \(\hat\tau^2_p\) from the
hyper-regression on \(b_{\cdot p}\)). See `R/Prior_Setup_GLMM.R`.

| Prior | Weight default | \(n_0=\dfrac{w}{1-w}\times(\text{units})\) | Shape / rate |
|---|---|---|---|
| Measurement \(\Omega\) (pooled) | `group.dispersion.pwt` \(=w_\Omega=0.01\) | \(n_0=\dfrac{w_\Omega}{1-w_\Omega}\,n\) | \(s_\Omega=\tfrac{n_0+1}{2}+\tfrac{p_{\mathrm{re}}}{2}\), \(r_\Omega=\hat\sigma^2\cdot\tfrac{n_0+p_{\mathrm{re}}-1}{2}\) |
| Measurement \(\Omega_j\) (per-group) | seed \(w_{\Omega,j}=0.01\); may rise under `group.alpha_target` \(=0.01\) | \(n_{0,j}=\dfrac{w_{\Omega,j}}{1-w_{\Omega,j}}\,n_j\) | same with \((n_{0,j},\hat\sigma^2_j,p_{\mathrm{re}})\) |
| Population \(\lambda_p=1/\tau_p^2\) | `pop.dispersion.pwt` \(=w_{\lambda,p}\) defaults to `pop.pwt` ( \(0.01\) if \(p_p<14\), else \(0.05\) ) | \(n_{0,p}=\dfrac{w_{\lambda,p}}{1-w_{\lambda,p}}\,J\) | \(s_{\lambda,p}=\tfrac{n_{0,p}+1}{2}+\tfrac{p_p}{2}\), \(r_{\lambda,p}=\hat\tau^2_p\cdot\tfrac{n_{0,p}+p_p-1}{2}\) |

#### Marginal centers (filled in)

**Per-group measurement.** For group \(j\), let \(\hat\beta_j\) be the
within-group OLS fit on \(D_j\), \(\mu_j\) the Block~1 prior mean for
\(\beta_j\), and \(\Sigma_j'= \Sigma_j+\Omega_j\) the Part-0 /
Part VI prior covariance (RE shrinkage plus hyper-mean uncertainty
through \(W_j\Sigma_\gamma W_j'\)). With
\(M_j=\bigl(\Sigma_j'+(D_j'D_j)^{-1}\bigr)^{-1}\),

\begin{align*}
S_{\mathrm{marg},j}
&=\mathrm{RSS}_{\mathrm{ols},j}
+(\hat\beta_j-\mu_j)'M_j(\hat\beta_j-\mu_j),
\\[4pt]
\hat\sigma^2_j
&=\frac{S_{\mathrm{marg},j}}{n_j-p_{\mathrm{re}}}.
\end{align*}

(`compute_gaussian_prior` / `.lmebayes_measurement_group_smarg_cal`.)

**Pooled measurement.** Aggregate the independent-group contributions:

\[
\hat\sigma^2
=\frac{\sum_{j=1}^{J} S_{\mathrm{marg},j}}{\sum_{j=1}^{J}(n_j-p_{\mathrm{re}})}
=\frac{\sum_{j=1}^{J} S_{\mathrm{marg},j}}{n-J\,p_{\mathrm{re}}}
\]

(when every group is full rank; singular within-group designs fall back
to intercept-only residual SS for that \(j\), with df \(n_j-1\)).

**Population \(\tau^2_p\).** Treat the reference group coefficients
\(b_{\cdot p}\) (length \(J\)) as the response in
\(b_{\cdot p}\sim N(W_p\gamma_p,\tau_p^2 I_J)\) with prior
\(\gamma_p\sim N(\mu_p,\Sigma_p)\). Let \(\hat\gamma_p\) be OLS of
\(b_{\cdot p}\) on \(W_p\) and

\begin{align*}
S_{\mathrm{marg},p}^{(\tau)}
&=\|b_{\cdot p}-W_p\hat\gamma_p\|^2
+(\hat\gamma_p-\mu_p)'
\bigl(\Sigma_p+(W_p'W_p)^{-1}\bigr)^{-1}
(\hat\gamma_p-\mu_p),
\\[4pt]
\hat\tau^2_p
&=\frac{S_{\mathrm{marg},p}^{(\tau)}}{J-p_p}.
\end{align*}

(`compute_gaussian_prior` on the hyper-regression /
`.lmebayes_compute_ing_prior_cal_tau2_hyper`.)

Truncation: `group.max_disp_perc` / `pop.max_disp_perc` default `0.99`.
By construction \(r/(s-1)\) equals the corresponding hat, so

\begin{align*}
\omega_\star
&=\frac1{\hat\sigma^2}
=\frac{n-J\,p_{\mathrm{re}}}{\sum_j S_{\mathrm{marg},j}}
\quad\text{(pooled)},
\\[4pt]
\omega_{\star,j}
&=\frac1{\hat\sigma^2_j}
=\frac{n_j-p_{\mathrm{re}}}{S_{\mathrm{marg},j}}
\quad\text{(per-group)},
\\[4pt]
\lambda_{\star,p}
&=\frac1{\hat\tau^2_p}
=\frac{J-p_p}{S_{\mathrm{marg},p}^{(\tau)}}.
\end{align*}

#### Default \(P\) (Normal / small-\(\mathrm{Var}_t\) reading)

With the filled-in `dGamma` plugs above, §1.5 specializes to
(pooled measurement; per-group replaces
\(\dfrac{n-J\,p_{\mathrm{re}}}{\sum_\ell S_{\mathrm{marg},\ell}}\) by
\(\dfrac{n_j-p_{\mathrm{re}}}{S_{\mathrm{marg},j}}\) in each \(j\)-block of
\(P_{22}\))

\begin{align*}
P_{11}
&=\sum_p \frac{J-p_p}{S_{\mathrm{marg},p}^{(\tau)}}\,W_p'W_p
+\mathrm{blockdiag}_p(V_p^{-1}),
\\[6pt]
P_{12}
&=\Biggl(-\frac{J-p_p}{S_{\mathrm{marg},p}^{(\tau)}}\,W_p'\Biggr)_p,
\qquad
P_{21}=P_{12}^{\top},
\\[6pt]
P_{22}
&=\mathrm{blockdiag}_j\Biggl(
\frac{n-J\,p_{\mathrm{re}}}{\sum_\ell S_{\mathrm{marg},\ell}}\,D_j'D_j
\Biggr)
+\mathrm{diag}\Biggl(
\biggl(\frac{J-p_p}{S_{\mathrm{marg},p}^{(\tau)}}\biggr)_p
\Biggr)
\quad\text{on each }\beta_j.
\end{align*}

Here \(V_p^{-1}\) is the Block~2 Normal prior precision on \(\gamma_p\)
(separate from the `dGamma` factor; under default scalar `pop.pwt`
\(=w_\gamma\), \(V_p^{-1}=\dfrac{w_\gamma}{1-w_\gamma}\,V_{\mathrm{fe},p}^{-1}\)).
The `dGamma` defaults enter \(P\) only through these \(S_{\mathrm{marg}}\)
ratios.

#### Default rate matrix \(A\)

\[
A
=P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2},
\qquad
\lambda^\star=\lambda_{\max}(A).
\]

Writing
\(\Lambda_\star=\mathrm{diag}\bigl((J-p_p)/S_{\mathrm{marg},p}^{(\tau)}\bigr)_p\)
in the \(\beta\)-ordering of the \(P_{12}\) strips and
\(H_\Omega^\star=\mathrm{blockdiag}_j\bigl(
\tfrac{n-J\,p_{\mathrm{re}}}{\sum_\ell S_{\mathrm{marg},\ell}}\,D_j'D_j\bigr)\)
(or per-group
\(\tfrac{n_j-p_{\mathrm{re}}}{S_{\mathrm{marg},j}}\)),

\[
P_{22}=H_\Omega^\star+\Lambda_\star,
\qquad
P_{12}=-W'\Lambda_\star,
\qquad
A
=P_{11}^{-1/2}\,
W'\Lambda_\star\,(H_\Omega^\star+\Lambda_\star)^{-1}\,\Lambda_\star W\,
P_{11}^{-1/2}.
\]

Whenever \(H_\Omega^\star\succ 0\) and \(P_{11}\succ 0\), §1.6 gives
\(\lambda^\star<1\). With the default weights
\(w_\Omega=w_{\lambda,p}=0.01\) (low-\(d\)) one sits at the **weak** end of
the §1.3 prior-weight continuum (\(n_0=\dfrac{w}{1-w}\,N\) small; \(K\)
finite). Truncation and (per-group) `alpha_target` keep mass inside the
ellipsoids without taking \(w\to 1\). Full truncated Hessians of §1.4
still apply when \(\mathrm{Var}_t>0\).

---

## 2. One Normal–Gamma factor, one pattern

Whenever a scalar precision \(\theta>0\) enters only through a Normal factor

\[
\tfrac12\,\theta\,\|r(x)\|^2 - \tfrac12\log\theta
\]

and \(\theta\) has a (possibly truncated) Gamma prior, integrating \(\theta\)
out yields a Student-\(t\)-type marginal in the residual map \(r(x)\) with
Hessian (in \(x\))

\[
H_r(x)
\;=\;
E_t[\theta;t(x)]\,I
\;-\;
\mathrm{Var}_t[\theta;t(x)]\,r(x)\,r(x)'.
\]

As \(\mathrm{Var}_t[\theta]\to 0\) with \(E_t[\theta]\to\theta_0\),
\(H_r(x)\to\theta_0 I\) (Normal / fixed-precision limit; residual
dependence drops out).

| Precision | Role | Residual map \(r(x)\) | Couples to |
|---|---|---|---|
| \(\Omega_j=1/\sigma_j^2\) | measurement | \(e_j(\beta_j)=y_j-D_j\beta_j\) | \(\beta_j\) only \(\to P_{22}\) |
| \(\lambda_p=1/\tau_p^2\) | population / RE | \(u_p(\beta_{\cdot p},\gamma_p)=\beta_{\cdot p}-W_p\gamma_p\) | \(\gamma_p\) and \(\beta_{\cdot p}\) \(\to P_{11},P_{12},P_{21},P_{22}\) |

---

## 3. Why \(\Omega_j\) hits only \(P_{22}\)

Chain rule through \(e_j(\beta_j)=y_j-D_j\beta_j\) gives
\(H_j^{(\Omega)}(\beta_j)\) in \(\beta_j\) alone. The likelihood never
involves \(\gamma\) in the centered parameterization, so \(P_{11}\) and
\(P_{12}/P_{21}\) are untouched by \(\Omega_j\).

---

## 4. Why \(\lambda_p\) hits \(P_{11}\), \(P_{12}/P_{21}\), and \(P_{22}\)

The RE prior links \(\gamma_p\) to \(\beta_{\cdot p}\) through
\(u_p(\beta_{\cdot p},\gamma_p)=\beta_{\cdot p}-W_p\gamma_p\).
Marginalizing \(\lambda_p\) leaves a \(t\)-prior on that link, so one
\(H_u^{(p)}(\beta_{\cdot p},\gamma_p)\) feeds the matched triple

\[
P_{11}\gets W_p'H_u^{(p)}(\beta_{\cdot p},\gamma_p)\,W_p,
\quad
P_{12}\gets -W_p'H_u^{(p)}(\beta_{\cdot p},\gamma_p),
\quad
P_{22}\gets H_u^{(p)}(\beta_{\cdot p},\gamma_p)
\]

(plus likelihood \(H_j^{(\Omega)}(\beta_j)\) in \(P_{22}\)).
Cross-component blocks (\(p\neq p'\)) stay zero under diagonal
\(\Psi^{-1}\) and independent ING priors.

---

## 5. Spectral roles: two different settings (do not conflate)

### 5.1 Setting A — Gaussian / fixed-precision (or joint) sub-chain

No marginal \(E_t-\mathrm{Var}_t\) Hessians: usual `.two_block_S_P11()` \(P\),
and Remark 8 always returns \(\lambda^\star\in[0,1)\).

- Raising \(\lambda_p\) pushes \(\lambda^\star\) **toward** 1, never above
  (\(H_u=\lambda_p I\succ 0\); Gaussian likelihood precision
  \(\omega D'D\succeq 0\) keeps \(P_{22}\) PD).
- Truncation on \(\lambda_p\) / \(\tau_p^2\) supplies a **gap** away from 1
  (§10 of `BLOCK_GIBBS_ERGODICITY_ING.md`). Without truncation there is no
  uniform \(\varepsilon>0\) with \(\lambda^\star\le 1-\varepsilon\).

### 5.2 Why matching helps — and why \(H_u\succeq 0\) is still needed

Do **not** write “\(L=H^{(\Omega)}\)” and then assume \(L\succ 0\).
The measurement block \(H^{(\Omega)}\) is itself a marginal Hessian and may
fail to be PD; that is a separate gate. Keep the two contributions
distinct:

\[
P_{22}=H^{(\Omega)}+H_u,
\qquad
P_{12}=-W'H_u,
\qquad
P_{11}=W'H_u W,
\qquad
S=W'H_u\,P_{22}^{-1}\,H_u\,W
\]

(one RE component; \(H_u\equiv H_u^{(p)}\) on the \(\beta_{\cdot p}\) slots).

**Claim (matched).** If \(H_u\succeq 0\), \(H^{(\Omega)}\succ 0\), and
\(P_{11}\succ 0\) (hence \(P_{22}=H^{(\Omega)}+H_u\succ 0\)), then
\(\lambda^\star<1\).

**Proof.** Let \(B:=P_{22}=H^{(\Omega)}+H_u\). For \(z=Wv\),

\begin{align*}
v'(P_{11}-S)v
&= z'H_u z - z'H_u B^{-1}H_u z
= z'H_u B^{-1}H^{(\Omega)} z.
\end{align*}

Under \(H_u\succeq 0\) and \(H^{(\Omega)}\succ 0\), one has
\(B\succeq H^{(\Omega)}\succ 0\), so
\((H^{(\Omega)})^{1/2}B^{-1}(H^{(\Omega)})^{1/2}\prec I\), and

\[
z'H_u B^{-1}H_u z - z'H_u z
= -\,z'(H^{(\Omega)})^{1/2}\Bigl(I-(H^{(\Omega)})^{1/2}B^{-1}(H^{(\Omega)})^{1/2}\Bigr)(H^{(\Omega)})^{1/2}z
\;<\; 0
\]

whenever \(H^{(\Omega)}z\neq 0\). Thus \(S\prec P_{11}\) and
\(\lambda_{\max}(A)<1\). Matching (\(H_u\) in \(P_{12}/P_{21}\) as well
as \(P_{11}/P_{22}\)) produces the identity; the hypothesis
\(H^{(\Omega)}\succ 0\) is stated on its own, not hidden inside a symbol
\(L\) that looks like a Gaussian likelihood precision. \(\square\)

**Matching alone is not enough if \(H_u\) is indefinite.** Numerically,
there exist \(H_u\not\succeq 0\), \(H^{(\Omega)}\succ 0\),
\(P_{22}\succ 0\), and \(P_{11}\succ 0\) with \(\lambda^\star\gg 1\).
So keep \(H_u\succeq 0\) in the sufficient condition.

**Unmatched \(\Omega_j\) blocks** (current §16 diagnostic; fixed RE
precision \(\Lambda\succ 0\), not \(H_u\)):

\[
P_{22}=H^{(\Omega)}+\Lambda,
\qquad
P_{12}=-W'\Lambda,
\qquad
P_{11}=W'\Lambda W,
\qquad
S=W'\Lambda\,P_{22}^{-1}\,\Lambda\,W.
\]

**Claim (unmatched).** \(P_{22}\succ 0\) is **not** enough for
\(\lambda^\star\le 1\). One needs \(H^{(\Omega)}\succeq 0\)
(equivalently \(P_{22}\succeq\Lambda\)).

**Scalar counterexample.** \(W=1\), \(\Lambda=1\), \(H^{(\Omega)}=-1/2\):
\(P_{22}=1/2\succ 0\), \(P_{11}=1\), but \(S=2>1\), so
\(\lambda^\star=2>1\).

| Gate | Matched \(\lambda_p\) (\(P_{12}=-W'H_u\)) | Unmatched \(\Omega_j\) (\(P_{12}=-W'\Lambda\)) |
|---|---|---|
| \(P_{22}\succ 0\) | follows from \(H_u\succeq 0\) and \(H^{(\Omega)}\succ 0\) | needed; **not** sufficient alone |
| \(P_{11}\succ 0\) | needed to form \(A\) | usually from \(\Lambda\succ 0\) |
| Extra | need \(H_u\succeq 0\) **and** \(H^{(\Omega)}\succ 0\) | need \(H^{(\Omega)}\succeq 0\) |

### 5.3 Summary

| Ingredient | Setting A (Gaussian \(P\)) | Setting B (marginal Hessian diagnostic) |
|---|---|---|
| \(\Omega_j\) | corner moves \(\lambda^\star\) in \([0,1)\) | unmatched: \(P_{22}\succ 0\) insufficient; need \(H_j^{(\Omega)}\succeq 0\) |
| \(\lambda_p\) | moves closeness to 1; never \(>1\); truncation \(\Rightarrow\) gap | matched: need \(H_u\succeq 0\), \(H^{(\Omega)}\succ 0\), \(P_{11}\succ 0\) |
| \(\gamma\) prior \(V_p^{-1}\) | moves closeness to 1 | strengthens \(P_{11}\) |

---

## 6. Relation to other notes

| Document | Role |
|---|---|
| `BLOCK_GIBBS_ERGODICITY_ING.md` §14–15 | Joint Hessian **keeping** \(\Omega,\lambda\) as \(x_1\) coordinates |
| `BLOCK_GIBBS_ERGODICITY_ING.md` §16 | Marginal \(\Omega\) only (feeds \(P_{22}\); implemented in `R/two_block_ergodicity_ing_marginal.R`) |
| `BLOCK_GIBBS_ERGODICITY_ING.md` §10 | Truncation corners for \(\lambda_p,\Omega\) give a certified gap from 1 |
| **This note** | Marginal \(\Omega\) **and** \(\lambda\): §1 table; default \(P\)/\(A\) in §1.7; spectral roles in §5 |
| `DERIVATION_sigma2_with_block2_v2.md` | Sampler partition that *samples* \(\Omega\) with Block~2 |

Implementation: measurement-marginal \(P_{22}\) path exists
(`.two_block_lambda_star_marginal_over_draws`). Population-marginal
\(H_u^{(p)}\) feeding \(P_{11}/P_{12}/P_{22}\) not yet wired.
