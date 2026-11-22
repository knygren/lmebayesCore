# Explicit Convergence Bounds for Two-Block Gibbs Samplers in Hierarchical Models with Log-Concave Likelihoods

*Draft.*

---

## Abstract

We consider the two-block Gibbs sampler for hierarchical models with a Gaussian population layer and log-concave group likelihoods, as in hierarchical generalized linear models. If each group's log-likelihood has curvature bounded below by a matrix $C_{i,\min}$, and possibly above by $C_{i,\max}$, the sampler's convergence can be compared with that of two Gaussian hierarchical models: those whose group likelihoods are Gaussian with precisions $C_{i,\min}$ and $C_{i,\max}$. We prove three results.

1. A coupling contracts at the exact rate $\rho^{\rm G}(C_{\min})$ of the first Gaussian model.
2. The total-variation distance to stationarity is bounded explicitly at that rate, from any start, with a simplified form for starts between the posterior mode and mean.
3. The $L^2$ geometric rate lies between $\rho^{\rm G}(C_{\max})$ and $\rho^{\rm G}(C_{\min})$.

The rate is below one for any proper prior on the population parameter. For common GLMs the global curvature floor is zero, and the rate then degrades toward one as the prior weakens. A ridge modification, mixing each log-likelihood with its quadratic approximation at weight $\omega$, restores a floor that does not depend on the prior. We bound the total-variation distance between the true and modified posteriors (linear in $\omega$ to first order), which gives an end-to-end bound for the distance between the modified sampler's output and the true posterior.

The comparison is between **rates**, not between distances. The total-variation bound is weaker than the exact Gaussian analysis of Nygren (2020), and a comparison of actual distances with the Gaussian sampler remains open.

---

## 1. Introduction

Explicit convergence bounds for Gibbs samplers usually rest on drift and minorization conditions (Rosenthal 1995), whose constants are hard to make sharp. For Gaussian hierarchical models, the two-block sampler is a Gaussian autoregression in the population parameter, and its distance from stationarity can be computed exactly (Nygren 2020). This paper extends what can be said, explicitly and with constants, to log-concave group likelihoods.

The key structural facts are:

1. **Each group-effect conditional is an exponential tilt** of a strongly log-concave measure, with the tilt parameter linear in the population parameter $\gamma$.
2. **The $\gamma$-update is exactly Gaussian.**
3. **Strong log-concavity gives contraction.** Through quantile and transport couplings, it gives a contraction of the $\gamma$-chain at the rate of a Gaussian model.
4. **The Gaussian final step turns contraction into total variation.** This is a one-shot coupling argument (Roberts and Rosenthal 2002; Madras and Sezer 2010), and no small set or minorization is needed.

**How the $k$-step kernel is avoided.** Without conjugacy, the law of $\gamma^{(k)}$ has no closed form, so the exact approach of the Gaussian case (compute $\mathcal L(\gamma^{(k)})$, then its TV distance) is unavailable. The argument never computes it.

1. **Two copies, one inequality.** Run the chain from $\gamma^{(0)}$ alongside a copy started at stationarity. Because $\pi_\gamma$ is invariant, the copy stays at stationarity, so its law at every step is $\pi_\gamma$, even though $\pi_\gamma$ is not known in closed form.
2. **Only a one-step sensitivity bound is needed.** For block 1 the argument uses no formula for the conditional $f_i(\cdot\mid\gamma)$, only how much it can move when $\gamma$ moves. Because it is an exponential tilt, in $\gamma$, of a strongly log-concave base, its quantiles (Lemma 3.4) or its $W_2$-optimal transport (Lemma 3.6) move by at most the Gaussian factor $(C_{i,\min}+\Psi_i^{-1})^{-1}$ per unit tilt. Combined with the linear block-2 mean, this gives a contraction inequality, $E\|\Delta\gamma'\|^2\le\rho^2\|\Delta\gamma\|^2$, that holds **for every pair of states**.
3. **Inequalities compose; kernels do not.** A bound valid for every pair of states can be applied at each step regardless of where the chains are, so after $k$ steps the distance is at most $\rho^k$ times the initial one. The $k$-step kernel is never written down.
4. **The only exact formula used is the last Gaussian step.** Given the final $\beta$'s, both copies' $\gamma^{(k)}$ are normal with common covariance $P_0^{-1}$. That converts the coupling distance into total variation (Lemma 3.9).
5. **The stationary law enters only through moments.** Strong log-concavity bounds the expected distance from the start to a stationary draw (Lemmas 3.7–3.8), without knowing $\pi_\gamma$.

The price of not computing the kernel is precision: the result is a bound at the Gaussian rate, not the exact distance. In particular it loses the $a^{2k}$ decay from the center that the exact Gaussian calculation gives.

Section 2 sets up the model. Section 3 collects the lemmas. Section 4 states and proves the main theorem, and §5 its corollaries. Section 6 treats the ridge modification and its distance from the true posterior. Section 7 specializes to the logit and Poisson families, §8 gives a worked example, and §9 discusses limitations and open problems.

---

## 2. Model, sampler and assumptions

### 2.1 Model

For $\gamma\in\mathbb R^p$ and groups $i=1,\dots,J$ with $\beta_i\in\mathbb R^q$:
$$\gamma\sim\mathcal N_p(\mu_\gamma,\Lambda^{-1}),\qquad\beta_i\mid\gamma\sim\mathcal N_q(H_i\gamma,\Psi_i),\qquad p(y_i\mid\beta_i)=e^{\ell_i(\beta_i)},$$
with $H_i\in\mathbb R^{q\times p}$, $\Psi_i\succ0$ and $\Lambda\succeq0$. Write
$$A_i=\Psi_i^{-1}H_i,\qquad P_0=\Lambda+\sum_iH_i^\top\Psi_i^{-1}H_i,$$
and assume $P_0\succ0$.

**What $P_0$ is.** $P_0$ is the precision of $\gamma$ **given all the group effects**: the precision of the block-2 conditional. It is neither the prior precision $\Lambda$ nor the posterior precision of $\gamma$, and it does not depend on the data.

### 2.2 Sampler

One iteration from $\gamma$:

1. **(Block 1)** draw $\beta_i\sim f_i(\cdot\mid\gamma)\propto e^{\ell_i(\beta)}\phi_q(\beta;H_i\gamma,\Psi_i)$, independently over $i$;
2. **(Block 2)** draw $\gamma'\sim\mathcal N_p\bigl(P_0^{-1}(\Lambda\mu_\gamma+\sum_iA_i^\top\beta_i),\ P_0^{-1}\bigr)$.

$K$ denotes the kernel of the $\gamma$-chain and $\pi_\gamma$ the marginal posterior of $\gamma$. $K$ is reversible with respect to $\pi_\gamma$, with non-negative spectrum (Liu, Wong and Kong 1994). $\mathcal L(\gamma^{(k)})$ is the law after $k$ iterations from a fixed $\gamma^{(0)}$.

### 2.3 Assumption

**Assumption LC (curvature bounds).** Each $\ell_i$ is concave on $\mathbb R^q$ with locally Lipschitz gradient, and
$$C_{i,\min}\preceq-\nabla^2\ell_i(\beta)\preceq C_{i,\max}\qquad\text{for almost every }\beta,$$
for symmetric $0\preceq C_{i,\min}\preceq C_{i,\max}$. The upper bound may be infinite, and is used only for lower bounds on the rate.

That is, each group's likelihood carries at least the information of a Gaussian observation with precision $C_{i,\min}$, and at most that of one with precision $C_{i,\max}$.

### 2.4 Rate quantities

**Definition 2.1.** For a collection $C=(C_i)$ of positive semidefinite matrices, let
$$\tilde A(C)=P_0^{-1/2}\Bigl[\sum_iA_i^\top\bigl(C_i+\Psi_i^{-1}\bigr)^{-1}A_i\Bigr]P_0^{-1/2},\qquad\rho^{\rm G}(C)=\lambda_{\max}\bigl(\tilde A(C)\bigr).$$
$\rho^{\rm G}(C)$ is the exact $L^2$ convergence rate of the two-block sampler in the **Gaussian** hierarchical model with group likelihoods $\mathcal N(\beta_i;\cdot,C_i^{-1})$, that is, group-likelihood precision $C_i$. Write
$$\rho_{\rm hi}=\rho^{\rm G}(C_{\min}),\qquad\rho_{\rm lo}=\lambda_{\min}\bigl(\tilde A(C_{\max})\bigr),\qquad K_i^{\rm lo}=C_{i,\min}+\Psi_i^{-1},\qquad K_i^{\rm hi}=C_{i,\max}+\Psi_i^{-1}.$$
For scalar $\gamma$ and $\beta_i$,
$$\rho_{\rm hi}=\frac1{P_0}\sum_i\frac{H_i^2}{\Psi_i(1+c_{i,\min}\Psi_i)}.$$

**Definition 2.2 (admissible reference set).** With $\bar\gamma$ and $\gamma^\star$ the mean and mode of $\pi_\gamma$,
$$\mathcal C=\Bigl\{c\in\mathbb R^p:\ E_{\pi_\gamma}\|\tilde\gamma-c\|^2_{P_0}\le\frac p{1-\rho_{\rm hi}}\Bigr\}=\Bigl\{c:\ \|c-\bar\gamma\|^2_{P_0}\le\frac p{1-\rho_{\rm hi}}-\operatorname{tr}\bigl(P_0\operatorname{Cov}_{\pi_\gamma}\bigr)\Bigr\}.$$
It is a closed ellipsoid in the $P_0$-metric, centered at the posterior mean, and it contains the segment $[\gamma^\star,\bar\gamma]$ (Lemma 3.8).

**Notation.** $\|v\|_M=(v^\top Mv)^{1/2}$; $\preceq$ is the Loewner order; $\varphi$ and $\Phi$ are the standard normal density and distribution function.

---

## 3. Lemmas

**Lemma 3.1 (tilt representation).** With $V_i(\beta)=-\ell_i(\beta)+\tfrac12\beta^\top\Psi_i^{-1}\beta$,
$$f_i(\beta\mid\gamma)\propto\exp\{\theta_i^\top\beta-V_i(\beta)\},\qquad\theta_i=A_i\gamma,$$
and $K_i^{\rm lo}\preceq\nabla^2V_i\preceq K_i^{\rm hi}$ almost everywhere. $\square$

**Lemma 3.2 (covariance sandwich).** $(K_i^{\rm hi})^{-1}\preceq\operatorname{Cov}(\beta_i\mid\gamma)\preceq(K_i^{\rm lo})^{-1}$.

*Proof.* *Upper bound:* Brascamp–Lieb. *Lower bound:* the Cramér–Rao inequality for the location family, $\operatorname{Cov}\succeq(E\nabla^2V_i)^{-1}$, which is valid because $\nabla V_i$ is locally Lipschitz. $\square$

**Lemma 3.3 (drift and curvature of $\pi_\gamma$).** With $d(\gamma)=E[\gamma'\mid\gamma]$ and $R(\gamma)=\sum_iA_i^\top\operatorname{Cov}(\beta_i\mid\gamma)A_i$,
$$\nabla d=P_0^{-1}R(\gamma),\qquad-\nabla^2\log\pi_\gamma=P_0-R(\gamma),\qquad R_{\rm lo}\preceq R(\gamma)\preceq R_{\rm hi},$$
where $R_{\rm lo/hi}=\sum_iA_i^\top(K_i^{\rm hi/lo})^{-1}A_i$. Consequently
$$-\nabla^2\log\pi_\gamma\succeq K_\gamma:=P_0-R_{\rm hi}\succeq(1-\rho_{\rm hi})P_0 .$$

*Proof.* For tilts, the derivative of the mean in $\theta$ is the covariance. Differentiating $\log\int e^{\ell_i}\phi$ twice in $\gamma$ gives the second identity (Tweedie). Then apply Lemma 3.2. $\square$

**Lemma 3.4 (quantile Lipschitz bound).** If $g_s\propto e^{sw-U(w)}$ on $\mathbb R$ with $U''\ge\kappa>0$ in the sense of distributions, the $u$-quantile $w_u(s)$ satisfies $0\le\partial_sw_u\le1/\kappa$.

*Proof.* With $G(w)=\int_{-\infty}^w(\mu_s-t)g_s(t)\,dt$, one has $\partial_sw_u=G(w_u)/g_s(w_u)\ge0$. Let $H=g_s/\kappa-G$. Then $H'=g_s\bigl(w-U'(w)/\kappa+s/\kappa-\mu_s\bigr)$, whose bracket is non-increasing. So $H'$ changes sign at most once, from $+$ to $-$. Since $H(\pm\infty)=0$, $H\ge0$. $\square$

**Lemma 3.5 (one-dimensional projections).** If $X$ has density $\propto e^{\theta^\top x-V(x)}$ with $\nabla^2V\succeq K\succ0$ and $a\in\mathbb R^q$, then $a^\top X$ has a log-concave density with curvature at least $(a^\top K^{-1}a)^{-1}$. If $\theta=sa$, the law of $a^\top X$ is a one-parameter tilt in $s$ of a fixed base.

*Proof.* Strong log-concavity is preserved under linear marginalization, with the projected constant (Saumard and Wellner 2014, §3, via Prékopa–Leindler). $\square$

**Lemma 3.6 (tilts in $W_2$).** In the setting of Lemma 3.5, $\mathrm{KL}(\mu_\theta\|\mu_{\theta'})\le\tfrac12\|\theta-\theta'\|^2_{K^{-1}}$, and
$$W_{2,K}(\mu_\theta,\mu_{\theta'})\le\|\theta-\theta'\|_{K^{-1}}.$$

*Proof.* The KL divergence is the Bregman divergence of the log-partition function, whose Hessian is $\operatorname{Cov}\preceq K^{-1}$. The $W_2$ bound is Talagrand's inequality in the coordinates $K^{1/2}x$. $\square$

**Lemma 3.7 (second moment about the mode).** If $-\nabla^2\log\pi\succeq K\succ0$ on $\mathbb R^p$, with mode $x^\star$ and mean $\bar x$, then $E\|X-x^\star\|^2_K\le p$ and $\operatorname{Cov}(X)\preceq K^{-1}$.

*Proof.* Integration by parts gives $E[\nabla U(X)^\top(X-x^\star)]=p$, and $\nabla U(x)^\top(x-x^\star)\ge\|x-x^\star\|^2_K$. The covariance bound is Brascamp–Lieb. $\square$

**Lemma 3.8 (the mode–mean segment is admissible).** Every $c\in[\gamma^\star,\bar\gamma]$ belongs to $\mathcal C$.

*Proof.* Let $K=(1-\rho_{\rm hi})P_0$, valid by Lemma 3.3. For $c$ on the segment, $\|c-\bar\gamma\|_K\le\|\gamma^\star-\bar\gamma\|_K$, so
$$E\|\tilde\gamma-c\|^2_K=\operatorname{tr}(K\operatorname{Cov})+\|c-\bar\gamma\|^2_K\le E\|\tilde\gamma-\gamma^\star\|^2_K\le p$$
by Lemma 3.7. Divide by $1-\rho_{\rm hi}$. $\square$

**Lemma 3.9 (Gaussian smoothing).**
$$\bigl\|\mathcal N(\mu,P_0^{-1})-\mathcal N(\tilde\mu,P_0^{-1})\bigr\|_{TV}=2\Phi\bigl(\tfrac12\|\mu-\tilde\mu\|_{P_0}\bigr)-1\le\frac{\|\mu-\tilde\mu\|_{P_0}}{\sqrt{2\pi}}.\qquad\square$$

**Lemma 3.10 (spectral bound from contraction).** If $W_1(K(x,\cdot),K(y,\cdot))\le\rho\|x-y\|_{P_0}$ for all $x,y$, then the $L^2(\pi_\gamma)$ spectral radius $\rho_\star$ of $K$ on mean-zero functions satisfies $\rho_\star\le\rho$ (Ollivier 2009, Prop. 30). $\square$

---

## 4. Main theorem

**Theorem 4.1.** Under Assumption LC:

1. **(Contraction.)** There is a coupling of one iteration from $\gamma$ and from $\tilde\gamma$ with
$$E\|\gamma'-\tilde\gamma'\|^2_{P_0}\le\rho_{\rm hi}^2\,\|\gamma-\tilde\gamma\|^2_{P_0}.$$
   If every $A_i$ has rank one, which is automatic when $p=1$ or $q=1$, the contraction holds pathwise: $\|\gamma'-\tilde\gamma'\|_{P_0}\le\rho_{\rm hi}\|\gamma-\tilde\gamma\|_{P_0}$.
2. **(Total variation.)** For every $k\ge1$ and every $c\in\mathcal C$,
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\frac{\rho_{\rm hi}^{\,k}}{\sqrt{2\pi}}\Bigl(\|\gamma^{(0)}-c\|_{P_0}+\sqrt{\frac p{1-\rho_{\rm hi}}}\Bigr),$$
   and the same holds for the joint law of $(\gamma^{(k)},\beta^{(k+1)})$. In bounded form, the right side may be replaced by $2\Phi\bigl(\tfrac12\rho_{\rm hi}^kD\bigr)-1$, with $D$ the bracket.
3. **(Rate sandwich.)** $\rho_{\rm lo}\le\rho_\star\le\rho_{\rm hi}$. The lower bound requires $C_{i,\max}<\infty$.

*Proof of 1.* Write $\Delta\gamma=\gamma-\tilde\gamma$.

- **Coupling for rank-one $A_i=a_iv_i^\top$.** The tilt is along $a_i$, with parameter $s_i=v_i^\top\gamma$ (Lemma 3.1). Draw $w_i=a_i^\top\beta_i$ by inversion with a common uniform, and $\beta_i$ given $w_i$ with common randomness. By Lemmas 3.4 and 3.5, $w_i-\tilde w_i=D_iv_i^\top\Delta\gamma$ with $0\le D_i\le a_i^\top(K_i^{\rm lo})^{-1}a_i$.
- **Coupling for other $A_i$.** Use a $W_{2,K_i^{\rm lo}}$-optimal coupling of $f_i(\cdot\mid\gamma)$ and $f_i(\cdot\mid\tilde\gamma)$. By Lemma 3.6, $E\|\Delta\beta_i\|^2_{K_i^{\rm lo}}\le\|A_i\Delta\gamma\|^2_{(K_i^{\rm lo})^{-1}}$. The rank-one coupling satisfies the same inequality pathwise.
- **Block 2** uses common normal noise, so
$$P_0^{1/2}(\gamma'-\tilde\gamma')=\sum_iP_0^{-1/2}A_i^\top\Delta\beta_i=Bv,$$
  with $v=\bigl((K_i^{\rm lo})^{1/2}\Delta\beta_i\bigr)_i$ and $B=\bigl[P_0^{-1/2}A_i^\top(K_i^{\rm lo})^{-1/2}\bigr]_i$.

Since $BB^\top=\tilde A(C_{\min})$, pathwise $\|Bv\|^2\le\rho_{\rm hi}\|v\|^2$. Taking expectations,
$$E\|\gamma'-\tilde\gamma'\|^2_{P_0}\le\rho_{\rm hi}\sum_i\|A_i\Delta\gamma\|^2_{(K_i^{\rm lo})^{-1}}=\rho_{\rm hi}\,\Delta\gamma^\top R_{\rm hi}\Delta\gamma\le\rho_{\rm hi}^2\|\Delta\gamma\|^2_{P_0}.$$
In the all-rank-one case, $P_0^{1/2}(\gamma'-\tilde\gamma')=P_0^{-1/2}\bigl[\sum_iv_iD_iv_i^\top\bigr]P_0^{-1/2}\,P_0^{1/2}\Delta\gamma$, with a symmetric matrix between $0$ and $\tilde A(C_{\min})$. This gives the pathwise statement. $\square$

*Proof of 2.* Couple the chain with a stationary copy $\tilde\gamma^{(0)}\sim\pi_\gamma$, using the coupling of part 1 at each step and common noise in the final $\gamma$-update. Let $\mathcal F$ be all randomness except the final normal draw.

1. By convexity of TV and Lemma 3.9,
$$\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le E\,\|\gamma^{(k)}-\tilde\gamma^{(k)}\|_{P_0}/\sqrt{2\pi}.$$
   (With common final noise, the difference of the conditional means equals $\gamma^{(k)}-\tilde\gamma^{(k)}$.)
2. Iterating part 1, $E\|\gamma^{(k)}-\tilde\gamma^{(k)}\|_{P_0}\le\rho_{\rm hi}^k\bigl(E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|^2_{P_0}\bigr)^{1/2}$.
3. By Minkowski's inequality and $c\in\mathcal C$, $\bigl(E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|^2_{P_0}\bigr)^{1/2}\le\|\gamma^{(0)}-c\|_{P_0}+\sqrt{p/(1-\rho_{\rm hi})}$.

The bounded form uses Jensen's inequality for the concave map $x\mapsto2\Phi(x/2)-1$. The joint statement holds because $\beta^{(k+1)}$ is drawn from the same conditional in both chains. $\square$

*Proof of 3.* *Upper bound:* $W_1\le W_2$, part 1, and Lemma 3.10. *Lower bound:* for $f$ linear, the autocorrelation lower-bounds $\rho_\star$. Summing over the coordinates of $x=P_0^{1/2}\gamma$, with $\tilde d$ the drift in these coordinates, Lemma 3.3 gives
$$\tfrac12E\bigl[(X-X')^\top(\tilde d(X)-\tilde d(X'))\bigr]\ge\rho_{\rm lo}\operatorname{tr}\operatorname{Cov}X .\qquad\square$$

**What Theorem 4.1 compares.** Parts 1 and 3 compare **rates**: $\rho^{\rm G}$ is the exact rate of a Gaussian model, so the true sampler's contraction and $L^2$ rate are bracketed by exact Gaussian rates. Part 2 is a **bound** whose rate is a Gaussian rate. It is not a comparison with the Gaussian sampler's actual distance from stationarity. It decays like $\rho^k$ with a stationary-coupling constant, whereas the exact Gaussian analysis gives $a^{2k}$ decay from the center (Nygren 2020).

**Why no small set is needed.** The contraction holds for every pair of states, because the curvature floor $C_{i,\min}+\Psi_i^{-1}$ is uniform, and the exactly Gaussian final step plays the role a minorization plays in drift–minorization arguments. The set $\mathcal C$ is not a small set. It only bounds the expected distance from the start to a stationary draw.

**Why no symmetry or variance monotonicity is needed.** Exact Gaussian analyses track the **law** of $\gamma^{(k)}$: its mean converges geometrically, and its variance increases monotonically toward the stationary variance. TV is then computed from the mean offset and the variance deficit. For non-Gaussian likelihoods, the variance of $\mathcal L(\gamma^{(k)})$ need not be monotone. From a start in a region where the conditionals of the $\beta_i$ are flat and wide, it can overshoot the stationary variance and then decrease.

The proof above never uses the variance of $\mathcal L(\gamma^{(k)})$. It tracks the **distance between two coupled copies**, and that distance shrinks by at least the factor $\rho_{\rm hi}$ at every step, monotonically, wherever the chains are. The worst case over all states, including the tails where the conditionals are widest, is already built into $\rho_{\rm hi}$ through the curvature floor: $\operatorname{Cov}(\beta_i\mid\gamma)\preceq(C_{i,\min}+\Psi_i^{-1})^{-1}$ everywhere. A start far in the tails enters only through the term $\|\gamma^{(0)}-c\|_{P_0}$. The price is that the variance-deficit structure, which gives $a^{2k}$ decay from the center in the Gaussian case, is not exploited.

---

## 5. Corollaries

**Corollary 5.1 (start in the admissible set; run length).**

1. For any start, the bound of Theorem 4.1(2) holds with $\|\gamma^{(0)}-c\|_{P_0}$ replaced by $\operatorname{dist}_{P_0}(\gamma^{(0)},\mathcal C)$.
2. If $\gamma^{(0)}\in\mathcal C$, in particular if $\gamma^{(0)}$ is the posterior mode, the posterior mean, or any point between them,
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\rho_{\rm hi}^{\,k}\sqrt{\frac p{2\pi(1-\rho_{\rm hi})}}\qquad(k\ge1).$$
3. From such a start, the TV distance is at most $\tau$ once
$$k\ge\max\Bigl\{1,\Bigl\lceil\frac{\log\bigl(\sqrt{p/(1-\rho_{\rm hi})}/z_\tau\bigr)}{\log(1/\rho_{\rm hi})}\Bigr\rceil\Bigr\},\qquad z_\tau=2\Phi^{-1}\Bigl(\frac{1+\tau}2\Bigr).$$

*Proof.* Minimize over $c\in\mathcal C$, and use Lemma 3.8 and the bounded form. $\square$

*Remark.* The bound's value is the same for every start in $\mathcal C$. This reflects the worst-case spread $p/(1-\rho_{\rm hi})$ used in the bound, not the chain. The sharper quantity $E\|\gamma^{(0)}-\tilde\gamma\|_{P_0}$ distinguishes starts: it is smaller at the posterior mean than at the mode when $\pi_\gamma$ is skewed.

**Corollary 5.2 (proper priors give a rate below one).** Since $C_{i,\min}\succeq0$,
$$\rho_{\rm hi}\le\rho^{\rm G}(0)=1-\lambda_{\min}\bigl(P_0^{-1/2}\Lambda P_0^{-1/2}\bigr),$$
which is strictly below $1$ whenever $\Lambda\succ0$, for every hierarchical model with log-concave group likelihoods. For scalar $\gamma$, $\rho^{\rm G}(0)=1-\Lambda/P_0$. $\square$

**Corollary 5.3 (zero curvature floor).** If $\inf_\beta\lambda_{\min}(-\nabla^2\ell_i(\beta))=0$, as for likelihoods whose curvature tends to $0$ in some direction, then $C_{i,\min}=0$ is the best global floor, and Theorem 4.1 gives only the prior-only rate $\rho^{\rm G}(0)$. That rate is close to $1$ when the prior on $\gamma$ is weak. $\square$

### 5.4 Starting values without the marginal mode

The marginal mode $\gamma^\star$ of $\pi_\gamma$ requires integrating out the $\beta_i$. The joint mode $(\hat\gamma,\hat\beta)$ is easier, since iterated conditional modes (ICM) need only one-dimensional or $q$-dimensional optimizations. This subsection bounds the starting constant from $\hat\gamma$, or from any other computed point. Write $K=(1-\rho_{\rm hi})P_0$, so that $-\nabla^2\log\pi_\gamma\succeq K$ (Lemma 3.3), and note $\|v\|^2_{P_0}=\|v\|^2_K/(1-\rho_{\rm hi})$.

**Lemma 5.4 (gradient certificate).** For any point $c\in\mathbb R^p$, let $G(c)=\|\nabla\log\pi_\gamma(c)\|_{K^{-1}}$. Then
$$E_{\pi_\gamma}\|\tilde\gamma-c\|^2_K\ \le\ t(c)^2,\qquad t(c)=\frac{G(c)+\sqrt{G(c)^2+4p}}2 .$$
At the mode $G=0$ and this is Lemma 3.7. In general $t(c)\le\sqrt p+G(c)$.

*Proof.* Let $U=-\log\pi_\gamma$. Integration by parts gives $E[\nabla U(X)^\top(X-c)]=p$ for **every** $c$, not only the mode. Strong convexity gives $(\nabla U(x)-\nabla U(c))^\top(x-c)\ge\|x-c\|^2_K$. Hence
$$E\|X-c\|^2_K\le p-\nabla U(c)^\top(\bar\gamma-c)\le p+G(c)\,\|\bar\gamma-c\|_K\le p+G(c)\,\bigl(E\|X-c\|^2_K\bigr)^{1/2},$$
using Jensen in the last step. Solve the quadratic inequality for $(E\|X-c\|^2_K)^{1/2}$. $\square$

The gradient is computable from conditional means alone (Tweedie):
$$\nabla\log\pi_\gamma(c)=-\Lambda(c-\mu_\gamma)+\sum_iA_i^\top\bigl(E[\beta_i\mid c]-H_ic\bigr).$$
For scalar $\beta_i$ these are $J$ one-dimensional integrals.

**Corollary 5.5 (start at any computed point).** For a fixed start $\gamma^{(0)}$ and $k\ge1$,
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le2\Phi\bigl(\tfrac12\rho_{\rm hi}^{\,k}D\bigr)-1,\qquad D=\frac{t(\gamma^{(0)})}{\sqrt{1-\rho_{\rm hi}}}.$$

*Proof.* The proof of Theorem 4.1(2), with $\bigl(E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|^2_{P_0}\bigr)^{1/2}\le D$ by Lemma 5.4. $\square$

No joint mode, ICM tolerance or reference point $c\in\mathcal C$ is needed: the certificate is evaluated at whatever point is used. Taking one or two Newton steps on $\log\pi_\gamma$ from $\hat\gamma$ drives $G$ toward $0$, so that $D\to\sqrt{p/(1-\rho_{\rm hi})}$, the constant for a start on the mode–mean segment. Each step uses $E[\beta_i\mid c]$ and $\operatorname{Cov}(\beta_i\mid c)$.

**Proposition 5.6 (a priori bound at the joint mode).** Let $(\hat\gamma,\hat\beta)$ be the joint mode of the posterior. Then
$$\|\nabla\log\pi_\gamma(\hat\gamma)\|^2_{P_0^{-1}}\le\rho_{\rm hi}\sum_iq_i,\qquad\text{so}\qquad G(\hat\gamma)\le\sqrt{\frac{\rho_{\rm hi}\sum_iq_i}{1-\rho_{\rm hi}}} .$$

*Proof.*
1. *Gradient at the joint mode.* Joint stationarity in $\gamma$ gives $-\Lambda(\hat\gamma-\mu_\gamma)+\sum_iA_i^\top(\hat\beta_i-H_i\hat\gamma)=0$, where $\hat\beta_i$ is the mode of $f_i(\cdot\mid\hat\gamma)$. Subtracting this from the gradient formula above, $\nabla\log\pi_\gamma(\hat\gamma)=\sum_iA_i^\top\delta_i$, with $\delta_i=E[\beta_i\mid\hat\gamma]-\hat\beta_i$, the mean minus the mode of group $i$'s conditional.
2. *Each mean–mode offset.* $f_i(\cdot\mid\hat\gamma)$ is $K_i^{\rm lo}$-strongly log-concave (Lemma 3.1). Lemma 3.7 in dimension $q_i$ and Jensen give $\|\delta_i\|_{K_i^{\rm lo}}\le\sqrt{q_i}$.
3. *Sum over groups.* $P_0^{-1/2}\nabla\log\pi_\gamma(\hat\gamma)=Bv$, with $v_i=(K_i^{\rm lo})^{1/2}\delta_i$ and $BB^\top=\tilde A(C_{\min})$. So $\|Bv\|^2\le\rho_{\rm hi}\sum_i\|\delta_i\|^2_{K_i^{\rm lo}}$. $\square$

*Remarks on Proposition 5.6.*

1. **The bound is valid but conservative.** Step 2 is nearly sharp in the worst case: a half-Gaussian attains $\|\delta\|_K=\sqrt{2/\pi}$ in one dimension. But that is a density with unbounded skewness. For GLM conditionals the offset is driven by the skewness; to first order, $\|\delta_i\|_{K}\approx|\varepsilon_i|/2$ (formal). The a priori bound therefore overstates the gradient by an order of magnitude or more, and through $D\approx\|g\|_{P_0^{-1}}/(1-\rho_{\rm hi})$ that error is multiplied by $1/(1-\rho_{\rm hi})$.
2. **Its $J$-dependence is real.** When the groups' skews have the same sign, $\hat\gamma-\gamma^\star$ is of order one in absolute terms, while the posterior s.d. of $\gamma$ shrinks like $J^{-1/2}$. So in posterior units the joint mode does drift away from the marginal mode as $J$ grows. The cost is only logarithmic in the run length.
3. **Upper curvature bounds sharpen step 2.** If $C_{i,\max}<\infty$, Cramér–Rao gives $\operatorname{Cov}\succeq(K_i^{\rm hi})^{-1}$, so $\|\delta_i\|^2_{K_i^{\rm lo}}\le q_i-\operatorname{tr}\bigl(K_i^{\rm lo}(K_i^{\rm hi})^{-1}\bigr)$. That is a modest gain.

**Proposition 5.6′ (sharper a priori offsets).** Let $f$ be a density on $\mathbb R^q$ with potential $U$, mode $m$ and mean $\mu$, and $K_{\rm lo}\preceq\nabla^2U\preceq K_{\rm hi}$. Let $\delta=\mu-m$, $\bar M=\tfrac12(K_{\rm lo}+K_{\rm hi})$ and $D=\tfrac12(K_{\rm hi}-K_{\rm lo})$.

(a) *Curvature-ratio bound.*
$$\|\delta\|_{\bar M}\le\theta\sqrt{c\,q},\qquad\theta=\lambda_{\max}\bigl(\bar M^{-1/2}D\bar M^{-1/2}\bigr),\quad c=\lambda_{\max}\bigl(K_{\rm lo}^{-1/2}\bar MK_{\rm lo}^{-1/2}\bigr).$$
For $q=1$ this reads $\sqrt{\kappa_{\rm lo}}\,|\delta|\le(\kappa_{\rm hi}-\kappa_{\rm lo})/(\kappa_{\rm hi}+\kappa_{\rm lo})$. The offset vanishes as the curvature range closes, as it must for Gaussians.

(b) *Third-derivative bound.* Let $A=\nabla^2U(m)$ and suppose $|D^3U(y)[u,u,v]|\le L\,\|u\|^2_{K_{\rm lo}}\|v\|_A$ for all $y,u,v$. Then $\|\delta\|_A\le Lq/2$. For $q=1$, with $L_3=\sup|U'''|$ and $a=U''(m)$,
$$|\delta|\le\frac{L_3}{2\,a\,\kappa_{\rm lo}} .$$

*Proof.* Both parts start from $E\nabla U(X)=0$ (integration by parts). Write $h=X-m$.

- (a) $\nabla U(X)=M(X)h$ with $M(X)=\int_0^1\nabla^2U(m+th)\,dt\in[K_{\rm lo},K_{\rm hi}]$. So $\bar M\delta=-E[(M(X)-\bar M)h]$. The matrix $N=\bar M^{-1/2}(M-\bar M)\bar M^{-1/2}$ lies between $\pm\bar M^{-1/2}D\bar M^{-1/2}$, so its spectral norm is at most $\theta$. Hence $\|\delta\|_{\bar M}\le\theta\,E\|h\|_{\bar M}\le\theta\,(c\,E\|h\|^2_{K_{\rm lo}})^{1/2}\le\theta\sqrt{cq}$, by Lemma 3.7.
- (b) $\nabla U(X)=Ah+R$, with $v^\top R=\int_0^1(1-s)\,D^3U(m+sh)[h,h,v]\,ds$, so $|v^\top R|\le\tfrac12L\|h\|^2_{K_{\rm lo}}\|v\|_A$. Then $A\delta=-ER$, and $\|\delta\|_A=\sup_{\|v\|_A\le1}v^\top A\delta\le\tfrac12L\,E\|h\|^2_{K_{\rm lo}}\le\tfrac12Lq$. $\square$

*Use at the ICM point.* Apply Proposition 5.6′ to each $f_i(\cdot\mid\hat\gamma)$, whose mode $\hat\beta_i$ is produced by ICM; take the smaller of (a) and (b) for each group. Then aggregate with the **triangle inequality**,
$$\|\nabla\log\pi_\gamma(\hat\gamma)\|_{P_0^{-1}}\le\sum_i\|P_0^{-1/2}A_i^\top\|\,b_i,$$
with $b_i$ the bound on $\|\delta_i\|$ in the matching norm, or with the Cauchy–Schwarz form of Proposition 5.6, whichever is smaller. Everything is evaluated at the ICM point, with no integrals.

*Logit constants* (ridge weight $\omega$, centred at $\hat\beta_i$, scalar $\beta_i$):

- $\kappa_{\rm lo}=\omega c_i^Q+\Psi_i^{-1}$;
- $\kappa_{\rm hi}=(1-\omega)n_i/4+\omega c_i^Q+\Psi_i^{-1}$;
- $a_i=n_iw(\hat\beta_i)+\Psi_i^{-1}$;
- $L_3=(1-\omega)\,n_i\sup_p|w(1-2p)|=(1-\omega)\,n_i/(6\sqrt3)$.

*What remains loose.*

1. In (b), $E h^2\le1/\kappa_{\rm lo}$ uses the global floor, which for small $\omega$ is close to the prior precision alone. The actual curvature near the mode is about twice that in the example.
2. The aggregation ignores signs, so there is no cancellation between groups.
3. The factor $1/(1-\rho_{\rm hi})$ from the global floor of $\pi_\gamma$ multiplies whatever gradient bound is used. Only the exact gradient (Lemma 5.4) avoids paying it on a loose bound.

**Check** (`start_values/apriori_bounds.py`; same model, data and ridge as the check below, $\omega=0.01$, $\tau=0.01$). Per-group offsets in $\sqrt{\kappa_{\rm lo}}$ units: actual maximum $0.039$; Cursor's bound $1$; ratio bound (a) $0.38$; third-derivative bound (b) $0.074$–$0.085$.

| a priori gradient bound at the ICM point | $H_i$ on $[-1,1]$: $\|g\|$ bound | $D$ | iterations | $H_i\equiv1$: iterations |
|---|---|---|---|---|
| Proposition 5.6 (Cursor: offsets $\le1$, Cauchy–Schwarz) | 4.445 | 373.6 | 803 | 823 |
| (a) ratio bound, Cauchy–Schwarz | 1.687 | 142.3 | 722 | 741 |
| (b) third-derivative bound, Cauchy–Schwarz | 0.347 | 31.8 | 597 | 613 |
| (b) third-derivative bound, triangle inequality | 0.301 | 28.3 | 587 | 613 |
| exact gradient (quadrature, not a priori) | 0.047 | 11.4 | 511 | 529 |

**Proposition 5.7 (random start).** Let $\gamma^{(0)}$ be random, independent of the stationary copy. Then $E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|^2_{P_0}=E\|\gamma^{(0)}-\bar\gamma\|^2_{P_0}+\operatorname{tr}(P_0\operatorname{Cov}_{\pi_\gamma})$, and
$$\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le2\Phi\Bigl(\tfrac12\rho_{\rm hi}^{\,k}\Bigl(E\|\gamma^{(0)}-\bar\gamma\|^2_{P_0}+\frac p{1-\rho_{\rm hi}}-\|\gamma^\star-\bar\gamma\|^2_{P_0}\Bigr)^{1/2}\Bigr)-1 .$$

*Proof.* The cross term vanishes because $E\tilde\gamma^{(0)}=\bar\gamma$ and the two starts are independent. By Lemma 3.7, $\operatorname{tr}(K\operatorname{Cov})+\|\gamma^\star-\bar\gamma\|^2_K\le p$. Then proceed as in Theorem 4.1(2), using the bounded form. $\square$

Dropping the last term gives a slightly weaker form, $E\|\gamma^{(0)}-\bar\gamma\|^2+p/(1-\rho_{\rm hi})$. That form combines Pythagoras with Brascamp–Lieb, which is sharper than the Minkowski step of Theorem 4.1(2) for starts outside $\mathcal C$. Keeping the last term recovers Corollary 5.1(2) at the mode.

*Remark (pilot chains).* Averaging $n_p$ independent pilot chains run for $k_p$ sweeps gives a random start with $E\|\bar X-\bar\gamma\|^2_{P_0}\le(\rho_{\rm hi}^{k_p}D_p)^2+V/n_p$, where $D_p$ is the pilot's starting constant and $V=(\rho_{\rm hi}^{k_p}D_p+\sqrt{p/(1-\rho_{\rm hi})})^2$; Proposition 5.7 then applies. This is valid, but at best it recovers the constant of a start on the mode–mean segment, at a cost of $n_pk_p$ sweeps. It pays only when many main chains share one pilot. For scalar $\beta_i$, a Newton step with Lemma 5.4 reaches the same constant deterministically.

**Proposition 5.8 (Monte Carlo certificate, for vector $\beta_i$).** At a fixed point $c$, draw $M$ independent exact samples from each $f_i(\cdot\mid c)$, with no chain and no burn-in, let $\bar\beta_i$ be their averages, and set $\hat g=-\Lambda(c-\mu_\gamma)+\sum_iA_i^\top(\bar\beta_i-H_ic)$. Then
$$E\|\hat g-\nabla\log\pi_\gamma(c)\|^2_{P_0^{-1}}\le\frac{\operatorname{tr}\tilde A(C_{\min})}M\le\frac{p\,\rho_{\rm hi}}M,$$
independently of $J$. Hence, with probability at least $1-\alpha$ over the draws,
$$G(c)\le\frac{\|\hat g\|_{P_0^{-1}}+\sqrt{p\rho_{\rm hi}/(\alpha M)}}{\sqrt{1-\rho_{\rm hi}}} .$$

*Proof.* $\hat g-\nabla\log\pi_\gamma(c)=\sum_iA_i^\top(\bar\beta_i-E[\beta_i\mid c])$, a sum of independent mean-zero terms with $\operatorname{Cov}(\bar\beta_i)\preceq(K_i^{\rm lo})^{-1}/M$ (Lemma 3.2). So the second moment is $\sum_i\operatorname{tr}\bigl(P_0^{-1}A_i^\top\operatorname{Cov}(\bar\beta_i)A_i\bigr)\le\operatorname{tr}\tilde A(C_{\min})/M$. Then apply Markov's inequality. $\square$

**Procedure 5.9 (run length without the posterior mean or mode).**

1. **Start.** Compute the joint mode $\hat\gamma$ by iterated conditional modes. This is deterministic and cheap: one-dimensional or $q$-dimensional optimizations.
2. **Certificate at the start.** Compute the gradient
$$g=\nabla\log\pi_\gamma(\hat\gamma)=-\Lambda(\hat\gamma-\mu_\gamma)+\sum_iA_i^\top\bigl(E[\beta_i\mid\hat\gamma]-H_i\hat\gamma\bigr).$$
   - *Scalar $\beta_i$:* $J$ one-dimensional integrals, by quadrature, or in closed form for piecewise-quadratic likelihoods. Inflate $\|g\|$ by the quadrature error bound.
   - *Vector $\beta_i$:* Monte Carlo with exact block-1 draws at $\hat\gamma$ (Proposition 5.8).
   - *No integrals at all:* the a priori bound of Proposition 5.6. It is valid, but loose.
   - *Optional:* one Newton step, $\hat\gamma\leftarrow\hat\gamma+\bigl(P_0-\sum_iA_i^\top\operatorname{Cov}(\beta_i\mid\hat\gamma)A_i\bigr)^{-1}g$, which uses the same integrals plus conditional variances. Then recompute $g$ at the new point.
3. **Starting constant.** Compute
$$G=\frac{\|g\|_{P_0^{-1}}}{\sqrt{1-\rho_{\rm hi}}},\qquad D=\frac{G+\sqrt{G^2+4p}}{2\sqrt{1-\rho_{\rm hi}}} .$$
4. **Run length.** Run for
$$k=\max\Bigl\{1,\Bigl\lceil\frac{\log(D/z_\tau)}{\log(1/\rho_{\rm hi})}\Bigr\rceil\Bigr\},\qquad z_\tau=2\Phi^{-1}\Bigl(\frac{1+\tau}2\Bigr),$$
   started at $\hat\gamma$. Then $\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le\tau$ (Corollary 5.5).

Neither $\bar\gamma$ nor $\gamma^\star$ is ever computed. The quantity the run length needs is $E\|\tilde\gamma-\hat\gamma\|^2_{P_0}$, the expected squared distance from the start to a posterior draw, and Lemma 5.4 bounds it from the gradient at the start alone. As a by-product, $\|\hat\gamma-\bar\gamma\|_{P_0}\le D$ by Jensen. For the ridge-modified model, $\rho_{\rm hi}$ is replaced by $\rho^{\rm G}(\omega C^Q)$ and the conditional expectations are taken under the modified conditionals.

**Check** (`start_values/icm_check.py`). Logit, $J=20$, $n_i=10$, $\Psi_i=0.5$, $\Lambda=0.01$, the simulated $y_i$ of §8. Ridge with $\omega=0.01$, centred at the ICM point (§6.1), so the ICM point is the joint mode of the modified model. $\pi_\gamma$, its mode and mean, and the exact gradient are computed by quadrature. $\rho_{\rm hi}=0.9881$, $\sqrt{p/(1-\rho_{\rm hi})}=9.16$, $\tau=0.01$.

| starting constant $D$ | $H_i$ on $[-1,1]$: $D$ | iterations | $H_i\equiv1$: $D$ | iterations |
|---|---|---|---|---|
| ICM point, Prop. 5.6 a priori gradient, via the mode–mean triangle inequality and Prop. 5.7 | 382.6 | 805 | 391.5 | 826 |
| ICM point, Prop. 5.6 a priori gradient, Lemma 5.4 | 373.6 | 803 | 382.3 | 823 |
| ICM point, **exact** gradient, Lemma 5.4 | 11.4 | 511 | 12.2 | 529 |
| one marginal Newton step from ICM, Lemma 5.4 | 9.165 | 493 | 9.271 | 506 |
| start on the mode–mean segment (Corollary 5.1) | 9.164 | 493 | 9.270 | 506 |
| actual $\sqrt{E\|\tilde\gamma-\hat\gamma\|^2_{P_0}}$ (not a bound) | 1.39 | — | 1.39 | — |

The exact gradient at the ICM point has $\|g\|_{P_0^{-1}}=0.047$ (and $0.060$ for $H_i\equiv1$), against the a priori bound $4.45$. The gap between the actual starting distance $1.39$ and the best available constant $9.16$ comes from the global curvature floor in Theorem 4.1, not from the start.

---

## 6. Ridge modification

### 6.1 A floor independent of the prior

**Proposition 6.1.** For $\omega\in(0,1]$ and reference points $m_i$, let
$$\ell_i^\omega(\beta)=(1-\omega)\,\ell_i(\beta)+\omega\,T_{m_i}(\beta),$$
where $T_{m_i}$ is the second-order Taylor polynomial of $\ell_i$ at $m_i$, and $C_i^Q=-\nabla^2\ell_i(m_i)$. Then $\ell_i^\omega$ satisfies Assumption LC with
$$C^\omega_{i,\min}=(1-\omega)C_{i,\min}+\omega C_i^Q\succeq\omega C_i^Q,\qquad C^\omega_{i,\max}=(1-\omega)C_{i,\max}+\omega C_i^Q .$$
So Theorem 4.1 applies to the modified model with rate at most $\rho^{\rm G}(\omega C^Q)$. $\square$

A natural deterministic choice of $m_i$ is the mode of $\beta_i$ given $\eta_i=H_i\gamma^\star$, with $\gamma^\star$ a mode of $\pi_\gamma$.

**Flat-prior limit.** For scalar $\gamma$ and $\beta_i$, as $\Lambda\to0$,
$$\rho^{\rm G}(\omega C^Q)\to\frac{\sum_i(H_i^2/\Psi_i)\,(1+\omega c_i^Q\Psi_i)^{-1}}{\sum_iH_i^2/\Psi_i}\le\max_i\frac1{1+\omega c_i^Q\Psi_i}<1 .$$

**Small $\omega$.**
$$1-\rho^{\rm G}(\omega C^Q)\approx\frac{\Lambda+\omega\sum_iH_i^2c_i^Q}{P_0}.$$
The ridge acts like adding prior precision equal to a fraction $\omega$ of the data's information about $\gamma$, evaluated at the reference points. Unlike $\Lambda$, it does not vanish when the prior becomes flat.

**Centring at the joint mode.** Since $\nabla T_{m_i}(m_i)=\nabla\ell_i(m_i)$, the ridge leaves the log-likelihood gradient unchanged at the centre. If the centres are the $\beta$-components of the joint mode of the original model, that point is also the joint mode of the modified model, and Proposition 5.6 applies to it. With the gradient certificate (Lemma 5.4) this requirement is not needed, but the conditional means must then be computed under the modified model.

**Sampling.** $e^{\ell_i^\omega}\phi$ is log-concave, so block 1 can be drawn exactly by likelihood-subgradient accept–reject (Nygren and Nygren 2006).

### 6.2 Distance between the true and modified posteriors

Let $d_i=T_{m_i}-\ell_i$ and $S(\beta)=\sum_id_i(\beta_i)$. Then exactly
$$\pi^\omega(\gamma,\beta)\propto\pi(\gamma,\beta)\,e^{\omega S(\beta)},\qquad\pi^\omega_\gamma(\gamma)\propto\pi_\gamma(\gamma)\,E_\pi[e^{\omega S}\mid\gamma].$$

**Proposition 6.2 (first-order characterization).** If $E_\pi e^{\omega_0S}<\infty$ for some $\omega_0>0$, then as $\omega\to0$
$$\|\pi-\pi^\omega\|_{TV}=\tfrac\omega2E_\pi|S-E_\pi S|+O(\omega^2),\qquad\|\pi_\gamma-\pi^\omega_\gamma\|_{TV}=\tfrac\omega2E_\pi\bigl|E_\pi[S\mid\gamma]-E_\pi S\bigr|+O(\omega^2).$$

*Proof.* Expand $e^{\omega S}/E_\pi e^{\omega S}$ to first order in $L^1(\pi)$. $\square$

**Proposition 6.3 (rigorous bounds).** Let $\Lambda_S(\omega)=\log E_\pi e^{\omega(S-E_\pi S)}$. Then $\|\pi_\gamma-\pi^\omega_\gamma\|_{TV}\le\|\pi-\pi^\omega\|_{TV}$, and:

(a) *In general,*
$$\|\pi-\pi^\omega\|_{TV}\le\frac\omega2\Bigl(\operatorname{Var}_\pi S+(\Lambda_S(\omega)/\omega)^2\Bigr)^{1/2}\Bigl(1+E_\pi e^{2\omega(S-E_\pi S)}\Bigr)^{1/2}.$$

(b) *If $S\le M$,* with $\Delta=M-E_\pi S$,
$$\|\pi-\pi^\omega\|_{TV}\le\frac\omega2e^{\omega\Delta}\Bigl(E_\pi|S-E_\pi S|+\frac\omega2e^{\omega\Delta}\operatorname{Var}_\pi S\Bigr).$$

*Proof.* Write the distance as $\tfrac12E_\pi|e^{\omega(S-c)}-1|$ with $c=\omega^{-1}\log E_\pi e^{\omega S}\ge E_\pi S$. Use $|e^x-1|\le|x|e^{x^+}$, then Cauchy–Schwarz for (a). For (b), use $\omega(S-c)\le\omega\Delta$ and $e^x\le1+x+\tfrac12x^2e^{x^+}$. $\square$

**Theorem 6.4 (end-to-end bound).** For the sampler targeting $\pi^\omega$, with exact block-1 draws and a start in the admissible set of the modified model,
$$\bigl\|\mathcal L^\omega(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\rho^{\rm G}(\omega C^Q)^{\,k}\sqrt{\frac p{2\pi\bigl(1-\rho^{\rm G}(\omega C^Q)\bigr)}}+\|\pi_\gamma-\pi^\omega_\gamma\|_{TV}.$$

*Proof.* Triangle inequality, Corollary 5.1 for the modified model, and Proposition 6.3. $\square$

The first term decreases in $\omega$ and the second increases, about linearly. This is a bound on the distance to the **true** posterior, given the bias term.

*Remark (reweighting).* Since $\pi\propto\pi^\omega e^{-\omega S}$, draws from the modified sampler can be reweighted by $e^{-\omega S}$ to estimate true-posterior expectations.

---

## 7. Families

### 7.1 Binomial logit

$\ell_i(\beta)=y_i\beta-n_i\log(1+e^\beta)$ and $-\ell_i''=n_iw(\beta)\in(0,n_i/4]$, with $w=p(1-p)$.

- **Global floor $0$:** Corollary 5.3 applies.
- **Upper bound $n_i/4$** (Böhning): $\rho_\star\ge\rho^{\rm G}(n/4)$. No logit hierarchical model mixes faster than the Gaussian model with group precisions $n_i/4$.
- **$d_i=T_{m_i}-\ell_i$ is bounded above.** The tangent quadratic exceeds $\ell_i$ only on a bounded interval, and $d_i\to-\infty$ quadratically. So Proposition 6.3(b) applies, with $M$ computed by one-dimensional maximization.

### 7.2 Poisson

$\ell_i(\beta)=y_i\beta-n_ie^\beta$ and $-\ell_i''=n_ie^\beta\in(0,\infty)$.

- **Global floor $0$, and no upper bound.** The lower end of the rate sandwich is not available globally.
- **$d_i$ is unbounded above** ($d_i\to+\infty$ as $\beta\to+\infty$), but $E_\pi e^{2\omega S}<\infty$ for $\omega<\tfrac12$. So Propositions 6.2 and 6.3(a) apply.

---

## 8. Example

**Setting.** Logit, $J=20$ groups, $n_i=10$, $H_i$ evenly spaced on $[-1,1]$, $\Psi_i=0.5$, $\Lambda=0.01$, so $P_0=14.75$. The counts are simulated from $\gamma=0.6$, $\beta_i=0.6H_i+e_i$, $e_i\sim\mathcal N(0,0.5)$, $y_i\sim\mathrm{Bin}(10,\operatorname{expit}\beta_i)$, with NumPy seed 3:
$$y=(10,3,4,2,5,4,1,5,1,9,6,3,6,7,3,5,8,5,7,7).$$
Ridge reference points $m_i$ are the modes of $\beta_i$ given $\eta_i=0.445H_i$, with $c_i^Q=10\,w(m_i)\in[1.96,2.50]$. Starts are in $\mathcal C$, and $\tau=0.01$.

| model | rate, $\Lambda=0.01$ | iterations | rate, $\Lambda=0$ | iterations | $\|\pi_\gamma-\pi^\omega_\gamma\|_{TV}$ (exact, quadrature) |
|---|---|---|---|---|---|
| true logit ($C_{\min}=0$) | 0.9993 | 10,813 | 1 | none | 0 |
| ridge $\omega=0.01$ | 0.98810 | 494 | 0.98877 | 526 | $3.1\times10^{-4}$ |
| ridge $\omega=0.05$ | 0.94565 | 93 | 0.94629 | 94 | $1.6\times10^{-3}$ |
| ridge $\omega=0.10$ | 0.8975 | 45 | 0.8981 | 45 | $3.1\times10^{-3}$ |

The first-order prediction of Proposition 6.2 matches the exact marginal distances to within 1.3%. Scripts: `positioning_rule/ridge_rate.py` and `positioning_rule/ridge_tv.py`.

---

## 9. Discussion and open problems

1. **Comparison of distances.** The natural strengthening of Theorem 4.1 is
$$\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le\|\mathcal L^{\rm G}(\gamma^{(k)})-\pi^{\rm G}_\gamma\|_{TV}$$
   for the Gaussian model with precisions $C_{\min}$ and a matched start. It would let the exact Gaussian formulas of Nygren (2020) be used directly for log-concave likelihoods. It is not established.
2. **$a^{2k}$ decay from the center.** The coupling bound loses the variance-deficit structure that gives $a^{2k}$ decay in the Gaussian case.
3. **Local floors.** For GLMs the curvature is positive wherever the posterior mass is but tends to $0$ in the tails. A theorem using a floor over a region of high posterior probability, with explicit control of excursions, would avoid both the prior-only rate and the ridge bias. Piecewise-quadratic approximations of the likelihood (a companion line of work) are one device for this.
4. **Directional rates** for $p\ge2$ without commuting rate matrices.
5. **Random variance components** $\Psi_i$ and $\Lambda$.
6. **Tighter bias bounds** for the ridge, without true-posterior moments of $S$.

## References

- Liu, J. S., Wong, W. H. and Kong, A. (1994). Covariance structure of the Gibbs sampler with applications to the comparisons of estimators and augmentation schemes. *Biometrika* 81, 27–40.
- Madras, N. and Sezer, D. (2010). Quantitative bounds for Markov chain convergence: Wasserstein and total variation distances. *Bernoulli* 16, 882–908.
- Nygren, K. N. and Nygren, L. M. (2006). Likelihood subgradient densities. *Journal of the American Statistical Association* 101. *(Page range to be checked.)*
- Nygren, K. N. (2020). Multivariate normal distances (manuscript dated 20 July 2020). *(Final title to be filled in.)*
- Ollivier, Y. (2009). Ricci curvature of Markov chains on metric spaces. *Journal of Functional Analysis* 256, 810–864.
- Roberts, G. O. and Rosenthal, J. S. (2002). One-shot coupling for certain stochastic recursive sequences. *Stochastic Processes and their Applications* 99, 195–208.
- Rosenthal, J. S. (1995). Minorization conditions and convergence rates for Markov chain Monte Carlo. *Journal of the American Statistical Association* 90, 558–566.
- Saumard, A. and Wellner, J. A. (2014). Log-concavity and strong log-concavity: a review. *Statistics Surveys* 8, 45–114.
