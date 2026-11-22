# Explicit Convergence Bounds for Two-Block Gibbs Samplers in Hierarchical Models with Log-Concave Likelihoods

*Draft.*

---

## Abstract

We consider the two-block Gibbs sampler for hierarchical models with a Gaussian population layer and log-concave group likelihoods, as in hierarchical generalized linear models. If each group's log-likelihood has curvature bounded below by a matrix $`C_{i,\min}`$, and possibly above by $`C_{i,\max}`$, the sampler's convergence can be compared with that of two Gaussian hierarchical models: those whose group likelihoods are Gaussian with precisions $`C_{i,\min}`$ and $`C_{i,\max}`$. We prove three results.

1. A coupling contracts at the exact rate $`\rho^{\rm G}(C_{\min})`$ of the first Gaussian model.
2. The total-variation distance to stationarity is bounded explicitly at that rate, from any start, with a simplified form for starts between the posterior mode and mean.
3. The $`L^2`$ geometric rate lies between $`\rho^{\rm G}(C_{\max})`$ and $`\rho^{\rm G}(C_{\min})`$.

The rate is below one for any proper prior on the population parameter. For common GLMs the global curvature floor is zero, and the rate then degrades toward one as the prior weakens. A ridge modification, mixing each log-likelihood with its quadratic approximation at weight $`\omega`$, restores a floor that does not depend on the prior. We bound the total-variation distance between the true and modified posteriors (linear in $`\omega`$ to first order), which gives an end-to-end bound for the distance between the modified sampler's output and the true posterior.

The comparison is between **rates**, not between distances. The total-variation bound is weaker than the exact Gaussian analysis of Nygren (2020), and a comparison of actual distances with the Gaussian sampler remains open.

---

## 1. Introduction

Explicit convergence bounds for Gibbs samplers usually rest on drift and minorization conditions (Rosenthal 1995), whose constants are hard to make sharp. For Gaussian hierarchical models, the two-block sampler is a Gaussian autoregression in the population parameter, and its distance from stationarity can be computed exactly (Nygren 2020). This paper extends what can be said, explicitly and with constants, to log-concave group likelihoods.

The key structural facts are:

1. **Each group-effect conditional is an exponential tilt** of a strongly log-concave measure, with the tilt parameter linear in the population parameter $`\gamma`$.
2. **The $`\gamma`$-update is exactly Gaussian.**
3. **Strong log-concavity gives contraction.** Through quantile and transport couplings, it gives a contraction of the $`\gamma`$-chain at the rate of a Gaussian model.
4. **The Gaussian final step turns contraction into total variation.** This is a one-shot coupling argument (Roberts and Rosenthal 2002; Madras and Sezer 2010), and no small set or minorization is needed.

**How the $`k`$-step kernel is avoided.** Without conjugacy, the law of $`\gamma^{(k)}`$ has no closed form, so the exact approach of the Gaussian case (compute $`\mathcal L(\gamma^{(k)})`$, then its TV distance) is unavailable. The argument never computes it.

1. **Two copies, one inequality.** Run the chain from $`\gamma^{(0)}`$ alongside a copy started at stationarity. Because $`\pi_\gamma`$ is invariant, the copy stays at stationarity, so its law at every step is $`\pi_\gamma`$, even though $`\pi_\gamma`$ is not known in closed form.
2. **Only a one-step sensitivity bound is needed.** For block 1 the argument uses no formula for the conditional $`f_i(\cdot\mid\gamma)`$, only how much it can move when $`\gamma`$ moves. Because it is an exponential tilt, in $`\gamma`$, of a strongly log-concave base, its quantiles (Lemma 3.4) or its $`W_2`$-optimal transport (Lemma 3.6) move by at most the Gaussian factor $`(C_{i,\min}+\Psi_i^{-1})^{-1}`$ per unit tilt. Combined with the linear block-2 mean, this gives a contraction inequality, $`E\|\Delta\gamma'\|^2\le\rho^2\|\Delta\gamma\|^2`$, that holds **for every pair of states**.
3. **Inequalities compose; kernels do not.** A bound valid for every pair of states can be applied at each step regardless of where the chains are, so after $`k`$ steps the distance is at most $`\rho^k`$ times the initial one. The $`k`$-step kernel is never written down.
4. **The only exact formula used is the last Gaussian step.** Given the final $`\beta`$'s, both copies' $`\gamma^{(k)}`$ are normal with common covariance $`P_0^{-1}`$. That converts the coupling distance into total variation (Lemma 3.9).
5. **The stationary law enters only through moments.** Strong log-concavity bounds the expected distance from the start to a stationary draw (Lemmas 3.7–3.8), without knowing $`\pi_\gamma`$.

The price of not computing the kernel is precision: the result is a bound at the Gaussian rate, not the exact distance. In particular it loses the $`a^{2k}`$ decay from the center that the exact Gaussian calculation gives.

Section 2 sets up the model. Section 3 collects the lemmas. Section 4 states and proves the main theorem, and §5 its corollaries. Section 6 treats the ridge modification and its distance from the true posterior. Section 7 specializes to the logit and Poisson families, §8 gives a worked example, and §9 discusses limitations and open problems.

---

## 2. Model, sampler and assumptions

### 2.1 Model

For $`\gamma\in\mathbb R^p`$ and groups $`i=1,\dots,J`$ with $`\beta_i\in\mathbb R^q`$:
```math
\gamma\sim\mathcal N_p(\mu_\gamma,\Lambda^{-1}),\qquad\beta_i\mid\gamma\sim\mathcal N_q(H_i\gamma,\Psi_i),\qquad p(y_i\mid\beta_i)=e^{\ell_i(\beta_i)},
```
with $`H_i\in\mathbb R^{q\times p}`$, $`\Psi_i\succ0`$ and $`\Lambda\succeq0`$. Write
```math
A_i=\Psi_i^{-1}H_i,\qquad P_0=\Lambda+\sum_iH_i^\top\Psi_i^{-1}H_i,
```
and assume $`P_0\succ0`$.

**What $`P_0`$ is.** $`P_0`$ is the precision of $`\gamma`$ **given all the group effects**: the precision of the block-2 conditional. It is neither the prior precision $`\Lambda`$ nor the posterior precision of $`\gamma`$, and it does not depend on the data.

### 2.2 Sampler

One iteration from $`\gamma`$:

1. **(Block 1)** draw $`\beta_i\sim f_i(\cdot\mid\gamma)\propto e^{\ell_i(\beta)}\phi_q(\beta;H_i\gamma,\Psi_i)`$, independently over $`i`$;
2. **(Block 2)** draw $`\gamma'\sim\mathcal N_p\bigl(P_0^{-1}(\Lambda\mu_\gamma+\sum_iA_i^\top\beta_i),\ P_0^{-1}\bigr)`$.

$`K`$ denotes the kernel of the $`\gamma`$-chain and $`\pi_\gamma`$ the marginal posterior of $`\gamma`$. $`K`$ is reversible with respect to $`\pi_\gamma`$, with non-negative spectrum (Liu, Wong and Kong 1994). $`\mathcal L(\gamma^{(k)})`$ is the law after $`k`$ iterations from a fixed $`\gamma^{(0)}`$.

### 2.3 Assumption

**Assumption LC (curvature bounds).** Each $`\ell_i`$ is concave on $`\mathbb R^q`$ with locally Lipschitz gradient, and
```math
C_{i,\min}\preceq-\nabla^2\ell_i(\beta)\preceq C_{i,\max}\qquad\text{for almost every }\beta,
```
for symmetric $`0\preceq C_{i,\min}\preceq C_{i,\max}`$. The upper bound may be infinite, and is used only for lower bounds on the rate.

That is, each group's likelihood carries at least the information of a Gaussian observation with precision $`C_{i,\min}`$, and at most that of one with precision $`C_{i,\max}`$.

### 2.4 Rate quantities

**Definition 2.1.** For a collection $`C=(C_i)`$ of positive semidefinite matrices, let
```math
\tilde A(C)=P_0^{-1/2}\Bigl[\sum_iA_i^\top\bigl(C_i+\Psi_i^{-1}\bigr)^{-1}A_i\Bigr]P_0^{-1/2},\qquad\rho^{\rm G}(C)=\lambda_{\max}\bigl(\tilde A(C)\bigr).
```
$`\rho^{\rm G}(C)`$ is the exact $`L^2`$ convergence rate of the two-block sampler in the **Gaussian** hierarchical model with group likelihoods $`\mathcal N(\beta_i;\cdot,C_i^{-1})`$, that is, group-likelihood precision $`C_i`$. Write
```math
\rho_{\rm hi}=\rho^{\rm G}(C_{\min}),\qquad\rho_{\rm lo}=\lambda_{\min}\bigl(\tilde A(C_{\max})\bigr),\qquad K_i^{\rm lo}=C_{i,\min}+\Psi_i^{-1},\qquad K_i^{\rm hi}=C_{i,\max}+\Psi_i^{-1}.
```
For scalar $`\gamma`$ and $`\beta_i`$,
```math
\rho_{\rm hi}=\frac1{P_0}\sum_i\frac{H_i^2}{\Psi_i(1+c_{i,\min}\Psi_i)}.
```

**Definition 2.2 (admissible reference set).** With $`\bar\gamma`$ and $`\gamma^\star`$ the mean and mode of $`\pi_\gamma`$,
```math
\mathcal C=\Bigl\{c\in\mathbb R^p:\ E_{\pi_\gamma}\|\tilde\gamma-c\|^2_{P_0}\le\frac p{1-\rho_{\rm hi}}\Bigr\}=\Bigl\{c:\ \|c-\bar\gamma\|^2_{P_0}\le\frac p{1-\rho_{\rm hi}}-\mathrm{tr}\bigl(P_0\mathrm{Cov}_{\pi_\gamma}\bigr)\Bigr\}.
```
It is a closed ellipsoid in the $`P_0`$-metric, centered at the posterior mean, and it contains the segment $`[\gamma^\star,\bar\gamma]`$ (Lemma 3.8).

**Notation.** $`\|v\|_M=(v^\top Mv)^{1/2}`$; $`\preceq`$ is the Loewner order; $`\varphi`$ and $`\Phi`$ are the standard normal density and distribution function.

---

## 3. Lemmas

**Conventions.** For a symmetric $`K\succ0`$, a density $`e^{-U}`$ on $`\mathbb R^d`$ is **$`K`$-strongly log-concave** if $`U`$ is convex and $`x\mapsto U(x)-\tfrac12\|x\|_K^2`$ is convex. No differentiability is required. In the notation of Saumard and Wellner (2014, Definition 2.9), these are exactly the densities $`p=g\,\varphi_{\Sigma}(\cdot-\mu)`$ with $`g`$ log-concave and $`\Sigma=K^{-1}`$, the class $`\mathrm{SLC}_2(\mu,\Sigma,d)`$. If $`U`$ has a locally Lipschitz gradient, $`K`$-strong log-concavity is equivalent to $`\nabla^2U\succeq K`$ almost everywhere, since the a.e. Hessian of a $`C^{1,1}_{\rm loc}`$ function is its distributional Hessian.

**Dependencies.** The numbering is kept for cross-reference; the logical order is 3.5 → 3.7 → 3.2 → 3.3 → 3.6. Lemmas 3.1, 3.4, 3.8, 3.9 and 3.10 use only what precedes them in this order. External results are named with exact theorem numbers and their hypotheses are checked; everything else is proved here. Appendix A collects the sources.

**Lemma 3.1 (tilt representation).** With $`V_i(\beta)=-\ell_i(\beta)+\tfrac12\beta^\top\Psi_i^{-1}\beta`$,
```math
f_i(\beta\mid\gamma)\propto\exp\{\theta_i^\top\beta-V_i(\beta)\},\qquad\theta_i=A_i\gamma,
```
and under Assumption LC, $`V_i`$ is $`K_i^{\rm lo}`$-strongly convex with $`\nabla^2V_i\preceq K_i^{\rm hi}`$ almost everywhere.

*Proof.* Expand $`-\tfrac12(\beta-H_i\gamma)^\top\Psi_i^{-1}(\beta-H_i\gamma)`$ and drop terms free of $`\beta`$. The curvature statements follow from Assumption LC and the conventions above. $`\square`$

**Lemma 3.5 (one-dimensional projections).** Let $`X`$ have a $`K`$-strongly log-concave density on $`\mathbb R^d`$, and let $`a\ne0`$. Then $`W=a^\top X`$ has a $`\kappa_a`$-strongly log-concave density on $`\mathbb R`$, with $`\kappa_a=(a^\top K^{-1}a)^{-1}`$. If the density of $`X`$ is $`\propto e^{s\,a^\top x-V(x)}`$, the law of $`W`$ is a one-parameter exponential tilt in $`s`$ of a fixed base law.

*Source.* Saumard and Wellner (2014), Theorem 3.8, attributed there to Brascamp and Lieb (1976), Theorem 4.3 (p. 380), and to Simon (2011), Theorem 13.13. Their hypothesis is $`p\in\mathrm{SLC}_2(0,\Sigma,m+n)`$; their conclusion is that the marginal over the last $`n`$ coordinates lies in $`\mathrm{SLC}_2(0,\Sigma_{11},m)`$. No smoothness is assumed. (A reading of the review flagged an inconsistent dimension label in its statement; the result used here is the standard one, as in Brascamp and Lieb 1976, Theorem 4.3.)

*Proof.* Choose an invertible $`T`$ whose first row is $`a^\top`$. If $`X`$ has density $`g\,\varphi_\Sigma(\cdot-\mu)`$ with $`\Sigma=K^{-1}`$, then $`Y=TX`$ has density $`|\det T|^{-1}g(T^{-1}y)\,\varphi_\Sigma(T^{-1}y-\mu)`$. Since $`\varphi_\Sigma(T^{-1}y-\mu)\propto\varphi_{T\Sigma T^\top}(y-T\mu)`$ and $`g\circ T^{-1}`$ is log-concave, $`Y\in\mathrm{SLC}_2(T\mu,T\Sigma T^\top,d)`$. By the cited theorem, after centering, the first-coordinate marginal is in $`\mathrm{SLC}_2\bigl((T\mu)_1,(T\Sigma T^\top)_{11},1\bigr)`$, and $`(T\Sigma T^\top)_{11}=a^\top K^{-1}a`$. For the tilt: $`s\,a^\top x=s\,w`$, so integrating out the other coordinates leaves $`e^{sw}`$ times a function of $`w`$ that does not depend on $`s`$. $`\square`$

**Lemma 3.7 (second moments about any point; covariance bound).** Let $`X`$ have a $`K`$-strongly log-concave density $`e^{-U}`$ on $`\mathbb R^d`$ with mean $`\bar x`$, and let $`U`$ be locally Lipschitz (true for every finite convex function). Then:

(i) $`E\bigl[\partial U(X)^\top(X-c)\bigr]=d`$ for every $`c\in\mathbb R^d`$, where $`\partial U`$ is the a.e. gradient;

(ii) if $`x^\star`$ is a minimizer of $`U`$ (the mode), then $`E\|X-x^\star\|_K^2\le d`$;

(iii) $`\mathrm{Cov}(X)\preceq K^{-1}`$.

*Proof.*

- (i) Fix a coordinate $`j`$ and the other coordinates $`x_{-j}`$. The map $`t\mapsto e^{-U}`$ along the line is locally Lipschitz, hence absolutely continuous, with derivative $`-\partial_jU\,e^{-U}`$ almost everywhere. One-dimensional integration by parts on $`[-R,R]`$ gives $`\int_{-R}^R(t-c_j)\,\partial_jU\,e^{-U}dt=\int_{-R}^Re^{-U}dt-\bigl[(t-c_j)e^{-U}\bigr]_{-R}^R`$. Strong convexity gives $`U(x)\ge U(x^\star)+\tfrac12\|x-x^\star\|_K^2`$, so the boundary term tends to $`0`$, and dominated convergence applies. The integrand is bounded below by an integrable function, since $`(x-c)^\top\partial U(x)\ge(x-c)^\top\partial U(c)+\|x-c\|^2_K`$ by monotonicity of $`\partial U-K\,\cdot`$. Integrate over $`x_{-j}`$ by Fubini, sum over $`j`$, and normalize.
- (ii) Since $`U-\tfrac12\|\cdot\|^2_K`$ is convex, its subdifferential is monotone, which gives $`(g_x-g_y)^\top(x-y)\ge\|x-y\|_K^2`$ for subgradients $`g_x\in\partial U(x)`$, $`g_y\in\partial U(y)`$. At the mode, $`0\in\partial U(x^\star)`$, so $`\partial U(x)^\top(x-x^\star)\ge\|x-x^\star\|^2_K`$ for almost every $`x`$. Take expectations and use (i) with $`c=x^\star`$.
- (iii) For $`a\ne0`$, $`W=a^\top X`$ is $`\kappa_a`$-strongly log-concave on $`\mathbb R`$ by Lemma 3.5. Apply (ii) with $`d=1`$ to $`W`$ and its mode $`w^\star`$: $`\mathrm{Var}(W)\le E(W-w^\star)^2\le1/\kappa_a=a^\top K^{-1}a`$. $`\square`$

*Remark.* Part (iii) is the linear-test-function case of the Brascamp–Lieb variance inequality (Brascamp and Lieb 1976, Theorem 4.1; stated with $`C^2`$, strictly convex potentials in Saumard and Wellner 2014, Proposition 10.1(a), and in Bobkov and Ledoux 2000, Proposition 2.1). The proof above avoids the $`C^2`$ hypothesis of those statements, which matters for potentials that are only $`C^{1,1}`$, such as piecewise-quadratic likelihoods.

**Lemma 3.2 (covariance sandwich).** Under Assumption LC, for every $`\gamma`$,
```math
(K_i^{\rm hi})^{-1}\preceq\mathrm{Cov}(\beta_i\mid\gamma)\preceq(K_i^{\rm lo})^{-1}.
```
The upper bound needs only $`K_i^{\rm lo}`$-strong log-concavity; the lower bound needs $`C_{i,\max}<\infty`$.

*Proof.* *Upper bound:* Lemma 3.1 and Lemma 3.7(iii). *Lower bound (Cramér–Rao for the location family, proved here):* write $`U=V_i-\theta_i^\top\beta`$, so that $`\nabla U`$ is locally Lipschitz with $`\nabla^2U\preceq K_i^{\rm hi}`$ a.e. Let $`S=\nabla U(X)`$.

1. By Lemma 3.7(i) applied coordinatewise (the same one-dimensional integration by parts with $`(t-c_j)`$ replaced by $`t`$ and $`e_k`$), $`E[S]=0`$ and $`E[(X-\mu)S^\top]=I`$.
2. The same argument applied to $`\partial_kU\,e^{-U}`$ gives $`E[SS^\top]=E[\nabla^2U(X)]\preceq K_i^{\rm hi}`$. This is valid because $`\partial_kU`$ is locally Lipschitz, so $`\partial_kU\,e^{-U}`$ is absolutely continuous on lines, and $`|\nabla U|`$ grows at most linearly, since $`\nabla^2U\preceq K_i^{\rm hi}`$.
3. For any $`u,v`$: $`(u^\top v)^2=\bigl(E[u^\top(X-\mu)\,S^\top v]\bigr)^2\le u^\top\mathrm{Cov}(X)u\cdot v^\top E[SS^\top]v`$. Taking $`v=E[SS^\top]^{-1}u`$ gives $`u^\top E[SS^\top]^{-1}u\le u^\top\mathrm{Cov}(X)u`$.

So $`\mathrm{Cov}(X)\succeq E[SS^\top]^{-1}\succeq(K_i^{\rm hi})^{-1}`$. $`\square`$

**Lemma 3.3 (drift and curvature of $`\pi_\gamma`$).** With $`d(\gamma)=E[\gamma'\mid\gamma]`$ and $`R(\gamma)=\sum_iA_i^\top\mathrm{Cov}(\beta_i\mid\gamma)A_i`$,
```math
\nabla d=P_0^{-1}R(\gamma),\qquad\nabla\log\pi_\gamma(\gamma)=-\Lambda(\gamma-\mu_\gamma)+\sum_iA_i^\top\bigl(E[\beta_i\mid\gamma]-H_i\gamma\bigr),\qquad-\nabla^2\log\pi_\gamma=P_0-R(\gamma),
```
and $`R_{\rm lo}\preceq R(\gamma)\preceq R_{\rm hi}`$, where $`R_{\rm lo/hi}=\sum_iA_i^\top(K_i^{\rm hi/lo})^{-1}A_i`$. Consequently $`-\nabla^2\log\pi_\gamma\succeq K_\gamma:=P_0-R_{\rm hi}\succeq(1-\rho_{\rm hi})P_0`$.

*Proof.* Let $`m_i(\gamma)=\int e^{\ell_i(\beta)}\phi_q(\beta;H_i\gamma,\Psi_i)\,d\beta`$, so that $`\pi_\gamma\propto\phi_p(\gamma;\mu_\gamma,\Lambda^{-1})\prod_im_i(\gamma)`$. Since $`\ell_i`$ is concave, $`e^{\ell_i(\beta)}\le e^{a+b^\top\beta}`$ for some $`a,b`$. So for $`\gamma`$ in a compact set, $`\beta\mapsto e^{\ell_i}\phi_q`$ and its first two $`\gamma`$-derivatives are dominated by an integrable function of $`\beta`$, and differentiation under the integral is justified by dominated convergence. This gives $`\nabla\log m_i=A_i^\top(E[\beta_i\mid\gamma]-H_i\gamma)`$, a multivariate form of Tweedie's formula (Efron 2011, who attributes it to Robbins 1956), and $`\nabla^2\log m_i=A_i^\top\mathrm{Cov}(\beta_i\mid\gamma)A_i-H_i^\top\Psi_i^{-1}H_i`$. Summing, and adding the prior, gives the two displayed identities for $`\pi_\gamma`$. For the drift, $`d(\gamma)=P_0^{-1}\bigl(\Lambda\mu_\gamma+\sum_iA_i^\top E[\beta_i\mid\gamma]\bigr)`$; by the same dominated differentiation of the tilt family, $`\nabla_\theta E_\theta[\beta]=\mathrm{Cov}_\theta(\beta)`$. The bounds on $`R`$ follow from Lemma 3.2. Finally, $`R_{\rm hi}=P_0^{1/2}\tilde A(C_{\min})P_0^{1/2}\preceq\rho_{\rm hi}P_0`$. $`\square`$

*(Efron states the scalar normal-means form; the identity used here is proved directly above, and the citation is for attribution only.)*

**Lemma 3.4 (quantile Lipschitz bound).** Let $`g_s\propto e^{sw-U(w)}`$ on $`\mathbb R`$, with $`U-\tfrac\kappa2w^2`$ convex, $`\kappa>0`$, and let $`w_u(s)`$ be its $`u`$-quantile, $`u\in(0,1)`$. Then $`0\le\partial_sw_u\le1/\kappa`$.

*Proof (ours).* Let $`\mu_s=E_sW`$ and $`G(w)=\int_{-\infty}^w(\mu_s-t)\,g_s(t)\,dt`$. Then $`G\ge0`$, because $`G(\pm\infty)=0`$ and $`G'=(\mu_s-w)g_s`$ changes sign once. Since $`\partial_sg_s=(w-\mu_s)g_s`$, one has $`\partial_sF_s(w)=-G(w)`$, hence $`\partial_sw_u=G(w_u)/g_s(w_u)\ge0`$. Let $`H=g_s/\kappa-G`$. Almost everywhere, $`H'=g_s\,\bigl(w-(U'(w)-s)/\kappa-\mu_s\bigr)`$, and $`w\mapsto w-(U'(w)-s)/\kappa`$ is non-increasing because $`U'-\kappa w`$ is non-decreasing. So $`H'`$ is positive and then negative, with at most one sign change. Since $`H(\pm\infty)=0`$, $`H\ge0`$, that is, $`G\le g_s/\kappa`$. Differentiation in $`s`$ under the integral is justified as in Lemma 3.3. $`\square`$

**Lemma 3.6 (tilts in $`W_2`$).** Let $`\mu_\theta`$ have density $`\propto e^{\theta^\top x-V(x)}`$ on $`\mathbb R^d`$, with $`V`$ $`K`$-strongly convex. Then for all $`\theta,\theta'`$:

(i) $`\mathrm{KL}(\mu_\theta\|\mu_{\theta'})\le\tfrac12\|\theta-\theta'\|^2_{K^{-1}}`$;

(ii) $`W_{2,K}(\mu_\theta,\mu_{\theta'})\le\|\theta-\theta'\|_{K^{-1}}`$, where $`W_{2,K}`$ is the quadratic Wasserstein distance for the norm $`\|\cdot\|_K`$;

(iii) an optimal coupling attaining $`W_{2,K}`$ exists. If $`\theta`$ and $`\theta'`$ depend measurably on a parameter, an optimal coupling can be chosen measurably in that parameter.

*Proof.*

- (i) Let $`\Lambda_V(\theta)=\log\int e^{\theta^\top x-V}`$. By dominated differentiation, $`\nabla\Lambda_V(\theta)=E_\theta X`$ and $`\nabla^2\Lambda_V(\theta)=\mathrm{Cov}_\theta X\preceq K^{-1}`$, by Lemma 3.7(iii). Directly, $`\mathrm{KL}(\mu_\theta\|\mu_{\theta'})=\Lambda_V(\theta')-\Lambda_V(\theta)-\nabla\Lambda_V(\theta)^\top(\theta'-\theta)`$, a Bregman divergence. By Taylor's theorem with integral remainder it is at most $`\tfrac12\|\theta'-\theta\|^2_{K^{-1}}`$.
- (ii) *Source:* Bobkov and Ledoux (2000), Proposition 4.2. On a finite-dimensional normed space $`(E,\|\cdot\|)`$, for $`\mu=e^{-V}dx`$ with $`V`$ convex satisfying their condition (3.6) with exponent $`p=2`$ and constant $`c>0`$, every $`\nu\ll\mu`$ with finite second moment satisfies $`W_2^2(\nu,\mu)\le\tfrac2c\,D(\nu\|\mu)`$. No $`C^2`$ assumption is required. *Verification:* with $`\|\cdot\|=\|\cdot\|_K`$ and $`t+s=1`$, $`K`$-strong convexity gives $`tV(x)+sV(y)-V(tx+sy)\ge\tfrac{ts}2\|x-y\|_K^2=\tfrac12(s-s^2)\|x-y\|^2_K`$. That is condition (3.6) with $`p=2`$ and $`c=1`$. The potential of $`\mu_{\theta'}`$ is $`V-\theta'^\top x`$ plus a constant, still $`K`$-strongly convex. Take $`\nu=\mu_\theta`$ and combine with (i).
- (iii) *Source:* Villani (2009), Theorem 4.1 (existence of an optimal coupling, for Polish spaces and a lower semicontinuous cost bounded below) and Corollary 5.22 (measurable selection of optimal plans, for a continuous cost bounded below and a measurable map $`\omega\mapsto(\mu_\omega,\nu_\omega)`$). The cost $`\|x-y\|^2_K`$ is continuous and non-negative, $`\mathbb R^d`$ is Polish, and $`\theta\mapsto\mu_\theta`$ is continuous in total variation by dominated convergence, hence measurable. $`\square`$

**Lemma 3.8 (the mode–mean segment is admissible).** Every $`c\in[\gamma^\star,\bar\gamma]`$ belongs to $`\mathcal C`$.

*Proof.* By Lemma 3.3, $`\pi_\gamma`$ is $`K`$-strongly log-concave with $`K=(1-\rho_{\rm hi})P_0`$. For $`c`$ on the segment, $`\|c-\bar\gamma\|_K\le\|\gamma^\star-\bar\gamma\|_K`$, so
```math
E\|\tilde\gamma-c\|^2_K=\mathrm{tr}(K\mathrm{Cov})+\|c-\bar\gamma\|^2_K\le E\|\tilde\gamma-\gamma^\star\|^2_K\le p,
```
by Lemma 3.7(ii). Divide by $`1-\rho_{\rm hi}`$. $`\square`$

**Lemma 3.9 (Gaussian smoothing).**
```math
\bigl\|\mathcal N(\mu,P_0^{-1})-\mathcal N(\tilde\mu,P_0^{-1})\bigr\|_{TV}=2\Phi\bigl(\tfrac12\|\mu-\tilde\mu\|_{P_0}\bigr)-1\le\frac{\|\mu-\tilde\mu\|_{P_0}}{\sqrt{2\pi}} .
```
Moreover, total variation is jointly convex: for Markov kernels $`P_z,Q_z`$ and a probability $`\lambda`$, $`\bigl\|\int P_z\,\lambda(dz)-\int Q_z\,\lambda(dz)\bigr\|_{TV}\le\int\|P_z-Q_z\|_{TV}\,\lambda(dz)`$.

*Proof.* TV is invariant under the bijection $`x\mapsto P_0^{1/2}(x-\tilde\mu)`$, which maps the pair to $`\mathcal N(\Delta,I)`$ and $`\mathcal N(0,I)`$ with $`\|\Delta\|=\|\mu-\tilde\mu\|_{P_0}`$. Rotating so that $`\Delta=\|\Delta\|e_1`$, both laws share the same factor in the other coordinates. So the TV equals that of $`\mathcal N(\|\Delta\|,1)`$ and $`\mathcal N(0,1)`$, whose densities cross once at $`\|\Delta\|/2`$, giving $`\Phi(\|\Delta\|/2)-\Phi(-\|\Delta\|/2)`$. The inequality is $`\Phi(x/2)-\Phi(-x/2)\le x\,\varphi(0)`$. Convexity: for any set $`B`$, $`|\int(P_z(B)-Q_z(B))\lambda(dz)|\le\int\sup_B|P_z(B)-Q_z(B)|\,\lambda(dz)`$, and the supremum over $`B`$ is measurable in $`z`$ (for kernels on a countably generated space). $`\square`$

**Lemma 3.10 (spectral bound from contraction).** Suppose that for all $`x,y`$, $`W_1\bigl(K(x,\cdot),K(y,\cdot)\bigr)\le\rho\,\|x-y\|_{P_0}`$ with $`\rho<1`$. Then the spectral radius $`\rho_\star`$ of $`K`$ on $`L^2(\pi_\gamma)`$ restricted to mean-zero functions satisfies $`\rho_\star\le\rho`$.

*Source.* Ollivier (2009), Proposition 29, as numbered in arXiv:math/0701886v4. That statement: for a metric space with random walk $`(X,d,m)`$ and invariant distribution $`\nu`$, if the coarse Ricci curvature is at least $`\kappa>0`$ and the average spread $`\sigma<\infty`$, and $`\nu`$ is reversible (or $`X`$ is finite), then the spectral radius of the averaging operator on $`L^2(X,\nu)/\{\mathrm{const}\}`$ is at most $`1-\kappa`$. Ollivier's standing assumptions (his Definitions 1–3) are: $`(X,d)`$ Polish, $`x\mapsto m_x`$ measurable, each $`m_x`$ with finite first moment, and $`\kappa(x,y)=1-W_1(m_x,m_y)/d(x,y)`$. The average spread is $`\sigma=\|\sigma(x)\|_{L^2(\nu)}`$, with $`\sigma(x)^2=\tfrac12\iint d(y,z)^2\,dm_x(y)\,dm_x(z)`$ (his Definition 18).

*Verification.* Take $`X=\mathbb R^p`$, $`d(x,y)=\|x-y\|_{P_0}`$, $`m_x=K(x,\cdot)`$ and $`\nu=\pi_\gamma`$.

- $`\mathbb R^p`$ is Polish, and $`x\mapsto K(x,\cdot)`$ is measurable.
- $`K(x,\cdot)`$ has finite moments of all orders, since it is a Gaussian mixture with bounded mean shift.
- The hypothesis gives $`\kappa\ge1-\rho>0`$.
- $`\sigma(x)^2=\mathrm{tr}\bigl(P_0\mathrm{Cov}(\gamma'\mid\gamma=x)\bigr)`$. Since $`\mathrm{Cov}(\gamma'\mid\gamma=x)=P_0^{-1}+P_0^{-1}R(x)P_0^{-1}\preceq P_0^{-1}(1+\rho_{\rm hi})`$, by Lemma 3.3, we have $`\sigma(x)^2\le p\,(1+\rho_{\rm hi})`$ and $`\sigma<\infty`$.
- Reversibility: write $`T:L^2(\pi_\beta)\to L^2(\pi_\gamma)`$, $`(Tf)(\gamma)=E[f(\beta)\mid\gamma]`$. Then $`K=TT^*`$, where $`(T^*g)(\beta)=E[g(\gamma)\mid\beta]`$. So $`K`$ is self-adjoint, hence $`\pi_\gamma`$ is reversible, and $`\langle g,Kg\rangle=\|T^*g\|^2\ge0`$, so the spectrum is non-negative. These are standard facts for two-block samplers (Liu, Wong and Kong 1994).

*(The numbering in the published version, J. Funct. Anal., DOI 10.1016/j.jfa.2008.11.001, has not been checked here and may differ.)* $`\square`$

---

## 4. Main theorem

**Theorem 4.1.** Under Assumption LC:

1. **(Contraction.)** There is a coupling of one iteration from $`\gamma`$ and from $`\tilde\gamma`$ with
```math
E\|\gamma'-\tilde\gamma'\|^2_{P_0}\le\rho_{\rm hi}^2\,\|\gamma-\tilde\gamma\|^2_{P_0}.
```
   If every $`A_i`$ has rank one, which is automatic when $`p=1`$ or $`q=1`$, the contraction holds pathwise: $`\|\gamma'-\tilde\gamma'\|_{P_0}\le\rho_{\rm hi}\|\gamma-\tilde\gamma\|_{P_0}`$.
2. **(Total variation.)** For every $`k\ge1`$ and every $`c\in\mathcal C`$,
```math
\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\frac{\rho_{\rm hi}^{\,k}}{\sqrt{2\pi}}\Bigl(\|\gamma^{(0)}-c\|_{P_0}+\sqrt{\frac p{1-\rho_{\rm hi}}}\Bigr),
```
   and the same holds for the joint law of $`(\gamma^{(k)},\beta^{(k+1)})`$. In bounded form, the right side may be replaced by $`2\Phi\bigl(\tfrac12\rho_{\rm hi}^kD\bigr)-1`$, with $`D`$ the bracket.
3. **(Rate sandwich.)** $`\rho_{\rm lo}\le\rho_\star\le\rho_{\rm hi}`$. The lower bound requires $`C_{i,\max}<\infty`$.

*Proof of 1.* Write $`\Delta\gamma=\gamma-\tilde\gamma`$.

- **Coupling for rank-one $`A_i=a_iv_i^\top`$.** The tilt is along $`a_i`$, with parameter $`s_i=v_i^\top\gamma`$ (Lemma 3.1). Draw $`w_i=a_i^\top\beta_i`$ by inversion with a common uniform, and $`\beta_i`$ given $`w_i`$ with common randomness. By Lemmas 3.4 and 3.5, $`w_i-\tilde w_i=D_iv_i^\top\Delta\gamma`$ with $`0\le D_i\le a_i^\top(K_i^{\rm lo})^{-1}a_i`$.
- **Coupling for other $`A_i`$.** Use a $`W_{2,K_i^{\rm lo}}`$-optimal coupling of $`f_i(\cdot\mid\gamma)`$ and $`f_i(\cdot\mid\tilde\gamma)`$, chosen measurably in $`(\gamma,\tilde\gamma)`$ by Lemma 3.6(iii), so that the coupled chain is a Markov chain on pairs. By Lemma 3.6(ii), $`E\|\Delta\beta_i\|^2_{K_i^{\rm lo}}\le\|A_i\Delta\gamma\|^2_{(K_i^{\rm lo})^{-1}}`$. The rank-one coupling is an explicit measurable map of a common uniform, and satisfies the same inequality pathwise, by Lemma 3.4.
- **Block 2** uses common normal noise, so
```math
P_0^{1/2}(\gamma'-\tilde\gamma')=\sum_iP_0^{-1/2}A_i^\top\Delta\beta_i=Bv,
```
  with $`v=\bigl((K_i^{\rm lo})^{1/2}\Delta\beta_i\bigr)_i`$ and $`B=\bigl[P_0^{-1/2}A_i^\top(K_i^{\rm lo})^{-1/2}\bigr]_i`$.

Since $`BB^\top=\tilde A(C_{\min})`$, pathwise $`\|Bv\|^2\le\rho_{\rm hi}\|v\|^2`$. Taking expectations,
```math
E\|\gamma'-\tilde\gamma'\|^2_{P_0}\le\rho_{\rm hi}\sum_i\|A_i\Delta\gamma\|^2_{(K_i^{\rm lo})^{-1}}=\rho_{\rm hi}\,\Delta\gamma^\top R_{\rm hi}\Delta\gamma\le\rho_{\rm hi}^2\|\Delta\gamma\|^2_{P_0}.
```
In the all-rank-one case, $`P_0^{1/2}(\gamma'-\tilde\gamma')=P_0^{-1/2}\bigl[\sum_iv_iD_iv_i^\top\bigr]P_0^{-1/2}\,P_0^{1/2}\Delta\gamma`$, with a symmetric matrix between $`0`$ and $`\tilde A(C_{\min})`$. This gives the pathwise statement. $`\square`$

*Proof of 2.* Couple the chain with a stationary copy $`\tilde\gamma^{(0)}\sim\pi_\gamma`$, using the coupling of part 1 at each step and common noise in the final $`\gamma`$-update. Let $`\mathcal F`$ be all randomness except the final normal draw.

1. By convexity of TV and Lemma 3.9,
```math
\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le E\,\|\gamma^{(k)}-\tilde\gamma^{(k)}\|_{P_0}/\sqrt{2\pi}.
```
   (With common final noise, the difference of the conditional means equals $`\gamma^{(k)}-\tilde\gamma^{(k)}`$.)
2. Iterating part 1, $`E\|\gamma^{(k)}-\tilde\gamma^{(k)}\|_{P_0}\le\rho_{\rm hi}^k\bigl(E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|^2_{P_0}\bigr)^{1/2}`$.
3. By Minkowski's inequality and $`c\in\mathcal C`$, $`\bigl(E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|^2_{P_0}\bigr)^{1/2}\le\|\gamma^{(0)}-c\|_{P_0}+\sqrt{p/(1-\rho_{\rm hi})}`$.

The bounded form uses Jensen's inequality for the concave map $`x\mapsto2\Phi(x/2)-1`$. The joint statement holds because $`\beta^{(k+1)}`$ is drawn from the same conditional in both chains. $`\square`$

*Proof of 3.* *Upper bound:* $`W_1\le W_2`$, part 1, and Lemma 3.10. *Lower bound:* for $`f`$ linear, the autocorrelation lower-bounds $`\rho_\star`$. Summing over the coordinates of $`x=P_0^{1/2}\gamma`$, with $`\tilde d`$ the drift in these coordinates, Lemma 3.3 gives
```math
\tfrac12E\bigl[(X-X')^\top(\tilde d(X)-\tilde d(X'))\bigr]\ge\rho_{\rm lo}\mathrm{tr}\mathrm{Cov}X .\qquad\square
```

**What Theorem 4.1 compares.** Parts 1 and 3 compare **rates**: $`\rho^{\rm G}`$ is the exact rate of a Gaussian model, so the true sampler's contraction and $`L^2`$ rate are bracketed by exact Gaussian rates. Part 2 is a **bound** whose rate is a Gaussian rate. It is not a comparison with the Gaussian sampler's actual distance from stationarity. It decays like $`\rho^k`$ with a stationary-coupling constant, whereas the exact Gaussian analysis gives $`a^{2k}`$ decay from the center (Nygren 2020).

**Why no small set is needed.** The contraction holds for every pair of states, because the curvature floor $`C_{i,\min}+\Psi_i^{-1}`$ is uniform, and the exactly Gaussian final step plays the role a minorization plays in drift–minorization arguments. The set $`\mathcal C`$ is not a small set. It only bounds the expected distance from the start to a stationary draw.

**Why no symmetry or variance monotonicity is needed.** Exact Gaussian analyses track the **law** of $`\gamma^{(k)}`$: its mean converges geometrically, and its variance increases monotonically toward the stationary variance. TV is then computed from the mean offset and the variance deficit. For non-Gaussian likelihoods, the variance of $`\mathcal L(\gamma^{(k)})`$ need not be monotone. From a start in a region where the conditionals of the $`\beta_i`$ are flat and wide, it can overshoot the stationary variance and then decrease.

The proof above never uses the variance of $`\mathcal L(\gamma^{(k)})`$. It tracks the **distance between two coupled copies**, and that distance shrinks by at least the factor $`\rho_{\rm hi}`$ at every step, monotonically, wherever the chains are. The worst case over all states, including the tails where the conditionals are widest, is already built into $`\rho_{\rm hi}`$ through the curvature floor: $`\mathrm{Cov}(\beta_i\mid\gamma)\preceq(C_{i,\min}+\Psi_i^{-1})^{-1}`$ everywhere. A start far in the tails enters only through the term $`\|\gamma^{(0)}-c\|_{P_0}`$. The price is that the variance-deficit structure, which gives $`a^{2k}`$ decay from the center in the Gaussian case, is not exploited.

---

## 5. Corollaries

**Corollary 5.1 (start in the admissible set; run length).**

1. For any start, the bound of Theorem 4.1(2) holds with $`\|\gamma^{(0)}-c\|_{P_0}`$ replaced by $`\mathrm{dist}_{P_0}(\gamma^{(0)},\mathcal C)`$.
2. If $`\gamma^{(0)}\in\mathcal C`$, in particular if $`\gamma^{(0)}`$ is the posterior mode, the posterior mean, or any point between them,
```math
\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\rho_{\rm hi}^{\,k}\sqrt{\frac p{2\pi(1-\rho_{\rm hi})}}\qquad(k\ge1).
```
3. From such a start, the TV distance is at most $`\tau`$ once
```math
k\ge\max\Bigl\{1,\Bigl\lceil\frac{\log\bigl(\sqrt{p/(1-\rho_{\rm hi})}/z_\tau\bigr)}{\log(1/\rho_{\rm hi})}\Bigr\rceil\Bigr\},\qquad z_\tau=2\Phi^{-1}\Bigl(\frac{1+\tau}2\Bigr).
```

*Proof.* Minimize over $`c\in\mathcal C`$, and use Lemma 3.8 and the bounded form. $`\square`$

*Remark.* The bound's value is the same for every start in $`\mathcal C`$. This reflects the worst-case spread $`p/(1-\rho_{\rm hi})`$ used in the bound, not the chain. The sharper quantity $`E\|\gamma^{(0)}-\tilde\gamma\|_{P_0}`$ distinguishes starts: it is smaller at the posterior mean than at the mode when $`\pi_\gamma`$ is skewed.

**Corollary 5.2 (proper priors give a rate below one).** Since $`C_{i,\min}\succeq0`$,
```math
\rho_{\rm hi}\le\rho^{\rm G}(0)=1-\lambda_{\min}\bigl(P_0^{-1/2}\Lambda P_0^{-1/2}\bigr),
```
which is strictly below $`1`$ whenever $`\Lambda\succ0`$, for every hierarchical model with log-concave group likelihoods. For scalar $`\gamma`$, $`\rho^{\rm G}(0)=1-\Lambda/P_0`$. $`\square`$

**Corollary 5.3 (zero curvature floor).** If $`\inf_\beta\lambda_{\min}(-\nabla^2\ell_i(\beta))=0`$, as for likelihoods whose curvature tends to $`0`$ in some direction, then $`C_{i,\min}=0`$ is the best global floor, and Theorem 4.1 gives only the prior-only rate $`\rho^{\rm G}(0)`$. That rate is close to $`1`$ when the prior on $`\gamma`$ is weak. $`\square`$

### 5.4 Starting values without the marginal mode

The marginal mode $`\gamma^\star`$ of $`\pi_\gamma`$ requires integrating out the $`\beta_i`$. The joint mode $`(\hat\gamma,\hat\beta)`$ is easier, since iterated conditional modes (ICM) need only one-dimensional or $`q`$-dimensional optimizations. This subsection bounds the starting constant from $`\hat\gamma`$, or from any other computed point. Write $`K=(1-\rho_{\rm hi})P_0`$, so that $`-\nabla^2\log\pi_\gamma\succeq K`$ (Lemma 3.3), and note $`\|v\|^2_{P_0}=\|v\|^2_K/(1-\rho_{\rm hi})`$.

**Lemma 5.4 (gradient certificate).** For any point $`c\in\mathbb R^p`$, let $`G(c)=\|\nabla\log\pi_\gamma(c)\|_{K^{-1}}`$. Then
```math
E_{\pi_\gamma}\|\tilde\gamma-c\|^2_K\ \le\ t(c)^2,\qquad t(c)=\frac{G(c)+\sqrt{G(c)^2+4p}}2 .
```
At the mode $`G=0`$ and this is Lemma 3.7. In general $`t(c)\le\sqrt p+G(c)`$.

*Proof.* Let $`U=-\log\pi_\gamma`$. Lemma 3.7(i) gives $`E[\nabla U(X)^\top(X-c)]=p`$ for **every** $`c`$, not only the mode. Strong convexity gives $`(\nabla U(x)-\nabla U(c))^\top(x-c)\ge\|x-c\|^2_K`$. Hence
```math
E\|X-c\|^2_K\le p-\nabla U(c)^\top(\bar\gamma-c)\le p+G(c)\,\|\bar\gamma-c\|_K\le p+G(c)\,\bigl(E\|X-c\|^2_K\bigr)^{1/2},
```
using Jensen in the last step. Solve the quadratic inequality for $`(E\|X-c\|^2_K)^{1/2}`$. $`\square`$

The gradient is computable from conditional means alone (Tweedie):
```math
\nabla\log\pi_\gamma(c)=-\Lambda(c-\mu_\gamma)+\sum_iA_i^\top\bigl(E[\beta_i\mid c]-H_ic\bigr).
```
For scalar $`\beta_i`$ these are $`J`$ one-dimensional integrals.

**Corollary 5.5 (start at any computed point).** For a fixed start $`\gamma^{(0)}`$ and $`k\ge1`$,
```math
\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le2\Phi\bigl(\tfrac12\rho_{\rm hi}^{\,k}D\bigr)-1,\qquad D=\frac{t(\gamma^{(0)})}{\sqrt{1-\rho_{\rm hi}}}.
```

*Proof.* The proof of Theorem 4.1(2), with $`\bigl(E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|^2_{P_0}\bigr)^{1/2}\le D`$ by Lemma 5.4. $`\square`$

No joint mode, ICM tolerance or reference point $`c\in\mathcal C`$ is needed: the certificate is evaluated at whatever point is used. Taking one or two Newton steps on $`\log\pi_\gamma`$ from $`\hat\gamma`$ drives $`G`$ toward $`0`$, so that $`D\to\sqrt{p/(1-\rho_{\rm hi})}`$, the constant for a start on the mode–mean segment. Each step uses $`E[\beta_i\mid c]`$ and $`\mathrm{Cov}(\beta_i\mid c)`$.

**Proposition 5.6 (a priori bound at the joint mode).** Let $`(\hat\gamma,\hat\beta)`$ be the joint mode of the posterior. Then
```math
\|\nabla\log\pi_\gamma(\hat\gamma)\|^2_{P_0^{-1}}\le\rho_{\rm hi}\sum_iq_i,\qquad\text{so}\qquad G(\hat\gamma)\le\sqrt{\frac{\rho_{\rm hi}\sum_iq_i}{1-\rho_{\rm hi}}} .
```

*Proof.*
1. *Gradient at the joint mode.* Joint stationarity in $`\gamma`$ gives $`-\Lambda(\hat\gamma-\mu_\gamma)+\sum_iA_i^\top(\hat\beta_i-H_i\hat\gamma)=0`$, where $`\hat\beta_i`$ is the mode of $`f_i(\cdot\mid\hat\gamma)`$. Subtracting this from the gradient formula above, $`\nabla\log\pi_\gamma(\hat\gamma)=\sum_iA_i^\top\delta_i`$, with $`\delta_i=E[\beta_i\mid\hat\gamma]-\hat\beta_i`$, the mean minus the mode of group $`i`$'s conditional.
2. *Each mean–mode offset.* $`f_i(\cdot\mid\hat\gamma)`$ is $`K_i^{\rm lo}`$-strongly log-concave (Lemma 3.1). Lemma 3.7 in dimension $`q_i`$ and Jensen give $`\|\delta_i\|_{K_i^{\rm lo}}\le\sqrt{q_i}`$.
3. *Sum over groups.* $`P_0^{-1/2}\nabla\log\pi_\gamma(\hat\gamma)=Bv`$, with $`v_i=(K_i^{\rm lo})^{1/2}\delta_i`$ and $`BB^\top=\tilde A(C_{\min})`$. So $`\|Bv\|^2\le\rho_{\rm hi}\sum_i\|\delta_i\|^2_{K_i^{\rm lo}}`$. $`\square`$

*Remarks on Proposition 5.6.*

1. **The bound is valid but conservative.** Step 2 is nearly sharp in the worst case: a half-Gaussian attains $`\|\delta\|_K=\sqrt{2/\pi}`$ in one dimension. But that is a density with unbounded skewness. For GLM conditionals the offset is driven by the skewness; to first order, $`\|\delta_i\|_{K}\approx|\varepsilon_i|/2`$ (formal). The a priori bound therefore overstates the gradient by an order of magnitude or more, and through $`D\approx\|g\|_{P_0^{-1}}/(1-\rho_{\rm hi})`$ that error is multiplied by $`1/(1-\rho_{\rm hi})`$.
2. **Its $`J`$-dependence is real.** When the groups' skews have the same sign, $`\hat\gamma-\gamma^\star`$ is of order one in absolute terms, while the posterior s.d. of $`\gamma`$ shrinks like $`J^{-1/2}`$. So in posterior units the joint mode does drift away from the marginal mode as $`J`$ grows. The cost is only logarithmic in the run length.
3. **Upper curvature bounds sharpen step 2.** If $`C_{i,\max}<\infty`$, Cramér–Rao gives $`\mathrm{Cov}\succeq(K_i^{\rm hi})^{-1}`$, so $`\|\delta_i\|^2_{K_i^{\rm lo}}\le q_i-\mathrm{tr}\bigl(K_i^{\rm lo}(K_i^{\rm hi})^{-1}\bigr)`$. That is a modest gain.

**Proposition 5.6′ (sharper a priori offsets).** Let $`f`$ be a density on $`\mathbb R^q`$ with potential $`U`$, mode $`m`$ and mean $`\mu`$, and $`K_{\rm lo}\preceq\nabla^2U\preceq K_{\rm hi}`$. Let $`\delta=\mu-m`$, $`\bar M=\tfrac12(K_{\rm lo}+K_{\rm hi})`$ and $`D=\tfrac12(K_{\rm hi}-K_{\rm lo})`$.

(a) *Curvature-ratio bound.*
```math
\|\delta\|_{\bar M}\le\theta\sqrt{c\,q},\qquad\theta=\lambda_{\max}\bigl(\bar M^{-1/2}D\bar M^{-1/2}\bigr),\quad c=\lambda_{\max}\bigl(K_{\rm lo}^{-1/2}\bar MK_{\rm lo}^{-1/2}\bigr).
```
For $`q=1`$ this reads $`\sqrt{\kappa_{\rm lo}}\,|\delta|\le(\kappa_{\rm hi}-\kappa_{\rm lo})/(\kappa_{\rm hi}+\kappa_{\rm lo})`$. The offset vanishes as the curvature range closes, as it must for Gaussians.

(b) *Third-derivative bound.* Suppose in addition that $`U\in C^3`$. Let $`A=\nabla^2U(m)`$ and suppose $`|D^3U(y)[u,u,v]|\le L\,\|u\|^2_{K_{\rm lo}}\|v\|_A`$ for all $`y,u,v`$. Then $`\|\delta\|_A\le Lq/2`$. For $`q=1`$, with $`L_3=\sup|U'''|`$ and $`a=U''(m)`$,
```math
|\delta|\le\frac{L_3}{2\,a\,\kappa_{\rm lo}} .
```

*Proof.* Both parts start from $`E\nabla U(X)=0`$, the one-dimensional integration by parts of Lemma 3.7(i) applied to $`\partial_kU\,e^{-U}`$. Write $`h=X-m`$.

- (a) $`\nabla U(X)=M(X)h`$ with $`M(X)=\int_0^1\nabla^2U(m+th)\,dt\in[K_{\rm lo},K_{\rm hi}]`$. So $`\bar M\delta=-E[(M(X)-\bar M)h]`$. The matrix $`N=\bar M^{-1/2}(M-\bar M)\bar M^{-1/2}`$ lies between $`\pm\bar M^{-1/2}D\bar M^{-1/2}`$, so its spectral norm is at most $`\theta`$. Hence $`\|\delta\|_{\bar M}\le\theta\,E\|h\|_{\bar M}\le\theta\,(c\,E\|h\|^2_{K_{\rm lo}})^{1/2}\le\theta\sqrt{cq}`$, by Lemma 3.7.
- (b) $`\nabla U(X)=Ah+R`$, with $`v^\top R=\int_0^1(1-s)\,D^3U(m+sh)[h,h,v]\,ds`$, so $`|v^\top R|\le\tfrac12L\|h\|^2_{K_{\rm lo}}\|v\|_A`$. Then $`A\delta=-ER`$, and $`\|\delta\|_A=\sup_{\|v\|_A\le1}v^\top A\delta\le\tfrac12L\,E\|h\|^2_{K_{\rm lo}}\le\tfrac12Lq`$. $`\square`$

*Use at the ICM point.* Apply Proposition 5.6′ to each $`f_i(\cdot\mid\hat\gamma)`$, whose mode $`\hat\beta_i`$ is produced by ICM; take the smaller of (a) and (b) for each group. Then aggregate with the **triangle inequality**,
```math
\|\nabla\log\pi_\gamma(\hat\gamma)\|_{P_0^{-1}}\le\sum_i\|P_0^{-1/2}A_i^\top\|\,b_i,
```
with $`b_i`$ the bound on $`\|\delta_i\|`$ in the matching norm, or with the Cauchy–Schwarz form of Proposition 5.6, whichever is smaller. Everything is evaluated at the ICM point, with no integrals.

*Logit constants* (ridge weight $`\omega`$, centred at $`\hat\beta_i`$, scalar $`\beta_i`$):

- $`\kappa_{\rm lo}=\omega c_i^Q+\Psi_i^{-1}`$;
- $`\kappa_{\rm hi}=(1-\omega)n_i/4+\omega c_i^Q+\Psi_i^{-1}`$;
- $`a_i=n_iw(\hat\beta_i)+\Psi_i^{-1}`$;
- $`L_3=(1-\omega)\,n_i\sup_p|w(1-2p)|=(1-\omega)\,n_i/(6\sqrt3)`$.

*What remains loose.*

1. In (b), $`E h^2\le1/\kappa_{\rm lo}`$ uses the global floor, which for small $`\omega`$ is close to the prior precision alone. The actual curvature near the mode is about twice that in the example.
2. The aggregation ignores signs, so there is no cancellation between groups.
3. The factor $`1/(1-\rho_{\rm hi})`$ from the global floor of $`\pi_\gamma`$ multiplies whatever gradient bound is used. Only the exact gradient (Lemma 5.4) avoids paying it on a loose bound.

**Check** (`start_values/apriori_bounds.py`; same model, data and ridge as the check below, $`\omega=0.01`$, $`\tau=0.01`$). Per-group offsets in $`\sqrt{\kappa_{\rm lo}}`$ units: actual maximum $`0.039`$; Cursor's bound $`1`$; ratio bound (a) $`0.38`$; third-derivative bound (b) $`0.074`$–$`0.085`$.

| a priori gradient bound at the ICM point | $`H_i`$ on $`[-1,1]`$: $`\Vert g\Vert `$ bound | $`D`$ | iterations | $`H_i\equiv1`$: iterations |
|---|---|---|---|---|
| Proposition 5.6 (Cursor: offsets $`\le1`$, Cauchy–Schwarz) | 4.445 | 373.6 | 803 | 823 |
| (a) ratio bound, Cauchy–Schwarz | 1.687 | 142.3 | 722 | 741 |
| (b) third-derivative bound, Cauchy–Schwarz | 0.347 | 31.8 | 597 | 613 |
| (b) third-derivative bound, triangle inequality | 0.301 | 28.3 | 587 | 613 |
| exact gradient (quadrature, not a priori) | 0.047 | 11.4 | 511 | 529 |

**Proposition 5.7 (random start).** Let $`\gamma^{(0)}`$ be random, independent of the stationary copy. Then $`E\|\gamma^{(0)}-\tilde\gamma^{(0)}\|^2_{P_0}=E\|\gamma^{(0)}-\bar\gamma\|^2_{P_0}+\mathrm{tr}(P_0\mathrm{Cov}_{\pi_\gamma})`$, and
```math
\bigl\|\mathcal L(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le2\Phi\Bigl(\tfrac12\rho_{\rm hi}^{\,k}\Bigl(E\|\gamma^{(0)}-\bar\gamma\|^2_{P_0}+\frac p{1-\rho_{\rm hi}}-\|\gamma^\star-\bar\gamma\|^2_{P_0}\Bigr)^{1/2}\Bigr)-1 .
```

*Proof.* The cross term vanishes because $`E\tilde\gamma^{(0)}=\bar\gamma`$ and the two starts are independent. By Lemma 3.7, $`\mathrm{tr}(K\mathrm{Cov})+\|\gamma^\star-\bar\gamma\|^2_K\le p`$. Then proceed as in Theorem 4.1(2), using the bounded form. $`\square`$

Dropping the last term gives a slightly weaker form, $`E\|\gamma^{(0)}-\bar\gamma\|^2+p/(1-\rho_{\rm hi})`$. That form combines Pythagoras with the covariance bound of Lemma 3.7(iii), which is sharper than the Minkowski step of Theorem 4.1(2) for starts outside $`\mathcal C`$. Keeping the last term recovers Corollary 5.1(2) at the mode.

*Remark (pilot chains).* Averaging $`n_p`$ independent pilot chains run for $`k_p`$ sweeps gives a random start with $`E\|\bar X-\bar\gamma\|^2_{P_0}\le(\rho_{\rm hi}^{k_p}D_p)^2+V/n_p`$, where $`D_p`$ is the pilot's starting constant and $`V=(\rho_{\rm hi}^{k_p}D_p+\sqrt{p/(1-\rho_{\rm hi})})^2`$; Proposition 5.7 then applies. This is valid, but at best it recovers the constant of a start on the mode–mean segment, at a cost of $`n_pk_p`$ sweeps. It pays only when many main chains share one pilot. For scalar $`\beta_i`$, a Newton step with Lemma 5.4 reaches the same constant deterministically.

**Proposition 5.8 (Monte Carlo certificate, for vector $`\beta_i`$).** At a fixed point $`c`$, draw $`M`$ independent exact samples from each $`f_i(\cdot\mid c)`$, with no chain and no burn-in, let $`\bar\beta_i`$ be their averages, and set $`\hat g=-\Lambda(c-\mu_\gamma)+\sum_iA_i^\top(\bar\beta_i-H_ic)`$. Then
```math
E\|\hat g-\nabla\log\pi_\gamma(c)\|^2_{P_0^{-1}}\le\frac{\mathrm{tr}\tilde A(C_{\min})}M\le\frac{p\,\rho_{\rm hi}}M,
```
independently of $`J`$. Hence, with probability at least $`1-\alpha`$ over the draws,
```math
G(c)\le\frac{\|\hat g\|_{P_0^{-1}}+\sqrt{p\rho_{\rm hi}/(\alpha M)}}{\sqrt{1-\rho_{\rm hi}}} .
```

*Proof.* $`\hat g-\nabla\log\pi_\gamma(c)=\sum_iA_i^\top(\bar\beta_i-E[\beta_i\mid c])`$, a sum of independent mean-zero terms with $`\mathrm{Cov}(\bar\beta_i)\preceq(K_i^{\rm lo})^{-1}/M`$ (Lemma 3.2). So the second moment is $`\sum_i\mathrm{tr}\bigl(P_0^{-1}A_i^\top\mathrm{Cov}(\bar\beta_i)A_i\bigr)\le\mathrm{tr}\tilde A(C_{\min})/M`$. Then apply Markov's inequality. $`\square`$

**Procedure 5.9 (run length without the posterior mean or mode).**

1. **Start.** Compute the joint mode $`\hat\gamma`$ by iterated conditional modes. This is deterministic and cheap: one-dimensional or $`q`$-dimensional optimizations.
2. **Certificate at the start.** Compute the gradient
```math
g=\nabla\log\pi_\gamma(\hat\gamma)=-\Lambda(\hat\gamma-\mu_\gamma)+\sum_iA_i^\top\bigl(E[\beta_i\mid\hat\gamma]-H_i\hat\gamma\bigr).
```
   - *Scalar $`\beta_i`$:* $`J`$ one-dimensional integrals, by quadrature, or in closed form for piecewise-quadratic likelihoods. Inflate $`\|g\|`$ by the quadrature error bound.
   - *Vector $`\beta_i`$:* Monte Carlo with exact block-1 draws at $`\hat\gamma`$ (Proposition 5.8).
   - *No integrals at all:* the a priori bound of Proposition 5.6. It is valid, but loose.
   - *Optional:* one Newton step, $`\hat\gamma\leftarrow\hat\gamma+\bigl(P_0-\sum_iA_i^\top\mathrm{Cov}(\beta_i\mid\hat\gamma)A_i\bigr)^{-1}g`$, which uses the same integrals plus conditional variances. Then recompute $`g`$ at the new point.
3. **Starting constant.** Compute
```math
G=\frac{\|g\|_{P_0^{-1}}}{\sqrt{1-\rho_{\rm hi}}},\qquad D=\frac{G+\sqrt{G^2+4p}}{2\sqrt{1-\rho_{\rm hi}}} .
```
4. **Run length.** Run for
```math
k=\max\Bigl\{1,\Bigl\lceil\frac{\log(D/z_\tau)}{\log(1/\rho_{\rm hi})}\Bigr\rceil\Bigr\},\qquad z_\tau=2\Phi^{-1}\Bigl(\frac{1+\tau}2\Bigr),
```
   started at $`\hat\gamma`$. Then $`\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le\tau`$ (Corollary 5.5).

Neither $`\bar\gamma`$ nor $`\gamma^\star`$ is ever computed. The quantity the run length needs is $`E\|\tilde\gamma-\hat\gamma\|^2_{P_0}`$, the expected squared distance from the start to a posterior draw, and Lemma 5.4 bounds it from the gradient at the start alone. As a by-product, $`\|\hat\gamma-\bar\gamma\|_{P_0}\le D`$ by Jensen. For the ridge-modified model, $`\rho_{\rm hi}`$ is replaced by $`\rho^{\rm G}(\omega C^Q)`$ and the conditional expectations are taken under the modified conditionals.

**Check** (`start_values/icm_check.py`). Logit, $`J=20`$, $`n_i=10`$, $`\Psi_i=0.5`$, $`\Lambda=0.01`$, the simulated $`y_i`$ of §8. Ridge with $`\omega=0.01`$, centred at the ICM point (§6.1), so the ICM point is the joint mode of the modified model. $`\pi_\gamma`$, its mode and mean, and the exact gradient are computed by quadrature. $`\rho_{\rm hi}=0.9881`$, $`\sqrt{p/(1-\rho_{\rm hi})}=9.16`$, $`\tau=0.01`$.

| starting constant $`D`$ | $`H_i`$ on $`[-1,1]`$: $`D`$ | iterations | $`H_i\equiv1`$: $`D`$ | iterations |
|---|---|---|---|---|
| ICM point, Prop. 5.6 a priori gradient, via the mode–mean triangle inequality and Prop. 5.7 | 382.6 | 805 | 391.5 | 826 |
| ICM point, Prop. 5.6 a priori gradient, Lemma 5.4 | 373.6 | 803 | 382.3 | 823 |
| ICM point, **exact** gradient, Lemma 5.4 | 11.4 | 511 | 12.2 | 529 |
| one marginal Newton step from ICM, Lemma 5.4 | 9.165 | 493 | 9.271 | 506 |
| start on the mode–mean segment (Corollary 5.1) | 9.164 | 493 | 9.270 | 506 |
| actual $`\sqrt{E\Vert \tilde\gamma-\hat\gamma\Vert ^2_{P_0}}`$ (not a bound) | 1.39 | — | 1.39 | — |

The exact gradient at the ICM point has $`\|g\|_{P_0^{-1}}=0.047`$ (and $`0.060`$ for $`H_i\equiv1`$), against the a priori bound $`4.45`$. The gap between the actual starting distance $`1.39`$ and the best available constant $`9.16`$ comes from the global curvature floor in Theorem 4.1, not from the start.

---

## 6. Ridge modification

### 6.1 A floor independent of the prior

**Proposition 6.1.** For $`\omega\in(0,1]`$ and reference points $`m_i`$, let
```math
\ell_i^\omega(\beta)=(1-\omega)\,\ell_i(\beta)+\omega\,T_{m_i}(\beta),
```
where $`T_{m_i}`$ is the second-order Taylor polynomial of $`\ell_i`$ at $`m_i`$, and $`C_i^Q=-\nabla^2\ell_i(m_i)`$. Then $`\ell_i^\omega`$ satisfies Assumption LC with
```math
C^\omega_{i,\min}=(1-\omega)C_{i,\min}+\omega C_i^Q\succeq\omega C_i^Q,\qquad C^\omega_{i,\max}=(1-\omega)C_{i,\max}+\omega C_i^Q .
```
So Theorem 4.1 applies to the modified model with rate at most $`\rho^{\rm G}(\omega C^Q)`$. $`\square`$

A natural deterministic choice of $`m_i`$ is the mode of $`\beta_i`$ given $`\eta_i=H_i\gamma^\star`$, with $`\gamma^\star`$ a mode of $`\pi_\gamma`$.

**Flat-prior limit.** For scalar $`\gamma`$ and $`\beta_i`$, as $`\Lambda\to0`$,
```math
\rho^{\rm G}(\omega C^Q)\to\frac{\sum_i(H_i^2/\Psi_i)\,(1+\omega c_i^Q\Psi_i)^{-1}}{\sum_iH_i^2/\Psi_i}\le\max_i\frac1{1+\omega c_i^Q\Psi_i}<1 .
```

**Small $`\omega`$.**
```math
1-\rho^{\rm G}(\omega C^Q)\approx\frac{\Lambda+\omega\sum_iH_i^2c_i^Q}{P_0}.
```
The ridge acts like adding prior precision equal to a fraction $`\omega`$ of the data's information about $`\gamma`$, evaluated at the reference points. Unlike $`\Lambda`$, it does not vanish when the prior becomes flat.

**Centring at the joint mode.** Since $`\nabla T_{m_i}(m_i)=\nabla\ell_i(m_i)`$, the ridge leaves the log-likelihood gradient unchanged at the centre. If the centres are the $`\beta`$-components of the joint mode of the original model, that point is also the joint mode of the modified model, and Proposition 5.6 applies to it. With the gradient certificate (Lemma 5.4) this requirement is not needed, but the conditional means must then be computed under the modified model.

**Sampling.** $`e^{\ell_i^\omega}\phi`$ is log-concave, so block 1 can be drawn exactly by likelihood-subgradient accept–reject (Nygren and Nygren 2006).

### 6.2 Distance between the true and modified posteriors

Let $`d_i=T_{m_i}-\ell_i`$ and $`S(\beta)=\sum_id_i(\beta_i)`$. Then exactly
```math
\pi^\omega(\gamma,\beta)\propto\pi(\gamma,\beta)\,e^{\omega S(\beta)},\qquad\pi^\omega_\gamma(\gamma)\propto\pi_\gamma(\gamma)\,E_\pi[e^{\omega S}\mid\gamma].
```

**Proposition 6.2 (first-order characterization).** If $`E_\pi e^{\omega_0S}<\infty`$ for some $`\omega_0>0`$, then as $`\omega\to0`$
```math
\|\pi-\pi^\omega\|_{TV}=\tfrac\omega2E_\pi|S-E_\pi S|+O(\omega^2),\qquad\|\pi_\gamma-\pi^\omega_\gamma\|_{TV}=\tfrac\omega2E_\pi\bigl|E_\pi[S\mid\gamma]-E_\pi S\bigr|+O(\omega^2).
```

*Proof.* Expand $`e^{\omega S}/E_\pi e^{\omega S}`$ to first order in $`L^1(\pi)`$. $`\square`$

**Proposition 6.3 (rigorous bounds).** Let $`\Lambda_S(\omega)=\log E_\pi e^{\omega(S-E_\pi S)}`$. Then $`\|\pi_\gamma-\pi^\omega_\gamma\|_{TV}\le\|\pi-\pi^\omega\|_{TV}`$, and:

(a) *In general,*
```math
\|\pi-\pi^\omega\|_{TV}\le\frac\omega2\Bigl(\mathrm{Var}_\pi S+(\Lambda_S(\omega)/\omega)^2\Bigr)^{1/2}\Bigl(1+E_\pi e^{2\omega(S-E_\pi S)}\Bigr)^{1/2}.
```

(b) *If $`S\le M`$,* with $`\Delta=M-E_\pi S`$,
```math
\|\pi-\pi^\omega\|_{TV}\le\frac\omega2e^{\omega\Delta}\Bigl(E_\pi|S-E_\pi S|+\frac\omega2e^{\omega\Delta}\mathrm{Var}_\pi S\Bigr).
```

*Proof.* Write the distance as $`\tfrac12E_\pi|e^{\omega(S-c)}-1|`$ with $`c=\omega^{-1}\log E_\pi e^{\omega S}\ge E_\pi S`$. Use $`|e^x-1|\le|x|e^{x^+}`$, then Cauchy–Schwarz for (a). For (b), use $`\omega(S-c)\le\omega\Delta`$ and $`e^x\le1+x+\tfrac12x^2e^{x^+}`$. $`\square`$

**Theorem 6.4 (end-to-end bound).** For the sampler targeting $`\pi^\omega`$, with exact block-1 draws and a start in the admissible set of the modified model,
```math
\bigl\|\mathcal L^\omega(\gamma^{(k)})-\pi_\gamma\bigr\|_{TV}\le\rho^{\rm G}(\omega C^Q)^{\,k}\sqrt{\frac p{2\pi\bigl(1-\rho^{\rm G}(\omega C^Q)\bigr)}}+\|\pi_\gamma-\pi^\omega_\gamma\|_{TV}.
```

*Proof.* Triangle inequality, Corollary 5.1 for the modified model, and Proposition 6.3. $`\square`$

The first term decreases in $`\omega`$ and the second increases, about linearly. This is a bound on the distance to the **true** posterior, given the bias term.

*Remark (reweighting).* Since $`\pi\propto\pi^\omega e^{-\omega S}`$, draws from the modified sampler can be reweighted by $`e^{-\omega S}`$ to estimate true-posterior expectations.

---

## 7. Families

### 7.1 Binomial logit

$`\ell_i(\beta)=y_i\beta-n_i\log(1+e^\beta)`$ and $`-\ell_i''=n_iw(\beta)\in(0,n_i/4]`$, with $`w=p(1-p)`$.

- **Global floor $`0`$:** Corollary 5.3 applies.
- **Upper bound $`n_i/4`$** (Böhning): $`\rho_\star\ge\rho^{\rm G}(n/4)`$. No logit hierarchical model mixes faster than the Gaussian model with group precisions $`n_i/4`$.
- **$`d_i=T_{m_i}-\ell_i`$ is bounded above.** The tangent quadratic exceeds $`\ell_i`$ only on a bounded interval, and $`d_i\to-\infty`$ quadratically. So Proposition 6.3(b) applies, with $`M`$ computed by one-dimensional maximization.

### 7.2 Poisson

$`\ell_i(\beta)=y_i\beta-n_ie^\beta`$ and $`-\ell_i''=n_ie^\beta\in(0,\infty)`$.

- **Global floor $`0`$, and no upper bound.** The lower end of the rate sandwich is not available globally.
- **$`d_i`$ is unbounded above** ($`d_i\to+\infty`$ as $`\beta\to+\infty`$), but $`E_\pi e^{2\omega S}<\infty`$ for $`\omega<\tfrac12`$. So Propositions 6.2 and 6.3(a) apply.

---

## 8. Example

**Setting.** Logit, $`J=20`$ groups, $`n_i=10`$, $`H_i`$ evenly spaced on $`[-1,1]`$, $`\Psi_i=0.5`$, $`\Lambda=0.01`$, so $`P_0=14.75`$. The counts are simulated from $`\gamma=0.6`$, $`\beta_i=0.6H_i+e_i`$, $`e_i\sim\mathcal N(0,0.5)`$, $`y_i\sim\mathrm{Bin}(10,\mathrm{expit}\beta_i)`$, with NumPy seed 3:
```math
y=(10,3,4,2,5,4,1,5,1,9,6,3,6,7,3,5,8,5,7,7).
```
Ridge reference points $`m_i`$ are the modes of $`\beta_i`$ given $`\eta_i=0.445H_i`$, with $`c_i^Q=10\,w(m_i)\in[1.96,2.50]`$. Starts are in $`\mathcal C`$, and $`\tau=0.01`$.

| model | rate, $`\Lambda=0.01`$ | iterations | rate, $`\Lambda=0`$ | iterations | $`\Vert \pi_\gamma-\pi^\omega_\gamma\Vert _{TV}`$ (exact, quadrature) |
|---|---|---|---|---|---|
| true logit ($`C_{\min}=0`$) | 0.9993 | 10,813 | 1 | none | 0 |
| ridge $`\omega=0.01`$ | 0.98810 | 494 | 0.98877 | 526 | $`3.1\times10^{-4}`$ |
| ridge $`\omega=0.05`$ | 0.94565 | 93 | 0.94629 | 94 | $`1.6\times10^{-3}`$ |
| ridge $`\omega=0.10`$ | 0.8975 | 45 | 0.8981 | 45 | $`3.1\times10^{-3}`$ |

The first-order prediction of Proposition 6.2 matches the exact marginal distances to within 1.3%. Scripts: `positioning_rule/ridge_rate.py` and `positioning_rule/ridge_tv.py`.

---

## 9. Discussion and open problems

1. **Comparison of distances.** The natural strengthening of Theorem 4.1 is
```math
\|\mathcal L(\gamma^{(k)})-\pi_\gamma\|_{TV}\le\|\mathcal L^{\rm G}(\gamma^{(k)})-\pi^{\rm G}_\gamma\|_{TV}
```
   for the Gaussian model with precisions $`C_{\min}`$ and a matched start. It would let the exact Gaussian formulas of Nygren (2020) be used directly for log-concave likelihoods. It is not established.
2. **$`a^{2k}`$ decay from the center.** The coupling bound loses the variance-deficit structure that gives $`a^{2k}`$ decay in the Gaussian case.
3. **Local floors.** For GLMs the curvature is positive wherever the posterior mass is but tends to $`0`$ in the tails. A theorem using a floor over a region of high posterior probability, with explicit control of excursions, would avoid both the prior-only rate and the ridge bias. Piecewise-quadratic approximations of the likelihood (a companion line of work) are one device for this.
4. **Directional rates** for $`p\ge2`$ without commuting rate matrices.
5. **Random variance components** $`\Psi_i`$ and $`\Lambda`$.
6. **Tighter bias bounds** for the ridge, without true-posterior moments of $`S`$.

## Appendix A. Sources of external results

Full statements of the cited results are in Appendix B. Each row gives the result as used, where it is used, the exact source, the source's hypotheses, and how our setting meets them. Results not listed are proved in the text, or are standard inequalities (Jensen, Minkowski, Cauchy–Schwarz, Markov, the triangle inequality).

| result used | where | exact source | hypotheses in the source | verification here |
|---|---|---|---|---|
| Marginals of strongly log-concave densities are strongly log-concave, with the marginal covariance parameter | Lemma 3.5, hence 3.7(iii), 3.2, 3.3, 3.6 | Saumard and Wellner (2014), Theorem 3.8; attributed there to Brascamp and Lieb (1976), Theorem 4.3, p. 380, and to Simon (2011), Theorem 13.13 | $`p\in\mathrm{SLC}_2(0,\Sigma,m+n)`$, i.e. $`p=g\varphi_\Sigma`$ with $`g`$ log-concave; no smoothness | $`K`$-strong log-concavity is $`\mathrm{SLC}_2`$ with $`\Sigma=K^{-1}`$ (Saumard and Wellner, Definition 2.9); the affine change of variables is in the proof of Lemma 3.5 |
| Transportation inequality $`W_2^2\le\tfrac2c\,D`$ for uniformly convex potentials, arbitrary norms | Lemma 3.6(ii), hence Theorem 4.1(1) when $`\mathrm{rank}A_i\ge2`$ | Bobkov and Ledoux (2000), Proposition 4.2 | $`\mu=e^{-V}dx`$ on a finite-dimensional normed space, $`V`$ convex, condition (3.6) with $`p\ge2`$, $`c>0`$; $`\nu\ll\mu`$ with finite $`p`$-th moment; no $`C^2`$ | (3.6) holds with $`p=2`$, $`c=1`$ in the norm $`\Vert \cdot\Vert _K`$, by $`K`$-strong convexity; tilts have Gaussian tails |
| Existence of optimal couplings | Lemma 3.6(iii) | Villani (2009), Theorem 4.1 | Polish spaces; cost lower semicontinuous and bounded below | $`\mathbb R^q`$; cost $`\Vert x-y\Vert ^2_K\ge0`$, continuous |
| Measurable selection of optimal plans | Lemma 3.6(iii), Theorem 4.1(1) | Villani (2009), Corollary 5.22 | Polish spaces; cost continuous and bounded below; $`\omega\mapsto(\mu_\omega,\nu_\omega)`$ measurable | $`(\gamma,\tilde\gamma)\mapsto(f_i(\cdot\mid\gamma),f_i(\cdot\mid\tilde\gamma))`$ is continuous in TV by dominated convergence |
| $`W_1`$-contraction implies an $`L^2`$ spectral-radius bound for reversible chains | Lemma 3.10, Theorem 4.1(3) | Ollivier (2009), Proposition 29 (arXiv:math/0701886v4 numbering); definitions in Definitions 1–3 and 18 | Polish space; $`x\mapsto m_x`$ measurable with finite first moments; coarse Ricci curvature $`\ge\kappa>0`$; average spread $`\sigma<\infty`$; $`\nu`$ reversible | Checked in the proof of Lemma 3.10: $`\kappa\ge1-\rho`$, $`\sigma^2\le p(1+\rho_{\rm hi})`$, reversibility from $`K=TT^*`$ |
| Brascamp–Lieb variance inequality | not needed (the linear case is Lemma 3.7(iii)) | Brascamp and Lieb (1976), Theorem 4.1; Saumard and Wellner (2014), Proposition 10.1(a); Bobkov and Ledoux (2000), Proposition 2.1 | $`V\in C^2`$, strictly convex | not used; Lemma 3.7(iii) proves the needed case without $`C^2`$ |
| Tweedie's formula (attribution only) | Lemma 3.3 | Efron (2011); Robbins (1956) | normal location model | the identity used is proved directly in Lemma 3.3 |
| Reversibility and non-negative spectrum of the two-block marginal chain | Lemma 3.10, Theorem 4.1(3) | Liu, Wong and Kong (1994) (context) | two-component Gibbs sampler | proved directly in Lemma 3.10 ($`K=TT^*`$) |

**Results proved in the text rather than cited:** the quantile Lipschitz bound (Lemma 3.4); second moments about the mode and the covariance bound (Lemma 3.7); the Cramér–Rao lower bound (Lemma 3.2); the Bregman form of the KL divergence and its bound (Lemma 3.6(i)); the Gaussian total-variation formula and convexity of total variation (Lemma 3.9); the gradient certificate (Lemma 5.4); the mean–mode bounds (Propositions 5.6 and 5.6′); and the ridge bounds (Propositions 6.2–6.3).

**Items still to verify against the original sources:**

1. Ollivier (2009): the proposition numbering in the journal version (only the arXiv v4 numbering was checked), and the journal volume and pages (only the DOI was checked).
2. Brascamp and Lieb (1976), Theorem 4.3: checked only through Saumard and Wellner's attribution; the dimension labelling in Saumard and Wellner's Theorem 3.8 should be compared with the original.
3. Simon (2011), Theorem 13.13: not checked.
4. Villani (2009), Theorem 4.1 and Corollary 5.22: checked only through statements quoted in other papers, not against the book (Appendix B.3).

## Appendix B. Statements of the cited results

The results the proofs rely on are restated below in full, with the original numbering, so that the proofs can be followed without the source papers. Each is restated in the source's own notation where that differs from ours, followed by how it is used. Each entry also records how closely the statement was checked.

### B.1 Saumard and Wellner (2014): strongly log-concave densities and their marginals

**Definition 2.9 ($`\mathrm{SLC}_2`$).** Let $`\mu\in\mathbb R^d`$ and let $`\Sigma`$ be a $`d\times d`$ positive definite matrix. A density $`p`$ on $`\mathbb R^d`$ belongs to $`\mathrm{SLC}_2(\mu,\Sigma,d)`$ if
```math
p(x)=g(x)\,\varphi_\Sigma(x-\mu)\qquad\text{for all }x\in\mathbb R^d,
```
for some log-concave function $`g`$, where $`\varphi_\Sigma`$ is the $`\mathcal N_d(0,\Sigma)`$ density. Equivalently, $`-\log p(x)-\tfrac12(x-\mu)^\top\Sigma^{-1}(x-\mu)`$ is convex. With $`\Sigma=\sigma^2I`$ this is the class $`\mathrm{SLC}_1(\sigma^2,d)`$ of their Definition 2.8.

**Theorem 3.8 (preservation of strong log-concavity under marginalization).** Let $`p\in\mathrm{SLC}_2(0,\Sigma,m+n)`$ be a density on $`\mathbb R^{m+n}`$, with $`x\in\mathbb R^m`$, $`y\in\mathbb R^n`$, and partition
```math
\Sigma=\begin{pmatrix}\Sigma_{11}&\Sigma_{12}\\\Sigma_{21}&\Sigma_{22}\end{pmatrix},\qquad\Sigma_{11}\in\mathbb R^{m\times m}.
```
Then the marginal $`q(x)=\int_{\mathbb R^n}p(x,y)\,dy`$ belongs to $`\mathrm{SLC}_2(0,\Sigma_{11},m)`$. No smoothness is assumed. Saumard and Wellner attribute the result to Brascamp and Lieb (1976), Theorem 4.3, p. 380, and to Simon (2011), Theorem 13.13.

*Use.* Lemma 3.5, with $`\Sigma=K^{-1}`$, after a translation (the class is translation-covariant) and the linear change of variables in the proof.

*Checked:* against the arXiv text of the review (arXiv:1404.5886). The review's statement labels the marginal's dimension inconsistently; the version above, with $`\Sigma_{11}`$ the block of the retained coordinates, is the one its proof supports.

### B.2 Bobkov and Ledoux (2000): transportation inequality

**Setting.** $`(E,\|\cdot\|)`$ is a finite-dimensional normed space with a fixed Lebesgue measure $`dx`$. For $`p\ge1`$ and Borel probability measures $`\mu,\nu`$ on $`E`$,
```math
W_p(\mu,\nu)=\inf_\pi\Bigl(\iint_{E\times E}\|x-y\|^p\,d\pi(x,y)\Bigr)^{1/p},
```
the infimum over all probability measures $`\pi`$ on $`E\times E`$ with marginals $`\mu`$ and $`\nu`$. For $`\nu\ll\mu`$, the relative entropy is $`D(\nu\|\mu)=\int\frac{d\nu}{d\mu}\log\frac{d\nu}{d\mu}\,d\mu`$.

**Condition (3.6).** $`V`$ is convex on $`E`$ and, for some $`p\ge2`$ and $`c>0`$,
```math
tV(x)+sV(y)-V(tx+sy)\ \ge\ \frac cp\,\bigl(s+o(s)\bigr)\,\|x-y\|^p
```
for all $`t,s>0`$ with $`t+s=1`$ and all $`x,y\in E`$, where $`o(s)/s\to0`$ as $`s\to0`$.

**Proposition 4.2.** Let $`\mu`$ be the probability measure $`d\mu=e^{-V}dx`$ on $`E`$, with $`V`$ convex and satisfying (3.6). Then for every probability measure $`\nu`$ on $`E`$ that is absolutely continuous with respect to $`\mu`$ and has a finite moment of order $`p`$,
```math
W_p^{\,p}(\mu,\nu)\ \le\ \frac pc\,D(\nu\|\mu).
```
No $`C^2`$ assumption is made.

*Use.* Lemma 3.6(ii), with $`E=\mathbb R^d`$, $`\|\cdot\|=\|\cdot\|_K`$, $`p=2`$ and $`c=1`$, which gives $`W_{2,K}^2(\mu,\nu)\le2D(\nu\|\mu)`$. Condition (3.6) is verified in that proof, with $`o(s)=-s^2`$.

*Checked:* against the authors' posted PDF of the published article.

### B.3 Villani (2009): optimal couplings

**Theorem 4.1 (existence of an optimal coupling).** Let $`(\mathcal X,\mu)`$ and $`(\mathcal Y,\nu)`$ be Polish probability spaces, and let $`c:\mathcal X\times\mathcal Y\to\mathbb R\cup\{+\infty\}`$ be a lower semicontinuous cost function such that $`c(x,y)\ge a(x)+b(y)`$ for some upper semicontinuous $`a\in L^1(\mu)`$ and $`b\in L^1(\nu)`$. Then there is a coupling $`\pi`$ of $`(\mu,\nu)`$ that minimizes $`\int c\,d\pi`$ among all couplings. In particular this holds when $`c`$ is lower semicontinuous and bounded below.

*Use.* Lemma 3.6(iii), with $`c(x,y)=\|x-y\|^2_K\ge0`$, continuous on $`\mathbb R^q\times\mathbb R^q`$.

**Corollary 5.22 (measurable selection of optimal plans).** Let $`\mathcal X,\mathcal Y`$ be Polish spaces and $`c:\mathcal X\times\mathcal Y\to\mathbb R`$ a continuous cost function with $`\inf c>-\infty`$. Let $`\Omega`$ be a measurable space and $`\omega\mapsto(\mu_\omega,\nu_\omega)`$ a measurable map $`\Omega\to\mathcal P(\mathcal X)\times\mathcal P(\mathcal Y)`$, where $`\mathcal P(\cdot)`$ carries the Borel $`\sigma`$-algebra of weak convergence. Then there is a measurable map $`\omega\mapsto\pi_\omega`$ such that, for each $`\omega`$, $`\pi_\omega`$ is an optimal transference plan between $`\mu_\omega`$ and $`\nu_\omega`$.

*Use.* Lemma 3.6(iii) and Theorem 4.1(1), with $`\Omega=\mathbb R^p\times\mathbb R^p\ni(\gamma,\tilde\gamma)`$ and $`(\mu_\omega,\nu_\omega)=(f_i(\cdot\mid\gamma),f_i(\cdot\mid\tilde\gamma))`$. This map is continuous, even in total variation, by dominated convergence.

*Checked:* only through statements quoted in several recent papers (for example arXiv:2311.05672, Propositions 2.4–2.5, and arXiv:2402.04159, Lemma 4.20), **not against the book itself**. The general lower bound $`c\ge a+b`$ in Theorem 4.1 is recalled from the book and should be checked. Only the special case "$`c`$ bounded below", which is the one quoted in those papers, is used here.

### B.4 Ollivier (2009): coarse Ricci curvature and the spectral gap

Numbering as in arXiv:math/0701886v4 (30 July 2007).

**Standing setting (Definitions 1–2).** $`(X,d)`$ is a Polish metric space with its Borel $`\sigma`$-algebra. A random walk $`m`$ on $`X`$ is a family $`(m_x)_{x\in X}`$ of probability measures on $`X`$ such that $`x\mapsto m_x`$ is measurable, and each $`m_x`$ has finite first moment: $`\int d(o,y)\,dm_x(y)<\infty`$ for some (hence every) $`o\in X`$. $`T_1`$ denotes the $`L^1`$ transportation distance,
```math
T_1(\mu,\nu)=\inf_\xi\iint d(x,y)\,d\xi(x,y),
```
over couplings $`\xi`$ of $`\mu`$ and $`\nu`$; this is $`W_1`$ in our notation. The averaging operator is $`Mf(x)=\int f\,dm_x`$.

**Definition 3 (coarse Ricci curvature).** For distinct $`x,y\in X`$,
```math
\kappa(x,y)=1-\frac{T_1(m_x,m_y)}{d(x,y)}.
```
The Ricci curvature of $`X`$ is at least $`\kappa`$ if $`\kappa(x,y)\ge\kappa`$ for all distinct $`x,y`$.

**Definition 18 (spread).** The spread of the walk at $`x`$ is
```math
\sigma(x)=\Bigl(\tfrac12\iint d(y,z)^2\,dm_x(y)\,dm_x(z)\Bigr)^{1/2}.
```
If $`\nu`$ is an invariant distribution, the average spread is $`\sigma=\|\sigma(\cdot)\|_{L^2(X,\nu)}`$.

**Proposition 29 (spectral radius).** Let $`(X,d,m)`$ be a metric space with random walk and invariant distribution $`\nu`$. Suppose that the Ricci curvature of $`X`$ is at least $`\kappa>0`$ and that $`\sigma<\infty`$. Suppose also that $`\nu`$ is reversible, or that $`X`$ is finite. Then the spectral radius of the averaging operator $`M`$ acting on $`L^2(X,\nu)/\{\mathrm{const}\}`$ is at most $`1-\kappa`$.

(Corollary 30 restates this as a lower bound $`\kappa`$ on the spectral gap of $`-\Delta=\mathrm{Id}-M`$, for an ergodic random walk with reversible $`\nu`$.)

*Use.* Lemma 3.10, with $`X=\mathbb R^p`$, $`d(x,y)=\|x-y\|_{P_0}`$, $`m_x=K(x,\cdot)`$, $`\nu=\pi_\gamma`$ and $`\kappa=1-\rho`$. The hypotheses are checked in that proof.

*Checked:* against the arXiv v4 text. The journal version's numbering was not checked.

### B.5 Brascamp and Lieb (1976), Theorem 4.1, for reference only

**Brascamp–Lieb variance inequality** (as stated in Saumard and Wellner 2014, Proposition 10.1(a)). Let $`X`$ have density $`e^{-\phi}`$ on $`\mathbb R^d`$, with $`\phi\in C^2(\mathbb R^d)`$ strictly convex, and let $`G\in C^1(\mathbb R^d)\cap L^2`$. Then
```math
\mathrm{Var}\bigl(G(X)\bigr)\le E\bigl[\nabla G(X)^\top(\nabla^2\phi(X))^{-1}\nabla G(X)\bigr].
```

*Not used.* Lemma 3.7(iii) proves the case $`G(x)=a^\top x`$ that the paper needs, with $`\nabla^2\phi\succeq K`$ replaced by $`K`$-strong convexity and no $`C^2`$ assumption.

### B.6 Tweedie's formula, for attribution only

If $`z\mid\mu\sim\mathcal N(\mu,\sigma^2)`$ and $`\mu`$ has prior $`g`$, with marginal density $`f`$ of $`z`$, then $`E[\mu\mid z]=z+\sigma^2\,\tfrac{d}{dz}\log f(z)`$ (Efron 2011, who attributes it to Robbins 1956).

*Use.* Lemma 3.3 proves the multivariate analogue, $`\nabla\log m_i=A_i^\top(E[\beta_i\mid\gamma]-H_i\gamma)`$, directly by differentiation under the integral. The formula is cited only for its name.

## References

- Bobkov, S. G. and Ledoux, M. (2000). From Brunn–Minkowski to Brascamp–Lieb and to logarithmic Sobolev inequalities. *Geometric and Functional Analysis* 10(5), 1028–1052. doi:10.1007/PL00001645.
- Brascamp, H. J. and Lieb, E. H. (1976). On extensions of the Brunn–Minkowski and Prékopa–Leindler theorems, including inequalities for log concave functions, and with an application to the diffusion equation. *Journal of Functional Analysis* 22(4), 366–389. doi:10.1016/0022-1236(76)90004-5.
- Efron, B. (2011). Tweedie's formula and selection bias. *Journal of the American Statistical Association* 106(496), 1602–1614. doi:10.1198/jasa.2011.tm11181.
- Liu, J. S., Wong, W. H. and Kong, A. (1994). Covariance structure of the Gibbs sampler with applications to the comparisons of estimators and augmentation schemes. *Biometrika* 81(1), 27–40. doi:10.1093/biomet/81.1.27.
- Madras, N. and Sezer, D. (2010). Quantitative bounds for Markov chain convergence: Wasserstein and total variation distances. *Bernoulli* 16(3), 882–908. doi:10.3150/09-BEJ238.
- Nygren, K. and Nygren, L. M. (2006). Likelihood subgradient densities. *Journal of the American Statistical Association* 101(475), 1144–1156. doi:10.1198/016214506000000357.
- Nygren, K. (2020). On the total variation distance between multivariate normal densities with applications to two-block Gibbs samplers. Manuscript, version dated 20 July 2020.
- Ollivier, Y. (2009). Ricci curvature of Markov chains on metric spaces. *Journal of Functional Analysis*. doi:10.1016/j.jfa.2008.11.001. Cited by the numbering of arXiv:math/0701886v4 (30 July 2007).
- Robbins, H. (1956). An empirical Bayes approach to statistics. *Proceedings of the Third Berkeley Symposium on Mathematical Statistics and Probability*, Vol. 1. *(Cited only through Efron 2011; details to be checked.)*
- Roberts, G. O. and Rosenthal, J. S. (2002). One-shot coupling for certain stochastic recursive sequences. *Stochastic Processes and their Applications* 99(2), 195–208. doi:10.1016/S0304-4149(02)00096-0.
- Rosenthal, J. S. (1995). Minorization conditions and convergence rates for Markov chain Monte Carlo. *Journal of the American Statistical Association* 90, 558–566.
- Saumard, A. and Wellner, J. A. (2014). Log-concavity and strong log-concavity: a review. *Statistics Surveys* 8, 45–114. doi:10.1214/14-SS107.
- Simon, B. (2011). *Convexity: An Analytic Viewpoint*. Cambridge Tracts in Mathematics 187. Cambridge University Press.
- Villani, C. (2009). *Optimal Transport: Old and New*. Grundlehren der mathematischen Wissenschaften 338. Springer, Berlin.
