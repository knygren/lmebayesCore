# Restricted Two-Block Gibbs: Minorization on $C_d$ and Distance to $\pi$

*Draft. This pass carries the full outline, the definitions needed to state the
main results, and the statements of Theorems 1 and 2. Supporting lemmas, closed
forms, and all proofs are left as stubs in §4–§7 and the Appendix, to be filled in
next.*

---

## 1. Introduction

### 1.1 Problem and setting

We study a two-block Gibbs sampler for the posterior of a Bayesian hierarchical
generalized linear mixed model (GLMM), and ask how quickly its distribution after
$n$ sweeps approaches the target $\pi$ in total variation. The sampler alternates
an exact draw of the random effects $\beta=(\beta_j)_{j\le J}$ given the population
parameter $\gamma$, and an exact (or near-exact) draw of $\gamma$ given $\beta$. The
$\gamma$-marginal of this chain is itself Markov, and it is on this marginal chain
that we establish geometric ergodicity.

The argument has two independent parts. Part I (§2) is a general fact about any
Markov chain that satisfies a Doeblin minorization condition on a subset $C_d$ of
its state space: the chain, run on $C_d$, converges geometrically to the
restriction of $\pi$ to $C_d$, and this transfers to a bound on distance from the
*full* target $\pi$ at the cost of the mass $\pi(C_d^c)$ left outside. Part II
(§3–§6) is specific to the hierarchical GLMM: it constructs a set $C_d$ and a
minorization constant explicitly, for every choice of tail-mass budget $\delta$,
under hypotheses satisfied by the package's supported GLM families.

### 1.2 Statement of the main results (informal)

- **Theorem 1** (§2.4): if a Markov kernel minorizes on $C_d$ with constant
  $\varepsilon_d$ and $\pi(\cdot\mid C_d)$ is stationary for the restricted chain,
  then the restricted chain converges to $\pi$ at rate $(1-\varepsilon_d)^n$, up to
  an additive error $\pi(C_d^c)$.
- **Theorem 2** (§3.4): for the $\gamma$-marginal of the two-block GLMM sampler,
  under hypotheses (H1)–(H3), such a $C_d$ and $\varepsilon_d$ exist explicitly for
  every tail-mass budget $\delta>0$ — with no bound on the minimum eigenvalue of the
  precision required, and with a nontrivial limit as the population prior becomes
  flat.

### 1.3 Relation to Rosenthal (1995)

Theorem 1 is an application of Rosenthal (1995), Proposition 2 (uniform ergodicity
from a Doeblin condition on the whole restricted state space); see `minor.pdf`,
p. 5. What is specific to this note is Theorem 2: an explicit, checkable
construction of the minorizing set and constant for a hierarchical GLMM posterior,
rather than an existence argument alone.

---

## 2. The Conditional Bound

### 2.1 Two-block Gibbs sampler and the restricted kernel

**Definition 1 (state space, restricted kernel).** Let $(\gamma,\beta)\in\mathcal X$
and let $\pi$ be the two-block Gibbs target on $\mathcal X$. Fix $C_d\subseteq
\mathcal X$ with $\pi(C_d)>0$. Write

$$
\pi(A\mid C_d):=\frac{\pi(A\cap C_d)}{\pi(C_d)},\qquad A\in\mathcal B(\mathcal X).
$$

Let $P_1(x,\cdot)$ be one Gibbs sweep on $\mathcal X$, and let $P_1(x,\cdot\mid C_d)$
be the kernel of the sampler **restricted to** $C_d$ — defined once, per
application, so that $P_1(x,C_d\mid C_d)=1$ for $x\in C_d$. Write
$P_n(x,\cdot\mid C_d):=\bigl(P_1(x,\cdot\mid C_d)\bigr)^n$.

**Definition 2 (Doeblin minorization on $C_d$).** There exist $\varepsilon_d\in(0,1]$
and a probability measure $Q$ on $C_d$ such that

$$
P_1(x,A\mid C_d)\;\ge\;\varepsilon_d\,Q(A)\qquad\forall\,x\in C_d,\ \forall\,\text{measurable }A\subseteq C_d.
$$

**Standing assumption.** $\pi(\cdot\mid C_d)$ is a stationary law for
$P_1(\cdot\mid C_d)$ on $C_d$. (Must be verified for the chosen restriction; not
automatic from stationarity of $\pi$ on $\mathcal X$.)

### 2.2 Total variation and L1 distance

**Definition 3.**

$$
\|\mu-\nu\|_{TV}:=\sup_{A\subseteq\mathcal X}\bigl|\mu(A)-\nu(A)\bigr|=\tfrac12\int|d\mu-d\nu|,
\qquad
\|\mu-\nu\|_1:=\int|d\mu-d\nu|=2\|\mu-\nu\|_{TV}.
$$

### 2.3 Truncation: $\pi(\cdot\mid C_d)$ vs. $\pi$

**Lemma 1 (truncation).** With $\mu:=\pi(\cdot\mid C_d)$ and $\nu:=\pi$,

$$
\|\pi(\cdot\mid C_d)-\pi\|_{TV}=\pi(C_d^c),
\qquad
\|\pi(\cdot\mid C_d)-\pi\|_1=2\,\pi(C_d^c).
$$

### 2.4 Theorem 1 — geometric convergence of the restricted chain

**Theorem 1.** Under Definitions 1–2 and the standing stationarity assumption, for
$x\in C_d$ and $n\ge1$,

$$
\boxed{\;\bigl\|P_n(x,\cdot\mid C_d)-\pi\bigr\|_{TV}\;\le\;(1-\varepsilon_d)^n+\pi(C_d^c).\;}
$$

Equivalently, $\bigl\|P_n(x,\cdot\mid C_d)-\pi\bigr\|_1\le2(1-\varepsilon_d)^n+2\,\pi(C_d^c)$.

*Proof outline (full proof in Appendix, §A.0 — carried over unchanged from the
working draft, not reproduced here).* Triangle inequality splits the bound into
(i) distance from the restricted chain to $\pi(\cdot\mid C_d)$, controlled by
Rosenthal's Proposition 2 applied on $C_d$ with kernel $P_1(\cdot\mid C_d)$; and
(ii) distance from $\pi(\cdot\mid C_d)$ to $\pi$, which is Lemma 1. $\qquad\blacksquare$

**Remark 1.1.** Nothing here is specific to the joint $(\gamma,\beta)$ state space;
Theorem 1 holds for any $(\mathcal X,\pi,C_d,P_1)$ satisfying Definitions 1–2. §3–§6
apply it with $\mathcal X=\mathbb R^q$, $\pi=\pi_\gamma$ (the $\gamma$-marginal
posterior), and $P_1=q(\cdot\mid\cdot\,;C)$.

**Remark 1.2 (Rosenthal / Meyn–Tweedie).** Theorem 1 uses **minorization only**
(Rosenthal Proposition 2 on $C_d$) plus **truncation** (Lemma 1). Full geometric
ergodicity with **both** minorization and Foster drift — Rosenthal (1995) Theorem 12
and the Meyn–Tweedie rate $\|P^k-\pi\|_{TV}\le M V(x)\rho^k$ — is documented in
`inst/GEOMETRIC_ERGODICITY_MT_ROSENTHAL.md`.

---

## 3. The Hierarchical GLMM Model

### 3.1 Model and hypotheses

Work in the two-block hierarchy of `notation.md` / `model_setup()`:

$$
y_j\mid\beta_j\ \text{GLM with canonical/noncanonical link},\qquad
\beta_j\mid\gamma\sim N(H_j\gamma,\,\Psi),\qquad
\gamma\sim N(\mu_0,\,\Lambda_\gamma^{-1}),
$$

with $P_b:=\Psi^{-1}$, $P_{11}^{\mathrm{RE}}:=\sum_j H_j^\top P_b H_j$, and
$P_{11}:=\Lambda_\gamma+P_{11}^{\mathrm{RE}}$ the exact conditional precision of
$\gamma$ given $\beta$.

- **(H1) Proper posterior.** $\pi(\gamma,\beta\mid y)$ is a probability measure on
  $\mathcal X=\mathbb R^q\times\prod_j\mathbb R^{p_{\mathrm{re},j}}$. Automatic when
  $\Lambda_\gamma\succ0$ and the likelihood is bounded; at $\Lambda_\gamma=0$ this is
  **derived** in §3.2 from (H2)+(H3).
- **(H2) Log-concave conditionals.** For each block draw, $\pi(\beta_j\mid\gamma,y_j)$
  is log-concave in $\beta_j$ (logit, probit, Poisson, Gamma-log, etc.).
- **(H3) Full rank / estimability.**
  - **(H3a)** Hyper-design full rank: $P_{11}^{\mathrm{RE}}\succ0$.
  - **(H3b)** Group-wise estimability: for each $j$, a finite MLE exists for the
    group GLM on $\beta_j$ (no complete separation, no all-zero Poisson intercept
    path, etc.).

(H1)–(H3) are the preflight gates under which the package runs production samplers.
Unlike the working draft, no hypothesis on the one-step kernel's positivity on
compacts (formerly (H4)) is carried as a standing assumption — the construction in
§4–§5 writes the minorizing density down explicitly, so that property is a
consequence, not an input.

### 3.2 Propriety of the flat-prior posterior

This subsection supplies **(H1) at $\Lambda_\gamma=0$**, needed because the flat
limit is not covered by the "bounded likelihood" argument available when
$\Lambda_\gamma\succ0$.

**Lemma 2 (propriety of the flat-$\gamma$ posterior).** Assume:

- **(P1)** each $\ell_j$ (group log-likelihood) is concave in $\beta_j$ with strictly
  positive IRLS weights at finite linear predictor — (H2);
- **(P2)** each group design $Z_j$ has full column rank;
- **(P3)** the stacked hyper-design has rank $q$, i.e. $P_{11}^{\mathrm{RE}}\succ0$ —
  (H3a);
- **(P4)** a finite group-wise MLE $\hat\beta_j$ exists for every $j$ — (H3b).

Then the flat-$\gamma$ posterior

$$
\pi_0(\gamma,\beta\mid y)\propto\exp\bigl(g(\gamma,\beta)\bigr),\qquad
g(\gamma,\beta):=\sum_j\ell_j(\beta_j)-\tfrac12\sum_j(\beta_j-H_j\gamma)^\top P_b(\beta_j-H_j\gamma),
$$

is **proper**: $\int_{\mathbb R^q}\int_{\mathbb R^{Jp_{\mathrm{re}}}}e^g\,d\beta\,d\gamma<\infty$.
Hence (H1) holds at $\Lambda_\gamma=0$.

*(Proof in Appendix A.1: $g$ is concave, bounded above by (P1)+(P4), and recedes to
$-\infty$ in every direction of $(\gamma,\beta)$-space — quadratically off the
$\gamma$-flat directions by strong concavity of the RE prior, and at least linearly
along them by (P2)+(P3) forcing every such direction to move some $\beta_j$ off its
finite maximum. (P4) is sufficient but not necessary; the sharp condition is stated
alongside the proof.)*

### 3.3 The mean map and the marginal $\gamma$-chain

In two-block Gibbs, the $\gamma$-marginal is itself Markov, with kernel

$$
q(\gamma'\mid\gamma)=\int\phi_q\bigl(\gamma';\,m(\beta),\,P_{11}^{-1}\bigr)\,\pi(\beta\mid\gamma,y)\,d\beta,
\qquad
m(\beta)=P_{11}^{-1}\Bigl(\Lambda_\gamma\mu_0+\sum_jH_j^\top P_b\,\beta_j\Bigr).
$$

**Definition 4 (mean map).** $M(\gamma):=P_{11}^{-1}\bigl(\Lambda_\gamma\mu_0+\sum_jH_j^\top P_b\,b_j(\gamma)\bigr)$,
where $b_j(\gamma):=E[\beta_j\mid\gamma,y]$.

**Lemma 3 (the mean map has a unique fixed point).** Assume (H2), (H3a), and either
$\Lambda_\gamma\succ0$ or the strict-curvature condition that the group likelihoods
contribute curvature ($-\nabla^2\ell_j\succeq\epsilon_jP_b$ uniformly, some
$\epsilon_j>0$). Then $M$ is a contraction in the norm $\|x\|_{P_{11}}=\|P_{11}^{1/2}x\|$
with modulus $\kappa<1$; it has a unique fixed point $\gamma^\star$, and the
iteration $\gamma^{(k+1)}=M(\gamma^{(k)})$ converges to it geometrically at rate
$\kappa$ from any start.

*(Proof in Appendix A.1: the Jacobian of $M$ is $J(\gamma)=P_{11}^{-1}\sum_jH_j^\top
P_bV_j(\gamma)P_bH_j$ with $V_j(\gamma)=\mathrm{Cov}(\beta_j\mid\gamma,y)$; a
Brascamp–Lieb bound under (H2) gives $V_j\preceq P_b^{-1}$, hence $\kappa<1$ either
from $\Lambda_\gamma\succ0$ directly or from the likelihood curvature margin in the
flat limit; Banach's fixed-point theorem finishes it.)*

$\gamma^\star$ is, equivalently, the mode of the marginal posterior $\pi(\gamma\mid
y)$ — this identification, and the objects it feeds into ($Q$ in §4, the
minorization profile in §5), are stated where first needed rather than collected
here.

### 3.4 Theorem 2 — existence of a certified safe set (statement)

**Theorem 2.** Under (H1), (H2), (H3a) — and (H3b) as well if $\Lambda_\gamma=0$ —
for every $\delta>0$ there exist:

1. a compact, convex set $C=\widetilde C_d\subseteq\mathbb R^q$, a superlevel set of
   a minorization profile $\varepsilon(\cdot)$ constructed in §5, with $d=d(\delta)$
   chosen so that $\pi_\gamma(C^c)<\delta$;
2. a restricted $\gamma$-chain kernel $q(\cdot\mid\cdot\,;C)$ — the two-block sweep
   confined to $C$;
3. a constant $\varepsilon=\varepsilon_d\,Q(C)\in(0,1]$ and the truncated Gaussian
   $Q_C=Q(\cdot\mid C)$, where $Q=N(\gamma^\star,\Sigma^\star)$, $\Sigma^\star=
   P_{11}^{-1}$, is the refresh measure of §4,

such that $\pi_\gamma(\cdot\mid C)$ is stationary for $q(\cdot\mid\cdot\,;C)$ and
$q(\gamma,A\mid C)\ge\varepsilon\,Q_C(A)$ for all $\gamma\in C$, $A\subseteq C$ —
hence, by Theorem 1, for $\gamma\in C$,

$$
\|q_n(\gamma,\cdot\mid C)-\pi_\gamma\|_{TV}\;\le\;(1-\varepsilon)^n+\delta.
$$

The refresh measure $Q$ is fixed **before** any set is chosen (§4); the certified
set is then the superlevel set of the resulting minorization profile (§5), so its
compactness, its exhaustion of $\mathbb R^q$, and the value of $\varepsilon_d$ are
all consequences of the same construction rather than separate searches. §6 proves
Theorem 2 by assembling §4–§5 and lifting the bound from the $\gamma$-marginal to
the joint $(\gamma,\beta)$ chain.

---

## 4. The Refresh Measure $Q$

The refresh measure is fixed **before** any set is chosen — $Q$ depends only on the
model, not on a tail-mass budget $\delta$ or a truncation level $d$. Write
$P_{11}=\Lambda_\gamma+P_{11}^{\mathrm{RE}}$, $P_{11}^{\mathrm{RE}}:=\sum_jH_j^\top
P_bH_j$, $P_b:=\Psi^{-1}$, as in §3.1. Two conditions pin $Q$ down.

**(Q1) Precision dominance.** $Q=N(\gamma^\star,\Sigma_Q)$ must satisfy
$\Lambda_Q:=\Sigma_Q^{-1}\succeq\Lambda_q(\gamma'\mid\gamma):=-\nabla^2_{\gamma'}
\log q(\gamma'\mid\gamma)$ for all $\gamma,\gamma'$ — the refresh measure must be at
least as concentrated as the transition kernel, or the log-ratio $\log q/q_Q$ is
unbounded below and no minorization constant exists. The choice $\Sigma_Q:=P_{11}^{-1}$
achieves this: since $q(\cdot\mid\gamma)$ is a mixture of Gaussians with common
covariance $P_{11}^{-1}$ and random mean $m(\beta)$, differentiating through the
mixture gives, for every GLM family,
$$
\Lambda_Q-\Lambda_q(\gamma'\mid\gamma)=\Lambda_Q\,\mathrm{Cov}\bigl(m(\beta)\mid\gamma,\gamma',y\bigr)\,\Lambda_Q\;\succeq\;0,
$$
the covariance taken over the bridge law $\pi(\beta\mid\gamma,\gamma',y)$ — an
*identity*, not a hypothesis, so (Q1) never has to be checked model by model.
(Minorizing the exact one-step Gaussian draw $\gamma'\mid\beta$ directly, without
integrating $\beta$ out, does not work: with $\beta$ held fixed the bridge
covariance above vanishes, the precision gap closes exactly, and the log-ratio
against any centered Gaussian is affine with infimum $-\infty$. Integrating $\beta$
out is what produces a strict gap.)

**(Q2) Mean.** $\gamma^\star=M(\gamma^\star)$, the fixed point of the mean map
$M$ of Definition 4 — equivalently (Cor. R2d-5, §3.3) the mode of the
$\gamma$-marginal posterior $\pi(\gamma\mid y)$.

Existence of $\Sigma^\star:=P_{11}^{-1}$ and $\gamma^\star$ — i.e. of $Q$ itself — is
what §4.1 and §4.2 establish, for a proper population prior and in the flat-prior
limit respectively.

### 4.1 Proper population prior ($\Lambda_\gamma\succ0$)

**Lemma 4 (the covariance of $Q$ exists).** Assume $\Lambda_\gamma\succ0$. Then
$P_{11}\succeq\Lambda_\gamma\succ0$, so
$$
\Sigma^\star:=P_{11}^{-1}\ \text{exists and}\ 0\prec\Sigma^\star\preceq\Lambda_\gamma^{-1}.
$$
No rank condition on the hyper-design is needed: a proper population prior supplies
invertibility on its own.

**Lemma 5 (the centre of $Q$ exists).** Assume (H2) and $\Lambda_\gamma\succ0$. The
negative log marginal posterior $\Phi(\gamma):=-\log\pi(\gamma\mid y)$ is
$\Lambda_\gamma$-strongly convex,
$$
\nabla^2\Phi(\gamma)\;=\;P_{11}-\sum_jH_j^\top P_bV_j(\gamma)P_bH_j\;\succeq\;\Lambda_\gamma\;\succ\;0
\qquad\text{for every }\gamma,
$$
where $V_j(\gamma):=\mathrm{Cov}(\beta_j\mid\gamma,y)$. Consequently $\Phi$ has a
unique, finite minimiser $\gamma^\star=\arg\min_\gamma\Phi(\gamma)=M(\gamma^\star)$ —
the mode of the $\gamma$-marginal posterior exists, is unique, is finite, and is the
fixed point of $M$. (Like Lemma 4, this uses no rank condition — $\Lambda_\gamma\succ0$
alone is enough. Uniqueness and existence can also be reached via Lemma 3, §3.3,
by Banach's fixed-point theorem; the two routes agree, but the variational route
here is the one that survives the flat-prior limit of §4.2, where Lemma 3's
contraction modulus $\kappa$ may equal $1$.)

**Corollary 6 (the refresh density $Q$ exists).** Under (H2) and $\Lambda_\gamma\succ0$,
$$
Q:=N(\gamma^\star,\Sigma^\star),\qquad
q_Q(\gamma')=(2\pi)^{-q/2}\det(P_{11})^{1/2}\exp\Bigl\{-\tfrac12(\gamma'-\gamma^\star)^\top P_{11}(\gamma'-\gamma^\star)\Bigr\},
$$
is a nondegenerate multivariate normal probability density, strictly positive and
continuous on all of $\mathbb R^q$, satisfying (Q1) (by the precision-dominance
identity above, which holds for every GLM family) and (Q2) (by Lemma 5).

### 4.2 Flat-prior limit ($\Lambda_\gamma\downarrow0$)

Write $\Lambda_\gamma=\tau\bar\Lambda$ for fixed $\bar\Lambda\succeq0$ and $\mu_0$,
and let $\tau\downarrow0$; write $P_{11}(\tau),\Sigma^\star(\tau),\gamma^\star(\tau),
Q_\tau$ for the corresponding objects. Both results of §4.1 degrade at $\tau=0$: the
prior no longer supplies invertibility (needed for Lemma 4) or strong convexity
(needed for Lemma 5). (H3a) replaces the first, (H3b) the second.

**Lemma 7 (the limiting covariance exists iff (H3a)).** As $\tau\downarrow0$,
$P_{11}(\tau)\downarrow P_{11}^{\mathrm{RE}}$ monotonically in the psd order. Under
(H3a) ($P_{11}^{\mathrm{RE}}\succ0$),
$$
\Sigma^\star(\tau)\ \uparrow\ \Sigma^{\star(0)}:=\bigl(P_{11}^{\mathrm{RE}}\bigr)^{-1}\;\prec\;\infty,
$$
monotonically in the psd order and entrywise. Without (H3a) the limit does not
exist: for $0\ne v\in\ker P_{11}^{\mathrm{RE}}$, $v^\top\Sigma^\star(\tau)v\to\infty$.
($\Sigma^{\star(0)}$ is the generalized-least-squares covariance for the
hyper-regression of the random effects on the hyper-design with weight $\Psi^{-1}$;
the prior term $\Lambda_\gamma\mu_0$ drops out of $M$ entirely in the limit, and
$\Sigma^\star(\tau)$ increases to $\Sigma^{\star(0)}$ — the *widest* admissible
refresh covariance.)

**Lemma 8 (the limiting centre is a finite vector, under (H3a)+(H3b)).** Assume
(H2), (H3a), (H3b). Let
$$
\Phi_0(\gamma)=\sum_j\mathcal D_j\bigl(P_b H_j\gamma\bigr)+\text{const}
$$
be the flat-limit negative log marginal posterior, where the group information
deficit $\mathcal D_j(\theta):=\tfrac12\theta^\top P_b^{-1}\theta+a_j-A_j(\theta)\ge0$
is convex (Brascamp–Lieb) with $\nabla^2\mathcal D_j=P_b^{-1}-V_j\succeq0$. Then
$\Phi_0$ is convex and coercive — coercivity following from propriety of the
flat-limit joint posterior (Lemma 2, §3.2) — so
$$
\gamma^{\star(0)}:=\arg\min_\gamma\Phi_0(\gamma)=\arg\max_\gamma\pi_0(\gamma\mid y)
$$
exists, is finite, and $\gamma^\star(\tau)\to\gamma^{\star(0)}$ as $\tau\downarrow0$.

*Open item carried over from removing the spectral-cutoff appendix:* the source
argument for **uniqueness** of $\gamma^{\star(0)}$ went through a strict-positivity
condition on the group Fisher information stated in that appendix (now dropped).
Convexity, coercivity, existence, and the convergence $\gamma^\star(\tau)\to
\gamma^{\star(0)}$ do not need it; uniqueness needs a replacement condition — most
likely $P_{11}^{\mathrm{RE}}\succ0$ itself (already assumed as (H3a)) together with
strict convexity of at least one $\mathcal D_j$ along every direction, but this needs
to be re-derived cleanly before Appendix A.2 is written, not carried over from the
old citation.

**Corollary 9 (the limiting refresh density exists).** Under (H2), (H3a), (H3b),
$$
Q_0:=N\bigl(\gamma^{\star(0)},\Sigma^{\star(0)}\bigr)
$$
is a nondegenerate multivariate normal, and $Q_\tau\to Q_0$ in total variation
(equivalently $L^1$, by Scheffé) as $\tau\downarrow0$.

**Remark 4.1 (the two hypotheses are sharp, and fail differently).**

| Failure | What breaks | Which result dies |
|---|---|---|
| (H3a) fails: stacked hyper-design rank-deficient | $P_{11}^{\mathrm{RE}}$ singular $\Rightarrow$ $\Sigma^{\star(0)}$ does not exist | Lemma 7 — and with it $q_{Q_0}$, so §5's flat-limit profile has no density to compare against |
| (H3b) fails: some group not estimable | $\Phi_0$ not coercive $\Rightarrow\ \pi_0$ improper | Lemma 8 — no finite mode, so no (Q2), and no probability measure for §5's exhaustion argument |

Both are visible directly in $\Phi_0$: if $0\ne v\in\bigcap_j\ker H_j$, then
$\Phi_0(\gamma+tv)=\Phi_0(\gamma)$ for all $t$ — the flat-limit marginal posterior is
constant along $v$ and has no mode, with $\tau\bar\Lambda$ having been the only
thing keeping $P_{11}$ invertible. If instead some group is separated, $\ell_j$
fails to decay along a direction, $\mathcal D_j$ stops growing, and $\Phi_0$
flattens at infinity along it. For every $\tau>0$ both failures are repaired by
$\Lambda_\gamma\succ0$ (Lemmas 4–5 used neither hypothesis), at the price of
constants that degrade as $\tau\downarrow0$ — which is exactly why the preflight
gates on hyper-design rank and group estimability are necessary, not merely
prudent, whenever the sampler runs with a flat population prior.

### 4.3 Summary: existence and nondegeneracy of $Q$

$Q=N(\gamma^\star,\Sigma^\star)$ exists and is nondegenerate under (H2)+(H3a), with
$\Lambda_\gamma\succ0$ needed only to reach it by the shorter argument of §4.1; the
same conclusion holds at $\Lambda_\gamma=0$ under the strictly stronger (H2)+(H3a)+(H3b)
via §4.2, and $Q_\tau\to Q_0$ continuously as the prior is relaxed. In particular
**the refresh measure does not degenerate as the population prior becomes diffuse** —
it converges to a specific, finite-variance limit rather than blowing up or ceasing
to exist. Precision dominance (Q1) holds unconditionally, in both regimes, as an
identity rather than a hypothesis. What §4 does *not* yet establish is that $Q$
actually minorizes the sampler kernel on a useful set with a useful constant — that
is §5.

---

## 5. The Minorization Constant and the Certified Set

With $Q$ fixed by §4, define the log-ratio $g(\gamma'\mid\gamma):=\log q(\gamma'\mid\gamma)-\log q_Q(\gamma')$ — finite and jointly continuous, since both densities are strictly positive and continuous.

**Definition 5 (profile, deficiency gap, certified set).**
$$
\varepsilon(\gamma):=\exp\Bigl\{\inf_{\gamma'\in\mathbb R^q}g(\gamma'\mid\gamma)\Bigr\}=\inf_{\gamma'}\frac{q(\gamma'\mid\gamma)}{q_Q(\gamma')},
\qquad
\mathcal D(\gamma):=\log\varepsilon(\gamma).
$$
For $d>0$,
$$
\Psi(\gamma):=\mathcal D(\gamma^\star)-\mathcal D(\gamma)\ge0,
\qquad
\boxed{\ \widetilde C_d:=\{\gamma:\varepsilon(\gamma)\ge\varepsilon_d\}=\{\gamma:\Psi(\gamma)\le d\},\qquad \varepsilon_d:=e^{-d}\varepsilon(\gamma^\star).\ }
$$
The certified set is a superlevel set of the very function being minimised for the constant, which is what makes the eventual minimum of $\varepsilon$ over $\widetilde C_d$ land on a *known* value ($\varepsilon_d$) rather than an unrelated number — that identification is §5.3.

### 5.1 The minorization profile $\varepsilon(\gamma)$

Assume (H2), (H3a), $\Lambda_\gamma\succ0$ throughout §5.1–§5.3; §5.4 relaxes to $\Lambda_\gamma=0$.

**Lemma 10 (uniqueness).** For fixed $\gamma$, $\gamma'\mapsto g(\gamma'\mid\gamma)$ is strictly convex,
$$
\nabla^2_{\gamma'}g(\gamma'\mid\gamma)=\Lambda_Q\,\mathrm{Cov}\bigl(m(\beta)\mid\gamma,\gamma',y\bigr)\,\Lambda_Q\succ0
$$
under (H3a), so it has at most one minimiser.

**Lemma 11 (existence and positivity).** For every $\gamma\in\mathbb R^q$ the infimum defining $\varepsilon(\gamma)$ is attained (at the unique point of Lemma 10), and $\varepsilon(\gamma)\in(0,1]$. Up to an affine change of variable, $g(\cdot\mid\gamma)$ is the cumulant generating function of the law of $m(\beta)$ under $\pi(\beta\mid\gamma,y)$, recentred at $\gamma^\star$; under (H3a) that law has full support $\mathbb R^q$, which forces coercivity along every ray and hence attainment. Coercivity here comes from *support*, not curvature, so this survives likelihoods with unbounded log-density curvature (e.g. Poisson-log).

**Lemma 12 (continuity).** $\varepsilon:\mathbb R^q\to(0,1]$ is continuous. A locally-uniform coercivity estimate, worked in $\beta$-space where the mixing law does not depend on $\Lambda_\gamma$, confines every minimiser over a compact $K\ni\gamma$ to one $K$-dependent (but $\gamma$-free) ball; Berge's maximum theorem then gives continuity on $K$. The same estimate is what lets §5.4 pass this result to $\Lambda_\gamma=0$.

**Remark 5.1 (range).** $\varepsilon(\gamma)\le1$ always: $q(\cdot\mid\gamma)\ge\varepsilon(\gamma)q_Q$ integrates to $1\ge\varepsilon(\gamma)$, with equality iff $q(\cdot\mid\gamma)=q_Q$.

### 5.2 Convexity, coercivity, and compactness of $\widetilde C_d$

**Lemma 13 (the deficiency is concave; closed form for both derivatives).** Under (H2), $\mathcal D\in C^2$ with
$$
\nabla\mathcal D(\gamma)=P_{11}\bigl(\gamma^\star-M(\gamma)\bigr),
\qquad
\nabla^2\mathcal D(\gamma)=-P_{11}J(\gamma)=-\sum_jH_j^\top P_bV_j(\gamma)P_bH_j\;\preceq\;0.
$$
Consequently: **(1)** $\mathcal D$ is concave, strictly so under (H3a), so every superlevel set $\{\varepsilon\ge\varepsilon_*\}$ is convex; **(2)** $\nabla\mathcal D(\gamma)=0\iff M(\gamma)=\gamma^\star\iff\gamma=\gamma^\star$ (Lemma 3, §3.3) — $\varepsilon$ attains its global maximum exactly at the refresh centre; **(3)** comparing with $\nabla\Phi(\gamma)=P_{11}(\gamma-M(\gamma))$ (the Hessian identity behind Lemma 5) gives the exact identity
$$
\boxed{\ \mathcal D(\gamma)-\mathcal D(\gamma^\star)=\bar\Phi(\gamma)-\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}},
\qquad\text{i.e.}\qquad
\Psi(\gamma)+\bar\Phi(\gamma)=\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}},\ }
$$
where $\bar\Phi:=\Phi-\min\Phi$ — the deficiency gap and the recentred negative-log-posterior partition the refresh quadratic exactly.

**Lemma 14 (the deficiency gap is coercive; the certified sets are compact).** Assume (H1), (H2), (H3a). With $S(\gamma):=\nabla^2\Psi(\gamma)=\sum_jH_j^\top P_bV_j(\gamma)P_bH_j$, there is an explicit $c>0$ with
$$
\boxed{\ \Psi(\gamma)\ge c\bigl(|\gamma-\gamma^\star|-1\bigr)\ \text{ for }|\gamma-\gamma^\star|\ge1,
\qquad\text{hence}\qquad
\widetilde C_d\subseteq\{\gamma:|\gamma-\gamma^\star|\le1+d/c\}.\ }
$$
So $\widetilde C_d$ is compact for every $d>0$ — closed (Lemma 12), convex (Lemma 13), now bounded. If in addition the group curvature is bounded above, $-\nabla^2\ell_j\preceq G_j$ (automatic when the GLM weights are bounded, e.g. binomial with any link), then Cramér–Rao gives $V_j\succeq(P_b+G_j)^{-1}$, and with $S_\flat:=\sum_jH_j^\top P_b(P_b+G_j)^{-1}P_bH_j\succ0$,
$$
\boxed{\ \{\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\le2d\}\ \subseteq\ \widetilde C_d\ \subseteq\ \{\|\gamma-\gamma^\star\|^2_{S_\flat}\le2d\},
\qquad 0\prec S_\flat\preceq P_{11}^{\mathrm{RE}}\preceq P_{11}\ }
$$
— two concentric ellipsoids sandwiching $\widetilde C_d$ for every $d$. This uses (H3a) alone and is **prior-free**: $S(\gamma)$ is built from random-effect conditional covariances, which never involve $\Lambda_\gamma$; only the centre $\gamma^\star$ moves with the prior, not the shape of $\widetilde C_d$.

**Corollary 15 (mass bounds).** For every $d>0$,
$$
Q(\widetilde C_d)\ \ge\ \Pr(\chi^2_q\le2d),
\qquad
\pi_\gamma(\widetilde C_d^{\,c})\ \le\ \pi_\gamma\bigl(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d\bigr),
$$
from $\bar\Phi\ge0$ (Lemma 13(2)) via the complementarity identity: the $Q$-ellipsoid of squared radius $2d$ is inscribed in $\widetilde C_d$. Both hold for every GLM family with no closure assumption, and both are tail probabilities of the same statistic $\|\gamma-\gamma^\star\|^2_{P_{11}}\sim\chi^2_q$ under $Q$ — the second bound is what §6 uses for exhaustion (Theorem 2, item 1).

### 5.3 The attained constant $\varepsilon_d$

**Lemma 16 (the constant on $\widetilde C_d$ is the attained minimum).** For every $d>0$,
$$
\boxed{\ \min_{\gamma\in\widetilde C_d}\varepsilon(\gamma)\;=\;\varepsilon_d\;=\;e^{-d}\varepsilon(\gamma^\star)\;\in\;(0,1),\quad\text{attained on }\partial\widetilde C_d=\{\Psi=d\}.\ }
$$
The constant is therefore **sharp** — no larger constant minorizes on $\widetilde C_d$ — and $\widetilde C_d$ is the largest set carrying it. (Weierstrass, from continuity (Lemma 12) and compactness (Lemma 14), gives attainment and $\ge\varepsilon_d$ immediately since every point of $\widetilde C_d$ satisfies $\varepsilon\ge\varepsilon_d$ by definition; the intermediate value theorem along a ray from $\gamma^\star$, using coercivity (Lemma 14), places a point exactly on the level $\Psi=d$, giving $\le\varepsilon_d$.)

**Lemma 17 (the minorization condition holds on $\widetilde C_d$).**

**(a) Unrestricted chain, untruncated refresh.** $q(\gamma,A)\ge\varepsilon_d\,Q(A)$ for all $\gamma\in\widetilde C_d$ and all $A\in\mathcal B(\mathbb R^q)$: $\widetilde C_d$ is a small set of the original $\gamma$-chain, with minorizing measure the untruncated Gaussian $Q$ and no mass discount.

**(b) Restricted chain, truncated refresh.** Let $q(\cdot,\cdot\mid\widetilde C_d)$ be the two-block sweep restricted to $\widetilde C_d$ (Definition 1, §2.1) and $Q_{\widetilde C_d}:=Q(\cdot\mid\widetilde C_d)$. Then
$$
\boxed{\ q\bigl(\gamma,A\mid\widetilde C_d\bigr)\ \ge\ \varepsilon\,Q_{\widetilde C_d}(A),
\qquad
\varepsilon:=\varepsilon_d\,Q(\widetilde C_d)\ \ge\ \varepsilon(\gamma^\star)\,e^{-d}\Pr(\chi^2_q\le2d),\ }
$$
for all $\gamma\in\widetilde C_d$, $A\subseteq\widetilde C_d$ — both factors now explicit functions of $d$, apart from the single scalar $\varepsilon(\gamma^\star)$ (Lemma 11). (Truncating the refresh Gaussian to $\widetilde C_d$ only *raises* the transition density on the retained region; the mass discount $Q(\widetilde C_d)$ accounts exactly for the $Q$-mass $\varepsilon_d\,Q(\widetilde C_d^{\,c})$ that (a) would have deposited outside $\widetilde C_d$, unavailable to a chain confined there.) This is the pair $(\widetilde C_d,\varepsilon)$ instantiating items 1 and 3 of Theorem 2.

### 5.4 Flat-prior limit ($\Lambda_\gamma\downarrow0$)

Write $\varepsilon(\cdot\mid\tau)$, $\widetilde C_d(\tau)$, $\varepsilon_d(\tau)$ for the objects at prior scale $\tau$ (notation of §4.2); assume (H2), (H3a), (H3b) throughout. Every dependence on $\tau$ in §5.1–§5.3 enters only through $Q_\tau$ — the mixing law $\pi(\beta\mid\gamma,y)$ does not involve $\Lambda_\gamma$ at all — which is what lets each result below carry over from §5.1–§5.3 essentially unchanged, read at $\tau=0$.

**Lemma 18 (limiting profile: well-posed, positive, continuous).** At $\Lambda_\gamma=0$, $\varepsilon(\cdot\mid0)$ is well defined, lies in $(0,1]$ pointwise, and is continuous — Lemmas 10–12 read at $\tau=0$, legitimate because their proofs consume only $\Lambda_Q=P_{11}^{\mathrm{RE}}\succ0$ (Lemma 7, under (H3a)) and properties of $\pi(\beta\mid\gamma,y)$ that never involve the prior. Moreover $\varepsilon(\cdot\mid\tau)\to\varepsilon(\cdot\mid0)$ uniformly on compacts as $\tau\downarrow0$.

**Lemma 19 (limiting set: convex, compact, uniformly in the prior).** $\widetilde C_d(0)=\{\Psi_0\le d\}$ is convex, compact, and nonempty with $\gamma^{\star(0)}$ (Lemma 8, §4.2) in its interior. The sandwich of Lemma 14 holds with the *same*, $\tau$-free constants at every $\tau\ge0$, so there is a single compact $K_d\supseteq\widetilde C_d(\tau)$ for all $\tau\le\tau_0$, and $\widetilde C_d(\tau)\to\widetilde C_d(0)$ in Hausdorff distance. Compactness needs only (H3a) here; it is (H3b), via Lemma 8, that places the centre $\gamma^{\star(0)}$ — the *only* place (H3b) enters §5.

**Lemma 20 (limiting constant; continuity in the prior).**
$$
\min_{\widetilde C_d(0)}\varepsilon(\cdot\mid0)\;=\;\varepsilon_d(0)\;=\;e^{-d}\varepsilon\bigl(\gamma^{\star(0)}\mid0\bigr)\;\in\;(0,1),\quad\text{attained on }\partial\widetilde C_d(0),
$$
and $\varepsilon_d(\tau)\to\varepsilon_d(0)$ as $\tau\downarrow0$ — **the certified constant does not collapse as the prior becomes diffuse.** This is the precise sense in which §1.2's claim holds.

**Lemma 21 (minorization holds in the limit).** Under (H2), (H3a), (H3b), at $\Lambda_\gamma=0$: $q_0(\gamma,A)\ge\varepsilon_d(0)\,Q_0(A)$ for $\gamma\in\widetilde C_d(0)$, $A\in\mathcal B(\mathbb R^q)$ (unrestricted form), and
$$
q_0\bigl(\gamma,A\mid\widetilde C_d(0)\bigr)\ge\varepsilon^{(0)}\,Q_{0,\widetilde C_d}(A),
\qquad
\varepsilon^{(0)}:=\varepsilon_d(0)\,Q_0\bigl(\widetilde C_d(0)\bigr)\ge\varepsilon\bigl(\gamma^{\star(0)}\mid0\bigr)e^{-d}\Pr(\chi^2_q\le2d),
$$
for $\gamma\in\widetilde C_d(0)$, $A\subseteq\widetilde C_d(0)$. **Theorem 2 holds at $\Lambda_\gamma=0$**, with the flat-limit $\gamma$-marginal in place of $\pi_\gamma$.

**Remark 5.2 (division of labour between §5.1–§5.3 and §5.4).** Positivity, uniqueness, continuity, convexity, and — crucially — compactness of the certified set need only (H2)+(H3a) and are identical in both regimes, because they rest entirely on $\nabla^2\Psi=\sum_jH_j^\top P_bV_j(\gamma)P_bH_j$, which never involves $\Lambda_\gamma$. (H3b) enters exactly once, to place the centre $\gamma^{\star(0)}$ (Lemma 8); everything downstream of the centre then follows the proper-prior argument verbatim.

### 5.5 Closed form under Gaussian closure

**Lemma 22 (closed form for the deficiency; LMM only).** Assume the hypotheses of §5.1–§5.3, and in addition Gaussian closure: $q(\cdot\mid\gamma)=N(M(\gamma),\Sigma)$ with $\Sigma$ not depending on $\gamma$. Then
$$
\boxed{\ \mathcal D(\gamma)=\tfrac12\log\det\bigl(\Sigma^\star\Sigma^{-1}\bigr)-\tfrac12\bigl\|\gamma^\star-M(\gamma)\bigr\|^2_{(\Sigma-\Sigma^\star)^{-1}},\ }
$$
so the one scalar left unevaluated by §5.1–§5.3 is
$$
\varepsilon(\gamma^\star)=\sqrt{\det(\Sigma^\star\Sigma^{-1})}=\det(I+\tilde J)^{-1/2}\ \ge\ (1+\kappa)^{-q/2}
$$
($\tilde J$, $\kappa$ as in Lemma 3, §3.3), and the certified sets are exactly the concentric ellipsoids
$$
\widetilde C_d=\bigl\{\gamma:\|\gamma^\star-M(\gamma)\|^2_{(\Sigma-\Sigma^\star)^{-1}}\le2d\bigr\},
\qquad
\varepsilon_d\ \ge\ (1+\kappa)^{-q/2}e^{-d}.
$$
(Existence, uniqueness, and attainment of the minimiser are Lemmas 10–11, which already apply; closure lets the minimisation be solved in closed form because $\Lambda_q=\Sigma^{-1}$ no longer depends on $\gamma'$.)

**Remark 5.3 (what closure buys, and what does not survive without it).** Only Lemma 22 is closure-dependent. Every result of §5.1–§5.4 holds for every GLM family and already gives $\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)$ exactly; closure supplies a formula for the one remaining scalar $\varepsilon(\gamma^\star)$, in place of a $q$-dimensional convex program solved numerically, and an exact ellipsoid for $\widetilde C_d$ in place of the two-sided sandwich of Lemma 14. Outside closure, $q(\cdot\mid\gamma)=\int\phi_q(\cdot;m(\beta),\Sigma^\star)\pi(\beta\mid\gamma,y)\,d\beta$ is a Gaussian *mixture* with no single covariance — which is exactly why (Q1) in §4 is stated on precisions rather than covariances — and $\Lambda_q(\gamma'\mid\gamma)$ genuinely varies with $\gamma'$, so $\mathcal D$ has no closed form. Every qualitative conclusion of §5.1–§5.4 — positivity, continuity, concavity, compactness, the attained constant $\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)$, and its nondegeneracy as $\Lambda_\gamma\downarrow0$ — survives regardless; what is lost is only a closed-form value for $\varepsilon(\gamma^\star)$ and the exact ellipsoidal shape.

---

## 6. Proof of Theorem 2

Two ingredients remain beyond §4–§5: that the certified sets $\widetilde C_d$ exhaust $\mathbb R^q$ as $d\to\infty$, so a set can actually be chosen to meet any tail-mass budget $\delta$; and that restricting the two-block sweep to a set $C$ produces a kernel for which $\pi_\gamma(\cdot\mid C)$ is genuinely stationary — the standing assumption Definition 1 (§2.1) required but did not prove.

### 6.1 Assembly of §4–§5

**Lemma 23 (exhaustion).** Under (H1) and the hypotheses of §5.1–§5.3 (or, at $\Lambda_\gamma=0$, of §5.4), the family $\widetilde C_d$ is nondecreasing in $d$ with $\bigcup_{d>0}\widetilde C_d=\mathbb R^q$. Consequently $\pi_\gamma(\widetilde C_d^{\,c})\downarrow0$ as $d\uparrow\infty$, and for every $\delta>0$ there is $d(\delta)<\infty$ with $\pi_\gamma\bigl(\widetilde C_{d(\delta)}^{\,c}\bigr)<\delta$. (Positivity of $\varepsilon$ pointwise — Lemma 11, or 18 at $\tau=0$ — gives $\Psi(\gamma)<\infty$ for every $\gamma$, hence every point eventually enters $\widetilde C_d$; continuity from above of the probability measure $\pi_\gamma$ finishes it. This is the same pointwise positivity that Lemma 11 already established — no separate tightness argument or compactness input is needed.)

**Lemma 24 (restricted kernel and stationarity).** Let $C\subseteq\mathbb R^q$ be measurable with $\pi_\gamma(C)>0$, and let $q(\cdot\mid\cdot\,;C)$ be the $\gamma$-marginal of the two-block Gibbs kernel *of the truncated target* $\pi(\cdot\mid\tilde C)$, $\tilde C:=C\times\mathbb R^{Jp_{\mathrm{re}}}$ — an exact draw $\beta\sim\pi(\beta\mid\gamma,y)$ followed by $\gamma'\sim N(m(\beta),\Sigma^\star)$ truncated to $C$. Then $q(\gamma,C\mid C)=1$ for $\gamma\in C$, and $\pi_\gamma(\cdot\mid C)$ is invariant for $q(\cdot\mid\cdot\,;C)$. (Each block is drawn from the exact conditional of $\pi(\cdot\mid\tilde C)$ given the other block, so the restricted sweep is Gibbs *for* the truncated target, not the unrestricted sweep with excursions repaired after the fact — a construction of the latter kind does not in general preserve $\pi(\cdot\mid\tilde C)$.) This discharges the standing assumption of Definition 1, for $C=\widetilde C_d$ or any measurable subset of it.

*Proof of Theorem 2.* Fix $\delta>0$.

1. $Q=N(\gamma^\star,\Sigma^\star)$ exists — Corollary 6 if $\Lambda_\gamma\succ0$, Corollary 9 under (H3a)+(H3b) if $\Lambda_\gamma=0$ — independent of $\delta$.
2. By Lemma 23, choose $d=d(\delta)$ with $\pi_\gamma(\widetilde C_d^{\,c})<\delta$, and set $C:=\widetilde C_d$.
3. By Lemma 14 (or 19 at $\tau=0$), $C$ is compact and convex; by Lemma 16 (or 20), $\min_C\varepsilon=\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\in(0,1)$, attained.
4. By Lemma 24, $\pi_\gamma(\cdot\mid C)$ is stationary for $q(\cdot\mid\cdot\,;C)$.
5. By Lemma 17(b) (or 21), $q(\gamma,A\mid C)\ge\varepsilon\,Q_C(A)$ for $\gamma\in C$, $A\subseteq C$, with $\varepsilon:=\varepsilon_d\,Q(C)\in(0,1]$ and $Q_C:=Q(\cdot\mid C)$.

Items 2–5 instantiate Theorem 2's items 1–3 exactly. Applying Theorem 1 (§2.4) with this $(C,\varepsilon)$: for $\gamma\in C$,
$$
\|q_n(\gamma,\cdot\mid C)-\pi_\gamma\|_{TV}\ \le\ (1-\varepsilon)^n+\pi_\gamma(C^c)\ <\ (1-\varepsilon)^n+\delta. \qquad\blacksquare
$$

**Remark 6.1 (how $\varepsilon$ degrades as $\delta$ shrinks).** Both terms competing in the choice of $d$ are now explicit: $\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)$ is decreasing in $d$, while $\pi_\gamma(\widetilde C_d^{\,c})$ is decreasing in $d$ toward $0$. Under (H2), $\pi_\gamma$ is log-concave, so with $\Sigma_\pi$ its (finite) covariance and $\rho:=\lambda_{\max}(\Sigma_\pi P_{11})\ge1$, the tail of $\Psi$ is sub-exponential and $\pi_\gamma(\widetilde C_d^{\,c})\lesssim e^{-d/\rho}$ (standard log-concave concentration; the exact constant is deferred to Appendix A.4). Choosing $d=\rho\log(1/\delta)$ to meet the budget then gives, for every GLM family,
$$
\boxed{\ \varepsilon\ \gtrsim\ \varepsilon(\gamma^\star)\,\delta^{\rho}.\ }
$$
The minorization constant degrades only **polynomially** in the tail-mass budget $\delta$, with exponent $\rho$ measuring the mismatch between the posterior's spread and the refresh metric $P_{11}$ — cheap whenever $\rho=O(1)$. The metric here is $P_{11}\succeq P_{11}^{\mathrm{RE}}\succ0$, not $\Lambda_\gamma$, so unlike an envelope built from the prior alone, this survives the flat-prior limit under (H3a). Under Gaussian closure, $\rho=(1-\kappa)^{-1}$ exactly (Lemma 22, §5.5).

### 6.2 Lifting to the joint $(\gamma,\beta)$ chain

Theorem 2 and its proof concern only the $\gamma$-marginal chain. The following closes the gap to the full two-block state space at no cost in rate.

**Lemma 25 (lifting to the joint chain).** The joint two-block chain is a de-initialization of the $\gamma$-chain: $(\gamma_n,\beta_n)$ is obtained from $\gamma_{n-1}$ by drawing $\beta_n\sim\pi(\beta\mid\gamma_{n-1},y)$ (exact) and then $\gamma_n\sim\pi(\gamma\mid\beta_n)$ (exact, or restricted to $C$ as in §6.1). Since the extra $\beta$-draw is from the exact conditional of the target, the data-processing inequality for total variation gives
$$
\|\mathcal L(\gamma_n,\beta_n)-\pi\|_{TV}\ \le\ \|\mathcal L(\gamma_{n-1})-\pi_\gamma\|_{TV}.
$$
Combined with Theorem 2, for $\gamma_0\in C$,
$$
\|\mathcal L(\gamma_n,\beta_n)-\pi\|_{TV}\ \le\ (1-\varepsilon)^{n-1}+\delta.
$$
Every bound proved for the $\gamma$-chain therefore transfers to the joint chain with the same constants and a one-step offset. Scoping Theorem 2 to the $\gamma$-marginal costs nothing at the end, and buys a construction on a state space of dimension $q$ rather than $q+Jp_{\mathrm{re}}$ throughout §4–§6 — which is why $\varepsilon$ carries no dependence on the number of groups $J$.

---

## 7. Scope and Extensions

*(Content to follow.)*

### 7.1 Symmetric case: sharper constants

### 7.2 What is and is not certified

### 7.3 Open problems

---

## Appendix — Proofs

### A.1 Proofs for §3.2 (propriety of the flat-prior posterior) and §3.3 (mean map)

**Proof of Lemma 2 (propriety of the flat-$\gamma$ posterior).**

*1. $g$ is concave and bounded above.* Each $\ell_j$ is concave by (P1) and the quadratic term is concave, so $g$ is concave on $\mathbb R^q\times\mathbb R^{Jp_{\mathrm{re}}}$. By (P4), $\ell_j\le\ell_j(\hat\beta_j)<\infty$ and the quadratic term is $\le0$, so $g\le\sum_j\ell_j(\hat\beta_j)<\infty$.

*2. Each $\ell_j$ has compact superlevel sets.* By (P2), $\beta_j\mapsto\eta_j=Z_j\beta_j$ is injective, and with the IRLS weights strictly positive (P1) the log-likelihood is strictly concave in $\beta_j$; by (P4) its maximum is attained at a finite $\hat\beta_j$. A strictly concave function with an attained maximum decays at least linearly along every ray: for unit $v$,
$$
\ell_j^\infty(v):=\lim_{t\to\infty}\frac{\ell_j(\hat\beta_j+tv)-\ell_j(\hat\beta_j)}{t}<0,
\qquad
c_j:=-\max_{\|v\|=1}\ell_j^\infty(v)>0
$$
(the max is attained since $\ell_j^\infty$ is concave, positively homogeneous, and u.s.c. on the compact unit sphere). This is exactly "finite, unique MLE with positive IRLS weights" — no separation, no flat direction.

*3. Every direction recedes.* Fix $(d_\gamma,d_\beta)\ne0$.
- *Case A: $d_{\beta,j}\ne H_jd_\gamma$ for some $j$.* The quadratic term $-\frac12\sum_j\|\beta_j+td_{\beta,j}-H_j(\gamma+td_\gamma)\|^2_{P_b}$ decreases quadratically in $t$ (since $P_b\succ0$), while $\sum_j\ell_j$ is bounded above by Step 1. So $g\to-\infty$ quadratically.
- *Case B: $d_{\beta,j}=H_jd_\gamma$ for all $j$.* The quadratic term is constant along the ray — the only direction propriety can fail. Since $d_\gamma\ne0$ here (else $d_\beta=0$, excluded), (P3) gives some $j_0$ with $H_{j_0}d_\gamma\ne0$; by Step 2, $\ell_{j_0}(\beta_{j_0}+tH_{j_0}d_\gamma)\to-\infty$ at least linearly with slope $\le-c_{j_0}\|H_{j_0}d_\gamma\|$, while every other $\ell_k$ stays bounded above. So $g\to-\infty$ at least linearly.

*4. Integrability.* $g$ is concave, bounded above, with strictly negative recession function in every direction; by compactness of the unit sphere there are $c>0$, $x_0=(\hat\gamma,\hat\beta)$ with $g(x)\le g(x_0)-c\|x-x_0\|$ for $\|x-x_0\|$ large, so
$$
\int e^g\ \le\ \mathrm{const}\cdot e^{g(x_0)}\int e^{-c\|x-x_0\|}dx\ <\ \infty. \qquad\blacksquare
$$

*Sharpness.* (P4) is sufficient, not necessary: the exact requirement in Case B is that the pooled profile $\gamma\mapsto\sum_j\ell_j(H_j\gamma)$ have compact superlevel sets — groups can borrow strength, and a single separated group need not destroy propriety if the others pin down $\gamma$. (P4) is therefore a conservative, per-group certificate: it certifies propriety but failing it does not prove impropriety. (Specializing to binomial with a flat prior recovers the classical fact that posterior propriety under a flat prior corresponds to existence of the MLE, i.e. absence of complete or quasi-complete separation.)

**Proof of Lemma 3 (the mean map has a unique fixed point).**

*(a) Jacobian via exponential tilt.* $\pi(\beta_j\mid\gamma,y)\propto\exp\{\ell_j(\beta_j)+(P_bH_j\gamma)^\top\beta_j-\frac12\beta_j^\top P_b\beta_j\}$ is an exponential family in natural parameter $\theta_j=P_bH_j\gamma$ with sufficient statistic $\beta_j$, so $\partial b_j/\partial\theta_j=V_j(\gamma)$ and, by the chain rule, $\partial b_j/\partial\gamma=V_jP_bH_j$. Hence
$$
J(\gamma):=\frac{\partial M}{\partial\gamma}=P_{11}^{-1}\sum_jH_j^\top P_bV_j(\gamma)P_bH_j=P_{11}^{-1/2}\tilde J(\gamma)P_{11}^{1/2},
\qquad
\tilde J:=P_{11}^{-1/2}\Bigl(\sum_jH_j^\top P_bV_jP_bH_j\Bigr)P_{11}^{-1/2}\succeq0.
$$

*(b) Brascamp–Lieb upper bound.* Under (H2), $-\nabla^2\log\pi(\beta_j\mid\gamma,y)=-\nabla^2\ell_j+P_b\succeq P_b$: the conditional is $P_b$-strongly log-concave, so Brascamp–Lieb gives $V_j\preceq P_b^{-1}$, hence $\sum_jH_j^\top P_bV_jP_bH_j\preceq\sum_jH_j^\top P_bH_j=P_{11}^{\mathrm{RE}}$, i.e. $\tilde J\preceq I-P_{11}^{-1/2}\Lambda_\gamma P_{11}^{-1/2}$.

*(c) Strictness.* If $\Lambda_\gamma\succ0$: $\kappa:=\sup_\gamma\lambda_{\max}(\tilde J(\gamma))\le1-\lambda_{\min}(P_{11}^{-1/2}\Lambda_\gamma P_{11}^{-1/2})<1$ unconditionally — the prior alone supplies the gap. At $\Lambda_\gamma=0$, one instead needs $V_j\prec P_b^{-1}$ strictly: if $-\nabla^2\ell_j\succeq\epsilon_jP_b$ uniformly, $V_j\preceq((1+\epsilon_j)P_b)^{-1}$ and $\kappa\le1/(1+\min_j\epsilon_j)<1$.

*(d) Contraction and Banach.* $M(\gamma_1)-M(\gamma_2)=\int_0^1J(\gamma_t)(\gamma_1-\gamma_2)\,dt$ along the segment $\gamma_t$, so $P_{11}^{1/2}(M(\gamma_1)-M(\gamma_2))=\int_0^1\tilde J(\gamma_t)P_{11}^{1/2}(\gamma_1-\gamma_2)\,dt$ and $\|M(\gamma_1)-M(\gamma_2)\|_{P_{11}}\le\kappa\|\gamma_1-\gamma_2\|_{P_{11}}$. $\mathbb R^q$ with $\|\cdot\|_{P_{11}}$ is complete, so Banach's fixed-point theorem gives existence, uniqueness, and geometric convergence of $\gamma^{(k+1)}=M(\gamma^{(k)})$. $\blacksquare$

---

### A.2 Proofs for §4

**Proof of Lemma 4 (the covariance of $Q$ exists).** $P_{11}^{\mathrm{RE}}$ is a sum of terms $H_j^\top P_bH_j\succeq0$, so $P_{11}\succeq\Lambda_\gamma\succ0$; inverting reverses the order. $\blacksquare$ (No rank condition is used: a proper prior supplies invertibility on its own.)

**Proof of Lemma 5 (the centre of $Q$ exists).** The Hessian formula follows from Lemma 3(a); Brascamp–Lieb (Lemma 3(b)) gives $V_j\preceq P_b^{-1}$, hence $\sum_jH_j^\top P_bV_jP_bH_j\preceq P_{11}^{\mathrm{RE}}$, which is the displayed bound $\nabla^2\Phi\succeq\Lambda_\gamma\succ0$. A $\Lambda_\gamma$-strongly convex function satisfies $\Phi(\gamma)\ge\Phi(\gamma_0)+\nabla\Phi(\gamma_0)^\top(\gamma-\gamma_0)+\frac12\|\gamma-\gamma_0\|^2_{\Lambda_\gamma}$, hence is coercive, so its minimum is attained; strict convexity makes it unique and finite. Stationarity $\gamma^\star=M(\gamma^\star)$ is the first-order condition. $\blacksquare$ (Lemma 3 also reaches $\gamma^\star$, via Banach; the two routes agree, but the variational route here is the one that survives $\Lambda_\gamma=0$, where $\kappa$ may equal $1$.)

**Proof of Corollary 6 (the refresh density $Q$ exists).** Combine Lemmas 4 and 5; positivity and continuity of $q_Q$ are immediate from the Gaussian formula, well posed since $\det P_{11}>0$. Precision dominance (Q1) is the identity derived in §4 (differentiate the mixture defining $q(\cdot\mid\gamma)$ twice in $\gamma'$: $\Lambda_Q-\Lambda_q(\gamma'\mid\gamma)=\Lambda_Q\,\mathrm{Cov}(m(\beta)\mid\gamma,\gamma',y)\,\Lambda_Q\succeq0$, the covariance under the bridge law $\pi(\beta\mid\gamma,\gamma',y)$), which holds for every GLM family with no further hypothesis. $\blacksquare$

**Proof of Lemma 7 (the limiting covariance exists iff (H3a)).** Monotonicity of $P_{11}(\tau)=\tau\bar\Lambda+P_{11}^{\mathrm{RE}}$ in $\tau$ is immediate, and inversion reverses the psd order. Under (H3a), $P_{11}^{\mathrm{RE}}$ is nonsingular and matrix inversion is continuous there, giving $\Sigma^\star(\tau)\to\Sigma^{\star(0)}$. Without (H3a): for $0\ne v\in\ker P_{11}^{\mathrm{RE}}$, Cauchy–Schwarz gives $(v^\top v)^2=(v^\top P_{11}^{1/2}P_{11}^{-1/2}v)^2\le(v^\top P_{11}v)(v^\top P_{11}^{-1}v)$, so $v^\top\Sigma^\star(\tau)v\ge(v^\top v)^2/(\tau\,v^\top\bar\Lambda v)\to\infty$. $\blacksquare$ ($\Sigma^{\star(0)}=(P_{11}^{\mathrm{RE}})^{-1}$ is the GLS covariance for the hyper-regression of the random effects on the hyper-design with weight $\Psi^{-1}$; the prior term drops out of $M$ entirely in the limit.)

**Proof of Lemma 8 (the limiting centre is a finite vector).**
*Coercivity.* Under (H2)+(H3a)+(H3b), Lemma 2 (§3.2/A.1) makes the flat-limit joint posterior proper, and its Step 4 supplies an exponential envelope $e^{-c\|x-x_0\|}$; marginalizing, $\Phi_0$ grows at least linearly and is coercive. A convex coercive function attains its minimum.
*Uniqueness.* $\nabla^2\Phi_0=P_{11}^{\mathrm{RE}}-\sum_jH_j^\top P_bV_jP_bH_j\succeq0$ always (Lemma 3(a)–(b)); strict positivity — needed for uniqueness — holds wherever the group conditional covariances $V_j$ sit strictly below their Brascamp–Lieb ceiling $P_b^{-1}$. This is supplied, with no measure-zero exclusion, by (P2) — $Z_j$ full column rank for every $j$, already part of Lemma 2's hypotheses (`groupef.rank`) — together with (H2)'s requirement that the GLM weights $w_{j,i}(\eta)>0$ at every finite $\eta$ (true for every canonical-link family, e.g. $\mu(1-\mu)>0$ for logit, $\mu=e^\eta>0$ for Poisson-log). With $V(\beta_j):=-\ell_j(\beta_j)+\frac12(\beta_j-H_j\gamma)^\top P_b(\beta_j-H_j\gamma)$, Brascamp–Lieb's variance formula gives, for any $a$, $a^\top V_j(\gamma)a=\mathrm{Var}_\pi(a^\top\beta_j)\le E_\pi[a^\top(\nabla^2V(\beta_j))^{-1}a]$ — requiring only $\nabla^2V\succ0$ a.e., not a uniform lower bound. Since $\nabla^2V(\beta_j)=Z_j^\top\mathrm{diag}(w(\eta_j))Z_j+P_b$ and $Z_j$ has full column rank with $w>0$ at every finite $\eta_j$, $\nabla^2V(\beta_j)\succ P_b$ — strictly, at *every* finite $\beta_j$, not almost every — so $(\nabla^2V(\beta_j))^{-1}\prec P_b^{-1}$ pointwise everywhere; taking $E_\pi$ of a pointwise-strict inequality under a measure of full support preserves strictness, giving $V_j(\gamma)\prec P_b^{-1}$ as a matrix inequality. For $v\ne0$, (H3a) gives some $j_0$ with $H_{j_0}v\ne0$, and $v^\top\nabla^2\Phi_0v\ge(P_bH_{j_0}v)^\top(P_b^{-1}-V_{j_0}(\gamma))(P_bH_{j_0}v)>0$. So $\nabla^2\Phi_0(\gamma)\succ0$ for every $\gamma$, and a strictly convex coercive function has a unique minimiser.
*Convergence of the modes.* $\Phi_\tau(\gamma)=\Phi_0(\gamma)+\frac\tau2(\gamma-\mu_0)^\top\bar\Lambda(\gamma-\mu_0)\ge\Phi_0(\gamma)$, so minimality gives $\Phi_0(\gamma^\star(\tau))\le\Phi_\tau(\gamma^\star(\tau))\le\Phi_\tau(\gamma^{\star(0)})=\Phi_0(\gamma^{\star(0)})+C\tau$ with $C=\frac12(\gamma^{\star(0)}-\mu_0)^\top\bar\Lambda(\gamma^{\star(0)}-\mu_0)$. Every $\gamma^\star(\tau)$ therefore lies in the sublevel set $\{\Phi_0\le\Phi_0(\gamma^{\star(0)})+C\tau\}$, bounded by coercivity; any limit point $\bar\gamma$ satisfies $\Phi_0(\bar\gamma)\le\Phi_0(\gamma^{\star(0)})$ by continuity, hence $\bar\gamma=\gamma^{\star(0)}$ by uniqueness, and a bounded family with one limit point converges. $\blacksquare$

*(Why Banach is unavailable here.* At $\tau=0$, Lemma 3(c) leaves only $\kappa\le1$, and for families whose GLM weights vanish in the tails — logit, probit at large $|\eta|$ — one has $\sup_\gamma\kappa(\gamma)=1$ exactly, so the contraction argument gives nothing globally; Lemma 8 is the variational replacement. If uniform group curvature $-\nabla^2\ell_j\succeq\epsilon_jP_b$ is assumed, Lemma 3(c) instead gives the prior-free ceiling $\kappa_0\le1/(1+\min_j\epsilon_j)<1$ and Banach is restored with a rate that does not degrade as $\tau\downarrow0$.)*

**Proof of Corollary 9 (the limiting refresh density exists).** Nondegeneracy is Lemma 7 plus finiteness of $\gamma^{\star(0)}$ (Lemma 8). Since $\gamma^\star(\tau)\to\gamma^{\star(0)}$ and $\Sigma^\star(\tau)\to\Sigma^{\star(0)}\succ0$, the Gaussian densities converge pointwise; both are probability densities, so Scheffé's lemma gives $L^1$ (equivalently TV) convergence. $\blacksquare$

---

### A.3 Proofs for §5

Throughout §5.1–§5.3, assume (H2), (H3a), $\Lambda_\gamma\succ0$, so $Q$ exists by Corollary 6.

**Proof of Lemma 10 (uniqueness).** Differentiating the mixture defining $q(\cdot\mid\gamma)$ twice in $\gamma'$ gives $\nabla^2_{\gamma'}g=\Lambda_Q-\Lambda_q(\gamma'\mid\gamma)=\Lambda_Q\,\mathrm{Cov}(m(\beta)\mid\gamma,\gamma',y)\,\Lambda_Q$, the covariance under the bridge law. With $m(\beta)=L\beta+a$, $L=P_{11}^{-1}[H_1^\top P_b,\dots,H_J^\top P_b]$ of full row rank $q$ under (H3a), and the bridge law of $\beta$ nondegenerate, $\mathrm{Cov}(m\mid\cdot)=L\,\mathrm{Cov}(\beta\mid\cdot)L^\top\succ0$; conjugating by $\Lambda_Q\succ0$ preserves positivity. A strictly convex function has at most one minimiser. $\blacksquare$

**Proof of Lemma 11 (existence and positivity).** Let $P_\gamma$ be the law of $m(\beta)$ under $\pi(\beta\mid\gamma,y)$. Dividing the two densities and cancelling the Gaussian normalisers and the $\gamma'$-quadratics, $g(\gamma'\mid\gamma)=\log\int e^{(m-\gamma^\star)^\top\Lambda_Q\gamma'}\,d\nu_\gamma(m)$ where $d\nu_\gamma(m)=e^{-\frac12m^\top\Lambda_Qm+\frac12{\gamma^\star}^\top\Lambda_Q\gamma^\star}dP_\gamma(m)$ — so in $t=\Lambda_Q\gamma'$, $g$ is the cumulant generating function of the finite positive measure $\nu_\gamma$, recentred at $\gamma^\star$: finite everywhere (the Gaussian factor dominates every exponential) and convex. For coercivity, restrict to a half-space: for unit $u$, $c>0$, $s>0$,
$$
g(s\Lambda_Q^{-1}u\mid\gamma)\ \ge\ sc+\log\nu_\gamma\bigl(\{m:(m-\gamma^\star)^\top u\ge c\}\bigr)\ \xrightarrow[s\to\infty]{}\ +\infty
$$
whenever that half-space carries positive $\nu_\gamma$-mass. Under (H3a), $m(\beta)$ is an affine map of $\beta$ with linear part of rank $q$, and $\pi(\beta\mid\gamma,y)>0$ on all of $\mathbb R^{Jp_{\mathrm{re}}}$, so $P_\gamma$ (hence $\nu_\gamma$) has full support $\mathbb R^q$ and every half-space qualifies. A finite convex coercive function attains its minimum, which is therefore $>-\infty$; exponentiating gives $\varepsilon(\gamma)>0$. $\blacksquare$

**Proof of Lemma 12 (continuity).** Joint continuity of $g$ gives upper semicontinuity of $\log\varepsilon=\inf_{\gamma'}g$ for free. For the reverse direction, fix compact $K\subset\mathbb R^q$; working in $\beta$-space (where $\pi(\beta\mid\gamma,y)>0$ everywhere), write $t=\Lambda_Q\gamma'$, $u=t/|t|$, $m(\beta)=L\beta+a$, so $(m(\beta)-\gamma^\star)^\top t=|t|(\beta^\top L^\top u+(a-\gamma^\star)^\top u)$. Since $L$ has full row rank under (H3a), $\sigma_0:=\sigma_{\min}(L)>0$ and $v(u):=L^\top u/|L^\top u|$ is a well-defined unit vector depending continuously on $u$; picking $\rho$ large enough that $c:=\sigma_0\rho/2-|a-\gamma^\star|>0$, on the cone-cap $B_v:=\{\beta:|\beta|\le\rho,\ \beta^\top v\ge\rho/2\}$ one has $(m(\beta)-\gamma^\star)^\top u\ge c$, so restricting the integral to $B_{v(u)}$,
$$
g(\gamma'\mid\gamma)\ \ge\ c|t|+\log\bigl(w_{\min}\,\pi(B_{v(u)}\mid\gamma,y)\bigr),
$$
where $w_{\min}>0$ lower-bounds the continuous positive Gaussian weight on $\{|\beta|\le\rho\}$. Since $(v,\gamma)\mapsto\pi(B_v\mid\gamma,y)$ is continuous (dominated convergence; the cone-cap boundary is Lebesgue-null) and strictly positive on the compact $\{|v|=1\}\times K$, $p_K:=\min\pi(B_v\mid\gamma,y)>0$ there; taking $\eta_K:=w_{\min}p_K$, $c_K:=c\,\lambda_{\min}(\Lambda_Q)$ gives $g(\gamma'\mid\gamma)\ge c_K|\gamma'|+\log\eta_K$ for all $\gamma\in K$, uniformly. Since $\log\varepsilon\le0$, every minimiser satisfies $|\gamma'_\star(\gamma)|\le R_K:=(\log(1/\eta_K))_+/c_K$, the same $R_K$ for all $\gamma\in K$; so on $K$, $\log\varepsilon(\gamma)=\min_{|\gamma'|\le R_K}g(\gamma'\mid\gamma)$ is a minimum of a jointly continuous function over a fixed compact set, and Berge's maximum theorem gives continuity. $K$ was arbitrary. $\blacksquare$

**Proof of Lemma 13 (concavity; closed-form derivatives).** Write $\pi(\beta\mid\gamma,y)=e^{L(\beta,\gamma)}/Z(\gamma)$, $L(\beta,\gamma)=\sum_j\ell_j(\beta_j)-\frac12\sum_j(\beta_j-H_j\gamma)^\top P_b(\beta_j-H_j\gamma)$.

*Score identity.* $\nabla_\gamma L=\sum_jH_j^\top P_b\beta_j-P_{11}^{\mathrm{RE}}\gamma=P_{11}m(\beta)-\Lambda_\gamma\mu_0-P_{11}^{\mathrm{RE}}\gamma$; since $\nabla_\gamma\log Z=E_{\beta\mid\gamma}[\nabla_\gamma L]$, subtracting gives $s(\beta,\gamma):=\nabla_\gamma\log\pi(\beta\mid\gamma,y)=P_{11}(m(\beta)-M(\gamma))$.

*Bridge representation.* Differentiating $q(\gamma'\mid\gamma)=\int\phi_q(\gamma';m(\beta),\Sigma^\star)\pi(\beta\mid\gamma,y)\,d\beta$ under the integral and dividing by $q$, $\nabla_\gamma\log q(\gamma'\mid\gamma)=E[s(\beta,\gamma)\mid\gamma,\gamma',y]=P_{11}(\bar m(\gamma,\gamma')-M(\gamma))$, the expectation under the bridge law $\pi(\beta\mid\gamma,\gamma',y)\propto\phi_q(\gamma';m(\beta),\Sigma^\star)\pi(\beta\mid\gamma,y)$, $\bar m:=E[m(\beta)\mid\gamma,\gamma',y]$. Likewise $\nabla_{\gamma'}\log q=P_{11}(\bar m-\gamma')$ and $\nabla_{\gamma'}\log q_Q=-P_{11}(\gamma'-\gamma^\star)$, so $\nabla_{\gamma'}g=P_{11}(\bar m(\gamma,\gamma')-\gamma^\star)$: the minimiser is characterised by $\bar m(\gamma,\gamma'_\star(\gamma))=\gamma^\star$.

*Envelope theorem.* The minimiser is unique (Lemma 10), interior, and locally bounded (Lemma 12), and $\nabla_\gamma g$ is continuous, so Danskin's theorem gives $\mathcal D\in C^1$ with $\nabla\mathcal D(\gamma)=\nabla_\gamma g(\gamma'\mid\gamma)|_{\gamma'=\gamma'_\star(\gamma)}$ — the term $-\log q_Q(\gamma')$ does not involve $\gamma$, and the chain-rule term through $\gamma'_\star(\gamma)$ vanishes since $\nabla_{\gamma'}g=0$ there. By the bridge representation at $\gamma'=\gamma'_\star(\gamma)$, where $\bar m=\gamma^\star$: $\nabla\mathcal D(\gamma)=P_{11}(\gamma^\star-M(\gamma))$.

*Second derivative.* Differentiating again, only $M$ depends on $\gamma$: $\nabla^2\mathcal D=-P_{11}\,\partial M/\partial\gamma=-P_{11}J(\gamma)=-\sum_jH_j^\top P_bV_jP_bH_j\preceq0$ by Lemma 3(a). $\blacksquare$

Item (1): concavity, strict under (H3a) with $V_j\succ0$, gives convex superlevel sets. Item (2): $\nabla\mathcal D=0\iff M(\gamma)=\gamma^\star\iff\gamma=\gamma^\star$ by uniqueness of the fixed point (Lemma 3). Item (3): comparing $\nabla\mathcal D=P_{11}(\gamma^\star-M(\gamma))$ with $\nabla\Phi=P_{11}(\gamma-M(\gamma))$ gives $\nabla(\Phi-\mathcal D)=P_{11}(\gamma-\gamma^\star)$, i.e. $\mathcal D(\gamma)-\mathcal D(\gamma^\star)=\bar\Phi(\gamma)-\frac12\|\gamma-\gamma^\star\|^2_{P_{11}}$ exactly, integrating from $\gamma^\star$ where both sides vanish.

**Proof of Lemma 14 (coercivity; compactness).** By Lemma 13, $\nabla\Psi(\gamma^\star)=0$, $\nabla^2\Psi=P_{11}J=S(\gamma):=\sum_jH_j^\top P_bV_jP_bH_j\succeq0$.

*$S(\gamma)\succ0$ everywhere.* Under (H1)+(H2), $\pi(\beta_j\mid\gamma,y)$ is a proper log-concave density, strictly positive on all of $\mathbb R^{p_{\mathrm{re}}}$, hence not supported in any hyperplane, with finite second moments, so $V_j(\gamma)\succ0$. For $u\ne0$, (H3a) gives $u^\top P_{11}^{\mathrm{RE}}u=\sum_j|H_ju|^2_{P_b}>0$, so $H_ju\ne0$ for some $j$; with $v:=P_bH_ju\ne0$, $u^\top S(\gamma)u\ge v^\top V_j(\gamma)v>0$.

*Linear growth along every ray.* Fix unit $u$, $\psi_u(s):=\Psi(\gamma^\star+su)$. Then $\psi_u'(s)=\int_0^su^\top S(\gamma^\star+tu)u\,dt$ is nondecreasing in $s$ (nonnegative integrand), so for $s\ge1$, $\psi_u'(s)\ge\psi_u'(1)\ge c$, giving $\psi_u(s)=\int_0^s\psi_u'\ge c(s-1)$, where $c:=\min_{|u|=1}\int_0^1u^\top S(\gamma^\star+tu)u\,dt$.

*$c>0$.* $V_j$ is the Hessian of a log-partition function, hence smooth in $\gamma$, so $(u,t)\mapsto u^\top S(\gamma^\star+tu)u$ is continuous and, by the first part, strictly positive on the compact $\{|u|=1\}\times[0,1]$; its $t$-integral is a continuous strictly positive function of $u$, whose minimum over the compact unit sphere is attained and positive.

Combining, $\Psi(\gamma)\ge c(|\gamma-\gamma^\star|-1)$ for $|\gamma-\gamma^\star|\ge1$, so $\widetilde C_d\subseteq\{|\gamma-\gamma^\star|\le1+d/c\}$: closed (Lemma 12), convex (Lemma 13), now bounded, hence compact.

*Ellipsoid sandwich.* Brascamp–Lieb gives $V_j\preceq P_b^{-1}$, so $S\preceq P_{11}^{\mathrm{RE}}$ and $\{\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\le2d\}\subseteq\widetilde C_d$. If additionally $-\nabla^2\ell_j\preceq G_j$ (automatic when GLM weights are bounded above: binomial with any link gives $G_j=\frac14Z_j^\top N_jZ_j$), Cramér–Rao for the location family $\pi(\beta_j\mid\gamma,y)$ gives $V_j\succeq(P_b+G_j)^{-1}$ (its negative-log-density Hessian is $\preceq G_j+P_b$, so its Fisher information for location satisfies $E[\nabla^2U]\preceq G_j+P_b$); with $S_\flat:=\sum_jH_j^\top P_b(P_b+G_j)^{-1}P_bH_j\succ0$, integrating the resulting constant Hessian bounds twice from $\gamma^\star$ (where $\Psi,\nabla\Psi$ vanish) gives $\frac12\|\gamma-\gamma^\star\|^2_{S_\flat}\le\Psi(\gamma)\le\frac12\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}$, i.e. the two-sided sandwich. $\blacksquare$

**Proof of Corollary 15 (mass bounds).** By Lemma 13(3), $\bar\Phi\ge0$ (since $\gamma^\star$ minimises $\Phi$) gives $\Psi(\gamma)\le\frac12\|\gamma-\gamma^\star\|^2_{P_{11}}$, so the $Q$-ellipsoid of squared radius $2d$ is inscribed in $\widetilde C_d$; since $\|\gamma-\gamma^\star\|^2_{P_{11}}\sim\chi^2_q$ under $Q$, $Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)$, and identically for $\pi_\gamma(\widetilde C_d^{\,c})$ using the same statistic under $\pi_\gamma$. $\blacksquare$

**Proof of Lemma 16 (attained minimum).** $\ge\varepsilon_d$ holds because every point of $\widetilde C_d$ satisfies $\varepsilon\ge\varepsilon_d$ by definition. For $\le$: fix unit $u$; $s\mapsto\Psi(\gamma^\star+su)$ is continuous, $0$ at $s=0$, $\to\infty$ as $s\to\infty$ (Lemma 14), so by the intermediate value theorem there is $s_d$ with $\Psi(\gamma^\star+s_du)=d$ exactly — that point lies in $\widetilde C_d$ and has $\varepsilon=\varepsilon_d$ there, so the minimum is $\le\varepsilon_d$. Existence and attainment of the minimum itself follow from Weierstrass: $\varepsilon$ is continuous (Lemma 12) and $\widetilde C_d$ compact (Lemma 14). $\blacksquare$

**Proof of Lemma 17 (minorization on $\widetilde C_d$).**

*(a)* By definition of the infimum, $q(\gamma'\mid\gamma)\ge\varepsilon(\gamma)q_Q(\gamma')$ for every pair, and on $\widetilde C_d$, $\varepsilon(\gamma)\ge\varepsilon_d$ (Lemma 16); the bound holds for all $\gamma'\in\mathbb R^q$ (nothing truncated on that side), so it integrates over any Borel $A$.

*(b)* Write the restricted step as: draw $\beta\sim\pi(\beta\mid\gamma,y)$ exactly, then $\gamma'\sim N(m(\beta),\Sigma^\star)$ truncated to $\widetilde C_d$. Its density is
$$
q(\gamma'\mid\gamma;\widetilde C_d)=\mathbf 1\{\gamma'\in\widetilde C_d\}\int\frac{\phi_q(\gamma';m(\beta),\Sigma^\star)}{N(\beta)}\pi(\beta\mid\gamma,y)\,d\beta,
\qquad
N(\beta):=\int_{\widetilde C_d}\phi_q(\gamma';m(\beta),\Sigma^\star)\,d\gamma'\le1,
$$
so the integrand is pointwise at least $\phi_q(\gamma';m(\beta),\Sigma^\star)$, whence $q(\gamma'\mid\gamma;\widetilde C_d)\ge q(\gamma'\mid\gamma)\ge\varepsilon_dq_Q(\gamma')$ for $\gamma,\gamma'\in\widetilde C_d$ — truncation only raises the density on the retained region. Integrating over $A\subseteq\widetilde C_d$ and using $Q(A)=Q(\widetilde C_d)Q_{\widetilde C_d}(A)$ gives $q(\gamma,A\mid\widetilde C_d)\ge\varepsilon_dQ(\widetilde C_d)Q_{\widetilde C_d}(A)=\varepsilon\,Q_{\widetilde C_d}(A)$. Finally $Q(\widetilde C_d)>0$: the inscribed ellipsoid $\{(\gamma-\gamma^\star)^\top P_{11}(\gamma-\gamma^\star)\le2d\}\subseteq\widetilde C_d$ (Corollary 15) carries $Q$-mass $\Pr(\chi^2_q\le2d)>0$. $\blacksquare$

**Flat-prior limit (§5.4).** By construction, $\nabla^2\Psi_\tau=\sum_jH_j^\top P_bV_j(\gamma)P_bH_j$ does not involve $\Lambda_\gamma$ at all — the mixing law $\pi(\beta\mid\gamma,y)$ is $\tau$-free — so every result of §5.1–§5.3 whose proof consumes only $\Lambda_Q\succ0$ (supplied at $\tau=0$ by Lemma 7 under (H3a)) and properties of $\pi(\beta\mid\gamma,y)$ transfers to $\tau=0$ verbatim. This is what licenses Lemmas 18–21 below as the $\tau=0$ instances of Lemmas 10–17 rather than independent arguments.

**Proof of Lemma 18.** Lemmas 10–12 read at $\tau=0$: legitimate for the reason above. For the convergence claim, $q_{Q_\tau}\to q_{Q_0}$ uniformly on compacts (explicit Gaussians with $\gamma^\star(\tau)\to\gamma^{\star(0)}$, $P_{11}(\tau)\to P_{11}^{\mathrm{RE}}\succ0$, Lemmas 7–8), and $q_\tau(\gamma'\mid\gamma)\to q_0(\gamma'\mid\gamma)$ locally uniformly by dominated convergence, so $g(\cdot\mid\cdot;\tau)\to g(\cdot\mid\cdot;0)$ uniformly on compacts. The uniform coercivity estimate underlying Lemma 12 holds uniformly in $\tau\in[0,\tau_0]$ as well — its ingredients ($\sigma_{\min}(L_\tau)\ge\sigma_0>0$, boundedness of $a_\tau-\gamma^\star(\tau)$, and the $\tau$-free $\pi(B_v\mid\gamma,y)$) are all $\tau$-uniform for small $\tau$ — so all minimisers lie in one $\tau$-free ball for $\gamma\in K$, giving $|\log\varepsilon(\gamma\mid\tau)-\log\varepsilon(\gamma\mid0)|\le\sup_{K\times B_R}|g(\cdot;\tau)-g(\cdot;0)|\to0$. $\blacksquare$

**Proof of Lemma 19.** Closedness is continuity of $\varepsilon(\cdot\mid0)$ (Lemma 18). Convexity and the compactness sandwich are Lemmas 13–14 read at $\tau=0$: since $\nabla^2\Psi_\tau$ is $\tau$-free, the entire two-sided sandwich $S_\flat\preceq\nabla^2\Psi_\tau\preceq P_{11}^{\mathrm{RE}}$ holds verbatim and uniformly in $\tau\in[0,\tau_0]$, and since $\gamma^\star(\tau)\to\gamma^{\star(0)}$ (Lemma 8), the ellipsoids all sit inside one compact $K_d$. Nonemptiness with $\gamma^{\star(0)}$ interior is Lemma 8 — the one place (H3b) enters. Exhaustion follows from pointwise positivity (Lemma 18). For Hausdorff convergence, $\sup_{K_d}|\Psi_\tau-\Psi_0|\to0$ (from Lemma 18 and $\gamma^\star(\tau)\to\gamma^{\star(0)}$), so $\{\Psi_0\le d-\epsilon_\tau\}\subseteq\widetilde C_d(\tau)\subseteq\{\Psi_0\le d+\epsilon_\tau\}$ with $\epsilon_\tau\to0$; since $\Psi_0$ is convex with $\Psi_0(\gamma^{\star(0)})=0<d$, $\{\Psi_0\le d\}=\overline{\{\Psi_0<d\}}$ (segment points toward the interior minimum stay strictly below $d$ and converge to any boundary point), so both bracketing families converge to $\widetilde C_d(0)$. $\blacksquare$

**Proof of Lemma 20.** The minimum exists and is attained by Weierstrass ($\varepsilon(\cdot\mid0)$ continuous, Lemma 18; $\widetilde C_d(0)$ compact, Lemma 19) and equals $\varepsilon_d(0)$ by the same two-sided argument as Lemma 16, read at $\tau=0$: $\ge$ by definition of the superlevel set, $\le$ by coercivity of $\Psi_0$ plus the intermediate value theorem along a ray from $\gamma^{\star(0)}$. For convergence, take a compact $K\supseteq\{\gamma^\star(\tau):\tau\le\tau_0\}\cup\{\gamma^{\star(0)}\}$ (exists since $\gamma^\star(\tau)\to\gamma^{\star(0)}$, Lemma 8); on $K$, $\varepsilon(\cdot\mid\tau)\to\varepsilon(\cdot\mid0)$ uniformly (Lemma 18) and $\varepsilon(\cdot\mid0)$ is continuous (Lemma 18), so $\varepsilon(\gamma^\star(\tau)\mid\tau)\to\varepsilon(\gamma^{\star(0)}\mid0)$, i.e. $\varepsilon_d(\tau)\to\varepsilon_d(0)$. $\blacksquare$

**Proof of Lemma 21.** Lemma 17(a)–(b) read at $\tau=0$, using Lemma 20 for the constant and Corollary 9 for $Q_0$ — the truncated-density step $N(\beta)\le1$ is $\tau$-free throughout. The refresh-mass bound is Corollary 15 read at $\tau=0$: $\bar\Phi_0\ge0$ gives $\Psi_0(\gamma)\le\frac12\|\gamma-\gamma^{\star(0)}\|^2_{P_{11}^{\mathrm{RE}}}$, so the corresponding ellipsoid is inscribed in $\widetilde C_d(0)$ and carries $Q_0$-mass $\Pr(\chi^2_q\le2d)$ since $\Lambda_{Q_0}=P_{11}^{\mathrm{RE}}$. $\blacksquare$

**Proof of Lemma 22 (Gaussian closure).** Existence, uniqueness, attainment of the minimiser are Lemmas 10–11, which already apply; what remains is to evaluate it. Closure makes $\Lambda_q=\Sigma^{-1}$ constant in $\gamma'$. With $a=M(\gamma)$, $u=\gamma^\star-a$,
$$
g(\gamma')=\tfrac12\log\det(\Sigma^\star\Sigma^{-1})+f(\gamma'),
\qquad
f(\gamma'):=\tfrac12(\gamma'-\gamma^\star)^\top\Lambda_Q(\gamma'-\gamma^\star)-\tfrac12(\gamma'-a)^\top\Lambda_q(\gamma'-a).
$$
Minimising $f$, with $\Delta:=\Lambda_Q-\Lambda_q\succ0$, the stationary point is $\gamma'_\star=\Delta^{-1}(\Lambda_Q\gamma^\star-\Lambda_qa)$, so $\gamma'_\star-\gamma^\star=\Delta^{-1}\Lambda_qu$, $\gamma'_\star-a=\Delta^{-1}\Lambda_Qu$, whence
$$
f(\gamma'_\star)=\tfrac12u^\top\bigl[\Lambda_q\Delta^{-1}\Lambda_Q\Delta^{-1}\Lambda_q-\Lambda_Q\Delta^{-1}\Lambda_q\Delta^{-1}\Lambda_Q\bigr]u=-\tfrac12u^\top\Lambda_q\Delta^{-1}\Lambda_Qu,
$$
expanding both triple products using $\Lambda_Q=\Lambda_q+\Delta$. Since $\Lambda_q(\Lambda_Q-\Lambda_q)^{-1}\Lambda_Q=(\Lambda_q^{-1}-\Lambda_Q^{-1})^{-1}=(\Sigma-\Sigma^\star)^{-1}$ (from $\Lambda_q^{-1}-\Lambda_Q^{-1}=\Lambda_q^{-1}(\Lambda_Q-\Lambda_q)\Lambda_Q^{-1}$), substituting gives the boxed $\mathcal D$. Since $\Sigma\succeq\Sigma^\star$, $\det(\Sigma^\star\Sigma^{-1})\le1$, so $\varepsilon_d\in(0,1]$. $\blacksquare$

---

### A.4 Proof of Theorem 2

**Proof of Lemma 23 (exhaustion).** $\varepsilon(\gamma)>0$ for every $\gamma$ (Lemma 11, or 18 at $\tau=0$), so $\Psi(\gamma)<\infty$ for every $\gamma$, and every $\gamma$ lies in $\widetilde C_d$ once $d\ge\Psi(\gamma)$ — this gives monotonicity and $\bigcup_d\widetilde C_d=\mathbb R^q$. Since $\pi_\gamma$ is a probability measure (H1) and $\widetilde C_d^{\,c}$ is nonincreasing with empty intersection, continuity from above gives $\pi_\gamma(\widetilde C_d^{\,c})\to0$. The explicit envelope $\pi_\gamma(\widetilde C_d^{\,c})\le\pi_\gamma(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d)$ is the complementarity identity ($\Psi\le\frac12\|\cdot-\gamma^\star\|^2_{P_{11}}$, i.e. $\bar\Phi\ge0$) of Lemma 13(3). $\blacksquare$

**Proof of Lemma 24 (restricted kernel and stationarity).** Both blocks of the restricted sweep are drawn from the exact conditional of $\pi(\cdot\mid\tilde C)$ given the other block, and inside the relevant slice of $\tilde C$ by construction — this gives $q(\gamma,C\mid C)=1$ for $\gamma\in C$ (slices are nonempty for $x\in C$ since $x$ itself lies in them), and each block update is therefore a genuine Gibbs update *for the target $\pi(\cdot\mid\tilde C)$*, which leaves its own target invariant; a composition of invariant kernels is invariant. $\blacksquare$ (This is why the kernel is defined as Gibbs for the truncated target rather than as "run the unrestricted sweep and repair excursions on exit" — the latter construction does not in general preserve $\pi(\cdot\mid C)$.)

**Derivation of the rate exponent $\rho$ in Remark 6.1.** Under (H2), $\pi_\gamma$ is log-concave (it is the marginal of the jointly log-concave $\pi(\gamma,\beta\mid y)$, and marginals of log-concave densities are log-concave). A log-concave density with covariance $\Sigma_\pi$ satisfies, for any positive-definite metric $A$, a sub-exponential tail bound for the quadratic form $\|\gamma-\gamma^\star\|^2_A$: writing $\rho:=\lambda_{\max}(\Sigma_\pi^{1/2}A\Sigma_\pi^{1/2})=\lambda_{\max}(\Sigma_\pi A)$, standard concentration for log-concave measures (Borell's lemma applied to the linear functionals diagonalizing $\Sigma_\pi^{1/2}A\Sigma_\pi^{1/2}$) gives $\Pr(\|\gamma-\gamma^\star\|^2_A>2d)\lesssim e^{-d/\rho}$ for $d$ large. Taking $A=P_{11}$ and combining with Lemma 23's envelope, $\pi_\gamma(\widetilde C_d^{\,c})\lesssim e^{-d/\rho}$; solving $e^{-d/\rho}=\delta$ for $d=\rho\log(1/\delta)$ and substituting into $\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)$ gives $\varepsilon\gtrsim\varepsilon(\gamma^\star)\delta^\rho$, since $\varepsilon=\varepsilon_dQ(\widetilde C_d)$ and $Q(\widetilde C_d)\to1$. Under Gaussian closure, $\Sigma_\pi=(P_{11}(I-J))^{-1}$ exactly, giving $\rho=\lambda_{\max}((I-\tilde J)^{-1})=(1-\kappa)^{-1}$. *(The multiplicative constant suppressed by $\lesssim$ is standard but not tracked here; sharpening it to an explicit constant is a numerical, not qualitative, refinement and does not affect the exponent $\rho$.)*

---

## References

- Rosenthal, J. S. (1995). Minorization conditions and convergence rates for Markov
  chain Monte Carlo. *JASA* 90, 558–566.
- `inst/MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` — symmetric GLM minorization, §7
- `inst/LOGIT_STATIC_TAIL_CERTIFICATION.md` §4.5–§4.6, §7.1 — tail mass / Approach A
- `inst/MINORIZATION_GAUSSIAN_REFRESH.md` — $Q$, $\varepsilon_d$
- `inst/PREFLIGHT_model_setup.md` — full rank, estimability (H3)
