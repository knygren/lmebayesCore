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

*(Content to follow. Profile $\varepsilon(\gamma):=\inf_{\gamma'}q(\gamma'\mid\gamma)/q_Q(\gamma')$,
a function of the source state $\gamma$ alone.)*

### 5.1 The minorization profile $\varepsilon(\gamma)$

- Well-posedness: unique minimiser, attained, strictly positive
- Continuity

### 5.2 Convexity, coercivity, and compactness of $\widetilde C_d$

- $\log\varepsilon$ concave (equivalently, the deficiency gap $\Psi$ convex)
- Coercivity of $\Psi$; the certified set trapped between concentric ellipsoids
- Compactness from (H3a) alone — prior-free

### 5.3 The attained constant $\varepsilon_d$

- $\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)$ as the attained minimum on
  $\widetilde C_d$ (Weierstrass, via §5.2)
- The minorization condition on $\widetilde C_d$, restricted and unrestricted forms

### 5.4 Flat-prior limit

- The profile, the certified set, and the constant as $\Lambda_\gamma\downarrow0$
- $\varepsilon_d(\Lambda_\gamma)\to\varepsilon_d(0)$ — the constant does not collapse
  as the prior becomes diffuse

### 5.5 Closed form under Gaussian closure

- Explicit $\mathcal D$ and $\varepsilon_d\ge(1+\kappa)^{-q/2}e^{-d}$ for LMMs
- What is lost for GLMMs (the kernel is a mixture, not Gaussian)

---

## 6. Proof of Theorem 2

*(Content to follow.)*

### 6.1 Assembly of §4–§5

### 6.2 Lifting to the joint $(\gamma,\beta)$ chain

---

## 7. Scope and Extensions

*(Content to follow.)*

### 7.1 Symmetric case: sharper constants

### 7.2 What is and is not certified

### 7.3 Open problems

---

## Appendix — Proofs

### A.1 Proof of §3.2 (propriety of the flat-prior posterior) and §3.3 (mean map)

### A.2 Proofs for §4

### A.3 Proofs for §5

### A.4 Proof of Theorem 2

---

## References

- Rosenthal, J. S. (1995). Minorization conditions and convergence rates for Markov
  chain Monte Carlo. *JASA* 90, 558–566.
- `inst/MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` — symmetric GLM minorization, §7
- `inst/LOGIT_STATIC_TAIL_CERTIFICATION.md` §4.5–§4.6, §7.1 — tail mass / Approach A
- `inst/MINORIZATION_GAUSSIAN_REFRESH.md` — $Q$, $\varepsilon_d$
- `inst/PREFLIGHT_model_setup.md` — full rank, estimability (H3)
