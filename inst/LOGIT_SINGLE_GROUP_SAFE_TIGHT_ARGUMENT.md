# Logit Single-Group Safe Region: A Tight Convergence Argument

A condensed rewrite of `LOGIT_SINGLE_GROUP_SAFE_REGION.md`, keeping only
the definitions, lemmas, and proofs needed to reach the final bound.
Companion note to `ELLIPSOID_TV_BOUND.md` (general safe/unsafe-region
machinery) and `SAFE_UNSAFE_TV_DECOMPOSITION.md` (the three-term TV split
this note instantiates). Discussion of motivation, alternatives, and
generalizations lives in `LOGIT_SINGLE_GROUP_SAFE_REGION.md`; this
document states the argument itself.

---

## 1. Setup

**Population parameter** $\theta\sim N(\mu_0,\sigma_0^2)$. **Single-group
random effect** $b\mid\theta \sim N(\theta,\tau^2)$. **Data** (one group,
$m\ge1$ binomial records sharing linear predictor $\eta=b$):

$$
y_i \sim \text{Binomial}(n_i,p(b)), \quad i=1,\dots,m, \quad 0<y_i<n_i,
\qquad p(b) := \frac{e^b}{1+e^b}
$$

($0<y_i<n_i$ rules out separation, so the group MLE is finite.) In
canonical-link form, $\ell(b)=\sum_i\big(y_ib-n_i\log(1+e^b)\big)$, so

$$
\ell'(b) = \sum_i(y_i-n_ip(b)), \qquad
\ell''(b) = -N\,p(b)(1-p(b)) \le 0, \qquad N:=\textstyle\sum_i n_i
$$

**Lemma 1 (bounded score).** $|\ell'(b)|\le L$ for every $b$, where
$L:=\sum_i\max(y_i,n_i-y_i)\le N$.

*Proof.* $y_i-n_ip(b)\in(y_i-n_i,y_i)$ since $p(b)\in(0,1)$, so
$|y_i-n_ip(b)|\le\max(y_i,n_i-y_i)$; sum and apply the triangle
inequality. $\blacksquare$

---

## 2. Safe region

$p(b)(1-p(b))$ is symmetric about $b=0$ and strictly decreasing in
$|b|$. Fix $b^\star>0$ and $\kappa:=N\,p(b^\star)(1-p(b^\star))$. Then

$$
A := \{b: N\,p(b)(1-p(b))\ge\kappa\} = [-b^\star,b^\star], \qquad
A^c = \{b:|b|>b^\star\}
$$

---

## 3. Closed-form bound on $\pi(A^c\mid\theta,y)$

**Lemma 2 (tangent bound).** For any $b_0$ and every $b$,
$\ell(b)\le\ell(b_0)+L|b-b_0|$ (concavity's tangent-line property, plus
Lemma 1).

**Lemma 3 (tilted Gaussian tail).** For $Z\sim N(\theta,\tau^2)$ with
density $\phi$, and $c\ge0$,

$$
\int_t^\infty e^{cz}\phi(z;\theta,\tau^2)\,dz
= \exp\!\Big(c\theta+\tfrac12c^2\tau^2\Big)\,\Phi\!\Big(-\tfrac{t-\theta-c\tau^2}{\tau}\Big)
$$

*Proof.* Complete the square: $e^{cz}\phi(z;\theta,\tau^2) =
e^{c\theta+c^2\tau^2/2}\phi(z;\theta+c\tau^2,\tau^2)$; integrate the
shifted Gaussian from $t$ to $\infty$. $\blacksquare$

**Lemma 4 (normalizer lower bound).** Writing $Z(\theta,y):=\int
e^{\ell(b)}\phi(b;\theta,\tau^2)\,db$,

$$Z(\theta,y) \ge \exp\!\big(\ell(\theta)-L\tau\sqrt{2/\pi}\big)$$

*Proof.* Jensen: $Z\ge e^{E_\phi[\ell(b)]}$. Lemma 1 gives $\ell(b)\ge
\ell(\theta)-L|b-\theta|$, and $E_\phi|b-\theta|=\tau\sqrt{2/\pi}$.
$\blacksquare$

**Proposition 1 (closed-form escape bound).**

$$
\pi(A^c\mid\theta,y) \;\le\; G(\theta) :=
\exp\!\Big(L\tau\sqrt{2/\pi}+\tfrac12L^2\tau^2\Big)
\left[\Phi\!\Big(-\tfrac{b^\star-\theta-L\tau^2}{\tau}\Big)+\Phi\!\Big(-\tfrac{b^\star+\theta-L\tau^2}{\tau}\Big)\right]
$$

*Proof.* By Lemma 2 ($b_0=\theta$), $e^{\ell(b)}\le e^{\ell(\theta)}e^{L|b-\theta|}$.
Split $A^c$ into its two tails and apply Lemma 3 with $c=L$ (and its
mirror image) to bound $\int_{A^c}e^{\ell(b)}\phi(b;\theta,\tau^2)\,db$;
divide by Lemma 4's lower bound on $Z(\theta,y)$ ($e^{\ell(\theta)}$
cancels). $\blacksquare$

This is fully explicit in $(\theta,\tau,b^\star,L)$ — no simulation.

---

## 4. Sampler order and the Foster–Lyapunov drift on $\theta$

**Convention** (`vignette("Chapter-C03")`, Definition 4): $b$ is updated
**first** within a sweep, $\theta$ **second** — matching the package's
actual sampler order. The chain is initialized at $\theta_0=x$ (no
"$b_0$"); sweep $n$ draws $b_n\mid\theta_{n-1}$, then $\theta_n\mid b_n$.
Following `Chapter-C03`, all drift/TV statements below are about
$\theta$, the block updated second.

**(a) Exact Gaussian update.** $\theta$'s likelihood-free conjugacy
(the group likelihood depends on $b$, not $\theta$, directly) gives

$$
\theta_n\mid b_n \sim N\big(wb_n+(1-w)\mu_0,\ \sigma_\theta^2\big), \qquad
w:=\frac{\tau^{-2}}{\sigma_0^{-2}+\tau^{-2}}\in(0,1), \quad
\sigma_\theta^2:=(\sigma_0^{-2}+\tau^{-2})^{-1}
$$

**Lemma 5 (Brascamp–Lieb variance bound).** $\operatorname{Var}(b\mid\theta,y)\le\tau^2$,
since $-\log\pi(b\mid\theta,y)$ has second derivative $\ge\tau^{-2}$
everywhere (log-concavity with a prior-precision floor).

**Lemma 6 (bounded mean-bias).** $|\delta(\theta)|\le L\tau^2$, where
$\delta(\theta):=E[b\mid\theta,y]-\theta$.

*Proof.* With $U(b):=-\ell(b)+(b-\theta)^2/(2\tau^2)$, integration by
parts gives $E_\pi[U'(b)]=0$, i.e. $E[b]-\theta=\tau^2E[\ell'(b)]$; Lemma
1 bounds $|E[\ell'(b)]|\le L$. $\blacksquare$

**Proposition 2 (drift condition).** Let $V(\theta):=(\theta-\mu_0)^2+1$.
Then for every $\theta_{n-1}\in\mathbb R$,

$$
\boxed{\;E\big[V(\theta_n)\mid\theta_{n-1}\big] \le \lambda V(\theta_{n-1}) + C\;}
\qquad
\lambda:=\frac{1+w^2}{2}<1
$$

$$
C := w^2(L\tau^2)^2\Big(1+\tfrac{2w^2}{1-w^2}\Big) + \sigma_\theta^2+w^2\tau^2+1-\lambda
$$

*Proof.* Composing $b_n\mid\theta_{n-1}$ (Lemmas 5–6) then $\theta_n\mid
b_n$ ((a)): $E[\theta_n\mid\theta_{n-1}]-\mu_0 = w(\theta_{n-1}-\mu_0)+\eta(\theta_{n-1})$
with $\eta:=w\delta(\theta_{n-1})$, $|\eta|\le wL\tau^2$. Young's
inequality ($2ab\le\epsilon a^2+b^2/\epsilon$, $\epsilon=(1-w^2)/(2w^2)$)
gives $\big(E[\theta_n\mid\theta_{n-1}]-\mu_0\big)^2 \le
\lambda(\theta_{n-1}-\mu_0)^2 + w^2(L\tau^2)^2(1+\tfrac{2w^2}{1-w^2})$.
For the variance, the law of total variance is exact here since
$\operatorname{Var}(\theta_n\mid b_n)=\sigma_\theta^2$ is constant:
$\operatorname{Var}(\theta_n\mid\theta_{n-1}) = \sigma_\theta^2 +
w^2\operatorname{Var}(b_n\mid\theta_{n-1}) \le \sigma_\theta^2+w^2\tau^2$
(Lemma 5). Add $E[(\theta_n-\mu_0)^2\mid\theta_{n-1}] =
(E[\theta_n\mid\theta_{n-1}]-\mu_0)^2 + \operatorname{Var}(\theta_n\mid\theta_{n-1})$
and the $+1$ in $V$. $\blacksquare$

---

## 5. From the $\theta$-drift to $P^n(x,A^c)$ and to $\pi(A^c)$

Let $K^\star := \min\big((b^\star-\mu_0)^2,(b^\star+\mu_0)^2\big)+1 = \inf_{b\in A^c}V(b)$,
and write $g(t):=\pi(A^c\mid\theta=t,y)\in[0,1]$ for the exact one-step
escape probability (so $b_n\mid\theta_{n-1}$ drawn from $\pi(\cdot\mid\theta_{n-1},y)$
gives $E[\mathbf 1_{A^c}(b_n)\mid\theta_{n-1}]=g(\theta_{n-1})$ pointwise).

**Lemma 7 (pointwise escape bound).** For every $t\in\mathbb R$,

$$
g(t) \;\le\; \Psi(t) := \frac{2V(t) + 2(L\tau^2)^2+\tau^2+1}{K^\star}
$$

*Proof.* By Lemmas 5–6, $E[(b_n-\mu_0)^2\mid\theta_{n-1}=t] =
(t-\mu_0+\delta(t))^2+\operatorname{Var}(b_n\mid\theta_{n-1}=t)
\le 2(t-\mu_0)^2+2(L\tau^2)^2+\tau^2 \le 2V(t)+2(L\tau^2)^2+\tau^2$ (using
$(t-\mu_0)^2\le V(t)$). Add $1$ and apply Markov's inequality to
$V(b_n)\ge0$ on the event $b_n\in A^c$, where $V(b)\ge K^\star$. $\blacksquare$

This one deterministic, $n$-free inequality now feeds two different
averages of $\theta_{n-1}$:

**Proposition 3** (started at $x$, iterated drift). Since
$P^n(x,A^c)=E[g(\theta_{n-1})\mid\theta_0=x]\le E[\Psi(\theta_{n-1})\mid\theta_0=x]$
and Proposition 2's drift iterates to $E[V(\theta_{n-1})\mid\theta_0=x]\le\lambda^{n-1}V(x)+C/(1-\lambda)$,

$$
P^n(x,A^c) \;\le\; \underbrace{\frac{2\lambda^{n-1}V(x) + 2C/(1-\lambda) + 2(L\tau^2)^2+\tau^2+1}{K^\star}}_{=:U_n(x)}
$$

**Proposition 3$'$** (stationary average). Since $\theta_{n-1}\sim\pi_\theta$
makes $b_n\sim\pi(\cdot\mid y)$ exactly (the $b$-marginal of the target),
$\pi(A^c) = E_{\pi_\theta}[g(\theta)] \le E_{\pi_\theta}[\Psi(\theta)]$, and
averaging Proposition 2's drift against its own stationary law gives the
standard identity $E_{\pi_\theta}[V]\le C/(1-\lambda)$ (Meyn and Tweedie,
1993, Ch. 14), so

$$
\pi(A^c) \;\le\; \underbrace{\frac{2C/(1-\lambda) + 2(L\tau^2)^2+\tau^2+1}{K^\star}}_{=:U_\infty} \;=\; \lim_{n\to\infty}U_n(x)
$$

Both $U_n(x)$ and $U_\infty$ shrink to $0$ as $b^\star\to\infty$ (since
$K^\star\to\infty$) — this closes Open Item 3 below (the fully marginal
bound on $\pi(A^c)$).

---

## 6. Bounding the gap $\big|\pi(A^c)-P^n(x,A^c)\big|$

The naive route bounds
$|\int g\,d\pi_\theta-\int g\,dP_\theta^{(n-1)}(x,\cdot)|$ by
$\|g\|_\infty\|P_\theta^{(n-1)}(x,\cdot)-\pi_\theta(\cdot)\|_{TV}$ — but
$g(t)\to1$ as $t\to\pm\infty$ (Proposition 1's tilted-Gaussian bound
saturates, it does not decay), so $\|g\|_\infty=1$ and this is exactly
the whole-support flaw again. Splitting $\theta$ into a core/tail (as in
an earlier draft of this argument) does not fix it either: bounding the
two pieces of the *gap* separately reduces to bounding $\pi(A^c)$ and
$P^n(x,A^c)$ on $B^c$ *separately*, which reintroduces a term
($\lambda^{n-1}V(x)+2C/(1-\lambda)$ over $K_B^\star$) that does not
vanish as $n\to\infty$ for any fixed $B$. What is required is a bound on
the **difference** that is small for two independent reasons — either
because the chains have already coupled ($n$ large) or because $A^c$
itself is rare ($b^\star$ large) — not a bound that needs both at once.

**Coupling construction.** Let $(\theta_n^x)_{n\ge0}$ be the
$\theta$-marginal chain started at $\theta_0^x=x$, and let
$(\theta_n^\pi)_{n\ge0}$ be a second copy started from $\theta_0^\pi\sim\pi_\theta$
(so $\theta_n^\pi\sim\pi_\theta$ for every $n$, by stationarity). Couple
them by the standard minorization construction underlying Theorem 1
below: run independently until both lie in the small set $C_d$
simultaneously, then merge with probability $\varepsilon$ at each such
joint visit, using the common measure $Q$ from the minorization; let $T$
be the resulting (random) coupling time, i.e. the first $n$ with
$\theta_n^x=\theta_n^\pi$. Once merged, run the two copies of the full
chain identically forever, so $\theta_{n'}^x=\theta_{n'}^\pi$ for all
$n'\ge T$. In particular, on $\{T\le n-1\}$, draw $b_n^x$ and $b_n^\pi$
from a single shared random draw of $\pi(\cdot\mid\theta_{n-1}^x,y)=\pi(\cdot\mid\theta_{n-1}^\pi,y)$
(valid since the two conditionals coincide), giving $b_n^x=b_n^\pi$ on
that event.

**Theorem 1 (Rosenthal, 1995, Theorem 12 / Theorem 5).** Fix any
$d>2C/(1-\lambda)$ and $C_d:=\{\theta:V(\theta)\le d\}$. If the one-step
$\theta$-kernel $q_\theta$ satisfies a minorization
$q_\theta(\theta,\cdot)\ge\varepsilon Q(\cdot)$ for $\theta\in C_d$
(existence: $q_\theta(\theta'\mid\theta)$ is jointly continuous and
strictly positive, hence bounded below on a compact product
$C_d\times[\theta_-',\theta_+']$ — the value of $\varepsilon$ itself is
not computed here), then for any $0<r<1$ the coupling time of the
construction above satisfies

$$
P(T>k) \;\le\;
(1-\varepsilon)^{rk} + \alpha^{-k}(\alpha\Lambda)^{rk}\Big[1+\frac{C}{1-\lambda}+V(x)\Big] \;=:\; R_k,
$$

$$
\alpha^{-1}:=\frac{1+2C+\lambda d}{1+d}<1, \qquad \Lambda:=1+2(\lambda d+C)
$$

(Rosenthal's Theorem 12 quotes the $\|P_\theta^{(k)}(x,\cdot)-\pi_\theta(\cdot)\|_{TV}\le R_k$
consequence of this; the bound on $P(T>k)$ itself, which is what his
proof — via renewal/regeneration on $C_d$ — actually establishes, is the
form needed here.)

**Proposition 4 (revised, Cauchy–Schwarz gap bound).**

$$
\big|\pi(A^c)-P^n(x,A^c)\big| \;\le\;
\Big(\sqrt{P^n(x,A^c)}+\sqrt{\pi(A^c)}\Big)\sqrt{P(T>n-1)}
\;\le\;
\Big(\sqrt{U_n(x)}+\sqrt{U_\infty}\Big)\sqrt{R_{n-1}}
$$

*Proof.* Using the coupling above and $b_n^x=b_n^\pi$ on $\{T\le n-1\}$,

$$
\pi(A^c)-P^n(x,A^c) = E\big[\mathbf 1_{A^c}(b_n^\pi)-\mathbf 1_{A^c}(b_n^x)\big]
= E\Big[\big(\mathbf 1_{A^c}(b_n^\pi)-\mathbf 1_{A^c}(b_n^x)\big)\mathbf 1_{T>n-1}\Big]
$$

since the summand is $0$ on $\{T\le n-1\}$. As both indicators are
nonnegative, $\big|\pi(A^c)-P^n(x,A^c)\big| \le E[\mathbf 1_{A^c}(b_n^x)\mathbf 1_{T>n-1}]+E[\mathbf 1_{A^c}(b_n^\pi)\mathbf 1_{T>n-1}]$.
Conditioning the first term on $\theta_{n-1}^x$ (whose $\sigma$-field
determines $\{T\le n-1\}$ up to independent coupling randomness, so
$b_n^x\mid\theta_{n-1}^x$ is unaffected) gives
$E[\mathbf 1_{A^c}(b_n^x)\mathbf 1_{T>n-1}] = E[g(\theta_{n-1}^x)\mathbf 1_{T>n-1}]$.
By Cauchy–Schwarz and $g^2\le g$ (as $g\in[0,1]$),

$$
E[g(\theta_{n-1}^x)\mathbf 1_{T>n-1}] \le \sqrt{E[g(\theta_{n-1}^x)^2]}\sqrt{P(T>n-1)}
\le \sqrt{E[g(\theta_{n-1}^x)]}\sqrt{P(T>n-1)} = \sqrt{P^n(x,A^c)}\sqrt{P(T>n-1)}
$$

and likewise $E[\mathbf 1_{A^c}(b_n^\pi)\mathbf 1_{T>n-1}]\le\sqrt{\pi(A^c)}\sqrt{P(T>n-1)}$.
Sum the two, then substitute Proposition 3/3$'$ and Theorem 1. $\blacksquare$

**Why this is both tail-aware and vanishing.** $R_{n-1}\to0$ as
$n\to\infty$ for fixed $b^\star$ (Theorem 1), while $\sqrt{U_n(x)}+\sqrt{U_\infty}$
stays bounded (it converges to $2\sqrt{U_\infty}\le2$), so the product
$\to0$ — the bound vanishes in $n$ with **no** additive floor, unlike the
core/tail split. Separately, $U_n(x)$ and $U_\infty\to0$ as
$b^\star\to\infty$ for fixed $n$ (Propositions 3/3$'$), so the bound also
shrinks as the safe region grows — it does not need $n\to\infty$ to see
$A^c$ shrink. The two properties come from two different factors in a
product, not from one term trying to do both jobs at once.

---

## 7. Assembled bound

$$
\big\|P^n(x,\cdot)-\pi(\cdot)\big\|_{TV} \;\le\;
\underbrace{\big\|P^n_A(x,\cdot)-\pi_A(\cdot)\big\|_{TV}}_{\text{safe/unsafe heuristic §7.1}}
\;+\; 4\,P^n(x,A^c) \;+\; 2\big|\pi(A^c)-P^n(x,A^c)\big|
$$

with the second term from Proposition 3, the third from the combined
Proposition 4/Theorem 1 bound above (using $\pi(A^c)\le U_\infty$ from
Proposition 3$'$ for the on-$A$ term's normalization correction as well).
Every constant is explicit given $(\theta_0,\tau,\mu_0,\sigma_0,b^\star,y,d,r)$
except $\varepsilon$.

---

## 8. Scope and open items

1. $\varepsilon$ (Theorem 1's minorization constant) is not computed
   numerically — the only unclosed input in the whole bound.
2. $r$ (Theorem 1) is a free parameter, not optimized.
3. ~~The fully marginal $\pi(A^c\mid y)$ is not derived.~~ Resolved by
   Proposition 3$'$ (§5): $\pi(A^c)\le U_\infty$, the $n\to\infty$ limit
   of Proposition 3's bound, via the same pointwise Lemma 7 averaged
   against $\pi_\theta$ instead of iterated from $x$.
4. Multiple groups / a random-effect design matrix with $p_j>1$ columns:
   $\theta$ stays scalar, so Proposition 2's drift should extend with
   Lemmas 5–6 applied group-by-group, but the aggregation and the
   resulting safe-region geometry (§2) are not worked out here. See
   `LOGIT_SINGLE_GROUP_SAFE_REGION.md` §6 for the fuller discussion.

---

## References

- `ELLIPSOID_TV_BOUND.md` — general safe/unsafe decomposition and its
  logit specialization, of which §2–§3 above are the scalar special case.
- `SAFE_UNSAFE_TV_DECOMPOSITION.md` — the three-term TV split instantiated
  in §7, and the computable on-$A$ heuristic (§7.1) its first term invokes.
- `vignette("Chapter-C03")` — Definition 4's update-order convention,
  used in §4 to anchor the drift on $\theta$.
- Rosenthal, J.S. (1995). *Minorization conditions and convergence rates
  for Markov chain Monte Carlo.* JASA 90, 558–566 — Theorem 12, applied
  in §6.
- Meyn, S.P. and Tweedie, R.L. (1993). *Markov Chains and Stochastic
  Stability* — the stationary-moment bound used in Proposition 3$'$.
- Brascamp, H.J. and Lieb, E.H. (1976). *On extensions of the
  Brascamp–Lieb–Luttinger inequality...* — Lemma 5's variance bound.
