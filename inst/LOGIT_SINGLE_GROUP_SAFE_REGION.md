# A Floor-Based Safe Region for a Single-Group Logit Random Effect

Companion note to `ELLIPSOID_TV_BOUND.md` (the general safe/unsafe-region
machinery, §§1–4, and its logit specialization, §6) and to
`vignette("Chapter-C05")` (the model-free Theorem 1 form of the three-term
TV bound). This note works a **concrete special case** all the way through
in closed form, as a template for later generalizing back into
`ELLIPSOID_TV_BOUND.md`'s §6.

## 0. Why this case, and how it differs from §6 of `ELLIPSOID_TV_BOUND.md`

`ELLIPSOID_TV_BOUND.md` §6 defines the logit safe region $A$ via a **pilot
eigenvalue ceiling** $M_j^{\text{pilot}}$ — safe means the *data* curvature
is no worse than what a pilot run already observed. That split is about
sharpening the *mixing rate on $A$* relative to the theoretical worst case
$M_j^{\text{theory}}=\tau_j^{-2}+\tfrac14\lambda_{\max}(D_j'D_j)$; it is
silent on $A^c$ beyond "the loose global rate still applies there," and
§6.2 explicitly gives up on a closed form for $\pi(A^c)$, resorting to an
empirical estimate (`ELLIPSOID_TV_BOUND.md` §6.2 *Status* remark).

Here the split runs the other way: $A$ is a **curvature-floor** region
(exactly Definition 1's original framing, mirrored for probit/Poisson in
§7–§8) — safe means the data curvature is *at least* a chosen level,
unsafe means it has decayed *below* that level. For logit's curvature
function $c(\eta)=p(\eta)(1-p(\eta))$, which is unimodal and *decreasing*
in $|\eta|$, this produces something §6 never has to deal with: a
**bounded** safe region (an interval), not a half-space or a
sublevel/superlevel set defined by a matrix eigenvalue. Because the
random effect here is scalar (one group, no covariates on the group
term), that interval is *exactly computable*, and so — unlike the general
multi-dimensional block case flagged as unresolved in §6.4 item 2 — is a
genuinely closed-form bound on $\pi(A^c)$. That is the payoff of working
the special case first.

---

## 1. Model and notation

**Population parameter** $\theta \in \mathbb{R}$ (block 1). **Single-group
random effect** $b \in \mathbb{R}$ (block 2), with prior centered at the
population parameter:

$$b \mid \theta \;\sim\; N(\theta, \tau^2)$$

**Data** (one group, $m \ge 1$ conditionally independent binomial records,
covariate-free so every record shares the linear predictor $\eta = b$):

$$
y_i \sim \text{Binomial}(n_i, p(b)), \quad i = 1,\dots,m, \quad n_i > 1,
\quad 0 < y_i < n_i, \qquad p(b) := \operatorname{expit}(b) = \frac{e^b}{1+e^b}
$$

The condition $0 < y_i < n_i$ for every $i$ is this model's instance of
Assumption 1 (`ELLIPSOID_TV_BOUND.md` §1): it rules out complete/quasi-
complete separation, so the group MLE $\hat b$ is finite and the machinery
below is not vacuous. Write $N := \sum_i n_i$ for the total trial count in
the group.

In canonical-link form, $\ell_i(b) = y_i b - n_i\log(1+e^b)$, so the
group's log-likelihood and its first two derivatives are

$$
\ell(b) = \sum_i \ell_i(b), \qquad
\ell'(b) = \sum_i \big(y_i - n_i\, p(b)\big), \qquad
\ell''(b) = -\sum_i n_i\, p(b)(1-p(b)) = -N\, p(b)(1-p(b))
$$

$\ell$ is concave ($\ell''\le0$ everywhere — the standard canonical-link
GLM fact). Writing $c(b) := -\ell''(b) = N\,p(b)(1-p(b))$ for the data
curvature (this is exactly `ELLIPSOID_TV_BOUND.md`'s $c(\eta)$ summed over
the group's $m$ records, since every record shares $\eta=b$), the
conditional negative-log-posterior Hessian for $b \mid \theta, y$ is

$$H(b) = \tau^{-2} + N\, p(b)(1-p(b))$$

and, as in §6.1 of the umbrella note, the **global sandwich**

$$\tau^{-2} \;\le\; H(b) \;\le\; \tau^{-2} + \tfrac{N}{4} \qquad \text{(for every finite } b\text{)}$$

holds unconditionally — the floor from the prior, the ceiling from
$p(1-p)\le\tfrac14$, attained only at $b=0$ ($p=1/2$).

---

## 2. The floor-based safe region

$p(b)(1-p(b)) = e^b/(1+e^b)^2$ is symmetric about $b=0$
($p(-b)(1-p(-b))=p(b)(1-p(b))$, since $p(-b)=1-p(b)$) and strictly
decreasing in $|b|$ (its derivative is $p(1-p)(1-2p)$, positive for $b<0$,
negative for $b>0$). Fix a threshold $b^\star > 0$ — the largest deviation
from $b=0$ (equivalently from $p=1/2$) you are willing to call "safe" —
and define the curvature floor at that threshold,

$$\kappa := N\, p(b^\star)\big(1-p(b^\star)\big)$$

Because of the symmetric/unimodal shape above, the safe region reduces to
an explicit, closed-form **interval**:

$$
A := \big\{ b : N\, p(b)(1-p(b)) \ge \kappa \big\} \;=\; [-b^\star, b^\star],
\qquad A^c = \{ b : |b| > b^\star \}
$$

On $A$: $\tau^{-2}+\kappa \le H(b) \le \tau^{-2}+N/4$ (the ceiling
is unchanged, since $b=0\in A$ always attains it; the floor is sharpened
from the bare prior precision to $\tau^{-2}+\kappa$). On $A^c$: only
$\tau^{-2}\le H(b) < \tau^{-2}+\kappa$ is guaranteed — curvature can decay
toward (but never below) the prior floor as $|b|\to\infty$, at the rate
$p(1-p)\sim e^{-|b|}$.

This is Definition 1 of `ELLIPSOID_TV_BOUND.md`, specialized: unlike §7's
probit or §8's Poisson (where the "good" region is an unbounded
half-space or exterior set because curvature is only bounded *below*),
here $A$ is a bounded interval around the point of maximal curvature,
because logit's curvature is bounded *above and below*, and is unimodal.

---

## 3. Bounding $\pi(A^c \mid \theta, y)$ in closed form

The key structural fact is that the logit score is globally bounded — the
same fact §6.1 uses for the curvature floor, one derivative up.

**Lemma L1 (bounded score).** For every $b \in \mathbb{R}$,

$$|\ell'(b)| \;\le\; L, \qquad \text{where} \qquad L := \sum_i \max(y_i,\, n_i - y_i) \;\le\; N$$

*Proof.* $y_i - n_i p(b) \in (y_i-n_i,\,y_i)$ since $p(b)\in(0,1)$
strictly, so $|y_i-n_ip(b)| \le \max(y_i,n_i-y_i)$ for every $b$ (a linear
function of $p\in[0,1]$ is maximized in absolute value at an endpoint).
Sum over $i$ and apply the triangle inequality. $\blacksquare$

**Lemma L2 (concavity gives a global tangent-line bound).** For any fixed
$b_0 \in \mathbb{R}$ and every $b$,

$$\ell(b) \;\le\; \ell(b_0) + \ell'(b_0)(b-b_0) \;\le\; \ell(b_0) + L|b-b_0|$$

*Proof.* First inequality: concavity ($\ell''\le0$), the standard
supporting-tangent-line property. Second: Lemma L1 applied at $b_0$, valid
for any $b_0$, any sign of $b-b_0$. $\blacksquare$

**Lemma L3 (exponentially tilted Gaussian tail — exact).** For
$Z \sim N(\theta,\tau^2)$ with density $\phi(\cdot;\theta,\tau^2)$, and any
$c \ge 0$,

$$
\int_t^\infty e^{cz}\,\phi(z;\theta,\tau^2)\,dz
\;=\; \exp\!\Big(c\theta + \tfrac12 c^2\tau^2\Big)\;
\Phi\!\Big(-\tfrac{t-\theta-c\tau^2}{\tau}\Big)
$$

*Proof.* Complete the square: $cz - \dfrac{(z-\theta)^2}{2\tau^2} =
c\theta + \tfrac12c^2\tau^2 - \dfrac{(z-\theta-c\tau^2)^2}{2\tau^2}$, so
$e^{cz}\phi(z;\theta,\tau^2) = e^{c\theta+c^2\tau^2/2}\,\phi(z;\theta+c\tau^2,\tau^2)$
— an exponentially tilted Gaussian density is itself Gaussian, shifted by
$c\tau^2$, same variance. Integrate the shifted density from $t$ to
$\infty$. $\blacksquare$

**Lemma L4 (Jensen lower bound on the normalizer).** Writing
$Z(\theta,y) := \displaystyle\int_{\mathbb{R}} e^{\ell(b)}\phi(b;\theta,\tau^2)\,db$
for the normalizing constant of $\pi(b\mid\theta,y) \propto
e^{\ell(b)}\phi(b;\theta,\tau^2)$,

$$Z(\theta,y) \;\ge\; \exp\!\big(\ell(\theta) - L\tau\sqrt{2/\pi}\big)$$

*Proof.* Jensen ($\exp$ convex): $Z = E_\phi[e^{\ell(b)}] \ge
e^{E_\phi[\ell(b)]}$. Lemma L1 gives the two-sided bound $\ell(b) \ge
\ell(\theta) - L|b-\theta|$ (not just the one-sided tangent bound of L2),
so $E_\phi[\ell(b)] \ge \ell(\theta) - L\,E_\phi|b-\theta|$; the mean
absolute deviation of $N(\theta,\tau^2)$ is $\tau\sqrt{2/\pi}$.
$\blacksquare$

**Proposition 1 (closed-form bound on the target escape probability).**

$$
\pi(A^c \mid \theta, y)
\;\le\;
\exp\!\Big(L\tau\sqrt{2/\pi} + \tfrac12 L^2\tau^2\Big)
\left[
\Phi\!\Big(-\tfrac{b^\star-\theta-L\tau^2}{\tau}\Big)
+
\Phi\!\Big(-\tfrac{b^\star+\theta-L\tau^2}{\tau}\Big)
\right]
$$

*Proof.* Apply Lemma L2 with $b_0=\theta$: $e^{\ell(b)} \le
e^{\ell(\theta)}e^{L|b-\theta|}$ for every $b$. Split $A^c$ into its two
tails and apply Lemma L3 with $c=L$ (upper tail, $b>b^\star$) and its
mirror image (lower tail, $b<-b^\star$) to bound

$$
\int_{A^c}e^{\ell(b)}\phi(b;\theta,\tau^2)\,db
\;\le\;
e^{\ell(\theta)}\,e^{L^2\tau^2/2}
\left[
\Phi\!\Big(-\tfrac{b^\star-\theta-L\tau^2}{\tau}\Big)
+
\Phi\!\Big(-\tfrac{b^\star+\theta-L\tau^2}{\tau}\Big)
\right].
$$

Divide by the Lemma L4 lower bound on $Z(\theta,y)$; the $e^{\ell(\theta)}$
factors cancel, leaving the $L\tau\sqrt{2/\pi}$ term from L4 in the
exponent. $\blacksquare$

This is fully explicit given $(\theta,\tau,b^\star)$ and the data-dependent
constant $L=\sum_i\max(y_i,n_i-y_i)$ — no simulation, no numerical
integration beyond evaluating $\Phi$. It answers, for this scalar case,
exactly the closed-form question `ELLIPSOID_TV_BOUND.md` §6.2 and §6.4
item 2 flag as open for the general (matrix-eigenvalue) logit block: there,
the safe-region indicator is a nonlinear eigenvalue functional with no
known closed-form tail; here, because the block is scalar, the "eigenvalue
functional" is just $p(b)(1-p(b))$ itself, and its super/sub-level sets
are the intervals worked out in §2 above.

**Remark (this is conditional on $\theta$, matching the note's existing
convention).** Exactly as `ELLIPSOID_TV_BOUND.md`'s Definition 1 and Lemma
2 state block-$j$'s safe region and tail bound conditional on "the rest"
($\pi(\beta_j\mid\text{rest})$), Proposition 1 bounds $\pi(A^c\mid\theta,y)$
for every fixed $\theta$. Recovering a fully marginal
$\pi(A^c\mid y) = E_{\theta\mid y}[\pi(A^c\mid\theta,y)]$ needs an outer
bound on the posterior of $\theta$ itself — exactly the "(analogous terms
for population/RE blocks)" placeholder in Claim 2 (§4) of the umbrella
note. If $\theta$ has a proper Gaussian (or flat, in the limit) prior,
$\theta\mid b,y$ is *exactly* Gaussian (the group likelihood depends on
$b$, not $\theta$, directly — $\theta$ only enters through the
conjugate-Gaussian prior link), so this outer step should itself be
closed-form; it is not carried out here (§7, open item).

---

## 4. Bounding the posterior-vs-sampler gap directly

**4.0 Why not the whole-space TV distance.** A tempting shortcut is
$|\pi(A^c)-P^n(x,A^c)| \le \|P^n(x,\cdot)-\pi(\cdot)\|_{TV} \le
M\rho_{\text{glob}}^{\,n}$, using the whole-space strong log-concavity of
§1 (`ELLIPSOID_TV_BOUND.md` §6.2's uniform-ergodicity argument). This is
*correct* but wasteful: $\|P^n(x,\cdot)-\pi(\cdot)\|_{TV}$ is a supremum
over **every** measurable set, and the resulting bound carries no memory
of the fact that $A^c$ is specifically a *far-tail* event — the bound
would be exactly the same size whether $b^\star$ is $2$ or $20$. What we
actually want is a bound that gets *smaller* as $A^c$ moves further into
the tail, and that targets the gap `\pi(A^c)-P^n(x,A^c)` itself rather
than bounding $\pi(A^c)$ and $P^n(x,A^c)$ separately and subtracting
(subtracting two separately-bounded quantities cannot exploit any
cancellation between them, and is no better than the whole-space route
once each is bounded that loosely). §4.5 below replaces the earlier
abstract Meyn–Tweedie-$f$-norm route with a sharper, fully explicit
argument specific to a two-block Gibbs sampler: the gap is bounded
directly by the total variation distance of a single scalar chain
(§4.5), for which Rosenthal's (1995) Theorem 12 supplies a concrete,
citable formula (§4.6) rather than a free constant $R$.

**4.1 Update order, and which variable is the recursive anchor.** As in
`vignette("Chapter-C03")` (Definition 4 there): $b$ is the block updated
**first** within a sweep, $\theta$ is the block updated **second**. The
package's actual sampler order — group parameters, then population
parameters — matches this exactly, with $b \leftrightarrow$
`Chapter-C03`'s $X_2$ and $\theta \leftrightarrow X_1$. Following
Definition 4's own indexing, the chain is initialized by specifying
$\theta_0 = x$ — there is no "$b_0$", exactly as `Chapter-C03` has no
$x_2^{(0)}$ — and sweep $n = 1, 2, \dots$ draws, in order,

$$
b_n \mid \theta_{n-1} \qquad\text{(first update of sweep $n$)}, \qquad
\theta_n \mid b_n \qquad\text{(second update of sweep $n$)}.
$$

`Chapter-C03` states every drift/mean/variance/TV result "about $X_1$,
the block updated second" (its remark immediately after Definition 4);
the same reasoning applies here, for the same structural reason —
$\theta$, not $b$, has a state defined at sweep $0$ and a one-step
transition whose *outer* half (the conjugate-Gaussian step below) is
exact in closed form. §4.3 below therefore states the Foster–Lyapunov
drift for $V(\theta)$, not $V(b)$; §4.4 then converts a tail statement
about $\theta_{n-1}$ into one about $b_n$ (one sweep-step ahead) using
the *other* building block; §4.5 pins down the posterior-vs-sampler gap
itself, which is about $b$, via a chain of exact conditional identities
that reduce it to a statement about the $\theta$-marginal chain alone.

**4.2 Two exact/explicit building blocks.**

*(a) The population block is exact conjugate Gaussian.* Give $\theta$ a
proper prior $\theta \sim N(\mu_0,\sigma_0^2)$ (the group likelihood never
depends on $\theta$ directly — it only enters through the Gaussian prior
link on $b$ — so this update is exact, no GLM nonlinearity). For the
*second* update of sweep $n$,

$$
\theta_n \mid b_n \;\sim\; N\big(w\,b_n + (1-w)\mu_0,\; \sigma_\theta^2\big),
\qquad
w := \frac{\tau^{-2}}{\sigma_0^{-2}+\tau^{-2}} \in (0,1), \quad
\sigma_\theta^2 := \big(\sigma_0^{-2}+\tau^{-2}\big)^{-1}
$$

*(b) The random-effect block's conditional mean and variance are both
bounded, uniformly over $\theta$.* These are properties of
$\pi(b\mid\theta,y)$ for any fixed $\theta$ — no sampler order enters —
so they apply directly to the *first* update of sweep $n$, $b_n \mid
\theta_{n-1}$. Lemma L5 is a Brascamp–Lieb variance bound; Lemma L6 plays
the same role for the *mean*.

**Lemma L5 (Brascamp–Lieb variance bound).** Since
$-\log\pi(b\mid\theta,y) = -\ell(b) + (b-\theta)^2/(2\tau^2) + \text{const}$
has second derivative $H(b) \ge \tau^{-2}$ for every $b$ (§1), the
Brascamp–Lieb inequality for log-concave densities gives

$$\operatorname{Var}(b \mid \theta, y) \;\le\; \tau^2$$

**Lemma L6 (bounded mean-bias, via Stein's identity).** Write
$\delta(\theta) := E[b\mid\theta,y] - \theta$. Then $|\delta(\theta)| \le
L\tau^2$ for every $\theta$.

*Proof.* With $U(b) := -\ell(b)+(b-\theta)^2/(2\tau^2)$ (so $\pi(b\mid
\theta,y)\propto e^{-U(b)}$, vanishing as $b\to\pm\infty$ since the
quadratic term dominates the at-most-linearly-growing $-\ell(b)$),
integration by parts gives the exact identity $E_\pi[U'(b)] =
-\int_{\mathbb R} \tfrac{d}{db}e^{-U(b)}\,db\big/Z = 0$. Since $U'(b) =
-\ell'(b) + (b-\theta)/\tau^2$, this reads $E[b]-\theta = \tau^2\,
E[\ell'(b)]$, and $|\ell'(b)|\le L$ pointwise (Lemma L1) gives
$|E[\ell'(b)]|\le L$. $\blacksquare$

**4.3 An explicit Foster–Lyapunov drift condition, on $\theta$.** Let
$V(\theta) := (\theta-\mu_0)^2+1$. Composing (b) then (a), **in that
order** — $b_n\mid\theta_{n-1}$ first, $\theta_n\mid b_n$ second, per
§4.1:

$$
E[\theta_n\mid\theta_{n-1}] - \mu_0 \;=\; w(\theta_{n-1}-\mu_0) + \eta(\theta_{n-1}), \qquad
|\eta(\theta_{n-1})| \le w L\tau^2 \le L\tau^2
$$

*Derivation.* $E[\theta_n\mid\theta_{n-1}] = E_{b_n\mid\theta_{n-1}}\big[w b_n+(1-w)\mu_0\big]
= w\,E[b_n\mid\theta_{n-1}]+(1-w)\mu_0 = w\big(\theta_{n-1}+\delta(\theta_{n-1})\big)+(1-w)\mu_0$
by (a) and Lemma L6, giving $\eta(\theta_{n-1}) := w\,\delta(\theta_{n-1})$.
Squaring and applying $2ab\le\epsilon a^2+b^2/\epsilon$ with
$\epsilon = (1-w^2)/(2w^2)$ (as before, chosen so the resulting quadratic
coefficient stays below 1 for *every* $w\in(0,1)$):

$$
\big(E[\theta_n\mid\theta_{n-1}]-\mu_0\big)^2 \;\le\; \frac{1+w^2}{2}(\theta_{n-1}-\mu_0)^2
\;+\; w^2(L\tau^2)^2\Big(1+\tfrac{2w^2}{1-w^2}\Big)
$$

For the variance, the law of total variance gives an **exact** (not
crudely split) decomposition, because $\operatorname{Var}(\theta_n\mid b_n)=\sigma_\theta^2$
is a constant from (a):

$$
\operatorname{Var}(\theta_n\mid\theta_{n-1})
= \underbrace{E\big[\operatorname{Var}(\theta_n\mid b_n)\mid\theta_{n-1}\big]}_{=\ \sigma_\theta^2}
+ \underbrace{\operatorname{Var}\big(E[\theta_n\mid b_n]\mid\theta_{n-1}\big)}_{=\ w^2\operatorname{Var}(b_n\mid\theta_{n-1})}
\;\le\; \sigma_\theta^2 + w^2\tau^2
$$

using Lemma L5 in the last step. (Anchoring on $\theta$ rather than $b$
avoids the crude $\operatorname{Var}(X+Y)\le2\operatorname{Var}X+2\operatorname{Var}Y$/Popoviciu split the
$b$-anchored composition needed for its second variance term — composing
*(b) then (a)* lets $\operatorname{Var}(b_n\mid\theta_{n-1})\le\tau^2$ enter directly via
L5, with no nonlinear $\delta(\cdot)$ inside a variance.) Adding the two
pieces and the $+1$ in $V$:

$$
\boxed{\;E\big[V(\theta_n) \mid \theta_{n-1}\big] \;\le\; \lambda\, V(\theta_{n-1}) \;+\; C\;}
\qquad
\lambda := \frac{1+w^2}{2} < 1
$$

$$
C := w^2(L\tau^2)^2\Big(1+\tfrac{2w^2}{1-w^2}\Big) + \sigma_\theta^2+w^2\tau^2 + 1-\lambda
$$

($\lambda$ is unchanged from the earlier $b$-anchored attempt — composing
the same exact-Gaussian step with the same bounded-bias/variance step in
either order contracts the mean gap at the same rate $w$ — but $C$ is
smaller here, since the variance term no longer needs the crude split.)
This drift condition holds for **every** $\theta_{n-1}\in\mathbb R$ (not
merely off a small set), because both building blocks of §4.2 are
uniform: (a) is exact for every $b_n$, and Lemmas L5–L6 hold for every
$\theta$ — a genuine, fully explicit Foster–Lyapunov pair on $\theta$,
with every constant computable from $(\tau,\sigma_0,\mu_0,L)$.

**4.4 From the $\theta$-drift to a bound on $P^n(x,A^c)$.** $A^c$ is a
statement about $b_n$, one sweep-step *ahead* of $\theta_{n-1}$ (§4.1).
By (b) (Lemmas L5–L6) and $(u+v)^2\le2u^2+2v^2$,

$$
E\big[(b_n-\mu_0)^2 \mid \theta_{n-1}\big]
= \big(\theta_{n-1}-\mu_0+\delta(\theta_{n-1})\big)^2 + \operatorname{Var}(b_n\mid\theta_{n-1})
\;\le\; 2(\theta_{n-1}-\mu_0)^2 + 2(L\tau^2)^2+\tau^2
$$

so, with $K^\star := \min\big((b^\star-\mu_0)^2,(b^\star+\mu_0)^2\big)+1 = \inf_{b\in A^c}V(b)$
as before and $(\theta_{n-1}-\mu_0)^2\le V(\theta_{n-1})$,

$$
E\big[V(b_n)\mid\theta_{n-1}\big] \;\le\; 2V(\theta_{n-1}) + 2(L\tau^2)^2+\tau^2+1,
\qquad
P(b_n\in A^c\mid\theta_{n-1}) \;\le\; \frac{E[V(b_n)\mid\theta_{n-1}]}{K^\star}
$$

(Markov's inequality, using $V(b)\ge K^\star$ for $b\in A^c$). Taking
$E[\cdot\mid\theta_0=x]$ and iterating the §4.3 drift,
$E[V(\theta_{n-1})\mid\theta_0=x]\le\lambda^{n-1}V(x)+C/(1-\lambda)$, gives
**fully explicitly**:

$$
P^n(x,A^c) \;\le\; \underbrace{\frac{2\lambda^{n-1}V(x) + 2C/(1-\lambda) + 2(L\tau^2)^2+\tau^2+1}{K^\star}}_{=:\,U_n(x)}
$$

**4.4$'$ The same bound, averaged against the stationary law instead of
iterated from $x$.** Write $g(t):=\pi(A^c\mid\theta=t,y)\in[0,1]$ for the
exact one-step conditional escape probability, so §4.4's derivation, up
to (and not including) "taking $E[\cdot\mid\theta_0=x]$", is really a
*pointwise, deterministic* statement: for every $t\in\mathbb R$,

$$
g(t) \;\le\; \Psi(t) \;:=\; \frac{2V(t) + 2(L\tau^2)^2+\tau^2+1}{K^\star}
$$

This single inequality feeds two different averages of $\theta$.
Averaging it against $P_\theta^{(n-1)}(x,\cdot)$ and iterating the §4.3
drift, as §4.4 does, gives $P^n(x,A^c)\le U_n(x)$ above. Averaging it
instead against the *stationary* law $\pi_\theta$ — using
$\theta_{n-1}\sim\pi_\theta \Rightarrow b_n\sim\pi(\cdot\mid y)$ exactly
(the $b$-marginal of the target itself, since $\pi_\theta$ is stationary
for the $(b,\theta)$-joint chain) and the standard stationary-moment
consequence of a Foster–Lyapunov drift, $E_{\pi_\theta}[V]\le C/(1-\lambda)$
(Meyn and Tweedie, 1993, Ch. 14: average the §4.3 drift against
$\pi_\theta$ itself, $E_{\pi_\theta}[V]=E_{\pi_\theta}\big[E[V(\theta_n)\mid\theta_{n-1}]\big]
\le\lambda E_{\pi_\theta}[V]+C$) — gives a genuine **marginal** bound on
$\pi(A^c)$ itself:

$$
\pi(A^c) \;=\; E_{\pi_\theta}[g(\theta)] \;\le\; E_{\pi_\theta}[\Psi(\theta)]
\;\le\; \underbrace{\frac{2C/(1-\lambda) + 2(L\tau^2)^2+\tau^2+1}{K^\star}}_{=:\,U_\infty} \;=\; \lim_{n\to\infty}U_n(x)
$$

This closes Open Item 1 below: unlike Proposition 1 (which bounds
$\pi(A^c\mid\theta,y)$ pointwise in $\theta$, leaving the outer average
over $\pi(\theta\mid y)$ undone), $U_\infty$ is a fully closed-form bound
on the actual marginal $\pi(A^c\mid y)$, via the same drift machinery
already built for $U_n(x)$ rather than a separate tail-integral
calculation. Both $U_n(x)$ and $U_\infty$ shrink to $0$ as $b^\star$
grows (since $K^\star\to\infty$), for any fixed $n$ (including
$n=\infty$) — this tail-awareness is inherited unchanged from §4.4 and
will matter in §4.5.

**4.5 Bounding the gap directly — why bounding the two pieces
separately (even the *rare-event* pieces) is still too crude, and what a
genuine coupling argument buys instead.** Because $b_n\mid\theta_{n-1}$
is drawn **exactly** from the target's own full conditional
$\pi(\cdot\mid\theta_{n-1},y)$ — that is what makes this a Gibbs sampler
rather than an approximation — both $\pi(A^c)=\int g\,d\pi_\theta$ and
$P^n(x,A^c)=\int g\,dP_\theta^{(n-1)}(x,\cdot)$ are the *same* functional
$g$ (§4.4$'$) integrated against, respectively, the true $\theta$-marginal
$\pi_\theta$ and the sampler's $(n-1)$-step $\theta$-marginal kernel
started at $x$ (the marginal chain $\{\theta_n\}$ is itself Markov, by
the same argument `Chapter-C03` uses for $X_1$: integrating the
intervening $b_n$ out of the composition $b_n\mid\theta_{n-1}$,
$\theta_n\mid b_n$).

The immediate bound $\big|\int g\,d\pi_\theta - \int g\,dP_\theta^{(n-1)}(x,\cdot)\big|
\le \|g\|_\infty\|P_\theta^{(n-1)}(x,\cdot)-\pi_\theta(\cdot)\|_{TV}$ uses
only $\|g\|_\infty\le1$ — but $\|g\|_\infty$ is genuinely close to $1$
here ($g(t)\to1$ as $t\to\pm\infty$: Proposition 1's tilted-Gaussian
bound *saturates* near its tail, it does not decay), reproducing exactly
the flaw §4.0 flags for the whole-space route.

An earlier version of this argument tried to fix this by splitting
$\theta$ into a bounded core $B$ (where $g$ is provably small, via
Proposition 1) and its complement $B^c$ (controlled by the §4.4/4.4$'$
drift machinery), bounding $\pi(A^c)$'s and $P^n(x,A^c)$'s restrictions
to $B^c$ **separately** and adding them. That fix is incomplete: bounding
two pieces of the *gap* separately, even when each piece is individually
small (a rare-event probability restricted to $B^c$), reduces to
$\pi_\theta(B^c)$ and $P_\theta^{(n-1)}(x,B^c)$ each appearing as
*additive, non-cancelling* terms of size roughly $C/(1-\lambda)$ divided
by $K_B^\star=\beta^2+1$ — a term that does **not** vanish as
$n\to\infty$ for any fixed $\beta$, because $\pi_\theta(B^c)$ alone does
not go away with $n$. It reintroduces, one level removed, exactly the
additive-floor problem the whole exercise is meant to avoid: the true
gap $\pi(A^c)-P^n(x,A^c)$ genuinely $\to0$ as $n\to\infty$ (both pieces
converge to the *same* limit $\pi(A^c)$), but a bound built by adding two
separately-controlled non-negative quantities can never see that
cancellation — it is structurally incapable of vanishing in $n$ unless
each piece is *individually* driven to $0$, which $\pi_\theta(B^c)$
(an $n$-free, fixed number) never is.

What is needed instead is a bound on the **difference** that can be
small for either of two independent reasons — the chains have run long
enough to have (with high probability) already coupled, *or* $A^c$
itself is a sufficiently rare event — rather than one that needs both
conditions to hold in the same term simultaneously. This is exactly what
a coupling argument on the $\theta$-marginal chain supplies.

**Coupling construction.** Let $(\theta_n^x)_{n\ge0}$ be the
$\theta$-marginal chain started at $\theta_0^x=x$, and let
$(\theta_n^\pi)_{n\ge0}$ be an independent copy started from
$\theta_0^\pi\sim\pi_\theta$ (so $\theta_n^\pi\sim\pi_\theta$ for every
$n$, by stationarity). Couple the two using the standard
minorization-based construction underlying §4.6's Rosenthal bound: run
them independently until both simultaneously lie in the small set $C_d$
of §4.6, then attempt to merge with probability $\varepsilon$ at each
such joint visit using the common measure $Q$ from the minorization; let
$T$ be the resulting random *coupling time*, the first sweep $n$ with
$\theta_n^x=\theta_n^\pi$. Once merged, run the two copies of the full
$(\theta,b)$ chain identically forever, so $\theta_{n'}^x=\theta_{n'}^\pi$
for all $n'\ge T$. In particular, on the event $\{T\le n-1\}$, draw
$b_n^x$ and $b_n^\pi$ from a single shared random draw of
$\pi(\cdot\mid\theta_{n-1}^x,y)=\pi(\cdot\mid\theta_{n-1}^\pi,y)$ — valid
since the two conditionals coincide exactly once $\theta_{n-1}^x=\theta_{n-1}^\pi$
— forcing $b_n^x=b_n^\pi$ on that event.

**Proposition 4$'$ (Cauchy–Schwarz gap bound).**

$$
\big|\pi(A^c)-P^n(x,A^c)\big| \;\le\;
\Big(\sqrt{P^n(x,A^c)}+\sqrt{\pi(A^c)}\Big)\sqrt{P(T>n-1)}
\;\le\;
\Big(\sqrt{U_n(x)}+\sqrt{U_\infty}\Big)\sqrt{R_{n-1}}
$$

where $R_{n-1}$ is §4.6's Rosenthal expression (below), reinterpreted as
a bound on $P(T>n-1)$ rather than directly on a TV distance.

*Proof.* Since $b_n^x=b_n^\pi$ on $\{T\le n-1\}$, the summand below
vanishes there, so

$$
\pi(A^c)-P^n(x,A^c) = E\big[\mathbf 1_{A^c}(b_n^\pi)-\mathbf 1_{A^c}(b_n^x)\big]
= E\Big[\big(\mathbf 1_{A^c}(b_n^\pi)-\mathbf 1_{A^c}(b_n^x)\big)\,\mathbf 1_{T>n-1}\Big]
$$

Both indicators are nonnegative, so
$\big|\pi(A^c)-P^n(x,A^c)\big| \le E[\mathbf 1_{A^c}(b_n^x)\mathbf 1_{T>n-1}]+E[\mathbf 1_{A^c}(b_n^\pi)\mathbf 1_{T>n-1}]$.
The event $\{T\le n-1\}$ is determined by $\theta_0^x,\dots,\theta_{n-1}^x$,
$\theta_0^\pi,\dots,\theta_{n-1}^\pi$, and independent coupling
randomness, all prior to the fresh draw of $b_n^x\mid\theta_{n-1}^x$; so
conditioning on $\theta_{n-1}^x$ leaves $\{T>n-1\}$ measurable and
$b_n^x\mid\theta_{n-1}^x$ unaffected, giving
$E[\mathbf 1_{A^c}(b_n^x)\mathbf 1_{T>n-1}]=E[g(\theta_{n-1}^x)\mathbf 1_{T>n-1}]$.
By Cauchy–Schwarz, then $g^2\le g$ (as $g\in[0,1]$),

$$
E[g(\theta_{n-1}^x)\mathbf 1_{T>n-1}]
\;\le\; \sqrt{E[g(\theta_{n-1}^x)^2]}\,\sqrt{P(T>n-1)}
\;\le\; \sqrt{E[g(\theta_{n-1}^x)]}\,\sqrt{P(T>n-1)}
\;=\; \sqrt{P^n(x,A^c)}\,\sqrt{P(T>n-1)},
$$

and likewise $E[\mathbf 1_{A^c}(b_n^\pi)\mathbf 1_{T>n-1}]\le\sqrt{\pi(A^c)}\sqrt{P(T>n-1)}$.
Summing the two and substituting §4.4/§4.4$'$'s $U_n(x)$, $U_\infty$
bounds gives the stated inequality. $\blacksquare$

**Why this genuinely satisfies both requirements at once.** As
$n\to\infty$ with $b^\star$ fixed, $R_{n-1}\to0$ (§4.6) while
$\sqrt{U_n(x)}+\sqrt{U_\infty}$ stays bounded (it converges to
$2\sqrt{U_\infty}\le2$) — so the *product* vanishes, with no additive,
non-vanishing floor of the kind the core/tail split left behind. As
$b^\star\to\infty$ with $n$ fixed, $U_n(x)\to0$ and $U_\infty\to0$ (§4.4,
§4.4$'$) while $\sqrt{R_{n-1}}$ is unaffected — so the bound also
shrinks toward $0$ purely from $A^c$ becoming rare, without needing
$n\to\infty$ at all. The two requirements are satisfied by two different
factors of a product rather than by one term trying to do both jobs, and
each factor is the *existing* machinery (§4.4/§4.4$'$'s drift bound,
§4.6's Rosenthal bound) — nothing new needed to be derived beyond the
coupling argument and the Cauchy–Schwarz split itself.

**4.6 Rosenthal's (1995) Theorem 12, applied to the $\theta$-chain, as a
coupling-time bound.** Given the §4.3 drift $(\lambda, C)$, fix any
$d > 2C/(1-\lambda)$ and take the small set $C_d := \{\theta : V(\theta)\le d\} =
[\mu_0-\sqrt{d-1},\mu_0+\sqrt{d-1}]$ — a compact interval, since sublevel
sets of the coercive $V$ are automatically compact. The one-step
$\theta$-kernel $q_\theta(\theta'\mid\theta) := \int_{\mathbb R}\pi(b\mid\theta,y)\,
\phi\big(\theta';wb+(1-w)\mu_0,\sigma_\theta^2\big)\,db$ is jointly
continuous in $(\theta,\theta')$ ($\pi(\cdot\mid\theta,y)$ and the
Gaussian kernel for $\theta'\mid b$ are each continuous, and Lemma L2's
uniform score bound lets dominated convergence pass continuity through
the $b$-integral) and strictly positive everywhere (neither factor has
compact support), so it attains a
positive minimum $\varepsilon>0$ on the compact product $C_d\times
[\theta'_-,\theta'_+]$, giving a minorization $q_\theta(\theta,\cdot)\ge
\varepsilon\,Q(\cdot)$ for $\theta\in C_d$, with $Q$ uniform on
$[\theta'_-,\theta'_+]$. With $(\lambda,C,\varepsilon,d)$ so defined, the
coupling construction above — the same one underlying Rosenthal (1995,
Theorem 12) — satisfies, for any $x$ and any $0<r<1$,

$$
P(T>k) \;\le\;
(1-\varepsilon)^{rk} \;+\; \alpha^{-k}(\alpha\Lambda)^{rk}
\left[1 + \frac{C}{1-\lambda} + V(x)\right] \;=:\; R_k,
$$

$$
\alpha^{-1} := \frac{1+2C+\lambda d}{1+d} \;<\; 1,
\qquad
\Lambda := 1+2(\lambda d+C).
$$

(Rosenthal's Theorem 12 quotes only the resulting
$\|P_\theta^{(k)}(x,\cdot)-\pi_\theta(\cdot)\|_{TV}\le R_k$ — the TV
consequence of $P(T>k)\le R_k$ via the standard coupling inequality — but
his proof, via renewal/regeneration on $C_d$, establishes the bound on
$P(T>k)$ itself, which is the form Proposition 4$'$ needs.) Combining
with Proposition 4$'$ ($k=n-1$):

$$
\big|\pi(A^c)-P^n(x,A^c)\big|
\;\le\;
\Big(\sqrt{U_n(x)}+\sqrt{U_\infty}\Big)
\sqrt{(1-\varepsilon)^{r(n-1)} \;+\; \alpha^{-(n-1)}(\alpha\Lambda)^{r(n-1)}
\left(1 + \frac{C}{1-\lambda} + V(x)\right)}
$$

Every constant here is an explicit function of $(\lambda,C,\varepsilon,d,r)$
— $(\lambda,C)$ from §4.3, $d$ any fixed value above $2C/(1-\lambda)$,
$r$ free to optimize — with **only $\varepsilon$** left as an unclosed
numeric quantity (Open Item 2 below).

---

## 5. Assembled bound

Combining §3 (Proposition 1), §4.4/§4.6, and Lemma 1/Theorem 1 of
`vignette("Chapter-C05")` for the on-$A$ term, with $x=\theta_0$ per
§4.1:

$$
\big\|P^n(x,\cdot) - \pi(\cdot)\big\|_{TV}
\;\le\;
\underbrace{\big\|P^n_A(x,\cdot)-\pi_A(\cdot)\big\|_{TV}}_{\text{on-}A\text{ gap, Chapter-C05 heuristic}}
\;+\;
4\,P^n(x,A^c)
\;+\;
2\big[\pi(A^c) - P^n(x,A^c)\big]
$$

with $4P^n(x,A^c)$ bounded fully explicitly by §4.4 (no abstract
constants), and the last term bounded by twice the Cauchy–Schwarz gap
bound of §4.5/§4.6 (Proposition 4$'$ combined with the Rosenthal
coupling-time bound) — every constant in it explicit given
$(\lambda,C,d,r)$, with $\varepsilon$ the sole remaining unclosed numeric
quantity (Open Item 2). $\pi(A^c)$ itself — needed for the on-$A$
normalization-correction machinery of `ELLIPSOID_TV_BOUND.md` §3, and as
a standalone sanity check on $4P^n(x,A^c)$'s bound — is available two
ways: pointwise in $\theta$ from Proposition 1 (§3), or as the fully
marginal, closed-form bound $U_\infty$ from §4.4$'$. All pieces are
computable from $(\theta_0,\tau,\mu_0,\sigma_0,b^\star,y,d,r)$, and every
constant is explicit except $\varepsilon$.

---

## 6. What generalizes, and what is special about this case

**Special to this case (scalar, single group, no covariates):**

- $A$ is an explicit interval $[-b^\star,b^\star]$ (§2) rather than an
  implicit level set of a matrix eigenvalue functional.
- $\pi(A^c\mid\theta,y)$ has a fully closed form (Proposition 1) because
  the "eigenvalue functional" $g(b)=N\,p(b)(1-p(b))$ here *is* the
  curvature itself, and its tail reduces to a one-dimensional
  tilted-Gaussian calculation.
- $\theta\mid b,y$ is exactly Gaussian (Remark, §3), since the group
  likelihood never depends on $\theta$ directly.

**What should generalize (candidate next steps, in order of how directly
this case's machinery extends):**

1. **Multiple groups, still scalar random effects, no covariates.** Each
   group $j$ gets its own $A_j=[-b_j^\star,b_j^\star]$, its own $L_j$ and
   Proposition-1 bound, and $A=\bigcap_j A_j$; Claim 2's union bound
   (§4 of the umbrella note) combines them. Because $\theta$ (the block
   updated second, per §4.1) is still scalar even with multiple groups —
   it is the *random effects* that become a vector $\mathbf b$, not the
   shared population parameter — the §4.3 drift argument should extend
   with $V(\theta)$ left scalar and Lemmas L5–L6 applied group-by-group
   (each $b_j\mid\theta$ contributing its own bounded mean-bias/variance
   to the composed update); the aggregation across groups happens inside
   $\eta(\theta_{n-1})$ and $\operatorname{Var}(b_n\mid\theta_{n-1})$ of §4.3,
   which needs checking before the single-group drift constants carry
   over unchanged.
2. **A random-effect design matrix $D_j$ with $p_j>1$ columns (e.g.
   random intercept *and* slope).** Lemmas L1–L2 (bounded/Lipschitz score,
   tangent-line bound) generalize immediately using the operator-norm
   Lipschitz constant of the score *vector*; Lemma L3's one-dimensional
   tilted-Gaussian tail becomes a multivariate Gaussian tilted-tail
   calculation (still closed form, via the multivariate normal MGF), but
   §2's interval $A$ becomes a genuine sublevel set of $\lambda_{\max}$
   or superlevel set of $\lambda_{\min}$ of
   $D_j'\operatorname{diag}(c(\eta_i))D_j$ — this is exactly where
   `ELLIPSOID_TV_BOUND.md` §6.4 item 2's gap reappears, now for a floor
   (not ceiling) split, and the union-bound sketch mentioned there (bound
   $g_j$ via $\max_i c(\eta_i)$, tail-bound each $c(\eta_i)$ via a
   Gaussian tail on $\eta_i$) is the natural next attempt.
3. **Feeding this back into `ELLIPSOID_TV_BOUND.md` §6.** Once step 2 is
   worked out, §6.1's ceiling-based split and this note's floor-based
   split could be merged into a single two-sided (floor-and-ceiling) safe
   region for the general logit block, closing open item 2 of §6.4.

---

## 7. Open items

1. ~~Marginal (not merely $\theta$-conditional) bound on $\pi(A^c\mid y)$~~
   (§3 Remark). **Resolved** by §4.4$'$'s $U_\infty$: averaging the same
   pointwise escape bound $g(t)\le\Psi(t)$ used for $P^n(x,A^c)$ against
   the stationary law $\pi_\theta$ (rather than iterating it from $x$)
   gives a closed-form bound on the true marginal $\pi(A^c\mid y)$,
   using only the §4.3 drift's stationary-moment consequence — no
   separate outer-tail calculation on $\pi(\theta\mid y)$ was needed.
2. **The minorization constant $\varepsilon$ in §4.6.** The small set
   $C_d$ is pinned down explicitly (any sublevel set of $V(\theta)$), and
   its existence as a genuine small set is justified by
   continuity/strict positivity of $q_\theta(\theta'\mid\theta)$ on the
   compact product $C_d\times[\theta'_-,\theta'_+]$ — but the actual
   minimum $\varepsilon$ of $q_\theta(\theta'\mid\theta)$ over that
   product is not computed here. Given §4.6's Rosenthal (1995) Theorem 12
   formula, $\varepsilon$ is now the *only* unclosed quantity in the
   entire assembled bound of §5 — a narrower, more tractable version of
   the same gap flagged for probit/Poisson in `ELLIPSOID_TV_BOUND.md` §9
   item 5, since both the drift constants and the small set are already
   explicit, leaving one numeric minimization over a compact product.
3. **Sharpening the Young's-inequality split in §4.3.** The choice
   $\epsilon=(1-w^2)/(2w^2)$ was made for universality (valid for every
   $w\in(0,1)$) rather than tightness; for a specific $w$, optimizing
   $\epsilon$ (or avoiding the split altogether via a more careful
   second-moment calculation) would tighten $\lambda$ and $C$.
4. **Extending §2's interval safe region to $p_j>1$ random-effect blocks**
   (§6, item 2) — the natural next case study before attempting to fold
   this back into the general §6 of `ELLIPSOID_TV_BOUND.md`.
5. **Tightness of the Lipschitz constant $L$.**
   $L=\sum_i\max(y_i,n_i-y_i)$ is the sharpest bound available from Lemma
   L1 alone; whether a tighter, still-closed-form bound exists by using
   the concavity of $\ell$ more fully (rather than only its global
   Lipschitz constant) is not investigated here.
6. **Extending the §4.3 drift argument to multiple groups** (§6, item 1)
   — checking that aggregating the per-group Lemma L5/L6 bounds inside
   the (still scalar, since $\theta$ stays scalar) drift condition
   preserves a bound of the same form with explicit constants.
7. **Optimizing $r\in(0,1)$ in §4.6's Rosenthal bound**, and comparing
   against the Roberts–Tweedie (1999) alternative form, once
   $\varepsilon$ (item 2) is available numerically — not attempted here.
8. **Sharpening the Cauchy–Schwarz split in §4.5.** Proposition 4$'$
   bounds $E[g(\theta_{n-1})^2]$ by the crude $E[g(\theta_{n-1})]$ (using
   only $g\le1$); using Proposition 1's sharper pointwise bound
   $g(t)\le G(t)$ on one factor of $g^2=g\cdot g$ instead — i.e.
   $E[g(\theta_{n-1})^2]\le E[g(\theta_{n-1})G(\theta_{n-1})]$ — would
   tighten the bound further, at the cost of a less clean closed form; not
   attempted here.

---

## References

- `ELLIPSOID_TV_BOUND.md` — the general safe/unsafe decomposition
  (§§1–4) and its logit specialization (§6), whose ceiling-based split
  and open closed-form question (§6.4 item 2) this note complements with
  a floor-based, closed-form special case.
- `vignette("Chapter-C05")` — the model-free statement of Lemma 1,
  Lemma 6/7, and the practical Theorem-1 form of the three-term TV bound
  assembled in §5 above.
- `vignette("Chapter-C03")` — Definition 4 and its remark ("the second
  block is updated first, and the first block is updated second... all
  the results below are statements about $X_1$, the block updated
  second"), the convention §4.1 aligns this note's $(\theta,b)$ ordering
  with, and which motivates anchoring §4.3's drift on $\theta$ rather
  than $b$.
- Mengersen, K.L. and Tweedie, R.L. (1996). *Rates of convergence of the
  Hastings and Metropolis algorithms.* Annals of Statistics — the
  whole-space uniform-ergodicity route deliberately set aside in §4.0 as
  too loose for a tail-specific gap.
- Rosenthal, J.S. (1995). *Minorization conditions and convergence rates
  for Markov chain Monte Carlo.* Journal of the American Statistical
  Association, 90, 558–566, Theorem 12 — the explicit drift-and-
  minorization total variation bound applied to the $\theta$-marginal
  chain in §4.6, replacing the earlier abstract Meyn–Tweedie-$f$-norm
  constant $R$ with a fully cited formula in $(\lambda,C,\varepsilon,d,r)$.
- Brascamp, H.J. and Lieb, E.H. (1976). *On extensions of the
  Brascamp–Lieb–Luttinger inequality...* — the variance bound for
  log-concave densities used in Lemma L5.
- Meyn, S.P. and Tweedie, R.L. (1993). *Markov Chains and Stochastic
  Stability* — the stationary-moment consequence of a Foster–Lyapunov
  drift, $E_{\pi_\theta}[V]\le C/(1-\lambda)$, used in §4.4$'$ to obtain
  $U_\infty$.
