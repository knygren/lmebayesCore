# Gaussian majorization on the escape set

Companion to `vignette("Chapter-C05")` (Total Variation Bounds for Restricted
Two-block Gibbs Samplers) and to `CHAPTER_C05_IMPLEMENTATION.md`. Chapter C05
establishes, for every tail budget \(\delta>0\), the existence of a \(d(\delta)\)
with \(\pi_\gamma(\widetilde C_d^{\,c})\le\delta\); it does not say how to
produce the number. This note supplies a route: **majorize the marginal
posterior by a scaled Gaussian on the escape set**, and read the escape
probability off that Gaussian.

> **Status: derived here, not in Chapter C05.** §1 is a one-line integration and
> is not in doubt. §3 onward rests on a chain of standard inequalities
> (Brascamp–Lieb, strong convexity) assembled here for the first time, checked
> numerically in §9 but not reviewed. Nothing below should be cited as proved.

**A word on the name.** "Majorization" here means *pointwise domination of one
density by a scaled other* — \(\pi_\gamma\le M q_2\) — not Schur majorization of
vectors. The term is chosen for the symmetry with the minorization
\(q(\gamma,\cdot\mid C)\ge\varepsilon Q_C(\cdot)\) that drives Theorem 1: the
certificate minorizes inside \(C\) and majorizes outside it.

**Notation.** As in `CHAPTER_C05_IMPLEMENTATION.md` §4A.0. The joint precision
blocks of \((\gamma,\beta)\) are
\(A=P_{11}=\Lambda_\gamma+\sum_jH_j^\top P_bH_j\),
\(B=P_{12}P_{22}^{-1}P_{21}=S=\sum_jH_j^\top P_bV_jP_bH_j\), and
\(\Pi=A-B\) is the marginal posterior precision;
\(V_j=\mathrm{Cov}(\beta_j\mid\gamma,y)\).
\(\Phi(\gamma)=-\log\pi_\gamma(\gamma)\),
\(\bar\Phi=\Phi-\Phi(\gamma^\star)\),
\(\Psi=\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}-\bar\Phi\) is the deficiency
gap, and \(\widetilde C_d=\{\Psi\le d\}\).

---

## 1. The statement

Let \(P_{22}^{\mathrm{LB}}\) be any symmetric matrix with
\(P_{22}\succeq P_{22}^{\mathrm{LB}}\), and set

\[
B^{\mathrm{UB}}:=P_{12}\bigl(P_{22}^{\mathrm{LB}}\bigr)^{-1}P_{21},
\qquad
\Pi^{\mathrm{LB}}:=P_{11}-B^{\mathrm{UB}} .
\]

Assume \(\Pi^{\mathrm{LB}}\succ0\) and let \(q_2\) be the density of
\(Q_2:=N\bigl(\gamma^\star,(\Pi^{\mathrm{LB}})^{-1}\bigr)\).

> **Proposition 1 (Gaussian majorization on the escape set).** Suppose there
> exists \(M\ge0\) with
> \[
> \pi_\gamma(\gamma)\ \le\ M\,q_2(\gamma)\qquad\text{for every }\gamma\in\widetilde C_d^{\,c}.
> \]
> Then
> \[
> \pi_\gamma\bigl(\widetilde C_d^{\,c}\bigr)\ \le\ M\,Q_2\bigl(\widetilde C_d^{\,c}\bigr).
> \]

*Proof.* Integrate the pointwise inequality over \(\widetilde C_d^{\,c}\). \(\blacksquare\)

Three things are worth saying about a statement this simple.

**The hypothesis is local to the escape set.** Domination is required only where
\(\Psi>d\) — that is, only in the tails, which is where a Gaussian envelope is
easiest to certify and where the mode's behaviour is irrelevant. §2 shows the
gain from this restriction is real, not cosmetic.

**The hypothesis is agnostic about provenance.** Proposition 1 does not care how
\(q_2\) or \(M\) were obtained. §3 constructs them from \(P_{22}^{\mathrm{LB}}\)
because that is the construction available in this package, but any argument
producing a dominating Gaussian on the tail qualifies. The \(P_{22}\) lower
bound is therefore a *sufficient device*, not part of the claim.

**It converts \(\delta\) from an existential into a number.** Theorem 1's
conclusion becomes

\[
\bigl\|q_n(\gamma,\cdot\mid C)-\pi_\gamma\bigr\|_{TV}
\ \le\ (1-\varepsilon)^n+M\,Q_2\bigl(C^c\bigr),
\qquad C=\widetilde C_d,
\]

with both terms explicit: the first a minorization constant, the second a
Gaussian measure. Chapter C05's Theorem 2 is then the statement that the second
term can be driven below any \(\delta\) by increasing \(d\), rather than an
existence claim about \(d(\delta)\).

**On naming.** \(M\) is written here rather than \(\varepsilon_2\), because
\(\varepsilon\) denotes the minorization constant throughout Chapter C05 and
satisfies \(\varepsilon\le1\), whereas \(M\ge1\) whenever the envelope is
global (§4). Reusing \(\varepsilon\) would invite exactly the wrong intuition
about its size.

---

## 2. The sharp constant, and why the escape-set restriction matters

The smallest admissible constant is

\[
M^\star:=\sup_{\gamma\in\widetilde C_d^{\,c}}\frac{\pi_\gamma(\gamma)}{q_2(\gamma)} .
\]

Write

\[
R(\gamma):=\tfrac12\|\gamma-\gamma^\star\|^2_{\Pi^{\mathrm{LB}}}-\bar\Phi(\gamma),
\qquad\text{so}\qquad
\frac{\pi_\gamma(\gamma)}{q_2(\gamma)}=c\,e^{R(\gamma)},
\quad
c:=\frac{\pi_\gamma(\gamma^\star)}{q_2(\gamma^\star)} .
\]

**\(R\) is concave with global maximum \(0\) at \(\gamma^\star\).** Concavity is
\(\nabla^2R=\Pi^{\mathrm{LB}}-\nabla^2\Phi=\Pi^{\mathrm{LB}}-\Pi(\gamma)\preceq0\),
which is §3's Lemma 2; stationarity is \(\nabla\Phi(\gamma^\star)=0\); and
\(R(\gamma^\star)=0\) by construction. Hence \(R\le0\) everywhere and the
density ratio is largest **at the mode**.

That is the point of restricting to \(\widetilde C_d^{\,c}\): \(\gamma^\star\) is
the one point guaranteed *not* to lie in the escape set. Since \(R\) is concave
with maximum at \(\gamma^\star\), it is nonincreasing along every ray out of
\(\gamma^\star\), so along each ray the supremum over the part lying in
\(\widetilde C_d^{\,c}\) is attained where the ray crosses
\(\partial\widetilde C_d\). Therefore

\[
M^\star=c\,\exp\Bigl(\max_{\partial\widetilde C_d}R\Bigr)\ \le\ c ,
\]

and the inequality is strict whenever \(R<0\) on the boundary.

**A closed form for the boundary maximum.** Complementarity gives
\(\bar\Phi=\tfrac12\|\cdot\|^2_{P_{11}}-\Psi\), so

\[
R(\gamma)=\Psi(\gamma)-\tfrac12\|\gamma-\gamma^\star\|^2_{B^{\mathrm{UB}}},
\]

and on \(\partial\widetilde C_d\) we have \(\Psi=d\). Hence

\[
\boxed{\ M^\star=c\,\exp\Bigl(d-\tfrac12 r_d\Bigr),
\qquad
r_d:=\min_{\partial\widetilde C_d}\|\gamma-\gamma^\star\|^2_{B^{\mathrm{UB}}}\ \ge\ 2d.\ }
\]

The inequality \(r_d\ge2d\) is \(\Psi\le\tfrac12\|\cdot\|^2_{B^{\mathrm{UB}}}\)
(Lemma 3), and it is what makes \(M^\star\le c\). Computing \(r_d\) is a
quadratic minimization over a convex level set — not free, but a standard convex
program with the same E-step ingredients as everything else.

**Sanity check.** When \(P_{22}^{\mathrm{LB}}=P_{22}\) and the model is Gaussian,
\(\Psi=\tfrac12\|\cdot\|^2_B\) exactly with \(B^{\mathrm{UB}}=B\), so
\(r_d=2d\) and \(M^\star=c=1\): the envelope is the posterior and Proposition 1
returns the exact escape probability. Unlike a Chernoff bound, this route is
*exact in the tight limit* — verified to round-off in §9.

**\(M^\star\) improves with \(d\).** Under closure,
\(r_d=2d\,\lambda_{\min}\bigl(B^{-1/2}B^{\mathrm{UB}}B^{-1/2}\bigr)\), so

\[
M^\star=c\,\exp\bigl(-d(\lambda_{\min}-1)\bigr),
\qquad \lambda_{\min}=\lambda_{\min}\bigl(B^{-1/2}B^{\mathrm{UB}}B^{-1/2}\bigr)\ge1 .
\]

The sharp constant therefore *decays exponentially in \(d\)* whenever the floor
is not tight, rather than sitting at the fixed value \(c\). §9 measures
\(M^\star/c\) ranging down to \(1.4\times10^{-3}\) across the test grid — nearly
three orders of magnitude recovered purely from restricting the supremum to the
escape set. This is the strongest single argument for stating the hypothesis on
\(\widetilde C_d^{\,c}\) rather than globally.

---

## 3. Constructing \(q_2\) from a lower bound on \(P_{22}\)

Three lemmas, each standard, assembled to verify Proposition 1's hypothesis
globally (hence in particular on \(\widetilde C_d^{\,c}\)).

> **Lemma 1 (curvature floor to covariance ceiling).** If
> \(P_{22}\succeq P_{22}^{\mathrm{LB}}\succ0\) then
> \(V_j=\mathrm{Cov}(\beta_j\mid\gamma,y)\preceq(P_{22}^{\mathrm{LB}})_j^{-1}\)
> for every \(\gamma\), hence \(S(\gamma)\preceq B^{\mathrm{UB}}\).

This is Brascamp–Lieb, the same inequality that gives \(V_j\preceq P_b^{-1}\) in
`CHAPTER_C05_IMPLEMENTATION.md` §4A.6. It deserves emphasis because **\(P_{22}\)
changes meaning across the step**: on the left it is the curvature block of the
joint negative log density, \(P_b+(-\nabla^2\ell_j)\); on the right the bound
constrains a *covariance*. The two coincide only when the \(\beta\)-conditional
is exactly Gaussian, and Brascamp–Lieb is what licenses the substitution in
general. Log-concavity of the \(\beta\)-conditional — (H2) — is the requirement,
and it survives truncation of \(\beta\) to a convex set (§7).

> **Lemma 2 (strong convexity).** Under Lemma 1,
> \(\nabla^2\Phi(\gamma)=P_{11}-S(\gamma)\succeq\Pi^{\mathrm{LB}}\) for every
> \(\gamma\); that is, \(\pi_\gamma\) is \(\Pi^{\mathrm{LB}}\)-strongly
> log-concave.

> **Lemma 3 (global envelope).** Under Lemma 2 and \(\nabla\Phi(\gamma^\star)=0\),
> \[
> \bar\Phi(\gamma)\ \ge\ \tfrac12\|\gamma-\gamma^\star\|^2_{\Pi^{\mathrm{LB}}},
> \qquad\text{equivalently}\qquad
> \pi_\gamma(\gamma)\ \le\ c\,q_2(\gamma)\ \ \text{for all }\gamma,
> \]
> and consequently
> \(\Psi(\gamma)\le\tfrac12\|\gamma-\gamma^\star\|^2_{B^{\mathrm{UB}}}\).

So Proposition 1's hypothesis holds with \(M=c\), and §2 sharpens that to
\(M^\star=c\,e^{d-r_d/2}\).

**The binding constraint is \(\Pi^{\mathrm{LB}}\succ0\).** A crude floor makes
\(B^{\mathrm{UB}}\) large enough that \(\Pi^{\mathrm{LB}}\) loses definiteness,
\(q_2\) ceases to exist, and the route is vacuous. This is the quantitative
limit on how weak the \(\beta\)-restriction of §7 may be, and it coincides with
the condition \(\kappa_0<1\) of `CHAPTER_C05_IMPLEMENTATION.md` §4B.3 — the
floor that makes the tilting exponent nonzero is the same floor that makes this
envelope exist.

---

## 3A. Existence of a certified floor

Lemma 1 assumes a floor \(P_{22}^{\mathrm{LB}}\). This section says when one
exists. It is the \(\beta\)-side counterpart of Chapter C05's Theorem 2: that
result produces a certified set for \(\gamma\), this one produces a certified set
for the group coefficients, and together they close the two-budget certificate
of §7.

Throughout, the certified \(\beta\)-set is written \(\widetilde B\), never \(C\):
\(\widetilde C_d\subset\mathbb R^q\) is the \(\gamma\)-set of Chapter C05 and
\(\widetilde B_r\subset\mathbb R^{Jp_{\mathrm{re}}}\) is the \(\beta\)-set built
here. The two are indexed by different levels of different functions and should
never be conflated.

### 3A.1 A necessary sharpening of the statement

The natural statement — *for every \(\delta>0\) there is a compact
\(\widetilde B\) and a positive definite \(P_{22}^{\mathrm{LB}}\) with
\(P_{22}(\beta)\succ P_{22}^{\mathrm{LB}}\) on \(\widetilde B\)* — is **true but
vacuous as written**.
Since \(P_{22}(\beta)=\mathrm{blockdiag}\bigl(P_b+G_j(\beta_j)\bigr)\) with
\(G_j=Z_j^\top W_jZ_j\succeq0\), one always has

\[
P_{22}(\beta)\ \succeq\ I_J\otimes P_b\ \succ\ 0
\qquad\text{for every }\beta\in\mathbb R^{Jp_{\mathrm{re}}},
\]

with no compactness and no assumptions. That floor is exactly the trivial one
that yields \(B^{\mathrm{UB}}=P_{11}^{\mathrm{RE}}\), \(\kappa_0=1-\mathtt{pop.pwt}\),
and the numerically vacuous certificate of
`CHAPTER_C05_IMPLEMENTATION.md` §4B.9.

The content is therefore not that \(P_{22}^{\mathrm{LB}}\succ0\), but that it can
be taken **strictly above the prior block** — that the *likelihood* contributes a
positive curvature margin. Theorem 2 is stated that way.

### 3A.2 The set: level sets of the unnormalized marginal

A ball, or a product of per-group boxes, is an arbitrary shape imposed from
outside. The natural object is already available:
`LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §4–§5 integrates \(\gamma\) out of the
prior exactly and gives the **unnormalized marginal on the group coefficients**

\[
\widetilde\pi(\beta\mid y)
\;=\;
\exp\Bigl(\sum_{j=1}^J\ell_j(\beta_j)\Bigr)\,
\exp\Bigl(-\tfrac12(\beta-\mu_\beta)^\top\Lambda_\beta(\beta-\mu_\beta)\Bigr),
\tag{M}
\]

with \(\mu_\beta=\mathcal W\mu_0\) and the **prior Schur precision**
\(\Lambda_\beta=\mathfrak P_{22}-\mathfrak P_{21}\mathfrak P_{11}^{-1}
\mathfrak P_{12}\), where \(\mathfrak P\) denotes the *prior-only* joint
precision of that note's §3 — not the posterior blocks \(P_{11},P_{22}\) used
elsewhere here. Note that \(\Lambda_\beta\) is **not** block diagonal:
integrating \(\gamma\) couples the groups.

Let \(\beta^\dagger:=\arg\max_\beta\widetilde\pi(\beta\mid y)\) and define the
**marginal deficiency**

\[
\Xi(\beta)
\;:=\;
\log\widetilde\pi(\beta^\dagger\mid y)-\log\widetilde\pi(\beta\mid y)
\;\ge\;0,
\qquad
\boxed{\ \widetilde B_r:=\{\beta:\Xi(\beta)\le r\}\ }
\tag{B}
\]

This is the exact \(\beta\)-side analogue of
\(\widetilde C_d=\{\gamma:\Psi(\gamma)\le d\}\), with \(\Xi\) playing the role of
the gap \(\Psi\), and it inherits the same four properties.

**Normalizer-free.** \(\Xi\) is a *difference* of log-densities, so the unknown
\(Z(y)=\int\widetilde\pi\) cancels identically. This is exactly why the
*unnormalized* marginal is the right object to level-set: the level is measured
against the mode, not against an intractable constant. `LOGIT_MARGINAL_INTEGRATE_GAMMA.md`
§8 and §12 record that \(Z(y)\) has no closed form for logit; (B) never needs it.

**Convex.** \(\Xi(\beta)=\sum_j\bigl(-\ell_j(\beta_j)\bigr)
+\tfrac12\|\beta-\mu_\beta\|^2_{\Lambda_\beta}+\text{const}\). Each \(-\ell_j\)
is convex (canonical link, log-concave likelihood) and \(\Lambda_\beta\succeq0\)
as a Schur complement of a positive semidefinite matrix, so \(\Xi\) is convex and
every \(\widetilde B_r\) is convex — which is what Lemma 1 needs for
Brascamp–Lieb to survive truncation (§7).

**Compact.** \(\Xi\) is coercive whenever \(\Lambda_\beta\succ0\), which holds iff
the population prior is proper (\(\Lambda_\gamma\succ0\)), since \(-\ell_j\) is
bounded below. In the flat limit \(\Lambda_\gamma=0\) the matrix \(\Lambda_\beta\)
is singular along the common-shift direction \(\mathcal W v\), and coercivity must
come from the likelihood instead — for logit, from non-separation of every group,
which is the same estimability condition (H3b) that Theorem 2 needs for a
different purpose.

**Exhausting.** \(\widetilde B_r\uparrow\mathbb R^{Jp_{\mathrm{re}}}\) as
\(r\uparrow\infty\), so \(\pi(\widetilde B_r^{\,c}\mid y)\downarrow0\) by
dominated convergence. Writing \(r(\delta)\) for the level attaining budget
\(\delta\) mirrors \(d(\delta)\) on the \(\gamma\) side.

The gain over a box is not only aesthetic. §3A.8 argues that the box should be
stretched where the marginal is diffuse and tightened where \(p\) is extreme, and
should carry unequal budget across groups; the level set does all of that
automatically, from one scalar \(r\), because it follows the contours of the
marginal itself.

### 3A.3 Statement

> **Theorem 2 (certified curvature floor on a \(\beta\)-set).** Assume
>
> - **(P)** the joint posterior \(\pi(\gamma,\beta\mid y)\) is proper;
> - **(W)** each family's GLM weight function \(w(\cdot)\) is continuous and
>   strictly positive on \(\mathbb R\);
> - **(H3b)** group-wise estimability: every \(Z_j\) has full column rank.
>
> Then for every \(\delta>0\) there exist a level \(r(\delta)<\infty\), with
> \(\widetilde B_{r(\delta)}\subset\mathbb R^{Jp_{\mathrm{re}}}\) the compact
> convex set (B), and matrices \(\Gamma_j(\delta)\succ0\) such that, with
> \(P_{22}^{\mathrm{LB}}(\delta):=\mathrm{blockdiag}\bigl(P_b+\Gamma_j(\delta)\bigr)\),
>
> 1. \(\pi\bigl(\widetilde B_{r(\delta)}^{\,c}\mid y\bigr)\le\delta\), and
> 2. \(P_{22}(\beta)\succeq P_{22}^{\mathrm{LB}}(\delta)\succ I_J\otimes P_b\)
>    for every \(\beta\in\widetilde B_{r(\delta)}\).

*Proof.* For (1): by (P), \(\pi(\cdot\mid y)\) is a probability measure, and by
§3A.2 the sets \(\widetilde B_r\) are compact and increase to
\(\mathbb R^{Jp_{\mathrm{re}}}\); choose \(r(\delta)\) with
\(\pi(\widetilde B_{r}^{\,c}\mid y)\le\delta\).

For (2): \(\widetilde B_r\) is compact, so the continuous map
\(\beta\mapsto\eta_{j,i}=o_{j,i}+d_{j,i}^\top\beta_j\) is bounded on it; write
\(\eta_{j,i}^{\max}:=\max_{\beta\in\widetilde B_r}|\eta_{j,i}(\beta)|<\infty\).
By (W), \(w\) attains a strictly positive minimum on the compact interval
\([-\eta_{j,i}^{\max},\eta_{j,i}^{\max}]\); write
\(\underline w_j:=\min_i\min_{|\eta|\le\eta_{j,i}^{\max}}w(\eta)>0\). Then for
every \(\beta\in\widetilde B_r\),

\[
G_j(\beta_j)=\sum_i w_{j,i}\,z_{j,i}z_{j,i}^\top
\ \succeq\ \underline w_j\sum_i z_{j,i}z_{j,i}^\top
=\underline w_j\,Z_j^\top Z_j\ =:\ \Gamma_j ,
\]

and \(\Gamma_j\succ0\) because \(\underline w_j>0\) and \(Z_j^\top Z_j\succ0\) by
(H3b). \(\blacksquare\)

The proof is pure existence: it uses compactness for (2) and tightness for (1),
and produces no numbers. §3A.4 replaces both steps with closed forms.

### 3A.4 Rays from the mode

Both halves of Theorem 2 become computable by working along rays
\(t\mapsto\beta^\dagger+tu\) out of the mode and continuing *past* the boundary
of \(\widetilde B_r\). Write

\[
\varphi_u(t):=\Xi(\beta^\dagger+tu),
\qquad
t_r(u):=\sup\{t:\varphi_u(t)\le r\},
\]

so \(\varphi_u\) is convex with \(\varphi_u(0)=0\) — the mode is the minimum —
and \(t_r(u)\) is where the ray crosses \(\partial\widetilde B_r\). Convexity
through the origin makes \(\varphi_u(t)/t\) nondecreasing, which is the whole
engine:

\[
\varphi_u(t)\ \le\ \frac{r\,t}{t_r(u)}\quad(t\le t_r),
\qquad
\varphi_u(t)\ \ge\ \frac{r\,t}{t_r(u)}\quad(t\ge t_r).
\tag{R}
\]

Past the boundary the curvature of \(\Xi\) is *declining* along every ray — the
logit weights vanish as \(|\eta|\to\infty\), so \(\nabla^2\Xi\to\Lambda_\beta\)
and the growth of \(\varphi_u\) is asymptotically linear. (R) is therefore not
merely a convenient relaxation: linear growth is the true behaviour, and any
bound assuming faster-than-linear growth outside \(\widetilde B_r\) is asserting
something false. §3A.4.3 records what that costs when ignored.

#### 3A.4.1 The mass bound

> **Proposition 2 (radial budget).** Let \(\Xi\) be convex on \(\mathbb R^n\),
> \(n=Jp_{\mathrm{re}}\), with a unique minimum \(\Xi(\beta^\dagger)=0\), and let
> \(\widetilde B_r=\{\Xi\le r\}\). Then
>
> \[
> \pi\bigl(\widetilde B_r^{\,c}\mid y\bigr)
> \ \le\
> \frac{\Gamma(n,r)}{\gamma(n,r)}
> \ =\
> \frac{\Pr\bigl(\mathrm{Gamma}(n,1)>r\bigr)}
>      {\Pr\bigl(\mathrm{Gamma}(n,1)\le r\bigr)} .
> \tag{P2}
> \]

*Proof.* In polar coordinates at \(\beta^\dagger\), for any direction \(u\) on the
unit sphere, the right half of (R) gives

\[
\int_{t_r}^{\infty}e^{-\varphi_u(t)}t^{n-1}\,dt
\ \le\
\int_{t_r}^{\infty}e^{-rt/t_r}t^{n-1}\,dt
=\Bigl(\tfrac{t_r}{r}\Bigr)^{\!n}\Gamma(n,r),
\]

and the left half gives
\(\int_0^{t_r}e^{-\varphi_u}t^{n-1}dt\ge(t_r/r)^n\gamma(n,r)\). Integrating both
over the sphere, the factor \(\int(t_r(u)/r)^n\sigma(du)\) is **common to
numerator and denominator and cancels**. Bounding the normalizer below by its
restriction to \(\widetilde B_r\) gives (P2). \(\blacksquare\)

The content of the cancellation is worth stating plainly: **the bound is free of
the geometry of \(\widetilde B_r\)**. No support function, no curvature, no
determinant, no eigenproblem — only the dimension and the level. That is exactly
what the failed Gaussian-majorant attempt could not achieve, and it is a direct
consequence of taking level sets of the log-density rather than any other shape.

The elementary argument above gives the ratio form; the sharper
\(\Gamma(n,r)/\Gamma(n)\) is standard for log-concave densities, and for \(r>n\)
the two differ negligibly (at \(n=20,r=30\): \(0.0223\) versus \(0.0218\)).

**Verified.** `data-raw/_chk_logit_levelset_bounds.R` Check 7 confirms (P2)
against \(2\times10^5\) exact iid draws from \(\pi(\beta\mid y)\) — obtained
without a sampler, by drawing \(\gamma\) from its grid posterior and then the
\(\beta_j\) independently from their conditionals — and separately verifies the
two chord inequalities (R) in 300 random directions with zero violations.

| \(r\) | exact mass | (P2) | old Gaussian-majorant bound |
|---|---|---|---|
| 20 | \(5.3\times10^{-3}\) | \(0.89\) | \(1\) |
| 30 | \(1.0\times10^{-5}\) | \(2.2\times10^{-2}\) | \(1\) |
| 40 | \(<10^{-5}\) | \(1.8\times10^{-4}\) | \(1\) |
| 50 | \(<10^{-5}\) | \(4.8\times10^{-7}\) | \(1\) |

#### 3A.4.2 The floor

The same rays give the floor, and give it *exactly* rather than through a
relaxation. The weight \(w_{j,i}\) depends on \(\beta\) only through the linear
functional \(\eta_{j,i}\), and for logit \(w=p(1-p)\) is **quasi-concave in
\(\eta\)**, peaking at \(\eta=0\). Hence along any line through the mode the
precision rises and then declines, and on a convex set the minimum of \(w\) is
attained where \(|\eta|\) is largest — that is, on \(\partial\widetilde B_r\),
with

\[
\min_{\beta\in\widetilde B_r}w\bigl(\eta_{j,i}(\beta)\bigr)
=\min\Bigl\{w\bigl(\eta^-_{j,i}\bigr),\,w\bigl(\eta^+_{j,i}\bigr)\Bigr\},
\qquad
\eta^{\pm}_{j,i}=\max_{\beta\in\widetilde B_r}\bigl(\pm\eta_{j,i}(\beta)\bigr).
\]

So the floor needs only the **support function of \(\widetilde B_r\)** in the
\(2\sum_jn_j\) directions \(\pm d_{j,i}\), each a one-dimensional root-find along
the corresponding ray. This is where "large enough that the precision is
declining in both directions" earns its keep: it is precisely the condition that
makes the endpoint evaluation correct, and quasi-concavity of \(w\) supplies it
for free on any convex set.

Relaxing the support function by an ellipsoid is possible but, done naively,
disastrous. From the Hessian

\[
\nabla^2\Xi(\beta)
=\Lambda_\beta+\mathrm{blockdiag}\bigl(Z_j^\top W_j(\beta_j)Z_j\bigr),
\]

and for **logit** the weights obey \(0\le w_{j,i}=n_{j,i}p(1-p)\le n_{j,i}/4\)
*globally*. Hence, with
\(\Lambda^{\mathrm{UB}}:=\Lambda_\beta+\tfrac14\,\mathrm{blockdiag}
\bigl(Z_j^\top N_jZ_j\bigr)\) and \(N_j=\mathrm{diag}(n_{j,i})\),

\[
\Lambda_\beta\ \preceq\ \nabla^2\Xi(\beta)\ \preceq\ \Lambda^{\mathrm{UB}}
\qquad\text{for every }\beta .
\tag{S}
\]

Both bounds are explicit matrices, requiring no evaluation of \(\beta^\dagger\)
beyond its location. Integrating (S) twice from the mode, where
\(\nabla\Xi(\beta^\dagger)=0\),

\[
\tfrac12\|\beta-\beta^\dagger\|^2_{\Lambda_\beta}
\ \le\ \Xi(\beta)\ \le\
\tfrac12\|\beta-\beta^\dagger\|^2_{\Lambda^{\mathrm{UB}}} .
\tag{S'}
\]

**Outer ellipsoid, hence the floor.** The left half of (S') gives
\(\widetilde B_r\subseteq\{\|\beta-\beta^\dagger\|^2_{\Lambda_\beta}\le2r\}\), so
by Cauchy–Schwarz in the \(\Lambda_\beta^{-1}\) inner product,

\[
\boxed{\
\eta_{j,i}^{\max}
\ \le\
\bigl|o_{j,i}+d_{j,i}^\top\beta_j^\dagger\bigr|
+\sqrt{2r\;d_{j,i}^\top\bigl(\Lambda_\beta^{-1}\bigr)_{jj}d_{j,i}}\ }
\tag{F}
\]

and \(\underline w_j=\min_i n_{j,i}\,p(\eta^{\max}_{j,i})
\bigl(1-p(\eta^{\max}_{j,i})\bigr)\). No root-finding is needed — but (F) uses
the *prior* Schur precision \(\Lambda_\beta\), discarding all data curvature, and
Check 2 of the verification script measures the damage: against the true support
function, (F) is loose by a factor of **4 to 20**, returning
\(\eta^{\max}\approx10\)–\(30\) where the truth is \(0.6\)–\(5.3\). The implied
weight floor is \(\sim10^{-9}\) and \(\omega_{\max}=1.000\): completely vacuous.

**The bootstrap, and its bistability.** (F) is repairable. A floor \(\Gamma_j\)
valid on \(\widetilde B_r\) upgrades the left side of (S) to
\(\Lambda^{\mathrm{LB}}=\Lambda_\beta+\mathrm{blockdiag}(Z_j^\top\Gamma_jZ_j)\),
shrinking the ellipsoid and raising the floor. The localisation is legitimate:
if \(\nabla^2\Xi\succeq\Lambda^{\mathrm{LB}}\) merely on
\(E=\{\|\cdot\|^2_{\Lambda^{\mathrm{LB}}}\le2r\}\), then any \(\beta\notin E\)
has its segment from \(\beta^\dagger\) leave \(E\) at some \(\bar\beta\) with
\(\Xi(\bar\beta)\ge r\), and \(\Xi\) increases along rays from the mode, so
\(\Xi(\beta)\ge r\) and \(\widetilde B_r\subseteq E\) as required.

The map \(T(\Lambda)=\Lambda_\beta+\mathrm{blockdiag}
\bigl(Z_j^\top\Gamma_j(\Lambda)Z_j\bigr)\) is **monotone increasing**, so by
Knaster–Tarski its fixed points form a lattice — and this matters, because
iterating from \(\Lambda_\beta\) converges to the *least* fixed point, which is
the trivial one \(\Lambda_\beta\) itself whenever the starting ellipsoid is so
wide that the weights underflow. Only the fixed point must satisfy the validity
condition, not the iterates, so one may start from
\(\Lambda_\beta+\tfrac14\mathrm{blockdiag}(Z_j^\top N_jZ_j)\) and descend to the
**greatest** fixed point. Check 4 shows the difference is the whole ballgame:

| \(r\) | true \(\eta^{\max}\) | plain (F) | from below | from above |
|---|---|---|---|---|
| 2 | 2.75 | 10.8 | 2.92 | 2.92 |
| 8 | 3.89 | 19.8 | 19.8 | **4.53** |
| 20 | 5.35 | 30.4 | 30.4 | **6.58** |
| 40 | 7.10 | 42.2 | 42.2 | 42.2 |

From above, the median conservatism falls from \(14.4\) to \(1.17\), and at
\(r=8\) the certified \(\omega_{\max}\) goes from \(1.000\) (vacuous) to
\(0.905\) against a true value of \(0.836\). But the last row is a **phase
transition**, not slow decay: past a critical level — here between \(r=28\) and
\(r=30\) — even the greatest fixed point is trivial and the floor vanishes
outright. This is the Theorem 2 trade-off appearing as a cliff rather than the
gentle \(\sqrt{\log}\) degradation of §3A.6.

#### 3A.4.3 What fails: the Gaussian majorant

The natural-looking companion to (F) is to use the right half of (S') to majorize
the escape mass, giving
\(\pi(\widetilde B_r^{\,c}\mid y)\le(\prod_i\theta_i)^{1/2}
\Pr(\sum_i\theta_iZ_i^2>2r)\) with
\(\theta_i=\mathrm{eig}(\Lambda_\beta^{-1}\Lambda^{\mathrm{UB}})\) — structurally
the twin of §5, with \((\prod\theta_i)^{1/2}\) in the role of \(c'\). It is
**valid and worthless**. In the worked example the constant is
\(1.8\times10^{11}\) and the bound is \(1\) at every level (Check 3).

The failure is structural, not a matter of tuning. The constant is
\(\prod_j(1+\text{data}_j/\text{prior}_j)^{1/2}\), *exponential in
\(Jp_{\mathrm{re}}\)*, because the only globally valid quadratic minorant of
\(\Xi\) is the prior one — logit weights vanish in the tails — while the
normalizer is set by the far more concentrated posterior. Open item 7 worried
that \(c'\) might degrade as \(q\) grows; on the \(\beta\) side the dimension is
\(Jp_{\mathrm{re}}\) rather than \(q\), and it degrades immediately.

Proposition 2 avoids this precisely by not pretending the tail is Gaussian. Past
the boundary the precision declines along every ray, \(\Xi\) grows linearly, and
(R) captures that exactly — which is why (P2) needs no constant at all.

#### 3A.4.4 Group-wise minima, and why the box is not a rival

The floor can be read off the set one group at a time. Since
\(P_{22}(\beta)=\mathrm{blockdiag}\bigl(P_b+G_j(\beta_j)\bigr)\) is **block
diagonal**, the positive-semidefinite ordering decomposes, so

\[
\Gamma_j:=\min_{\beta\in\widetilde B}G_j(\beta_j)
\qquad\Longrightarrow\qquad
P_{22}(\beta)\succeq\mathrm{blockdiag}(P_b+\Gamma_j)\ \ \forall\beta\in\widetilde B ,
\]

and **no single \(\beta\in\widetilde B\) need attain all the minima
simultaneously** — each block sees only its own. When \(p_{\mathrm{re}}=1\) each
\(\Gamma_j\) is an ordinary minimum of a continuous scalar function on a compact
set, hence attained; for \(p_{\mathrm{re}}>1\) it must be replaced by
\(\underline w_jZ_j^\top Z_j\), a valid PSD lower bound that is generally not
attained by any point. \(G_j\) depends on \(\beta\) only through \(\beta_j\), so
\(\Gamma_j\) is a minimum over the **projection** of \(\widetilde B\) onto
coordinate \(j\) — which is the support-function computation of §3A.4.2.

**A correction.** An earlier version of this section reported that a per-group
box beats \(\widetilde B_r\) at matched budget (\(\omega_{\max}=0.801\) versus
\(0.950\) at \(\delta_\beta=0.01\)) and concluded that the shape was wrong. The
comparison was not admissible. That box was sized from the exact per-group
marginals \(\pi(\beta_j\mid y)\), which for \(J>1\) have **no closed form**
(`LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §8) and were obtained in the script only by
grid-integrating \(\gamma\) out — the very computation the certificate exists to
avoid. It is an oracle, not a constructible set.

The constructible box is the **projection box** of the level set,
\(PB:=\prod_j\bigl[\min_{\widetilde B_r}\beta_j,\ \max_{\widetilde B_r}\beta_j\bigr]\).
Since \(PB\supseteq\widetilde B_r\), its budget is inherited from (P2), and since
\(\omega_{\max}\) depends on a set only through its coordinate projections, its
floor is **identical** to the level set's. So the level set and the only box we
can actually build are the same certificate, and the shape question does not
arise. Check 8:

| certified \(\delta_\beta\) | \(r\) | \(\omega_{\max}\) certified | \(\omega_{\max}\) oracle box | exact mass outside \(PB\) |
|---|---|---|---|---|
| \(0.5\) | 21.6 | 0.962 | 0.675 | \(<5\times10^{-6}\) |
| \(0.2\) | 24.3 | 0.970 | 0.710 | \(<5\times10^{-6}\) |
| \(0.05\) | 28.0 | 0.979 | 0.757 | \(<5\times10^{-6}\) |

**Where the remaining gap actually lives.** Against the oracle it splits in two,
and the larger piece is not the shape:

- *Conservatism of (P2)*, \(\omega_{\max}\ 0.866\to0.962\). Certifying
  \(\delta_\beta=0.5\) forces \(r=21.6\) when the true median of \(\Xi\) is
  \(9.7\). The last column shows how far this overshoots: the certified set is so
  large that not one of \(2\times10^5\) exact draws falls outside it.
- *Shape*, \(\omega_{\max}\ 0.675\to0.866\) at the exact level \(r=9.7\). A
  superlevel set of a sum does permit single-coordinate excursions that a product
  set forbids; that effect is real, but it is the smaller half and it is not
  actionable without marginals we do not have.

So the target is the dimension dependence of (P2), which needs \(r\gtrsim n\).
The one-dimensional radial bound is \(\Gamma(1,s)/\gamma(1,s)=e^{-s}/(1-e^{-s})\)
— incomparably better — so a *per-coordinate* budget would close most of the gap.
The obstruction is exactly the missing marginals: a coordinate-wise radial
argument needs the deficiency of \(\pi(\beta_j\mid y)\), and profiling the other
\(n-1\) coordinates out reintroduces a determinant ratio in dimension \(n-1\),
which is the failure of §3A.4.3 again. Open item 6.

#### 3A.4.5 The map \(\beta_j\mapsto G_j(\beta_j)\), and the floor when \(p_{\mathrm{re}}>1\)

Written out, the data precision of group \(j\) is

\[
\eta_{ij}(\beta_j)=o_{ij}+z_{ij}^\top\beta_j,
\qquad
w_{ij}(\beta_j)=n_{ij}\,p_{ij}(1-p_{ij}),
\qquad
G_j(\beta_j)=\sum_i w_{ij}(\beta_j)\,z_{ij}z_{ij}^\top=Z_j^\top W_j(\beta_j)Z_j .
\tag{W}
\]

Two properties of this map decide everything downstream.

**(i) Each weight sees only one linear functional.** \(w_{ij}\) depends on
\(\beta_j\) solely through the scalar \(\eta_{ij}\). So minimising it over a
convex set is a *one-dimensional* question no matter what \(p_{\mathrm{re}}\) is.

**(ii) \(w\) is strictly log-concave in \(\eta\), but *not* concave.** Writing
\(p'=p(1-p)=w\),

\[
w'=w\,(1-2p),
\qquad
w''=w(1-2p)^2-2w^2=w\,(1-6w),
\qquad
(\log w)''=-2w<0 .
\tag{C}
\]

Log-concavity is global. Concavity is not: since \((1-2p)^2=1-4w\), the sign of
\(w''\) is the sign of \(1-6w\), so

\[
w\ \text{is concave}\iff w>\tfrac16
\iff |\eta|<\log(2+\sqrt3)\approx1.3170
\iff p\in(0.211,\,0.789),
\]

and \(w\) is **convex** outside that band — as it must be, since
\(w\approx e^{-|\eta|}\) in the tails. The distinction is recorded because it
decides §3A.4.5's last paragraph.

What the floor needs is weaker than either: only that \(w\) have **no interior
local minimum**, so that

\[
\min_{\eta\in[a,b]}w(\eta)=\min\bigl\{w(a),w(b)\bigr\} .
\tag{EP}
\]

Log-concavity gives quasi-concavity, hence unimodality with peak at \(\eta=0\),
hence (EP) — the minimum sitting at whichever endpoint has larger \(|\eta|\).
(EP) is what §3A.4.2 uses, and it is not specific to logit: for Poisson with the
canonical log link \(w=e^\eta\) is strictly *increasing*, so (EP) holds trivially
with the minimum always at the left endpoint. Probit and cloglog are unimodal
like logit.

Together: the range of \(\eta_{ij}\) over \(\widetilde B_r\) is
\([\,-h(-z_{ij}),\,h(z_{ij})\,]+o_{ij}\) with \(h\) the support function, so

\[
\underline w_{ij}
=\min\Bigl\{w\bigl(\eta^-_{ij}\bigr),\,w\bigl(\eta^+_{ij}\bigr)\Bigr\},
\qquad
\boxed{\ \Gamma_j:=\sum_i \underline w_{ij}\,z_{ij}z_{ij}^\top\ \preceq\ G_j(\beta_j)
\quad\forall\beta\in\widetilde B_r\ }
\tag{FL}
\]

valid termwise, since \(w_{ij}\ge\underline w_{ij}\) and \(z z^\top\succeq0\).
**The weights are still found from \(2n_j\) support-function evaluations
regardless of \(p_{\mathrm{re}}\); only the assembly changes.** \(\Gamma_j\succ0\)
iff \(Z_j\) has full column rank — (H3b) again.

There is no "minimum of \(G_j\)" to take: the positive-semidefinite order is
partial and \(\{G_j(\beta_j):\beta_j\in B_j\}\) has no least element in general.
(FL) is a *construction*, and like the group-wise minima of §3A.4.4 it is attained
by no single \(\beta_j\) — different observations reach their worst weights at
different points. Validity is all that is required, and the same argument that
lets the floor be assembled group by group lets it be assembled observation by
observation within a group.

The cruder alternative \(\Gamma_j^{\mathrm{scalar}}=
\bigl(\min_i\underline w_{ij}\bigr)Z_j^\top Z_j\) is also valid and is what
§3A.2–§3A.3 used, but it applies the single worst weight to every observation.
`data-raw/_chk_logit_matrix_floor_2d.R` measures the cost on a random-intercept-
and-slope logit (\(p_{\mathrm{re}}=2\), \(J=6\), \(r=6\)): both floors pass the
PSD check at 600 points of \(\widetilde B_r\) with zero violations,
\(\Gamma_j\succeq\Gamma_j^{\mathrm{scalar}}\) for every group, and

\[
\omega_{\max}:\quad 0.929\ \text{(scalar)}\ \longrightarrow\ 0.849\ \text{(per-observation)} .
\]

The within-group spread \(\max_i\underline w_{ij}/\min_i\underline w_{ij}\) runs
to \(9.6\) here, and that spread is precisely what the scalar version discards.

**Anisotropy is the new phenomenon.** For \(p_{\mathrm{re}}=1\) the certificate
was set by the worst group; now it is set by the worst *direction* in the worst
group, through \(\Omega_j=(I+P_b^{-1/2}G_jP_b^{-1/2})^{-1}\). The floor can be
**more** anisotropic than the design — in the example \(\mathrm{cond}(\Gamma_j)\)
reaches \(5.2\) against \(\mathrm{cond}(Z_j^\top Z_j)=2.0\) — because observations
differ in how far their own \(\eta_{ij}\) can run inside \(\widetilde B_r\), so
the weights re-tilt the Gram matrix rather than merely rescaling it. A group whose
covariate values are nearly constant has a near-singular \(Z_j^\top Z_j\), hence a
\(\Gamma_j\) with a small eigenvalue and \(\omega_j\approx1\) in that direction,
no matter how many observations it contributes.

**What is left on the table, and when it is reachable.** The tightest
constant-matrix floor is the largest \(\Gamma\) with
\(v^\top\Gamma v\le g_j(v):=\min_{\beta\in B_j}v^\top G_j(\beta)v\) for every unit
\(v\); the condition is sufficient as well as necessary, and in 2-D the
constraints are indexed by a single angle, so the outer problem is a small
semidefinite program. The difficulty is entirely in \(g_j\).

Here (C) matters. Log-concavity is *not* preserved by summation, so
\(v^\top G_jv=\sum_iw_i(\beta)(z_i^\top v)^2\) is in general neither convex nor
quasi-concave, and its minimum need not sit at an extreme point: two weights
centred at \(\eta=\pm2.5\) put both terms in their convex tails near \(\eta=0\)
and produce an interior dip, so (EP) fails for the sum. Concavity, by contrast,
*is* preserved by summation. So on the region where every
\(|\eta_{ij}|\le\log(2+\sqrt3)\), each \(w_i\) is concave, \(v^\top G_jv\) is
concave, and \(g_j\) becomes a **concave minimisation over a convex set — attained
at an extreme point**, hence tractable.

That band is narrow (\(p_{ij}\in(0.211,0.789)\)) and the worked example runs to
\(|\eta|\approx4\text{–}5\), so it does not apply there. But it locates the
difficulty precisely: the directional bound is hard only once some linear
predictor leaves \((-1.317,\,1.317)\), where the objective turns **indefinite** —
concave in the central band, convex in the tails — so neither convex programming
nor "concave minimum at an extreme point" applies.

**What separates this from (FL) is reducibility, not convexity.** Both objectives
are non-convex: \(w\) is log-concave and unimodal, so even
\(\beta\mapsto w(c_{ij}^\top\beta)\) in (FL) is a non-convex function being
minimised over a convex set. Convexity enters only in two *other* places — the
constraint function \(\Xi\), globally strictly convex because
\(\nabla^2\Xi=\Lambda_\beta+\mathrm{blkdiag}(G_j)\succeq\Lambda_\beta\succ0\), and
the support-function subproblems, whose objectives are linear. What makes (FL)
exact is instead the two-step reduction

\[
\min_{\beta\in\widetilde B_r}w(c^\top\beta)
=\min_{\eta\in[\eta^-,\eta^+]}w(\eta)
=\min\{w(\eta^-),w(\eta^+)\},
\]

the first equality because the objective factors through a **single scalar**
linear functional, the second by unimodality (EP). The non-convex minimisation is
discharged by a change of variable and a theorem; the only optimisation performed
is the pair of support functions.

The interval here is the image of the **joint** set,
\([\eta^-,\eta^+]=c^\top\widetilde B_r=\bigcup_{\beta_{-j}}[\ell(\beta_{-j}),
u(\beta_{-j})]\) — an interval because \(\widetilde B_r\) is convex and compact —
so the other groups have already been maximised over, and the first equality is a
tautology (\(\min_{\beta\in K}g(c^\top\beta)=\min_{\eta\in c^\top K}g(\eta)\) for
*any* \(g\)) rather than a relaxation. Nothing is conditioned, so nothing is left
depending on \(\beta_{-j}\); computationally the KKT system
\(\nabla\Xi(\beta)=t\,c\) moves all blocks at once, so \(\Xi^{\mathrm{prof}}_j\)
of (PR) never has to be formed and there is no nested optimisation.

The sum keeps neither leg. It depends on \(\beta_j\) through \(n_j\) linear
functionals at once, so there is no scalar image to reduce to, and summation
destroys quasi-concavity, so interior local minima appear. On top of that, feasibility itself is an optimisation, since
\(B_j\) is described by \(\Xi^{\mathrm{prof}}_j\) of (PR). Decisively: a locally
optimal answer is not merely imprecise but **invalid**, because overestimating
the minimum yields a \(\Gamma\) that is not a floor, with no signal — the unsafe,
silent failure of §3A.4.7. What is required is a *certified global lower bound*,
which is why the discretised SDP (finite angle grid, finite inner search) failed
the PSD check outright: it relaxes in the wrong direction.

**A route that would certify.** The objective is globally Lipschitz with a clean
constant: from (C), \(|w'|=n\,p(1-p)|1-2p|\), and substituting \(q=p-\tfrac12\)
gives \(|w'|/n=2|q|(\tfrac14-q^2)\), maximised at \(q=1/(2\sqrt3)\), so

\[
\|w'\|_\infty=\frac{n}{6\sqrt3}\approx0.0962\,n,
\qquad\text{attained where } p(1-p)=\tfrac16,
\tag{L}
\]

necessarily the same point as the concavity inflection, since \(w''=0\) is the
first-order condition for extremising \(w'\). Hence
\(\beta_j\mapsto\sum_iw_i(\beta)(z_i^\top v)^2\) is Lipschitz with constant
\(L=\tfrac{n}{6\sqrt3}\sum_i\|z_{ij}\|(z_{ij}^\top v)^2\), and branch-and-bound
over \(B_j\) — pruning cells on (centre value) \(-\,L\times\)(radius) — returns a
*certified* lower bound.

**This is available in principle, not cheap.** Certifying a gap \(\epsilon\)
needs cell radius \(\epsilon/L\), hence
\((\mathrm{diam}\cdot L/\epsilon)^{p_{\mathrm{re}}}\) cells — exponential in
\(p_{\mathrm{re}}\), with \(L\) growing in \(n\) — and that is the cost for a
*single* direction \(v\), while the SDP constrains a continuum of them, for every
group. Nor is low dimension a refuge. With
\(\eta_{ij}=o_{ij}+z_{ij}^\top\beta_j\), each term peaks where its own linear
predictor vanishes; only if the offsets \(o_{ij}\) are absent do those peaks
coincide, making every term decreasing in \(|\beta_j|\) and the sum unimodal. Any
fixed-effect contribution or varying exposure scatters them, and a scalar
\(\beta_j\) with \(k\) separated peaks already has \(k-1\) interior local minima
plus two endpoints. Multiple local minima are therefore present at
\(p_{\mathrm{re}}=1\).

**The real dividing line is where a climbing algorithm is *sufficient*.** The
support-function solve is convex — a linear functional over a convex body — so
KKT is sufficient and Newton on \(\nabla\Xi(\beta)=t\,c\) returns the global
answer with a certificate. The inner minimisation is non-convex, so the same
algorithm returns a stationary point with no gap bound, and a wrong basin yields
an invalid \(\Gamma\) silently. Compounding this, \(B_j\) is defined implicitly
through \(\Xi^{\mathrm{prof}}_j\), so each feasibility test is itself an inner
minimisation over \(\beta_{-j}\).

So (FL) is not a cheaper attack on the same problem; it **replaces search by
theorem**. Each \(\underline w_{ij}\) costs two support functions (convex,
KKT-sufficient) and a comparison of two numbers (justified by unimodality), with
no step at which a local optimum could pass for a global one. That, rather than
cost, is why it is the recommendation; the price is the safe slack of
non-simultaneity.

#### 3A.4.6 The feasible set for \(\beta_j\) depends on the other groups

(FL) minimises \(w_{ij}\) over \(\widetilde B_r\), and \(\widetilde B_r\) is a set
in \(\mathbb R^{Jp_{\mathrm{re}}}\), not a product. So although \(\eta_{ij}\)
involves only \(\beta_j\), **the constraint does not**: \(\Xi\) couples the groups
through \(\Lambda_\beta\), which is not block diagonal because integrating
\(\gamma\) out couples them (`LOGIT_MARGINAL_INTEGRATE_GAMMA.md` §4). The correct
object is the **infimal projection**

\[
B_j=\bigl\{\beta_j:\Xi^{\mathrm{prof}}_j(\beta_j)\le r\bigr\},
\qquad
\Xi^{\mathrm{prof}}_j(\beta_j):=\min_{\beta_{-j}}\Xi(\beta_j,\beta_{-j}),
\tag{PR}
\]

convex because partial minimisation of a jointly convex function is convex.
Maximising \(z_{ij}^\top\beta_j\) over \(B_j\) is therefore a constrained
optimisation over the **full** vector, with \(\beta_{-j}\) free to move so as to
spend as little \(\Xi\)-budget as possible.

**Where the other groups enter the support-function solve.** Splitting
\(\nabla\Xi(\beta)=t\,c\) by block, with \(c\) supported on block \(j\),

\[
\nabla_j\Xi(\beta)=t\,z_{ij},
\qquad
\nabla_k\Xi(\beta)=0\quad(k\neq j).
\]

The second set of equations says \(\beta_{-j}\) sits at its **conditional mode
given \(\beta_j\)** — the profile minimiser — so (PR) is computed implicitly and
never has to be formed. The coupling channel is exactly the off-diagonal of
\(\Lambda_\beta\): the likelihood score is block separable
(\(Z_k^\top(y_k-n\,p(\eta_k))\) involves \(\beta_k\) alone), so only
\(\Lambda_\beta(\beta-\mu)\) links blocks, through \((\Lambda_\beta)_{kj}\). That
off-diagonal exists *because* \(\gamma\) was integrated out; conditioning on
\(\gamma\) instead would make \(\Lambda_\beta\) block diagonal and decouple the
problem entirely.

Eliminating \(\beta_{-j}\) leaves the profile Hessian as a Schur complement of
\(H=\nabla^2\Xi\):

\[
\nabla^2\Xi^{\mathrm{prof}}_j
=\underbrace{(\Lambda_\beta)_{jj}+G_j}_{H_{jj}}
-(\Lambda_\beta)_{j,-j}
\Bigl[(\Lambda_\beta)_{-j,-j}+\mathrm{blkdiag}(G_k)\Bigr]^{-1}
(\Lambda_\beta)_{-j,j}
\ \preceq\ H_{jj}.
\tag{SC}
\]

The correction is PSD, so the other groups **reduce** the curvature along
\(\beta_j\): the level set reaches further, \([\eta^-,\eta^+]\) widens, and by
(FL) the floor falls. That is the safe direction, and it is why freezing
\(\beta_{-j}\) — which substitutes \(H_{jj}\) for the Schur complement — overstates
curvature and hence the floor. Since \(S^{-1}=(H^{-1})_{jj}\), the same dependence
appears in the closed-form fallback as \(c^\top H^{-1}c=z_{ij}^\top(H^{-1})_{jj}
z_{ij}\), the **marginal** block of the inverse rather than the inverse of the
block.

**The shortcut to avoid.** Freezing \(\beta_{-j}\) at the mode and using the
*conditional* deficiency \(\Xi(\cdot,\beta^\dagger_{-j})\) gives a smaller set,
since \(\Xi^{\mathrm{prof}}_j\le\Xi(\cdot,\beta^\dagger_{-j})\), hence a larger
floor: anti-conservative, the same failure mode as simulation. Measured on the
\(p_{\mathrm{re}}=2\) example, it reports \(\omega_{\max}=0.816\) against the
correct \(0.849\).

**How much do the other groups actually contribute?** Here, little: profile
extents exceed conditional ones by only \(0.7\%\). The reason is that the
likelihood curvature is block diagonal and dominates the rank-\(q\) prior
coupling, so there is not much for the other groups to do. That is a statement
about this design, not a general licence — and note that a \(0.7\%\) error in
extent still moved \(\omega_{\max}\) by \(4\%\), because \(w\) decays
exponentially in \(|\eta|\) and \(\omega\) is stiff near \(1\).

**Which endpoint is active, and why it does not matter.** Since \(w\) is symmetric
(\(w(-\eta)=w(\eta)\), because \(p(-\eta)=1-p(\eta)\)) and strictly decreasing in
\(|\eta|\), write the certified range as \([\ell,u]\) with midpoint
\(m=(\ell+u)/2\) and half-width \(\rho=(u-\ell)/2\). From
\((m+\rho)^2-(m-\rho)^2=4m\rho\),

\[
\underline w=w\bigl(|m|+\rho\bigr),
\qquad
\arg\min=\begin{cases}u,&m>0\\ \ell,&m<0\end{cases}
\qquad\text{tie iff } m=0 .
\tag{EN}
\]

So the active endpoint is decided by the **sign of the midpoint**. This flips as
the conditioning moves, and the flip is real. In the symmetric two-group scalar
model (\(J=2\), \(y_j=n/2\), prior mean \(0\), \(\Lambda_\beta=aI+b\mathbf 1\mathbf
1^\top\) with \(b=-\lambda_b^2/(\lambda_\gamma+2\lambda_b)<0\)): freezing
\(\beta_1=0\) makes \(\Xi(0,\cdot)\) exactly even, so the slice is \([-\rho,\rho]\)
and **both** endpoints attain the minimum. Moving \(\beta_1=t\) shifts the
conditional mode by \(\beta_2^{*\prime}(0)=-b/(w_0+a+b)>0\), and a cubic expansion
of \(\Xi(t,\beta_2^*+s)-\Xi_{\min}(t)=\tfrac12\kappa_2s^2+\tfrac16\kappa_3s^3\)
with \(\kappa_2=w(\beta_2^*)+a+b\), \(\kappa_3=w'(\beta_2^*)\) displaces *both*
roots equally by \(-\kappa_3\rho_0^2/(6\kappa_2)\), giving

\[
m(t)\approx\beta_2^*(t)-\frac{w'(\beta_2^*(t))\,r_t}{3\kappa_2^2},
\qquad
m'(0)=\beta_2^{*\prime}(0)\Bigl[1-\frac{w''(0)\,r}{3\kappa_2^2}\Bigr]>0,
\]

the bracket exceeding \(1\) because \(w''(0)=w_0(1-6w_0)<0\) by (C) whenever
\(n>2/3\). So \(t>0\) activates \(u\), \(t<0\) activates \(\ell\), and the flip at
\(t=0\) is transversal.

**The certificate is immune because it takes the union.** The static range is the
projection of the *joint* set,
\([\eta^-,\eta^+]=\bigcup_{\beta_{-j}}[\ell(\beta_{-j}),u(\beta_{-j})]\), an
interval because it is the image of a convex compact set under a linear
functional. Hence \(|m|+\rho=\max(|\eta^-|,|\eta^+|)\) is a max over **both**
endpoints and over **all** conditionings: the branches are enveloped, and
\(\underline w_{ij}=\min_{\beta\in\widetilde B_r}w_{ij}(\beta)\) holds with
equality — (FL) is exact per observation, not a relaxation. In the symmetric model
oddness gives \(\eta^-=-\eta^+\), so the static answer sits at the tie.

**What the flip does damage is any *conditional* floor.** Near a tie \(|m|\) grows
linearly in the displacement while \(\rho\) moves only at second order (the
residual budget is \(r_t=r-\tfrac12\tilde\kappa t^2\) and \(w'(0)=0\)), so
\(\beta_{-j}\mapsto w(|m|+\rho)\) has a **downward kink** on the locus
\(\{m_{ij}=0\}\), with gradient jump \(2w'(\rho)\,\nabla m\neq0\). The tie is a
local *maximum* of the conditional floor, so displacing the other groups is pure
loss — which is why the static certificate lands at the extremes. Two corollaries:

- \(\underline w=w(|m|+\rho)\le w(\rho)\), equality iff \(m=0\). **Eccentricity of
  the certified \(\eta\)-range about the origin is pure loss**, so centering that
  symmetrises the ranges about \(0\) maximises the floor. This is the quantitative
  content of the Centering remark in §7.3.3.
- The §3A.4.2 bootstrap iterates through these minima, so its map is a min of
  smooth functions: continuous, with Jacobian discontinuous across each
  \(\{m_{ij}=0\}\). Piecewise-smooth maps are the standard setting for coexisting
  attractors and border-collision transitions, a plausible — but here unverified —
  mechanism for the observed bistability and cliff.

**Where it is not small at all is the ellipsoid relaxation.** (F) sees only
\(\Lambda_\beta\), where the coupling is essentially the whole story, and it uses
the **marginal** block \((\Lambda_\beta^{-1})_{jj}=\Psi+\mathcal W_j
\Lambda_\gamma^{-1}\mathcal W_j^\top\) rather than
\(\bigl((\Lambda_\beta)_{jj}\bigr)^{-1}\approx\Psi\). Those differ by a factor of
\(6.3\) in the \(\eta\) scale in the example — a factor of \(40\) in variance —
because the marginal block carries the population-level uncertainty in
\(\gamma\), which is large when \(\Lambda_\gamma\) is weak. Any Gaussian
shortcut that conditions rather than marginalises is badly wrong in the unsafe
direction; (F) as written is correct.

#### 3A.4.7 Why the minima are not found by simulating from \(\widetilde B\)

It is natural to expect that, in general, the minima in (FL) are hard and must be
approximated by drawing from \(\widetilde B\) and minimising across draws. Two
things need separating.

**The minima in (FL) are not hard.** Because \(w_{ij}\) depends on \(\beta\)
through the single linear functional \(\eta_{ij}\), and \(w\) is monotone in
\(|\eta|\) by (ii), computing \(\underline w_{ij}\) reduces to **maximising a
linear function over a convex set** — the support function \(h(\pm z_{ij})\).
That is a convex program with a unique solution, obtained by a one-dimensional
root-find in the multiplier with a Newton solve inside. The cost is
\(2\sum_jn_j\) such solves, *linear in the data size and independent of
\(Jp_{\mathrm{re}}\)*. Nothing here needs simulating.

What *is* hard is the tighter directional bound of §3A.4.5, whose inner minimum
is over a sum of log-concave functions and so is a genuine global optimisation.
That bound is not used.

**Simulation errs in the unsafe direction, and the assembly slack does not cover
it.** The minimum over draws satisfies \(\hat{\underline w}\ge\underline w\), so
a simulated \(\widehat\Gamma_j\succeq\Gamma_j\) is **too large** — not a lower
bound. Against this stands a genuine conservatism: (FL) lets each observation hit
its own worst weight, which cannot happen simultaneously, so \(\Gamma_j\) sits
strictly below the tightest valid floor. Whether that slack absorbs the
simulation error is decidable, and the verification script tests it
adversarially, at the exact extremal points where a simulated floor should fail:

| draws \(M\) | mean range recovered | \(\widehat{\underline w}/\underline w\) | PSD violations | \(\omega_{\max}\) |
|---|---|---|---|---|
| \(200\) | 0.62 | 1.70 | 58 | 0.724 |
| \(2\,000\) | 0.76 | 1.38 | 48 | 0.765 |
| \(20\,000\) | 0.83 | 1.30 | 45 | 0.784 |
| exact | 1 | 1 | **0** | 0.849 |

The slack does **not** absorb the error. Violations persist at a hundredfold
increase in draws, and the reported \(\omega_{\max}=0.784\) is *better* than the
true \(0.849\) — the certificate would claim a rate it has not earned, which is
the one failure mode that must not be tolerated.

**And more draws will not rescue it.** The obstruction is geometric: recovering
the extent of a convex body in one *fixed* direction from random draws is
exponentially inefficient in dimension. For uniform draws in the unit ball, the
fraction of the true extent recovered is

| \(n\) | 2 | 5 | 10 | 20 | 50 | 100 |
|---|---|---|---|---|---|---|
| \(M=10^3\) | 1.00 | 0.91 | 0.69 | 0.60 | 0.40 | 0.36 |
| \(M=10^5\) | 1.00 | 0.97 | 0.93 | 0.75 | 0.61 | 0.46 |

with the cost of an extra digit growing like \(M^{2/(n+1)}\). Since
\(n=Jp_{\mathrm{re}}\), this is worst exactly in the regime the certificate is
for. Note the contrast with Proposition 2, which is *helped* by nothing and hurt
by dimension only through \(r\gtrsim n\): there the radial cancellation removed
the geometry, whereas here the geometry in a specific direction is the whole
question.

**The rigorous fallback, if a support function is ever unavailable** — a family
whose weight is not monotone in \(|\eta|\), or a set not presented as a smooth
convex sublevel set — is the bootstrapped outer ellipsoid of §3A.4.2, which errs
in the safe direction and is merely loose. Simulation is appropriate only for
*exploring* how much the directional bound would buy, and any number it produces
must be labelled a target rather than a floor.

### 3A.5 Remarks

**(H3b) is indispensable here.** Chapter C05 treats group-wise estimability as
mattering only for quantitative margins and flat-limit propriety. In Theorem 2 it
is what makes \(\Gamma_j\) nonsingular: without full column rank,
\(Z_j^\top Z_j\) is singular, \(\Gamma_j\not\succ0\), and the conclusion
degenerates to the trivial floor. A single rank-deficient group is enough to
destroy the margin *in the directions it fails to identify*, which is the precise
sense in which adding non-estimable groups cannot help.

**(W) holds for the families of interest, and trivially for Gaussian.** Logit,
probit, Poisson-log, cloglog and Gamma-log all have continuous, strictly positive
weight functions. For Gaussian, \(w\equiv1/\phi\) is *constant*, so
\(\underline w_j\) is independent of \(r\), \(\widetilde B_r\) can be taken to be
all of \(\mathbb R^{Jp_{\mathrm{re}}}\), and no restriction is needed —
consistent with the fact that the Gaussian case is exactly the closure case where
none of this machinery is required.

**The two conclusions pull against each other.** This is the substance of the
theorem and should not be hidden by the existential quantifier. Raising \(r\)
makes (1) easier and (2) *worse*: for every family whose weights vanish in the
tails, \(\underline w_j\downarrow0\) as \(r\uparrow\infty\) — visibly so in (F),
where \(\eta^{\max}_{j,i}\) grows like \(\sqrt{2r}\). Hence

\[
\Gamma_j(\delta)\ \longrightarrow\ 0
\qquad\text{as }\delta\downarrow0 ,
\]

so the *quality* of the floor — and with it \(\omega_{\max}\), \(\kappa_0\), and
the exponent \(\rho=1/g_{\min}\) of §6 — degrades as the \(\beta\)-budget is
tightened. Theorem 2 asserts that both conditions hold simultaneously for every
\(\delta\); it does not assert that the resulting certificate stays useful as
\(\delta\to0\). That is a quantitative question, and §3A.6 is the reassuring
part of the answer.

**Convexity is not incidental.** \(\widetilde B_r\) is convex because \(\Xi\) is
(§3A.2), which is what Lemma 1 needs for Brascamp–Lieb to survive truncation
(§7). Any other shape must be convex for the same reason — one of the several
respects in which level sets of a log-concave marginal are the right family
rather than a convenient one.

### 3A.6 How fast does the floor degrade?

The composition is gentler than either factor suggests, because the two rates
are of different orders.

For a log-concave posterior with sub-Gaussian tails — which is what (H2) plus a
proper \(P_b\) delivers — \(\pi(\|\beta_j\|>b\mid y)\lesssim e^{-cb^2}\), so
meeting a budget \(\delta\) needs only

\[
b\ \sim\ \sqrt{\log(1/\delta)/c}\, .
\]

For logit, \(w(\eta)=p(\eta)(1-p(\eta))\approx e^{-|\eta|}\) for large
\(|\eta|\), so \(\underline w\sim e^{-b}\) and therefore

\[
\underline w(\delta)\ \sim\ \exp\Bigl(-C\sqrt{\log(1/\delta)}\Bigr),
\]

decaying in \(\log(1/\delta)\) only at the square-root rate rather than
geometrically. Concretely, at \(\delta=10^{-6}\) one has
\(\sqrt{\log(1/\delta)}\approx3.7\), so the weight floor is suppressed by a
factor of order \(e^{-3.7}\approx0.025\) rather than by \(10^{-6}\). With
\(N_j=50\) trials and \(\lambda_b=1\) that leaves \(g\approx1.2\) and
\(\rho\approx0.8\) — still a usable exponent. §9's check 9 tabulates this
trade-off.

The structural reason is worth stating: **the posterior tails are Gaussian while
the weights decay only exponentially**, so the box need grow only like
\(\sqrt{\log(1/\delta)}\) while the penalty is exponential in the box radius. Had
the weights decayed like \(e^{-\eta^2}\), the composition would be polynomial in
\(\delta\) and this route would be far weaker.

This is an *asymptotic* argument. §3A.7 measures the actual degradation on a
worked example and finds it mild, but also finds that over the practically
relevant range of \(\delta\) the data do not discriminate between this form and a
weak power law. The magnitude is what should be relied on, not the functional
form.

### 3A.7 A worked logit example

`data-raw/_ex_logit_majorization_floor.R` builds a binomial random-intercept
GLMM with \(J=20\) groups, \(N_j\in\{20,60,200\}\) trials, \(\tau^2=0.5\)
(\(\lambda_b=2\)), and a weak population prior \(\lambda_\gamma=0.05\). Since
\(q=1\) and \(p_{\mathrm{re}}=1\), every quantity is scalar and everything is
computed by direct 1-D quadrature — no sampler, no Laplace — so the *true*
\(\kappa\) is exact and the capped \(\kappa^{\mathrm{UB}}\) can be compared
against it rather than against an approximation.

The group data precision is the one asked for above,
\(G_j(\beta_j)=N_j\,p(\beta_j)\bigl(1-p(\beta_j)\bigr)\).

**The example uses the box construction, not \(\widetilde B_r\).** It takes
per-group equal-tailed central intervals of the marginal
\(\pi(\beta_j\mid y)\), each carrying mass \(1-\delta_\beta/J\), rather than the
level set (B). Theorem 2 holds for either — the proof needs only compactness and
tightness — and the box is what the script implements. The numbers below are
therefore an *upper* bound on the cost of the restriction: §3A.8 shows the box is
demonstrably the wrong shape here, so a level-set version would report smaller
\(\kappa^{\mathrm{UB}}\) at the same \(\delta_\beta\). Re-running the example on
\(\widetilde B_r\) is open item 6.

**The trade-off, measured.** At \(\gamma^\star=-1.254\) the true values are
\(\kappa=0.196\), \(\rho=0.244\).

| \(\delta_\beta\) | box half-width | \(\min_j\underline w_j\) | \(\min_j\Gamma_j\) | \(\omega_{\max}\) | \(\kappa^{\mathrm{UB}}\) | \(\rho\) | \(c'\) |
|---|---|---|---|---|---|---|---|
| \(0.5\) | 0.66 | 0.048 | 0.97 | 0.674 | 0.281 | 0.391 | 1.18 |
| \(10^{-2}\) | 1.04 | 0.025 | 0.50 | 0.799 | 0.345 | 0.527 | 1.24 |
| \(10^{-3}\) | 1.21 | 0.018 | 0.36 | 0.846 | 0.376 | 0.604 | 1.27 |
| \(10^{-4}\) | 1.37 | 0.013 | 0.27 | 0.882 | 0.404 | 0.679 | 1.30 |
| \(10^{-6}\) | 1.66 | 0.008 | 0.16 | 0.927 | 0.452 | 0.826 | 1.35 |

Three readings.

**The cap is affordable.** Driving \(\delta_\beta\) from \(0.5\) to \(10^{-6}\)
— six orders of magnitude — only moves \(\rho\) from \(0.391\) to \(0.826\), a
factor of \(2.1\). Compare the population-prior route, where the package default
\(\mathtt{pop.pwt}=0.01\) gives \(\rho=99\). Capping the linear predictor is
cheap; relying on the prior is not.

**The floor survives with room to spare.** \(\Pi^{\mathrm{LB}}\succ0\) at every
budget, and \(c'\) stays between \(1.18\) and \(1.35\) — so for \(q=1\) the
majorization constant is essentially free and the bound is close to the exact
weighted-\(\chi^2\) tail.

**The cap costs a factor of about two in \(\kappa\), not an order of magnitude.**
True \(\kappa=0.196\) against \(\kappa^{\mathrm{UB}}=0.28\)–\(0.45\). That gap is
the entire price of replacing exact conditional variances by the worst-case
weight inside the box.

### 3A.8 What actually needs capping

The per-group table at \(\delta_\beta=10^{-3}\) makes the binding constraint
plain. Pooling weights \(\omega_j=\lambda_b/(\lambda_b+\Gamma_j)\) range over
more than an order of magnitude:

| \(j\) | \(N_j\) | \(\hat p_j\) | box | \(\underline w_j\) | \(\Gamma_j\) | \(\omega_j\) |
|---|---|---|---|---|---|---|
| 3 | 20 | 0.10 | \([-3.97,\,-0.03]\) | 0.018 | 0.36 | **0.846** |
| 1 | 20 | 0.15 | \([-3.66,\,0.14]\) | 0.024 | 0.49 | 0.804 |
| 7 | 60 | 0.50 | \([-1.17,\,0.85]\) | 0.181 | 10.8 | 0.156 |
| 18 | 200 | 0.38 | \([-1.12,\,0.05]\) | 0.185 | 37.1 | 0.051 |

Since \(\kappa^{\mathrm{UB}}\) is driven by \(\omega_{\max}\), **the certificate
is set by the single least informative group, not by the average or the total**.
Here that is group 3: \(N_j=20\) trials at \(\hat p=0.10\), whose box must stretch
to \(\beta\approx-4\) and whose weight floor is therefore \(0.018\) rather than
the \(0.25\) available near \(p=1/2\).

Two design consequences follow.

**Small \(N_j\) and extreme \(\hat p_j\) compound.** The floor is
\(N_j\underline w_j\), and a small group is also a *wide* group — its marginal is
diffuse, so its box is long, so \(\underline w_j\) is small at both ends. The two
effects multiply rather than offset. A group with \(N_j=200\) at
\(\hat p=0.38\) contributes \(\Gamma_j=37\); one with \(N_j=20\) at
\(\hat p=0.10\) contributes \(0.36\), a hundredfold difference.

**Capping asymmetrically is worth it.** Group 3's box \([-3.97,-0.03]\) is
dominated by its lower end: \(p(1-p)\) at \(-3.97\) is \(0.018\), at \(-0.03\) it
is \(0.250\). The floor is set entirely by the left endpoint, so shrinking only
that side buys nearly all the available improvement at a fraction of the escape
cost. Per-group asymmetric margins — which
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0.4 already permits — are therefore the right
shape for the box, not a common symmetric \(b^\star\).

The practical rule this suggests: **cap on the linear-predictor scale, per group,
asymmetrically, and budget \(\delta_\beta\) unequally across groups** — spending
more of the escape budget on the groups whose weights are already healthy, and
keeping the boxes tight where \(p\) is extreme. The example does none of this
(it splits \(\delta_\beta/J\) equally), so the numbers above are a floor on what
a tuned version would achieve.

**The level set does all of it, and is the only constructible option.** The
asymmetry comes free because \(\Xi\) is skewed whenever \(\hat p_j\ne1/2\), and
the cross-group allocation comes free because the set is not a product — both
from the single scalar \(r\), with no per-group tuning. The tuned box that would
beat it on \(\omega_{\max}\) requires per-group marginals that do not exist in
closed form for \(J>1\), so it cannot be built without the sampling this whole
development is designed to avoid; §3A.4.4 works through the comparison and its
projection-box resolution.

### 3A.9 What Theorem 2 buys

Combined with Proposition 1 and Lemmas 1–3, every ingredient of the escape bound
is now certified rather than assumed:

\[
\pi_\gamma\bigl(\widetilde C_d^{\,c}\bigr)\ \le\ M\,
\Pr\Bigl(\sum_i\nu_iZ_i^2>2d\Bigr),
\qquad
\nu_i=\frac{\kappa_i^{\mathrm{UB}}}{1-\kappa_i^{\mathrm{UB}}},
\]

with \(\kappa^{\mathrm{UB}}\) built from \(P_{22}^{\mathrm{LB}}(\delta_\beta)\)
and \(M\le c'\) computable from eigenvalues alone. The remaining budget
arithmetic — how to split a total tolerance between \(\delta_\gamma\),
\(\delta_\beta\) and the geometric term \((1-\varepsilon)^n\) — is §7 and §10(4),
and is not settled here.

---

## 4. The constant \(c\), and how to avoid needing it

\(c=\pi_\gamma(\gamma^\star)/q_2(\gamma^\star)\) involves the posterior
normalizing constant. Write \(I:=\int e^{-\bar\Phi(\gamma)}\,d\gamma\), so that

\[
c=\frac{(2\pi)^{q/2}\det(\Pi^{\mathrm{LB}})^{-1/2}}{I} .
\]

Two remarks make this cheaper than "the marginal likelihood is required".

**\(I\) is a \(q\)-dimensional integral of a computable integrand.** \(\bar\Phi\)
is available pointwise from the ray quadrature of
`CHAPTER_C05_IMPLEMENTATION.md` §7, which needs only conditional means from the
E-step. The \(\beta\) dimension never enters, and no Laplace approximation of
the full marginal is required. Cost scales with \(q\), not with
\(Jp_{\mathrm{re}}\).

**\(c\) is sandwiched for free.** The two curvature bounds on \(\Phi\) give
two-sided control of \(I\):

\[
\underbrace{\nabla^2\Phi\succeq\Pi^{\mathrm{LB}}}_{\text{Lemma 2}}
\ \Longrightarrow\ I\le(2\pi)^{q/2}\det(\Pi^{\mathrm{LB}})^{-1/2},
\qquad
\underbrace{\nabla^2\Phi\preceq P_{11}}_{S\succeq0,\ \text{always}}
\ \Longrightarrow\ I\ge(2\pi)^{q/2}\det(P_{11})^{-1/2},
\]

hence

\[
\boxed{\ 1\ \le\ c\ \le\ c':=\sqrt{\frac{\det P_{11}}{\det\Pi^{\mathrm{LB}}}}
=\prod_i\bigl(1-\kappa_i^{\mathrm{UB}}\bigr)^{-1/2},
\qquad
\kappa_i^{\mathrm{UB}}:=\mathrm{eig}\bigl(P_{11}^{-1}B^{\mathrm{UB}}\bigr).\ }
\]

The lower end \(c\ge1\) also follows from \(1=\int\pi_\gamma\le c\int q_2=c\).
So an implementation can run the certificate immediately at \(c'\), which needs
only eigenvalues already computed, and refine toward \(c\) by quadrature only if
the resulting sweep count is unacceptable.

**\(c'\) is not a new object.** It equals \(M_Q(1)=\mathbb E_Q[e^{\Psi}]\) from
the tilting identity of `CHAPTER_C05_IMPLEMENTATION.md` §4B.0, evaluated at the
upper-bound spectrum. The two routes are therefore the same construction: §4B
applies Chernoff and pays a polynomial factor to remove the tail shape, while
Proposition 1 keeps the tail shape and pays \(M\). §8 compares them.

---

## 5. The computable form

\(Q_2(\widetilde C_d^{\,c})\) is a Gaussian measure of the complement of a
convex set, which has no closed form. Lemma 3's second conclusion supplies the
reduction:

\[
\widetilde C_d^{\,c}\ \subseteq\ \bigl\{\|\gamma-\gamma^\star\|^2_{B^{\mathrm{UB}}}>2d\bigr\},
\]

and under \(Q_2\) that quadratic form is a weighted \(\chi^2_1\) sum with weights
\(\mathrm{eig}\bigl((\Pi^{\mathrm{LB}})^{-1}B^{\mathrm{UB}}\bigr)\). Since
\(B^{\mathrm{UB}}=P_{11}-\Pi^{\mathrm{LB}}\), those weights are

\[
\nu_i=\frac{\kappa_i^{\mathrm{UB}}}{1-\kappa_i^{\mathrm{UB}}},
\]

the same \(\kappa\mapsto\kappa/(1-\kappa)\) transform as
`CHAPTER_C05_IMPLEMENTATION.md` §4A.0. Combining with Proposition 1:

\[
\boxed{\ \pi_\gamma\bigl(\widetilde C_d^{\,c}\bigr)\ \le\
M\;\Pr\Bigl(\sum_i\nu_iZ_i^2>2d\Bigr),
\qquad M\in\{M^\star,\,c,\,c'\}\ \text{in increasing conservatism}.\ }
\]

This is §4A.3's *exact* closure formula promoted to a bound outside closure: same
shape, with \(\kappa_i\) replaced by upper bounds and one multiplicative
constant. Evaluate the tail by Davies or Imhof, as route 4(a) already does.

### 5.1 Inverting for the budget

\[
d(\delta)=\tfrac12\,Q_\nu\bigl(1-\delta/M\bigr),
\qquad Q_\nu=\text{weighted-}\chi^2\text{ quantile},
\]

by bisection. **The constant enters only through \(\delta\mapsto\delta/M\)**, and
since the tail decays like \(e^{-d/\nu_{\max}}\),

\[
d(\delta)\ \approx\ \nu_{\max}\bigl[\log(1/\delta)+\log M\bigr] .
\]

So a loose \(M\) costs an *additive* shift \(\nu_{\max}\log M\) in \(d\) — the
same structural form as the \(t_\delta/2\) shift of §4A.4, and the reason \(c'\)
is usable even when it is orders of magnitude above \(M^\star\).

There is **no side condition**. Because \(M\ge1\) and \(\delta<1\) always leave
\(\delta/M<1\), the inversion is solvable for every budget, unlike §4B.3's
vacuity threshold \(d>q\kappa_0/(2(1-\kappa_0))\).

### 5.2 The certified rate

Feeding \(d(\delta)\) into \(\varepsilon=\varepsilon_d\,Q(\widetilde C_d)\) with
\(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\):

\[
\varepsilon\ \gtrsim\ \varepsilon(\gamma^\star)\,(\delta/M)^{\nu_{\max}},
\qquad
n\ \approx\ \frac{\log\bigl(1/(\mathrm{tol}-\delta)\bigr)}{\varepsilon(\gamma^\star)}
\Bigl(\frac{M}{\delta}\Bigr)^{\nu_{\max}} .
\]

The polynomial exponent is \(\rho=\nu_{\max}=\mu_{\max}\), one unit better than
Remark 6.1's \(\rho=1+\mu_{\max}\), for the reason identified in §4B.2: the
exponent here is the exact metric-\(B\) rate rather than the \(P_{11}\)-envelope
rate.

---

## 6. Where \(P_{22}^{\mathrm{LB}}\) comes from: pooling weights

\(S\) is the pooling-weighted \(P_{11}^{\mathrm{RE}}\). With
\(G_j:=Z_j^\top W_jZ_j\) and

\[
\Omega_j:=P_b^{1/2}(P_b+G_j)^{-1}P_b^{1/2}=\bigl(I+P_b^{-1/2}G_jP_b^{-1/2}\bigr)^{-1},
\]

one has

\[
S=\sum_jH_j^\top P_b^{1/2}\,\Omega_j\,P_b^{1/2}H_j,
\qquad
P_{11}^{\mathrm{RE}}=\sum_jH_j^\top P_b^{1/2}\,I\,P_b^{1/2}H_j ,
\]

so the trivial \(\Omega_j\preceq I\) is exactly what yields
\(S\preceq P_{11}^{\mathrm{RE}}\) and hence \(\kappa_0=1-\mathtt{pop.pwt}\). In
the scalar-RE case \(\Omega_j\) reduces to the shrinkage weight
\(\omega_j=\lambda_b/(\lambda_b+W_j)\) of `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md`
§0.2 — the "pooling weight" in the usual sense.

A floor \(W_j\succeq\) something positive therefore buys

\[
\kappa_0\ \le\ \omega_{\max}\bigl(1-\mathtt{pop.pwt}\bigr),
\qquad
\omega_{\max}=\max_j\frac{1}{1+g_j},
\quad
g_j:=\lambda_{\min}\bigl(P_b^{-1/2}G_jP_b^{-1/2}\bigr),
\]

and in the flat-prior limit

\[
\rho=\nu_{\max}\ \le\ \frac{\omega_{\max}}{1-\omega_{\max}}=\frac{1}{g_{\min}} :
\]

the reciprocal of the worst group's information relative to the random-effect
prior precision. This is the quantity that decides whether the certificate is
usable, and it is prior-free.

**(H3b) becomes load-bearing.** \(g_j>0\) requires \(Z_j\) full column rank.
Chapter C05 treats group-wise estimability as mattering only for quantitative
margins and flat-limit propriety; on this route it is what makes the exponent
nonzero at all.

**Why this matters more than the prior route.** `CHAPTER_C05_IMPLEMENTATION.md`
§4B.9 records that at the package default \(\mathtt{pop.pwt}=0.01\) one gets
\(\kappa_0=0.99\), hence \(\rho=99\) and a numerically vacuous certificate. A
logit group with \(N_j=50\), \(b^\star=1\), \(\lambda_b=1\) gives
\(g_j\approx9.8\), hence \(\rho\approx0.10\) and \(\delta^\rho\approx0.63\) at
\(\delta=0.01\) — essentially no polynomial degradation. The pooling-weight
floor is the difference between a theorem and a number.

---

## 7. What the \(\beta\)-restriction costs

A global \(P_{22}^{\mathrm{LB}}\succ0\) **does not exist** for logit, probit,
Poisson-log or cloglog: the GLM weights vanish as \(|\eta|\to\infty\), so
\(\inf_\beta W_j(\beta)=0\). The floor is therefore available only on a
restricted \(\beta\)-region, and that has three consequences.

**A second budget.** The certificate covers the sweep restricted to
\(\widetilde C_d\times B\) — Gibbs for the doubly truncated target, by the same
argument as Lemma 24 — so the escape budget splits into \(\delta_\gamma\) for
\(\widetilde C_d\) and \(\delta_\beta=\pi(\beta\notin B)\) for the
\(\beta\)-block. The second is what
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §1–§6 bounds per likelihood, and
`SAFE_UNSAFE_TV_DECOMPOSITION.md` assembles into a three-term total variation
bound. The two lines of work in this package meet here.

**\(B\) must be convex.** Lemma 1 needs Brascamp–Lieb for the *truncated*
conditional. Truncating a log-concave density to a convex set preserves both
log-concavity and the Hessian floor, so boxes and ellipsoids are admissible and
nothing else obviously is. The product-of-intervals construction of
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0.4 already satisfies this.

**Chapter C05's framing changes.** The chapter's selling point is a *static,
drift-free* certificate that restricts \(\gamma\) only. Restricting \(\beta\)
reintroduces the safe/unsafe machinery the chapter was built to avoid. That is a
genuine trade, not a refinement: a usable exponent in exchange for an extra
term and an extra budget. Any write-up should say so plainly.

Per-likelihood floors are already tabulated: logit
\(\kappa_j=N_jp(b^\star)(1-p(b^\star))\), Poisson \(n_js_j\), and the probit,
cloglog and Gamma analogues, in `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` §0.4 and
§1–§6.

---

## 8. Comparison with the tilting route

`CHAPTER_C05_IMPLEMENTATION.md` §4B bounds the same quantity by exponentially
tilting \(Q\); the relationship is close.

| | tilting (§4B) | majorization (here) |
|---|---|---|
| identity used | \(\pi_\gamma/q_Q\propto e^{\Psi}\) | \(\pi_\gamma\le Mq_2\) on \(C^c\) |
| reference Gaussian | \(Q=N(\gamma^\star,P_{11}^{-1})\) | \(Q_2=N(\gamma^\star,(\Pi^{\mathrm{LB}})^{-1})\) |
| tail shape | Chernoff, exponential | exact weighted \(\chi^2\) |
| constant | \(M_Q(1)\ge1\), i.e. \(c'\) | \(M^\star\le c\le c'\) |
| exact in tight limit | no | yes |
| side condition | \(d>q\kappa_0/(2(1-\kappa_0))\) | none |
| needs \(\beta\)-restriction | no (works from \(\Lambda_\gamma\)) | yes, for \(\Pi^{\mathrm{LB}}\succ0\) |

They use the same constant in the worst case, since \(c'=M_Q(1)\) (§4). The
majorization route is tighter — it keeps the tail shape the Chernoff step
discards — and it is exact when the floor is tight. Indicative numbers at
\(q=3,\kappa=0.5\): constant \(2^{3/2}\approx2.8\) here against §4B.9's measured
conservatism of about \(10\); at \(q=8,\kappa=0.95\), \(1.6\times10^5\) against
\(5.9\times10^5\).

The honest division of labour: **§4B where only a population prior is
available**, since it needs no \(\beta\)-restriction and degrades gracefully;
**this note where a pooling-weight floor can be certified**, since it is both
tighter and the only one of the two whose exponent is usable at weak
\(\mathtt{pop.pwt}\).

---

## 9. Numerical verification

`data-raw/_chk_majorization_envelope.R` checks the identities of §2, §4 and §5
against Gaussian-closure cases where every quantity has an independently known
value. The floor is made loose by a single scalar,
\(B^{\mathrm{UB}}=B+t(P_{11}-B)\) for \(t\in[0,1)\), which gives
\(\Pi^{\mathrm{LB}}=(1-t)\Pi\) and \(t=0\) as the tight case. Weighted-\(\chi^2\)
tails use Monte Carlo with common random numbers across the bound and the exact
probability, so the tight case must agree identically rather than up to
simulation error.

Over 70 cases with \(q\in\{2,\dots,5\}\):

| claim | result |
|---|---|
| \(\nu_i=\kappa_i^{\mathrm{UB}}/(1-\kappa_i^{\mathrm{UB}})\) | \(6\times10^{-12}\) |
| \(c'\): determinant form \(=\) product form | \(3\times10^{-14}\) relative |
| sandwich \(1\le c\le c'\) | holds, all cases |
| \(r_d\ge2d\) | holds, all cases |
| \(M^\star=c=1\) when \(t=0\) | holds, all tight cases |
| bound \(=\) exact when \(t=0\) | agree to \(4\times10^{-14}\) |
| bound dominates exact | 70/70, no violations |
| \(M^\star=c\,e^{d-r_d/2}\) is the sup over \(\widetilde C_d^{\,c}\) | matches a \(q=2\) grid scan to grid resolution (\(3\times10^{-4}\) to \(8\times10^{-3}\)) |
| \(c'=M_Q(1)=\mathbb E_Q[e^{\Psi}]\) | \(2\times10^{-2}\), Monte Carlo limited |

Observed magnitudes from the same run, which matter for how the routes should be
presented:

- \(M^\star/c\in[1.4\times10^{-3},\,0.95]\) — the escape-set restriction of §2
  recovers up to three orders of magnitude over the global constant.
- \(c'/c\in[1.4,\,14.8]\) — the price of taking the free constant instead of
  doing the \(q\)-dimensional quadrature for \(c\).
- bound/exact \(\in[1.16,\,276]\) at \(d=3\) across the loose cases, so the route
  is within an order of magnitude of exact when the floor is decent and degrades
  gracefully when it is not.

Separately, over 200 random hierarchical configurations with \(q\le4\),
\(p_{\mathrm{re}}\le3\), \(J\le6\), for §6:

| claim | result |
|---|---|
| \(S=\sum_jH_j^\top P_b^{1/2}\Omega_jP_b^{1/2}H_j\) | \(5\times10^{-15}\) relative |
| \(\kappa_0\le\omega_{\max}\,\lambda_{\max}(P_{11}^{-1}P_{11}^{\mathrm{RE}})\) | holds, 200/200; worst slack \(5\times10^{-4}\) |

The near-zero worst-case slack in the second row says the \(\omega_{\max}\)
factorization is close to tight, so §6's exponent \(1/g_{\min}\) is not a lossy
detour — it is very nearly \(\kappa_0\) itself.

Two caveats on the checks themselves. The \(c'=M_Q(1)\) agreement is limited by
Monte Carlo, not by the identity: \(e^{\Psi}\) has infinite variance under \(Q\)
once \(\kappa>1/2\), the same degeneracy recorded in
`CHAPTER_C05_IMPLEMENTATION.md` §4B.5. And everything here is a closure check —
it validates the algebra, not Lemmas 1–3, which are what carry the argument
outside closure and are checkable only by proof.

---

## 10. Open items

1. **Proposition 1 is trivial; the construction is not.** Lemmas 1–3 chain
   Brascamp–Lieb with strong convexity in a way that is standard step by step
   but has not been reviewed as a whole. In particular Lemma 1's
   curvature-to-covariance step under truncation deserves a written proof.
2. **\(r_d\) is a convex program, not a formula.** §2's sharp constant needs
   \(\min_{\partial\widetilde C_d}\|\cdot\|^2_{B^{\mathrm{UB}}}\). Whether a
   useful closed-form lower bound on \(r_d\) exists — which is all that is
   needed, since \(r_d\) appears with a negative sign — is open.
3. **\(\varepsilon(\gamma^\star)\) is untouched.** This note addresses only the
   set-sizing stage. Outside Gaussian closure there is still no general route to
   \(\varepsilon(\gamma^\star)\), which after this note is the only remaining
   stage of the pipeline without one.
4. **The two-budget assembly is not written down.** §7 identifies the pieces —
   \(\delta_\gamma\) here, \(\delta_\beta\) from
   `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md`, three-term structure from
   `SAFE_UNSAFE_TV_DECOMPOSITION.md` — but not the combined theorem, and the
   interaction between truncating \(\beta\) and the tilting identity's target
   needs care.
5. **The target changes.** With both blocks restricted, the certificate concerns
   convergence to \(\pi\) truncated to \(\widetilde C_d\times B\), not \(\pi\).
   Chapter C05 §7.2 ("What is and is not certified") is an empty heading; this
   note enlarges what belongs in it.
6. **(P2)'s dimension dependence is the binding constraint.** It needs
   \(r\gtrsim n=Jp_{\mathrm{re}}\), and §3A.4.4 shows this, not the shape of the
   set, is the larger source of slack — the certified set is so oversized that no
   draw out of \(2\times10^5\) lands outside it. A per-coordinate budget would
   behave like \(n=1\) and close most of the gap, but requires the deficiency of
   \(\pi(\beta_j\mid y)\), and profiling out the other coordinates reintroduces
   the dimension-\((n-1)\) determinant ratio that sinks §3A.4.3. Whether the
   sharper log-concave form of (P2), a curvature-aware refinement of (R) inside
   the boundary, or a genuinely per-coordinate argument is the way through is
   open. Underneath it the extremal problem \(\min\max_j\omega_j\) subject to
   \(\pi(\text{set}^c\mid y)\le\delta_\beta\) remains unposed.

   **The level as a function of dimension.** (Not to be confused with the
   bootstrap cliff \(r^\star\) of open item 8, which is unrelated and
   uncharacterised.) For escape budget \(\epsilon\),

   \[
   r^\star_{\mathrm{P2}}(n,\epsilon)=q_{1-\epsilon/(1+\epsilon)}\bigl(\mathrm{Gamma}(n,1)\bigr),
   \qquad
   r^\star_{\mathrm{Gauss}}(n,\epsilon)=q_{1-\epsilon}\bigl(\mathrm{Gamma}(n/2,1)\bigr)
   =\tfrac12\chi^2_{n,1-\epsilon},
   \]

   the second being the 99%-HPD level of \(N(0,I_n)\), whose density threshold is
   \(e^{-r^\star}\) of the mode. Hence the exact identity
   \(r^\star_{\mathrm{P2}}(n)=r^\star_{\mathrm{Gauss}}(2n)\): **(P2) charges the
   level appropriate to twice the dimension.** Cornish–Fisher on
   \(\mathrm{Gamma}(k,1)\) gives \(q_{1-\epsilon}\approx k+z\sqrt k+(z^2-1)/3\),
   so \(r^\star_{\mathrm{P2}}\approx n+z\sqrt n\) and
   \(r^\star_{\mathrm{Gauss}}\approx n/2+z\sqrt{n/2}\) — accurate to three digits
   at \(n=30\), \(\epsilon=0.01\) (\(44.21\) vs \(44.216\); \(25.48\) vs
   \(25.446\)).

   Equivalently, \(r^\star_{\mathrm{Gauss}}(n)\) is the level at which
   \(\Pr[\pi(\beta\mid 0,I_n)/\pi(0\mid 0,I_n)\le e^{-r}]=0.01\), since that
   ratio is \(e^{-\|\beta\|^2/2}\) and \(\|\beta\|^2/2\sim\mathrm{Gamma}(n/2,1)\):

   | \(n\) | 1 | 5 | 10 | 30 | 60 | 100 | 500 |
   |---|---|---|---|---|---|---|---|
   | \(r(n)\) | 3.32 | 7.54 | 11.61 | 25.45 | 44.19 | 67.90 | 288.2 |
   | \(r(n)/n\) | 3.32 | 1.51 | 1.16 | 0.85 | 0.74 | 0.68 | 0.58 |
   | \(e^{-r(n)}\) | 3.6e-2 | 5.3e-4 | 9.1e-6 | 1.1e-11 | 6.4e-20 | 3.2e-30 | 6.5e-126 |

   So \(r(n)/n\downarrow\tfrac12\), the excess over \(n/2\) grows only like
   \(\sqrt n\), and the density ratio defining a 99% region collapses
   geometrically — which is why a certified set sitting at \(10^{-11}\) of the
   modal density at \(n=30\) is not an artefact of the bound.

   The decisive decomposition is \(r(n)=\mathrm{med}(\Xi)+[r(n)-\mathrm{med}(\Xi)]\).
   The first term is \(\approx n/2\) and is where the mass already sits: a
   *typical* draw at \(n=30\) has density ratio \(4.3\times10^{-7}\), with no
   tail behaviour involved. The second term is only \(O(\sqrt n)\) — the set has
   to reach a further factor \(e^{r-\mathrm{med}}\approx 4.8\times10^{4}\) at
   \(n=30\). Hence the tiny absolute level is inherited from the location of the
   mass, not from the size of the escape budget, and shrinking \(\epsilon\)
   cannot recover it.

   **Report the level in radius units.** Because \(\Xi\) is a squared quantity,
   \(r\) is the wrong scale to quote. Reparameterise the same nested family by

   \[
   \widetilde B_r=\Bigl\{\beta:\ \frac{\widetilde\pi(\beta\mid y)}{\widetilde\pi(\beta^\dagger\mid y)}\ \ge\ e^{-s^2/2}\Bigr\},
   \qquad s=\sqrt{2r},
   \]

   which is a monotone change of units, not a new bound. For \(N(0,I_n)\), \(s\)
   is \(\|\beta\|\) in sd units, and \(s(1)=2.576\) recovers the univariate 99%
   \(z\)-value — the reason for \(e^{-s^2/2}\) rather than \(e^{-s^2}\), which
   would report \(2.576/\sqrt2\). Then

   | \(n\) | 1 | 5 | 10 | 30 | 100 | 1000 |
   |---|---|---|---|---|---|---|
   | \(s_{\mathrm{Gauss}}(n)\) | 2.58 | 3.88 | 4.82 | 7.13 | 11.65 | 33.27 |
   | \(\sqrt n+z/\sqrt2\) | 2.65 | 3.88 | 4.81 | 7.12 | 11.65 | 33.27 |
   | \(s_{\mathrm{P2}}/s_{\mathrm{Gauss}}\) | 1.18 | 1.24 | 1.27 | 1.32 | 1.36 | 1.39 |

   so \(s(n)\approx\sqrt n+z_{1-\epsilon}/\sqrt2\) to better than 0.2% for
   \(n\ge5\). Three things become legible at once. The radius grows like
   \(\sqrt n\), not like \(n\), so a 99% region at \(n=30\) is a ball of 7.13 sd
   — unremarkable, and the apparent collapse of \(e^{-r}\) was only \(e^{-s^2/2}\)
   at \(s\approx\sqrt n\). Proposition 2's factor of two in \(r\) is a factor
   \(\sqrt2=1.414\) in \(s\) (1.32 at \(n=30\)). And the certificate is *linear*
   in these units,

   \[
   \eta^\star\approx|\eta^\dagger|+s/\sqrt\kappa,
   \qquad \underline w\sim e^{-s/\sqrt\kappa},
   \qquad 1-\omega\sim e^{-s/\sqrt\kappa},
   \]

   so \(s\) is directly comparable to the per-direction budget
   \(z_{1-\epsilon/2m}\approx\sqrt{2\log(2m/\epsilon)}\), and the ratio
   \(s(n)/z_{1-\epsilon/2m}\) is exactly the exponent of the loss. **Adopt \(s\)
   as the reporting convention throughout**; \(r\) remains the natural variable
   for the Gamma tail in Proposition 2 only.

   Propagating: \(\eta^\star\approx|\eta^\dagger|+\sqrt{2r^\star/\kappa}\) and
   \(1-\omega\approx\lambda_{\min}(P_b^{-1/2}\Gamma P_b^{-1/2})\propto\underline w
   \sim e^{-\eta^\star}\), so

   \[
   1-\omega_{\max}\ \sim\ \exp\bigl(-c\sqrt n\bigr),\qquad n=Jp_{\mathrm{re}},
   \]

   and (P2)'s factor of two costs \(\sqrt2\) in that exponent — it raises the
   floor to the power \(1/\sqrt2\), which is why correcting it moves
   \(\omega_{\max}\) so little.

   **The prize from a per-direction budget, quantified**
   (`_chk_hpd_density_level.R`). The HPD region is a ball of radius
   \(\sqrt{2r^\star_{\mathrm{Gauss}}}\approx\sqrt n+z/\sqrt2\), but a single unit
   functional needs only \(z_{1-\epsilon/2m}\approx\sqrt{2\log(2m/\epsilon)}\)
   after a union bound over \(m\) of them: \(\sqrt n\) against \(\sqrt{\log m}\).
   Since \(w\sim e^{-|\eta|}\) the difference exponentiates into the floor:

   | \(n\) | 5 | 10 | 30 | 60 | 100 |
   |---|---|---|---|---|---|
   | ball radius | 3.88 | 4.82 | 7.13 | 9.40 | 11.65 |
   | per-direction, \(m=100\) | 3.89 | 3.89 | 3.89 | 3.89 | 3.89 |
   | gain \(e^{\Delta}\) | 0.99 | 2.5 | 25.6 | 247 | 2352 |

   with a crossover near \(n\approx5\).

   **How much slack, and where it comes from** (`_ex_logit_floor_J10_p3.R` §8).
   (P2) assumes only convexity, whose extremal case is the **cone**
   \(\Xi(\beta)=\|\beta\|\): linear growth along rays, so \(s\propto v\) and the
   polar Jacobian \(s^{n-1}ds\) becomes \(v^{n-1}dv\), giving \(\mathrm{Gamma}(n)\).
   A near-Gaussian \(\Xi\) grows quadratically, \(s\propto\sqrt v\), and the
   reference is \(\mathrm{Gamma}(n/2)\) — a factor of two in the shape parameter.
   The measured cost, at escape \(0.01\):

   | \(n\) | 6 | 12 | 20 | 30 | 60 |
   |---|---|---|---|---|---|
   | \(r\), convexity only | 13.1 | 21.5 | 31.9 | 44.2 | 79.5 |
   | \(r\), Gaussian reference | 8.4 | 13.1 | 18.8 | 25.5 | 44.2 |
   | ratio | 1.56 | 1.64 | 1.70 | 1.74 | 1.80 |

   climbing toward \(2\) as predicted. At \(n=30\), (P2) charges \(r=44.2\) for an
   escape claim of \(10^{-2}\) where a Gaussian posterior at that \(r\) escapes with
   probability \(1.1\times10^{-7}\).

   Two caveats keep this from being a fix. **It is not what blocks the hard
   cases**: since \(\eta^\star\sim\sqrt r\) and \(w\sim e^{-|\eta|}\), the floor is
   \(e^{-C\sqrt r}\), and in the weak-data \(n=30\) example moving to the tighter
   level changes \(\omega_{\max}\) only from \(0.99996\) to \(0.99966\) — both
   vacuous. **And convexity alone cannot be improved**, since the cone attains
   \(\mathrm{Gamma}(n)\); an extra hypothesis is needed. Strong convexity is
   available here (\(\nabla^2\Xi\succeq\Lambda_\beta\succ0\) globally), but the
   direct route fails: it yields \(s_r\le\sqrt{2r/\mu}\) whereas the ray
   substitution needs \(s_r^2\ge2r/\mu\), the wrong direction. So the tightening is
   plausible and the hypothesis is in hand, but unproved.

   **The level, read as a density ratio, shows the set is the wrong object.**
   Since \(\delta/\widetilde\pi(\beta^\dagger\mid y)=e^{-r}\), the \(n=30\)
   certificate admits everything down to \(6\times10^{-18}\)% of the modal
   density. That is not profligacy in (P2): it is concentration of measure and it
   binds first. For a Gaussian posterior \(\Xi\sim\mathrm{Gamma}(n/2,1)\), so the
   *median* draw already sits \(4\times10^{-5}\)% of the mode at \(n=30\), and 99%
   coverage needs \(8.9\times10^{-10}\)% — irreducibly, however the mass bound is
   proved. Any level set anchored at the mode pays that gap in full
   (`_chk_level_as_density_ratio.R`).

   **The mismatch to fix.** The certificate never needs a set in
   \(\mathbb R^{n}\); by property (i) of §3A.4.5 its entire input is the
   \(\sum_jn_j\) **scalar** ranges \([\eta^-_{ij},\eta^+_{ij}]\). Proving the mass
   bound on a joint set costs \(r\gtrsim n\) where controlling that many scalars
   needs only a union bound, \(r\sim\log(\sum_jn_j)\) — linear in dimension versus
   logarithmic in data size. The joint set is a convenience of (P2)'s radial
   argument, not a requirement of the statement being proved.

   By **Prékopa**, every linear image of a log-concave measure is log-concave, so
   each \(\eta_{ij}\) has a log-concave one-dimensional marginal and the \(n=1\)
   radial bound \(e^{-s}/(1-e^{-s})\) applies with no dimension penalty. For the
   scale, \(\Pr(|\eta_{ij}-\mathbb E\eta_{ij}|\ge t\sigma_{ij})\le e^{1-t}\) with
   \(\sigma^2_{ij}=z_{ij}^\top\mathrm{Cov}(\beta_j\mid y)z_{ij}\), and
   **Brascamp–Lieb** bounds the covariance: \(\nabla^2\Xi\succeq\Lambda\succ0\)
   implies \(\mathrm{Cov}\preceq\Lambda^{-1}\).

   This is *not* yet a fix, for two reasons. Brascamp–Lieb needs a curvature floor,
   and the only globally valid one is \(\Lambda_\beta\) — precisely what made (F)
   vacuous — so the argument is circular and returns to the §3A.4.2 bootstrap. And
   the log-concave constant is crude: at two-sided \(10^{-4}\), \(e^{1-t}\) needs
   \(t\approx10.2\) against \(3.9\) for a Gaussian, a factor \(2.6\) in \(\sigma\)
   units. On the \(J=10\) example, \(\sigma\approx0.97\) from the Schur complement
   gives \(\pm10\) against a half-extent of \(5.3\) from the exact joint solve — so
   at these constants the replacement *loses*. The reformulation is right; the
   constants are the open problem, and they have moved from log-concavity (settled
   by Prékopa) to the tail and covariance inequalities.

   **A better route: per-group marginals, which avoid both inequalities**
   (`_chk_group_marginal_bound.R`). The step above went all the way down to
   \(d=1\) and therefore had to pay for a *scale* — hence Brascamp–Lieb, hence
   the circularity. Stopping at \(d=p_{\mathrm{re}}\) avoids that entirely,
   because (P2) applied to the group marginal is already scale-free. By Prékopa
   the marginal \(\widetilde\pi_j(\beta_j)=\int\widetilde\pi(\beta\mid y)\,
   d\beta_{-j}\) is log-concave on \(\mathbb R^{p_{\mathrm{re}}}\), so (P2)
   applies in dimension \(p_{\mathrm{re}}\) with budget \(\epsilon/J\), and a
   union bound over groups costs only \(\log J\).

   The marginal level set is not computable, but it can be sandwiched by one
   that is. Writing \(\Xi_j^{\mathrm{prof}}(\beta_j)=\min_{\beta_{-j}}\Xi(\beta)\),

   \[
   \frac{\widetilde\pi_j(\beta_j)}{\widetilde\pi(\beta^\dagger\mid y)}
   =e^{-\Xi_j^{\mathrm{prof}}(\beta_j)}\,V(\beta_j),
   \quad
   V(\beta_j)=\!\int\! e^{-[\Xi(\beta_j,\beta_{-j})-\Xi_j^{\mathrm{prof}}(\beta_j)]}d\beta_{-j}.
   \]

   The exponent is convex in \(\beta_{-j}\) with minimum \(0\) and Hessian
   between \(\Lambda_\beta^{-j,-j}\) and \(\Lambda_\beta^{-j,-j}+\bar G_{-j}\)
   (\(\bar G\) the maximal-weight data precision, \(n_{ij}/4\) for logit), so
   \(V\) lies between two Gaussian volumes and \(|\log V(a)-\log V(b)|\le
   \kappa_j:=\tfrac12\log\det(I+(\Lambda_\beta^{-j,-j})^{-1}\bar G_{-j})\).
   Applying this twice — once for the level, once because \(\beta_j^{\mathrm{marg}}
   \ne\beta_j^\dagger\), using \(\Xi_j^{\mathrm{prof}}(\beta_j^{\mathrm{marg}})
   \le\kappa_j\) from maximality — gives the certified inclusion

   \[
   \{\Xi_j^{\mathrm{marg}}\le r\}\ \subseteq\ \{\Xi_j^{\mathrm{prof}}\le r+2\kappa_j\}.
   \tag{M}
   \]

   Nothing new has to be built: \(\max\{c^\top\beta:\Xi(\beta)\le R\}\) for \(c\)
   supported on group \(j\) *is* the support function of
   \(\{\Xi_j^{\mathrm{prof}}\le R\}\), since the KKT system already places the
   other groups at their conditional mode. The scheme is purely a change of
   **level**, from \(r_{\mathrm{joint}}(n)\) to \(r_{\mathrm{grp}}(p_{\mathrm{re}})+2\kappa_j\).

   On the \(J=10\), \(p_{\mathrm{re}}=3\) example at \(\epsilon=0.01\):

   | scheme | \(d\) | sets | budget | \(r\) | \(s\) |
   |---|---|---|---|---|---|
   | joint | 30 | 1 | \(10^{-2}\) | 44.22 | 9.40 |
   | per group | 3 | 10 | \(10^{-3}\) | 11.23 | 4.74 |
   | per functional | 1 | 100 | \(10^{-4}\) | 9.21 | 4.29 |

   so almost all of the available gain is captured at \(d=p_{\mathrm{re}}\);
   going down to \(d=1\) adds only \(4.74\to4.29\) in \(s\) and is what forces
   the Brascamp–Lieb detour.

   The two levers are **not symmetric**, which is what makes the scheme work.
   For fixed \(d\) a Gamma quantile grows like \(\log(1/\epsilon)\), so dividing
   the budget by \(J\) is *additive* — at \(d=3\), \(r:8.42\to11.23\), an increase
   of \(2.81\) against \(\log J=2.30\) — while dividing the dimension by \(J\) is
   *multiplicative*, \(r:44.22\to8.42\). The union bound is cheap and the
   dimension is dear. Consequently the total budget can be tightened almost for
   free: \(\epsilon=10^{-2}\to10^{-3}\) costs 36% of \(1-\omega_{\max}\)
   (\(2.07\to1.32\times10^{-3}\)) for a tenfold safer certificate, and
   \(\epsilon=10^{-4}\) costs 57%. **Everything now rests on \(\kappa_j\)**, and the
   decisive fact is that \(\kappa_j=0\) *exactly* when \(\Xi\) is quadratic,
   because profile and marginal are then the same Schur complement. \(\kappa_j\)
   prices non-Gaussianity alone.

   | | crude sandwich | Laplace truth | break-even |
   |---|---|---|---|
   | \(\kappa_j\) (max over \(j\)) | 18.37 | 0.89 | 16.49 |

   The uniform Hessian sandwich charges \(\tfrac{m}{2}\log(1+\mathrm{SNR})\) with
   \(m=n-p_{\mathrm{re}}=27\) — it prices the *absolute* size of the other
   groups' information — and so misses break-even by 11%. The quantity it is
   bounding, the variation of \(\tfrac12\log\det H_{-j}\) along the extremal
   path, is \(0.89\): an overcharge of a factor \(45\), with \(19\times\)
   headroom to spare. The payoff is flat in \(\kappa\) up to about \(2\):

   | \(\kappa\) | 0 | 1 | 2 | 4 | 8 | 16 |
   |---|---|---|---|---|---|---|
   | \(1-\omega_{\max}\) | 2.97e-3 | 2.07e-3 | 1.48e-3 | 7.99e-4 | 2.73e-4 | 4.62e-5 |
   | gain over joint | 71× | 49× | 35× | 19× | 6.5× | 1.1× |

   against \(1-\omega_{\max}=4.19\times10^{-5}\) for the joint scheme. With
   informative data (25 trials/obs) the difference is qualitative rather than
   cosmetic: \(\omega_{\max}\) falls from \(0.9853\) — effectively no certificate
   — to \(0.8716\), a usable geometric contraction rate, with \(\eta^\star_{\max}\)
   dropping \(13.97\to8.78\) and \(\underline w\) rising \(2.1\times10^{-5}\to
   3.8\times10^{-3}\).

   **Conditional vs profile vs marginal, and which way the error points.**
   Three curvatures for \(\beta_j\), in increasing variance: conditional
   \(H_{jj}\) (others frozen at \(\beta_{-j}^\dagger\)); profile
   \(S_j=H_{jj}-H_{j,-j}H_{-j,-j}^{-1}H_{-j,j}\preceq H_{jj}\) (others relax to
   their conditional mode); and marginal (others integrated out), which equals
   the profile exactly for quadratic \(\Xi\). Measured as the sd of
   \(\eta_{1j}=z^\top\beta_j\):

   | | profile/conditional | marginal/profile |
   |---|---|---|
   | weak data (3 trials/obs) | 1.016 | 1.0008 |
   | informative (25 trials/obs) | 1.0006 | 1.0000 |

   The ordering conditional \(\le\) profile \(\le\) marginal is the safety
   direction: **freezing the other groups is anti-conservative**, shrinking the
   set, raising the weight floor and certifying a contraction the sampler has
   not earned. Note also that the marginal being *wider* is a **cost**, not a
   benefit — at a fixed level a wider set gives a lower floor. The benefit of
   going marginal is entirely the dimension in (P2).

   The gaps are small here for a structural reason worth recording. All group
   coupling sits in the off-diagonal blocks \(-M\) of \(\Lambda_\beta\) with
   \(M=\Lambda_b(\Lambda_\gamma+J\Lambda_b)^{-1}\Lambda_b\), which is monotone
   decreasing in \(\Lambda_\gamma\) and therefore **maximised at
   \(\Lambda_\gamma\to0\), where \(M\to\Lambda_b/J\)**. Coupling can never
   exceed \(O(1/J)\) however weak the population prior, so with \(J=10\) the
   three curvatures must agree to \(O(1/J)\) — confirmed numerically at
   \(\|M\|/\|\Lambda_b\|=0.0999\) against the bound \(0.1\). Small \(J\) is
   where the distinction would bite. The accounting is thus favourable for a
   sharp reason: **the width cost is \(O(1/J)\) while the level gain is a factor
   of four.**

   So the open problem has moved again, and to a much better place. It is no
   longer "beat the dimension" but the single well-posed question: **bound the
   variation of \(\tfrac12\log\det H_{-j}\) along the profile path by something
   closer to \(0.9\) than to \(18\)** — a quantity that vanishes identically in
   the Gaussian case, so any bound respecting that structure (a determinant
   *ratio* rather than a determinant *range*) should succeed. Candidate routes:
   \(|\log\det H(a)-\log\det H(b)|\le\mathrm{tr}(H^{-1}\Delta G)\) by concavity
   of \(\log\det\), combined with \(|w'|\le w\) for the logit, which converts
   \(\Delta G\) into the certified movement of the *other* groups' linear
   predictors rather than their absolute weights.
7. **The example is \(q=1\).** The constant \(c'\) is nearly free there
   (\(1.18\)–\(1.35\)) because it is a single eigenvalue factor. Whether it stays
   benign as \(q\) grows — \(c'=\prod_i(1-\kappa_i^{\mathrm{UB}})^{-1/2}\) over
   \(q\) terms — has not been measured.
8. **The bootstrap cliff is unexplained.** §3A.4.2 finds that the greatest fixed
   point of \(T\) is nontrivial only below a critical level (here \(r^\star\)
   between \(28\) and \(30\)) and collapses to \(\Lambda_\beta\) above it. No
   characterisation of \(r^\star\) in terms of the design exists, so there is no
   way to predict whether a given problem admits a usable floor without running
   the iteration. Relatedly, (P2) is \(10^2\)–\(10^3\) conservative in the
   example, and whether the sharper log-concave form or a curvature-aware
   refinement of (R) inside the boundary would close enough of that gap to move
   \(r(\delta)\) below \(r^\star\) is unknown — that, not the shape question, is
   what currently blocks certifying \(\delta_\beta=0.01\).
9. **What is logit-specific is narrower than it first appeared.** Proposition 2
   needs only convexity of \(\Xi\), which holds for every family checked — logit,
   probit, cloglog, Poisson-log and Gamma-log all have convex per-observation
   negative log-likelihood in \(\eta\). The endpoint rule (EP) also survives
   broadly, from unimodality for logit/probit/cloglog and from monotonicity for
   Poisson-log. What genuinely fails outside logit-like families is the **global
   weight ceiling** \(w\le n/4\) behind \(\Lambda^{\mathrm{UB}}\): for
   Poisson-log \(w=e^\eta\) is unbounded above, so (S), (F) and the §3A.4.2
   bootstrap need a different upper envelope, while the floor machinery itself
   carries over unchanged.

10. **Production implementation.** `group_precision_floor()` in `R/` runs the
   four-step per-group certificate (integrate \(\gamma\) → mode → profile level
   with \(\kappa\) → support-function floors). v1: `gaussian` / `binomial(logit)`
   / `poisson(log)`, fixed `dNormal` \(\tau^2\), `kappa_method = "crude"|"laplace"|"none"`.
   **Theoretical basis:** `inst/LAPLACE_PROFILE_GROUP_MARGINAL_BOUND.md`
   (Proposition L: Gaussian in \(\beta_{-j}\), not in \(\beta_j\); Prop 2 on
   \(\Xi_j^{\mathrm{prof}}\) only; \(\kappa\) prices \(\tfrac12\log\det H_{-j}\)).
   Validation: `data-raw/_chk_group_precision_floor_pkg.R`. Remaining gap: certified
   \(\kappa\) still uses the crude sandwich; tightening to the Laplace-scale
   determinant-ratio bound is open item 6 above.

---

## References

- `vignette("Chapter-C05")` — Theorem 1, Theorem 2, Lemmas 5/13/14/17/24,
  Corollary 15, Remark 6.1 as cited above.
- `CHAPTER_C05_IMPLEMENTATION.md` — §4A (precision-matrix form, closure
  comparison), §4B (tilting route), §7 (identities, ray quadrature for
  \(\bar\Phi\)).
- `RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` — §0.2 pooling weights, §0.4 inner-box
  construction, §1–§6 per-likelihood escape bounds.
- `SAFE_UNSAFE_TV_DECOMPOSITION.md` — the three-term decomposition the
  two-budget version reduces to.
- `HIERARCHICAL_GENERALIZED_LINEAR_MODEL_NOTATION.md` — \(H_j\), \(P_b\),
  \(\Lambda_\gamma\), \(W_j\) conventions.
- `data-raw/_chk_majorization_envelope.R` — the checks of §9.
- `LOGIT_MARGINAL_INTEGRATE_GAMMA.md` — §3–§5 for the unnormalized marginal
  \(\widetilde\pi(\beta\mid y)\), the prior Schur precision \(\Lambda_\beta\) and
  its cross-group coupling; §8, §12 for the absence of a closed-form \(Z(y)\).
  This is the object \(\widetilde B_r\) is a level set of.
- `data-raw/_ex_logit_majorization_floor.R` — the worked logit example of
  §3A.7–§3A.8.
- `data-raw/_chk_logit_levelset_bounds.R` — verification of §3A.4: the sandwich
  (S'), the support-function floor and its bootstrap, the failure of the Gaussian
  majorant, Proposition 2, and the level-set/box comparison.
- `data-raw/_chk_logit_matrix_floor_2d.R` — the map (W), log-concavity of \(w\),
  and the matrix floor (FL) for \(p_{\mathrm{re}}=2\).
- `data-raw/_ex_logit_floor_J10_p3.R` — the full pipeline at \(J=10\),
  \(p_{\mathrm{re}}=3\) in the heavy-shrinkage regime: the KKT path and the budget
  decomposition \(\Xi=L_j+L_{-j}+Q\), free-versus-frozen conditioning, the
  identity \(\omega_j=\) shrinkage weight \(=\) contraction factor, and the (P2)
  slack measurement of open item 6.
