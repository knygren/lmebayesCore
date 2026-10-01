# Hierarchical Linear Models in Standard Form

**Status:** derivation note. Shows how to re-parameterize a two-level hierarchical
linear model so that the population-effect precision becomes the identity matrix and
the sampler's own convergent matrix becomes exactly diagonal — the "standard form"
used throughout the companion theoretical documents (`foundational_lemmas.md`) when
comparing total variation distances via eigenvalues.

## 1. The model, in its natural (design-matrix) form

Following the notation of the `lmebayesCore` package documentation (Chapter C01):
write $j=1,\dots,J$ for the groups. Group $j$ contributes $n_j$ observations $y_j$
with its own design matrix $D_j$ and its own coefficient vector $\beta_j$ (the
**group effects**):
$$
y_j \mid \beta_j \sim N\big(D_j\beta_j,\ \sigma^2 I_{n_j}\big).
$$
The group effects are drawn from a common distribution, itself a regression on
group-level characteristics via a design matrix $\mathcal W_j$:
$$
\beta_j \mid \gamma \sim N\big(\mathcal W_j\gamma,\ \Psi\big),
$$
where $\gamma$ holds the **population effects**, and $\Psi$ is the between-group
variance. Finally $\gamma$ receives its own prior:
$$
\gamma \sim N(\mu_0, V).
$$

**The precision-matrix blocks this induces.** Writing $P_{11}$ for $\gamma$'s own
effective precision (prior plus every group's contribution), $P_{12}$ for the
cross-block linking $\gamma$ to the stacked group effects, and $P_{22}$ for the
(block-diagonal, given conditional independence across groups) total precision of
$\beta:=(\beta_1,\dots,\beta_J)$:
$$
P_{11} = V^{-1} + \sum_{j=1}^J \mathcal W_j^\top\Psi^{-1}\mathcal W_j, \qquad
P_{12} = \big[\mathcal W_1^\top\Psi^{-1},\ \dots,\ \mathcal W_J^\top\Psi^{-1}\big],
$$
$$
P_{22} = \mathrm{blockdiag}\big(\Psi^{-1}+\kappa_{\mathrm{true},1},\ \dots,\ \Psi^{-1}
+\kappa_{\mathrm{true},J}\big), \qquad \kappa_{\mathrm{true},j} :=
\frac{D_j^\top D_j}{\sigma^2}
$$
(the last using the exact Gaussian likelihood case; $\kappa_{\mathrm{true},j}$
generalizes to the local curvature of a non-Gaussian group likelihood elsewhere in
this series).

## 2. Why a "standard form" is wanted

The two-block Gibbs sampler's own convergence rate, and the total-variation-distance
comparisons used to bound it (Nygren, *Total variation distance between multivariate
normal densities...*), are governed entirely by the eigenvalues of
$$
A(P) := P_{11}^{-1/2}P_{12}P_{22}^{-1}P_{21}P_{11}^{-1/2}, \qquad P:=(P_{11},P_{12},
P_{22}),
$$
the **convergent matrix**, written here explicitly as a function of $P$ — this is
deliberate, not decorative, and is explained in §2.1 below. Working directly with
$P_{11},P_{12},P_{22}$ in their natural, design-matrix-derived form obscures this —
the eigenvalue structure is present but not visible. The standard form makes it
explicit: a single re-parameterization of $\gamma$ turns $P_{11}$ into the identity
and $A$ into an exactly diagonal matrix of its own eigenvalues, with no loss of
generality and no approximation.

### 2.1 $C$ is a function of $P$ — which blocks it depends on, and which it doesn't

**Three independent blocks determine $C$, not four.** $P_{21}$ is never a free input
— it is forced to equal $P_{12}^\top$ by symmetry of the joint precision matrix, so
it carries no independent information. Writing
$$
X(P) := P_{11}^{-1/2}P_{12}P_{22}^{-1/2},
$$
$A(P)=X(P)X(P)^\top$, so **$C(P)$ is exactly the left-singular-vector matrix of
$X(P)$**, and $\Lambda(P)$ (the eigenvalues $a_k$) are $X(P)$'s squared singular values —
verified directly. This makes the dependence complete and precise:
$$
\boxed{C = C(P_{11},P_{12},P_{22}).}
$$
Perturbing *any one* of $P_{11},P_{12},P_{22}$ alone, holding the others fixed,
changes $C$ (verified numerically, $d=3,p=4$) — so all three genuinely matter, and
$C$ depends on each one's *full matrix content* (not merely its eigenvectors or its
eigenvalues alone: $P_{11}^{-1/2}$ and $P_{22}^{-1/2}$ each require the whole matrix).

**Consequence: $C$ is not a fixed, universal object for a given model — it depends on
which likelihood is being diagonalized against.** Since $P_{22}=P_{22}(t)=\Psi^{-1}+
\kappa_t$ varies along the floor-to-true homotopy while $P_{11},P_{12}$ stay fixed,
$$
C = C\big(P_{11},P_{12},P_{22}(t)\big) =: C(t)
$$
genuinely moves with $t$ in general. "Standard form" should be read as *standard form
relative to a chosen $P_{22}$* (equivalently, a chosen $t$) — not a single coordinate
system that simultaneously diagonalizes the sampler at every point along the
homotopy. Where a fixed, $t$-independent $C$ is wanted (e.g. to compare across $t$),
it must be chosen deliberately — e.g. $C(1)$, diagonalizing against the true model's
own $P_{22}$ — rather than assumed to exist for free.

### 2.2 $\Lambda$ as a function of $P$, and a monotonicity result

Write $\Lambda(P):=\Lambda(P_{11},P_{12},P_{22})$ for the (sorted) vector of eigenvalues of
$A(P)$, in parallel with $C(P)$ above.

**Key identity.** $A(P)$'s eigenvalues are exactly the *generalized* eigenvalues of
the pair $(S,P_{11})$, where $S:=P_{12}P_{22}^{-1}P_{21}$ (PSD, fixed once $P_{12},
P_{22}$ are fixed) — i.e. the solutions of $Sv=\lambda P_{11}v$. This holds because
$A(P)=P_{11}^{-1/2}SP_{11}^{-1/2}$ is similar to $P_{11}^{-1}S$ (same nonzero
eigenvalues, standard fact for $XY$ vs. $YX$), and $P_{11}^{-1}Sv=\lambda v
\iff Sv=\lambda P_{11}v$. *(Verified directly: $A(P_{11})$'s eigenvalues match the
generalized eigenvalues of $(S,P_{11})$ exactly.)*

**Monotonicity in $P_{11}$, with $P_{12},P_{22}$ fixed.** By the Courant–Fischer
min-max characterization for generalized eigenvalues,
$$
\lambda_k(S,P_{11}) = \min_{\substack{V\subseteq\mathbb R^d\\ \dim V = d-k+1}}\
\max_{v\in V\setminus\{0\}}\ \frac{v^\top Sv}{v^\top P_{11}v}.
$$
If $P_{11}^{(1)}\preceq P_{11}^{(2)}$ (Loewner order), then $v^\top P_{11}^{(1)}v \le
v^\top P_{11}^{(2)}v$ for every $v$, so (since $S\succeq0$ makes the numerator
non-negative) every Rayleigh quotient is no larger under $P_{11}^{(2)}$ than under
$P_{11}^{(1)}$ — hence, by the min-max formula, **every** generalized eigenvalue
$\lambda_k(S,P_{11})$ is monotonically non-increasing as $P_{11}$ increases:
$$
\boxed{P_{11}^{(1)}\preceq P_{11}^{(2)} \implies \Lambda\big(P_{11}^{(2)},P_{12},P_{22}\big)
\le \Lambda\big(P_{11}^{(1)},P_{12},P_{22}\big) \text{ componentwise.}}
$$
*(Verified: 0/200 random trials violated this, for either the full eigenvalue vector
or the maximal eigenvalue alone.)*

**Monotonicity in $P_{22}$, with $P_{11},P_{12}$ fixed.** This is exactly Theorem 2's
own Loewner-chain result, restated in the current notation: $P_{22}$ increasing makes
$P_{22}^{-1}$ decrease, congruence by the fixed $P_{11}^{-1/2}P_{12}$ preserves the
order, so
$$
P_{22}^{(1)}\preceq P_{22}^{(2)} \implies \Lambda\big(P_{11},P_{12},P_{22}^{(2)}\big) \le
\Lambda\big(P_{11},P_{12},P_{22}^{(1)}\big) \text{ componentwise.}
$$

### 2.3 Along the floor-to-true homotopy: $\tilde P_{22}(t)$, $\tilde C(t)$, $\tilde\Lambda(t)$

Since $P_{11},P_{12}$ do not depend on $t$ — only the total $\beta$-precision does,
$\tilde P_{22}(t):=\Psi^{-1}+\kappa_t=(1-t)\underbrace{(\Psi^{-1}+L)}_{=:\Xi_0}+t
\underbrace{(\Psi^{-1}+\kappa_{\mathrm{true}})}_{=:\Xi_1}$, exactly linear in $t$ —
every $t$-dependence in this whole construction factors through this single object.
Define
$$
\tilde C(t) := C\big(P_{11},P_{12},\tilde P_{22}(t)\big), \qquad
\tilde\Lambda(t) := \Lambda\big(P_{11},P_{12},\tilde P_{22}(t)\big),
$$
the eigenvector and eigenvalue functions restricted to this one-parameter family.

**$\tilde\Lambda(t)$ is componentwise strictly decreasing.** $\tilde P_{22}(t)$ is
Loewner-increasing in $t$ ($\tilde P_{22}'(t)=\Xi_1-\Xi_0=\delta\succeq0$, and
strictly so whenever $\delta\succ0$), so §2.2's $P_{22}$-monotonicity result applies
directly along this specific path:
$$
\boxed{t_1 < t_2 \implies \tilde D(t_2) \le \tilde D(t_1) \text{ componentwise (strict
if } \delta\succ0).}
$$
This is exactly Theorem 2 (`foundational_lemmas.md`), now recovered as the special
case of §2.2's general two-block monotonicity result restricted to the specific
linear path $t\mapsto\tilde P_{22}(t)$ that this whole homotopy construction uses.

**$\tilde C(t)$ has no corresponding monotonicity — it can rotate arbitrarily.**
Eigenvectors are not ordered the way eigenvalues are; there is no Loewner-style
statement available for $\tilde C(t)$ itself. Two consequences worth being explicit
about: (i) $\tilde C(0)$ (diagonalizing the floor) and $\tilde C(1)$ (diagonalizing
the true model) are generally *different* orthogonal matrices — the "standard form"
coordinate system is not shared across the homotopy's endpoints, consistent with
§2.1's warning; (ii) whether $\tilde C(t)$ moves continuously in $t$ depends on
$\tilde\Lambda(t)$ having no repeated eigenvalues along the path (crossing eigenvalues can
make the eigenvector assignment discontinuous) — generic for a random $\delta$, but
not guaranteed, and not checked here.

**None of this threatens $\phi(t)$ itself — the total variation distance depends on
$\Lambda$ alone, never on $C$.** Theorem 1's result is stated purely in terms of sorted
eigenvalues; $C$ never enters the distance formula. And eigenvalues of a matrix
varying continuously (here, linearly) in $t$ are themselves always continuous
functions of $t$ — a standard fact that holds regardless of whether the
*eigenvectors* behave continuously at crossings. So $\phi(t)$ is guaranteed
continuous (indeed monotone, by the boxed result above) even in the case where
$\tilde C(t)$ itself jumps. The caveats above matter only if $\tilde C(t)$ is used
for something beyond computing $\phi(t)$ — e.g. tracking the $u$-basis coordinates
themselves, as the $J_\pi,J_Q$ covariance-sum machinery in `foundational_lemmas.md`
does — not for the distance bound this document's whole construction is aimed at.

**Together**: $\Lambda(P_{11},P_{12},P_{22})$ is componentwise non-increasing in *both*
$P_{11}$ and $P_{22}$ separately (Loewner order), with $P_{12}$ held fixed — not
merely the maximal eigenvalue, but every eigenvalue simultaneously, in both cases.

## 3. The re-parameterization

**Diagonalize $A$.** Since $A$ is symmetric (by construction), it admits an ordinary
orthogonal eigendecomposition
$$
A = C\,\Lambda\,C^\top, \qquad \Lambda=\mathrm{diag}(a_1,\dots,a_d), \qquad C^\top C = CC^\top = I,
$$
$d:=\dim(\gamma)$.

**Define the new population parameter.**
$$
u := C^\top P_{11}^{1/2}\gamma, \qquad \gamma = P_{11}^{-1/2}Cu \quad(\text{exact
inverse, since } C \text{ orthogonal}).
$$

**What this does to the precision blocks.** Substituting $\gamma=P_{11}^{-1/2}Cu$
into the joint quadratic form $\gamma^\top P_{11}\gamma - 2\gamma^\top P_{12}\beta +
\beta^\top P_{22}\beta$ and reading off the new coefficients:
$$
\tilde P_{11} = I, \qquad \tilde P_{12} = C^\top P_{11}^{-1/2}P_{12}, \qquad
\tilde P_{21} = \tilde P_{12}^\top, \qquad \tilde P_{22} = P_{22}\ \text{(unchanged)}.
$$

**The convergent matrix in standard form.**
$$
\tilde A := \tilde P_{11}^{-1/2}\tilde P_{12}\tilde P_{22}^{-1}\tilde P_{21}
\tilde P_{11}^{-1/2} = \tilde P_{12}\tilde P_{22}^{-1}\tilde P_{21} = C^\top AC = \Lambda,
$$
using $\tilde P_{11}=I$ and the defining property of $C$. **This is exact, not
approximate**: $\tilde A$ is a similarity transform of $A$ by the very matrix built to
diagonalize it.

*(Verified numerically, $d=3$: $\tilde P_{11}=I$ to machine precision; $\tilde A=
\mathrm{diag}(0.1344,0.2630,0.7588)$, matching $A$'s own eigenvalues exactly.)*

## 4. What transforms, and what doesn't, in design-matrix terms

**The group-level design matrices $D_j$ — untouched.** The re-parameterization acts
on $\gamma$ alone; $\beta$ is never transformed, so $D_j$ (mapping $\beta_j\to y_j$)
stays exactly as calibrated from the data.

**The population-level design matrices $\mathcal W_j$ — transform on the right:**
$$
\boxed{\tilde{\mathcal W}_j := \mathcal W_j\,P_{11}^{-1/2}C.}
$$
This follows directly from substituting $\gamma=P_{11}^{-1/2}Cu$ into
$\beta_j\mid\gamma\sim N(\mathcal W_j\gamma,\Psi)$: the conditional mean becomes
$\mathcal W_jP_{11}^{-1/2}Cu$, so $\tilde{\mathcal W}_j$ is simply the *effective*
design matrix once $u$ replaces $\gamma$ as the parameter.
*(Verified: $P_{12}$ rebuilt from $\{\tilde{\mathcal W}_j\}$ matches
$C^\top P_{11}^{-1/2}P_{12}$ exactly.)*

**The prior — both mean and covariance transform, by the ordinary linear-Gaussian
rule:**
$$
\tilde\mu_0 := C^\top P_{11}^{1/2}\mu_0, \qquad
\tilde V := C^\top P_{11}^{1/2}\,V\,P_{11}^{1/2}C.
$$

## 5. A consequence worth flagging explicitly

**A diagonal, independent-component prior on $\gamma$ generally becomes a fully
correlated prior on $u$.** Verified directly: a diagonal $V=\mathrm{diag}(10,8,12)$
transforms into a dense $\tilde V$ with off-diagonal entries on the order of the
diagonal ones ($\pm10$). This is the price of buying a diagonal $\tilde P_{11}$ and
$\tilde A$ — the standard form is the natural coordinate system for the *sampler's*
convergence properties, not for *prior elicitation*, and the two should not be
conflated. A prior specified independently component-by-component (as
`Prior_Setup_GLMM()` naturally does, matching $\mathcal W_j$'s own block structure
from the model formula) is not independent in $u$-coordinates, and vice versa.

## 6. Why this is the form used throughout the theoretical documents

Once $\tilde P_{11}=I$ and $\tilde A=D$, the target and $l$-step iterate covariances
$\Sigma_t,\Sigma_t^{(l)}$ (in the original $\gamma$-coordinates) become, in
$u$-coordinates, *simultaneously diagonal for every $l$* — this is the exact content
of Nygren's Claim 3/Remark 13, and it is what lets Theorem 1's eigenvalue-only TV
distance result apply directly, coordinate by coordinate, with no separate
diagonalization needed for each $l$. This is the coordinate system implicitly assumed
whenever `foundational_lemmas.md` writes quantities "per eigen-direction $k$" — this
document makes explicit how to get there from an ordinary, design-matrix-specified
hierarchical model, and what that costs in terms of the prior's own interpretability.

**What does not automatically follow.** Diagonalizing $\Sigma_t,\Sigma_t^{(l)}$ this
way says nothing about whether the *non-Gaussian* tilt $\log g(\beta)=
\ell_{\mathrm{true}}(\beta)-\ell_{\mathrm{floor}}(\beta)$ — a $\beta$-space object —
also decomposes cleanly in a basis related to $u$. A natural candidate transform on
$\beta$ (mirroring $A$'s own structure via $P_{22}^{-1/2}P_{21}P_{11}^{-1}P_{12}
P_{22}^{-1/2}$, which shares $A$'s eigenvalues by a standard linear-algebra identity)
does **not** diagonalize $\delta:=\kappa_{\mathrm{true}}-L$ in general — checked
directly, with clearly nonzero off-diagonal terms — because that candidate transform
is built from $P_{22}(t)$, which is $t$-dependent, while $\delta$ is fixed. Whether
*some* fixed transform on $\beta$ aligns $\delta$ with $u$'s own coordinates remains
open; this document only resolves the $\gamma$-side (Gaussian, $\Sigma_t,\Sigma_t^{(l)}$)
half of the picture.

## 7. The posterior marginal precision and variance for $u$

Under the standardized model, $u$'s own posterior marginal precision and variance
have a clean, diagonal closed form, directly in terms of $\Lambda(P_{11},P_{22},P_{12})$
— keeping the functional dependence explicit, as in §2.2–2.3.

**The result.**
$$
\boxed{\tilde\Sigma_u^{-1}(P_{11},P_{22},P_{12}) = I - \Lambda(P_{11},P_{22},P_{12}),}
$$
$$
\boxed{\tilde\Sigma_u(P_{11},P_{22},P_{12}) = \big(I - \Lambda(P_{11},P_{22},P_{12})\big)^{-1}
= \mathrm{diag}\!\left(\frac{1}{1-a_1},\ \dots,\ \frac{1}{1-a_d}\right).}
$$

**Why.** The Schur-complement formula for $\gamma$'s own marginal precision,
$\Sigma_\gamma^{-1}=P_{11}-P_{12}P_{22}^{-1}P_{21}$, applies identically to the
transformed blocks:
$$
\tilde\Sigma_u^{-1} = \tilde P_{11} - \tilde P_{12}\tilde P_{22}^{-1}\tilde P_{21}.
$$
Since $\tilde P_{11}=I$, and
$$
\tilde P_{12}\tilde P_{22}^{-1}\tilde P_{21}
= C^\top P_{11}^{-1/2}P_{12}P_{22}^{-1}P_{21}P_{11}^{-1/2}C
$$
$$
= C^\top A(P_{11},P_{22},P_{12})\,C = \Lambda(P_{11},P_{22},P_{12}),
$$
this collapses to $\tilde\Sigma_u^{-1}=I-\Lambda(P_{11},P_{22},P_{12})$ directly.
*(Verified: matches exactly, both as $I-D$ computed directly and as $\gamma$'s own
Schur-complement precision transformed via $C^\top P_{11}^{-1/2}(\cdot)P_{11}^{-1/2}C$.)*

**Why this is clean, and had to come out this way.** The posterior marginal for $u$
is diagonal — its $d$ components are *posteriorly independent*, each with its own
precision $1-a_k$ — purely as a consequence of the standard-form construction itself,
with no extra assumption needed. This is the same scalar formula
$(\Sigma_t)_{u,kk}=1/(1-a_k)$ used throughout §§2–3 for the target covariance,
now recognized as simply the posterior marginal precision of the standardized
population parameter, stated once in general — as a function of $P_{11},P_{22},
P_{12}$ — rather than re-derived case by case.

**A sanity check on the formula's shape.** Since $a_k\in(0,1)$ (eigenvalues of the
convergent matrix), $1-a_k\in(0,1)$, so $\tilde\Sigma_u(P_{11},P_{22},P_{12})$ is
always finite and positive — and a slow-mixing direction ($a_k$ close to $1$)
corresponds to a
*large* posterior variance in that same direction: a direction the sampler struggles
to move in is also one the data constrain only weakly relative to the prior and
random-effects layer, which is the right intuition to have going in either
direction.

## 8. The posterior mode $u^\star(P)$, by substitution — keeping $C(P)$, $\Lambda(P)$ explicit

$\gamma^\star$, in the model's own design-matrix/prior/data notation, is
$$
\gamma^\star = \Sigma_\gamma\left[V^{-1}\mu_0 + \sum_{j=1}^J\mathcal W_j^\top\Psi^{-1}
\left(\Psi^{-1}+\frac{D_j^\top D_j}{\sigma^2}\right)^{-1}\frac{D_j^\top y_j}{\sigma^2}
\right].
$$
Substituting into $u^\star(P)=C(P)^\top P_{11}^{1/2}\gamma^\star$ — and writing every
transformed object with its $P$-dependence intact, exactly as $C(P)$ and $\Lambda(P)$
themselves carry it — gives:
$$
\boxed{u^\star(P) = \big(I-\Lambda(P)\big)^{-1}\left[\tilde V(P)^{-1}\tilde\mu_0(P) +
\sum_{j=1}^J\tilde{\mathcal W}_j(P)^\top\Psi^{-1}\left(\Psi^{-1}+\frac{D_j^\top D_j}
{\sigma^2}\right)^{-1}\frac{D_j^\top y_j}{\sigma^2}\right],}
$$
where
$$
\tilde{\mathcal W}_j(P) := \mathcal W_j\,P_{11}^{-1/2}C(P), \qquad
\tilde\mu_0(P) := C(P)^\top P_{11}^{1/2}\mu_0, \qquad
\tilde V(P) := C(P)^\top P_{11}^{1/2}\,V\,P_{11}^{1/2}C(P),
$$
and $P:=(P_{11},P_{12},P_{22})$ throughout, as in §§2.1–2.2. *(Verified: matches the
direct transform of $\gamma^\star$ to machine precision.)*

**Nothing here is a fixed object — every tilde quantity, and $u^\star(P)$ itself, is a
function of $P$, inherited entirely from $C(P)$ and $\Lambda(P)$.** In particular,
along the floor-to-true homotopy (§2.3), this becomes $u^\star(t):=\big(I-\tilde
\Lambda(t)\big)^{-1}[\cdots]$ built from $\tilde C(t)$ — tracking both the moving
mean $\gamma^\star(t)$ and the moving coordinate system $C(t)$ together, exactly as
anticipated when $u^\star(t)$ was first introduced as the standardized-space
counterpart of $\gamma^\star(t)$. The formula's *shape* does not change with $t$ —
only $C(t)$ and $\Lambda(t)$, hence every tilde object built from them, do.

**The one part of the formula that never gets a tilde**: $D_j^\top D_j/\sigma^2$ and
$D_j^\top y_j$, consistent with §4 — the group-level design matrices and the data
itself live in $\beta$-space, which this whole re-parameterization never touches.

## 9. The sampler's variance after $k$ iterations, in standard form

Nygren's own $l$-step result (Claim 3), $\Sigma_{11}^{(l)}=P_{11}^{-1/2}\big[\sum_{i=1}^{2l}
A^{i-1}\big]P_{11}^{-1/2}$, specializes immediately once $\tilde P_{11}=I$ and
$\tilde A=\Lambda(P)$: the $P_{11}^{-1/2}$ factors drop out entirely, leaving a pure
geometric matrix series in $\Lambda(P)$, which sums entrywise since $\Lambda(P)$ is
diagonal.

**The result.**
$$
\boxed{\tilde\Sigma_u^{(k)}(P) = \sum_{i=1}^{2k}\Lambda(P)^{i-1}}
$$
$$
\boxed{= \big(I-\Lambda(P)\big)^{-1}\big(I-\Lambda(P)^{2k}\big)
= \tilde\Sigma_u(P)\cdot\big(I-\Lambda(P)^{2k}\big),}
$$
componentwise
$$
\big(\tilde\Sigma_u^{(k)}(P)\big)_{jj} = \sum_{i=1}^{2k}a_j^{i-1} = \frac{1-a_j^{2k}}{1-a_j}.
$$
The first line is Nygren's own formula, $\Sigma_{11}^{(l)}=P_{11}^{-1/2}\big[\sum_{i=1}^{2l}
A^{i-1}\big]P_{11}^{-1/2}$, with the $P_{11}^{-1/2}$ factors simply dropping out once
$\tilde P_{11}=I$ and $A\to\Lambda(P)$ — kept explicit here, ahead of the summed
closed form, since **it is this sum-of-powers form, not the simplified geometric
closed form, that is likely to generalize** beyond the fully Gaussian case (e.g. to
settings where the per-step matrix is not exactly $\Lambda(P)$ raised to a power, but
the telescoping sum structure itself still holds).
*(Verified: matches the corrected recursion run directly in $u$-space, and matches
$\Sigma_\gamma^{(k)}$ — built via the same corrected recursion in the original
$\gamma$-space — transformed via $C(P)^\top P_{11}^{1/2}(\cdot)P_{11}^{1/2}C(P)$,
exactly.)*

**Reading it.** This is Nygren's Remark 11 ($k_i^{(l)}=1/(1-a_i^{2l})$) made fully
explicit in the standardized parameterization: each coordinate's sampler variance
starts at $\tilde\Sigma_u^{(0)}(P)=0$ and grows geometrically toward the target
$\tilde\Sigma_u(P)$ (§7) at rate $a_j^2$ per iteration, with the gap to the target
exactly $\tilde\Sigma_u(P)\Lambda(P)^{2k}$ — the same $a_k^{2l}$ geometric decay that
has recurred throughout this whole investigation, now stated as the sampler's own
$k$-step variance rather than derived as a byproduct of something else.

**How the sum itself changes with $P_{22}$ — the scalar case.** Writing
$a=a(P_{22}):=P_{12}P_{21}/(P_{11}P_{22})$ for the univariate convergent "matrix,"
differentiating $\tilde\Sigma_u^{(k)}(P_{22})=\sum_{i=1}^{2k}a(P_{22})^{i-1}$ term by
term drops the $i=1$ term (constant, zero derivative), so the sum naturally starts at
$i=2$:
$$
\frac{d}{dP_{22}}\tilde\Sigma_u^{(k)}(P_{22}) = \frac{da}{dP_{22}}\sum_{i=2}^{2k}(i-1)\,a^{i-2}.
$$
Using $da/dP_{22}=-a/P_{22}$ (immediate from $a\propto P_{22}^{-1}$):
$$
\boxed{\frac{d}{dP_{22}}\tilde\Sigma_u^{(k)}(P_{22}) = -\frac{1}{P_{22}}\sum_{i=2}^{2k}(i-1)\,a(P_{22})^{i-1}.}
$$
*(Verified against direct numerical differentiation to 10 decimal places.)* Every
term in this sum is non-negative ($a\in(0,1)$, $i-1\ge0$), so the whole derivative is
$\le0$ — the same $P_{22}$-monotonicity established in §2.2 for $\Lambda(P)$ itself,
now shown to propagate cleanly through the entire $k$-step sum term by term, rather
than only through the simplified closed form. This is a second, independent reason
to keep the sum-of-powers form visible: differentiating it directly reproduces the
monotonicity result, without needing to first collapse to the geometric closed form.

**The multivariate case, decomposed along a homotopy in $t$.** Write $P_{22}(t) =
S_1 + PD_{LB} + t(PD-PD_{LB})$, where $S_1$ is the prior part linked to the group
effects ($\Psi^{-1}$), $PD_{LB}$ a data-precision lower bound, $PD$ the true data
precision, and $\delta:=PD-PD_{LB}$ (PSD by construction, $PD\succeq PD_{LB}$) —
exactly $\tilde P_{22}(t)$ from §2.3.

*Step 1: the eigenvalue derivative, via the envelope theorem.* Since only the
current eigenvector $c_j(t)$ is needed even though it is itself changing,
$$
\boxed{a_j'(t) = -c_j(t)^\top P_{11}^{-1/2}P_{12}\,P_{22}(t)^{-1}\,\delta\,
P_{22}(t)^{-1}P_{21}\,P_{11}^{-1/2}\,c_j(t).}
$$
Every $a_j'(t)\le0$, since $\delta\succeq0$ makes this a non-positive quadratic form
— exactly Theorem 2's monotonicity, now at the level of each individual eigenvalue's
derivative rather than the Loewner-chain argument. *(Verified: matches direct
numerical differentiation of $a_j(t)$ exactly, $d=3$.)*

*Step 2: the sum's derivative, in the same form as the scalar case.* Since
$\Lambda(t)$ and $\Lambda'(t)$ are both diagonal (hence commute),
$$
\frac{d}{dt}\tilde\Sigma_u^{(k)}(t) = \sum_{i=2}^{2k}(i-1)\,\Lambda(t)^{i-2}\,\Lambda'(t)
$$
holds entrywise, per eigen-direction $j$:
$$
\boxed{\frac{d}{dt}\big(\tilde\Sigma_u^{(k)}(t)\big)_{jj} = a_j'(t)\sum_{i=2}^{2k}(i-1)\,a_j(t)^{i-2}.}
$$
*(Verified: matches direct numerical differentiation of the sum exactly, all $d=3$
eigen-directions.)* Every factor is now explicit and signed: $a_j'(t)\le0$ (data
precision increasing along the homotopy) times a manifestly non-negative sum, so the
whole derivative is $\le0$ for every $j$, every $k$ — confirming §2.3's monotonicity
propagates through the finite-$k$ sum exactly as it did for the scalar
$P_{22}$-derivative above, now identifying *which* piece of the model (the
data-precision gap $\delta$, via $a_j'(t)$) drives the decrease, separately from the
sum's own combinatorial structure.

## 10. The distance equality — un-normalized sampler, standard-form sampler

Since $u=C(P)^\top P_{11}^{1/2}\gamma$ is an invertible linear map, and total
variation distance is exactly preserved under any invertible linear transformation,
the sampler's distance to its target computed in the original, un-normalized
$\gamma$-space equals the same distance computed entirely in standardized $u$-space:

$$
\boxed{\Big\|N\big(\gamma^\star,\ \Sigma_\gamma^{(k)}\big) - N\big(\gamma^\star,\ \Sigma_\gamma\big)\Big\|_{TV}}
$$
$$
\boxed{=\ \Big\|N\big(u^\star(P),\ \tilde\Sigma_u^{(k)}(P)\big) - N\big(u^\star(P),\ \tilde\Sigma_u(P)\big)\Big\|_{TV}}
$$
$$
\boxed{\le\ \sum_{j=1}^d\left[\frac{a_j(P)^{2k}}{\sqrt{1-a_j(P)^2}\sqrt{1-a_{j-1}(P)^2}}\right]Z_{d-j}}
$$
$$
\boxed{\xrightarrow[k\to\infty]{}\ 0,\quad\text{asymptotic per-step ratio}\ \to\ \big(\max_j a_j(P)\big)^2,}
$$

with $a_0:=0$ and
$$
Z_{d-j} := \frac{e^{-(d-j)/2}\big[(d-j)/2\big]^{(d-j)/2}}{\displaystyle\int_0^\infty e^{-u^2}u^{d-j}\,du}
$$
Nygren's Remark 2 constant — the slope bound on the multivariate error function
$\mathrm{erf}_n$, evaluated at the dimension $d-j$ remaining once $j$ eigen-directions
have already been accounted for in Lemma 2's telescoping construction. It depends
only on dimension, not on $k$ or $t$. This third row is Nygren's own Corollary 1
bound, specialized to this setup: since the sampler is started at the mode,
$\gamma^{(0)}=\gamma^\star$, the mean-mismatch term in Nygren's original bound
vanishes exactly, leaving only the same-mean piece shown here.

**This rate is asymptotic ($k\to\infty$), not the rate at every $k$ — the transient
is generally faster.** Each term decays at its own rate $a_j^{2k}\le(\max_ja_j)^{2k}$,
so at any finite $k$ the sum is a blend of several geometric rates, all at least as
fast as $(\max_ja_j)^2$. Early on, the faster-decaying terms (smaller $a_j$) still
contribute meaningfully, and the *sum's* own decline is correspondingly faster than
$(\max_ja_j)^2$ per step; only once $k$ is large enough for the slowest term to
dominate every other term in the sum does the per-step ratio settle at
$(\max_ja_j)^2$. *(Verified directly: with $a=(0.15,0.35,0.85)$, so
$(\max a_j)^2=0.7225$, the per-step ratio starts at $0.677$ at $k=2$ — genuinely
faster than $0.7225$ — and rises monotonically, reaching $0.7225$ to six decimal
places only by $k\approx9$; using the more closely-spaced $a$'s from the running
numerical example, this transient happens to be short enough that the asymptotic
rate is nearly reached already by $k=2$, which is a property of that particular
example, not of the bound in general.)*

**Why the asymptotic rate is $(\max_ja_j)^2$, specifically.** For large enough $k$
the term with the *largest* $a_j$ eventually decays the slowest and so eventually
dominates every other term in the sum. Since $\max_ja_j(P)$ is exactly the top
eigenvalue of the convergent matrix $A(P)$ — the slowest-mixing direction identified
back in Theorem 2 — the sampler's *ultimate* convergence rate is governed by that
single number, though how quickly that limit is actually reached depends on the full
spread of $\{a_j(P)\}$ and their respective coefficients, not on $\max_ja_j(P)$
alone. This is Nygren's geometric ergodicity claim (Definition 7, Corollary 1), made
concrete here in the standard-form parameterization — as an asymptotic statement,
not a uniform per-iteration rate.

Both sides use the *same mean* for sampler and target — $\gamma^\star$ (equivalently
$u^\star(P)$, §8) — because starting the sampler at the mode gives exact
mean-cancellation, $\mathbb E[\gamma^{(k)}\mid\gamma^{(0)}=\gamma^\star]=\gamma^\star$
for every $k$, which transforms to the same fact for $u^\star(P)$.

**Why this is the useful direction to state it in.** The right-hand side is where
all the earlier work pays off: since $\tilde\Sigma_u(P)$ and $\tilde\Sigma_u^{(k)}(P)$
are *both diagonal in the same basis* (§7, §9), Theorem 1's eigenvalue-only result
applies immediately — the distance depends purely on $\Lambda(P)$, through the ratios
$\kappa_j^{(k)}=(\tilde\Sigma_u)_{jj}/(\tilde\Sigma_u^{(k)})_{jj}=1/(1-a_j^{2k})$ —
with no need to ever diagonalize the *original* $\Sigma_\gamma,\Sigma_\gamma^{(k)}$
pair directly, which do not share eigenvectors with each other in $\gamma$-space at
all. The standardization does not merely simplify bookkeeping — it is what makes the
distance computable via Theorem 1 in the first place.

*(Verified: $\Sigma_\gamma^{(k)}$, built via the corrected $\gamma$-space recursion,
transforms to $\tilde\Sigma_u^{(k)}(P)$ under $C(P)^\top P_{11}^{1/2}(\cdot)
P_{11}^{1/2}C(P)$ exactly matching the closed form of §9 — confirming the equality's
two sides are genuinely the same object viewed in two coordinate systems, not merely
similarly-shaped formulas.)*


