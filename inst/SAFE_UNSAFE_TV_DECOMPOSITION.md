# A Three-Term Total Variation Bound via Safe/Unsafe Decomposition

*Retired from `vignette("Chapter-C05")`, which now covers restricted two-block
Gibbs minorization instead. This note preserves that material unchanged.*

`vignette("Chapter-C03")` reproduces Nygren (2020)'s total variation bound for
a two-block Gibbs sampler whose target is *exactly* multivariate normal.
That argument leans on the target's precision matrix being fixed and known,
so that the \(l\)-step kernel is itself Gaussian (Remark 7) and Theorem 2's
eigenvalue machinery applies to it directly. Once a block's local curvature
is allowed to depend on the current state — a per-group dispersion with its
own prior (`Chapter-C07`), or a non-Gaussian likelihood — the \(l\)-step
kernel is no longer Gaussian, and Theorem 2 no longer applies as it stands.

This note develops the general-purpose device that stands in for it: a
partition of the state space into a *safe* region \(A\), where the local
curvature is bounded well enough to run a comparison-density argument, and
its complement \(A^c\), where only a cruder escape-probability estimate is
available. Sections 1–3 restate Sections 1–3 of `ELLIPSOID_TV_BOUND.md`
model-free: everything there concerns a generic transition kernel
\(P^n(x,\cdot)\), its target \(\pi(\cdot)\), and a measurable partition of
the state space. Sections 6–7 below are specific to this note and have no
counterpart in `ELLIPSOID_TV_BOUND.md`. `Chapter-C06`/`Chapter-C07` and that
note's §5–§8 supply the concrete \(A\), curvature floor and likelihood that
turn each general lemma below into a number.

Throughout, `x = (gamma, beta)` denotes the full state of the sampler
(population effects \(\gamma\), group effects \(\beta\)), \(P^n(x,\cdot)\)
is the \(n\)-step transition kernel starting from \(x\), and \(\pi(\cdot)\)
is the target density. The total variation distance \(\|\cdot\|_{TV}\) is
Definition 2 of `Chapter-C03`.

---

## 1. General setup and notation

For a generic block \(\beta_j\) (dimension \(p_j\), design matrix \(D_j\),
group size \(n_j\)), write its contribution to the negative log-target as
\[
-\log \pi(\beta_j \mid \text{rest}) =
-\sum_i \log L(y_i ; \eta_i) + \tfrac{1}{2}\tau_j^{-2}\|\beta_j\|^2 + \text{const},
\qquad \eta_i = D_j[i,\,]\,\beta_j,
\]
with local Hessian
\[
H_j(\beta_j) = \tau_j^{-2} I_{p_j} + D_j' \operatorname{diag}\!\big(c(\eta_i)\big) D_j,
\qquad
c(\eta) := -\frac{d^2}{d\eta^2}\log L(y;\eta),
\]
where \(c(\cdot)\) is the likelihood's curvature function — equal to the
precision \(\Omega_j\) in the Gaussian/t case, to \(p(1-p)\) for logit, and
so on. \(\tau_j^{-2}\) is the prior precision on \(\beta_j\), which is why it
always appears as a floor: however uninformative the data become, the prior
still supplies it.

> **Assumption 1 (full-rank / estimability).** \(D_j\) has full column rank
> \(p_j\), and (for non-Gaussian likelihoods) the group-level MLE
> \(\hat\beta_j\) exists and is finite — i.e. the group's data do not exhibit
> complete or quasi-complete separation and do not have a degenerate
> sufficient statistic. This is the same full-rank / estimability screen the
> package already runs on \(D_j\) before attempting a group-level fit, and it
> is what guarantees the anchor point used throughout (\(\hat\beta_j\) for
> GLM likelihoods, the OLS estimate for the Gaussian/t case) is well defined.

> **Definition 1 (safe region via a precision floor).** Fix a threshold
> \(\kappa_j > 0\). Block \(j\) is *safe* at \(\beta_j\) if
> \[
> \lambda_{\min}\!\big(H_j(\beta_j)\big) \;\ge\; \kappa_j.
> \]
> The **safe region** is
> \(A := \bigcap_j \{\beta_j : \lambda_{\min}(H_j(\beta_j)) \ge \kappa_j\}\)
> (intersected, where relevant, with analogous conditions for
> population/random-effect blocks); \(A^c\) is its complement.

Definition 1 is deliberately abstract. `Chapter-C03`'s setting — fixed,
known precision matrices — sits inside it as the degenerate case
\(A = \) the whole space, \(A^c = \varnothing\): the curvature never varies,
so there is nothing to be unsafe about, and every bound below reduces to
zero extra terms. The bounds in this note are needed only once curvature
can vary with the state, which is exactly what happens once a variance
component or dispersion parameter gets a prior instead of a fixed value.

---

## 2. The safe/unsafe decomposition

The starting point is a triangle-inequality split across any partition of
the state space into \(A\) and \(A^c\), with no assumption yet on what \(A\)
is or on the likelihood generating \(P^n\), \(\pi\).

> **Lemma 1 (total variation splits across a partition).** For any
> measurable partition \(\Omega = A \cup A^c\),
> \[
> \|P^n(x,\cdot) - \pi(\cdot)\|_{TV} \;\le\; d_A(x) \;+\; P^n(x,A^c) \;+\; \pi(A^c),
> \]
> where \(d_A(x) := \sup_{B \subseteq A} \big|P^n(x,B) - \pi(B)\big|\).

*Proof.* For any measurable \(B\), write
\(B = (B \cap A) \cup (B \cap A^c)\), so
\[
P^n(x,B) - \pi(B) =
\big[P^n(x, B\cap A) - \pi(B \cap A)\big] + \big[P^n(x, B \cap A^c) - \pi(B \cap A^c)\big].
\]
Take absolute values and the supremum over \(B\). The supremum of the first
bracket over subsets of \(A\) is \(d_A(x)\) by definition. For the second,
since \(B \cap A^c \subseteq A^c\), \(P^n(x, B\cap A^c) \le P^n(x,A^c)\) and
\(\pi(B \cap A^c) \le \pi(A^c)\), so its contribution is at most
\(P^n(x,A^c) + \pi(A^c)\). \(\square\)

Lemma 1 is entirely likelihood-agnostic: \(A\) and \(A^c\) only need to
partition the state space. It replaces a single hard-to-control quantity —
the raw total variation distance — with three separately tractable
quantities: an "on-\(A\)" gap where the comparison-density machinery of
`Chapter-C03` can be brought to bear, and two escape probabilities that
measure how much of the sampler's (respectively the target's) mass ever
leaves the region where that machinery applies.

### 2.1 A sharper split: converting the escape-mass term into an estimable quantity

The two escape terms \(P^n(x,A^c)\) and \(\pi(A^c)\) play different roles.
\(P^n(x,A^c)\) is, by definition, the probability that the chain — started
at \(x\), run for \(n\) sweeps — lands in the unsafe region; that is exactly
what a running MCMC chain lets you estimate directly, as the empirical
proportion of post-burn-in draws falling in \(A^c\). \(\pi(A^c)\), by
contrast, is a property of the *target* alone, and is typically available
only through a closed-form (and often loose) tail bound.

The following identity — true for any two real numbers, so it needs no
hypotheses at all — lets a bound built from both quantities be re-expressed
so that the directly-estimable one carries all of the weight and the other
survives only as a correction:
\[
2\,P^n(x,A^c) + 2\,\pi(A^c)
\;=\;
4\,P^n(x,A^c) + 2\big[\pi(A^c) - P^n(x,A^c)\big].
\]
This is a no-op unless the bracketed gap \(\pi(A^c) - P^n(x,A^c)\) can be
argued to be small — which it is, under exactly the geometric-ergodicity
hypotheses (a Foster–Lyapunov drift/minorization pair) that the package's
convergence machinery already establishes for the likelihoods it supports;
see `ELLIPSOID_TV_BOUND.md` §2.1 and §5–§8 for the argument in each
case. That refinement is outside the scope of this note, which stops at
the algebraic identity above; the point of recording it here is that it is
what turns the Corollary of the next section into the practical, empirically
estimable inequality this note is building toward.

---

## 3. Bounding \(d_A\): a normalization-corrected version

### 3.1 Why the naive approach is invalid

\(P^n(x,\cdot)\) and \(\pi(\cdot)\) restricted to \(B \subseteq A\) are
*sub-probability* measures, integrating to \(p := P^n(x,A) \le 1\) and
\(q := \pi(A) \le 1\) respectively, generally with \(p \ne q\). A theorem
stated for two properly normalized densities — such as `Chapter-C03`'s
Theorem 2, or any comparable comparison device for a GLM likelihood — cannot
be applied directly to these unnormalized restrictions: doing so silently
compares two measures of different total mass.

The fix is to compare the *conditional* laws on \(A\) instead, and to
account separately for the normalization each one drops.

> **Definitions.** For \(B \subseteq A\),
> \(P^n_A(x,B) := P^n(x,B)/p\) and \(\pi_A(B) := \pi(B)/q\) — the properly
> normalized conditional laws on \(A\).

> **Lemma 5 (three-term decomposition).** For \(B \subseteq A\),
> \[
> P^n(x,B) - \pi(B)
> = \big[P^n_A(x,B) - \pi_A(B)\big]
> + \big[P^n(x,B) - P^n_A(x,B)\big]
> - \big[\pi(B) - \pi_A(B)\big].
> \]

*Proof.* Add and subtract \(P^n_A(x,B)\) and \(\pi_A(B)\); the resulting
identity is pure algebra, independent of the likelihood generating
\(P^n\), \(\pi\). \(\square\)

> **Lemma 6 (bound on \(d_A\)).**
> \[
> d_A(x) \;\le\; \|P^n_A(x,\cdot) - \pi_A(\cdot)\|_{TV} \;+\; P^n(x,A^c) \;+\; \pi(A^c).
> \]

*Proof.* Take \(\sup_{B \subseteq A}\) of Lemma 5 and apply the triangle
inequality across its three terms. For the first, the supremum over
\(B \subseteq A\) is \(\|P^n_A(x,\cdot) - \pi_A(\cdot)\|_{TV}\) by
Definition 2 (`Chapter-C03`). For the second and third, the identity
\(P^n(x,B) = p \cdot P^n_A(x,B)\) (and similarly \(\pi(B) = q \cdot \pi_A(B)\))
holds for any likelihood, so
\(\sup_B |P^n(x,B) - P^n_A(x,B)| = |1-p|\,\sup_B P^n_A(x,B) \le (1-p) = P^n(x,A^c)\),
and likewise \(\sup_B|\pi(B) - \pi_A(B)| \le (1-q) = \pi(A^c)\). \(\square\)

> **Corollary (cost of the fix).** Combining Lemma 1 and Lemma 6,
> \[
> \|P^n(x,\cdot) - \pi(\cdot)\|_{TV}
> \;\le\;
> \|P^n_A(x,\cdot) - \pi_A(\cdot)\|_{TV}
> \;+\; 2\,P^n(x,A^c) \;+\; 2\,\pi(A^c).
> \]

Substituting \(d_A(x)\)'s bound from Lemma 6 into Lemma 1's escape terms is
what doubles the coefficient on \(P^n(x,A^c)\) and \(\pi(A^c)\): the "cost"
of correctly normalizing the on-\(A\) comparison is that each escape term
now has to account for both the outer split (Lemma 1) and the inner one
(Lemma 6). This is the price paid for working with the exactly-normalized
conditional laws \(P^n_A, \pi_A\) — the only laws a Gaussian-comparison
theorem (or any comparable GLM device) can legitimately be applied to.

### 3.2 Bounding the conditional term \(\|P^n_A - \pi_A\|_{TV}\)

\(P^n_A(x,\cdot)\) and \(\pi_A(\cdot)\) are the sampler and target laws
truncated to \(A\) and renormalized — not, in general, literally Gaussian
(or any other tractable closed form) even when the untruncated density is.
Two routes, likelihood-agnostic in statement though not in execution:

- **Route 1 — Pinsker + truncated KL.**
  \(\|P^n_A - \pi_A\|_{TV} \le \sqrt{0.5\,D_{KL}(P^n_A \| \pi_A)}\)
  (Theorem 1 of `Chapter-C03`), computing the Kullback–Leibler divergence
  between the two densities truncated to the common region \(A\).
- **Route 2 — domination by an untruncated comparison density.** Bound the
  true truncated target on \(A\) above and below by an untruncated
  comparison density built from the curvature floor \(\kappa_j\) — a
  Gaussian in every case considered here, since \(H_j(\beta_j) \succeq
  \kappa_j I\) on \(A\) by Definition 1 — then apply a Gaussian-comparison
  theorem (`Chapter-C03`'s Theorem 2) to that comparison pair, absorbing the
  domination gap as additional slack.

Neither route is carried to a closed form here; that is deliberately left
to the model- and likelihood-specific notes, since what the comparison
density looks like depends on which curvature floor \(\kappa_j\) applies.

---

## 4. The assembled bound

Putting the Corollary above together with the algebraic identity of §2.1
gives the headline result of this note.

> **Theorem 1 (safe/unsafe total variation bound, practical form).**
> \[
> \|P^n(x,\cdot) - \pi(\cdot)\|_{TV}
> \;\le\;
> \underbrace{\|P^n_A(x,\cdot) - \pi_A(\cdot)\|_{TV}}_{\text{on-}A\text{ gap}}
> \;+\;
> \underbrace{4\,P^n(x,A^c)}_{\text{estimable}}
> \;+\;
> \underbrace{2\big[\pi(A^c) - P^n(x,A^c)\big]}_{\text{correction}}.
> \]

*Proof.* By the Corollary,
\(\|P^n(x,\cdot)-\pi(\cdot)\|_{TV} \le \|P^n_A(x,\cdot)-\pi_A(\cdot)\|_{TV} + 2P^n(x,A^c) + 2\pi(A^c)\).
Apply the identity
\(2P^n(x,A^c)+2\pi(A^c) = 4P^n(x,A^c) + 2[\pi(A^c)-P^n(x,A^c)]\)
to the last two terms. \(\square\)

Each of the three terms plays a distinct role, and is controlled by a
different piece of machinery:

1. **The on-\(A\) gap**, \(\|P^n_A(x,\cdot)-\pi_A(\cdot)\|_{TV}\), is where
   the comparison-density theory of `Chapter-C03` (or its GLM analogues) is
   applied — on the safe region only, where the curvature floor of
   Definition 1 gives a usable sandwich on \(H_j\).
2. **The escape term**, \(4\,P^n(x,A^c)\), is the probability that the
   *sampler itself* has, after \(n\) sweeps, landed outside the safe region.
   Unlike a closed-form tail bound, this is exactly what a stored chain lets
   you estimate directly: the empirical proportion of post-burn-in draws in
   \(A^c\), with a Monte Carlo standard error or a Clopper–Pearson interval
   if a rigorous band is wanted.
3. **The correction term**, \(2[\pi(A^c) - P^n(x,A^c)]\), is the gap between
   the *stationary* escape probability and the sampler's, after \(n\)
   sweeps. It is a no-op unless this gap can be shown to be small; under a
   geometric-ergodicity (Foster–Lyapunov) argument it decays like \(\rho^n\)
   for some \(\rho \in (0,1)\) (Meyn & Tweedie, 1993, Thm 16.0.1), which is
   exactly the kind of argument the package's likelihood-specific
   convergence certificates are built to supply. Making \((\rho, R)\)
   explicit for a given likelihood is beyond Sections 1–3 of
   `ELLIPSOID_TV_BOUND.md`, and is left to that note's later sections.

The practical upshot is a bound whose *dominant*, empirically checkable term
is \(4\,P^n(x,A^c)\) — four times the fraction of stored draws that ever
left the safe region — rather than a closed-form but possibly loose
analytic bound on \(\pi(A^c)\) alone.

---

## 5. Connection to the package's convergence diagnostics

Theorem 1 is the properly normalized counterpart of a simpler diagnostic
already implemented (as a `known_vcov`-engine, empirical augmentation of
the calibrated `tv_tol`) in `BLOCK_GIBBS_ERGODICITY_ING.md` §17. That
diagnostic applies Lemma 1 directly — without Lemma 6's normalization
correction — and approximates *both* escape terms, \(P^n(x,A^c)\) and
\(\pi(A^c)\), by the same empirical fraction \(\hat p\) of main-stage draws
outside the safe region, giving
\[
\mathrm{TV}_{\mathrm{revised}} = \texttt{tv\_tol} + 2\hat p.
\]
Theorem 1 shows what changes if the on-\(A\) term is instead the properly
normalized \(\|P^n_A - \pi_A\|_{TV}\) rather than the raw \(d_A(x)\) that
Lemma 1 leaves unbounded: the escape terms double, to \(4\hat p\) under the
same same-\(\hat p\) heuristic, because the normalization fix (Lemma 6)
introduces its own copy of each escape term on top of the one Lemma 1
already contributes. The §17 diagnostic and Theorem 1 are therefore two
different points on the same trade-off — a smaller multiplier on \(\hat p\)
paired with an unnormalized, harder-to-interpret on-\(A\) term, versus a
larger multiplier paired with a properly normalized one — rather than
competing derivations of the same quantity.

---

## 6. Open items

Sections 1–4 assemble Theorem 1 but deliberately leave two of its three
terms unfinished — bounded in the abstract, but not yet reduced to a
concrete, checkable condition. Both are open questions this note does
not resolve.

1. **When does the on-\(A\) gap inherit `Chapter-C03`'s bound?**
   `Chapter-C03`'s Theorem 2 bounds the total variation distance between two
   multivariate normal densities on the *whole* space. Route 2 of §3.2
   proposes reaching that bound on the safe region by *dominating* the true
   conditional laws \(P^n_A(x,\cdot)\), \(\pi_A(\cdot)\) with an untruncated
   comparison Gaussian built from the curvature floor \(\kappa_j\), but does
   not say when the domination gap is small enough for the resulting bound
   to be useful rather than vacuous. Concretely: what conditions on \(A\)
   (its shape, how \(\kappa_j\) compares to the floor/ceiling of \(H_j\) on
   \(A\), how far \(A\) sits from the anchor \(\hat\beta_j\)) guarantee that
   \[
   \|P^n_A(x,\cdot) - \pi_A(\cdot)\|_{TV}
   \]
   is bounded by `Chapter-C03`'s Theorem 2 applied directly to the
   *untruncated* comparison pair — with no extra slack — versus requiring an
   explicit additive correction for the truncation itself? A partial answer
   would identify a class of \(A\) (e.g. a sublevel set of the quadratic
   form \(q_j\), as in the ellipsoid case of `ELLIPSOID_TV_BOUND.md`
   §5.1) for which the truncated and untruncated comparisons coincide, or
   differ by a correction that itself vanishes as \(A\) grows to fill the
   space. §7 below gives one partial route: it trades the need for a
   comparison theorem on the *renormalized* conditional laws for a
   comparison theorem on the *unnormalized* restrictions \(P^n(x,\cdot)|_A\),
   \(\pi(\cdot)|_A\), at the cost of a term that is controlled by the same
   escape-probability gap as item 2 below.

2. **When does the correction term shrink quickly, and for which
   likelihoods/priors?** Theorem 1's last term,
   \(2[\pi(A^c) - P^n(x,A^c)]\), is a no-op unless it can be shown to decay
   in \(n\). §2.1 notes that a Foster–Lyapunov drift/minorization pair gives
   geometric decay in the abstract (Meyn & Tweedie, 1993, Thm 16.0.1;
   `ELLIPSOID_TV_BOUND.md` Claim 5), but the resulting rate \(\rho\)
   and constant \(R\) are explicit only for the multivariate-t/ellipsoid
   case (§5.3 of that note, via the exact corner-covariance rate) and remain
   abstract — asserted to exist, not computed — for probit and Poisson
   (§7.4, §8.4), while logit gets a sharper *uniform*-ergodicity shortcut
   that avoids the \(V\)-norm argument entirely (§6.3). The open question is
   which combinations of log-likelihood and prior admit an explicit,
   checkable \((\rho, R)\) or drift pair \((\lambda, b, \epsilon)\) — beyond
   the multivariate-t and logit cases already worked out — and what the
   resulting rate looks like as a function of the curvature floor
   \(\kappa_j\) and the prior precision \(\tau_j^{-2}\), so that "shrinks
   quickly" becomes a number rather than a qualitative claim.

---

## 7. A direct bound on the on-\(A\) gap via the escape-probability gap

This section records a candidate route toward Open Item 1, obtained by
add-and-subtract on the renormalized densities themselves rather than by a
comparison-density argument. As in the discussion above, write \(f, g\) for
the densities of \(P^n(x,\cdot)\), \(\pi(\cdot)\) restricted to \(A\), and
suppose (as in the worked example above) that
\(p := P^n(x,A) > \pi(A) =: q\).

For \(y \in A\), add and subtract \(g(y)/p\) — that is, \(\pi(y)\) rescaled
by the *sampler's* mass on \(A\) rather than by its own:
\[
f_A(y) - g_A(y) \;=\; \frac{f(y)}{p} - \frac{g(y)}{q}
\;=\; \underbrace{\frac{f(y)-g(y)}{p}}_{\text{(i)}}
\;-\; \underbrace{g(y)\,\frac{p-q}{pq}}_{\text{(ii)}}.
\]
Term (i) is the *unnormalized* restricted difference, divided by the
constant \(p\) — no renormalization of either density is needed to make
sense of it, unlike \(f_A\) or \(g_A\) individually. Term (ii) is
\(\pi\)'s density on \(A\), scaled by a single scalar built from the two
escape masses.

> **Lemma 7 (on-\(A\) gap via the unnormalized restriction).** With \(p, q\)
> as above,
> \[
> \big\|P^n_A(x,\cdot) - \pi_A(\cdot)\big\|_{TV}
> \;\le\;
> \frac{d_A(x) + |p-q|}{p}.
> \]

*Proof.* Integrating the display above over \(B \subseteq A\),
\[
P^n_A(x,B) - \pi_A(B) = \frac{P^n(x,B)-\pi(B)}{p} - \frac{p-q}{pq}\,\pi(B).
\]
Taking absolute values, applying the triangle inequality, and then the
supremum over \(B \subseteq A\),
\[
\big\|P^n_A(x,\cdot)-\pi_A(\cdot)\big\|_{TV}
\;\le\;
\frac{1}{p}\sup_{B\subseteq A}\big|P^n(x,B)-\pi(B)\big|
\;+\;
\frac{|p-q|}{pq}\sup_{B\subseteq A}\pi(B).
\]
The first supremum is \(d_A(x)\) by definition (Lemma 1). The second is
\(\pi(A) = q\), attained at \(B=A\) itself, since \(\pi\) is a positive
measure and therefore monotone in \(B\). Substituting and simplifying
\(\frac{|p-q|}{pq}\cdot q = \frac{|p-q|}{p}\) gives the claim. \(\square\)

Two things are worth noting about Lemma 7. First, since
\(p - q = \pi(A^c) - P^n(x,A^c)\) exactly, \(|p-q|\) is the same
escape-probability gap that appears in Theorem 1's correction term — so
Lemma 7 makes the connection between Open Items 1 and 2 concrete rather than
merely thematic: whatever argument controls how quickly Theorem 1's last
term shrinks *also* controls the size of the extra slack Lemma 7 pays for
routing around a renormalized comparison. Second, and more usefully, Lemma 7
replaces the problem of bounding \(\|P^n_A - \pi_A\|_{TV}\) — a distance
between two *properly normalized* conditional laws, for which Route 1/2 of
§3.2 need machinery that does not yet exist for truncated, renormalized
densities — with the problem of bounding \(d_A(x)\), a distance between the
*unnormalized* restrictions of \(P^n(x,\cdot)\) and \(\pi(\cdot)\) to \(A\).
The latter does not have a mismatched-total-mass problem to begin with
(nothing has been divided by \(p\) or \(q\) inside the supremum that defines
it), so a domination argument in the style of Route 2 — sandwiching \(f\)
and \(g\) on \(A\) between untruncated comparison Gaussians built from the
curvature floor \(\kappa_j\) — can be applied to \(f\) and \(g\) directly,
without first having to justify normalizing two sub-probability measures of
different total mass. Carrying that domination step out explicitly, for a
concrete class of safe regions \(A\), is what would turn Lemma 7 into a full
answer to Open Item 1; it is not carried out here.

### 7.1 A computable heuristic: Corollary 1 in place of \(d_A(x)\)

**The underlying idea, stated globally first.** Suppose the sampler's
precision were such that the two-block Gibbs contraction matrix's
eigenvalues were bounded *above, everywhere*, by values strictly less than
1 — i.e. a single eigenvalue bound valid on the whole state space, not just
on \(A\). Then it is natural to expect that the total variation distance
for the two-block Gibbs sampler on a multivariate normal density
*consistent with that eigenvalue bound* — the comparison chain
`Chapter-C03`'s Corollary 1 is a proven statement about — should dominate
(be an upper bound for) the total variation distance of the *actual*
sampler, whose true, state-dependent contraction is no worse than that
bound anywhere. This is the domination principle behind every "worst-case
rate" already used elsewhere in the package (`two_block_l_for_tv()`'s
sweep budget, the ING marginal safeguard); it is a well-motivated
expectation, not a theorem proved here — Corollary 1 itself only ever
establishes the bound for the literally-fixed-precision Gaussian chain, not
for a chain dominated by it in this eigenvalue sense.

**Restricting the eigenvalue bound to a safe region.** If that eigenvalue
bound can only be established on the safe region \(A\) (Definition 1's
floor \(\kappa_j\), or a pilot ceiling as in
`ELLIPSOID_TV_BOUND.md` §6.1) rather than everywhere, the global
picture above does not apply as stated — but Theorem 1 has already done
the work of separating exactly the two pieces this requires. The part of
the state space *outside* \(A\), where no such bound is available, is
handled by Theorem 1's escape terms, \(4P^n(x,A^c)+2[\pi(A^c)-P^n(x,A^c)]\),
not by this section. What remains is the on-\(A\) term itself, and that is
where the eigenvalue bound — even though only established on \(A\) — is
put to use: Lemma 7,
\[
\big\|P^n_A(x,\cdot) - \pi_A(\cdot)\big\|_{TV}
\;\le\;
\frac{d_A(x) + |p-q|}{p},
\]
already isolates the on-\(A\) gap as a *scaled* version (by \(1/p\)) of
\(d_A(x)\), plus one additional, exact term \(|p-q|/p\) that is the price
of that scaling (its proof, above, needs nothing beyond the triangle
inequality). \(d_A(x)\) — the *unscaled* distance on the safe region,
before Lemma 7's \(1/p\) rescaling — is not directly computable once the
curvature on \(A\) is state-dependent; this is exactly Open Item 1. The
heuristic step is a direct substitution for \(d_A(x)\) itself: in its
place, use the number Corollary 1 already returns for the fixed-precision
comparison chain over the *full* region, evaluated at the safe region's
bounding precision — exactly what `two_block_rate()`/`two_block_tv_bound()`
already compute. That is, the full-region bound is used as an upper bound
for the unscaled on-\(A\) distance \(d_A(x)\), not the other way around,
and not by claiming the eigenvalue bound itself literally holds outside
\(A\). Spelled out, Corollary 1's bound is
\[
\begin{aligned}
\widehat d \;:=\;&
\frac{0.5(\lambda^{*})^{l}\sqrt{\left(x_1^{(0)}-\mu_1\right)^{T}\Sigma_{11}^{-1}\left(x_1^{(0)}-\mu_1\right)}}{\sqrt{2}}
\left[\frac{e^{-(n-i)/2}\left[\sqrt{(n-i)/2}\right]^{n-i}}{\int_0^{\infty}e^{-u^2}u^{\,n-i}du}\right] \\
&\;+\;
\sum_{i=1}^{n}\left[\frac{a_i^{2l}}{\sqrt{1-a_i^{2}}\sqrt{1-a_{i-1}^{2}}}\right]
\left[\frac{e^{-(n-i)/2}\left[\sqrt{(n-i)/2}\right]^{n-i}}{\int_0^{\infty}e^{-u^2}u^{\,n-i}du}\right],
\end{aligned}
\]
with \(\lambda^{*}\) and \(a_i\) computed from \(P_{11}, P_{12}, P_{22}\)
evaluated at the safe region's bounding precision. Substituting
\(\widehat d\) for \(d_A(x)\) in Lemma 7 gives the practical, computable
heuristic:
\[
\big\|P^n_A(x,\cdot) - \pi_A(\cdot)\big\|_{TV}
\;\le\;
\frac{d_A(x) + |p-q|}{p}
\;\;\overset{\text{heuristic}}{\lesssim}\;\;
\frac{\widehat d}{p} + \frac{|p-q|}{p}.
\]
The first inequality is exact — it is Lemma 7. The second,
\(d_A(x)\overset{\text{heuristic}}{\lesssim}\widehat d\), is the one
approximation in the chain: it uses \(\widehat d\) — a bound established
for a fixed-precision Gaussian chain over the *entire* state space — as an
upper bound for \(d_A(x)\), the unscaled distance defined only over the
safe region \(A\) for the true, state-dependent sampler. This is the
global domination principle above, applied with the *only* eigenvalue
bound actually in hand — the one certified on \(A\) — pressed into service
as a stand-in for the whole-space bound the domination principle actually
calls for. It is not a claim that the on-\(A\) eigenvalue bound itself
extends to \(A^c\); \(A^c\) is handled entirely by Theorem 1's escape
terms. Nothing else in the substitution is approximate: the \(1/p\)
scaling and the additional \(|p-q|/p\) term both come from Lemma 7's exact
algebra, not from the heuristic step.

The corresponding caution is the same one already attached to every other
worst-case rate in the package: because \(\widehat d\) is computed at a
single fixed precision rather than the true, varying one, it can
*understate* the on-\(A\) gap whenever the state-dependent curvature spends
a non-trivial fraction of its time near the edge of the sandwich rather
than comfortably inside it — the same caveat `ELLIPSOID_TV_BOUND.md`
§6.1 attaches to treating \(M_j^{\text{pilot}}\) as a fixed threshold
rather than an estimated one.

---

## References

- `vignette("Chapter-C03")` — Nygren (2020)'s total variation theory for the
  exactly-Gaussian case, which is the tool Theorem 1's on-\(A\) term
  reaches for once the safe region has a usable curvature sandwich, and the
  source of Definition 2, Theorem 2 and Corollary 1 as cited above.
- `vignette("Chapter-C04")` — the two-block Gibbs sampler this note's
  \(P^n(x,\cdot)\) refers to, in the fixed-precision case where
  \(A^c = \varnothing\) and every term but the on-\(A\) gap vanishes.
- `vignette("Chapter-C05")` — the restricted two-block Gibbs minorization
  route to a total variation bound, which replaced this material in that
  chapter.
- `ELLIPSOID_TV_BOUND.md` — the full note, including the
  escape-probability bounds (§4), the drift/Lyapunov arguments that control
  \(P^n(x,A^c)\) explicitly, and the likelihood-specific specializations
  (§5–§8: multivariate-t, logit, probit, Poisson) that Sections 1–3 above
  deliberately factor out.
- `BLOCK_GIBBS_ERGODICITY_ING.md` §17 — the simpler, already-implemented
  diagnostic this note's Theorem 1 refines.
- `BLOCK_GIBBS_ERGODICITY.md`, `P_MATRIX_MARGINAL_PRECISIONS.md` — the
  spectral-radius and ellipsoid-score machinery Definition 1's safe region
  specializes to in the t-likelihood case.
