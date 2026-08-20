# Chapter C05 — Implementation Notes

Companion to `vignette("Chapter-C05")` (Total Variation Bounds for Restricted
Two-block Gibbs Samplers). That chapter proves the certificate exists; this
note works out what has to be *computed* to produce one, in what order, and
where the theory currently stops short of a number.

**Scope.** Design only. No function signatures are proposed here and no code
is drafted — §9 lists computational responsibilities, not an API. Section
references of the form "Lemma 13", "§5.2" point into Chapter C05.

---

## 0. What a certificate consists of

Theorem 2 delivers, for a tail-mass budget \(\delta\), a triple
\((C,\varepsilon,Q)\) with \(C=\widetilde C_d\), \(Q=N(\gamma^\star,\Sigma^\star)\),
and \(\varepsilon=\varepsilon_d\,Q(C)\), such that for \(\gamma\in C\)

\[
\bigl\|q_n(\gamma,\cdot\mid C)-\pi_\gamma\bigr\|_{TV}\ \le\ (1-\varepsilon)^n+\delta .
\]

Turning that into software means producing five numbers — \(\gamma^\star\),
\(\varepsilon(\gamma^\star)\), \(d\), \(\varepsilon\), \(n\) — from the model
objects the package already builds (`model_setup()`, `pfamily_list()`, the
hyper-designs \(H_j\), the random-effect precision \(P_b=\Psi^{-1}\), and the
population prior \((\mu_0,\Lambda_\gamma)\)).

Standing structural assumption, inherited from §3.1: \(\Psi\) is **known and
fixed**, so \(P_b\) and hence
\(P_{11}=\Lambda_\gamma+\sum_jH_j^\top P_bH_j\) are constant matrices. If
variance components are sampled, \(P_{11}\) moves and the entire construction
must be re-derived; nothing below survives that change unexamined.

---

## 1. The computational spine

In dependency order:

1. Assemble \(P_{11}^{\mathrm{RE}}=\sum_jH_j^\top P_bH_j\) and
   \(P_{11}=\Lambda_\gamma+P_{11}^{\mathrm{RE}}\); factor once. Gives
   \(\Sigma^\star=P_{11}^{-1}\) immediately, and the ceiling \(\kappa_0\) from
   the pencil \((\Lambda_\gamma,P_{11})\) (§4B.4).
2. Find \(\gamma^\star\), the fixed point of the mean map \(M\) (§2 below).
   Byproducts: \(\tilde J\), \(\kappa\), \(\Sigma_\pi\), \(\rho\).
3. Evaluate the scalar \(\varepsilon(\gamma^\star)\) (§3).
4. Choose \(d\) meeting the budget \(\pi_\gamma(\widetilde C_d^{\,c})\le\delta\)
   (§4; closed form in §4B.3 outside closure).
5. Evaluate \(\varepsilon=\varepsilon_d\,Q(\widetilde C_d)\) (§5).
6. Invert Theorem 1 for \(n\), optimizing over \(d\) (§6).

**Alternative set calibration (MVN reference).** `MVN_CALIBRATED_MINORIZATION_SET.md`
defines tail budget **\(\delta\)**, radius \(r(\delta)\), level \(d=r(\delta)^2/2\),
and ratio cut \(\widetilde C_d\); minorization uses \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\)
and \(\varepsilon=\varepsilon_d Q(\widetilde C_d)\). See that note §0 for notation.

Steps 1–2 are the bulk of the work and produce most of the inputs to 3–6.
Step 3 is the only place where a genuinely hard numerical problem remains
outside the Gaussian case; step 4, which used to be the other one, reduces to a
scalar bisection once \(\kappa_0\) is in hand.

---

## 2. Stage 1 — the refresh measure \(Q\)

### 2.1 The covariance is free

\(\Sigma^\star=P_{11}^{-1}\) is fixed by (Q1) and does not depend on
\(\gamma\), on \(d\), on \(\delta\), or on the likelihood. It is one Cholesky
of a matrix assembled directly from the hyper-designs and the priors. Half of
\(Q\) costs essentially nothing, and the factor is reused throughout (the
\(M\) solve, the \(P_{11}\)-norm in \(\Psi\), the \(\chi^2\) mass bound).

Note the asymmetry this creates: everything expensive is concentrated in the
*centre* of \(Q\).

### 2.2 The centre: \(\gamma^\star\) is an EM fixed point

\(\gamma^\star\) is the fixed point of the **conditional-mean** map, not of the
transition kernel (a kernel's fixed point would be \(\pi_\gamma\), a
distribution). With \(b_j(\gamma):=\mathbb E[\beta_j\mid\gamma,y]\),

\[
M(\gamma)=P_{11}^{-1}\Bigl(\Lambda_\gamma\mu_0+\sum_jH_j^\top P_b\,b_j(\gamma)\Bigr).
\]

Taking the expectation of the complete-data log posterior over
\(\beta\mid\gamma^{(k)},y\) and maximizing in \(\gamma\) gives first-order
condition

\[
\Bigl(\Lambda_\gamma+\sum_jH_j^\top P_bH_j\Bigr)\gamma
=\Lambda_\gamma\mu_0+\sum_jH_j^\top P_b\,b_j(\gamma^{(k)}),
\]

i.e. \(\gamma^{(k+1)}=M(\gamma^{(k)})\). **The Banach iteration of Lemma 3 is
exactly EM for the marginal posterior mode of \(\gamma\), with \(\beta\) as
missing data.** The E-step is the conditional means \(b_j\); the M-step is a
single solve against the pre-factored \(P_{11}\).

This is the organizing fact for the implementation, for four reasons.

**(a) Monotone ascent, hence convergence where Banach fails.** Appendix A.2
records that at \(\tau=0\), for families whose GLM weights vanish in the tails
(logit, probit at large \(|\eta|\)), \(\sup_\gamma\kappa(\gamma)=1\) exactly,
so "the contraction argument gives nothing globally" and Lemma 8 is brought in
as a variational replacement. EM's ascent property converges in precisely that
regime without needing a global contraction, so the natural implementation is
on firmer footing than Lemma 3's route. No line search is required.

**(b) \(\kappa\) is the fraction of missing information.** The Jacobian
\(J(\gamma)=P_{11}^{-1}\sum_jH_j^\top P_bV_j(\gamma)P_bH_j\), with
\(V_j=\mathrm{Cov}(\beta_j\mid\gamma,y)\), is \(I_{\mathrm{com}}^{-1}I_{\mathrm{mis}}\)
in the standard EM sense, since \(P_{11}\) is the complete-data information for
\(\gamma\). So \(\kappa=\lambda_{\max}(\tilde J)\) is the classical EM linear
rate and is a byproduct of the terminal E-step.

**(c) \(\Sigma_\pi\) is the inverse observed information.**
\(I_{\mathrm{obs}}=I_{\mathrm{com}}-I_{\mathrm{mis}}=P_{11}(I-J)\), which
reproduces the \(\Sigma_\pi=(P_{11}(I-J))^{-1}\) asserted under closure in
Appendix A.4. No pilot chain is needed to estimate it.

**(d) One eigendecomposition serves everything.** From the spectrum of
\(\tilde J\) at \(\gamma^\star\): \(\kappa=\lambda_{\max}\),
\(\rho=(1-\kappa)^{-1}\), \(\varepsilon(\gamma^\star)=\det(I+\tilde J)^{-1/2}\)
under closure, and the weights for the \(\pi_\gamma\)-tail computation of §4.

### 2.3 The mean-versus-mode trap

\(b_j\) is a conditional **mean**. Substituting the conditional **mode**
computes the joint mode of \(\pi(\gamma,\beta\mid y)\), which is a different
point: it is the PQL-versus-Laplace gap, exact only in the Gaussian case and
materially wrong for GLMMs with small \(n_j\), especially binary responses.

This is the single most dangerous error available here, because
\(\gamma^\star\) is simultaneously the centre of \(Q\), the maximizer of
\(\varepsilon\), and the zero of \(\Psi\). A mode-based \(b_j\) still produces
a convergent-looking iteration and a plausible certificate; every downstream
number is quietly displaced. The risk is concrete in this codebase, where
`two_block_mode_weights()` already takes a `b_mode` argument that must **not**
be reused for this purpose.

### 2.4 E-step tiers

By (H2) each \(\pi(\beta_j\mid\gamma,y)\) is log-concave, and the problem is
group-separable: \(J\) independent \(p_{\mathrm{re}}\)-dimensional
computations. Three tiers are worth supporting.

- **Exact (Gaussian).** \(b_j\) and \(V_j\) are closed-form solves. This tier
  should exist if only as the reference implementation for validating the
  others.
- **Adaptive Gauss–Hermite quadrature** (small \(p_{\mathrm{re}}\)).
  Deterministic, so the EM stopping rule stays a clean tolerance on
  \(\|\gamma^{(k+1)}-\gamma^{(k)}\|\); returns \(V_j\) in the same pass, which
  §2.2(b)–(d) need anyway.
- **Monte Carlo (MCEM).** Nearly free to build, because Block 1 of the
  existing two-block sampler already draws *exactly* from
  \(\pi(\beta_j\mid\gamma,y)\) (`rglmb`, `rNormalGLM_reg_group`). Averaging
  those draws is an E-step using tested machinery. The cost is that the
  stopping rule must then account for Monte Carlo error rather than treating
  successive iterates as exact — a ramping sample size across iterations
  (Booth–Hobert style) rather than a fixed tolerance.

The MC tier is what makes the general GLMM case tractable at all, and it
reuses code that already exists; it should probably be the default for
non-Gaussian families rather than an afterthought.

### 2.5 Acceleration

The regime that matters is the slow one. \(\kappa\to1\) means weakly
informative groups (the centred-parameterization problem), which is
simultaneously when EM crawls and when \(\rho=(1-\kappa)^{-1}\) inflates the
certified constant. Plain EM is not adequate there. Two options, both cheap:

- **SQUAREM / Aitken** applied to the EM sequence directly. Derivative-free,
  drops onto \(M\) with no extra model code.
- **Quasi-Newton on \(\Phi\).** Fisher's identity gives
  \(\nabla\Phi(\gamma)=P_{11}(\gamma-M(\gamma))\), so the same E-step that
  produces \(M\) produces the *exact* gradient of the marginal log-posterior.
  L-BFGS is therefore available at no additional cost per iteration.

The second is attractive because the gradient is exact rather than
approximated, but it forfeits the monotone-ascent guarantee of §2.2(a); the
safest default is accelerated EM with a monotonicity fallback.

### 2.6 Edge cases to handle explicitly

- **Flat population prior** (\(\Lambda_\gamma=0\)). \(P_{11}=P_{11}^{\mathrm{RE}}\),
  still \(\succ0\) under (H3a), so \(\Sigma^\star\) is fine and the \(\mu_0\)
  term drops from \(M\). Existence of \(\gamma^\star\) is Lemma 8 and needs
  (H3b) — group-wise estimability — so the preflight estimability gate is a
  hard prerequisite here, not advisory.
- **\(\kappa\) at or near 1.** Detect and report rather than iterate blindly;
  it is both a convergence problem and a signal that the eventual certificate
  will be weak.
- **(H3a) failure.** \(P_{11}^{\mathrm{RE}}\) singular means no \(\Sigma^\star\)
  and no certificate. This must be checked before anything else runs.

---

## 3. Stage 2 — the scalar \(\varepsilon(\gamma^\star)\)

This is the one hard quantity. Definition 5 makes it a \(q\)-dimensional
convex program, \(\varepsilon(\gamma)=\exp\{\inf_{\gamma'}g(\gamma'\mid\gamma)\}\),
over a Gaussian *mixture* density — so each objective evaluation itself needs
integration over \(\beta\).

- **Under Gaussian closure (LMM):** Lemma 22 gives
  \(\varepsilon(\gamma^\star)=\det(I+\tilde J)^{-1/2}\), free from the
  eigendecomposition of §2.2(d).
- **Outside closure:** the convex program must be solved numerically, with a
  Monte Carlo inner objective. Expensive and noisy. It is solved **once**, at
  \(\gamma^\star\) only — never on a grid — which is what keeps the overall
  scheme feasible.

**Do not reuse \((1+\kappa)^{-q/2}\) outside closure.** That inequality appears
*inside* Lemma 22 and is a closure result. Using it as a cheap GLM fallback
would produce an uncertified number wearing the same name.

---

## 4. Stage 3 — the set size \(d(\delta)\)

The requirement is \(\pi_\gamma(\widetilde C_d^{\,c})\le\delta\).

**The central difficulty, stated plainly.** This is a tail probability under
\(\pi_\gamma\), and outside the Gaussian case \(\pi_\gamma\) is exactly the
object we do not have. Corollary 15 and Lemma 14 simplify the **set** — they
replace \(\widetilde C_d\) by an ellipsoid — but they do **not** touch the
**measure**. An ellipsoid tail under an unknown \(\pi_\gamma\) is no more
computable than the original. So neither is a computational route on its own;
each must be paired with something that controls \(\pi_\gamma\) itself.

The way out is to **change the measure rather than the set**. §4B.0 shows
\(\pi_\gamma\) is the exponential tilt of the refresh Gaussian \(Q\) by
\(e^{\Psi}\), which converts the tail into an expectation under a known law.
Route (d) below is the consequence and is the recommended default outside
closure; §4A.7 supplies the one scalar it needs.

**(a) Closure — exact, no Monte Carlo.** Under closure \(\pi_\gamma=N(\gamma^\star,\Sigma_\pi)\)
with \(\Sigma_\pi=(P_{11}(I-J))^{-1}\), so the envelope of Corollary 15,

\[
\pi_\gamma(\widetilde C_d^{\,c})\ \le\ \pi_\gamma\bigl(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d\bigr),
\]

is the tail of a **weighted sum of \(\chi^2_1\) variates**, with weights the
eigenvalues of \((I-J)^{-1}\) — already in hand from §2.2(d), and with largest
weight exactly \(\rho=(1-\kappa)^{-1}\). Evaluate by Davies or Imhof and
root-find in \(d\). This gives a fully closed path with no sampling anywhere
and should be the first implementation target. §4A quantifies how much this
envelope gives away, and §4A.6 shows that \(P_{11}^{\mathrm{RE}}\) should be
used in place of \(P_{11}\) here.

**Package implementation (v1).** `deficiency_spectrum()` reports
\(\kappa_i=\mathrm{eig}(A)\) and \(w_i=\kappa_i/(1-\kappa_i)=\mathrm{eig}(B)\)
from the terminal `population_mode()` E-step (\(S=P_{11}\tilde J\),
\(B=A(I-A)^{-1}\) spectrally). `deficiency_calibrate()` inverts
\(\Pr(\sum_i w_i Z_i^2 > r^2)=\delta\) and sets \(d=r^2/2\). Limits:
\(\kappa_i\to0\Rightarrow w_i\to0\Rightarrow d\to0\);
\(\kappa_i\to1\Rightarrow w_i\to\infty\Rightarrow d\to\infty\).
Weighted-\(\chi^2\) tails use \pkg{CompQuadForm} or \pkg{mgcv} when installed.
The legacy `mvn_calibrate(q, delta)` equal-weight shortcut is retained for
reference only.

**(b) Empirical.** Evaluate \(\Psi\) at stored pilot draws and take the
empirical fraction with \(\Psi>d\). Requires \(\Psi\), which by the
complementarity identity (§7) needs the marginal log-posterior \(\bar\Phi\) —
for GLMMs that means a Laplace or AGHQ evaluation per draw, which is an
approximation and therefore degrades the certificate's status from proved to
estimated. This mirrors the empirical \(\hat p\) pattern already used in
`BLOCK_GIBBS_ERGODICITY_ING.md` §17 and inherits the same caveats, plus Monte
Carlo error on the fraction itself.

**(c) Sandwich.** Lemma 14 gives
\(\{\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\le2d\}\subseteq\widetilde C_d
\subseteq\{\|\gamma-\gamma^\star\|^2_{S_\flat}\le2d\}\) with
\(S_\flat=\sum_jH_j^\top P_b(P_b+G_j)^{-1}P_bH_j\). Prior-free, avoids
\(\bar\Phi\) entirely, conservative. Requires group curvature ceilings
\(-\nabla^2\ell_j\preceq G_j\), automatic for binomial with any link but **not**
available for Poisson-log — so this route is family-dependent and must refuse
to run where \(G_j\) does not exist.

*This is a set reduction, not a computation.* It removes the need for
\(\bar\Phi\) in describing \(\widetilde C_d\), which is worth having, but the
resulting ellipsoid tail is still under \(\pi_\gamma\). Pair it with route (d)
to obtain a number.

**(d) Exponential tilt — closed form, no sampling, any family.** Developed in
§4B. Given a uniform ceiling \(S(\gamma)\preceq\kappa_0P_{11}\),

\[
\pi_\gamma(\widetilde C_d^{\,c})\ \le\
\Bigl(\frac{2ed}{q\,\kappa_0}\Bigr)^{q/2}\exp\Bigl(-\frac{d(1-\kappa_0)}{\kappa_0}\Bigr),
\qquad d>\frac{q\,\kappa_0}{2(1-\kappa_0)} ,
\]

inverted for \(\delta\) by scalar bisection. \(\kappa_0\) is one generalized
eigenvalue when \(\Lambda_\gamma\succ0\) (§4B.4). This requires no
\(\bar\Phi\), no pilot chain, no quadrature and no closure, and its exponent is
the *exact* tail rate rather than the envelope rate (§4B.2). **It should be the
default route outside closure**, with §4B.5 available when the extra sharpness of
evaluating \(M_Q\) directly is worth the quadrature.

**(e) Gaussian majorization — tighter than (d), needs a \(\beta\)-restriction.**
Developed in `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md`. Where a lower bound
\(P_{22}\succeq P_{22}^{\mathrm{LB}}\) can be certified — which for logit,
probit and Poisson-log requires restricting \(\beta\) to a convex region — the
posterior is majorized on the escape set by a scaled
\(N(\gamma^\star,(\Pi^{\mathrm{LB}})^{-1})\) with
\(\Pi^{\mathrm{LB}}=P_{11}-P_{12}(P_{22}^{\mathrm{LB}})^{-1}P_{21}\), giving

\[
\pi_\gamma(\widetilde C_d^{\,c})\ \le\ M\,\Pr\Bigl(\sum_i\nu_iZ_i^2>2d\Bigr),
\qquad \nu_i=\frac{\kappa_i^{\mathrm{UB}}}{1-\kappa_i^{\mathrm{UB}}} .
\]

This keeps the exact weighted-\(\chi^2\) tail shape that route (d) Chernoffs
away, has no side condition, and is *exact* when the floor is tight. It costs a
second budget \(\pi(\beta\notin B)\), which is what
`RATE_Ac_BOUNDS_BY_LIKELIHOOD.md` bounds per likelihood. Prefer (e) over (d)
wherever a pooling-weight floor is available; (d) remains the route when only a
population prior is.

**What is not usable as stated:** Remark 6.1's
\(\pi_\gamma(\widetilde C_d^{\,c})\lesssim e^{-d/\rho}\) cannot certify a numeric
\(\delta\), since Appendix A.4 states plainly that "the multiplicative constant
suppressed by \(\lesssim\) is standard but not tracked here". It fixes the
*exponent* \(\rho\), which is the right object for the scaling claim
\(\varepsilon\gtrsim\varepsilon(\gamma^\star)\delta^\rho\) and for reporting, but
it is a heuristic, not a certificate. Any implementation offering it must label
it as such. Route (d) is the same statement with the constant tracked and a
sharper exponent, so it should be preferred wherever Remark 6.1 would have been
invoked.

---

## 4A. How good is the escape bound? Exact comparison under closure

Route 4(a) uses Corollary 15's envelope rather than the true escape
probability. Under Gaussian closure both are available in closed form, so the
slack can be quantified exactly rather than guessed at.

### 4A.0 Precision-matrix form

Everything below is cleanest in the joint precision blocks of the two-block
state \((\gamma,\beta)\), which is also the `Chapter-C03` idiom. Write

\[
A:=P_{11},\qquad
B:=P_{12}P_{22}^{-1}P_{21},\qquad
\Pi:=P_{11}-P_{12}P_{22}^{-1}P_{21},
\]

so \(A\) is the **conditional** precision of \(\gamma\) given \(\beta\), \(\Pi\)
is the **marginal (posterior)** precision of \(\gamma\), and \(B\) is the Schur
correction relating them. Trivially but importantly,

\[
\boxed{\ A=\Pi+B\ }
\]

— the Schur complement identity, rearranged.

**\(B\) is exactly the \(S\) of Lemma 14.** With
\(P_{12}=[-H_1^\top P_b,\dots,-H_J^\top P_b]\) and \(P_{22}\) block-diagonal
with blocks \(P_b+Z_j^\top W_jZ_j\), the signs cancel and the block structure
turns the Schur product into a sum over groups:

\[
B=P_{12}P_{22}^{-1}P_{21}=\sum_jH_j^\top P_bV_jP_bH_j=S,
\qquad V_j=\bigl(P_b+Z_j^\top W_jZ_j\bigr)^{-1}=\mathrm{Cov}(\beta_j\mid\gamma,y).
\]

So the information decomposition of §4A.1 is nothing but the Schur identity:
\(I_{\mathrm{com}}=A\), \(I_{\mathrm{obs}}=\Pi\), \(I_{\mathrm{mis}}=B\).

**The three precisions have three distinct roles:**

| object | precision | role |
|---|---|---|
| refresh measure \(Q\) | \(A=P_{11}\) | \(\Sigma^\star=A^{-1}\) |
| posterior \(\pi_\gamma\) | \(\Pi=P_{11}-P_{12}P_{22}^{-1}P_{21}\) | law the tail is taken under |
| certified set \(\widetilde C_d\) | \(B=P_{12}P_{22}^{-1}P_{21}\) | \(\Psi=\tfrac12\|\gamma-\gamma^\star\|^2_B\) |

**Answering the question directly: the envelope uses \(A=P_{11}\) where the
exact set uses \(B=P_{12}P_{22}^{-1}P_{21}\).** Compared with the posterior
precision it is larger by exactly the Schur term, \(A=\Pi+B\succeq\Pi\), so the
envelope is a smaller (conservative) ellipsoid — correct in direction, and
loose by precisely \(B\).

The relation to §4A.2–§4A.4's \(\kappa\) is a change of reference matrix:
\(\kappa_i=\mathrm{eig}(A^{-1}B)\) measures missing information against
*complete*, while

\[
\mu_i:=\mathrm{eig}(\Pi^{-1}B)=\frac{\kappa_i}{1-\kappa_i}
\]

measures it against *observed*. The \(\mu_i\) are the natural coordinates here,
because \(\Pi^{-1}A=I+\Pi^{-1}B\) gives the weights of the two quadratic forms
under \(\gamma-\gamma^\star\sim N(0,\Pi^{-1})\) as

\[
\underbrace{\pi_\gamma(\widetilde C_d^{\,c})=\Pr\Bigl(\sum_i\mu_i Z_i^2>2d\Bigr)}_{\text{exact, metric }B}
\qquad\le\qquad
\underbrace{\Pr\Bigl(\sum_i(1+\mu_i)Z_i^2>2d\Bigr)}_{\text{envelope, metric }A} .
\]

**The envelope adds exactly \(1\) to every weight** — that is the entire
approximation, and it is immediate from \(A=\Pi+B\). Two consequences drop out
without further work. First, \(\rho=1+\mu_{\max}\) is the envelope's largest
weight, which is what Remark 6.1's exponent is. Second, since each weight
enters \(d\) linearly, adding \(1\) to all of them costs

\[
d_{\text{env}}-d_{\text{exact}}=\frac{t_\delta}{2},\qquad t_\delta=\chi^2_{q,1-\delta},
\]

which is why the additive penalty of §4A.4 is free of \(\kappa\) — a fact that
looks like a coincidence in the \(\kappa\) parameterization and is obvious in
this one. In the univariate case \(A,B,\Pi\) are scalars,
\(B=P_{12}^2/P_{22}\), \(\mu=B/\Pi\), and the two tails are
\(2\Phi(-\sqrt{2d/\mu})\) and \(2\Phi(-\sqrt{2d/(1+\mu)})\).

Finally, §4A.6's improvement reads: replace \(A=P_{11}\) by
\(A-\Lambda_\gamma=P_{11}^{\mathrm{RE}}\), legitimate because
\(B\preceq P_{11}^{\mathrm{RE}}\) by Brascamp–Lieb. In precision terms, strip
the population prior's contribution out of the conditional precision before
using it as the envelope metric.

### 4A.1 Where the slack comes from

**First, a framing caution.** Both the exact and the envelope probability are
taken under the *same* measure, the marginal posterior \(\pi_\gamma\) with
precision \(\Pi\); what changes between them is the **set**, i.e. the metric
(\(B\) versus \(P_{11}\)). It is not a comparison of the conditional against
the marginal law.

The conditional-versus-marginal comparison instead defines the *set*. Since
\(Q=N(\gamma^\star,P_{11}^{-1})\) carries exactly the conditional precision,
\(\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}\) is the recentred negative
log-density under \(Q\) and \(\bar\Phi\) the recentred negative log-density
under \(\pi_\gamma\), so Lemma 13(3) read as densities gives

\[
\boxed{\ \Psi(\gamma)=\log\frac{\pi_\gamma(\gamma)}{q_Q(\gamma)}
-\log\frac{\pi_\gamma(\gamma^\star)}{q_Q(\gamma^\star)}\ }
\]

— **general, not a closure result**. \(\widetilde C_d\) is a sublevel set of
the log-ratio of the marginal posterior to the conditional-precision Gaussian,
normalized at \(\gamma^\star\); since \(\Pi\preceq P_{11}\) the posterior is
the more diffuse of the two, so certification is lost exactly where the
posterior has outrun the refresh measure by a factor \(e^{d}\).

The other measure does appear, but for the other factor: Corollary 15's
\(Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)\) feeds the mass discount in
\(\varepsilon=\varepsilon_d\,Q(\widetilde C_d)\). So \(\pi_\gamma\) sets the
tail budget \(\delta\) and \(Q\) sets the discount. With
\(\mu_i=\mathrm{eig}(\Pi^{-1}B)\), the three combinations that arise are:

| measure | metric | weights of \(\sum_iw_iZ_i^2\) | used for |
|---|---|---|---|
| \(Q\) | \(P_{11}\) | \(1\) | mass discount \(Q(\widetilde C_d)\) |
| \(\pi_\gamma\) | \(B\) | \(\mu_i\) | exact escape probability |
| \(\pi_\gamma\) | \(P_{11}\) | \(1+\mu_i\) | envelope |

The envelope row is the sum of the other two — the clearest statement of what
the bound gives away.

**Now the algebra.** The envelope is the complementarity identity (§7) with
\(\bar\Phi\) discarded:

\[
\Psi(\gamma)=\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}-\bar\Phi(\gamma)
\ \le\ \tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}},
\]

so the discarded quantity is the *entire* recentred negative log-posterior.
Under closure \(M\) is affine with \(M(\gamma)-\gamma^\star=J(\gamma-\gamma^\star)\)
and \(\bar\Phi(\gamma)=\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}(I-J)}\), so
\(\Psi\) is exactly quadratic:

\[
\boxed{\ \Psi(\gamma)=\tfrac12\|\gamma-\gamma^\star\|^2_{S},
\qquad S=P_{11}J=\sum_jH_j^\top P_bV_jP_bH_j\ }
\]

where \(V_j(\gamma)=\mathrm{Cov}(\beta_j\mid\gamma,y)\) and \(S=\nabla^2\Psi\)
is Lemma 14's symbol. Under closure \(V_j\) is constant, so \(S\) is a constant
matrix and \(\Psi\) is exactly quadratic.

**The three metrics in play form one decomposition.** In the EM reading of
§2.2, \(S\) is the *missing information* matrix, and

\[
\underbrace{P_{11}}_{I_{\mathrm{com}}}
=\underbrace{\Sigma_\pi^{-1}}_{I_{\mathrm{obs}}}
+\underbrace{S}_{I_{\mathrm{mis}}} ,
\]

with \(P_{11}\) the metric of the refresh measure \(Q\), \(\Sigma_\pi^{-1}\)
the metric of the posterior \(\pi_\gamma\), and \(S\) the metric of the
certified set \(\widetilde C_d\). The complementarity identity is exactly this
decomposition read as quadratic forms: \(\Psi\) is the missing-information
quadratic, \(\bar\Phi\) the observed-information quadratic. It also identifies
the \(\kappa_i\) as the generalized eigenvalues of \((S,P_{11})\) — the
fraction of information missing along each direction, necessarily in
\([0,1)\) — and makes \(S\preceq P_{11}^{\mathrm{RE}}\preceq P_{11}\) (§4A.6)
the statement that missing information cannot exceed complete information.

So

\[
\widetilde C_d=\{\|\gamma-\gamma^\star\|^2_{S}\le 2d\}
\quad\text{(exact)},
\qquad
\{\|\gamma-\gamma^\star\|^2_{P_{11}}\le 2d\}\subseteq\widetilde C_d
\quad\text{(envelope)} .
\]

**The envelope replaces \(S=P_{11}J\) by \(P_{11}\), i.e. replaces \(J\) by
\(I\): it assumes the fraction of missing information is \(1\).** That is the
whole content of the approximation, and it predicts the behaviour below —
tight as \(\kappa\to1\), useless as \(\kappa\to0\).

### 4A.2 Univariate case

With \(q=1\), \(J=\kappa\), and \(\pi_\gamma=N(\gamma^\star,\sigma_\pi^2)\),
\(\sigma_\pi^2=\bigl(P_{11}(1-\kappa)\bigr)^{-1}\).

Both sets here are **symmetric intervals** about \(\gamma^\star\): \(\Psi\) is
an even quadratic in \(\gamma-\gamma^\star\), so

\[
\widetilde C_d=\Bigl\{|\gamma-\gamma^\star|\le\sqrt{2d/S}\Bigr\},\quad S=P_{11}\kappa,
\qquad
\text{envelope}=\Bigl\{|\gamma-\gamma^\star|\le\sqrt{2d/P_{11}}\Bigr\},
\]

and since \(\pi_\gamma\) is Gaussian centred at \(\gamma^\star\), each escape
probability is a **two-sided normal tail** with equal halves. In standardized
units the radii are \(\sqrt{2d(1-\kappa)/\kappa}\) and \(\sqrt{2d(1-\kappa)}\)
posterior standard deviations respectively, so

\[
\underbrace{\pi_\gamma(\widetilde C_d^{\,c})
=2\,\Phi\!\Bigl(-\sqrt{\tfrac{2d(1-\kappa)}{\kappa}}\Bigr)}_{\text{exact}}
\qquad\le\qquad
\underbrace{2\,\Phi\!\bigl(-\sqrt{2d(1-\kappa)}\bigr)}_{\text{envelope}} .
\]

Equivalently, using \(\chi^2_1=Z^2\) and \(\Pr(\chi^2_1>t)=2\Phi(-\sqrt t)\),

\[
\underbrace{\Pr\!\Bigl(\chi^2_1>\tfrac{2d(1-\kappa)}{\kappa}\Bigr)}_{\text{exact}}
\qquad\le\qquad
\underbrace{\Pr\bigl(\chi^2_1>2d(1-\kappa)\bigr)}_{\text{envelope}} ,
\]

which is the form that generalizes in §4A.3. The two differ only in the
argument, by the factor \(1/\kappa\) — equivalently, the envelope shrinks the
certified interval by a factor \(\sqrt\kappa\) in posterior standard
deviations.

In this notation the budget calculation reads directly off the two-sided
normal critical value \(z:=z_{1-\delta/2}\): the exact and envelope radii
must reach \(z\) standard deviations, giving \(d_{\text{exact}}=\kappa
z^2/(2(1-\kappa))\), \(d_{\text{env}}=z^2/(2(1-\kappa))\), and the constant
gap \(z^2/2\) of §4A.4 (with \(z^2=\chi^2_{1,1-\delta}=t_\delta\)).

Note the two-sided reading is special to \(q=1\), where an ellipsoid
degenerates to an interval; for \(q>1\) the complement is the exterior of an
ellipsoid and has no directional decomposition.

### 4A.3 Multivariate case

Let \(\kappa_1,\dots,\kappa_q\) be the eigenvalues of \(\tilde J\) (so
\(\kappa=\max_i\kappa_i\)) and \(Z_i\) iid \(N(0,1)\). With
\(\Sigma_\pi=(P_{11}(I-J))^{-1}\), the exact statistic has weights
\(\mathrm{eig}(\Sigma_\pi S)=\mathrm{eig}\bigl((I-J)^{-1}J\bigr)\) and the
envelope's has \(\mathrm{eig}(\Sigma_\pi P_{11})=\mathrm{eig}\bigl((I-J)^{-1}\bigr)\):

\[
\underbrace{\pi_\gamma(\widetilde C_d^{\,c})
=\Pr\Bigl(\sum_{i=1}^q\frac{\kappa_i}{1-\kappa_i}Z_i^2>2d\Bigr)}_{\text{exact}}
\qquad\le\qquad
\underbrace{\Pr\Bigl(\sum_{i=1}^q\frac{1}{1-\kappa_i}Z_i^2>2d\Bigr)}_{\text{envelope}} .
\]

Both are weighted-\(\chi^2\) tails on the *same* eigenbasis, evaluable by
Imhof or Davies. The relation between them is exactly

\[
w_i^{\text{env}}=\frac{w_i^{\text{exact}}}{\kappa_i},
\]

each weight inflated by the reciprocal of its own eigenvalue. The envelope's
largest weight is \(\rho=(1-\kappa)^{-1}\), which is where the \(\rho\) of
Remark 6.1 comes from.

### 4A.4 How the two compare

Two different questions, with two different answers; both are worth reporting.

**At fixed \(d\), the envelope is exponentially loose.** From the \(\chi^2_1\)
tail asymptotic,

\[
\frac{\text{envelope}}{\text{exact}}\ \approx\
\kappa^{-1/2}\exp\Bigl(d\,\frac{(1-\kappa)^2}{\kappa}\Bigr),
\]

growing exponentially in \(d\) at rate \((1-\kappa)^2/\kappa\).

**At fixed \(\delta\), the cost is a constant.** Solving each tail \(=\delta\)
for \(d\), with \(t_\delta:=\chi^2_{q,1-\delta}\) (and \(q=1\), or all
\(\kappa_i\) equal):

\[
d_{\text{exact}}=\frac{\kappa\,t_\delta}{2(1-\kappa)},
\qquad
d_{\text{env}}=\frac{t_\delta}{2(1-\kappa)},
\qquad
\frac{d_{\text{env}}}{d_{\text{exact}}}=\frac1\kappa,
\]

so the ratio blows up as \(\kappa\to0\) — but the *difference* is

\[
\boxed{\ d_{\text{env}}-d_{\text{exact}}=\frac{t_\delta}{2}
\quad\text{— independent of }\kappa .\ }
\]

Since \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\), the envelope therefore
costs a **fixed multiplicative factor** in the minorization constant,

\[
\frac{\varepsilon_d^{\text{env}}}{\varepsilon_d^{\text{exact}}}
=e^{-t_\delta/2}=e^{-\chi^2_{q,1-\delta}/2},
\]

with no \(\kappa\) dependence at all (slightly offset in the assembled
\(\varepsilon=\varepsilon_dQ(\widetilde C_d)\), since the larger \(d\) pushes
\(Q(\widetilde C_d)\) closer to \(1\)).

The two statements reconcile because \(d_{\text{exact}}\) itself scales like
\(\kappa/(1-\kappa)\): the exponential looseness in the *probability* is
absorbed into a constant additive shift in \(d\).

### 4A.5 Reading

The envelope is **tightest exactly where mixing is worst**. As \(\kappa\to1\)
the rate \((1-\kappa)^2/\kappa\to0\) and \(S\to P_{11}\), so the envelope
becomes exact — but that is the regime where \(\rho=(1-\kappa)^{-1}\) is
already inflating everything. As \(\kappa\to0\) the envelope is catastrophic
in relative terms: the true \(\widetilde C_d\) is enormous (\(\Psi\approx0\)
everywhere) and the true escape probability nearly zero, while the envelope
remains a fixed \(\chi^2\) tail. So easy problems get a needlessly conservative
\(d\), and the certificate is worse than it needs to be precisely where the
sampler is good.

### 4A.6 A strictly better envelope is already available

Corollary 15's \(P_{11}\) envelope is **not** the tightest inner ellipsoid the
chapter supplies. Lemma 14's inner set uses \(P_{11}^{\mathrm{RE}}\), and
Brascamp–Lieb (\(V_j\preceq P_b^{-1}\)) gives

\[
S=\sum_jH_j^\top P_bV_jP_bH_j\ \preceq\ \sum_jH_j^\top P_bH_j=P_{11}^{\mathrm{RE}}
\ \preceq\ P_{11}=\Lambda_\gamma+P_{11}^{\mathrm{RE}} .
\]

A smaller metric means a larger ellipsoid, hence a better inner approximation
to \(\widetilde C_d\) and a tighter escape bound. So

\[
\pi_\gamma(\widetilde C_d^{\,c})\ \le\
\pi_\gamma\bigl(\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}>2d\bigr)\ \le\
\pi_\gamma\bigl(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d\bigr),
\]

with equality throughout when \(\Lambda_\gamma=0\), and strict improvement
whenever the population prior is proper. The \(P_{11}^{\mathrm{RE}}\) version
is the same weighted-\(\chi^2\) computation with weights
\(\mathrm{eig}(\Sigma_\pi P_{11}^{\mathrm{RE}})\); it has no single-parameter
closed form, but costs nothing extra to evaluate. **Route 4(a) should use
\(P_{11}^{\mathrm{RE}}\), not \(P_{11}\).**

### 4A.7 What survives outside closure: strong log-concavity

Everything in §4A.1–§4A.6 evaluates a tail under \(\pi_\gamma\), which is only
available in the Gaussian case. The general case needs a handle on
\(\pi_\gamma\) itself, and Lemma 5 supplies one that is often overlooked:

\[
\nabla^2\Phi(\gamma)=P_{11}-\sum_jH_j^\top P_bV_j(\gamma)P_bH_j=P_{11}-S(\gamma)=\Pi(\gamma)\ \succeq\ \Lambda_\gamma
\quad\text{for every }\gamma,
\]

so **\(\pi_\gamma\) is \(\Lambda_\gamma\)-strongly log-concave**, for every GLM
family, with no closure assumption. Two consequences:

- **Brascamp–Lieb:** \(\Sigma_\pi\preceq\Lambda_\gamma^{-1}\), so
  \(\rho=\lambda_{\max}(\Sigma_\pi P_{11})\le\lambda_{\max}(\Lambda_\gamma^{-1}P_{11})\)
  is computable without knowing \(\pi_\gamma\).
- **Explicit tails:** strongly log-concave measures admit sub-Gaussian
  concentration for quadratic forms with *explicit* constants, unlike the
  general log-concave setting (Borell) that Remark 6.1 appeals to.

This is the substantive refinement of §8(1): the untracked constant there is
untracked because Remark 6.1 invokes general log-concavity. Under (H2) the
measure is *strongly* log-concave, where the relevant constants are standard.
The gap is more closable than the chapter implies.

**Prior-free version.** The above degenerates as \(\Lambda_\gamma\downarrow0\),
which is precisely the limit the chapter is built to survive. The repair is the
uniform curvature margin already used by Lemma 3(c) and Appendix A.2: if
\(-\nabla^2\ell_j\succeq\epsilon_jP_b\) uniformly, then
\(\kappa_0\le(1+\min_j\epsilon_j)^{-1}<1\) prior-free, hence
\(S(\gamma)\preceq\kappa_0P_{11}\) and

\[
\boxed{\ \nabla^2\Phi(\gamma)\succeq(1-\kappa_0)P_{11}\succ0
\quad\text{uniformly in }\gamma,\ \text{free of }\Lambda_\gamma.\ }
\]

So \(\Sigma_\pi\preceq\bigl((1-\kappa_0)P_{11}\bigr)^{-1}\) and
\(\rho\le(1-\kappa_0)^{-1}\).

**The practical upshot.** The closure formulas of §4A.2–§4A.4 become valid
*bounds* in general on replacing the spectrum \(\{\kappa_i\}\) by the single
uniform ceiling \(\kappa_0\) — e.g. the envelope tail becomes
\(\Pr\bigl(\chi^2_q>2d(1-\kappa_0)\bigr)\), and the \(d\)-for-\(\delta\)
arithmetic goes through unchanged with \(\kappa\mapsto\kappa_0\). That is a
computable, certified, conservative route for non-Gaussian families, and it
does not require \(\bar\Phi\), a pilot chain, or knowledge of \(\pi_\gamma\).

**What remains open.** The step from the covariance bound to the tail carries a
concentration constant that is explicit for strongly log-concave measures but
is not written down anywhere in Chapter C05, and a strongly log-concave measure
is not Gaussian, so the \(\chi^2\) form above is the right shape rather than a
finished inequality. Pinning that constant is the single highest-value piece of
theory work for making this pipeline certified outside the Gaussian case — and
it is a much better-posed problem than the general log-concave version.

But see §4B: the concentration step can be bypassed entirely. The ceiling
\(\kappa_0\) introduced here is exactly the input §4B needs, and routed through
the tilting identity it yields a closed-form tail bound with every constant
explicit, rather than a \(\chi^2\) "right shape". §4B supersedes the closing
paragraph above.

---

## 4B. A sampling-free bound: the exponential tilt

> **Status: derived here, not in Chapter C05.** §4B.0–§4B.1 follow from Lemma 13
> by elementary steps; the algebra of §4B.3, the closure MGF of §4B.2, the
> tilting identity of §4B.0 and the domination claim have all been checked
> numerically (`data-raw/_chk_tilt_bound.R`, results in §4B.9). None of it has
> been through the chapter's proof discipline, so §4B.3 remains the claim to
> review first.

### 4B.0 The posterior is an exponential tilt of \(Q\)

Complementarity (Lemma 13(3), restated in §7) says
\(\Psi=\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}-\bar\Phi\). Both terms are
log-density ratios: the first is *exactly* the recentred negative log-density of
the refresh measure,
\(\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}=\log\bigl(q_Q(\gamma^\star)/q_Q(\gamma)\bigr)\),
and the second is
\(\bar\Phi(\gamma)=\log\bigl(\pi_\gamma(\gamma^\star)/\pi_\gamma(\gamma)\bigr)\).
Subtracting,

\[
\boxed{\ \frac{\pi_\gamma(\gamma)}{q_Q(\gamma)}
=\frac{\pi_\gamma(\gamma^\star)}{q_Q(\gamma^\star)}\;e^{\Psi(\gamma)}\ }
\]

— \(\pi_\gamma\) is \(Q\) **exponentially tilted by \(e^{\Psi}\)**. No closure,
no Gaussian assumption, no family restriction, and the posterior normalizing
constant cancels between the two ratios. Consequently, for measurable \(A\),

\[
\pi_\gamma(A)=\frac{\mathbb E_Q\bigl[e^{\Psi}\mathbf 1_A\bigr]}{\mathbb E_Q\bigl[e^{\Psi}\bigr]},
\qquad\text{in particular}\qquad
\pi_\gamma(\widetilde C_d^{\,c})
=\frac{\mathbb E_Q\bigl[e^{\Psi}\mathbf 1\{\Psi>d\}\bigr]}{\mathbb E_Q\bigl[e^{\Psi}\bigr]} .
\]

This is the reduction §4 was missing. The escape probability is no longer an
expectation under the unknown \(\pi_\gamma\) but under
\(Q=N(\gamma^\star,P_{11}^{-1})\), which is known in closed form and already in
hand from Stage 1. The unknown-ness has been pushed entirely into the scalar
function \(\Psi\).

Write \(M_Q(\lambda):=\mathbb E_Q\bigl[e^{\lambda\Psi}\bigr]\). Two facts are
free: \(M_Q(1)=q_Q(\gamma^\star)/\pi_\gamma(\gamma^\star)<\infty\), so the
denominator above is never degenerate; and \(\Psi\ge0\) everywhere (Lemma 13(2):
\(\varepsilon\) attains its global maximum at \(\gamma^\star\)), so
\(M_Q(\lambda)\ge1\) for every \(\lambda\ge0\).

### 4B.1 A rigorous Chernoff bound

For \(\lambda>1\), the elementary inequality
\(\mathbf 1\{\Psi>d\}\le e^{(\lambda-1)(\Psi-d)}\) applied inside the numerator
gives \(\mathbb E_Q[e^{\Psi}\mathbf 1\{\Psi>d\}]\le e^{-(\lambda-1)d}M_Q(\lambda)\),
hence

\[
\boxed{\ \pi_\gamma(\widetilde C_d^{\,c})\ \le\
\inf_{\lambda>1}\ e^{-(\lambda-1)d}\,\frac{M_Q(\lambda)}{M_Q(1)}
\ \le\ \inf_{\lambda>1}\ e^{-(\lambda-1)d}\,M_Q(\lambda)\ }
\]

the second form using \(M_Q(1)\ge1\) to drop a denominator that costs a
normalizing constant to evaluate. Inverting the second form for the budget,

\[
d(\delta)\ =\ \inf_{\lambda>1}\ \frac{\log M_Q(\lambda)+\log(1/\delta)}{\lambda-1}.
\]

\(M_Q\) is the moment generating function of \(\Psi\) under a **known Gaussian**.
No posterior draws, no marginal likelihood, no untracked constant. Note the cost
structure: \(M_Q(\lambda)\) does not depend on \(d\) or \(\delta\), so it is
computed once on a \(\lambda\)-grid and reused for every budget.

### 4B.2 The bound is rate-exact

Under closure \(\Psi=\tfrac12\|\gamma-\gamma^\star\|^2_B\) with \(B=S\) constant,
and under \(Q\) the form \(\|\gamma-\gamma^\star\|^2_B\) is a weighted \(\chi^2_1\)
sum with weights \(\kappa_i=\mathrm{eig}(A^{-1}B)\) (§4A.0). So

\[
M_Q(\lambda)=\prod_i(1-\lambda\kappa_i)^{-1/2},\qquad \lambda<1/\kappa,
\]

and letting \(\lambda\uparrow1/\kappa\) the bound decays at rate
\(e^{-d(1-\kappa)/\kappa}=e^{-d/\mu_{\max}}\) — **exactly** the true tail rate of
the weighted-\(\chi^2\) escape probability of §4A.3. The bound is therefore
lossy only in the polynomial prefactor, not in the exponent.

Worth contrasting with Remark 6.1, whose exponent is \(d/\rho=d/(1+\mu_{\max})\)
— the *envelope* rate of §4A.0, obtained after adding \(1\) to every weight.
The tilting bound achieves the exact metric-\(B\) rate instead, so it is both
sharper in the exponent than Remark 6.1 and explicit where Remark 6.1 is not.

### 4B.3 Fully explicit form via a uniform ceiling

Suppose \(S(\gamma)\preceq\kappa_0P_{11}\) uniformly in \(\gamma\) — the ceiling
of §4A.7. Then \(\nabla^2\Phi\succeq(1-\kappa_0)P_{11}\), and since
\(\nabla\Phi(\gamma^\star)=0\), Taylor gives
\(\bar\Phi(\gamma)\ge\tfrac12(1-\kappa_0)\|\gamma-\gamma^\star\|^2_{P_{11}}\).
Substituting into complementarity yields a **pointwise ceiling on \(\Psi\)
itself**:

\[
\Psi(\gamma)\ \le\ \tfrac12\kappa_0\|\gamma-\gamma^\star\|^2_{P_{11}}
\qquad\text{for every }\gamma. \tag{*}
\]

Under \(Q\) the right side is \(\tfrac12\kappa_0\chi^2_q\), so
\(M_Q(\lambda)\le(1-\lambda\kappa_0)^{-q/2}\) for \(\lambda<1/\kappa_0\), and the
\(\lambda\)-optimization is elementary — the optimum is
\(\lambda^\star=\kappa_0^{-1}-q/(2d)\), giving

\[
\boxed{\ \pi_\gamma(\widetilde C_d^{\,c})\ \le\
\Bigl(\frac{2ed}{q\,\kappa_0}\Bigr)^{q/2}\exp\Bigl(-\frac{d(1-\kappa_0)}{\kappa_0}\Bigr),
\qquad d>\frac{q\,\kappa_0}{2(1-\kappa_0)}\ }
\]

with the side condition being exactly \(\lambda^\star>1\); for smaller \(d\) the
bound is vacuous, which is honest — small sets carry no certificate.

This is the object §4 has been missing: a closed-form, certified, sampling-free
bound on the escape probability, valid for any GLM family, requiring **one
scalar** \(\kappa_0\) and nothing else. It needs no \(\bar\Phi\), no pilot
chain, no quadrature, and no knowledge of \(\pi_\gamma\). It also settles §8(1)
— the previously untracked multiplicative constant is now written down, as the
explicit polynomial factor \((2ed/q\kappa_0)^{q/2}\).

Root-find in \(d\) for a target \(\delta\); the right side is decreasing in
\(d\) beyond the side condition, so the inversion is a scalar bisection.

### 4B.4 Where \(\kappa_0\) comes from

Two suppliers, in order of preference:

1. **Proper population prior (universal).** \(S\preceq P_{11}^{\mathrm{RE}}
   =P_{11}-\Lambda_\gamma\) by Brascamp–Lieb (§4A.6), so
   \[
   \kappa_0\ \le\ 1-\lambda_{\min}\bigl(P_{11}^{-1}\Lambda_\gamma\bigr)\ <\ 1
   \qquad\text{whenever }\Lambda_\gamma\succ0,
   \]
   a single generalized eigenvalue of the pencil \((\Lambda_\gamma,P_{11})\),
   already factored at Stage 1. **No family restriction and no curvature
   margin.** This is the default route.
2. **Uniform curvature margin (family-dependent).** If
   \(-\nabla^2\ell_j\succeq\epsilon_jP_b\) uniformly then
   \(\kappa_0\le(1+\min_j\epsilon_j)^{-1}\), prior-free. As §4A.7 notes, this is
   unavailable for logit, probit and Poisson-log, whose GLM weights vanish in
   the tails, so it must not be the primary route.

Under closure \(S\) is constant and \(\kappa_0=\kappa\) exactly, recovering the
rate-exactness of §4B.2. As \(\Lambda_\gamma\downarrow0\) route 1 degrades
gracefully — \(\kappa_0\uparrow1\) and the exponent \((1-\kappa_0)/\kappa_0\to0\)
— rather than failing outright, which is the correct behaviour: the flat limit
genuinely has no certified rate absent a curvature margin.

### 4B.5 Sharper: evaluate \(M_Q\) instead of bounding it

\((\ast)\) is a worst-case ceiling on \(\Psi\), so §4B.3 gives away whatever
\(\Psi\) is smaller than \(\tfrac12\kappa_0\|\cdot\|^2_{P_{11}}\). Evaluating
\(M_Q(\lambda)=\mathbb E_Q[e^{\lambda\Psi}]\) directly recovers that. Two ways,
neither of which touches the posterior:

- **Deterministic quadrature.** Gauss–Hermite in the \(Q\) metric via the
  Cholesky factor of \(P_{11}\). Tensor grids are fine for small \(q\); by
  \(q\approx5\) a sparse grid is needed. Yields a deterministic bound.
- **Independent \(Q\)-draws.** Draw iid \(\gamma^{(m)}\sim Q\) and average
  \(e^{\lambda\Psi}\). This is *not* posterior sampling: iid, exact, no
  convergence question, and with a bounded-below integrand a Bernstein interval
  applies. But the result is a confidence statement, not a proved bound, and the
  variance is finite only when \(M_Q(2\lambda)<\infty\) — under closure
  \(\lambda<1/(2\kappa)\), which bites in the same regime everything else does.

Either way \(\Psi\) is needed at each node, which is where §7's third identity
earns its place: \(\bar\Phi\) comes from a one-dimensional integral of the exact
gradient, so no Laplace approximation of the full marginal is required. And the
node values of \(\Psi\) are reused across the whole \(\lambda\)-grid, so the
\(\lambda\)-optimization and the sweep over \(\delta\) are both free once the
nodes are evaluated.

### 4B.6 Byproduct: a better \(\chi^2\) floor

\((\ast)\) rearranges to
\(\{\|\gamma-\gamma^\star\|^2_{P_{11}}\le2d/\kappa_0\}\subseteq\widetilde C_d\),
an inscribed ellipsoid larger by a factor \(1/\kappa_0\) in squared radius than
Corollary 15's. So §5's floor improves to

\[
Q(\widetilde C_d)\ \ge\ \Pr\bigl(\chi^2_q\le2d/\kappa_0\bigr)\ \ge\ \Pr\bigl(\chi^2_q\le2d\bigr),
\]

which feeds straight through to a larger certified \(\varepsilon\) at no cost
once \(\kappa_0\) is computed. Stage 4 should use it.

### 4B.7 Consistency check: the tilt is where \(\mu=\kappa/(1-\kappa)\) comes from

The identity explains a transform that looks unmotivated in §4A.0. Tilting a
standard normal coordinate by \(e^{\kappa_iz^2/2}\) produces
\(N\bigl(0,(1-\kappa_i)^{-1}\bigr)\); applying this to \(\Psi=\tfrac12\sum_i\kappa_iZ_i^2\)
under \(Q\) converts the weights from \(\kappa_i\) to
\(\kappa_i/(1-\kappa_i)=\mu_i\) under \(\pi_\gamma\). So the \(\kappa\mapsto\mu\)
change of reference matrix in §4A.0 *is* the exponential tilt, and the two
sections are describing the same object from opposite ends.

### 4B.8 What still needs checking

- **Finiteness of \(M_Q(\lambda)\) for \(\lambda>1\)** is the crux, and is not
  automatic: \(M_Q(\lambda)\propto\int\pi_\gamma^\lambda q_Q^{1-\lambda}\) is a
  Rényi divergence of order \(\lambda\), finite only if \(\pi_\gamma\)'s tails
  are not too much heavier than \(Q\)'s — and \(\Pi\preceq A\) says they are
  heavier. At \(\lambda=1\) it is finite but the bound is vacuous. §4B.3 is
  precisely the certificate that closes this: the ceiling \(\kappa_0\) gives
  finiteness on \(\lambda<1/\kappa_0\) with an explicit majorant, so §4B.3 is
  self-contained while §4B.5 depends on a finiteness check the implementation
  must perform rather than assume.
- The Taylor step behind \((\ast)\) uses \(\nabla\Phi(\gamma^\star)=0\), which is
  the defining property of \(\gamma^\star\) as EM fixed point (§2.2) — worth
  confirming it is the *global* minimizer and not merely stationary, though
  convexity of \(\Phi\) under (H2) should settle that. This is the one step in
  §4B.3 that numerics cannot check, since the check assumes closure.
- Whether the polynomial prefactor can be improved. §4B.9 shows it costs one to
  three orders of magnitude in probability, which is a real but bounded price.

### 4B.9 Numerical verification

`data-raw/_chk_tilt_bound.R` checks four things over a grid in
\((q,\kappa_0,d)\), and all four hold:

| claim | result |
|---|---|
| \(\lambda^\star=\kappa_0^{-1}-q/(2d)\) is the argmin | matches numeric argmin to \(2\times10^{-5}\) |
| closed form of §4B.3 equals \(\inf_\lambda\) | relative gap \(<4\times10^{-6}\) |
| side condition \(\Leftrightarrow d>q\kappa_0/(2(1-\kappa_0))\) | exact agreement, all 100 cases |
| \(M_Q(\lambda)=\prod_i(1-\lambda\kappa_i)^{-1/2}\) under closure | matches Monte Carlo |
| tilting identity reproduces §4A.3's exact tail | matches Monte Carlo |
| §4B.3 dominates the exact tail | 472/472 cases, flat and non-flat spectra |

Note the objective is *convex* in \(\lambda\) and diverges as
\(\lambda\uparrow1/\kappa_0\), so the interior stationary point is a minimum and
the side condition is exactly the condition for it to be interior. An
implementation that optimizes in the wrong direction gets no warning from the
formula, so this is worth an assertion in code.

**Conservatism.** The ratio of the §4B.3 bound to the exact closure tail has
median \(\approx1.3\times10^2\) over the grid, ranging from \(\approx10\) at
small \(q,d\) to \(\approx6\times10^5\) at \(q=8,\kappa_0=0.95,d=100\). In the
currency that matters this is an additive shift in \(d\) of
\(\log(\text{ratio})\cdot\kappa_0/(1-\kappa_0)\) — a few units of \(d\) for
well-conditioned problems, tens when \(\kappa_0\) approaches 1. So the bound is
usable but not sharp, and §4B.5 is worth having where the extra \(d\) propagates
into an expensive \(n\).

**Confirmed caveat.** Check C's agreement degrades from \(6\times10^{-3}\) at
\(d=1\) to \(1.3\times10^{-1}\) at \(d=6\) with \(\kappa_{\max}=0.7\) — exactly
the infinite-weight-variance regime \(\kappa>1/2\) flagged in §4B.5. That
degradation is a property of the Monte Carlo estimator, not of the identity, and
it is the concrete reason to prefer quadrature over \(Q\)-draws when
\(\kappa\) is large.

---

## 5. Stage 4 — the constant \(\varepsilon\) on the set

Given \(d\) and \(\varepsilon(\gamma^\star)\), Lemma 17(b) collapses to
arithmetic:

\[
\varepsilon=\varepsilon_d\,Q(\widetilde C_d),\qquad
\varepsilon_d=e^{-d}\varepsilon(\gamma^\star),\qquad
\varepsilon\ \ge\ \varepsilon(\gamma^\star)\,e^{-d}\Pr(\chi^2_q\le 2d).
\]

The \(\chi^2\) floor is what removes the need to compute a Gaussian measure of
a non-ellipsoidal convex set (§7). Everything here is elementary given the two
inputs; there is no reason for this stage to be anything but a closed-form
evaluation, and it should also report \(\varepsilon_d\) and \(Q(\widetilde C_d)\)
separately so the two sources of degradation are visible.

Where a uniform ceiling \(\kappa_0\) is available — which by §4B.4 is whenever
the population prior is proper — §4B.6 sharpens the floor to
\(\Pr(\chi^2_q\le2d/\kappa_0)\) at no additional cost, and this stage should use
that form.

---

## 6. Stage 5 — sweeps \(n\), and the \(d\) trade-off

Theorem 1 requires \((1-\varepsilon)^n+\delta\le\mathrm{tol}\), hence first
\(\delta<\mathrm{tol}\), and then

\[
n\ \ge\ \frac{\log(\mathrm{tol}-\delta)}{\log(1-\varepsilon)} .
\]

**\(d\) must be optimized, not supplied.** Raising \(d\) shrinks \(\delta\),
enlarging the numerator's budget, but shrinks \(\varepsilon\) like \(e^{-d}\),
hurting the denominator. So \(n(d)\) has a genuine interior minimum and the
implementation should locate it rather than pushing a non-obvious
one-dimensional optimization onto the user. This is structurally the same
trade-off `two_block_optimize_pilot_cost()` already solves for pilot length,
and should follow that precedent.

Two reporting points. The bound never drops below \(\delta\) however large
\(n\) is, so a requested `tol` at or below the achievable \(\delta\) must fail
loudly rather than return a large \(n\). And the bound is **uniform over
starting values in \(\widetilde C_d\)** — starting at \(\gamma^\star\) earns
nothing from Theorem 1 as stated, so no start-dependent discount should be
offered (see §8).

---

## 7. Three identities that make this cheap

Worth stating separately because they are what turn an apparently intractable
scheme into arithmetic.

**Complementarity (Lemma 13(3)).**
\(\Psi(\gamma)=\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}-\bar\Phi(\gamma)\).
Set membership therefore needs one scalar evaluation of the marginal negative
log-posterior — no inner convex program per point. The gradient comes free
too: \(\nabla\Psi(\gamma)=P_{11}(M(\gamma)-\gamma^\star)\), reusing the same
E-step as \(\gamma^\star\).

**\(\bar\Phi\) from exact gradients (Fisher's identity).** The gradient above
rearranges to \(\nabla\Phi(\gamma)=P_{11}(\gamma-M(\gamma))\), which is *exact*.
Integrating along the ray \(\gamma(t)=\gamma^\star+t(\gamma-\gamma^\star)\),

\[
\bar\Phi(\gamma)=\int_0^1(\gamma-\gamma^\star)^\top P_{11}
\bigl(\gamma(t)-M(\gamma(t))\bigr)\,dt ,
\]

with the integrand vanishing at \(t=0\) because \(M(\gamma^\star)=\gamma^\star\).
So \(\bar\Phi\) — and hence \(\Psi\) — is obtainable by a one-dimensional
quadrature over E-step evaluations, with no Laplace approximation of the
marginal.

Be precise about what this buys. It does not make \(\bar\Phi\) exact: the E-step
still needs \(\mathbb E[\beta_j\mid\gamma,y]\), which for GLMMs is itself a
quadrature. What it does is trade a Laplace approximation of a
\(Jp_{\mathrm{re}}\)-dimensional integral for \(J\) independent
\(p_{\mathrm{re}}\)-dimensional AGHQ problems, which is a far better-conditioned
approximation at controllable accuracy — not the removal of one. This revises the
pessimism of route 4(b) and §8(3): the approximation there is better than
"a Laplace evaluation per draw" suggests, though still an approximation.

**The \(\chi^2\) floor (Corollary 15).** \(Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)\),
because the \(Q\)-ellipsoid of squared radius \(2d\) is inscribed in
\(\widetilde C_d\) and \(\|\gamma-\gamma^\star\|^2_{P_{11}}\sim\chi^2_q\) under
\(Q=N(\gamma^\star,P_{11}^{-1})\). Replaces a Gaussian measure of a convex set
with a `pchisq` call.

---

## 8. Where the theory stops short of a number

To be carried into any user-facing documentation, since these determine which
computations may be labelled *certified*.

1. **Remark 6.1's concentration constant is untracked** (Appendix A.4), so the
   \(\rho\)-route of §4 is heuristic *as stated in the chapter*. §4B.3 closes
   this: the tilting bound carries the explicit factor
   \((2ed/q\kappa_0)^{q/2}\) and a sharper exponent, so nothing needs to rest on
   Remark 6.1. What remains open is only whether §4B survives review — see
   §4B.8.
1b. **Outside closure, \(d(\delta)\) needs a change of measure, not a change of
   set.** Corollary 15 and Lemma 14 reduce the *set* to an ellipsoid but leave
   the *measure* untouched, so neither certifies a number on its own. Route 4(d)
   (§4B) supplies the missing step by tilting to \(Q\); route 4(b) remains an
   estimate rather than a proof. The residual gap here is narrow: 4(d) is
   certified modulo the verification listed in §4B.8.
2. **\(\varepsilon(\gamma^\star)\) has no closed form outside Gaussian
   closure**, and the closure bound \((1+\kappa)^{-q/2}\) does not transfer.
   This is now the *only* stage of the pipeline with no general route, and is
   therefore the highest-value remaining piece of theory work.
3. **\(\bar\Phi\) for GLMMs is an approximation**, so the empirical route 4(b)
   is estimated rather than proved — though §7's third identity reduces this to
   per-group AGHQ rather than a full-marginal Laplace, which is a materially
   better approximation than the original phrasing implied. Route 4(d) is the
   certified escape hatch and, unlike 4(c), carries no family restriction.
4. **Theorem 1's proof is not in the vignette.** Its outline points to
   "Appendix, §A.0", but the appendix runs A.1–A.4. Worth resolving before
   this machinery is documented as user-facing.
5. **The certificate covers the restricted chain**, which is Gibbs for the
   truncated target (Lemma 24), *not* the production sampler. §7.2 of Chapter
   C05 ("What is and is not certified") is currently an empty heading; until
   it is written, any reported number should carry that caveat explicitly.

---

## 9. Computational responsibilities

Named units, not an API — no signatures are proposed and nothing here is
settled. Public-facing names require sign-off before anything lands in `R/`;
prototypes belong in `data-raw/` first.

**Model assembly.** Build \(P_{11}^{\mathrm{RE}}\), \(P_{11}\), and its
factorization from the hyper-designs, \(P_b\), and \((\mu_0,\Lambda_\gamma)\);
check (H3a). Produces \(\Sigma^\star\) with no further work.

**E-step.** Given \(\gamma\), return \(\{b_j\}\) and optionally \(\{V_j\}\),
with the exact / AGHQ / MC tiers of §2.4 behind one switch. Group-separable,
and the natural parallelization point.

**Mean map.** One solve against the pre-factored \(P_{11}\) given \(\{b_j\}\).
Trivial once the E-step exists.

**Fixed-point solver.** Accelerated EM to \(\gamma^\star\), returning
\(\gamma^\star\), \(\tilde J\) and its spectrum, \(\kappa\), \(\rho\),
\(\Sigma_\pi\), and \(\varepsilon(\gamma^\star)\) when closure holds. These
should be returned *together*: they all come from the terminal E-step, so
splitting them across separate entry points forces redundant computation.

**Refresh measure.** Package \(N(\gamma^\star,\Sigma^\star)\) — mostly
assembly once the solver has run.

**Profile evaluation.** \(\Psi(\gamma)\) via complementarity, with \(\bar\Phi\)
from the ray quadrature of §7, plus the membership test \(\Psi\le d\) and the
cheap sandwich variants of §4(c).

**Curvature ceiling.** \(\kappa_0\) from the pencil \((\Lambda_\gamma,P_{11})\),
or from curvature margins where the family supplies them (§4B.4). One
generalized eigenvalue, and the sole input to the default sizing route.

**Tilt MGF.** \(M_Q(\lambda)\) on a \(\lambda\)-grid by quadrature or \(Q\)-draws
(§4B.5), evaluating \(\Psi\) once per node and reusing it across \(\lambda\).
Optional — only needed when the closed form of §4B.3 is too conservative.

**Set sizing.** \(d\) from \(\delta\), with the closure / empirical / sandwich /
tilt routes of §4 as explicit alternatives that report which one was used.
Default to the tilt route (d) outside closure.

**Constant and budget.** \(\varepsilon\) from \(d\); the TV bound at given
\((d,n)\); and the sweep count for a tolerance with the \(d\)-optimization of
§6 folded in.

---

## 10. Validation plan

- **Stationarity residual.** \(\|P_{11}(\gamma^\star-M(\gamma^\star))\|\approx0\)
  as the solver's own convergence check — this is \(\nabla\Phi(\gamma^\star)\).
- **Gaussian reference.** In the LMM case every quantity has a closed form;
  the exact tier validates the AGHQ and MC tiers end to end.
- **Mean-versus-mode guard.** Compare \(\gamma^\star\) against the mode of a
  long pilot \(\gamma\)-chain. A mode-based E-step (§2.3) lands visibly off
  the marginal's peak, so this is the specific test that catches the trap.
- **EM rate check.** The observed linear convergence rate of the iteration
  should match \(\kappa=\lambda_{\max}(\tilde J)\); disagreement indicates an
  inconsistent \(V_j\).
- **Monotonicity.** \(\bar\Phi\) must decrease across EM iterations wherever it
  is evaluable; a violation localizes an E-step error immediately.
- **Certificate self-consistency.** \(\varepsilon\in(0,1]\),
  \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\), and the inner/outer
  sandwich sets must bracket the empirical \(\widetilde C_d\).
- **Tilting identity.** Under closure, \(M_Q(\lambda)\) evaluated by quadrature
  must match \(\prod_i(1-\lambda\kappa_i)^{-1/2}\) (§4B.2), and the §4B.3 closed
  form must dominate the exact weighted-\(\chi^2\) tail of §4A.3 at every \(d\)
  above the side condition. A crossing indicates an algebra error in \(\lambda^\star\).
- **Tilt ceiling.** \(\Psi(\gamma)\le\tfrac12\kappa_0\|\gamma-\gamma^\star\|^2_{P_{11}}\)
  should be spot-checked at sampled \(\gamma\); a violation falsifies either
  \(\kappa_0\) or the E-step's \(V_j\).
- **Simulation check.** For a small LMM, compare the certified \(n\) against an
  independently estimated mixing time; the certificate should be conservative
  but not absurd.

---

## 11. Open questions

1. Is the MC (MCEM) tier accurate enough for \(\varepsilon(\gamma^\star)\)'s
   convex program (§3), or does that need a separate deterministic route?
2. Does §4B survive review? Specifically the \(\lambda^\star\) algebra of §4B.3
   and the Rényi finiteness condition of §4B.8. If it does, §8(1) and §8(1b)
   both close and the \(\rho\)-route becomes redundant rather than heuristic.
3. Is there a usable \(G_j\) for Poisson-log — a restricted region, or a
   different argument — so route 4(c) covers it?
4. Should the certificate be recomputed as \(\gamma\) moves, or fixed once at
   \(\gamma^\star\)? The theory is static; a moving version would change what
   is being certified.
5. How should this interact with the existing preflight gates, given that
   (H3a) and (H3b) are already checked there?

---

## References

- `vignette("Chapter-C05")` — the theory this note implements.
- `vignette("Chapter-C03")` — `two_block_rate()` / `two_block_tv_bound()` /
  `two_block_l_for_tv()`, the precedent this pipeline should mirror.
- `notation.md`, `HIERARCHICAL_GENERALIZED_LINEAR_MODEL_NOTATION.md` — the
  \(H_j\), \(P_b\), \(\Lambda_\gamma\) conventions used above.
- `BLOCK_GIBBS_ERGODICITY_ING.md` §17 — the empirical escape-fraction pattern
  reused in §4(b).
- `PREFLIGHT_model_setup.md` — where (H3a)/(H3b) are already enforced.
- `data-raw/_chk_tilt_bound.R` — the numerical checks behind §4B.9.
- `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` — route 4(e); the escape bound via a
  scaled Gaussian envelope on \(\widetilde C_d^{\,c}\), which supersedes §4B
  wherever a \(P_{22}\) floor is available.
