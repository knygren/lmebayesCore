# Restricted Two-Block Gibbs: Minorization on $C_d$ and Distance to $\pi$

Bound on the **full** target $\pi$ when a two-block Gibbs sampler is **restricted**
to a subset $C_d$ of $(\gamma,\beta)$-space and satisfies Doeblin minorization on
$C_d$. Supports **Approach A** in `LOGIT_STATIC_TAIL_CERTIFICATION.md` §7.1.

**References.** Rosenthal (1995), Proposition 2 — `lmebayes/references/minor.pdf` (p. 5);
`MINORIZATION_GAUSSIAN_REFRESH.md`; `LOGIT_STATIC_TAIL_CERTIFICATION.md` §7.1.

**How this document is organized.** Part I proves the conditional bound (Proposition
R1): *given* a certified safe set $C_d$ with a minorization constant, the restricted
chain converges geometrically to $\pi$ up to a truncation error. Part II states the
existence result (Proposition R2) that supplies such a set for the hierarchical GLMM
posterior, together with the model hypotheses and the definitions it needs. Part II
states results only — every proof is relegated to a lettered appendix, one appendix
per existence ingredient, so that no single proof runs more than a few pages and each
can be read independently of the others. Appendix F assembles the ingredients of
Appendices A–E into the proof of Proposition R2 itself. A full index of every lemma,
corollary, and remark in the appendices, with a one-line statement and its location,
is given at the end of Part II.

*A note on cross-references carried over from the working draft.* Citations of the
form "§7.3.2A" or "Lemma 2A.3-G" below refer to the lemma **labels**, which are
unchanged from the draft this note was assembled from — only their **physical
location** has moved into the lettered appendices. A pointer from label to appendix
is in the Index of Results (§10). Renumbering the labels themselves (e.g. `2A.3-G` →
`B.7`) is a mechanical follow-up, not done here to avoid introducing transcription
errors into a proof this long; flag if you'd like that pass done as well.

---

# Part I — The Conditional Bound

## 1. Setup


**State space.** \((\gamma,\beta) \in \mathcal X\). Fix \(C_d \subseteq \mathcal X\) with
\(\pi(C_d) > 0\). Write \(\pi\) for the two-block Gibbs **target** on \(\mathcal X\).

**Conditional target.**

\[
\pi(A \mid C_d) := \frac{\pi(A \cap C_d)}{\pi(C_d)},
\qquad A \in \mathcal B(\mathcal X).
\]

**Restricted kernel.** \(P_1(x,\cdot)\) is one Gibbs sweep on \(\mathcal X\).
\(P_1(x,\cdot \mid C_d)\) is the Markov kernel of the sampler **restricted to**
\(C_d\) (defined once for the application — reject, renormalize, or confine so that
\(P_1(x,C_d \mid C_d)=1\) for \(x \in C_d\)). Let
\(P_n(x,\cdot \mid C_d) := (P_1(x,\cdot \mid C_d))^n\).

**Minorization (Doeblin on \(C_d\)).** There exist \(\varepsilon_d \in (0,1]\) and a
probability measure \(Q\) on \(C_d\) such that

\[
\boxed{
P_1(x,A \mid C_d) \;\ge\; \varepsilon_d\, Q(A)
\qquad
\forall\, x \in C_d,\ \forall\,\text{measurable } A \subseteq C_d.
}
\]

**Standing assumption (stationarity).** \(\pi(\cdot \mid C_d)\) is a stationary law
for \(P_1(\cdot \mid C_d)\) on \(C_d\). *(Must be verified for the chosen restriction;
not automatic from stationarity of \(\pi\) on \(\mathcal X\).)*

---

## 2. Distance conventions


**Total variation (standard).**

\[
\|\mu - \nu\|_{TV}
\;:=\;
\sup_{A \subseteq \mathcal X} \bigl|\mu(A) - \nu(A)\bigr|
\;=\;
\tfrac12 \int \bigl|d\mu - d\nu\bigr|.
\]

**L1 mass of the difference** (no factor \(\tfrac12\)):

\[
\|\mu - \nu\|_{1}
\;:=\;
\int \bigl|d\mu - d\nu\bigr|
\;=\;
2\,\|\mu - \nu\|_{TV}.
\]

---

## 3. Truncation lemma: $\pi(\cdot\mid C_d)$ vs $\pi$

**Lemma 1 (truncation).**


Write \(\mu := \pi(\cdot \mid C_d)\) and \(\nu := \pi\). On the partition
\(\mathcal X = C_d \cup C_d^c\):

| Piece | \(\mu\) | \(\nu\) | \(|\mu-\nu|\) on the piece |
|---|---|---|---|
| \(C_d^c\) | \(0\) | \(\pi(C_d^c)\) | \(\pi(C_d^c)\) |
| \(C_d\) | renormalized \(\pi\) | \(\pi\) | \(1 - \pi(C_d) = \pi(C_d^c)\) |

**Two effects (same size, different roles).**

1. **Complement:** \(\nu\) puts mass \(\pi(C_d^c)\) on \(C_d^c\); \(\mu\) puts none.
2. **Renormalization on \(C_d\):** \(\mu(A) = \pi(A)/\pi(C_d) > \pi(A)\) for
   \(A \subseteq C_d\), \(\pi(A)>0\); at \(A=C_d\) the gap is \(1-\pi(C_d)=\pi(C_d^c)\).

**Partition sum (L1).** The two pieces **add**:

\[
\boxed{
\|\pi(\cdot \mid C_d) - \pi\|_{1}
\;=\;
\pi(C_d^c) + \bigl(1-\pi(C_d)\bigr)
\;=\;
2\,\pi(C_d^c).
}
\]

**Sup over single events (TV).** One event \(A\) at a time; the maximum is **not**
attained at \(A = C_d \cup C_d^c\) (that gives \(0\)), but at \(A=C_d\) or
\(A=C_d^c\):

\[
\boxed{
\|\pi(\cdot \mid C_d) - \pi\|_{TV}
\;=\;
\pi(C_d^c).
}
\]

*(Same magnitude \(\pi(C_d^c)\) from either piece; \(\sup_A|\cdot|\) does not add
them.)*


## 4. Main theorem

**Theorem 1 (Proposition R1).**


**Proposition R1.** For \(x \in C_d\), \(n \ge 1\),

\[
\boxed{
\bigl\| P_n(x,\cdot \mid C_d) - \pi \bigr\|_{TV}
\;\le\;
(1-\varepsilon_d)^{n}
\;+\;
\pi(C_d^c).
}
\]

Equivalently, in L1 notation,

\[
\bigl\| P_n(x,\cdot \mid C_d) - \pi \bigr\|_{1}
\;\le\;
2(1-\varepsilon_d)^{n}
\;+\;
2\,\pi(C_d^c).
\]

*Proof.* Triangle inequality on \(\mathcal X\):

\[
\|P_n - \pi\|_{TV}
\;\le\;
\underbrace{\|P_n - \pi(\cdot \mid C_d)\|_{TV}}_{\text{Step A}}
\;+\;
\underbrace{\|\pi(\cdot \mid C_d) - \pi\|_{TV}}_{\text{Step B}}.
\]

**Step A.** Rosenthal (1995), Proposition 2 (`minor.pdf`, p. 5), applied on state
space \(C_d\) with kernel \(P_1(\cdot \mid C_d)\), minorization constant
\(\varepsilon_d\), and stationary law \(\pi(\cdot \mid C_d)\):

\[
\|P_n(x,\cdot \mid C_d) - \pi(\cdot \mid C_d)\|_{TV}
\;\le\;
(1-\varepsilon_d)^{n}.
\]

**Step B.** Section 3: \(\|\pi(\cdot \mid C_d) - \pi\|_{TV} = \pi(C_d^c)\). \(\square\)

*Nothing here is specific to the joint state space.* R1 holds for any \((\mathcal
X,\pi,C_d,P_1)\) with the stated minorization and stationarity, so §7.2 applies it
with \(\mathcal X=\mathbb R^q\), \(\pi=\pi_\gamma\) and \(P_1=q(\cdot\mid\cdot\,;C)\).

---

**Remark 1.1 (L1 form).** The same bound in L1, using
$\|\mu-\nu\|_1=2\|\mu-\nu\|_{TV}$ throughout (§2): the proof needs no change, only the
constant-2 dictionary of §2 applied to each term —
$\|P_n-\pi\|_1\le2(1-\varepsilon_d)^n+2\pi(C_d^c)$. *(Consolidates what were two
separately-proved copies of the same bound in the working draft — the second, "in
partition notation," added no new content.)*

**Remark 1.2 (nothing here is specific to the joint state space).** R1 holds for any
$(\mathcal X,\pi,C_d,P_1)$ with the stated minorization and stationarity; Part II
applies it with $\mathcal X=\mathbb R^q$, $\pi=\pi_\gamma$, $P_1=q(\cdot\mid\cdot\,;C)$.

---

# Part II — Existence of a Certified Safe Set (Proposition R2)

Part I is conditional: given $C_d$, a minorization constant, and stationarity of
$\pi(\cdot\mid C_d)$, Theorem 1 follows. This part states the model hypotheses under
which such a $C_d$ *provably exists* for the hierarchical GLMM posterior, defines the
objects the construction needs, and states the existence proposition. Every proof is
in Appendices A–F.


## 6. Model hypotheses


Work in the **two-block** hierarchy of `notation.md` / `model_setup()`:

\[
y_j \mid \beta_j \ \text{GLM with canonical/noncanonical link},\quad
\beta_j \mid \gamma \sim N(\mathcal W_j \gamma,\, \Psi),\quad
\gamma \sim N(\mu_0,\, \Lambda_\gamma^{-1}),
\]

with the package's **full-rank** and **estimability** diagnostics satisfied
(`check_identifiability()` / `groupef.estimable`, `popef.rank_ok` — see
`PREFLIGHT_model_setup.md`):

- **(H1) Proper posterior.** \(\pi(\gamma,\beta \mid y)\) is a probability measure on
  \(\mathcal X = \mathbb R^q \times \prod_j \mathbb R^{p_{\mathrm{re},j}}\) (no
  improper tail from rank deficiency). *Automatic when \(\Lambda_\gamma\succ0\) and the
  likelihood is bounded; in the flat limit \(\Lambda_\gamma=0\) it is **derived** from
  (H2)+(H3a)+(H3b)+`groupef.rank` — Lemma R2a‴, **Appendix A.3**.*
- **(H2) Log-concave conditionals.** For each block draw in two-block Gibbs,
  \(\pi(\beta_j \mid \gamma, y_j)\) is **log-concave** in \(\beta_j\) (logit, probit,
  Poisson, Gamma log, etc. — package GLM families with conjugate Gaussian RE prior).
- **(H3) Full rank / estimability.**
  - **(H3a)** Hyper-design full rank (`popef.rank_ok`): \(P_{11}^{\mathrm{RE}}\succ 0\).
  - **(H3b) Group-wise estimability:** for each group \(j\), a **finite maximum
    likelihood estimate** exists for the group GLM on \(\beta_j\) (given \(D_j\),
    link, and \(y_j\)). Equivalently: no complete separation, no all-zero Poisson
    intercept path, etc. — what `groupef.estimable` encodes via per-group `glm()`
    (`check_identifiability.R`, `R/check_identifiability.R`).
  - *With \(\Lambda_\gamma\succ0\) these are quality gates for the **existence of
    \(Q\)**: Lemmas 1A.1–1A.2 use neither. (H3a) is still needed for the
    \(\varepsilon\)-machinery of §7.3.2A, which rests on full support of \(m(\beta)\).
    **In the flat limit \(\Lambda_\gamma=0\) both become load-bearing** — (H3a) is
    what makes \(\Sigma^{\star(0)}=(\sum_j\mathcal W_j^\top\Psi^{-1}\mathcal W_j)^{-1}\)
    exist at all (Lemma 1B.1), and (H3b) is what makes the refresh centre
    \(\gamma^{\star(0)}\) exist (Lemma 1B.2). See **Remark 1B.4** (§7.3.1B) for how each
    failure propagates.*
- **(H4) Continuous positive kernel on compacts.** The one-step two-block kernel
  \(P_1(x,\cdot)\) has a **strictly positive continuous** density on compact
  subsets of \(\mathcal X\) (holds for Gaussian block updates with log-concave
  conditionals). ***No longer assumed.*** The Q-first construction of §7.2–§7.3
  writes the minorizing density down explicitly, so (H4) is a **consequence** of
  (H2)+(H3a) rather than a hypothesis. It is retained here only because the earlier
  C-first argument, and Approach B, still refer to it.

*(H1)–(H3) are the preflight gates under which the package runs production
samplers.*

## 7. Definitions

*(Reference section: nothing here is proved. The third column points to the appendix
result that establishes existence or well-posedness.)*


*(Reference section. Nothing here is proved; the third column points at the result
that establishes existence or well-posedness. Objects used only in §7 are defined
here anyway so that the definitions are not scattered across the proofs. Labels of
the form 1A.x, 1B.x, 2A.x, 2B.x are the lemmas of §7.3; the subsections below are
left unnumbered to avoid colliding with them.)*

Throughout, \(\gamma\in\mathbb R^q\) is the population parameter, \(\beta=(\beta_j)_{j\le J}\)
the random effects, and \(y\) the data. The hyper-design is written \(H_j\) or
\(\mathcal W_j\) interchangeably.

> **The star convention.** A superscript \(\star\) marks the *canonical* choice of an
> object, the one the Q-first construction is obliged to make: \(\gamma^\star\) is the
> centre of the refresh measure and \(\Sigma^\star=P_{11}^{-1}\) its covariance, with
> \(\lambda^\star\) the spectral cut. Do not read \(\gamma^\star\) as "some optimum" —
> it is pinned down three equivalent ways at once (mode, minimiser of \(\Phi\), fixed
> point of \(M\)), which is exactly what makes it computable. Note that \(\mu_0\)
> remains the **prior** mean and is a different object; \(\mu\) no longer denotes any
> parameter of the model, and survives only as a generic measure in §2–§3 and as the
> GLM mean function in Appendix A.

### Model primitives

These are **inputs**: they are read off `model_setup()` and require no existence
argument beyond the hypotheses of §7.1.

| Object | Definition | Existence / source |
|---|---|---|
| \(\ell_j(\beta_j)\) | group \(j\) GLM log-likelihood | model input |
| \(\mathcal W_j=H_j\), \(\Psi\), \(P_b:=\Psi^{-1}\) | hyper-design and RE covariance in \(\beta_j\mid\gamma\sim N(\mathcal W_j\gamma,\Psi)\) | model input |
| \(\Lambda_\gamma,\ \mu_0\) | population prior precision and mean, \(\gamma\sim N(\mu_0,\Lambda_\gamma^{-1})\) | model input |
| \(P_{11}^{\mathrm{RE}}\) | \(\sum_{j=1}^J H_j^\top P_b H_j\succeq0\) | **(H3a)** \(\iff P_{11}^{\mathrm{RE}}\succ0\) |
| \(P_{11}\) | \(\Lambda_\gamma+P_{11}^{\mathrm{RE}}\) — the **exact** conditional precision of \(\gamma\) given \(\beta\) | Lemma 1A.1 (\(\succ0\)) |
| \(m(\beta)\) | \(P_{11}^{-1}\bigl(\Lambda_\gamma\mu_0+\sum_j H_j^\top P_b\beta_j\bigr)\) — conditional mean of \(\gamma\) given \(\beta\) | needs \(P_{11}\succ0\) |
| \(\pi,\ \pi_\gamma\) | joint posterior on \(\mathcal X\); its \(\gamma\)-marginal | **(H1)**; in the flat limit, Lemma R2a‴ (App. A.3) |
| \(q(\gamma'\mid\gamma)\) | \(\displaystyle\int\phi_q\bigl(\gamma';m(\beta),P_{11}^{-1}\bigr)\pi(\beta\mid\gamma,y)\,d\beta\) — one-sweep \(\gamma\)-transition density | §7.2; Markov because \(\gamma\)-marginal of two-block Gibbs is Markov |

### The mean map, its curvature, and the rate

| Object | Definition | Existence / source |
|---|---|---|
| \(b_j(\gamma)\) | \(E[\beta_j\mid\gamma,y]\) | **(H1)**, **(H2)** |
| \(V_j(\gamma)\) | \(\mathrm{Cov}(\beta_j\mid\gamma,y)\); \(V_j\preceq P_b^{-1}\) | Brascamp–Lieb under **(H2)** — Lemma R2d(b) |
| \(M(\gamma)\) | \(E[\gamma'\mid\gamma]=P_{11}^{-1}\bigl(\Lambda_\gamma\mu_0+\sum_j H_j^\top P_b\,b_j(\gamma)\bigr)\) | Lemma R2d |
| \(J(\gamma)\), \(\tilde J(\gamma)\) | \(J=\partial M/\partial\gamma=P_{11}^{-1}\sum_j H_j^\top P_bV_jP_bH_j\); \(\tilde J=P_{11}^{-1/2}(\sum_j H_j^\top P_bV_jP_bH_j)P_{11}^{-1/2}\) symmetric psd | Lemma R2d(a) (exponential tilt) |
| \(\kappa(\gamma)\) | \(\lambda_{\max}\bigl(\tilde J(\gamma)\bigr)\in[0,1]\) — pointwise contraction rate | Lemma R2d(b): \(\kappa\le1\) |
| \(E_{\lambda^\star}\), \(\lambda^\star\) | \(\{\gamma:\kappa(\gamma)\le\lambda^\star\}\), closed; \(\lambda^\star\in(0,1)\) chosen so \(\pi_\gamma(E_{\lambda^\star}^c)<\delta_1\) | **Cor. R2d-4**; \(E_{\lambda^\star}=\mathbb R^q\) and \(\delta_1=0\) when \(\Lambda_\gamma\succ0\) |

### The posterior potential, its minimiser, and the refresh measure \(Q\)

| Object | Definition | Existence / source |
|---|---|---|
| \(A_j(\theta)\) | \(\log\int\exp\{\ell_j(\beta)-\tfrac12\beta^\top P_b\beta+\theta^\top\beta\}\,d\beta\) — log-partition of the tilted group conditional | finite under **(H1)**+**(H2)**; **Cor. R2d-5** |
| \(\Phi(\gamma)\) | \(\tfrac12(\gamma-\mu_0)^\top\Lambda_\gamma(\gamma-\mu_0)+\tfrac12\gamma^\top P_{11}^{\mathrm{RE}}\gamma-\sum_jA_j(P_bH_j\gamma)\), so \(\log\pi(\gamma\mid y)=-\Phi(\gamma)+\mathrm{const}\) | **Cor. R2d-5(1)**; convex by **R2d-5(3)** |
| \(\gamma^\star\) | \(\arg\min_\gamma\Phi(\gamma)\) \(=\) **mode** of \(\pi(\gamma\mid y)\) \(=\) **fixed point of the kernel's mean map**, \(M(\gamma^\star)=\gamma^\star\) | **Lemma 1A.2** (existence/uniqueness); **Lemma R2d** (contraction); **Cor. R2d-5(2)** (the three descriptions agree) |
| \(\Sigma^\star\) | \(P_{11}^{-1}\) — hypothesis **(Q1)** | **Lemma 1A.1** |
| \(Q\), \(q_Q\) | \(Q:=N(\gamma^\star,\Sigma^\star)\), density \(q_Q=\phi_q(\cdot;\gamma^\star,\Sigma^\star)\) — the **refresh measure**, fixed before any set | **Cor. 1A.3** |

### The certified set and the minorization constant

| Object | Definition | Existence / source |
|---|---|---|
| \(g(\gamma'\mid\gamma)\) | \(\log q(\gamma'\mid\gamma)-\log q_Q(\gamma')\), minimised over the **destination** \(\gamma'\) | finite and jointly continuous; strictly convex in \(\gamma'\) by **Lemma 2A.1** |
| \(\gamma'_\star(\gamma)\) | \(\arg\min_{\gamma'}g(\gamma'\mid\gamma)\) — the **unique worst destination** from \(\gamma\); the \(\gamma'\) at which the kernel is hardest to minorize | exists (**2A.2**), unique (**2A.1**), and \(\lvert\gamma'_\star(\gamma)\rvert\le R_K\) uniformly for \(\gamma\) in a compact \(K\) (**2A.3**) |
| \(\mathcal D(\gamma)\) | \(g\bigl(\gamma'_\star(\gamma)\mid\gamma\bigr)=\min_{\gamma'}g(\gamma'\mid\gamma)\in(-\infty,0]\) — the **deficiency** | \(C^1\) by the envelope identity (**Remark 2A.3-E**); \(C^2\) and **concave**, \(\nabla^2\mathcal D=-P_{11}J\), by **Lemma 2A.3-F** |
| \(\varepsilon(\gamma)\) | \(e^{\mathcal D(\gamma)}=\inf_{\gamma'}q(\gamma'\mid\gamma)/q_Q(\gamma')\in(0,1]\) — the **profile**, a function of the **source** \(\gamma\); carries no \(d\) | attained and \(>0\) by **Lemma 2A.2**; unique minimiser by **2A.1**; continuous by **2A.3**; log-concave with peak at \(\gamma^\star\) by **2A.3-F** |
| \(\bar\Phi(\gamma)\) | \(\Phi(\gamma)-\Phi(\gamma^\star)\ge0\), \(=0\) only at \(\gamma^\star\) — the **HPD potential** | normalisation is **essential**: \(\Phi\) is defined only up to the unknown marginal-likelihood constant, \(\bar\Phi\) is not |
| \(\Psi(\gamma)\) | \(\mathcal D(\gamma^\star)-\mathcal D(\gamma)=\log\frac{\varepsilon(\gamma^\star)}{\varepsilon(\gamma)}\ge0\) — the **deficiency gap**, the function that indexes the certified family | convex, \(\nabla^2\Psi=P_{11}J\succeq0\); **complementarity** \(\bar\Phi+\Psi=\frac12\|\gamma-\gamma^\star\|^2_{P_{11}}\) — **Lemma 2A.3-F(3)**; **coercive**, \(S_\flat\preceq\nabla^2\Psi\preceq P_{11}^{\mathrm{RE}}\) — **Lemma 2A.3-G** |
| \(S_\flat\) | \(\sum_jH_j^\top P_b(P_b+G_j)^{-1}P_bH_j\succ0\), where \(\mathcal G_j(\beta_j)\preceq G_j\) is a curvature ceiling — the **circumscribing metric** | \(\succ0\) under **(H3a)**; Cramér–Rao floor on \(V_j(\gamma)\), **Lemma 2A.3-G**. **Prior-free**; needs bounded GLM weights, else use the radius \(1+d/c\) |
| \(d>0\) | radius parameter of the certified family, in **nats of deficiency** (not a distance) | chosen in §7.2 Step 4 via **Lemma R2a** |
| \(\varepsilon_d\) | \(e^{-d}\varepsilon(\gamma^\star)\in(0,1)\) — the **attained minimum** \(\min_{\widetilde C_d}\varepsilon\), reached on \(\partial\widetilde C_d\) | **Lemma 2A.4**, by continuity (2A.3) + compactness (2A.3-G); positivity is \(\varepsilon(\gamma^\star)>0\) (2A.2). **No closed form** for that scalar outside Gaussian closure — Lemma R2b′ |
| \(\widetilde C_d\) | \(\{\gamma:\varepsilon(\gamma)\ge\varepsilon_d\}=\{\gamma:\Psi(\gamma)\le d\}\) — the **largest** set carrying the constant \(\varepsilon_d\) | closed by **2A.3**, convex by **2A.3-F**, **compact** by **2A.3-G** with \(\widetilde C_d\subseteq\{\|\gamma-\gamma^\star\|^2_{S_\flat}\le2d\}\); \(\pi_\gamma(\widetilde C_d^{\,c})\to0\) by **Lemma R2a** |
| \(C\) | \(\widetilde C_d\cap E_{\lambda^\star}\) — the **certified set**; **compact**, and convex when the spectral truncation is vacuous | §7.2 Step 5; \(\pi_\gamma(C^c)<\delta_1+\delta_2\le\delta\) |
| \(Q_C\) | \(Q(\cdot\mid C)\) — \(Q\) **truncated** to \(C\); a truncated normal, not normal | \(Q(C)>0\); \(Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)\) by the inscribed-ellipsoid bound of §7.3.2 |
| \(\varepsilon\) | \(\varepsilon_d\,Q(C)\in(0,1]\) — the **Doeblin constant of the restricted chain** | **Lemma 2A.5(b)**; the factor \(Q(C)\) is the price of confining the sampler |
| \(\delta,\delta_1,\delta_2\) | total mass budget and its split, \(\delta_1+\delta_2\le\delta\): \(\delta_1\) for the spectral set, \(\delta_2\) for the deficiency tail | §7.2 Step 4 |

### Flat-limit (\(\Lambda_\gamma\downarrow0\)) counterparts

Write \(\Lambda_\gamma=\tau\bar\Lambda\) with \(\bar\Lambda\succeq0\) and \(\mu_0\)
fixed, \(\tau\downarrow0\). Each object below is the \(\tau=0\) reading of its
counterpart in the two preceding tables; **(H3a)** and **(H3b)** replace what the
prior used to supply.

| Object | Definition | Existence / source |
|---|---|---|
| \(\Phi_\tau,\ \Phi_0\) | \(\Phi_\tau(\gamma)=\Phi_0(\gamma)+\tfrac{\tau}{2}(\gamma-\mu_0)^\top\bar\Lambda(\gamma-\mu_0)\); \(\Phi_0\) is \(\Phi\) with the prior term deleted | quadratic increment because the population prior is Gaussian |
| \(\Sigma^{\star(0)}\) | \((P_{11}^{\mathrm{RE}})^{-1}=\bigl(\sum_j\mathcal W_j^\top\Psi^{-1}\mathcal W_j\bigr)^{-1}\) | **Lemma 1B.1** — exists **iff (H3a)** |
| \(\gamma^{\star(0)}\) | \(\arg\min\Phi_0\) = mode of the flat-limit \(\gamma\)-marginal | **Lemma 1B.2** — needs **(H3a)+(H3b)** for coercivity; \(\gamma^\star(\tau)\to\gamma^{\star(0)}\) |
| \(Q_0\) | \(N\bigl(\gamma^{\star(0)},\Sigma^{\star(0)}\bigr)\); \(Q_\tau\to Q_0\) in total variation | **Cor. 1B.3** |
| \(\varepsilon(\gamma\mid\Lambda_\gamma)\), \(\varepsilon(\gamma\mid0)\) | profile at prior scale \(\tau\), and its limit | **Lemma 2B.1** (well defined, finite, uniform on compacts), **2B.2** (continuous) |
| \(\Psi_\tau,\ \Psi_0\) | deficiency gap at prior scale \(\tau\), and its limit; \(\Psi_0=\frac12\|\cdot-\gamma^{\star(0)}\|^2_{P_{11}^{\mathrm{RE}}}-\bar\Phi_0\) | **Lemma 2B.3** — convex and coercive, with the **\(\tau\)-free** curvature sandwich \(S_\flat\preceq\nabla^2\Psi_\tau\preceq P_{11}^{\mathrm{RE}}\); \(\Psi_\tau\to\Psi_0\) uniformly on compacts |
| \(\widetilde C_d(\Lambda_\gamma)\), \(\widetilde C_d(0)\) | \(\{\varepsilon(\cdot\mid\tau)\ge\varepsilon_d(\tau)\}=\{\Psi_\tau\le d\}\), and its limit | **Lemma 2B.3** — convex, **compact**, all inside one \(\tau\)-free compact \(K_d\), Hausdorff-convergent. Compactness needs **(H3a)** only |
| \(\varepsilon_d(0)\) | \(e^{-d}\,\varepsilon\bigl(\gamma^{\star(0)}\mid0\bigr)>0\) | **Lemma 2B.4**; also \(\varepsilon_d(\Lambda_\gamma)\to\varepsilon_d(0)\), so the constant does not collapse |

### Three distinctions worth keeping straight

1. **The set and the constant come from the same function.** \(\widetilde C_d\) is a
   **superlevel set of \(\varepsilon\) itself**, equivalently a sublevel set of the
   convex deficiency gap \(\Psi\). The constant is then the **attained minimum**
   \(\min_{\widetilde C_d}\varepsilon=\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\),
   produced exactly as one would expect — continuity plus compactness, i.e. Weierstrass
   (Lemma 2A.4) — but with a **known value**, because the set is indexed by the very
   function being minimised, so the minimum sits on the level that names the set. The
   pair is a fixed point of the correspondence between sets and constants:
   \(\min_{\widetilde C_d}\varepsilon=\varepsilon_d\) and
   \(\{\varepsilon\ge\varepsilon_d\}=\widetilde C_d\). Indexing by \(\bar\Phi\) instead
   — the HPD regions — is mass-optimal but breaks both identities, leaving a minimum
   with no evaluable value and no visible \(d\)-dependence. Both families are convex
   and compact; they differ in the *price* of compactness, which here is **(H3a)** and
   prior-free (Lemma 2A.3-G, from nondegeneracy of the random-effect conditional
   covariances) and there is coercivity of \(\bar\Phi\), costing **(H3b)** in the flat
   limit. The two are related by the complementarity identity
   \(\bar\Phi+\Psi=\frac12\|\gamma-\gamma^\star\|^2_{P_{11}}\), which also shows why
   the curvatures split the way they do — \(\nabla^2\Psi\) prior-free,
   \(\nabla^2\bar\Phi=P_{11}-\nabla^2\Psi\) carrying \(\Lambda_\gamma\). Compared in
   Remark R2b-0.
2. **\(g\) is minimised over \(\gamma'\); \(\varepsilon\) is a function of \(\gamma\).**
   The destination variable is integrated out of the picture by the infimum, leaving a
   function of the source state only. Both appear in the same displays, so the
   distinction is easy to lose.
3. **\(\varepsilon_d\) and \(\varepsilon\) are different constants.**
   \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\) minorizes the
   **unrestricted** chain against the **untruncated** \(Q\) (Lemma 2A.5(a));
   \(\varepsilon=\varepsilon_d\,Q(C)\) minorizes the **restricted** chain against the
   **truncated** \(Q_C\) (Lemma 2A.5(b)). Only the second appears in Proposition R2.

---

## 8. Existence proposition

**Proposition R2 (existence of a certified safe set for the \(\gamma\)-chain).**
Under **(H1)–(H3)**, for every \(\delta > 0\) there exist:

1. a **compact** set \(C \subseteq \mathbb R^q\) — explicitly a **superlevel set of the
   minorization profile itself**,
   \(\widetilde C_d=\{\gamma:\varepsilon(\gamma)\ge\varepsilon_d\}\), which is convex
   and compact, intersected with the closed spectral set of Cor. R2d-4,
2. a **restricted** \(\gamma\)-chain kernel \(q(\cdot \mid \cdot\,;C)\) supported on
   \(C\),
3. a constant \(\varepsilon \in (0,1]\) and a probability measure \(Q_C\) on \(C\)
   which is the **restriction to \(C\) of a multivariate normal** — explicitly,
   \(Q_C=Q(\cdot\mid C)\) with
   \(Q=N\bigl(\gamma^\star,\Sigma^\star\bigr)\), \(\Sigma^\star=P_{11}^{-1}\) and \(\gamma^\star\) the
   fixed point of the mean map — equivalently, the **mode of the \(\gamma\)-marginal
   posterior** \(\pi(\gamma\mid y)\) (**Lemma 1A.2**, **Cor. R2d-5**). So \(Q_C\) is a
   **truncated multivariate normal** centred at the posterior mode of the population
   parameters, with the exact conditional precision \(P_{11}\) as its precision,

such that

\[
\pi_\gamma(C^c) < \delta,
\qquad
\pi_\gamma(\cdot \mid C)\ \text{is stationary for}\ q(\cdot\mid\cdot\,;C),
\qquad
q(\gamma,A \mid C) \ge \varepsilon\, Q_C(A)\ \ \forall \gamma \in C,\ A \subseteq C,
\]

and hence **Proposition R1** gives, for \(\gamma \in C\),

\[
\|q_n(\gamma,\cdot \mid C) - \pi_\gamma\|_{TV}
\;\le\;
(1-\varepsilon)^{n} + \pi_\gamma(C^c)
\;<\;
(1-\varepsilon)^{n} + \delta.
\]

## 9. Certification roadmap

Proposition R2 needs three existence ingredients and one assembly step. This table
is the single place that tracks what each depends on and where it is proved; it
replaces three overlapping lists in the working draft (an "open items" list, a
"Status" table, and "What R2 does not claim").

| Ingredient | What is needed | Appendix | Status |
|---|---|---|---|
| Refresh measure $Q$ | $\Sigma^\star=P_{11}^{-1}\succ0$ (Lemma 1A.1 / 1B.1) and centre $\gamma^\star$, the fixed point of the mean map / mode of $\pi(\gamma\mid y)$ (Lemma 1A.2 / 1B.2, mechanism Lemma R2d) | **A** | Proved under (H2); proper prior needs no rank condition, flat limit needs (H3a)+(H3b) |
| Minorization constant $\varepsilon_d$ and certified set $\widetilde C_d$ | Profile $\varepsilon(\gamma)>0$, unique minimiser, continuous, $\log\varepsilon$ concave (Lemma 2A.3-F) and coercive (Lemma 2A.3-G) $\Rightarrow$ $\widetilde C_d=\{\varepsilon\ge\varepsilon_d\}$ compact, $\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)$ attained (Lemma 2A.4) | **B** | Proved for every GLM family under (H3a); no closed-form value for the scalar $\varepsilon(\gamma^\star)$ outside Gaussian closure |
| Gaussian-closure specialization | Explicit $\mathcal D$, ellipsoidal $\widetilde C_d$, constant $\ge(1+\kappa)^{-q/2}e^{-d}$ | **C** | Closed form; LMM only |
| Spectral cutoff $\lambda^\star$ | $\exists\lambda^\star\in(0,1)$ with $\pi(\kappa>\lambda^\star)<\delta_1$; attained, continuous | **E** | Proved (Lemma R2a′, Cor. R2a″); vacuous when $\Lambda_\gamma\succ0$ |
| Exhaustion + restricted-kernel stationarity | $\pi_\gamma(\widetilde C_d^{\,c})\downarrow0$ (Lemma R2a); $\pi(\cdot\mid C)$ stationary for the restricted kernel (Lemma R2c) | **D** | Proved |
| Assembly | Combine A + B (or C) + D + E into Proposition R2; lift $\gamma$-chain bound to the joint chain (Remark R2-J) | **F** | Proved |
| Scope / symmetric-case constants | What is and isn't certified; sharper constants under S1–S4 | **G** | Reference only, no new proofs |

**Outstanding beyond Proposition R2** (tracked in Appendix G, not restated here):
a numerical value for $\varepsilon(\gamma^\star)$ outside Gaussian closure; a
computable exponent $\rho$ in $\varepsilon_d\gtrsim\varepsilon(\gamma^\star)\delta_2^\rho$;
a margin for $\kappa$ bounded away from 1 in the flat limit; and the relationship
between the *restricted* kernel R2 certifies and the *unrestricted* kernel the
package actually runs (Approach B territory, not this note).


## 10. Index of results (Appendices A–F)

Every named lemma, corollary, and remark used in the appendices, in one table.
Statement text is abbreviated to a phrase; full statements and proofs are at the
cited location.

| Label | One-line content | Location |
|---|---|---|
| Lemma R2d | Mean map $M$ has a unique fixed point $\gamma^\star$ (contraction) | A.0 |
| Remark R2d-1 | Computability: each iterate needs only $J$ group-conditional solves | A.0 |
| Cor. R2d-2 | Identification with the rate matrix / Laplace analogue | A.0 |
| Cor. R2d-3 | When $\gamma^\star$ is the posterior mean | A.0 |
| Cor. R2d-4 | Spectral tail set $E_{\lambda^\star}$ for the $\gamma$-chain | A.0 |
| Cor. R2d-5 | $\gamma^\star$ is the marginal posterior mode (variational form) | A.0 |
| Lemma 1A.1 | Covariance of $Q$ exists, $\Lambda_\gamma\succ0$ (no rank condition) | A.1 |
| Lemma 1A.2 | Centre $\gamma^\star$ of $Q$ exists, $\Lambda_\gamma\succ0$ | A.1 |
| Cor. 1A.3 | Refresh density $Q$ exists, $\Lambda_\gamma\succ0$ | A.1 |
| Lemma 1B.1 | Limiting covariance $\Sigma^{\star(0)}$ exists iff (H3a) | A.2 |
| Lemma 1B.2 | Limiting centre $\gamma^{\star(0)}$ finite under (H3a)+(H3b) | A.2 |
| Cor. 1B.3 | Limiting refresh density $Q_0$ exists | A.2 |
| Remark 1B.4 | (H3a)/(H3b) are sharp and fail differently | A.2 |
| Lemma 2A.1 | Strict convexity: minimiser of $g(\cdot\mid\gamma)$ unique | B.1 |
| Lemma 2A.2 | Minimum exists for every $\gamma$; $\varepsilon(\gamma)>0$ | B.1 |
| Lemma 2A.3 | Profile $\varepsilon$ continuous | B.1 |
| Remark 2A.3-E | Envelope identity for the slope of $\mathcal D$ | B.1 |
| Lemma 2A.3-F | Deficiency $\mathcal D$ concave; closed form for both derivatives; complementarity identity | B.1 |
| Lemma 2A.3-G | Deficiency gap $\Psi$ coercive; $\widetilde C_d$ compact between concentric ellipsoids | B.1 |
| Lemma 2A.4 | $\varepsilon_d=\min_{\widetilde C_d}\varepsilon$ is the attained minimum | B.1 |
| Lemma 2A.5 | Minorization condition holds on $\widetilde C_d$ | B.1 |
| Lemma 2B.1 | Limiting profile well defined and finite at $\Lambda_\gamma=0$ | B.2 |
| Lemma 2B.2 | Limiting profile continuous | B.2 |
| Lemma 2B.3 | Limiting set convex/compact, uniform in the prior | B.2 |
| Lemma 2B.4 | Limiting constant; continuity in the prior scale $\tau$ | B.2 |
| Lemma 2B.5 | Minorization condition holds in the flat limit | B.2 |
| Remark 2B.6 | What 2B adds beyond pointwise existence at $\tau=0$ | B.2 |
| Lemma R2b$'$ | Closed form for $\mathcal D$ and $\varepsilon_d$ under Gaussian closure | C |
| Remark R2b-0 | What's available in general: existence, not a formula | C |
| Remark R2b-1 | Precision gap is an exact bridge covariance (identity, not hypothesis) | C |
| Remark R2b-2 | LMM level sets are ellipsoids about $\gamma^\star$ | C |
| Remark R2b-3 | GLMM kernel is a mixture, not Gaussian — what's lost | C |
| Remark R2b-4 | Certified constants under closure | C |
| Remark R2b-5 | Division of labour between A and B subsections | C |
| Lemma R2a | Exhaustion: $\pi_\gamma(\widetilde C_d^{\,c})\downarrow0$ | D |
| Remark R2a-1 | What changed from the earlier (HPD) construction | D |
| Remark R2a-2 | The Pareto tension between $\delta$ and $\varepsilon_d$, quantified | D |
| Lemma R2c | Restricted kernel well defined; $\pi(\cdot\mid C)$ stationary | D |
| Remark R2c-1 | Why this definition of the restricted kernel | D |
| Remark R2c-2 | Relation to the sampler actually run by the package | D |
| Remark R2c-3 | The $\gamma$-chain instance of Lemma R2c | D |
| Lemma R2a$'$ | $A^{(0)}$ well defined; $0\le T\le1$; $\exists\lambda^\star$ with $\pi(T>\lambda^\star)<\delta_1$ | E.1 |
| Cor. R2a$''$ | $T$ continuous, attained eigenvalue; limits $\Lambda_\gamma\downarrow0,\ \delta_1\downarrow0$ | E.2 |
| Lemma R2a$'''$ | Flat-limit posterior $\pi_0$ proper (supplies (H1) at $\Lambda_\gamma=0$) | E.3 |
| A.4 Example | Separation degrades $\lambda^\star$; what actually fails | E.4 |
| Remark R2-J | Lifting the $\gamma$-chain bound to the joint $(\gamma,\beta)$ chain | F |

---

# Appendices

Each appendix proves one existence ingredient of Proposition R2 and is readable on
its own given Part II's definitions (§7) and hypotheses (§6). Appendix F assembles
A–E into the proof of Proposition R2 itself.

## Appendix A — Existence of the refresh measure $Q$

### A.0 Machinery: the mean map, the marginal posterior, and the rate


##### Lemma R2d (the mean map has a unique fixed point)

**Statement.** Assume **(H2)**, **(H3a)**, and either \(\Lambda_\gamma\succ0\) or the
strict-curvature condition (c) below. Then the mean map

\[
M(\gamma)
\;=\;
P_{11}^{-1}\Bigl(\Lambda_\gamma\mu_0 + \sum_j H_j^\top P_b\, b_j(\gamma)\Bigr),
\qquad
b_j(\gamma)=E[\beta_j\mid\gamma,y],
\]

is a **contraction** in the norm \(\|x\|_{P_{11}}=\|P_{11}^{1/2}x\|\) with modulus
\(\kappa<1\); it has a **unique** fixed point \(\gamma^\star\), and the iteration
\(\gamma^{(k+1)}=M(\gamma^{(k)})\) converges to it geometrically at rate \(\kappa\) from
any start.

*Proof.*

**(a) Jacobian via exponential tilt.** The RE conditional is

\[
\pi(\beta_j\mid\gamma,y)
\;\propto\;
\exp\Bigl\{\ell_j(\beta_j) + (P_b H_j\gamma)^\top\beta_j
- \tfrac12\beta_j^\top P_b\beta_j\Bigr\},
\]

an exponential family in the natural parameter \(\theta_j=P_bH_j\gamma\) with
sufficient statistic \(\beta_j\). Hence \(\partial b_j/\partial\theta_j
= V_j(\gamma):=\mathrm{Cov}(\beta_j\mid\gamma,y)\), and by the chain rule
\(\partial b_j/\partial\gamma = V_j P_b H_j\). Therefore

\[
\boxed{\;
J(\gamma)
\;:=\;
\frac{\partial M}{\partial\gamma}
\;=\;
P_{11}^{-1}\sum_j H_j^\top P_b\, V_j(\gamma)\, P_b H_j
\;=\;
P_{11}^{-1/2}\,\tilde J(\gamma)\,P_{11}^{1/2},
\;}
\qquad
\tilde J := P_{11}^{-1/2}\Bigl(\sum_j H_j^\top P_b V_j P_b H_j\Bigr)P_{11}^{-1/2},
\]

with \(\tilde J\) **symmetric positive semidefinite**; so \(J\) has real
nonnegative spectrum, and \(J\) is \(\succ0\) under **(H3a)** whenever \(V_j\succ0\).

**(b) Brascamp–Lieb upper bound.** Under **(H2)** each \(\ell_j\) is concave, so
\(-\nabla^2\log\pi(\beta_j\mid\gamma,y) = -\nabla^2\ell_j + P_b \succeq P_b\): the
conditional is \(P_b\)-strongly log-concave. Brascamp–Lieb (applied to linear
functionals) gives \(V_j\preceq P_b^{-1}\), hence

\[
\sum_j H_j^\top P_b V_j P_b H_j
\;\preceq\;
\sum_j H_j^\top P_b H_j
\;=\;
P_{11}^{\mathrm{RE}}
\qquad\Longrightarrow\qquad
\tilde J
\;\preceq\;
I - P_{11}^{-1/2}\Lambda_\gamma P_{11}^{-1/2}.
\]

**(c) Strictness.** If \(\Lambda_\gamma\succ0\), the display gives
\(\kappa:=\sup_\gamma\lambda_{\max}(\tilde J(\gamma))
\le 1-\lambda_{\min}(P_{11}^{-1/2}\Lambda_\gamma P_{11}^{-1/2})<1\)
**unconditionally** — the population prior alone supplies the gap. In the flat
limit \(\Lambda_\gamma=0\) one needs \(V_j\prec P_b^{-1}\) strictly, i.e. the group
**likelihood** must contribute curvature: if \(-\nabla^2\ell_j\succeq\epsilon_j P_b\)
uniformly then \(V_j\preceq\bigl((1+\epsilon_j)P_b\bigr)^{-1}\) and

\[
\kappa\;\le\;\frac{1}{1+\min_j\epsilon_j}\;<\;1 ,
\]

the same relative group information \(\epsilon_j\) as **Appendix A.1–A.2**, and the
same prior-free ceiling as **Cor. R2a″(b2)**.

**(d) Contraction and Banach.** For any \(\gamma_1,\gamma_2\),
\(M(\gamma_1)-M(\gamma_2)=\int_0^1 J(\gamma_t)(\gamma_1-\gamma_2)\,dt\) with
\(\gamma_t\) the segment, so
\(P_{11}^{1/2}(M(\gamma_1)-M(\gamma_2))
=\int_0^1\tilde J(\gamma_t)\,P_{11}^{1/2}(\gamma_1-\gamma_2)\,dt\) and

\[
\|M(\gamma_1)-M(\gamma_2)\|_{P_{11}}\;\le\;\kappa\,\|\gamma_1-\gamma_2\|_{P_{11}} .
\]

\(\mathbb R^q\) with \(\|\cdot\|_{P_{11}}\) is complete, so Banach's fixed-point
theorem gives existence, uniqueness, and geometric convergence of
\(\gamma^{(k+1)}=M(\gamma^{(k)})\). \(\square\)

**Remark R2d-1 (computability).** Each iterate needs only \(J\) group-conditional
means \(b_j(\gamma^{(k)})\) — the same quantity the sampler's Block-1 update already
computes. So \(\gamma^\star\) is **computable offline** at geometric cost, unlike
\(E[\gamma\mid y]\), which is not available in closed form.

**Corollary R2d-2 (identification with the rate matrix).** With
\(S(\beta)=\sum_j H_j^\top P_b B_j(\beta_j)^{-1} P_b H_j\) and
\(B_j = P_b + G_j(\beta_j)\) the Laplace precision of Appendix A.0, the rate matrix
is \(A = P_{11}^{-1/2}S\,P_{11}^{-1/2}\). Comparing with \(\tilde J\): the two
differ **only** by replacing the exact conditional covariance \(V_j(\gamma)\) with
its Laplace approximation \(B_j^{-1}\). Hence

> "\(M\) is a contraction" and "\(\lambda_{\max}(A)<1\)" are the **same statement**,
> modulo Laplace-vs-exact covariance.

In particular \(\kappa\) and \(T=\lambda_{\max}(A^{(0)})\) are controlled by the same
\(\epsilon_j\), and Appendix A.1's \(\lambda^\star\) doubles as a bound on \(\kappa\)
on the certified set.

**Corollary R2d-3 (when \(\gamma^\star\) is the posterior mean).** Let \(\pi_\gamma\) be the
\(\gamma\)-marginal, stationary for \(q\). Then
\(E_{\pi_\gamma}[\gamma]=E_{\pi_\gamma}\bigl[E[\gamma'\mid\gamma]\bigr]
=E_{\pi_\gamma}[M(\gamma)]\).

- If \(M\) is **affine** — the LMM / Gaussian-closure case, where \(b_j(\gamma)\) is
  affine in \(\gamma\) — then \(E_{\pi_\gamma}[M(\gamma)]=M(E_{\pi_\gamma}[\gamma])\),
  so \(E[\gamma\mid y]\) is a fixed point, and by uniqueness
  \(\gamma^\star=E[\gamma\mid y]\). **Hypothesis 2 of `MINORIZATION_GAUSSIAN_REFRESH.md` is
  therefore derived, not assumed.**
- For non-Gaussian likelihoods \(M\) is nonlinear and the two differ by a Jensen
  gap \(E_{\pi_\gamma}[M(\gamma)]-M(E_{\pi_\gamma}[\gamma])\). The **fixed point**,
  not the posterior mean, is the correct centre (Centering remark below).

**Corollary R2d-4 (spectral set for the \(\gamma\)-chain).** Define the
**exact-covariance rate**

\[
\kappa(\gamma)\;:=\;\lambda_{\max}\bigl(\tilde J(\gamma)\bigr),
\qquad
E_{\lambda^\star}\;:=\;\{\gamma:\ \kappa(\gamma)\le\lambda^\star\}\subseteq\mathbb R^q ,
\]

so that the contraction modulus of Lemma R2d is \(\kappa=\sup_\gamma\kappa(\gamma)\).
Then:

- **(a) Proper population prior.** If \(\Lambda_\gamma\succ0\), Lemma R2d(c) gives
  \(\kappa(\gamma)\le1-\lambda_{\min}(P_{11}^{-1/2}\Lambda_\gamma P_{11}^{-1/2})<1\)
  **uniformly in \(\gamma\)**. So \(E_{\lambda^\star}=\mathbb R^q\) for that
  \(\lambda^\star\): **no spectral truncation is needed**, and \(\delta_1=0\).
- **(b) Flat limit.** If \(\Lambda_\gamma=0\), then \(\kappa\le1\) everywhere
  (Lemma R2d(b)) and \(\kappa(\gamma)<1\) at every finite \(\gamma\) whenever each
  \(\mathcal G_j\succ0\) (group full rank, positive GLM weights). The argument of
  **Lemma R2a′ step (iv)** (Appendix A.1) applies **verbatim** with \(T\) replaced by
  \(\kappa\) and \(\pi\) by \(\pi_\gamma\): the decreasing events \(\{\kappa>t\}\)
  shrink to \(\varnothing\), so by **(H1)** \(\pi_\gamma(\kappa>t)\to0\) as
  \(t\uparrow1\) and for every \(\delta_1>0\) there is \(\lambda^\star\in(0,1)\) with
  \(\pi_\gamma(E_{\lambda^\star}^c)<\delta_1\). \(E_{\lambda^\star}\) is closed since
  \(\kappa\) is continuous.
- **(c) Relation to Appendix A.** \(\kappa\) is the exact-covariance analogue of
  \(T=\lambda_{\max}(A^{(0)})\) (Cor. R2d-2): the two coincide when the Laplace
  approximation \(B_j^{-1}\) is exact (Gaussian, and asymptotically in group size),
  and both obey the prior-free ceiling \(1/(1+\min_j\epsilon_j)\).

**Corollary R2d-5 (variational form: \(\gamma^\star\) is the marginal posterior *mode*).** For
each group let

\[
A_j(\theta)\;:=\;\log\int
\exp\bigl\{\ell_j(\beta)-\tfrac12\beta^\top P_b\beta+\theta^\top\beta\bigr\}\,d\beta
\]

be the log-partition function of the tilted group conditional, so that
\(\nabla A_j(P_bH_j\gamma)=b_j(\gamma)\) and \(\nabla^2A_j(P_bH_j\gamma)=V_j(\gamma)\)
(this is the exponential-family structure of part (a)). Define

\[
\boxed{\;
\Phi(\gamma)\;:=\;\tfrac12(\gamma-\mu_0)^\top\Lambda_\gamma(\gamma-\mu_0)
+\tfrac12\gamma^\top P_{11}^{\mathrm{RE}}\gamma
-\sum_j A_j\bigl(P_bH_j\gamma\bigr).\;}
\]

Then:

1. **\(\Phi\) is the negative log marginal posterior of \(\gamma\).** Completing the
   square inside \(\int e^{\ell_j(\beta_j)}\phi(\beta_j;H_j\gamma,\Psi)\,d\beta_j
   =\text{const}\cdot e^{-\frac12\gamma^\top H_j^\top P_bH_j\gamma+A_j(P_bH_j\gamma)}\)
   and multiplying over \(j\) against the \(\gamma\)-prior gives
   \(\log\pi(\gamma\mid y)=-\Phi(\gamma)+\text{const}\).
2. **Stationarity is the fixed-point equation.**
   \(\nabla\Phi(\gamma)=P_{11}\bigl(\gamma-M(\gamma)\bigr)\), so — since
   \(P_{11}\succ0\) —
   \[
   \gamma^\star=M(\gamma^\star)\iff\nabla\Phi(\gamma^\star)=0 .
   \]
3. **Curvature is the contraction modulus.**
   \(\nabla^2\Phi(\gamma)=P_{11}\bigl(I-J(\gamma)\bigr)
   =P_{11}^{1/2}\bigl(I-\tilde J(\gamma)\bigr)P_{11}^{1/2}\), hence
   \[
   \nabla^2\Phi(\gamma)\succeq0\iff\kappa(\gamma)\le1,
   \qquad
   \nabla^2\Phi(\gamma)\succ0\iff\kappa(\gamma)<1 .
   \]
   So Lemma R2d(b) — Brascamp–Lieb, \(\tilde J\preceq I\) — is exactly the statement
   that the \(\gamma\)-**marginal posterior is log-concave**, and "\(M\) is a
   contraction near \(\gamma\)" is exactly "\(\pi(\gamma\mid y)\) is strictly
   log-concave at \(\gamma\)".

Consequently \(\gamma^\star=\arg\max_\gamma\pi(\gamma\mid y)\) whenever the maximiser exists.
This sharpens Cor. R2d-3: the fixed point is always the **mode**, and Cor. R2d-3 is
the observation that mode \(=\) mean under Gaussian closure. It also gives a second
route to \(\gamma^\star\) — minimise \(\Phi\) — which is available when the Banach iteration is
not (Lemma 1B.2).

---


### A.1 Proper population prior ($\Lambda_\gamma\succ0$)


The refresh measure is \(Q=N(\gamma^\star,\Sigma^\star)\) with \(\Sigma^\star\) fixed by
**(Q1)** and \(\gamma^\star\) by **(Q2)**. For \(Q\) to be a genuine nondegenerate multivariate
normal density two things must exist: a **positive definite covariance**, and a
**finite centre**. Subsection 1A does this for a proper population prior, 1B in the
flat limit. Throughout, \(P_{11}=\Lambda_\gamma+P_{11}^{\mathrm{RE}}\) with
\(P_{11}^{\mathrm{RE}}=\sum_j\mathcal W_j^\top\Psi^{-1}\mathcal W_j\succeq0\).


**Lemma 1A.1 (the covariance of \(Q\) exists).** Assume \(\Lambda_\gamma\succ0\). Then

\[
P_{11}=\Lambda_\gamma+P_{11}^{\mathrm{RE}}\;\succeq\;\Lambda_\gamma\;\succ\;0,
\qquad\text{so}\qquad
\Sigma^\star:=P_{11}^{-1}\ \text{exists and }\ 0\prec\Sigma^\star\preceq\Lambda_\gamma^{-1}.
\]

*Proof.* \(P_{11}^{\mathrm{RE}}\) is a sum of terms
\(\mathcal W_j^\top\Psi^{-1}\mathcal W_j\succeq0\), so \(P_{11}\succeq\Lambda_\gamma\succ0\);
inverting reverses the order. \(\square\)

> **No rank condition is used.** Lemma 1A.1 does **not** need **(H3a)**: a proper
> population prior supplies invertibility on its own, and the hyper-design may be
> rank-deficient. This is the first of the two places where 1B will have to work
> harder.

**Lemma 1A.2 (the centre of \(Q\) exists: the marginal posterior mode).** Assume
**(H2)** and \(\Lambda_\gamma\succ0\). Then the negative log marginal posterior
\(\Phi\) of Cor. R2d-5 is **\(\Lambda_\gamma\)-strongly convex**:

\[
\boxed{\;\nabla^2\Phi(\gamma)\;=\;P_{11}-\sum_j\mathcal W_j^\top P_bV_j(\gamma)P_b\mathcal W_j
\;\succeq\;P_{11}-P_{11}^{\mathrm{RE}}\;=\;\Lambda_\gamma\;\succ\;0
\quad\text{for every }\gamma.\;}
\]

Consequently \(\Phi\) has a **unique finite minimiser**

\[
\gamma^\star\;=\;\arg\min_{\gamma\in\mathbb R^q}\Phi(\gamma)
\;=\;\arg\max_{\gamma\in\mathbb R^q}\pi(\gamma\mid y),
\qquad
\gamma^\star=M(\gamma^\star),
\]

i.e. the mode of the \(\gamma\)-marginal posterior exists, is unique, is finite, and is
the fixed point of the mean map.

*Proof.* The Hessian formula is Cor. R2d-5(3). Brascamp–Lieb under **(H2)**
(Lemma R2d(b)) gives \(V_j\preceq P_b^{-1}\), hence
\(\sum_j\mathcal W_j^\top P_bV_jP_b\mathcal W_j\preceq\sum_j\mathcal W_j^\top P_b\mathcal W_j
=P_{11}^{\mathrm{RE}}\), which is the displayed bound. A \(\Lambda_\gamma\)-strongly
convex function satisfies
\(\Phi(\gamma)\ge\Phi(\gamma_0)+\nabla\Phi(\gamma_0)^\top(\gamma-\gamma_0)
+\tfrac12\|\gamma-\gamma_0\|^2_{\Lambda_\gamma}\), hence is **coercive**, so its
minimum is attained; strict convexity makes the minimiser unique and finite.
Stationarity is the fixed-point equation by Cor. R2d-5(2). \(\square\)

> **Two independent routes, and why the variational one is recorded.** Lemma R2d(c)
> also produces \(\gamma^\star\), by Banach: \(\Lambda_\gamma\succ0\) forces
> \(\kappa\le1-\lambda_{\min}(P_{11}^{-1/2}\Lambda_\gamma P_{11}^{-1/2})<1\), so \(M\)
> is a contraction and \(\gamma^{(k+1)}=M(\gamma^{(k)})\) converges geometrically. The two
> routes agree by Cor. R2d-5(2). The contraction route is the better **algorithm**;
> the variational route is the better **existence proof**, because it is the one that
> survives 1B, where \(\kappa\) may equal \(1\). Note also that Lemma 1A.2, like 1A.1,
> uses no rank condition.

**Corollary 1A.3 (the refresh density \(Q\) exists).** Under **(H2)** and
\(\Lambda_\gamma\succ0\),

\[
Q:=N\bigl(\gamma^\star,\Sigma^\star\bigr),
\qquad
q_Q(\gamma')=(2\pi)^{-q/2}\det(P_{11})^{1/2}
\exp\bigl\{-\tfrac12(\gamma'-\gamma^\star)^\top P_{11}(\gamma'-\gamma^\star)\bigr\},
\]

is a **nondegenerate multivariate normal probability density**, strictly positive and
continuous on all of \(\mathbb R^q\). It satisfies **(Q1)** — precision dominance
\(\Lambda_q\preceq\Lambda_Q=P_{11}\), which is the identity of Remark R2b-1 and holds
for every GLM family — and **(Q2)** by Lemma 1A.2.

*Proof.* Combine Lemmas 1A.1 and 1A.2; positivity and continuity are immediate from
the displayed formula, which is well posed because \(\det P_{11}>0\). \(\square\)


### A.2 Flat-prior limit ($\Lambda_\gamma\downarrow0$)


Write \(\Lambda_\gamma=\tau\bar\Lambda\) with \(\bar\Lambda\succeq0\) fixed and
\(\mu_0\) fixed, and let \(\tau\downarrow0\). Objects carrying the prior are written
\(P_{11}(\tau),\Sigma^\star(\tau),\gamma^\star(\tau),Q_\tau\). Both lemmas of 1A degrade at
\(\tau=0\): the prior no longer supplies invertibility, and it no longer supplies
strong convexity. **(H3a)** replaces the first, **(H3b)** the second.

**Lemma 1B.1 (the limiting covariance exists iff (H3a)).** As \(\tau\downarrow0\),
\(P_{11}(\tau)=\tau\bar\Lambda+P_{11}^{\mathrm{RE}}\downarrow P_{11}^{\mathrm{RE}}\)
monotonically in the psd order. Then:

- **Under (H3a)** (\(P_{11}^{\mathrm{RE}}\succ0\), i.e. the stacked hyper-design
  \([\mathcal W_1;\dots;\mathcal W_J]\) has rank \(q\)),
  \[
  \boxed{\;
  \Sigma^\star(\tau)\;\uparrow\;
  \Sigma^{\star(0)}:=\bigl(P_{11}^{\mathrm{RE}}\bigr)^{-1}
  =\Bigl(\sum_j\mathcal W_j^\top\Psi^{-1}\mathcal W_j\Bigr)^{-1}\;\prec\;\infty ,
  \;}
  \]
  the convergence being monotone in the psd order and entrywise.
- **Without (H3a)** the limit does not exist: for \(0\ne v\in\ker P_{11}^{\mathrm{RE}}\),
  \(v^\top\Sigma^\star(\tau)v\ge(v^\top v)^2/(v^\top P_{11}(\tau)v)
  =(v^\top v)^2/(\tau\,v^\top\bar\Lambda v)\to\infty\).

**Interpretation.** \(\Sigma^{\star(0)}\) is the **generalised least squares**
covariance for the hyper-regression of the random effects on the hyper-design with
weight \(\Psi^{-1}\), and correspondingly the limiting mean function
\(m_0(\beta)=\Sigma^{\star(0)}\sum_j\mathcal W_j^\top\Psi^{-1}\beta_j\) is the GLS
estimate of \(\gamma\) from the \(\beta_j\). The prior terms \(\Lambda_\gamma\mu_0\)
drop out of \(m\) and \(M\) entirely. So the flat-limit refresh measure is the sampling
law of the GLS fit — and it is the **widest** admissible refresh measure, since
\(\Sigma^\star(\tau)\) increases to it.

*Proof.* Monotonicity of \(P_{11}(\tau)\) is immediate and inversion reverses it. Under
(H3a) the limit \(P_{11}^{\mathrm{RE}}\) is nonsingular, and matrix inversion is
continuous at nonsingular matrices, giving \(\Sigma^\star(\tau)\to\Sigma^{\star(0)}\).
The divergence claim is Cauchy–Schwarz:
\((v^\top v)^2=(v^\top P_{11}^{1/2}P_{11}^{-1/2}v)^2\le
(v^\top P_{11}v)(v^\top P_{11}^{-1}v)\). \(\square\)

**Lemma 1B.2 (the limiting centre is a finite vector, under (H3a)+(H3b)).** Assume
**(H2)**, **(H3a)**, **(H3b)**. Let

\[
\Phi_0(\gamma)=\tfrac12\gamma^\top P_{11}^{\mathrm{RE}}\gamma
-\sum_jA_j\bigl(P_b\mathcal W_j\gamma\bigr)
\;=\;\sum_j\mathcal D_j\bigl(P_b\mathcal W_j\gamma\bigr)+\text{const},
\qquad P_b=\Psi^{-1},
\]

be the flat-limit negative log marginal posterior, where the **group information
deficit**

\[
\mathcal D_j(\theta):=\tfrac12\theta^\top P_b^{-1}\theta+a_j-A_j(\theta)
\quad\text{satisfies}\quad
\mathcal D_j\ge0,
\qquad
\nabla^2\mathcal D_j(\theta)=P_b^{-1}-V_j\succeq0
\]

(nonnegativity from \(\ell_j\le\ell_j^{\max}\) with
\(a_j:=\ell_j^{\max}+\tfrac12\log\det(2\pi P_b^{-1})\); convexity from Brascamp–Lieb).
Then \(\Phi_0\) is convex and **coercive**, so

\[
\gamma^{\star(0)}:=\arg\min_\gamma\Phi_0(\gamma)=\arg\max_\gamma\pi_0(\gamma\mid y)
\]

exists, is **unique**, and is a **finite** vector. Moreover the modes converge:
\(\gamma^\star(\tau)\to\gamma^{\star(0)}\) as \(\tau\downarrow0\).

*Proof.* **Coercivity.** Under **(H2)+(H3a)+(H3b)** Lemma R2a‴ (**Appendix A.3**) makes
the flat-limit joint posterior proper, hence \(\pi_0(\gamma\mid y)\) proper; its step 4
supplies an exponential envelope \(e^{-c\|x-x_0\|}\), so \(\Phi_0\) grows at least
linearly and is coercive. A convex coercive function attains its minimum.
**Uniqueness.** \(\nabla^2\Phi_0=P_{11}^{\mathrm{RE}}\bigl(I-J_0\bigr)\succ0\) wherever
\(\mathcal G_j\succ0\) (Cor. R2d-4(b)), so \(\Phi_0\) is strictly convex.
**Convergence of the modes.** \(\Phi_\tau(\gamma)=\Phi_0(\gamma)
+\tfrac{\tau}{2}(\gamma-\mu_0)^\top\bar\Lambda(\gamma-\mu_0)\ge\Phi_0(\gamma)\).
Minimality gives
\(\Phi_0(\gamma^\star(\tau))\le\Phi_\tau(\gamma^\star(\tau))\le\Phi_\tau(\gamma^{\star(0)})
=\Phi_0(\gamma^{\star(0)})+C\tau\) with
\(C=\tfrac12(\gamma^{\star(0)}-\mu_0)^\top\bar\Lambda(\gamma^{\star(0)}-\mu_0)\). So every \(\gamma^\star(\tau)\)
lies in the sublevel set \(\{\Phi_0\le\Phi_0(\gamma^{\star(0)})+C\tau\}\), which is bounded by
coercivity; any limit point \(\bar\gamma\) satisfies
\(\Phi_0(\bar\gamma)\le\Phi_0(\gamma^{\star(0)})\) by continuity, hence \(\bar\gamma=\gamma^{\star(0)}\) by
uniqueness. A bounded family with a single limit point converges. \(\square\)

> **Banach is unavailable here, and that is the point.** At \(\tau=0\) Lemma R2d(c)
> leaves only \(\kappa\le1\), and for families whose GLM weights vanish in the tails
> (logit, probit at large \(|\eta|\)) one has \(\sup_\gamma\kappa(\gamma)=1\) exactly,
> so the contraction argument yields nothing globally. Lemma 1B.2 is precisely the
> variational replacement. *If* one additionally has uniform group curvature
> \(-\nabla^2\ell_j\succeq\epsilon_jP_b\), then Lemma R2d(c) gives the **prior-free
> ceiling** \(\kappa_0\le1/(1+\min_j\epsilon_j)<1\) and the Banach iteration is
> restored with a rate that does not degrade as \(\tau\downarrow0\).

**Corollary 1B.3 (the limiting refresh density exists).** Under **(H2)**, **(H3a)**,
**(H3b)**,

\[
Q_0:=N\bigl(\gamma^{\star(0)},\Sigma^{\star(0)}\bigr)
\]

is a nondegenerate multivariate normal, and \(Q_\tau\to Q_0\) in **total variation** as
\(\tau\downarrow0\).

*Proof.* Nondegeneracy is Lemma 1B.1 plus finiteness of \(\gamma^{\star(0)}\) (Lemma 1B.2).
Since \(\gamma^\star(\tau)\to\gamma^{\star(0)}\) and \(\Sigma^\star(\tau)\to\Sigma^{\star(0)}\succ0\), the
Gaussian densities converge pointwise; both are probability densities, so Scheffé gives
\(L^1\), i.e. TV, convergence. \(\square\)

**Remark 1B.4 (the two hypotheses are sharp, and they fail differently).**

| Failure | What breaks | Which lemma dies |
|---|---|---|
| **(H3a)** stacked \(\mathcal W\) rank-deficient | \(P_{11}^{\mathrm{RE}}\) singular ⟹ \(\Sigma^{\star(0)}\) does not exist | **1B.1** — and with it \(q_{Q_0}\), so §7.3.2B has no density to compare against |
| **(H3b)** some group not estimable | \(\Phi_0\) not coercive ⟹ \(\pi_0\) improper | **1B.2** — no finite mode, so no **(Q2)**; and no probability measure for Lemma R2a |

Both are visible in \(\Phi_0\). If \(0\ne v\in\bigcap_j\ker\mathcal W_j\) then
\(\Phi_0(\gamma+tv)=\Phi_0(\gamma)\) for all \(t\) — the flat-limit marginal posterior
is **constant** along \(v\), hence improper, with no mode; the term
\(\tau\bar\Lambda\) was the only thing keeping \(P_{11}\) invertible. If instead some
group is separated, \(\ell_j\) fails to decay along a direction, \(A_j\) saturates its
Gaussian ceiling there, \(\mathcal D_j\) stops growing, and \(\Phi_0\) flattens at
infinity — **Appendix A.4** is exactly this example. The degenerate no-data case is the
extreme: \(\ell_j\equiv\text{const}\) gives
\(A_j(\theta)=\text{const}+\tfrac12\theta^\top P_b^{-1}\theta\), so
\(\mathcal D_j\equiv0\) and \(\Phi_0\equiv\text{const}\), as it must be — a flat prior
with no information is improper.

For every \(\tau>0\) both failures are **repaired** by \(\Lambda_\gamma\succ0\)
(Lemmas 1A.1–1A.2 used neither hypothesis), at the price of constants that degrade as
\(\tau\downarrow0\). This is the same trade-off as Cor. R2a″(b3), and it is why the
preflight gates `popef.rank_ok` and `groupef.estimable` are **necessary**, not merely
prudent, whenever the package runs with a flat population prior.

---


## Appendix B — Existence of the minorization constant $\varepsilon_d$ and the certified set $\widetilde C_d$


With \(Q\) in hand from Section 1, fix the **log-ratio**

\[
\boxed{\;
g(\gamma'\mid\gamma)\;:=\;\log q(\gamma'\mid\gamma)-\log q_Q(\gamma'),
\qquad \gamma,\gamma'\in\mathbb R^q,
\;}
\]

where \(q(\cdot\mid\gamma)\) is the one-sweep \(\gamma\)-transition density of §7.2 and
\(q_Q=q_{N(\gamma^\star,\Sigma^\star)}\). Both are strictly positive and jointly continuous, so
\(g\) is finite and jointly continuous. Define the **minorization profile**

\[
\varepsilon(\gamma)\;:=\;\exp\Bigl\{\inf_{\gamma'\in\mathbb R^q}g(\gamma'\mid\gamma)\Bigr\}
\;=\;\inf_{\gamma'\in\mathbb R^q}\frac{q(\gamma'\mid\gamma)}{q_Q(\gamma')} .
\]

*(The profile does not depend on \(d\); the subscript in \(\varepsilon_d\) belongs to
the **constant** that indexes the certified family, not to the function.)*

**The certified family: superlevel sets of the profile itself.** The certified set is
the set on which the minorization is asked to hold, so take it to be **exactly** that
set and nothing else. Measure the profile against its own maximum value, which by
Lemma 2A.3-F(2) is attained at the refresh centre \(\gamma^\star\), and define the
**deficiency gap**

\[
\Psi(\gamma)\;:=\;\mathcal D(\gamma^\star)-\mathcal D(\gamma)
\;=\;\log\frac{\varepsilon(\gamma^\star)}{\varepsilon(\gamma)}\;\ge\;0,
\qquad
\boxed{\;
\widetilde C_d\;:=\;\bigl\{\gamma\in\mathbb R^q:\ \varepsilon(\gamma)\ge\varepsilon_d\bigr\}
\;=\;\bigl\{\gamma:\ \Psi(\gamma)\le d\bigr\},
\quad
\varepsilon_d:=e^{-d}\varepsilon(\gamma^\star),
\quad d>0 .
\;}
\]

The one thing this buys, and it is the whole point, is that the Weierstrass minimum
**has a value**. The minorization constant is still produced the standard way — the
minimum of \(\varepsilon\) over the certified set, existing and attained because
\(\varepsilon\) is continuous (Lemma 2A.3) and \(\widetilde C_d\) is compact
(Lemma 2A.3-G) — but because the set is a superlevel set *of the very function being
minimised*, the minimum lands on the level that names the set:

\[
\min_{\widetilde C_d}\varepsilon\;=\;\varepsilon_d\;=\;e^{-d}\varepsilon(\gamma^\star),
\qquad\text{attained on }\partial\widetilde C_d=\{\Psi=d\}
\]

(Lemma 2A.4). So the whole \(d\)-dependence of the certificate is explicit, the constant
is **sharp** rather than merely valid, and the only quantity left to evaluate is the
single number \(\varepsilon(\gamma^\star)\), the profile at the mode. Index the family
by \(\bar\Phi\) instead and the same Weierstrass argument runs, but its output is a
number with no known relation to \(\varepsilon(\gamma^\star)\) or to \(d\).

> **The complementarity identity.** \(\Psi\) is not a new function. By Lemma 2A.3-F(3),
> \[
> \boxed{\;
> \Psi(\gamma)\;=\;\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}-\bar\Phi(\gamma),
> \qquad\text{i.e.}\qquad
> \bar\Phi(\gamma)+\Psi(\gamma)\;=\;\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}}
> \;=\;-\log\frac{q_Q(\gamma)}{q_Q(\gamma^\star)} .
> \;}
> \]
> The HPD potential and the deficiency gap **partition the refresh quadratic**.
> Exponentiating, this is the exact density identity
> \[
> \pi(\gamma\mid y)\,\varepsilon(\gamma)\;=\;c\,q_Q(\gamma),
> \qquad
> c:=\frac{\pi(\gamma^\star\mid y)\,\varepsilon(\gamma^\star)}{q_Q(\gamma^\star)}\in(0,1],
> \]
> valid for **every** GLM family. (That \(c\le1\): integrating
> \(\pi=c\,q_Q/\varepsilon\) and using \(\varepsilon\le1\) gives
> \(1=c\int q_Q/\varepsilon\ge c\). Since \(\varepsilon\le1\) the identity also says
> \(\pi_\gamma\) is everywhere at least \(c\,q_Q\) — the posterior is the heavier of
> the two, as it must be for a Gaussian refresh to minorize the kernel at all.)
> So the certified set is an **importance-ratio** sublevel set,
> \[
> \widetilde C_d=\Bigl\{\gamma:\ \frac{\pi(\gamma\mid y)}{q_Q(\gamma)}
> \ \le\ e^{d}\,\frac{\pi(\gamma^\star\mid y)}{q_Q(\gamma^\star)}\Bigr\}
> \]
> — the region where the posterior is no more than \(e^{d}\) times as heavy as the
> refresh Gaussian, relative to that ratio at the mode. That is precisely the right
> object to certify a minorization *against \(q_Q\)*, and it is why the HPD region was
> never quite the natural set here: \(\bar\Phi\) measures the posterior against a flat
> reference, \(\Psi\) measures it against the reference actually being used.
> Membership is testable from \(\pi(\cdot\mid y)\) up to its normalising constant.

Every property needed downstream follows, and each is one line:

| Property | Reason |
|---|---|
| **\(\Psi\ge0\)**, so \(\varepsilon_d\le\varepsilon(\gamma^\star)\le1\) | \(\gamma^\star\) maximises \(\mathcal D\) — Lemma 2A.3-F(2). Equivalently \(\bar\Phi\le\frac12\|\cdot\|^2_{P_{11}}\) |
| **the constant holds on the set** | by definition of a superlevel set |
| **and is the attained minimum**, \(\min_{\widetilde C_d}\varepsilon=\varepsilon_d\), on \(\partial\widetilde C_d\) | continuity (2A.3) + compactness (2A.3-G): Weierstrass gives attainment, coercivity + IVT gives that the level \(d\) is reached. So the constant is **sharp** — Lemma 2A.4 |
| **convex** | \(\mathcal D\) is concave — Lemma 2A.3-F(1); equivalently \(\nabla^2\Psi=P_{11}J\succeq0\) |
| **closed**, and \(\gamma^\star\) is an interior point | \(\varepsilon\) is continuous (Lemma 2A.3) and \(\Psi(\gamma^\star)=0<d\) |
| **bounded**, hence **compact** | \(\nabla^2\Psi\succ0\) because the RE conditional covariances \(V_j(\gamma)\) are nondegenerate — Lemma 2A.3-G, which needs only **(H3a)** and is **prior-free** |
| **nested**, \(\bigcup_{d>0}\widetilde C_d=\mathbb R^q\) | \(\varepsilon>0\) pointwise (Lemma 2A.2), i.e. \(\Psi<\infty\) everywhere |
| **carries all the mass in the limit** | continuity from above of \(\pi_\gamma\) — Lemma R2a |
| **maximal** | \(\widetilde C_d\) is the **largest** set on which \(\varepsilon\ge\varepsilon_d\) holds; any other certified set at that constant is a subset of it |

> **Both mass bounds come from the same scalar.** The inequality \(\bar\Phi\ge0\) —
> which is nothing but \(\gamma^\star\) minimising \(\Phi\) — turns the identity into
> \[
> 0\;\le\;\Psi(\gamma)\;\le\;\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}} ,
> \]
> so the \(Q\)-ellipsoid of squared radius \(2d\) is **inscribed** in
> \(\widetilde C_d\) and, since
> \(\|\gamma-\gamma^\star\|^2_{P_{11}}\sim\chi^2_q\) under \(Q\),
> \[
> \boxed{\;
> Q\bigl(\widetilde C_d\bigr)\;\ge\;\Pr\bigl(\chi^2_q\le 2d\bigr),
> \qquad
> \pi_\gamma\bigl(\widetilde C_d^{\,c}\bigr)\;\le\;
> \pi_\gamma\bigl(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d\bigr) .
> \;}
> \]
> Both hold for **every** GLM family with no closure assumption, and both are tails of
> the **same** scalar statistic \(\|\gamma-\gamma^\star\|^2_{P_{11}}\) — under \(Q\) for
> the refresh mass, under \(\pi_\gamma\) for the posterior mass. Note the second bound
> is not available for the HPD family, where the corresponding inequality
> \(\bar\Phi\ge\frac12\|\cdot\|^2_{\Lambda_\gamma}\) degrades with the prior and is
> vacuous at \(\Lambda_\gamma=0\); here the metric is \(P_{11}\succeq P_{11}^{\mathrm{RE}}\),
> which survives the flat limit under **(H3a)**. Remark R2a-2 turns this into a rate.

> **And the set is trapped between two concentric ellipsoids.** The upper bound above
> is one half of Lemma 2A.3-G; the other half is that
> \(\nabla^2\Psi=\sum_jH_j^\top P_bV_j(\gamma)P_bH_j\) is bounded **away from zero**,
> because the conditional covariances \(V_j(\gamma)=\mathrm{Cov}(\beta_j\mid\gamma,y)\)
> of the random-effect block are nondegenerate — a positive log-concave density on
> \(\mathbb R^{p_{\mathrm{re}}}\) cannot have singular covariance. With
> \(S_\flat:=\sum_jH_j^\top P_b(P_b+G_j)^{-1}P_bH_j\succ0\) from the Cramér–Rao bound
> (\(G_j\) any ceiling on the group curvature \(\mathcal G_j(\beta_j)\);
> \(G_j=\frac14Z_j^\top N_jZ_j\) for binomial links),
> \[
> \boxed{\;
> \bigl\{\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\le2d\bigr\}
> \;\subseteq\;\widetilde C_d\;\subseteq\;
> \bigl\{\|\gamma-\gamma^\star\|^2_{S_\flat}\le2d\bigr\},
> \qquad 0\prec S_\flat\preceq P_{11}^{\mathrm{RE}}\preceq P_{11},
> \;}
> \]
> so \(\widetilde C_d\) is **compact**, with an explicit containment ellipsoid, for
> every \(d\). Two features of this are worth flagging now: it consumes **(H3a)** and
> nothing else — in particular **not (H3b)** — and \(\nabla^2\Psi\) does not involve
> \(\Lambda_\gamma\) at all, so the **shape of the certified set is prior-free** and
> only its centre moves with the prior. Both fail for the HPD family, whose curvature
> \(\nabla^2\bar\Phi=P_{11}-\nabla^2\Psi\) carries the prior explicitly and whose
> coercivity in the flat limit needs **(H3b)**. See Lemma 2A.3-G.

> **Range.** \(\varepsilon(\gamma)\in[0,1]\) always. The **upper** bound is automatic
> and has nothing to do with the model: \(q(\cdot\mid\gamma)\ge\varepsilon(\gamma)q_Q\)
> integrates to \(1\ge\varepsilon(\gamma)\), with equality iff
> \(q(\cdot\mid\gamma)=q_Q\). Everything below concerns the **lower** bound.


### B.1 Proper population prior ($\Lambda_\gamma\succ0$)


Throughout 2A assume **(H2)**, **(H3a)** and \(\Lambda_\gamma\succ0\), so that \(Q\)
exists by Cor. 1A.3. The order of business is: the inner program is well posed
(2A.1–2A.2), the profile is regular (2A.3), the profile is **concave in the log**
(2A.3-F) and **coercive in the log** (2A.3-G) — which is what makes its superlevel sets
usable as certified sets, convex and compact — and the minorization then holds on
\(\widetilde C_d\) (2A.4–2A.5).

**Lemma 2A.1 (strict convexity: any minimum is unique).** For each fixed \(\gamma\),
the map \(\gamma'\mapsto g(\gamma'\mid\gamma)\) is **strictly convex**:

\[
\nabla^2_{\gamma'}g(\gamma'\mid\gamma)
=\Lambda_Q-\Lambda_q(\gamma'\mid\gamma)
=\Lambda_Q\,\mathrm{Cov}\bigl(m(\beta)\mid\gamma,\gamma',y\bigr)\,\Lambda_Q
\;\succ\;0 .
\]

Hence \(g(\cdot\mid\gamma)\) has **at most one** minimiser.

*Proof.* The middle identity is Remark R2b-1: differentiating the mixture twice gives
\(\Lambda_q=\Lambda_Q-\Lambda_Q\mathrm{Cov}(m\mid\gamma,\gamma',y)\Lambda_Q\), where the
covariance is under the **bridge law** of \(\beta\) given both endpoints. Strictness:
\(m(\beta)=L\beta+a\) with \(L=P_{11}^{-1}[\mathcal W_1^\top P_b,\dots,\mathcal W_J^\top P_b]\)
of full row rank \(q\) under **(H3a)**, and the bridge law of \(\beta\) is
nondegenerate, so \(\mathrm{Cov}(m\mid\cdot)=L\,\mathrm{Cov}(\beta\mid\cdot)L^\top\succ0\);
conjugating by \(\Lambda_Q\succ0\) preserves this. A strictly convex function has at
most one minimiser. \(\square\)

**Lemma 2A.2 (a minimum exists, for every \(\gamma\)).** For every
\(\gamma\in\mathbb R^q\) the infimum defining \(\varepsilon(\gamma)\) is **attained**
at a (unique, by 2A.1) point \(\gamma'_\star(\gamma)\), and

\[
\varepsilon(\gamma)=e^{\,g(\gamma'_\star(\gamma)\mid\gamma)}\;\in\;(0,1] .
\]

*Proof.* Let \(P_\gamma\) be the law of \(m(\beta)\) under \(\pi(\beta\mid\gamma,y)\).
Dividing the two densities, the Gaussian normalisers and the \(\gamma'\)-quadratics
cancel:

\[
g(\gamma'\mid\gamma)=\log\int\exp\bigl\{(m-\gamma^\star)^\top\Lambda_Q\gamma'\bigr\}\,d\nu_\gamma(m),
\qquad
d\nu_\gamma(m)=e^{-\frac12m^\top\Lambda_Qm+\frac12{\gamma^\star}^\top\Lambda_Q\gamma^\star}\,dP_\gamma(m),
\]

so in the variable \(t=\Lambda_Q\gamma'\), **\(g\) is the cumulant generating function
of the finite positive measure \(\nu_\gamma\) recentred at \(\gamma^\star\)**. It is finite
everywhere (the Gaussian factor in \(\nu_\gamma\) dominates every exponential) and
convex. For coercivity, restrict the integral to a half-space: for any unit \(u\), any
\(c>0\), any \(s>0\),

\[
g\bigl(s\Lambda_Q^{-1}u\mid\gamma\bigr)\;\ge\;sc+\log\nu_\gamma\bigl(\{m:(m-\gamma^\star)^\top u\ge c\}\bigr)
\;\xrightarrow[s\to\infty]{}\;+\infty
\]

whenever that half-space carries positive \(\nu_\gamma\)-mass. Under **(H3a)**,
\(m(\beta)\) is an affine map of \(\beta\) whose linear part has rank \(q\), and
\(\pi(\beta\mid\gamma,y)>0\) on all of \(\mathbb R^{Jp_{\mathrm{re}}}\), so \(P_\gamma\)
— hence \(\nu_\gamma\) — has **full support \(\mathbb R^q\)** and every half-space
qualifies. A finite convex coercive function attains its minimum, which is therefore
\(>-\infty\); exponentiating gives \(\varepsilon(\gamma)>0\). \(\square\)

> **What this does and does not give.** It certifies that the minimum is **attained**;
> it produces **no value**. Outside Gaussian closure (§7.3.3A) there is no closed form
> for \(\varepsilon(\gamma)\) — it is the value of a strictly convex \(q\)-dimensional
> program, evaluated numerically. Remark R2b-0 discusses this.

**Lemma 2A.3 (\(\varepsilon\) is continuous).** \(\varepsilon:\mathbb R^q\to(0,1]\) is
**continuous**; in particular it is continuous on \(\widetilde C\).

*Proof.* Joint continuity of \(g\) gives **upper** semicontinuity of
\(\log\varepsilon(\gamma)=\inf_{\gamma'}g(\gamma'\mid\gamma)\) for free — an infimum
of continuous functions is upper semicontinuous. The other direction
needs the minimisers not to escape as \(\gamma\) varies, i.e. **locally uniform**
coercivity, which is the estimate below.

> **Uniform coercivity estimate (used again in 2B).** Let \(K\subset\mathbb R^q\) be
> compact. There are constants \(c_K,\eta_K>0\) with
> \[
> g(\gamma'\mid\gamma)\;\ge\;c_K\,|\gamma'|+\log\eta_K
> \qquad\text{for all }\gamma\in K,\ \gamma'\in\mathbb R^q .
> \]
> *Proof.* Work in \(\beta\)-space, where the mixing law \(\pi(\beta\mid\gamma,y)\) is
> a probability density that is **strictly positive everywhere**. Write
> \(t=\Lambda_Q\gamma'\), \(u=t/|t|\), \(m(\beta)=L\beta+a\), so that
> \((m(\beta)-\gamma^\star)^\top t=|t|\bigl(\beta^\top L^\top u+(a-\gamma^\star)^\top u\bigr)\).
>
> Since \(L\) has full row rank \(q\) under **(H3a)**, \(\sigma_0:=\sigma_{\min}(L)>0\),
> so \(|L^\top u|\ge\sigma_0\) for every unit \(u\) and \(v(u):=L^\top u/|L^\top u|\) is
> a well-defined unit vector in \(\beta\)-space depending continuously on \(u\). Put
> \(A_0:=|a-\gamma^\star|\) and fix \(\rho\) large enough that \(c:=\sigma_0\rho/2-A_0>0\). On the
> **cone-cap** \(B_v:=\{\beta:|\beta|\le\rho,\ \beta^\top v\ge\rho/2\}\) we have
> \(\beta^\top L^\top u=|L^\top u|\,\beta^\top v(u)\ge\sigma_0\rho/2\), hence
> \((m(\beta)-\gamma^\star)^\top u\ge c\). Restricting the integral to \(B_{v(u)}\),
> \[
> g(\gamma'\mid\gamma)\;\ge\;c\,|t|+\log\bigl(w_{\min}\,\pi(B_{v(u)}\mid\gamma,y)\bigr),
> \]
> where \(w_{\min}>0\) lower-bounds the continuous positive weight
> \(e^{-\frac12m(\beta)^\top\Lambda_Qm(\beta)+\frac12{\gamma^\star}^\top\Lambda_Q\gamma^\star}\) on the
> bounded set \(\{|\beta|\le\rho\}\).
>
> Finally \((v,\gamma)\mapsto\pi(B_v\mid\gamma,y)\) is continuous (dominated
> convergence; the cone-cap boundary is Lebesgue-null) and strictly positive
> (\(B_v\) has positive Lebesgue measure) on the **compact** set
> \(\{|v|=1\}\times K\), so
> \(p_K:=\min_{|v|=1,\gamma\in K}\pi(B_v\mid\gamma,y)>0\). Take
> \(\eta_K:=w_{\min}p_K\) and \(c_K:=c\,\lambda_{\min}(\Lambda_Q)\), using
> \(|t|\ge\lambda_{\min}(\Lambda_Q)|\gamma'|\). \(\;\square\)
>
> *(This supersedes the total-variation argument of an earlier draft: no continuity of
> \(\gamma\mapsto\nu_\gamma\) in TV is needed, and — crucially for 2B — the same
> estimate holds uniformly in a varying prior, since \(\pi(\beta\mid\gamma,y)\) does not
> depend on \(\Lambda_\gamma\) at all.)*

Since \(\log\varepsilon\le0\) always, the estimate forces every minimiser to satisfy
\(|\gamma'_\star(\gamma)|\le R_K:=\bigl(\log(1/\eta_K)\bigr)_+/c_K\) — **the same
\(R_K\) for all \(\gamma\in K\)**. Hence on \(K\),
\(\log\varepsilon(\gamma)=\min_{|\gamma'|\le R_K}g(\gamma'\mid\gamma)\) is a minimum of
a jointly continuous function over a **fixed compact** set, and Berge's maximum theorem
(constant, compact constraint correspondence) makes it continuous. As \(K\) was
arbitrary, \(\varepsilon\in C(\mathbb R^q)\). \(\square\)

> **By-product.** The radius \(R_K\) is explicit, so evaluating \(\varepsilon(\gamma)\)
> numerically is a **bounded** optimisation.

**Remark 2A.3-E (envelope identity: the slope of \(\mathcal D\) is the slope of
\(\log q\)).** Write \(\mathcal D(\gamma)=g(\gamma'_\star(\gamma)\mid\gamma)
=\log\varepsilon(\gamma)\). The minimiser is unique (2A.1), interior, and locally
bounded (2A.3), and \(\nabla_\gamma g\) is continuous, so **Danskin's envelope
theorem** applies: \(\mathcal D\in C^1\) and

\[
\boxed{\;
\nabla\mathcal D(\gamma)
\;=\;\nabla_\gamma g\bigl(\gamma'\mid\gamma\bigr)\Big|_{\gamma'=\gamma'_\star(\gamma)}
\;=\;\nabla_\gamma\log q\bigl(\gamma'\mid\gamma\bigr)\Big|_{\gamma'=\gamma'_\star(\gamma)} .
\;}
\]

Two things make this work. The term \(-\log q_Q(\gamma')\) does not involve \(\gamma\)
at all, so it contributes nothing to \(\nabla_\gamma\). And the chain-rule term through
\(\gamma'_\star(\gamma)\) — which would otherwise require differentiating the
minimiser, an object with no formula — vanishes because
\(\nabla_{\gamma'}g=0\) there. **One differentiates as if the minimiser were frozen.**

*Second order.* Differentiating again (Danskin, with \(g_{\gamma'\gamma'}\succ0\) by
2A.1),

\[
\nabla^2\mathcal D(\gamma)
=\underbrace{g_{\gamma\gamma}}_{=\;\nabla^2_\gamma\log q(\gamma'\mid\gamma)}
-\;\underbrace{g_{\gamma\gamma'}\,g_{\gamma'\gamma'}^{-1}\,g_{\gamma'\gamma}}_{\succeq\,0}
\;\preceq\;\nabla^2_\gamma\log q\bigl(\gamma'_\star(\gamma)\mid\gamma\bigr),
\]

all derivatives evaluated at \(\gamma'=\gamma'_\star(\gamma)\). This gives a clean
**sufficient condition for the superlevel sets of \(\varepsilon\) to be convex**:

> \(\mathcal D\) is concave — hence \(\{\varepsilon\ge\varepsilon_\ast\}\) is convex for
> every \(\varepsilon_\ast\) — **whenever the transition density
> \(\gamma\mapsto q(\gamma'\mid\gamma)\) is log-concave in its *source* for each fixed
> destination \(\gamma'\).**

Whether that holds is settled — affirmatively, for every family — in **Lemma 2A.3-F**,
which is what licenses the certified family of the preamble.

**Lemma 2A.3-F (the deficiency is concave; closed form for both derivatives).** Under
**(H2)**, \(\mathcal D\in C^2\) with

\[
\boxed{\;
\nabla\mathcal D(\gamma)=P_{11}\bigl(\gamma^\star-M(\gamma)\bigr),
\qquad
\nabla^2\mathcal D(\gamma)=-P_{11}J(\gamma)
=-\sum_{j}H_j^\top P_b\,V_j(\gamma)\,P_bH_j\;\preceq\;0 .
\;}
\]

Hence:

1. \(\mathcal D\) is **concave** on \(\mathbb R^q\), and strictly so wherever
   \(\sum_jH_j^\top P_bV_jP_bH_j\succ0\) — in particular under **(H3a)** with
   \(V_j\succ0\). Therefore **every superlevel set
   \(\{\gamma:\varepsilon(\gamma)\ge\varepsilon_\ast\}\) is convex.**
2. \(\nabla\mathcal D(\gamma)=0\iff M(\gamma)=\gamma^\star\iff\gamma=\gamma^\star\)
   (uniqueness of the fixed point, Lemma R2d). So \(\varepsilon\) attains its **global
   maximum exactly at the refresh centre** — the Centering remark's item 1, which was
   previously argued only under closure, is a theorem.
3. Comparing with \(\nabla\Phi(\gamma)=P_{11}(\gamma-M(\gamma))\) (Cor. R2d-5(2)) gives
   \(\nabla(\Phi-\mathcal D)(\gamma)=P_{11}(\gamma-\gamma^\star)\), i.e. the **exact
   identity**
   \[
   \boxed{\;
   \mathcal D(\gamma)-\mathcal D(\gamma^\star)
   \;=\;\bar\Phi(\gamma)-\tfrac12\|\gamma-\gamma^\star\|^2_{P_{11}} .
   \;}
   \]
   Equivalently \(\nabla^2\mathcal D=\nabla^2\Phi-P_{11}\). The inscribed-ellipsoid
   bound \(\bar\Phi\le\frac12\|\gamma-\gamma^\star\|^2_{P_{11}}\) used for
   \(Q(\widetilde C_d)\) is **exactly** the statement
   \(\mathcal D(\gamma)\le\mathcal D(\gamma^\star)\) of item 2.

*Proof.* Write \(\pi(\beta\mid\gamma,y)=e^{L(\beta,\gamma)}/Z(\gamma)\) with
\(L(\beta,\gamma)=\sum_j\ell_j(\beta_j)-\frac12\sum_j(\beta_j-H_j\gamma)^\top
P_b(\beta_j-H_j\gamma)\).

**(i) Score.** \(\nabla_\gamma L=\sum_jH_j^\top P_b\beta_j-P_{11}^{\mathrm{RE}}\gamma
=P_{11}m(\beta)-\Lambda_\gamma\mu_0-P_{11}^{\mathrm{RE}}\gamma\), by definition of
\(m\). Since \(\nabla_\gamma\log Z=E_{\beta\mid\gamma}[\nabla_\gamma L]\), subtracting
gives the **score identity**
\[
s(\beta,\gamma):=\nabla_\gamma\log\pi(\beta\mid\gamma,y)=P_{11}\bigl(m(\beta)-M(\gamma)\bigr).
\]

**(ii) Bridge representation.** Differentiating
\(q(\gamma'\mid\gamma)=\int\phi_q(\gamma';m(\beta),\Sigma^\star)\pi(\beta\mid\gamma,y)d\beta\)
under the integral and dividing by \(q\),
\[
\nabla_\gamma\log q(\gamma'\mid\gamma)=E\bigl[s(\beta,\gamma)\bigm|\gamma,\gamma',y\bigr]
=P_{11}\bigl(\bar m(\gamma,\gamma')-M(\gamma)\bigr),
\]
where the expectation is under the **bridge law**
\(\pi(\beta\mid\gamma,\gamma',y)\propto\phi_q(\gamma';m(\beta),\Sigma^\star)
\pi(\beta\mid\gamma,y)\) of Remark R2b-1 and
\(\bar m(\gamma,\gamma'):=E[m(\beta)\mid\gamma,\gamma',y]\).

**(iii) Stationarity in \(\gamma'\).** Likewise
\(\nabla_{\gamma'}\log q=P_{11}(\bar m-\gamma')\) and
\(\nabla_{\gamma'}\log q_Q=-P_{11}(\gamma'-\gamma^\star)\), so
\[
\nabla_{\gamma'}g(\gamma'\mid\gamma)=P_{11}\bigl(\bar m(\gamma,\gamma')-\gamma^\star\bigr),
\]
and the minimiser is characterised by the strikingly simple condition
\(\bar m\bigl(\gamma,\gamma'_\star(\gamma)\bigr)=\gamma^\star\): **the worst destination
is the one that makes the bridge mean land on the refresh centre.**

**(iv) First derivative.** By the envelope identity 2A.3-E and (ii)–(iii),
\(\nabla\mathcal D(\gamma)=P_{11}(\bar m-M(\gamma))\) evaluated at
\(\gamma'=\gamma'_\star(\gamma)\), where \(\bar m=\gamma^\star\). Hence
\(\nabla\mathcal D(\gamma)=P_{11}(\gamma^\star-M(\gamma))\).

**(v) Second derivative.** Differentiate (iv): only \(M\) depends on \(\gamma\), so
\(\nabla^2\mathcal D=-P_{11}\,\partial M/\partial\gamma=-P_{11}J(\gamma)\), and
\(P_{11}J=\sum_jH_j^\top P_bV_jP_bH_j\succeq0\) by Lemma R2d(a). \(\square\)

> **Why the earlier "difference of concave functions" reasoning was inconclusive.**
> Bounding \(\nabla^2\mathcal D\preceq\nabla^2_\gamma\log q\) (2A.3-E) and then bounding
> \(\nabla^2_\gamma\log q\) by Prékopa throws away the Schur complement
> \(g_{\gamma\gamma'}g_{\gamma'\gamma'}^{-1}g_{\gamma'\gamma}\). Here that term is not a
> harmless correction: one computes
> \(g_{\gamma\gamma'}=g_{\gamma'\gamma'}=P_{11}\mathrm{Cov}(m\mid\gamma,\gamma',y)P_{11}\),
> so the Schur complement **exactly cancels** the indefinite part of
> \(g_{\gamma\gamma}\), leaving \(-P_{11}J\). Equivalently
> \(\partial\gamma'_\star/\partial\gamma=-(g_{\gamma'\gamma'})^{-1}g_{\gamma'\gamma}=-I\):
> the bridge law depends on \((\gamma,\gamma')\) only through \(\gamma+\gamma'\), so the
> worst destination reflects through a fixed point, \(\gamma'_\star(\gamma)+\gamma\)
> being constant. Under Gaussian closure this reads \(\gamma'_\star(\gamma)=2\gamma^\star-\gamma\).

**The Prékopa route, for contrast, stops short.** Write the joint exponent

\[
L(\beta,\gamma):=\sum_j\ell_j(\beta_j)
-\tfrac12\sum_j(\beta_j-H_j\gamma)^\top P_b(\beta_j-H_j\gamma),
\]

which is **jointly concave** in \((\beta,\gamma)\) under **(H2)** — each \(\ell_j\) is
concave and the quadratic is a negative semidefinite form in the linear map
\(\beta_j-H_j\gamma\). Since \(\pi(\beta\mid\gamma,y)=e^{L(\beta,\gamma)}/Z(\gamma)\),

\[
\log q(\gamma'\mid\gamma)=\log N(\gamma',\gamma)-\log Z(\gamma),
\qquad
\begin{aligned}
N(\gamma',\gamma)&:=\int\phi_q\bigl(\gamma';m(\beta),\Sigma^\star\bigr)e^{L(\beta,\gamma)}d\beta,\\
Z(\gamma)&:=\int e^{L(\beta,\gamma)}d\beta .
\end{aligned}
\]

Both are log-concave in \(\gamma\) by **Prékopa** — the integrand of \(N\) is jointly
log-concave in \((\beta,\gamma)\) because \(m\) is affine, and \(Z\) is the same
statement without the Gaussian factor. So

\[
\log q(\gamma'\mid\gamma)=\underbrace{\log N(\gamma',\gamma)}_{\text{concave in }\gamma}
-\underbrace{\log Z(\gamma)}_{\text{concave in }\gamma},
\]

a **difference of concave functions**, which a curvature bound alone cannot resolve.
Taking the minimum over \(\gamma'\) (an infimum of a family concave in \(\gamma\), hence
concave) isolates the whole difficulty in a single term:

\[
\boxed{\;
\mathcal D(\gamma)\;=\;\underbrace{\mathcal C(\gamma)}_{\text{concave}}
\;+\;\underbrace{\Phi(\gamma)}_{\text{convex}}
\;-\;\tfrac12\|\gamma-\mu_0\|^2_{\Lambda_\gamma}\;+\;\mathrm{const},
\qquad
\mathcal C(\gamma):=\min_{\gamma'}\bigl[\log N(\gamma',\gamma)-\log q_Q(\gamma')\bigr],
\;}
\]

using \(-\log Z=\Phi-\tfrac12\|\gamma-\mu_0\|^2_{\Lambda_\gamma}+\mathrm{const}\)
(Cor. R2d-5(1)). Read on its own this looks inconclusive — \(\mathcal C\) is concave
but \(\Phi\) is convex — and that is as far as log-concavity arguments get.

**Lemma 2A.3-F resolves it, and identifies \(\mathcal C\) exactly.** Matching
\(\nabla^2\mathcal D=\nabla^2\Phi-P_{11}\) against the display gives
\(\nabla^2\mathcal C=-P_{11}+\Lambda_\gamma=-P_{11}^{\mathrm{RE}}\): the concave part
is a **pure quadratic**, \(\mathcal C(\gamma)=-\frac12\gamma^\top
P_{11}^{\mathrm{RE}}\gamma+\text{affine}\), with no dependence on the likelihood
curvature at all. So the resolution is not that \(\Phi\)'s convexity cancels, but that
it is **capped**:

\[
\nabla^2\mathcal D=\nabla^2\Phi-P_{11}\preceq0
\iff
\nabla^2\Phi\preceq P_{11}
\iff
J(\gamma)\succeq0
\iff
V_j(\gamma)\succeq0 ,
\]

the last being automatic since \(V_j\) is a covariance. The posterior potential is
convex, but never *more* convex than the refresh precision — and \(\mathcal D\)
subtracts exactly \(P_{11}\). Note this is the **same** inequality
\(\nabla^2\Phi\preceq P_{11}\) that inscribes the \(Q\)-ellipsoid in
\(\widetilde C_d\); by item 3 of 2A.3-F the two statements are literally equivalent.

**Under Gaussian closure the closed form agrees.** Lemma R2b′ gives
\(\mathcal D(\gamma)=\tfrac12\log\det(\Sigma^\star\Sigma^{-1})
-\tfrac12\bigl(M(\gamma)-\gamma^\star\bigr)^\top(\Sigma-\Sigma^\star)^{-1}
\bigl(M(\gamma)-\gamma^\star\bigr)\), a concave quadratic with
\(\nabla^2\mathcal D=-J^\top(\Sigma-\Sigma^\star)^{-1}J=-P_{11}J\), matching 2A.3-F —
the last equality because \(\Sigma-\Sigma^\star=JP_{11}^{-1}\) under closure. The
superlevel sets \(\{\varepsilon\ge\varepsilon_\ast\}\) are then ellipsoids in the metric
\(P_{11}J\) — a **third** metric, distinct from both \(\nabla^2\Phi\) (which shapes
\(\widetilde C_d\)) and \(P_{11}\) (which shapes \(Q\)). See Remark R2b-2.

**Lemma 2A.3-G (the deficiency gap is coercive; the certified sets are compact).**
Assume **(H1)**, **(H2)**, **(H3a)**. Write

\[
S(\gamma)\;:=\;\nabla^2\Psi(\gamma)\;=\;P_{11}J(\gamma)
\;=\;\sum_j H_j^\top P_b\,V_j(\gamma)\,P_bH_j ,
\qquad
c\;:=\;\min_{|u|=1}\int_0^1 u^\top S\bigl(\gamma^\star+tu\bigr)u\;dt .
\]

Then \(c>0\) and

\[
\boxed{\;
\Psi(\gamma)\;\ge\;c\bigl(|\gamma-\gamma^\star|-1\bigr)
\quad\text{for }|\gamma-\gamma^\star|\ge1,
\qquad\text{hence}\qquad
\widetilde C_d\;\subseteq\;\bigl\{\gamma:\ |\gamma-\gamma^\star|\le 1+d/c\bigr\} .
\;}
\]

So \(\Psi\) is **coercive** and \(\widetilde C_d\) is **compact** — closed by 2A.3,
convex by 2A.3-F, bounded here — for every \(d>0\), with an explicit containment
radius. Moreover:

1. **A two-sided ellipsoid sandwich.** Brascamp–Lieb gives \(V_j\preceq P_b^{-1}\),
   hence \(S(\gamma)\preceq P_{11}^{\mathrm{RE}}\) and
   \[
   \bigl\{\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\le2d\bigr\}
   \;\subseteq\;\widetilde C_d .
   \]
   If in addition the group curvature is bounded, \(-\nabla^2\ell_j=\mathcal G_j(\beta_j)
   \preceq G_j\) for some fixed \(G_j\) (automatic whenever the GLM weights are bounded
   above: binomial with any link gives \(G_j=\frac14Z_j^\top N_jZ_j\), Gaussian gives
   \(G_j=\mathcal G_j\) exactly), then Cramér–Rao gives
   \(V_j\succeq(P_b+G_j)^{-1}\), so with
   \(S_\flat:=\sum_jH_j^\top P_b(P_b+G_j)^{-1}P_bH_j\succ0\),
   \[
   \boxed{\;
   \bigl\{\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\le2d\bigr\}
   \;\subseteq\;\widetilde C_d\;\subseteq\;
   \bigl\{\|\gamma-\gamma^\star\|^2_{S_\flat}\le2d\bigr\},
   \;}
   \]
   two **concentric ellipsoids** with \(0\prec S_\flat\preceq S(\gamma)\preceq
   P_{11}^{\mathrm{RE}}\preceq P_{11}\). The certified set is squeezed between them,
   and \(\Psi\) is \(S_\flat\)-strongly convex.
2. **The shape is prior-free.** \(V_j(\gamma)=\mathrm{Cov}(\beta_j\mid\gamma,y)\) does
   not involve \(\Lambda_\gamma\) — the population prior never enters the
   \(\beta\)-block — so \(S(\gamma)=\nabla^2\Psi(\gamma)\), and with it \(c\),
   \(S_\flat\) and the whole sandwich, are **independent of the population prior**.
   Only the centre \(\gamma^\star\) moves with it. Contrast
   \(\nabla^2\bar\Phi=P_{11}-S\), which carries \(\Lambda_\gamma\) explicitly.
3. **The inner ellipsoid is Lemma 1A.2 in disguise.** By complementarity,
   \(S\preceq P_{11}^{\mathrm{RE}}\) is the *same* inequality as
   \(\bar\Phi\ge\frac12\|\gamma-\gamma^\star\|^2_{\Lambda_\gamma}\), the
   \(\Lambda_\gamma\)-strong convexity of \(\Phi\). The coarser bound
   \(S\preceq P_{11}\) is what gives the \(\chi^2\) refresh-mass bound of the preamble;
   \(P_{11}^{\mathrm{RE}}\) sharpens it.

*Proof.* By 2A.3-F, \(\nabla\Psi(\gamma^\star)=0\) and
\(\nabla^2\Psi=P_{11}J=\sum_jH_j^\top P_bV_jP_bH_j=S\succeq0\).

**(i) \(S(\gamma)\succ0\) for every \(\gamma\).** Under **(H1)**+**(H2)** the law
\(\pi(\beta_j\mid\gamma,y)\propto e^{\ell_j(\beta_j)}\phi(\beta_j;H_j\gamma,\Psi)\) is a
proper log-concave density that is strictly positive on all of
\(\mathbb R^{p_{\mathrm{re}}}\); it is therefore not supported in any hyperplane and has
finite second moments, so \(V_j(\gamma)\succ0\). Given \(u\ne0\), **(H3a)** makes
\(u^\top P_{11}^{\mathrm{RE}}u=\sum_j|H_ju|^2_{P_b}>0\), so \(H_ju\ne0\) for some \(j\);
then \(v:=P_bH_ju\ne0\) and \(u^\top S(\gamma)u\ge v^\top V_j(\gamma)v>0\).

**(ii) Linear growth along every ray.** Fix a unit \(u\) and put
\(\psi_u(s):=\Psi(\gamma^\star+su)\). Then \(\psi_u'(s)=\int_0^su^\top S(\gamma^\star+tu)u\,dt\),
which is **nondecreasing** in \(s\) because the integrand is \(\ge0\). Hence for
\(s\ge1\), \(\psi_u'(s)\ge\psi_u'(1)\ge c\), and \(\psi_u(s)=\int_0^s\psi_u'\ge c(s-1)\).

**(iii) \(c>0\).** \(V_j\) is the Hessian of a log-partition function (the exponential
tilt of Lemma R2d(a)), hence smooth in \(\gamma\); so
\((u,t)\mapsto u^\top S(\gamma^\star+tu)u\) is continuous and, by (i), strictly positive
on the compact \(\{|u|=1\}\times[0,1]\). Its \(t\)-integral is therefore a continuous
strictly positive function of \(u\), whose minimum over the compact unit sphere is
attained and positive.

Combining, \(\Psi(\gamma)\ge c(|\gamma-\gamma^\star|-1)\), so \(\Psi(\gamma)\le d\)
forces \(|\gamma-\gamma^\star|\le1+d/c\). For item 1, \(V_j\preceq P_b^{-1}\) is
Lemma R2d(b) and \(V_j\succeq(P_b+G_j)^{-1}\) is the Cramér–Rao bound for a location
family: the negative log-density of \(\pi(\beta_j\mid\gamma,y)\) has Hessian
\(-\nabla^2\ell_j+P_b\preceq G_j+P_b\), and its Fisher information for location is
\(E[\nabla U\nabla U^\top]=E[\nabla^2U]\preceq G_j+P_b\), so
\(V_j\succeq(G_j+P_b)^{-1}\). Integrating the resulting constant Hessian bounds twice
from \(\gamma^\star\), where \(\Psi\) and \(\nabla\Psi\) both vanish, gives
\(\frac12\|\gamma-\gamma^\star\|^2_{S_\flat}\le\Psi(\gamma)\le
\frac12\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\). \(\square\)

> **Why this is easy here and was hard for the HPD family.** Coercivity of \(\Psi\)
> reduces, via 2A.3-F, to \(S=\sum_jH_j^\top P_bV_jP_bH_j\) being bounded away from
> \(0\) — a statement about **conditional covariances of the random-effect block**,
> which are nondegenerate for the trivial reason that a positive log-concave density on
> \(\mathbb R^{p_{\mathrm{re}}}\) has a nondegenerate covariance. Coercivity of
> \(\bar\Phi\) is the opposite kind of statement: it needs \(\nabla^2\bar\Phi=P_{11}-S\)
> bounded **below**, i.e. the conditional covariances to be bounded **away from their
> Brascamp–Lieb ceiling** \(P_b^{-1}\), which the prior supplies when
> \(\Lambda_\gamma\succ0\) and which otherwise requires **(H3b)** — no separated group
> may flatten the potential at infinity. That asymmetry is why boundedness of the
> certified set costs **(H3b)** on the HPD route and nothing here.

**What this settles.** Concavity of \(\mathcal D\) holds for every family, so the
superlevel sets \(\{\varepsilon\ge\varepsilon_\ast\}\) are convex, and by 2A.3-F(3)
each is the sublevel set \(\{\Psi\le\mathcal D(\gamma^\star)-\log\varepsilon_\ast\}\)
of the convex nonnegative \(\Psi=\frac12\|\cdot-\gamma^\star\|^2_{P_{11}}-\bar\Phi\),
which 2A.3-G shows is coercive. So the certified family is **convex and compact**, with
an explicit containment ellipsoid, and Weierstrass therefore applies to it: Lemma 2A.4
takes the minimum of \(\varepsilon\) over \(\widetilde C_d\) and finds it equal to
\(e^{-d}\varepsilon(\gamma^\star)\), attained on the boundary. The alternative of
indexing by \(\bar\Phi\) — HPD regions — is a *different* convex compact family,
discussed for contrast in Remark R2b-0; the same Weierstrass argument applies there but
returns an unevaluated number, and its compactness needs **(H3b)** in the flat limit.

**Lemma 2A.4 (the constant on \(\widetilde C_d\) is the attained minimum).** For every
\(d>0\), with \(\varepsilon_d:=e^{-d}\varepsilon(\gamma^\star)\),

\[
\boxed{\;
\min_{\gamma\in\widetilde C_d}\varepsilon(\gamma)
\;=\;\varepsilon_d\;=\;e^{-d}\varepsilon(\gamma^\star)\;\in\;(0,1),
\;}
\]

the minimum being **attained**, on the boundary
\(\partial\widetilde C_d=\{\gamma:\Psi(\gamma)=d\}\ne\varnothing\). Consequently
\(\varepsilon(\gamma)\ge\varepsilon_d\) for all \(\gamma\in\widetilde C_d\), the
constant is **sharp** — no larger one minorizes on \(\widetilde C_d\) — and
\(\widetilde C_d\) is the **largest** set carrying it. Also \(d\mapsto\varepsilon_d\)
is strictly decreasing and \(d\mapsto\widetilde C_d\) nondecreasing.

*Proof.* **Existence of the minimum: continuity plus compactness.** \(\varepsilon\) is
continuous (Lemma 2A.3) and \(\widetilde C_d\) is compact (Lemma 2A.3-G), so by
Weierstrass the minimum is attained at some \(\gamma_\varepsilon\in\widetilde C_d\).

**Its value.** \(\ge\varepsilon_d\) holds because every point of \(\widetilde C_d\)
satisfies \(\varepsilon\ge\varepsilon_d\) by definition of the superlevel set. For
\(\le\), the boundary is nonempty: fix a unit \(u\); the map
\(s\mapsto\Psi(\gamma^\star+su)\) is continuous with value \(0\) at \(s=0\) and
\(\to\infty\) as \(s\to\infty\) (coercivity, Lemma 2A.3-G), so by the intermediate value
theorem there is \(s_d\) with \(\Psi(\gamma^\star+s_du)=d\) exactly, i.e.
\(\varepsilon=\varepsilon_d\) there. That point lies in \(\widetilde C_d\), so the
minimum is \(\le\varepsilon_d\). Hence \(\min=\varepsilon_d\), attained on
\(\partial\widetilde C_d\), and the minorization constant cannot be improved.

**Positivity and range.** \(\varepsilon(\gamma^\star)>0\) is Lemma 2A.2 at the single
point \(\gamma^\star\), and \(\varepsilon(\gamma^\star)\le1\) is the Range box, so
\(\varepsilon_d\le e^{-d}<1\). Maximality of the set is the definition of a superlevel
set. \(\square\)

> **What compactness buys, stated precisely.** An earlier draft of this note claimed
> the superlevel construction "removes Weierstrass". That is wrong, and the correction
> matters. \(\varepsilon_d\) **is** the minimum of \(\varepsilon\) over
> \(\widetilde C_d\), and continuity (2A.3) plus compactness (2A.3-G) is exactly the
> proof that the minimum exists and is attained. What the superlevel construction
> changes is not whether a minimum is taken but whether its **value is known**:
> \[
> \min_{\widetilde C_d}\varepsilon=e^{-d}\varepsilon(\gamma^\star)
> \qquad\text{versus}\qquad
> \min_{\widetilde C^{\mathrm{HPD}}_d}\varepsilon=\;?
> \]
> On the HPD route Weierstrass produces a number with no known relation to anything
> else in the construction, and the entire \(d\)-dependence of the certificate is
> hidden inside it. On the superlevel route Weierstrass produces
> \(e^{-d}\varepsilon(\gamma^\star)\), exactly, because the set was indexed by the very
> function being minimised — the minimum sits on the level that names the set. The
> \(d\)-dependence is then explicit and the unknown collapses to **one scalar at one
> point**, \(\varepsilon(\gamma^\star)\), still supplied by Lemma 2A.2.
>
> So compactness is load-bearing here, in three separate places: it gives the minimum
> (this lemma), it makes the rate \(\lambda_{C_d}=\max_{C_d}T\) an attained maximum
> (Appendix A.2), and it bounds the numerical search for \(\varepsilon(\gamma^\star)\)
> (Remark R2b-0). What *did* get cheaper is its **price**: Lemma 2A.3-G proves it from
> nondegeneracy of the random-effect conditional covariances, which needs **(H3a)**
> alone and is prior-free, where coercivity of \(\bar\Phi\) needs **(H3b)** in the flat
> limit.

**Lemma 2A.5 (the minorization condition holds on \(\widetilde C_d\)).** Let
\(\varepsilon_d\) be as in 2A.4. Then:

**(a) Unrestricted chain, untruncated refresh.** The density bound

\[
q(\gamma'\mid\gamma)\;\ge\;\varepsilon_d\,q_Q(\gamma')
\qquad\text{for all }\gamma\in\widetilde C_d\ \text{and all }\gamma'\in\mathbb R^q
\]

holds, and integrating it over an arbitrary Borel \(A\subseteq\mathbb R^q\) gives

\[
\boxed{\;
q(\gamma,A)\;\ge\;\varepsilon_d\,Q(A)
\qquad\forall\,\gamma\in\widetilde C_d,\ \forall\,A\in\mathcal B(\mathbb R^q).
\;}
\]

That is, \(\widetilde C_d\) is a **small set of the original \(\gamma\)-chain**, with
minorizing constant \(\varepsilon_d\) and minorizing measure the **untruncated**
Gaussian \(Q\). Nothing is truncated on either side and there is **no mass factor**.

**(b) Restricted chain, truncated refresh.** Let \(q(\cdot,\cdot\mid\widetilde C_d)\)
be the restricted kernel of **Lemma R2c** — the two-block sweep of the target
truncated to \(\widetilde C_d\), which by construction never leaves
\(\widetilde C_d\) — and let \(Q_{\widetilde C_d}(A):=Q(A\cap\widetilde
C_d)/Q(\widetilde C_d)\) be \(Q\) conditioned on \(\widetilde C_d\). Then

\[
\boxed{\;
q\bigl(\gamma,A\mid\widetilde C_d\bigr)\;\ge\;\varepsilon\,Q_{\widetilde C_d}(A),
\qquad
\varepsilon\;:=\;\varepsilon_d\,Q(\widetilde C_d),
\qquad\forall\,\gamma\in\widetilde C_d,\ A\subseteq\widetilde C_d .
\;}
\]

Both sides now live on the **same** state space \(\widetilde C_d\): this is the
Doeblin condition in the form Rosenthal's theorem requires, with \(\varepsilon\in(0,1]\)
and refresh measure the truncated normal \(Q_{\widetilde C_d}\). Both factors are now
explicit in \(d\), the mass by the inscribed-ellipsoid bound of the preamble:

\[
\varepsilon\;=\;\varepsilon_d\,Q\bigl(\widetilde C_d\bigr)
\;\ge\;\varepsilon(\gamma^\star)\;e^{-d}\;\Pr\bigl(\chi^2_q\le 2d\bigr) .
\]

Everything on the right is a computable function of \(d\) except the single number
\(\varepsilon(\gamma^\star)\) — one \(q\)-dimensional convex program, solved once, at
one point. That is the quantitative content of certifying the superlevel set: the
unevaluated object is no longer a minimum over a set but a value at a point.

The proof of (b) uses nothing about \(\widetilde C_d\) beyond the density bound of (a)
holding on it, so it applies verbatim to **any** measurable \(C\subseteq\widetilde C_d\)
with \(Q(C)>0\), giving \(q(\gamma,A\mid C)\ge\varepsilon_d Q(C)\,Q_C(A)\) for
\(\gamma\in C\), \(A\subseteq C\). This is the form §7.2 Step 7 invokes, where \(C\) is
the intersection of \(\widetilde C_d\) with the spectral set \(E_{\lambda^\star}\).

> **Why (a) and (b) are different statements, and why the factor appears only in (b).**
> Truncating \(Q\) is not a free reformulation of (a): a minorizing measure must be a
> probability measure **on the state space of the chain being minorized**. In (a) that
> space is \(\mathbb R^q\) and \(Q\) already qualifies, so the constant is the full
> \(\varepsilon_d\). In (b) the sampler has been confined to \(\widetilde C_d\), so the
> refresh must be too; renormalising \(Q\) onto \(\widetilde C_d\) *divides* by
> \(Q(\widetilde C_d)\le1\), and multiplying the constant by the same
> \(Q(\widetilde C_d)\) is what pays for it. Concretely, the \(Q\)-mass
> \(\varepsilon_d\,Q(\widetilde C_d^{\,c})\) that (a) would have deposited outside
> \(\widetilde C_d\) is simply **discarded** — it is unavailable to a chain that cannot
> go there. That discarded mass *is* the factor. So (b) is the weaker statement, and it
> is the one Proposition R2 uses, because Proposition R2 certifies the **restricted**
> sampler.

*Proof.* **(a)** By definition of the infimum, \(q(\gamma'\mid\gamma)\ge
\varepsilon(\gamma)q_Q(\gamma')\) for **every** pair \((\gamma,\gamma')\); on
\(\widetilde C_d\) we have \(\varepsilon(\gamma)\ge\varepsilon_d\) by 2A.4. The bound
holds for **all** \(\gamma'\in\mathbb R^q\), not merely
\(\gamma'\in\widetilde C_d\) — nothing is truncated on the \(\gamma'\) side — so it may
be integrated over any Borel \(A\).

**(b)** Write the restricted \(\gamma\)-step as: draw \(\beta\sim\pi(\beta\mid\gamma,y)\)
exactly, then \(\gamma'\sim N(m(\beta),\Sigma^\star)\) **truncated to**
\(\widetilde C_d\). Its transition density is

\[
q\bigl(\gamma'\mid\gamma;\widetilde C_d\bigr)
=\mathbf 1\{\gamma'\in\widetilde C_d\}
\int\frac{\phi_q\bigl(\gamma';m(\beta),\Sigma^\star\bigr)}{N(\beta)}\,
\pi(\beta\mid\gamma,y)\,d\beta,
\qquad
N(\beta):=\!\!\int_{\widetilde C_d}\!\!\phi_q\bigl(\gamma';m(\beta),\Sigma^\star\bigr)d\gamma'.
\]

Since \(N(\beta)\le1\) for every \(\beta\), the integrand is pointwise at least
\(\phi_q(\gamma';m(\beta),\Sigma^\star)\), whence

\[
q\bigl(\gamma'\mid\gamma;\widetilde C_d\bigr)\;\ge\;q(\gamma'\mid\gamma)
\;\ge\;\varepsilon_d\,q_Q(\gamma')
\qquad\text{for }\gamma\in\widetilde C_d,\ \gamma'\in\widetilde C_d
\]

— truncation only *raises* the density on the retained region. Integrating over
\(A\subseteq\widetilde C_d\) gives \(q(\gamma,A\mid\widetilde C_d)\ge\varepsilon_d Q(A)
=\varepsilon_d Q(\widetilde C_d)\,Q_{\widetilde C_d}(A)\), using
\(Q(A)=Q(\widetilde C_d)Q_{\widetilde C_d}(A)\) for \(A\subseteq\widetilde C_d\).
Finally \(Q(\widetilde C_d)>0\) because \(\widetilde C_d\) contains a ball around
\(\gamma^\star\) and \(q_Q>0\); quantitatively, the inscribed ellipsoid
\(\{(\gamma-\gamma^\star)^\top P_{11}(\gamma-\gamma^\star)\le2d\}\subseteq\widetilde C_d\) of the
preamble carries \(Q\)-mass \(\Pr(\chi^2_q\le2d)\), since
\((\gamma-\gamma^\star)^\top P_{11}(\gamma-\gamma^\star)\sim\chi^2_q\) under \(Q\). \(\square\)

> **A free sharpening of (b).** The step \(N(\beta)\le1\) is wasteful. Keeping the
> normaliser gives \(\varepsilon=\varepsilon_d\,Q(\widetilde C_d)\big/
> \sup_{\beta}N(\beta)\ge\varepsilon_d\,Q(\widetilde C_d)\), and
> \(\sup_\beta N(\beta)<1\) strictly whenever no \(m(\beta)\) places essentially all
> its Gaussian mass inside \(\widetilde C_d\). We do not use this, since Proposition R2
> claims no formula for \(\varepsilon\); it is recorded to show that the
> \(Q(\widetilde C_d)\) factor is an upper bound on the true cost of restricting, not
> an exact accounting of it.


### B.2 Flat-prior limit ($\Lambda_\gamma\downarrow0$)


Now let the population prior vary. Everything in 2A carries a hidden dependence on
\(\Lambda_\gamma\) through \(Q\); make it explicit by writing

\[
g\bigl(\gamma'\mid\gamma;\Lambda_\gamma\bigr)
:=\log q_\tau(\gamma'\mid\gamma)-\log q_{Q_\tau}(\gamma'),
\qquad
\varepsilon\bigl(\gamma\mid\Lambda_\gamma\bigr):=
\exp\Bigl\{\min_{\gamma'}g(\gamma'\mid\gamma;\Lambda_\gamma)\Bigr\},
\]

\[
\varepsilon_d\bigl(\Lambda_\gamma\bigr)
:=e^{-d}\,\varepsilon\bigl(\gamma^\star(\tau)\mid\Lambda_\gamma\bigr),
\qquad
\widetilde C_d\bigl(\Lambda_\gamma\bigr)
:=\bigl\{\gamma:\ \varepsilon(\gamma\mid\Lambda_\gamma)\ge\varepsilon_d(\Lambda_\gamma)\bigr\}
=\bigl\{\gamma:\ \Psi_\tau(\gamma)\le d\bigr\},
\]

with \(\Psi_\tau:=\mathcal D(\gamma^\star(\tau)\mid\Lambda_\gamma)-
\mathcal D(\cdot\mid\Lambda_\gamma)\) and, by the complementarity identity,
\(\Psi_\tau=\frac12\|\cdot-\gamma^\star(\tau)\|^2_{P_{11}(\tau)}-\bar\Phi_\tau\), where
\(\Phi_\tau(\gamma)=\Phi_0(\gamma)+\frac{\tau}{2}(\gamma-\mu_0)^\top
\bar\Lambda(\gamma-\mu_0)\) is the potential of Cor. R2d-5 at prior scale \(\tau\).

Assume **(H2)**, **(H3a)**, **(H3b)** throughout 2B, so that Section 1B applies.

> **What 2B has to prove is now much less.** In the sublevel-set formulation the
> flat-limit work was dominated by *compactness* of the limiting set, which needed
> **(H3b)** through coercivity of \(\Phi_0\), and by a Hausdorff-convergence argument
> whose only purpose was to pass a minimum to the limit. Compactness is still proved —
> it is wanted downstream — but it is now the **prior-free** Lemma 2A.3-G, which needs
> only **(H3a)** and transfers to \(\tau=0\) verbatim, with a containment ellipsoid
> that does not depend on \(\tau\) at all. And no minimum has to be passed to the
> limit: only the **single scalar** \(\varepsilon(\gamma^\star(\tau)\mid\tau)\)
> (2B.4). **(H3b)** survives in 2B for one purpose only — placing the centre
> \(\gamma^{\star(0)}\) (Lemma 1B.2).

**Lemma 2B.1 (the limiting profile is well defined and finite).** At \(\Lambda_\gamma=0\)
the objects \(q_0(\cdot\mid\gamma)\) and \(q_{Q_0}\) are well defined and strictly
positive (Cor. 1B.3), so \(g(\cdot\mid\gamma;0)\) is finite and jointly continuous, and

1. \(g(\cdot\mid\gamma;0)\) is **strictly convex**, so its minimum is unique;
2. that minimum is **attained**, so
   \(\varepsilon(\gamma\mid0)\in(0,1]\) for every \(\gamma\) — in particular it is
   **finite and strictly positive**, i.e. \(\log\varepsilon(\gamma\mid0)>-\infty\).

Moreover \(\varepsilon(\gamma\mid\Lambda_\gamma)\to\varepsilon(\gamma\mid0)\) as
\(\tau\downarrow0\), uniformly for \(\gamma\) in compact sets.

*Proof.* Claims 1–2 are Lemmas 2A.1 and 2A.2 read at \(\tau=0\). This is legitimate
because their proofs consume only: \(\Lambda_Q\succ0\) — Lemma 1B.1 under **(H3a)**;
full row rank of the linear part
\(L_0=\Sigma^{\star(0)}[\mathcal W_1^\top P_b,\dots,\mathcal W_J^\top P_b]\) — again
**(H3a)**; full support and propriety of \(\pi(\beta\mid\gamma,y)\) — which is
\(\tau\)-**free**, since the population prior does not enter the \(\beta\)-block at all
(\(\pi(\beta_j\mid\gamma,y)\propto e^{\ell_j(\beta_j)}\phi(\beta_j;\mathcal W_j\gamma,\Psi)\));
and the bridge identity of Remark R2b-1, which is an identity for every family. **The
entire dependence of 2A on a proper prior is the existence of \(Q\)**, which 1B supplies
in the limit.

For the convergence, note \(q_{Q_\tau}\to q_{Q_0}\) uniformly on compacts (explicit
Gaussians with \(\gamma^\star(\tau)\to\gamma^{\star(0)}\), \(P_{11}(\tau)\to P_{11}^{\mathrm{RE}}\succ0\)),
and \(q_\tau(\gamma'\mid\gamma)\to q_0(\gamma'\mid\gamma)\) locally uniformly by
dominated convergence — the mixing law \(\pi(\beta\mid\gamma,y)\) is \(\tau\)-free,
\(m_\tau(\beta)\to m_0(\beta)\), \(\Sigma^\star(\tau)\to\Sigma^{\star(0)}\succ0\), and the
Gaussian integrand is bounded by \((2\pi)^{-q/2}\det(P_{11}(\tau))^{1/2}\le C\) for
small \(\tau\). Hence \(g(\cdot\mid\cdot;\tau)\to g(\cdot\mid\cdot;0)\) uniformly on
compacts. The uniform coercivity estimate of 2A.3 holds **uniformly in
\(\tau\in[0,\tau_0]\)** as well: its ingredients are \(\sigma_{\min}(L_\tau)\ge\sigma_0>0\)
(true for small \(\tau\) since \(L_\tau\to L_0\) of rank \(q\)), boundedness of
\(a_\tau-\gamma^\star(\tau)\) (Lemmas 1B.1–1B.2), and the same \(\tau\)-free
\(\pi(B_v\mid\gamma,y)\). So all minimisers lie in one ball \(|\gamma'|\le R\) for all
\(\tau\le\tau_0\) and \(\gamma\in K\), whence
\(\bigl|\log\varepsilon(\gamma\mid\tau)-\log\varepsilon(\gamma\mid0)\bigr|
\le\sup_{K\times B_R}|g(\cdot;\tau)-g(\cdot;0)|\to0\). \(\square\)

**Lemma 2B.2 (the limiting profile is continuous).**
\(\varepsilon(\cdot\mid0):\mathbb R^q\to(0,1]\) is continuous.

*Proof.* Lemma 2A.3 verbatim at \(\tau=0\): joint continuity of \(g(\cdot\mid\cdot;0)\)
gives upper semicontinuity, the uniform coercivity estimate (valid at \(\tau=0\), as
checked in 2B.1) confines the minimisation to a fixed compact ball locally uniformly in
\(\gamma\), and Berge concludes. (Alternatively: uniform-on-compacts convergence of
continuous functions, Lemma 2B.1, preserves continuity.) \(\square\)

**Lemma 2B.3 (the limiting set is convex and compact, uniformly in the prior).**
The set

\[
\widetilde C_d(0)\;=\;\bigl\{\gamma:\ \varepsilon(\gamma\mid0)\ge\varepsilon_d(0)\bigr\}
\;=\;\bigl\{\gamma:\ \Psi_0(\gamma)\le d\bigr\},
\qquad
\Psi_0=\tfrac12\|\cdot-\gamma^{\star(0)}\|^2_{P_{11}^{\mathrm{RE}}}-\bar\Phi_0,
\]

is **convex, compact and nonempty** with \(\gamma^{\star(0)}\) in its interior, and
\(\bigcup_{d>0}\widetilde C_d(0)=\mathbb R^q\). The containment bound of Lemma 2A.3-G
applies **with the same constants** at every \(\tau\), including \(\tau=0\):

\[
\widetilde C_d(\Lambda_\gamma)\;\subseteq\;
\bigl\{\|\gamma-\gamma^\star(\tau)\|^2_{S_\flat}\le2d\bigr\}
\quad\text{(bounded curvature)},
\qquad
\widetilde C_d(\Lambda_\gamma)\;\subseteq\;
\bigl\{|\gamma-\gamma^\star(\tau)|\le1+d/c_0\bigr\}
\quad\text{(in general)},
\]

with \(S_\flat=\sum_jH_j^\top P_b(P_b+G_j)^{-1}P_bH_j\succ0\) and \(c_0>0\) a
\(\tau\)-free constant. So there is \(\tau_0>0\) and a **single compact** \(K_d\)
containing \(\widetilde C_d(\Lambda_\gamma)\) for all \(\tau\le\tau_0\). Finally
\(\Psi_\tau\to\Psi_0\) uniformly on compacts and
\(\widetilde C_d(\Lambda_\gamma)\to\widetilde C_d(0)\) in **Hausdorff** distance.

*Proof.* **Closed:** \(\varepsilon(\cdot\mid0)\) is continuous (2B.2). **Convex and
compact:** Lemmas 2A.3-F and 2A.3-G read at \(\tau=0\). This is legitimate, and needs
no separate argument, because \(\nabla^2\Psi_\tau=\sum_jH_j^\top P_bV_j(\gamma)P_bH_j\)
**does not involve \(\Lambda_\gamma\)**: the mixing law
\(\pi(\beta_j\mid\gamma,y)\propto e^{\ell_j(\beta_j)}\phi(\beta_j;H_j\gamma,\Psi)\) is
prior-free, so the entire two-sided sandwich \(S_\flat\preceq\nabla^2\Psi_\tau\preceq
P_{11}^{\mathrm{RE}}\) holds verbatim and uniformly in \(\tau\in[0,\tau_0]\). Only the
centre moves, and \(\gamma^\star(\tau)\to\gamma^{\star(0)}\) by Lemma 1B.2, so the
ellipsoids all sit inside one compact \(K_d\). Without a curvature ceiling the same
conclusion follows from the ray bound: since \(\{\gamma^\star(\tau):\tau\le\tau_0\}\)
lies in a compact neighbourhood \(N\) of \(\gamma^{\star(0)}\), the constant
\(c_0:=\min_{|u|=1}\min_{\gamma_0\in N}\int_0^1u^\top S(\gamma_0+tu)u\,dt\) is positive
by the same continuity-and-compactness argument, and is \(\tau\)-free.
**Nonempty with interior point:** \(\Psi_0(\gamma^{\star(0)})=0<d\), and
\(\gamma^{\star(0)}\) exists and is finite by Lemma 1B.2 — this is where **(H3b)**
enters, and it is the only place in 2B that it does. **Exhausting:**
\(\varepsilon(\cdot\mid0)>0\) pointwise (2B.1), so \(\Psi_0<\infty\) everywhere.
**Hausdorff convergence:** \(\Psi_\tau=\log\varepsilon(\gamma^\star(\tau)\mid\tau)
-\log\varepsilon(\cdot\mid\tau)\), and \(\varepsilon(\cdot\mid\tau)\to\varepsilon(\cdot\mid0)\)
uniformly on compacts (2B.1) with \(\gamma^\star(\tau)\to\gamma^{\star(0)}\) (1B.2), so
\(\epsilon_\tau:=\sup_{K_d}|\Psi_\tau-\Psi_0|\to0\) and

\[
\{\Psi_0\le d-\epsilon_\tau\}\;\subseteq\;\widetilde C_d(\Lambda_\gamma)
\;\subseteq\;\{\Psi_0\le d+\epsilon_\tau\} .
\]

Both bracketing families converge to \(\widetilde C_d(0)\), by the standard fact that a
**convex** \(f\) attaining a value strictly below \(d\) satisfies
\(\{f\le d\}=\overline{\{f<d\}}\): here \(\Psi_0(\gamma^{\star(0)})=0<d\), so for any
\(\gamma\) with \(\Psi_0(\gamma)=d\) the segment points
\(\gamma_t=\gamma^{\star(0)}+t(\gamma-\gamma^{\star(0)})\), \(t<1\), have
\(\Psi_0(\gamma_t)\le td<d\) and converge to \(\gamma\). So no sublevel set is separated
from its interior and the sandwich closes. \(\square\)

> **What each hypothesis is doing, and how that has changed.**
> **(H3a)** remains indispensable and now does **more** work, not less: without it
> \(P_{11}^{\mathrm{RE}}\) is singular, so at \(\Lambda_\gamma=0\) there is no
> \(\Sigma^{\star(0)}\), no refresh density \(q_{Q_0}\) and no \(\varepsilon(\cdot\mid0)\)
> (Remark 1B.4) — *and* \(S_\flat\) is singular, so \(\Psi_0\) is flat along
> \(\bigcap_j\ker\mathcal W_j\) and \(\widetilde C_d(0)\) is an unbounded cylinder.
> **(H3b)** is used **only** to place the centre, via coercivity of \(\Phi_0\)
> (Lemma 1B.2). That is a real change: on the HPD route **(H3b)** was also what made
> the certified set bounded, since \(\{\bar\Phi_0\le d\}\) is unbounded whenever a
> separated group flattens \(\Phi_0\) at infinity. Here compactness comes from
> \(\nabla^2\Psi\succeq S_\flat\), which is prior-free and estimability-free, so a
> separated group costs the **centre** but not the **shape**. Even the residual use is
> soft: by the Centering remark any centre \(\gamma_c\) gives a valid profile and a
> valid compact superlevel family, at a fixed multiplicative cost.

**Lemma 2B.4 (the limiting constant, and continuity in the prior).**

\[
\min_{\widetilde C_d(0)}\varepsilon(\cdot\mid0)
\;=\;\varepsilon_d(0)\;:=\;e^{-d}\,\varepsilon\bigl(\gamma^{\star(0)}\bigm|0\bigr)\;\in\;(0,1),
\quad\text{attained on }\partial\widetilde C_d(0),
\]

and \(\varepsilon_d(\Lambda_\gamma)\to\varepsilon_d(0)\) as \(\tau\downarrow0\): **the
certified constant does not collapse in the limit**.

*Proof.* The minimum exists and is attained by Weierstrass — \(\varepsilon(\cdot\mid0)\)
is continuous (2B.2) and \(\widetilde C_d(0)\) is compact (2B.3) — and equals
\(\varepsilon_d(0)\) by the same two-sided argument as in 2A.4: \(\ge\) is the
definition of the superlevel set, and \(\le\) holds because coercivity of \(\Psi_0\)
plus the intermediate value theorem along a ray from \(\gamma^{\star(0)}\) puts a point
on the level \(\Psi_0=d\) exactly. Positivity is
\(\varepsilon(\gamma^{\star(0)}\mid0)>0\), which is 2B.1 at the single point
\(\gamma^{\star(0)}\). For the convergence there is now just one scalar to track:
\(\varepsilon_d(\tau)/\varepsilon_d(0)
=\varepsilon(\gamma^\star(\tau)\mid\tau)/\varepsilon(\gamma^{\star(0)}\mid0)\). Take a
compact \(K\) containing \(\{\gamma^\star(\tau):\tau\le\tau_0\}\cup\{\gamma^{\star(0)}\}\),
which exists since \(\gamma^\star(\tau)\to\gamma^{\star(0)}\) (1B.2). On \(K\),
\(\varepsilon(\cdot\mid\tau)\to\varepsilon(\cdot\mid0)\) uniformly (2B.1) and
\(\varepsilon(\cdot\mid0)\) is continuous (2B.2), so
\(\varepsilon(\gamma^\star(\tau)\mid\tau)\to\varepsilon(\gamma^{\star(0)}\mid0)\).
\(\square\)

> **Compare the work this used to take.** In both formulations the constant is a
> minimum over a compact set, produced by Weierstrass. What differs is what has to be
> passed to the limit. On the sublevel (HPD) route it is *the minimum itself*, which
> forces Hausdorff convergence of the sets and a common compact, and those in turn force
> coercivity of \(\Phi_0\) and hence **(H3b)**. Here the minimum has a known value,
> \(e^{-d}\varepsilon(\gamma^\star(\tau)\mid\tau)\), so what has to converge is a
> **scalar at a moving point** and the argument is "uniform convergence on a compact
> plus a convergent point". Hausdorff convergence and the common compact are still
> proved (2B.3) and still used — but they now come from the prior-free bound
> \(\nabla^2\Psi\succeq S_\flat\), so they no longer drag **(H3b)** in with them.

**Lemma 2B.5 (the minorization condition holds in the limit).** Under **(H2)**,
**(H3a)**, **(H3b)**, at \(\Lambda_\gamma=0\):

\[
q_0(\gamma'\mid\gamma)\;\ge\;\varepsilon_d(0)\,q_{Q_0}(\gamma')
\qquad\text{for all }\gamma\in\widetilde C_d(0),\ \gamma'\in\mathbb R^q ,
\]

so that \(\widetilde C_d(0)\) is a small set of the **unrestricted** limiting
\(\gamma\)-chain, \(q_0(\gamma,A)\ge\varepsilon_d(0)\,Q_0(A)\) for every Borel \(A\);
and, for the **restricted** limiting kernel \(q_0(\cdot,\cdot\mid\widetilde C_d(0))\)
of Lemma R2c, with \(Q_{0,\widetilde C_d}:=Q_0(\cdot\mid\widetilde C_d(0))\),

\[
q_0\bigl(\gamma,A\mid\widetilde C_d(0)\bigr)\;\ge\;\varepsilon^{(0)}\;Q_{0,\widetilde C_d}(A)
\qquad\forall\,\gamma\in\widetilde C_d(0),\ A\subseteq\widetilde C_d(0),
\qquad
\varepsilon^{(0)}:=\varepsilon_d(0)\,Q_0\bigl(\widetilde C_d(0)\bigr),
\]

with \(\varepsilon^{(0)}\in(0,1]\) and \(Q_0(\widetilde C_d(0))\ge\Pr(\chi^2_q\le2d)\).
So **Proposition R2 holds at \(\Lambda_\gamma=0\)**, with the flat-limit
\(\gamma\)-marginal in place of \(\pi_\gamma\) — note that the spectral truncation of
§7.2 Step 5 now does real work, since Cor. R2d-4(a) is vacuous only when
\(\Lambda_\gamma\succ0\).

*Proof.* Lemma 2A.5(a) and (b) verbatim, using 2B.4 for the constant and Cor. 1B.3 for
\(Q_0\); in particular the truncated-density step \(N(\beta)\le1\) is prior-free.
The refresh-mass bound is the inscribed-ellipsoid argument of the §7.3.2 preamble read
at \(\tau=0\): \(\gamma^{\star(0)}\) minimises \(\Phi_0\), so \(\bar\Phi_0\ge0\), and
the complementarity identity gives
\(\Psi_0(\gamma)\le\frac12(\gamma-\gamma^{\star(0)})^\top P_{11}^{\mathrm{RE}}
(\gamma-\gamma^{\star(0)})\). Hence \(\{(\gamma-\gamma^{\star(0)})^\top P_{11}^{\mathrm{RE}}
(\gamma-\gamma^{\star(0)})\le2d\}\subseteq\widetilde C_d(0)\), and that ellipsoid carries
\(Q_0\)-mass \(\Pr(\chi^2_q\le2d)\) since \(\Lambda_{Q_0}=P_{11}^{\mathrm{RE}}\).
\(\square\)

**Remark 2B.6 (what 2B adds beyond "the objects exist at \(\tau=0\)").** Two distinct
statements are proved, and it is worth separating them.

1. **Well-posedness at the limit.** Every object of Section 2A — \(g\),
   \(\varepsilon(\cdot)\), the certified set, the constant, the minorization — is
   defined and behaves identically at \(\Lambda_\gamma=0\), because the
   \(\varepsilon\)-machinery never touches the population prior. All the prior does is
   deliver \(Q\), and Section 1B delivers it.
2. **Continuity in the prior.** \(\varepsilon(\cdot\mid\Lambda_\gamma)\to
   \varepsilon(\cdot\mid0)\) uniformly on compacts (2B.1),
   \(\widetilde C_d(\Lambda_\gamma)\to\widetilde C_d(0)\) in **Hausdorff** distance,
   all inside one \(\tau\)-free compact (2B.3), and therefore
   \(\varepsilon_d(\Lambda_\gamma)\to\varepsilon_d(0)\) (2B.4).
   The certified constant does **not** degenerate as the prior is weakened.

Claim 2 is what an earlier draft flagged as unproven, on the grounds that the
continuity argument then in use rested on total-variation continuity of a *fixed*
pushforward measure and would have had to be redone for a varying one. The cone-cap
estimate of 2A.3 removes that obstruction: it is stated in \(\beta\)-space, where the
mixing law is prior-free, so adjoining \(\tau\) to the compact parameter set costs
nothing.

---


## Appendix C — Gaussian-closure closed forms


### C.1 Closed form (Gaussian closure only)

**Lemma R2b′ (closed form for the deficiency, under Gaussian closure).** Assume the
hypotheses of §7.3.2A, and in addition **Gaussian closure**:
\(q(\cdot\mid\gamma)=N(M(\gamma),\Sigma)\) with \(\Sigma\) **not depending on**
\(\gamma\) (the LMM case). Write \(\mathcal D(\gamma):=\log\varepsilon(\gamma)\) for
the **deficiency**. Then the program of Lemma 2A.2 can be solved explicitly:

\[
\boxed{\;
\mathcal D(\gamma)
=
\tfrac12\log\det\bigl(\Sigma^\star\Sigma^{-1}\bigr)
-\tfrac12\bigl\|\gamma^\star-M(\gamma)\bigr\|^2_{(\Sigma-\Sigma^\star)^{-1}},
\;}
\]

Since \(M(\gamma^\star)=\gamma^\star\), the value at the centre — the one number the
general construction leaves unevaluated — is

\[
\varepsilon(\gamma^\star)=\sqrt{\det\bigl(\Sigma^\star\Sigma^{-1}\bigr)}
=\det\bigl(I+\tilde J\bigr)^{-1/2}\;\ge\;(1+\kappa)^{-q/2},
\]

and the certified sets of §7.3.2 are **ellipsoids**,

\[
\widetilde C_d=\{\varepsilon\ge\varepsilon_d\}
=\bigl\{\gamma:\|\gamma^\star-M(\gamma)\|^2_{(\Sigma-\Sigma^\star)^{-1}}\le 2d\bigr\},
\qquad
\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)
=\sqrt{\det(\Sigma^\star\Sigma^{-1})}\,e^{-d}\in(0,1) .
\]

Both terms of \(\mathcal D\) are \(\le0\), as they must be:
\(\Sigma\succeq\Sigma^\star\) makes the first \(\le0\), and the second is a squared
norm.

> **Scope.** *Only* this lemma is closure-dependent. Lemmas 2A.1–2A.5 and 2B.1–2B.5
> hold for every GLM family and produce \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\)
> exactly; what closure adds is a **formula for the one remaining scalar**
> \(\varepsilon(\gamma^\star)\), together with the ellipsoidal description of the sets.
> Outside closure that scalar is the value of a strictly convex \(q\)-dimensional
> program at a single point (Remark R2b-0). *(The radius convention differs from
> earlier drafts by a factor of two: \(d\) is now measured in **nats of deficiency**,
> \(\Psi\le d\), rather than as a squared radius, which is why the constant reads
> \(e^{-d}\) and the ellipsoid \(\le2d\).)*

**Both constants are precision objects.** The proof derives
\(\det(\Sigma^\star\Sigma^{-1})=\det(\Lambda_q\Lambda_Q^{-1})\) and
\((\Sigma-\Sigma^\star)^{-1}=\Lambda_q(\Lambda_Q-\Lambda_q)^{-1}\Lambda_Q\), so the
boxed display is

\[
\mathcal D(\gamma)
=\tfrac12\log\det\bigl(\Lambda_q\Lambda_Q^{-1}\bigr)
-\tfrac12\bigl\|\gamma^\star-M(\gamma)\bigr\|^2_{\Lambda_q(\Lambda_Q-\Lambda_q)^{-1}\Lambda_Q}.
\]

Everything is driven by the **gap** \(\Lambda_Q-\Lambda_q\); the covariance symbols
are shorthand available only under closure.

**Both constants are precision objects.** The proof below derives
\(\det(\Sigma^\star\Sigma^{-1})=\det(\Lambda_q\Lambda_Q^{-1})\) and
\((\Sigma-\Sigma^\star)^{-1}=\Lambda_q(\Lambda_Q-\Lambda_q)^{-1}\Lambda_Q\), so the
boxed display is

\[
\mathcal D(\gamma)
=\tfrac12\log\det\bigl(\Lambda_q\Lambda_Q^{-1}\bigr)
-\tfrac12\bigl\|\gamma^\star-M(\gamma)\bigr\|^2_{\Lambda_q(\Lambda_Q-\Lambda_q)^{-1}\Lambda_Q}.
\]

Everything is driven by the **gap** \(\Lambda_Q-\Lambda_q\); the covariance symbols
are shorthand available only under closure.

*Proof.* Existence, uniqueness and attainment of the minimiser are Lemmas 2A.1–2A.2,
which already apply; what remains is to **evaluate** it. Gaussian closure makes
\(\Lambda_q=\Sigma^{-1}\)
**constant** in \(\gamma'\) — this, and only this, is what the following computation
needs. With \(a=M(\gamma)\) and \(u=\gamma^\star-a\),
\[
g(\gamma')
=\tfrac12\log\det(\Sigma^\star\Sigma^{-1}) + f(\gamma'),
\qquad
f(\gamma'):=\tfrac12(\gamma'-\gamma^\star)^\top\Lambda_Q(\gamma'-\gamma^\star)
-\tfrac12(\gamma'-a)^\top\Lambda_q(\gamma'-a).
\]
Minimising \(f\), with \(\Delta:=\Lambda_Q-\Lambda_q\succ0\), gives the stationary
point \(\gamma'_\star=\Delta^{-1}(\Lambda_Q\gamma^\star-\Lambda_q a)\),
so \(\gamma'_\star-\gamma^\star=\Delta^{-1}\Lambda_q u\) and \(\gamma'_\star-a=\Delta^{-1}\Lambda_Q u\),
whence
\[
f(\gamma'_\star)
=\tfrac12 u^\top\bigl[\Lambda_q\Delta^{-1}\Lambda_Q\Delta^{-1}\Lambda_q
-\Lambda_Q\Delta^{-1}\Lambda_q\Delta^{-1}\Lambda_Q\bigr]u
=-\tfrac12 u^\top\Lambda_q\Delta^{-1}\Lambda_Q\,u ,
\]
using \(\Lambda_Q=\Lambda_q+\Delta\) to expand both triple products. Finally
\[
\Lambda_q(\Lambda_Q-\Lambda_q)^{-1}\Lambda_Q
=\bigl(\Lambda_q^{-1}-\Lambda_Q^{-1}\bigr)^{-1}
=\bigl(\Sigma-\Sigma^\star\bigr)^{-1},
\]
since \(\Lambda_q^{-1}-\Lambda_Q^{-1}=\Lambda_q^{-1}(\Lambda_Q-\Lambda_q)\Lambda_Q^{-1}\).
Substituting gives the boxed \(\mathcal D\). Since \(\Sigma\succeq\Sigma^\star\) we have
\(\det(\Sigma^\star\Sigma^{-1})\le1\), so \(\varepsilon_d\in(0,1]\). \(\square\)

### C.2 Remarks on the construction

**Remark R2b-0 (what is actually available: existence, not a formula).** It is worth
being blunt about the logical status of each object, because the closed form is
seductive and applies only to LMMs.

| Object | General GLMM (§7.3.2) | Gaussian closure (Lemma R2b′) |
|---|---|---|
| \(\varepsilon(\gamma)=e^{\min_{\gamma'}g(\gamma'\mid\gamma)}\) | **exists** in \((0,1]\), **continuous**, **log-concave**; value of a strictly convex \(q\)-dimensional program (2A.1–2A.3, 2A.3-F) | closed form (boxed) |
| minimiser \(\gamma'_\star(\gamma)\) | exists, unique, confined a priori to \(|\gamma'|\le R_K\); no formula | \(2\gamma^\star-\gamma\) |
| \(\min_{\widetilde C_d}\varepsilon\) | attained, \(=e^{-d}\varepsilon(\gamma^\star)\) **exactly** (2A.4) | same, with \(\varepsilon(\gamma^\star)=\det(I+\tilde J)^{-1/2}\) |
| the certified set | \(\widetilde C_d=\{\varepsilon\ge\varepsilon_d\}\) — convex and compact, of no particular shape, but sandwiched between the \(P_{11}^{\mathrm{RE}}\)- and \(S_\flat\)-ellipsoids of radius \(\sqrt{2d}\) | exactly an ellipsoid in the metric \(P_{11}J\), which the sandwich then brackets |
| a numerical value for \(\varepsilon_d\) | one convex program, at the single point \(\gamma^\star\) | \(\ge(1+\kappa)^{-q/2}e^{-d}\) |

The point to keep hold of: **Proposition R2 never needs a numerical value.** The
Weierstrass minimum over \(\widetilde C_d\) is taken and is attained, but its value is
known in terms of \(d\) and the single scalar \(\varepsilon(\gamma^\star)\); the only
other thing proved about the model is that \(\varepsilon>0\) pointwise (2A.2), which is
what makes the family exhaust \(\mathbb R^q\). Everything quantitative beyond that is
optional refinement, and what refinement now costs is **one evaluation at one point**
rather than a minimisation whose answer varies with \(d\).

**Two routes, and why the superlevel one is used here.** Both families are indexed by
level sets of a convex function, both are convex (Lemma 2A.3-F and Cor. R2d-5(3)), both
are compact, and on **both** the minorization constant is an attained minimum of
\(\varepsilon\) obtained by continuity plus compactness. Nothing in that list
discriminates. What discriminates is **which** function indexes the family, and hence
whether the minimum has a **value**.

| | Superlevel route: \(\{\varepsilon\ge\varepsilon_d\}\) — **used in §7.3.2** | Sublevel (HPD) route: \(\{\bar\Phi\le d\}\) |
|---|---|---|
| indexing function | deficiency gap \(\Psi\), convex, \(\ge0\) | potential \(\bar\Phi\), convex, \(\ge0\) |
| the constant | \(\min_{\widetilde C_d}\varepsilon\), attained by Weierstrass | \(\min_{\widetilde C^{\mathrm{HPD}}_d}\varepsilon\), attained by Weierstrass |
| **its value** | \(=e^{-d}\varepsilon(\gamma^\star)\) **exactly** — the set is indexed by the function being minimised, so the minimum sits on the level that names the set | **unknown** — no relation to \(\varepsilon(\gamma^\star)\), and the whole \(d\)-dependence of the certificate is hidden inside it |
| what remains unknown | one scalar at one point, \(\varepsilon(\gamma^\star)\): Lemma 2A.2 | a minimum over a set, for each \(d\) |
| compactness of the set | holds, from \(\nabla^2\Psi\succeq S_\flat\succ0\): **(H3a)** only, **prior-free**, one lemma for both regimes (2A.3-G) | holds, from coercivity of \(\bar\Phi\): needs **(H3b)** in the flat limit, and must be proved separately per regime |
| posterior mass | \(\pi_\gamma(\widetilde C_d^{\,c})\le\Pr(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d\mid y)\) | \(\to0\) by exhaustion; envelope only via \(\Lambda_\gamma\), vacuous in the flat limit |
| mass optimality | not optimal — an HPD region of the same Lebesgue measure carries more | **optimal** by construction |
| size at a given constant | **maximal** — the largest set carrying \(\varepsilon_d\) | strictly smaller: \(\widetilde C^{\mathrm{HPD}}_d\subseteq\{\varepsilon\ge\min_{\widetilde C^{\mathrm{HPD}}_d}\varepsilon\}\) |

Two deciding considerations, then. **First, the value of the minimum.** Indexing by
\(\varepsilon\) itself makes the pair \((\widetilde C_d,\varepsilon_d)\) a fixed point of
the correspondence "set \(\mapsto\) its minimum" and "constant \(\mapsto\) its superlevel
set": \(\min_{\widetilde C_d}\varepsilon=\varepsilon_d\) and
\(\{\varepsilon\ge\varepsilon_d\}=\widetilde C_d\), each exactly. Indexing by \(\bar\Phi\)
breaks both identities, and the resulting minimum is a number one can neither evaluate
nor relate to \(d\). **Second, the price of compactness.** On the sublevel route it comes
from coercivity of \(\bar\Phi\), whose curvature
\(\nabla^2\bar\Phi=P_{11}-\nabla^2\Psi\) carries \(\Lambda_\gamma\) explicitly, so when
the prior is withdrawn the burden shifts to **(H3b)**. On the superlevel route it comes
from \(\nabla^2\Psi=\sum_jH_j^\top P_bV_j(\gamma)P_bH_j\succeq S_\flat\succ0\), which is
prior-free and estimability-free — it needs only that the random-effect conditional
covariances be nondegenerate, which is automatic, plus **(H3a)**. One lemma serves both
regimes and the containment ellipsoid does not move as \(\Lambda_\gamma\downarrow0\).

What is genuinely given up is mass *optimality*: at equal Lebesgue measure the HPD
region captures more posterior mass, so the superlevel route may need a slightly larger
\(d\) to meet the same \(\delta_2\). Since \(d\) enters the constant only through
\(e^{-d}\), that is a bounded and quantified loss (Remark R2a-2), not a structural one.

**Explicit constants require more than existence.** A *numerical* lower bound on
\(\mathcal D(\gamma)\) needs a **uniform lower** bound on the curvature — coercivity
alone is not enough. If \(\Lambda_Q-\Lambda_q(\gamma'\mid\gamma)\succeq cI\) for all
\(\gamma'\), then strong convexity gives, from any reference point \(\gamma'_0\),
\[
\mathcal D(\gamma)\;\ge\;g(\gamma'_0)-\tfrac1{2c}\bigl\|\nabla g(\gamma'_0)\bigr\|^2,
\qquad
\nabla g(\gamma'_0)=\Lambda_Q\bigl(\bar m(\gamma'_0)-\gamma^\star\bigr),
\]
with \(\bar m\) the bridge mean of Remark R2b-1 — computable by the same machinery as
Lemma R2d. This is the honest role of a uniform precision gap: **not** needed for
Proposition R2, **needed** for a certified number.

**Remark R2b-1 (the precision gap is an exact bridge covariance).** This is the
identity behind the whole construction, and it holds for **every** GLM family with no
Gaussian-closure assumption.

*Derivation.* \(q(\gamma'\mid\gamma)=\int\phi_q(\gamma';m(\beta),\Sigma^\star)\,
\pi(\beta\mid\gamma,y)\,d\beta\), so
\(\nabla_{\gamma'}\log q=-\Lambda_Q\bigl(\gamma'-\bar m(\gamma')\bigr)\) with
\(\bar m(\gamma')=E[m(\beta)\mid\gamma,\gamma',y]\) taken under the **bridge law**

\[
\pi(\beta\mid\gamma,\gamma',y)\;\propto\;\pi(\beta\mid\gamma,y)\,
\exp\Bigl\{m(\beta)^\top\Lambda_Q\gamma'-\tfrac12 m(\beta)^\top\Lambda_Q m(\beta)\Bigr\}
\]

— the law of the RE block given **both** endpoints of the step. This is again an
exponential tilt (natural parameter \(\theta=\Lambda_Q\gamma'\), statistic
\(m(\beta)\)), exactly as in Lemma R2d(a), so
\(\partial\bar m/\partial\theta=\mathrm{Cov}(m(\beta)\mid\gamma,\gamma',y)\) and
\(\partial\bar m/\partial\gamma'=\mathrm{Cov}(\cdot)\Lambda_Q\). Hence

\[
\boxed{\;
\Lambda_Q-\Lambda_q(\gamma'\mid\gamma)
\;=\;
\Lambda_Q\,\mathrm{Cov}\bigl(m(\beta)\mid\gamma,\gamma',y\bigr)\,\Lambda_Q
\;=:\;\Lambda_Q^{1/2}\,\tilde C(\gamma,\gamma')\,\Lambda_Q^{1/2},
\;}
\qquad
\tilde C:=\Lambda_Q^{1/2}\mathrm{Cov}(m\mid\gamma,\gamma',y)\Lambda_Q^{1/2}.
\]

Three consequences.

1. **Dominance is free, and no uniform version of it is needed.**
   \(\tilde C\succeq0\) always, so \(\Lambda_q\preceq\Lambda_Q=P_{11}\)
   **unconditionally** — this is why **(Q1)** never has to be verified. It is tempting
   to demand a *uniform* gap \(\tilde C\succeq cI\) to force coercivity of \(g\), but
   that is the wrong condition and it is not needed: Lemma 2A.2 gets coercivity from
   the **support** of \(P_\gamma\), not from curvature. The sharp
   statement is
   \[
   \varepsilon(\gamma)>0
   \iff
   \gamma^\star\in\mathrm{int}\,\mathrm{conv}\,\mathrm{supp}\,P_\gamma ,
   \]
   i.e. no hyperplane through \(\gamma^\star\) has all of the mass of \(m(\beta)\) on one side.
   Under **(H3a)** \(P_\gamma\) has full support \(\mathbb R^q\) and this is
   automatic — for every family, including those with unbounded log-likelihood
   curvature. The degenerate case of §7.2 is exactly the failure of this condition:
   holding \(\beta\) fixed makes \(P_\gamma\) a point mass at \(m(\beta)\), which lies
   in a supporting hyperplane, and \(\varepsilon\equiv0\).
   *(Consistency check: \(\Lambda_q\succeq0\) by Prékopa under **(H2)**, so the
   identity forces \(\tilde C\preceq I\) automatically.)*
2. **Gaussian closure: the covariance shadow, and the link to Lemma R2d.** When
   \(q(\cdot\mid\gamma)=N(M(\gamma),\Sigma)\), conditional independence of the
   \(\beta_j\) gives
   \[
   \Sigma-\Sigma^\star=\mathrm{Cov}\bigl(m(\beta)\mid\gamma,y\bigr)
   =P_{11}^{-1}\Bigl(\sum_j H_j^\top P_b V_j P_b H_j\Bigr)P_{11}^{-1}
   =P_{11}^{-1/2}\,\tilde J\,P_{11}^{-1/2},
   \]
   **the same \(\tilde J\) as Lemma R2d** — the excess *variance* is the Jacobian of
   the mean map. Equivalently \(\tilde C=\tilde J(I+\tilde J)^{-1}\): the **bridge**
   covariance is a shrunken version of the **marginal** one, as it must be. So
   \(\Sigma\succ\Sigma^\star\), \(\tilde J\succ0\) and \(\tilde C\succ0\) all say the
   same thing here — but only \(\tilde C\) survives outside closure.
3. **A uniform curvature sandwich — and what it does and does not buy.** Since
   \(m\) is **affine** in \(\beta\), \(\mathrm{Cov}(m\mid\cdot)=L\,\mathrm{Cov}
   (\beta\mid\cdot)\,L^\top\) with \(L=\partial m/\partial\beta
   =P_{11}^{-1}[\,H_1^\top P_b,\dots,H_J^\top P_b\,]\), and the bridge log-density has
   curvature \(\succeq D+L^\top\Lambda_Q L\), \(D:=\mathrm{blockdiag}_j\bigl((1+\epsilon_j)P_b\bigr)\)
   — the second term is the *extra* concentration from conditioning on \(\gamma'\).
   Brascamp–Lieb then gives \(\mathrm{Cov}(\beta\mid\text{bridge})\preceq
   (D+L^\top\Lambda_QL)^{-1}\), and Woodbury collapses the sandwich:
   \[
   \mathrm{Cov}(m\mid\gamma,\gamma',y)\;\preceq\;L\bigl(D+L^\top\Lambda_QL\bigr)^{-1}L^\top
   \;=\;\bigl(W^{-1}+\Lambda_Q\bigr)^{-1},
   \qquad
   W:=LD^{-1}L^\top=P_{11}^{-1/2}\tilde J_b P_{11}^{-1/2},
   \]
   \(\tilde J_b:=P_{11}^{-1/2}\bigl(\sum_j(1+\epsilon_j)^{-1}H_j^\top P_bH_j\bigr)P_{11}^{-1/2}
   \preceq\kappa I\). Conjugating by \(\Lambda_Q^{1/2}\),
   \[
   \tilde C\;\preceq\;\tilde J_b\bigl(I+\tilde J_b\bigr)^{-1}\;\preceq\;
   \frac{\kappa}{1+\kappa}I,
   \qquad\text{i.e.}\qquad
   \frac{\Lambda_Q}{1+\kappa}\;\preceq\;\Lambda_q(\gamma'\mid\gamma)\;\preceq\;\Lambda_Q
   \quad\text{for all }\gamma,\gamma' ,
   \]
   valid for **every** GLM family. Two careful readings:
   - **Under Gaussian closure this gives the constant.** There \(\Lambda_q\) is the
     constant \(\Sigma^{-1}\) and part (b) makes the leading factor exactly
     \(\sqrt{\det(\Lambda_q\Lambda_Q^{-1})}=\sqrt{\det(I-\tilde C)}
     =\det(I+\tilde J)^{-1/2}\), so
     \(\varepsilon(\gamma^\star)\ge(1+\kappa)^{-q/2}\) and hence
     \(\varepsilon_d\ge(1+\kappa)^{-q/2}e^{-d}\). The **same** spectral quantity that
     makes the chain contract lower-bounds the minorization constant, and the bound
     does not degenerate as \(\kappa\uparrow1\).
   - **Outside closure it does *not*.** There is no \(d\), and
     \(\mathcal D(\gamma)\) is not \(\tfrac12\log\det(I-\tilde C)\) minus a quadratic —
     that factorisation is a property of the Gaussian program. What the sandwich gives
     is an **upper** bound on \(\nabla^2 g=\Lambda_Q-\Lambda_q\), which is the wrong
     direction for bounding \(\min g\) below. A certified number needs the *opposite*
     inequality, \(\tilde C\succeq\tilde c I\) — a **lower** bound on the bridge
     covariance, i.e. Cramér–Rao from an **upper** bound on group curvature — fed into
     the strong-convexity estimate of Remark R2b-0.

The structural point behind **(Q1)** survives in either case: \(\tilde C\), the bridge
covariance, is **zero** exactly when \(\beta\) is held fixed. That is the degenerate
case of §7.2 — no integration, no precision gap, no minorization.

**Remark R2b-2 (LMM level sets are ellipsoids about \(\gamma^\star\)).** Under Gaussian
closure \(M\) is affine with \(M(\gamma)-\gamma^\star=J(\gamma-\gamma^\star)\) and \(\tilde J\)
constant. Using \(J=P_{11}^{-1/2}\tilde JP_{11}^{1/2}\) and
\((\Sigma-\Sigma^\star)^{-1}=P_{11}^{1/2}\tilde J^{-1}P_{11}^{1/2}\),
\[
\|\gamma^\star-M(\gamma)\|^2_{(\Sigma-\Sigma^\star)^{-1}}
=(\gamma-\gamma^\star)^\top P_{11}^{1/2}\tilde J P_{11}^{1/2}(\gamma-\gamma^\star),
\qquad\text{so}\qquad
C_d=\bigl\{\gamma:\|\gamma-\gamma^\star\|^2_{\Lambda_C}\le d\bigr\},
\quad
\Lambda_C:=P_{11}^{1/2}\tilde JP_{11}^{1/2}.
\]
So under closure the **certified sets** \(\widetilde C_d=\{\varepsilon\ge\varepsilon_d\}\)
are a nested family of concentric ellipsoids centred at \(\gamma^\star\), **bounded**
exactly when \(\tilde J\succ0\), i.e. under **(H3a)** — the closure instance of
Lemma 2A.3-G, where the two-sided sandwich degenerates to the single exact metric
\(\Lambda_C\). This is the level-set family of
`MINORIZATION_GAUSSIAN_REFRESH.md` §6.

**Three concentric ellipsoids, and an exact additive split.** It is worth writing the
metrics side by side, because they are easy to conflate:

| Family | Metric \(\Lambda\), with \(\{\|\gamma-\gamma^\star\|^2_\Lambda\le2d\}\) | Role |
|---|---|---|
| certified set \(\widetilde C_d=\{\Psi\le d\}\) | \(P_{11}J=P_{11}^{1/2}\tilde JP_{11}^{1/2}\) | carries the **constant** \(\varepsilon_d\) |
| HPD region \(\{\bar\Phi\le d\}\) | \(P_{11}(I-J)=P_{11}^{1/2}(I-\tilde J)P_{11}^{1/2}\) | carries the **mass**, optimally |
| refresh measure \(Q\) | \(\Lambda_Q=P_{11}\) | carries the refresh **mass** \(Q(\widetilde C_d)\) |

The three are not independent: the first two metrics **add to the third**, which is the
complementarity identity \(\Psi+\bar\Phi=\frac12\|\cdot-\gamma^\star\|^2_{P_{11}}\) in
closure notation. Since \(0\preceq\tilde J\preceq I\), each of the first two is
dominated by the third, so each of the first two sets *contains* the \(Q\)-ellipsoid of
squared radius \(2d\) — this is the inscribed-ellipsoid bound
\(Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)\) of §7.3.2, and it holds for the HPD family
for the same reason. The two families are complementary in a literal sense:
\(\tilde J\) large means a strong mean map, hence a wide certified set and a narrow
HPD region, and conversely. Outside closure the third is still an ellipsoid and the
first two are general convex bodies (Lemma 2A.3-F and Cor. R2d-5(3)); the additive
identity survives verbatim, since it was proved there without closure.

**Remark R2b-3 (GLMM: the kernel is a mixture, not a Gaussian).** For non-Gaussian
likelihoods \(q(\cdot\mid\gamma)=\int\phi_q(\cdot;m(\beta),\Sigma^\star)\,
\pi(\beta\mid\gamma,y)\,d\beta\) is a **Gaussian mixture**. What survives verbatim
and what needs work:

| Ingredient | GLMM status |
|---|---|
| \(q(\cdot\mid\gamma)\) log-concave | **Holds** — Prékopa: the joint in \((\gamma',\beta)\) is log-concave under **(H2)**, marginalise |
| **Precision dominance** \(\Lambda_q\preceq\Lambda_Q\) | **Holds exactly** — it is the identity of Remark R2b-1, valid for every family; nothing to verify |
| \(\varepsilon(\gamma)>0\) pointwise | **Holds** under **(H3a)** alone (Lemma 2A.2) — coercivity of the CGF comes from full support of \(m(\beta)\), not from curvature, so unbounded-curvature families (Poisson log) are fine |
| \(\varepsilon\) **continuous** (indeed \(C^1\)) | **Holds** (Lemma 2A.3) — locally uniform coercivity from the cone-cap estimate, then Berge/Danskin |
| \(\min_{\widetilde C_d}\varepsilon=\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)>0\), **attained** | **Holds** (Lemma 2A.4) — Weierstrass, from continuity (2A.3) + compactness (2A.3-G); the value is known because the set is indexed by \(\varepsilon\) |
| \(\mathcal D=\log\varepsilon\) **concave**, so \(\widetilde C_d\) **convex** | **Holds** (Lemma 2A.3-F) — \(\nabla^2\mathcal D=-P_{11}J\preceq0\), from the bridge law |
| certified set \(\widetilde C_d\) **closed and convex** | **Holds** — closed by continuity of \(\varepsilon\) (2A.3), convex by 2A.3-F |
| \(\widetilde C_d\) **compact** | **Holds** (Lemma 2A.3-G) — \(\nabla^2\Psi=\sum_jH_j^\top P_bV_jP_bH_j\succ0\) because the RE conditional covariances are nondegenerate; needs **(H3a)** only, and is prior-free |
| \(\widetilde C_d\) an **ellipsoid** | **Fails** — a general convex body; but it is **trapped between two concentric ellipsoids**, \(\{\|\cdot\|^2_{P_{11}^{\mathrm{RE}}}\le2d\}\subseteq\widetilde C_d\subseteq\{\|\cdot\|^2_{S_\flat}\le2d\}\) (2A.3-G), and the inscribed \(Q\)-ellipsoid gives \(Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)\) |
| complementarity \(\bar\Phi+\Psi=\frac12\|\cdot-\gamma^\star\|^2_{P_{11}}\) | **Holds** — Lemma 2A.3-F(3), proved without closure |
| Curvature sandwich \(\Lambda_Q/(1+\kappa)\preceq\Lambda_q\preceq\Lambda_Q\) | **Holds** — Remark R2b-1(3) used only affineness of \(m\), Brascamp–Lieb and Woodbury |
| \(d\)-dependence \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\) | **Holds** — exact for every family (Lemma 2A.4); \(d\) is the deficiency level, and the minimum lands on it by construction |
| Leading constant \(\varepsilon(\gamma^\star)\ge(1+\kappa)^{-q/2}\) | **Fails** — the sandwich of Remark R2b-1(3) bounds \(\nabla^2g\) from *above*, the wrong direction for bounding \(\min g\) below. A number needs Cramér–Rao plus Remark R2b-0 |
| A single \(\Sigma\) (covariance dominance) | **Meaningless** — the mixture has no one covariance; this is exactly why the hypothesis is stated on precisions |
| \(\Lambda_q\) **constant** | **Fails** — \(\Lambda_q(\gamma'\mid\gamma)=\Lambda_Q-\Lambda_Q\mathrm{Cov}(m\mid\gamma,\gamma',y)\Lambda_Q\) varies with the step |
| Closed-form \(\mathcal D\) | **Fails** — the boxed formula needs \(\Lambda_q\) constant; in general \(\mathcal D\) is the value of a smooth strictly convex \(q\)-dimensional program |

The genuine losses are losses of **closed form**, not of validity: the minorization,
the compact convex certified set with its explicit containment ellipsoid, and the
positivity of \(\varepsilon_d\) all survive intact — that is precisely why §7.3.2 was
stated without any formula. The
residue is a **single scalar**, \(\varepsilon(\gamma^\star)\): everything else in
\(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\) is explicit. To obtain a number one
either solves the \(q\)-dimensional convex program for \(\mathcal D(\gamma^\star)\)
numerically — a bounded search, by the radius \(R_K\) of Lemma 2A.3 — or substitutes
the uniform envelope
\(\bar\Lambda_q:=\Lambda_Q-\Lambda_Q\inf_{\gamma'}\mathrm{Cov}(m\mid\gamma,\gamma',y)\Lambda_Q
\preceq\Lambda_Q\) — note this needs an **infimum** of the bridge covariance
(Cramér–Rao), not the Brascamp–Lieb supremum of Remark R2b-1(3) — which restores
strong convexity and hence the boxed formula as a **lower** bound on
\(\mathcal D\). See §6.1 and **open item 2** of
`MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md`.

**Remark R2b-4 (certified constants).** Under Gaussian closure the bound
\(\varepsilon_d\ge(1+\kappa)^{-q/2}e^{-d}\) is explicit but dimension-sensitive
through \(q\) (the *population* dimension, typically small — not \(Jp_{\mathrm{re}}\)).
Outside closure there is no analogue without the extra curvature input above; the
statement of Proposition R2 does not depend on one. Sharper values come from
`MINORIZATION_GAUSSIAN_REFRESH.md` §5–§6 (exact inner problem) and
`MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` §6–§8 (symmetric GLM).

**Remark R2b-5 (division of labour between the A and B subsections).** Subsections 1A/2A
and 1B/2B prove the same things under different priors, and it is worth being exact
about *why* two versions are needed, because most of the machinery is prior-free.

| Ingredient | Really needs | Supplied in A by | Supplied in B by |
|---|---|---|---|
| \(\Lambda_Q=P_{11}\succ0\), so \(\Sigma^\star\) exists | \(\Lambda_\gamma\succ0\) **or** **(H3a)** | the prior (1A.1) | **(H3a)** (1B.1) |
| full support of \(m(\beta)\), hence \(\varepsilon>0\) | **(H3a)** alone | 2A.2 | 2B.1 — same argument |
| unique minimiser, continuity of \(\varepsilon\) | **(H2)**+**(H3a)** | 2A.1, 2A.3 | 2B.1–2B.2 — same argument |
| **a centre \(\gamma^\star\) satisfying (Q2)** | \(\Lambda_\gamma\succ0\) **or** **(H3b)** | strong convexity of \(\Phi\) (1A.2) | coercivity of \(\Phi_0\) (1B.2) |
| closed, convex certified set | **(H2)**+**(H3a)** — prior-free | 2A.3 + 2A.3-F | 2B.2 + 2B.3 |
| **compact** certified set | **(H2)**+**(H3a)** — prior-free | 2A.3-G | 2A.3-G verbatim (2B.3) |
| \(\min_{\widetilde C_d}\varepsilon\) attained (Weierstrass) | continuity + compactness | 2A.3 + 2A.3-G \(\Rightarrow\) 2A.4 | same, via 2B.2 + 2B.3 \(\Rightarrow\) 2B.4 |
| **(H1)** for Lemma R2a (exhaustion) | automatic, **or** **(H3b)** | proper prior | Lemma R2a‴ |

The compactness row is the one that changed. On the sublevel (HPD) route it read
"\(\Lambda_\gamma\succ0\) **or** **(H3a)**+**(H3b)**", and it was the single largest
consumer of **(H3b)** in the flat limit, because \(\{\bar\Phi_0\le d\}\) is bounded only
if \(\Phi_0\) is coercive. Certifying the superlevel set moves the burden onto
\(\nabla^2\Psi=\sum_jH_j^\top P_bV_j(\gamma)P_bH_j\), which does not involve
\(\Lambda_\gamma\) at all, so **one** lemma covers both columns and estimability drops
out. Note what did *not* change: the row beneath it. The constant is still an attained
minimum obtained by Weierstrass, so compactness remains a prerequisite for it — it has
simply become cheap, and prior-free, and its output has become computable. **(H3b)**
now enters Section 2 only through the centre \(\gamma^\star\) (via 1B.2), which the
Centering remark shows is a matter of a fixed multiplicative constant rather than of
validity.

Two consequences worth stating plainly:

1. **\(\varepsilon>0\) does not need (H3b).** By the Centering remark, *any* centre
   gives a valid minorization; the profile built from an arbitrary centre is strictly
   positive and continuous under (H2)+(H3a) whatever the prior. What (H3b) buys in the
   flat limit is the **canonical** centre — the one that sits inside every level set —
   and, separately, propriety of \(\pi_0\) so that Lemma R2a has a probability measure
   to work with.
2. **(H3a) is different in kind.** It is not a quality gate but a well-posedness
   condition: without it \(P_{11}^{\mathrm{RE}}\) is singular, so at
   \(\Lambda_\gamma=0\) there is no \(\Sigma^{\star(0)}\), no density \(q_{Q_0}\), and
   \(\varepsilon(\cdot\mid0)\) is not even defined. See Remark 1B.4.

##### Centering remark (why \(\gamma^\star=M(\gamma^\star)\) is the right centre)

Write \(\gamma_c\) for an **arbitrary candidate centre**, so that the refresh measure
under test is \(N(\gamma_c,\Sigma^\star)\); the claim of **(Q2)** is that
\(\gamma_c=\gamma^\star\) is the right choice. Any \(\gamma_c\) gives a **valid**
minorization — nothing in Lemmas 2A.1–2A.5 used **(Q2)**; they need only *some*
nondegenerate Gaussian to compare against. What the fixed point buys is that the
family is *usable*. Items 1 and 3 are now **theorems for every family**, because
Lemma 2A.3-F computes \(\nabla\mathcal D\) without any closure assumption; item 2 is
still stated in the closure notation of Lemma R2b′, where the penalty is visible as a
number.

1. **The refresh centre and the peak of the profile coincide only at the fixed
   point.** Building the profile from \(Q_c=N(\gamma_c,\Sigma^\star)\) and repeating
   Lemma 2A.3-F gives \(\nabla\mathcal D_c(\gamma)=P_{11}(\gamma_c-M(\gamma))\), so
   \(\mathcal D_c\) — still concave — peaks at \(M^{-1}(\gamma_c)\), not at
   \(\gamma_c\). The certified family \(\{\mathcal D_c\ge\mathcal D_c(M^{-1}(\gamma_c))-d\}\)
   is therefore a nest of convex sets centred at \(M^{-1}(\gamma_c)\) while the refresh
   measure sits at \(\gamma_c\): the two **drift apart** unless
   \(M(\gamma_c)=\gamma_c\), i.e. unless \(\gamma_c=\gamma^\star\). (The map \(M\) is
   injective, so the fixed point is the only such centre:
   \(\langle P_{11}^{1/2}(M(\gamma_1)-M(\gamma_2)),P_{11}^{1/2}(\gamma_1-\gamma_2)\rangle
   =\int_0^1\|\gamma_1-\gamma_2\|^2_{P_{11}^{1/2}\tilde J(\gamma_t)P_{11}^{1/2}}\,dt>0\)
   under **(H3a)**.) This was previously a closure heuristic; it is now a consequence
   of the envelope computation.
2. **Any other centre costs a fixed factor.** For an arbitrary \(\gamma_c\), let
   \(d_{\min}:=\inf_\gamma\|\gamma_c-M(\gamma)\|^2_{(\Sigma-\Sigma^\star)^{-1}}\), which is
   \(>0\) unless \(\gamma_c\) lies in the closure of the range of \(M\). Then every constant
   is scaled: \(\varepsilon_d\) at a given *effective* radius \(d-d_{\min}\) carries the
   extra factor \(e^{-d_{\min}/2}\). So a bad centre costs a **fixed multiplicative
   penalty, never validity** — the same robustness as Hypothesis 2 of the refresh
   note. Nor does a bad centre cost compactness: \(\nabla\mathcal D_c(\gamma)
   =P_{11}(\gamma_c-M(\gamma))\) differs from \(\nabla\mathcal D\) by a constant, so
   \(\nabla^2\Psi_c=\nabla^2\Psi\) and Lemma 2A.3-G applies verbatim with the same
   \(c\) and \(S_\flat\) — only the centre of the containment ellipsoid moves. This is
   why **(H3b)**, whose only remaining job in Section 2 is to supply
   \(\gamma^{\star(0)}\) in the flat limit (Lemma 1B.2, used in 2B.3), is a quality
   condition rather than a well-posedness one.
3. **It is what makes the refresh mass bound work.** Once the sampler is restricted to
   \(\widetilde C_d\) (§7.2 Step 7), the surviving constant is
   \(\varepsilon_d\,Q(\widetilde C_d)\). The inscribed-ellipsoid bound
   \(Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)\) is the complementarity identity plus
   \(\bar\Phi\ge0\), and **both** halves of that require \(\gamma^\star\): the identity
   is normalised at the peak of \(\mathcal D\), and \(\bar\Phi\ge0\) says the same point
   minimises \(\Phi\). Off-centre the two normalisation points separate, the quadratic
   picks up a linear term, the inscribed ellipsoid shrinks and its \(\chi^2\) becomes
   noncentral: mass is lost on **both** factors at once. The coincidence of the two
   points is Cor. R2d-5 — \(\gamma^\star\) is at once the fixed point of \(M\) and the
   mode of \(\pi(\gamma\mid y)\) — and it is the reason the two families of
   Remark R2b-2 are concentric.


## Appendix D — Exhaustion and the restricted kernel


### D.1 Lemma R2a (exhaustion: the certified sets carry all the mass)

**Statement.** Under **(H1)** and the hypotheses of §7.3.2A (or, at
\(\Lambda_\gamma=0\), of §7.3.2B), the family
\(\widetilde C_d=\{\varepsilon\ge e^{-d}\varepsilon(\gamma^\star)\}=\{\Psi\le d\}\) is
**nondecreasing** in \(d\) with \(\bigcup_{d>0}\widetilde C_d=\mathbb R^q\).
Consequently \(\pi_\gamma(\widetilde C_d^{\,c})\downarrow0\) as \(d\uparrow\infty\), and
for every \(\delta_2>0\) there is \(d(\delta_2)<\infty\) with
\(\pi_\gamma\bigl(\widetilde C_{d(\delta_2)}^{\,c}\bigr)<\delta_2\). Explicitly, and
with an envelope that needs no further hypothesis,

\[
\pi_\gamma\bigl(\widetilde C_d^{\,c}\bigr)
\;=\;\Pr\bigl(\Psi(\gamma)>d\bigm| y\bigr)
\;\le\;
\Pr\bigl(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d\bigm| y\bigr).
\]

*Proof.* Monotonicity and the union: \(\varepsilon(\gamma)>0\) for every \(\gamma\)
(Lemma 2A.2, or 2B.1 at \(\tau=0\)), so \(\Psi(\gamma)<\infty\) and every \(\gamma\)
lies in \(\widetilde C_d\) once \(d\ge\Psi(\gamma)\). Since \(\pi_\gamma\) is a
probability measure (**H1**) and \(\widetilde C_d^{\,c}\) is nonincreasing with empty
intersection, continuity from above gives \(\pi_\gamma(\widetilde C_d^{\,c})\to0\). The
envelope is \(\Psi\le\frac12\|\cdot-\gamma^\star\|^2_{P_{11}}\), which by the
complementarity identity is the statement \(\bar\Phi\ge0\), i.e. that \(\gamma^\star\)
minimises \(\Phi\). \(\square\)

> **Pointwise positivity is doing the work, and it is enough.** Exhaustion is exactly
> the statement that the certified family is not stuck: it rests on \(\varepsilon>0\)
> at *every* \(\gamma\) (Lemma 2A.2) — the one genuinely model-dependent input of
> §7.3.2, coming from full support of \(m(\beta)\) under **(H3a)**. Note what is *not*
> needed **for this lemma specifically**: no tightness argument, no appeal to
> compactness of \(\widetilde C_d\) (which does hold, by 2A.3-G, and which Lemma 2A.4
> does use), and no separate tail estimate for \(\bar\Phi\). The set family and the
> constant come from the same function, so
> "the constant holds on the set" (2A.4) and "the sets exhaust the space" (here) are
> two readings of the same two facts, \(\varepsilon\le\varepsilon(\gamma^\star)\) and
> \(\varepsilon>0\).

**Remark R2a-1 (what changed).** Two earlier versions of this lemma have been
discarded. The first invoked tightness of a Borel measure on a Polish space to produce
*some* compact \(K\); that is true but useless, since \(K\) had no relation to the set
on which minorization holds, so the two had to be intersected and the minorization
re-derived. The second indexed the family by the posterior potential \(\bar\Phi\) —
HPD regions — which made exhaustion trivial but left the constant as a Weierstrass
minimum of **unknown value**, and required the set to be proved compact separately for
each prior regime, in the flat limit at the cost of **(H3b)**. Here **one** nested
family carries **both** properties, indexed by a single knob \(d\). The Weierstrass
minimum is still taken (Lemma 2A.4) and compactness is still needed for it — but the
minimum now *evaluates* to \(e^{-d}\varepsilon(\gamma^\star)\), and compactness comes
from the prior-free Lemma 2A.3-G.

**Remark R2a-2 (the Pareto tension, quantitatively).** The two requirements pull in
opposite directions along the *same* axis \(d\), and both sides are now explicit:

\[
\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\quad\text{(decreasing in }d),
\qquad
\pi_\gamma(\widetilde C_d^{\,c})\quad\text{(decreasing in }d),
\qquad
Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)\ \uparrow1 ,
\]

the certified Doeblin constant being \(\varepsilon=\varepsilon_d\,Q(\widetilde C_d)\).

The parametrisation makes the mass side a **one-dimensional** tail — that of the scalar
\(\Psi(\gamma)\) under the posterior — and Lemma R2a bounds it by the tail of
\(\|\gamma-\gamma^\star\|^2_{P_{11}}\), the *same* statistic whose \(Q\)-tail gives the
refresh mass. Under **(H2)**, \(\pi_\gamma\) is log-concave (Cor. R2d-5(3)), so with

\[
\rho\;:=\;\lambda_{\max}\bigl(\Sigma_\pi^{1/2}P_{11}\Sigma_\pi^{1/2}\bigr)
\;=\;\lambda_{\max}\bigl(\Sigma_\pi P_{11}\bigr)\;\ge\;1
\]

the envelope is sub-exponential, \(\pi_\gamma(\widetilde C_d^{\,c})\lesssim e^{-d/\rho}\).
Choosing \(d=\rho\log(1/\delta_2)\) to meet the mass budget then gives, for **every**
GLM family and with no closure assumption,

\[
\boxed{\;
\varepsilon_d\;\gtrsim\;\varepsilon(\gamma^\star)\,\delta_2^{\,\rho} .
\;}
\]

So the minorization constant degrades only **polynomially** in the tail budget
\(\delta_2\), with exponent \(\rho\) — the mismatch between the posterior spread and
the refresh metric \(P_{11}\). This is the §4.6.3 tension of
`LOGIT_STATIC_TAIL_CERTIFICATION.md` made explicit, and it says the certification is
cheap when \(\rho=O(1)\).

Two things are worth noting about this display, because they are what the superlevel
parametrisation bought. First, the metric is \(P_{11}\), not \(\Lambda_\gamma\): the
old HPD envelope came from \(\bar\Phi\ge\frac12\|\cdot\|^2_{\Lambda_\gamma}\) and
therefore degraded with the prior and was **vacuous at \(\Lambda_\gamma=0\)**, whereas
\(P_{11}\succeq P_{11}^{\mathrm{RE}}\succ0\) survives the flat limit under **(H3a)**.
Second, the relation holds for every family, not only under closure, because
\(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\) is exact by construction; the only
unevaluated quantities left are the scalar \(\varepsilon(\gamma^\star)\) and the spread
constant \(\rho\).

**Under Gaussian closure everything closes.** There \(\Sigma_\pi=(P_{11}(I-J))^{-1}\),
so \(\rho=\lambda_{\max}((I-\tilde J)^{-1})=(1-\kappa)^{-1}\) and
\(\|\gamma-\gamma^\star\|^2_{P_{11}}\) is a weighted \(\chi^2_q\) with weights
\((1-\lambda_i(\tilde J))^{-1}\); also \(\varepsilon(\gamma^\star)=\det(I+\tilde J)^{-1/2}
\ge(1+\kappa)^{-q/2}\) by Lemma R2b′. Hence taking
\(d=\tfrac12(1-\kappa)^{-1}\chi^2_{q,1-\delta_2}\),

\[
\pi_\gamma\bigl(\widetilde C_d^{\,c}\bigr)\le\delta_2,
\qquad
\varepsilon_d\;\ge\;(1+\kappa)^{-q/2}
\exp\Bigl\{-\frac{\chi^2_{q,1-\delta_2}}{2(1-\kappa)}\Bigr\},
\qquad
Q\bigl(\widetilde C_d\bigr)\ge\Pr\bigl(\chi^2_q\le2d\bigr),
\]

a complete certification in terms of \(q\), \(\kappa\) and \(\delta_2\) alone. Outside
closure the same three displays hold with \(\varepsilon(\gamma^\star)\) and \(\rho\)
left as the two numbers to supply.

### D.2 Lemma R2c (restricted kernel and stationarity)

**Statement.** Let \(C\subseteq\mathcal X\) be measurable with \(\pi(C)>0\). Define
the **truncated target** \(\pi(A\mid C):=\pi(A\cap C)/\pi(C)\) and let
\(P_1(\cdot\mid C)\) be the two-block Gibbs kernel **of that target**: one sweep
drawing each block from the conditional of \(\pi(\cdot\mid C)\) given the other
block, i.e. from the full conditional **restricted to the corresponding slice of
\(C\)** and renormalized. Then

1. \(P_1(x,C\mid C)=1\) for every \(x\in C\), and
2. \(\pi(\cdot\mid C)\) is invariant for \(P_1(\cdot\mid C)\).

*Proof.* Both blocks are drawn inside the relevant slice of \(C\) by construction,
giving (1); the slices are nonempty for \(x\in C\) since \(x\) itself lies in them.
For (2), each block update is a Gibbs update for the target \(\pi(\cdot\mid C)\) —
it draws from the exact conditional of that target — and a Gibbs update leaves its
own target invariant; a composition of invariant kernels is invariant. \(\square\)

**Remark R2c-1 (why this definition).** Stationarity is immediate **because** the
kernel is defined as Gibbs *for the truncated target*, not as "run the original
sweep and repair excursions". Constructions that run the unrestricted sweep and then
project or resample on exit generally do **not** preserve \(\pi(\cdot\mid C)\) and
were the source of the earlier open item here.

**Remark R2c-2 (relation to the sampler actually run).** The package runs the
**unrestricted** two-block Gibbs kernel, which is *not* \(P_1(\cdot\mid C)\).
Proposition R2 therefore certifies a **companion** chain, not the production one;
closing that gap is Approach B (escape control for the full kernel) — see §7.5.

**Remark R2c-3 (the \(\gamma\)-chain instance).** Apply Lemma R2c with the
**cylinder** \(\tilde C=C\times\mathbb R^{Jp_{\mathrm{re}}}\), \(C\subseteq\mathbb R^q\).
Three things follow.

1. The \(\beta\)-slice of \(\tilde C\) is all of \(\mathbb R^{Jp_{\mathrm{re}}}\), so
   **Block 1 is unrestricted**: the sweep is an exact draw
   \(\beta\sim\pi(\beta\mid\gamma,y)\) followed by
   \(\gamma'\sim N(m(\beta),\Sigma^\star)\) **truncated to \(C\)**. Only the cheap,
   low-dimensional block is truncated.
2. \(\pi(\cdot\mid\tilde C)\) has \(\gamma\)-marginal \(\pi_\gamma(\cdot\mid C)\),
   and the \(\gamma\)-component of the restricted sweep is Markov with kernel
   \(q(\cdot\mid\cdot\,;C)\); Lemma R2c gives its invariance.
3. Truncating the \(\gamma\)-draw divides by a slice mass
   \(\int_C\phi_q(\gamma';m(\beta),\Sigma^\star)\,d\gamma'\le1\), so the restricted
   density **dominates** the unrestricted one on \(C\). The minorization of §7.2
   Step 6 is therefore inherited without change.


## Appendix E — Spectral truncation: existence of $\lambda^\star$

Supplies the spectral-cutoff ingredient of Appendix F: a rate threshold
$\lambda^\star<1$ whose exceedance set has small posterior mass. Nothing here is
used outside Appendix F except through the four results below (E.1–E.4).


Supplies **Step 5** of the §7.2 proof: a rate threshold \(\lambda^\star<1\) whose
exceedance set has small posterior mass. Nothing here is used elsewhere in the note
except through the three results below.

The appendix works with the **Laplace** rate \(T(\beta)=\lambda_{\max}(A^{(0)}(\beta))\)
on the \(\beta\)-block. The \(\gamma\)-chain of §7.2 uses the **exact-covariance**
analogue \(\kappa(\gamma)=\lambda_{\max}(\tilde J(\gamma))\) of **Cor. R2d-4**, which
inherits step (iv) below verbatim; Cor. R2d-2 records that the two objects differ
only by Laplace-vs-exact conditional covariance.

| Result | Statement | Used in |
|---|---|---|
| **A.1 Lemma R2a′** | \(A^{(0)}\) well defined; \(0\le T\le1\); \(\exists\lambda^\star\in(0,1)\) with \(\pi(T>\lambda^\star)<\delta_1\); weak-prior dominance | §7.2 Step 5, via **Cor. R2d-4** |
| **A.2 Corollary R2a″** | \(T\) continuous and an **attained** eigenvalue; limits \(\Lambda_\gamma\downarrow0\), \(\delta_1\downarrow0\) | §7.2 Step 5 |
| **A.3 Lemma R2a‴** | flat-limit posterior \(\pi_0\) **proper** — supplies **(H1)** at \(\Lambda_\gamma=0\) | §7.1 (H1) |
| **A.4 Example** | separation: what degrades, and what actually fails | — |


### E.0 Definition (rate matrix and rate sublevel set)

Throughout, **\(A\)** means the two-block Gibbs **rate matrix** as a function of
the RE block \(\beta\) (through \(P_{22}(\beta)\), \(S(\beta)\), and the
\(\beta\)-dependent blocks of \(P_{12},P_{21}\)), in the sense of
`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §I.3–§I.4 (A.1 below for
notation):

\[
\boxed{
A(\beta)
\;:=\;
P_{11}^{-1/2}\,
P_{12}\,
P_{22}(\beta)^{-1}\,
P_{21}\,
P_{11}^{-1/2},
\qquad
\lambda^\star(\beta)
\;:=\;
\lambda_{\max}\bigl(A(\beta)\bigr).
}
\]

*(Population-prior scaling enters \(P_{11}\); with fixed hyperparameters the
state dependence of \(A\) is through \(\beta\). If \(\gamma\) also enters
\(P_{12},P_{11}\) at the current sweep, write \(A(\gamma,\beta)\) — the spectral
event below is still \(\lambda_{\max}(A) > \lambda^\star\).)*

Given \(\delta_1>0\), Lemma R2a′ (A.1, claim 3) supplies \(\lambda^\star\in(0,1)\)
**sufficiently large** (close to \(1\) from below) with

\[
\pi\bigl(\lambda_{\max}\bigl(A^{(0)}(\beta)\bigr) > \lambda^\star\bigr) \;<\; \delta_1 .
\]

Define the **rate sublevel set**

\[
\boxed{
E_{\lambda^\star}
\;:=\;
\bigl\{
(\gamma,\beta) :
\lambda_{\max}\bigl(A^{(0)}(\beta)\bigr) \le \lambda^\star
\bigr\}.
}
\]

*(With \(\lambda^\star<1\), this is a **genuine** truncation — not the whole space.)*

### E.1 Lemma R2a′ (existence of \(\lambda^\star\); weak population prior)

**Setup (notation).** Use the package two-block precision blocks
(`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §I.3–§I.4):
\[
P_{11}(\Lambda_\gamma)
=
\underbrace{\sum_{j=1}^{J} H_j^\top P_b H_j}_{P_{11}^{\mathrm{RE}}}
+
\underbrace{\mathrm{blockdiag}_k(V_k^{-1})}_{\text{population prior on }\gamma},
\]
\[
S(\beta)
=
P_{12}(\beta)\,P_{22}(\beta)^{-1}\,P_{21}(\beta)
=
\sum_{j=1}^{J} H_j^\top P_b\,\bigl(P_b+\mathcal{G}_j(\beta)\bigr)^{-1} P_b\, H_j,
\]
\[
A(\beta;\Lambda_\gamma)
=
P_{11}(\Lambda_\gamma)^{-1/2}\,S(\beta)\,P_{11}(\Lambda_\gamma)^{-1/2},
\qquad
A^{(0)}(\beta)
=
{P_{11}^{\mathrm{RE}}}^{-1/2}\,S(\beta)\,{P_{11}^{\mathrm{RE}}}^{-1/2}.
\]
The **flat / weak limit** is \(\Lambda_\gamma \downarrow 0\) (equivalently
\(V_k^{-1}\to 0\) Loewner), so \(P_{11}(\Lambda_\gamma)\to P_{11}^{\mathrm{RE}}\)
with \(S(\beta)\) unchanged. Write
\[
T(\beta) := \lambda_{\max}\bigl(A^{(0)}(\beta)\bigr).
\]

**Claim.** Under **(H1)–(H3)**:

1. \(P_{11}^{\mathrm{RE}} \succ 0\) and \(A^{(0)}(\beta)\) is well defined — **(H3a)**,
   step (i).
2. \(0 \le T(\beta) \le 1\) for all \(\beta\) — step (iii); matrix algebra only.
3. **Spectral tail** (step (iv)): for every \(\delta_1 > 0\),
   \(\exists\,\lambda^\star \in (0,1)\) with \(\pi(T>\lambda^\star) < \delta_1\).
   Existence uses **(H3a)** + group full rank + positive GLM weights + **(H1)**;
   **(H2)** and **(H3b)** give the **quantitative** version (a \(\lambda^\star\)
   bounded away from \(1\) by a computable margin).
4. Weak-prior dominance — step (v); then combine with (3).
5. **Attainment, including the weak-prior limit** — Corollary R2a″: \(T\) is an
   attained eigenvalue; \(\lambda_{C_d}=\sup_{C_d}T\le\lambda^\star<1\) on the
   certified set, and is a maximum whenever that set is bounded; as
   \(\Lambda_\gamma\downarrow0\) the eigenvalue
   \(\lambda_{\max}(A(\beta;\Lambda_\gamma))\) increases to \(T(\beta)\), stays
   attained, and obeys the **prior-free** ceiling \(1/(1+\epsilon(\beta))<1\); and
   \(\delta_1=0\) is admissible exactly when the family's weights are bounded below.

**(H3b) in one line.** For each group \(j\), the group GLM \(\ell_j(\beta_j)\) has a
**finite** maximizer (finite `glm` solution). The package tests this as
`groupef.estimable` (`check_identifiability.R`).

**What (H3b) does and does not do for claim (3).** Pointwise, \(T(\beta)<1\) at every
finite \(\beta\) as soon as each \(\mathcal G_j(\beta_j)\succ0\), which follows from
full column rank of \(Z_j\) and strict positivity of the GLM weights — no estimability
needed. Since \(\pi\) is a probability measure on a finite-dimensional space **(H1)**,
the decreasing events \(\{T>t\}\) shrink to \(\varnothing\) and \(F(t)\to0\). What
**(H3b)** controls is the **rate**: without a finite group-wise MLE the posterior is
pushed up a likelihood recession direction where \(\mathcal G_j\to0\), forcing
\(\lambda^\star\) so close to \(1\) that the certificate is useless (and, in the flat-\(\gamma\)
limit, destroying **(H1)** itself).

**Where each hypothesis enters.**

| Hypothesis | Step | Role |
|---|---|---|
| **(H3a)** `popef.rank_ok` | **(i), (iv)-2** | \(P_{11}^{\mathrm{RE}}\succ 0\); \(A^{(0)}\), \(T\) defined |
| full-rank \(Z_j\) (`groupef.rank`) | **(iv)-3** | \(\mathcal G_j(\beta_j)\succ0\) \(\Rightarrow T(\beta)<1\) |
| **(H1)** proper \(\pi\) | **(iv)-4** | \(\{T>t\}\downarrow\varnothing\Rightarrow F(t)\to0\) |
| **(H2)** log-concavity | **(iv)** *(quantitative)* | Sub-Gaussian \(\beta_j\) tails \(\Rightarrow\) explicit \(\lambda^\star(\delta_1)\) |
| **(H3b)** finite group-wise MLE | **(ii)**, **(iv)** *(quantitative)* | Data-determined center; keeps \(1-\lambda^\star\) usable |
| *(none)* | **(iii), (v)** | Matrix bounds; no \(\pi\) |

*Proof.*

**(i) Identifiability \(\Rightarrow P_{11}^{\mathrm{RE}} \succ 0\).**
`check_identifiability()` requires `popef.rank_ok`: the hyper-design \(W\)
has **full column rank \(q\)** on estimable population coordinates. With
\(P_b \succ 0\), the Gram matrix
\[
P_{11}^{\mathrm{RE}}
=
\sum_{j=1}^{J} H_j^\top P_b H_j
\]
is positive definite (standard block-Gram argument: full rank of the stacked
\(H_j\) on the estimable subspace). Hence \(P_{11}^{\mathrm{RE}}\) is invertible
and \(A^{(0)}\) is defined in the flat-\(\gamma\) limit \(\Lambda_\gamma=0\).

**(ii) \(S(\beta)\) is well defined unconditionally.** GLM weights are nonnegative, so
\(\mathcal{G}_j(\beta_j)\succeq 0\) and \(B_j(\beta_j)=P_b+\mathcal{G}_j(\beta_j)\succ 0\)
whenever \(P_b\succ 0\) — no rank or estimability condition is needed here
(`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §I.5: \(P_{22}\succ0\) for every
likelihood). Full column rank of \(Z_j\) (`groupef.rank`) is used later, in **(iv)-3**,
to upgrade \(\mathcal G_j\succeq 0\) to \(\mathcal G_j\succ 0\); finite group-wise MLE
**(H3b)** is used only for the quantitative margin in **(iv)**.

**(iii) Spectral bound \(T(\beta)\le 1\) (multivariate) — no \(\pi\).**
For each group \(j\) and each \(u\in\mathbb R^q\), with \(t:=P_b^{1/2}H_j u\),
\[
u^\top H_j^\top P_b B_j^{-1} P_b H_j u
=
t^\top B_j^{-1} t
\;\le\;
t^\top P_b^{-1} t
=
u^\top H_j^\top P_b H_j u,
\]
because \(B_j\succeq P_b\) implies \(B_j^{-1}\preceq P_b^{-1}\). Summing over
\(j\),
\[
u^\top S(\beta)\, u \;\le\; u^\top P_{11}^{\mathrm{RE}} u
\qquad \forall u\in\mathbb R^q.
\]
By the generalized Rayleigh quotient identity
(`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §I.4),
\[
T(\beta)
=
\lambda_{\max}\bigl(A^{(0)}(\beta)\bigr)
=
\max_{u\neq 0}\frac{u^\top S(\beta)\, u}{u^\top P_{11}^{\mathrm{RE}} u}
\;\le\; 1.
\]
Nonnegativity is immediate since \(S(\beta)\succeq 0\).

**(iv) Spectral tail — proof of claim (3).**

*What we prove.* Define the **tail function**
\[
F(t)\;:=\;\pi\bigl(T(\beta) > t\bigr),
\qquad t\in[0,1).
\]
Then \(F\) is **non-increasing** in \(t\). Claim (3) is equivalent to
\[
\lim_{t\uparrow 1} F(t) = 0.
\]
Given that limit, for any \(\delta_1>0\) choose
\[
\lambda^\star
\;:=\;
\inf\bigl\{ t \in [0,1) : F(t) \le \delta_1 \bigr\}.
\]
Then \(\pi(T>\lambda^\star)=F(\lambda^\star)\le \delta_1\) and \(\lambda^\star<1\).
So the entire problem is: **as \(\lambda^\star\) is raised toward \(1\), show the
posterior mass on \(\{T>\lambda^\star\}\) vanishes.**

*Two separate questions.* (a) **Existence** of \(\lambda^\star<1\) — settled by
Steps 1–5 below from rank/positivity and **(H1)**. (b) **Size** of \(1-\lambda^\star\)
— this is where estimability **(H3b)** and log-concavity **(H2)** matter; see the
quantitative refinement after Step 5 and the degenerate example below.

---

**Step 1 — rewrite \(T\) as a Fisher-information gap (uses (iii) only).**

**Group data precision — general definition (every family, every \(p_{\mathrm{re}}\)).**
\(\mathcal{G}_j(\beta_j)\) denotes the **group data-precision Gram matrix** of
`CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §I.2 (written \(\mathcal W_j\) there;
renamed here because \(\mathcal W_j\) is the hyper design in §7.1):
\[
\mathcal{G}_j(\beta_j)
\;=\;
Z_j^\top W_j(\beta_j)\, Z_j
\;=\;
\sum_{i=1}^{n_j} w_{j,i}(\eta_{j,i})\, z_{j,i} z_{j,i}^\top
\;\in\; \mathbb R^{p_{\mathrm{re}}\times p_{\mathrm{re}}},
\qquad
w_{j,i}
=
\frac{\mathrm{wt}_{j,i}\,\bigl[\mu'(\eta_{j,i})\bigr]^2}{V(\mu_{j,i})\,\phi},
\]
\(\eta_{j,i}=\mathrm{offset}_{j,i}+z_{j,i}^\top\beta_j\). This is the **generic**
IRLS / Fisher weight (Gaussian, binomial with any link, Poisson, Gamma, …) and
\(p_{\mathrm{re}}\ge 1\) is arbitrary: nothing in (iv) is link-specific or univariate.
Recall \(B_j(\beta_j) = P_b + \mathcal{G}_j(\beta_j)\).

From (iii), \(S(\beta)\preceq P_{11}^{\mathrm{RE}}\). Write the **information gap**
\[
\Delta(\beta)
\;:=\;
P_{11}^{\mathrm{RE}} - S(\beta)
\;=\;
\sum_{j=1}^{J} H_j^\top K_j(\beta_j)\, H_j
\;\succeq\; 0,
\qquad
K_j \;:=\; P_b - P_b B_j^{-1} P_b .
\]
By the **parallel-sum (harmonic-mean) identity** applied to \(B_j = P_b+\mathcal G_j\),
\[
\boxed{
K_j(\beta_j)
\;=\;
P_b\bigl(P_b+\mathcal{G}_j\bigr)^{-1}\mathcal{G}_j
\;=\;
\bigl(P_b^{-1} + \mathcal{G}_j^{-1}\bigr)^{-1}
\;=:\;
P_b : \mathcal{G}_j ,
}
\]
symmetric PSD in all cases (the inverse forms when \(\mathcal G_j\succ0\)), and
**Loewner-monotone increasing** in \(\mathcal G_j\).
Set \(M(\beta):= {P_{11}^{\mathrm{RE}}}^{-1/2}\,\Delta(\beta)\,{P_{11}^{\mathrm{RE}}}^{-1/2}\).
For every \(\beta\),
\[
T(\beta)
=
\max_{u\neq 0}\frac{u^\top S(\beta)\,u}{u^\top P_{11}^{\mathrm{RE}} u}
=
1 - \min_{u\neq 0}\frac{u^\top \Delta(\beta)\,u}{u^\top P_{11}^{\mathrm{RE}} u}
=
1 - \lambda_{\min}\bigl(M(\beta)\bigr).
\]
Hence, for \(\eta\in(0,1]\),
\[
\boxed{
\{\,T(\beta) > 1-\eta\,\}
\;=\;
\{\,\lambda_{\min}(M(\beta)) < \eta\,\}.
}
\]
*Contribution:* pure matrix algebra — no \(\pi\), no estimability yet. **Raising
\(\lambda^\star\)** is the same as **shrinking the Fisher-gap threshold \(\eta\).**

---

**Step 2 — an explicit bound: relative group information controls \(T\) (uses (H3a)).**

Define the **relative group information** (dimensionless; invariant under
reparametrization of \(\beta_j\))
\[
\epsilon_j(\beta_j)
\;:=\;
\lambda_{\min}\Bigl(P_b^{-1/2}\,\mathcal{G}_j(\beta_j)\,P_b^{-1/2}\Bigr),
\qquad
\epsilon(\beta) \;:=\; \min_{1\le j\le J}\epsilon_j(\beta_j)
\]
— how much data precision group \(j\) carries **relative to** its RE prior precision,
in the worst direction.

If \(\mathcal G_j \succeq \epsilon P_b\), monotonicity of the parallel sum gives
\[
K_j
=
\bigl(P_b^{-1}+\mathcal G_j^{-1}\bigr)^{-1}
\;\succeq\;
\bigl(P_b^{-1}+\epsilon^{-1}P_b^{-1}\bigr)^{-1}
=
\frac{\epsilon}{1+\epsilon}\,P_b .
\]
Summing over \(j\) and using \(P_{11}^{\mathrm{RE}}=\sum_j H_j^\top P_b H_j\),
\[
\Delta(\beta) \;\succeq\; \frac{\epsilon(\beta)}{1+\epsilon(\beta)}\,P_{11}^{\mathrm{RE}}
\quad\Longrightarrow\quad
\lambda_{\min}\bigl(M(\beta)\bigr) \;\ge\; \frac{\epsilon(\beta)}{1+\epsilon(\beta)},
\]
so by Step 1
\[
\boxed{
T(\beta) \;\le\; \frac{1}{1+\epsilon(\beta)},
\qquad\text{hence}\qquad
\{\,T > 1-\eta\,\}
\;\subseteq\;
\bigcup_{j=1}^{J}
\Bigl\{\,\epsilon_j(\beta_j) < \tfrac{\eta}{1-\eta}\,\Bigr\}.
}
\]
*Contribution:* **(H3a)** gives \(P_{11}^{\mathrm{RE}}\succ 0\), so \(M\) and \(T\)
are defined. Still no \(\pi\).

*Why not decompose \(M=\sum_j M_j\) directly?* A per-group term
\(M_j = {P_{11}^{\mathrm{RE}}}^{-1/2}H_j^\top K_j H_j {P_{11}^{\mathrm{RE}}}^{-1/2}\)
has rank \(\le p_{\mathrm{re}} < q\) in general, so \(\lambda_{\min}(M_j)=0\) for every
\(\beta_j\) whenever \(H_j\) has a null space — a "one group at a time" spectral
argument on \(M_j\) is **vacuous**. Going through \(\epsilon_j\) (a statement about
\(\mathcal G_j\) versus \(P_b\) in the \(p_{\mathrm{re}}\)-dimensional group space)
avoids the hyper-design null space entirely.

---

**Step 3 — \(T(\beta)<1\) at every finite \(\beta\) (uses group full rank + weight positivity).**

Two family-generic facts:

1. **Positive weights at finite linear predictors.** For every package family the
   weight \(w(\eta)=\mathrm{wt}\,[\mu'(\eta)]^2/(V(\mu)\phi)\) is **continuous and
   strictly positive** at finite \(\eta\): the link is a diffeomorphism onto the
   interior of the mean domain, where \(\mu'(\eta)\neq 0\), \(0<V(\mu)<\infty\),
   \(0<\phi<\infty\). Weight degeneracy \(w\to 0\) can occur **only** as
   \(|\eta|\to\infty\) — and for some families (Gaussian; Gamma with log link) not
   even then.
2. **Full column rank \(Z_j\)** on estimable coordinates (`groupef.rank`, part of
   **(H3)**): for finite \(\beta_j\) and any \(v\ne 0\),
   \[
   v^\top \mathcal{G}_j(\beta_j)\, v
   =\sum_{i=1}^{n_j} w_{j,i}(\eta_{j,i})\,(z_{j,i}^\top v)^2 \;>\; 0 .
   \]

Hence \(\mathcal{G}_j(\beta_j)\succ 0\) and \(\epsilon_j(\beta_j)>0\) at every finite
\(\beta_j\), so \(\epsilon(\beta)>0\) and Step 2 yields the **strict** bound
\[
\boxed{\,
T(\beta)\;\le\;\frac{1}{1+\epsilon(\beta)}\;<\;1
\qquad\text{for every }\beta\in\mathbb R^{Jp_{\mathrm{re}}} .}
\]
Both \(\beta\mapsto\epsilon(\beta)\) and \(\beta\mapsto T(\beta)\) are continuous, so
on every compact \(K\): \(\inf_K\epsilon>0\) and \(\sup_K T<1\). **The ceiling
\(T=1\) is approached only as \(\|\beta\|\to\infty\).**

*Contribution:* group full rank + positive GLM weights. No estimability, no
log-concavity, no link-specific structure.

---

**Step 4 — the tail vanishes (uses (H1)).**

The events \(\{T>t\}\) **decrease** as \(t\uparrow 1\), and by Step 3
\[
\bigcap_{t<1}\{\,T>t\,\}
\;=\;
\{\,T=1\,\}
\;=\;
\varnothing
\]
as a subset of \(\mathcal X=\mathbb R^q\times\mathbb R^{Jp_{\mathrm{re}}}\): no
**finite** state attains the ceiling.

**(H1)** says \(\pi\) is a probability measure on \(\mathcal X\). Continuity of a
finite measure from above gives
\[
\lim_{t\uparrow 1}F(t)
=
\pi\Bigl(\bigcap_{t<1}\{T>t\}\Bigr)
=
\pi(\varnothing)
=
0 .
\]
Equivalently via Step 2:
\(F(1-\eta)\le\sum_j \pi\bigl(\epsilon_j(\beta_j)<\eta/(1-\eta)\bigr)\), and each
event decreases to \(\varnothing\) as \(\eta\downarrow 0\).

*Contribution:* **(H1)** is the entire probabilistic input — a proper posterior
puts no mass "at infinity", which is the only place \(T\) reaches \(1\).

---

**Step 5 — the quantile.**

For any \(\delta_1>0\), Step 4 makes \(\{t\in[0,1):F(t)\le\delta_1\}\) nonempty, so
\[
\lambda^\star \;:=\; \inf\bigl\{t\in[0,1):F(t)\le\delta_1\bigr\}
\;\in\;(0,1),
\qquad
\pi(T>\lambda^\star)\;\le\;\delta_1 .
\]
\(\square\)

---

**Summary — what each assumption buys (existence).**

| Step | Fact used | Hypothesis |
|---|---|---|
| 1 | \(T=1-\lambda_{\min}(M)\); \(K_j=P_b:\mathcal G_j\) | (iii) / matrix algebra |
| 2 | \(T\le 1/(1+\epsilon)\) with \(\epsilon=\min_j\lambda_{\min}(P_b^{-1/2}\mathcal G_jP_b^{-1/2})\) | **(H3a)** |
| 3 | \(\mathcal G_j(\beta_j)\succ0\) at finite \(\beta_j\) \(\Rightarrow T<1\) | full-rank \(Z_j\) + \(w(\eta)>0\) |
| 4 | \(\{T>t\}\downarrow\varnothing\), continuity from above | **(H1)** |
| 5 | Quantile | none |

**Honest accounting.** Existence of \(\lambda^\star\in(0,1)\) needs **only** (H3a),
group full rank, positive GLM weights, and **(H1)**. It does **not** need **(H2)**,
**(H3b)**, symmetry, drift, or a "\(\pi(T=1)=0\)" axiom. What **(H2)** and **(H3b)**
buy is the **quantitative** statement below — a \(\lambda^\star\) bounded away from
\(1\) by a computable amount, which is what makes \(\varepsilon_d\) usable.

---

**Quantitative refinement — where (H2) and (H3b) actually matter.**

Existence alone is weak: if \(\lambda^\star = 1-10^{-12}\), the certified set
\(E_{\lambda^\star}\) is useless in practice (Proposition R2 gives a valid but
astronomically slow bound, and \(\varepsilon_d\) degrades). A **usable**
\(\lambda^\star\) needs a rate for \(F\). Chain the same steps quantitatively:

\[
F(1-\eta)
\;\le\;
\sum_{j=1}^{J}
\pi\Bigl(\epsilon_j(\beta_j)<\tfrac{\eta}{1-\eta}\Bigr)
\;\le\;
\sum_{j=1}^{J}
\pi\bigl(\|\beta_j - m_j\| > R_j(\eta)\bigr)
\;\le\;
\sum_{j=1}^{J} C_j\, e^{-\alpha_j R_j(\eta)^2}.
\]

- **Middle inequality — (H3b).** Requires a **deterministic** implication
  "\(\epsilon_j\) small \(\Rightarrow\) \(\beta_j\) far from a fixed center \(m_j\)",
  with \(R_j(\eta)\to\infty\) as \(\eta\downarrow0\). Finite group-wise MLE
  \(\hat\beta_j\) supplies a data-determined center \(m_j=\hat\beta_j\) at which
  \(\epsilon_j\) is bounded below by an **observable** quantity (the fitted
  IRLS weights), so \(R_j(\eta)\) is computable from the fit. Without a finite MLE,
  the natural center is prior-driven and \(R_j\) degrades to prior scale.
- **Right inequality — (H2).** Log-concavity of \(\pi(\beta_j\mid\gamma,y_j)\) plus
  the Gaussian RE penalty gives sub-Gaussian concentration (Brascamp–Lieb), i.e.
  explicit \(C_j,\alpha_j\).

So: **(H1) + (H3a) + full rank ⇒ \(\lambda^\star\) exists; (H2) + (H3b) ⇒
\(\lambda^\star\) is computable and not absurdly close to \(1\).**

**(v) Weak population prior (multivariate).**
Since \(P_{11}(\Lambda_\gamma)=P_{11}^{\mathrm{RE}}+\Lambda_\gamma\succeq
P_{11}^{\mathrm{RE}}\),
\[
\frac{u^\top S(\beta)\, u}{u^\top P_{11}(\Lambda_\gamma)\, u}
\;\le\;
\frac{u^\top S(\beta)\, u}{u^\top P_{11}^{\mathrm{RE}}\, u}
\qquad \forall u\neq 0,
\]
so \(\lambda_{\max}(A(\beta;\Lambda_\gamma))\le T(\beta)\), with equality at
\(\Lambda_\gamma=0\). \(\square\)

### E.2 Corollary R2a″ (the rate is an attained eigenvalue, including in the limits)

Lemma R2a′ produces a threshold \(\lambda^\star\) per \(\delta_1\). For the
certificate to deliver a **geometric** rate one also needs the rate to be an
**attained** eigenvalue of a well-defined matrix — on the certified set, in the
weak-prior limit, and (where possible) in the \(\delta_1\downarrow0\) limit.

**(a) Attainment on the certified set (Weierstrass).**
\(A^{(0)}(\beta)\) is a real symmetric \(q\times q\) matrix for every finite
\(\beta\), so by the spectral theorem \(T(\beta)=\lambda_{\max}(A^{(0)}(\beta))\) is
attained by a maximizing unit vector \(u(\beta)\) — it is a genuine eigenvalue, not a
supremum. Moreover \(\beta\mapsto T(\beta)\) is continuous (weights are continuous;
eigenvalues are continuous functions of a symmetric matrix). Hence on the certified set
\(C_d=\widetilde C_d\cap E_{\lambda^\star}\) of §7.2, which is **compact** by
Lemma 2A.3-G,
\[
\boxed{
\lambda_{C_d} \;:=\; \max_{(\gamma,\beta)\in C_d} T(\beta)
\quad\text{is attained, and}\quad
\lambda_{C_d} \;\le\; \lambda^\star \;<\; 1 .
}
\]
So the restricted kernel has a **uniform, attained** contraction factor
\(\lambda_{C_d}\), not merely a pointwise bound — this is what the geometric factor
\((1-\varepsilon_d)^n\) of Proposition R1 is paired with.

**(b) Weak population prior \(\Lambda_\gamma\downarrow0\): a valid eigenvalue survives the limit.**

This is the limit in which a certified rate is most at risk, since \(P_{11}\) loses the
prior term that anchors the \(\gamma\)-block. Three statements, in order.

**(b1) The limit exists and is attained.** By (v),
\(\Lambda_\gamma\mapsto\lambda_{\max}(A(\beta;\Lambda_\gamma))\) is non-increasing in
\(\Lambda_\gamma\) (Loewner) and bounded above by \(T(\beta)\), so the limit exists and
\[
\lambda_{\max}\bigl(A(\beta;\Lambda_\gamma)\bigr)\;\uparrow\;T(\beta)
\qquad\text{as }\Lambda_\gamma\downarrow0 .
\]
The limiting value is an **attained eigenvalue** of the symmetric matrix
\(A^{(0)}(\beta) = {P_{11}^{\mathrm{RE}}}^{-1/2}S(\beta){P_{11}^{\mathrm{RE}}}^{-1/2}\),
which is well defined at \(\Lambda_\gamma=0\) because **(H3a)** gives
\(P_{11}^{\mathrm{RE}}\succ0\) without any prior contribution (step (i)). No
degeneracy, no \(0/0\): the flat limit is a legitimate eigenvalue problem, and
\(T=\sup_{\Lambda_\gamma\succeq0}\lambda_{\max}(A(\cdot;\Lambda_\gamma))\) is attained
at \(\Lambda_\gamma=0\).

**(b2) The ceiling is prior-free, hence uniform in \(\Lambda_\gamma\).** The bound of
Step 2 involves only the **group data precision** and the **RE prior precision**,
\[
\epsilon_j(\beta_j)=\lambda_{\min}\bigl(P_b^{-1/2}\mathcal{G}_j(\beta_j)P_b^{-1/2}\bigr),
\]
and contains **no \(\Lambda_\gamma\)** whatsoever. Therefore, for **every**
\(\Lambda_\gamma\succeq0\) including the flat limit,
\[
\boxed{
\lambda_{\max}\bigl(A(\beta;\Lambda_\gamma)\bigr)
\;\le\;
T(\beta)
\;\le\;
\frac{1}{1+\epsilon(\beta)}
\;<\;1 ,
}
\]
and on the certified set,
\(\sup_{\Lambda_\gamma\succeq0}\sup_{C_d}\lambda_{\max}(A(\cdot;\Lambda_\gamma))
=\lambda_{C_d}\le\lambda^\star<1\), both suprema attained — the outer at
\(\Lambda_\gamma=0\), the inner because \(C_d\) is **compact** (Lemma 2A.3-G). *(Only
the bound \(\le\lambda^\star\) is used in §7.2 Step 5; attainment is what makes
\(\lambda_{C_d}\) a genuine rate rather than an infimum of valid rates.)*
**The spectral half of the certificate does not weaken as the population prior
weakens** — the rate ceiling is carried entirely by \((\mathcal G_j, P_b)\), i.e. by
the data and the RE-level prior. Weakening \(\Lambda_\gamma\) moves
\(\lambda_{\max}(A)\) up to \(T\), never past it.

**(b3) What *does* depend on the limit: the measure, not the matrix.** \(T\) is a
fixed function of \(\beta\), but the tail mass is taken under
\(\pi_{\Lambda_\gamma}\), which changes with the prior. Claim (3) applies to each
proper \(\pi_{\Lambda_\gamma}\) separately; to keep **one** \(\lambda^\star\) valid
along the whole limit one needs **uniform tightness**:
\[
\lim_{R\to\infty}\ \sup_{0\preceq\Lambda_\gamma\preceq\Lambda_0}
\pi_{\Lambda_\gamma}\bigl(\|\beta\|>R\bigr)\;=\;0
\quad\Longrightarrow\quad
\exists\,\lambda^\star<1:\ \sup_{0\preceq\Lambda_\gamma\preceq\Lambda_0}
\pi_{\Lambda_\gamma}(T>\lambda^\star)<\delta_1 .
\]
Sufficient condition: the flat-limit posterior \(\pi_0\) is itself **proper**
(**(H1)** at \(\Lambda_\gamma=0\)); then \(\pi_{\Lambda_\gamma}\Rightarrow\pi_0\)
weakly and the family is tight, so \(\lambda^\star\) can be chosen once and reused.
If \(\pi_0\) is **improper** — the separation case below, where mass escapes along
\(\beta\sim\gamma\to\infty\) — tightness fails, \(\lambda^\star(\delta_1)\uparrow1\)
as \(\Lambda_\gamma\downarrow0\), and no uniform certificate exists.
**Propriety of \(\pi_0\) is not an extra axiom**: Lemma R2a‴ below derives it from
log-concavity + full rank + estimability, i.e. from **(H2)+(H3a)+(H3b)** and
`groupef.rank`.

**Summary of (b).** As the population prior gets weak: the eigenvalue itself stays
well defined, attained, and bounded by the prior-free ceiling
\(1/(1+\epsilon(\beta))<1\); only the **mass** statement can fail, and it fails
exactly when the flat-limit posterior is improper. *(Separately, weakening
\(\Lambda_\gamma\) degrades the minorization constant \(\varepsilon_d\) — the
Pareto tension of §7.5 — which is a different mechanism from the rate ceiling.)*

**(c) Limit \(\delta_1\downarrow0\): when truncation can be removed.**
\(\lambda^\star(\delta_1)\) is non-increasing in \(\delta_1\), so
\[
\lambda^\star(\delta_1)\;\uparrow\;\lambda_{\mathrm{ess}}
\;:=\;\operatorname*{ess\,sup}_{\pi} T
\;=\;\inf\{t:\pi(T>t)=0\}
\qquad\text{as }\delta_1\downarrow0 .
\]
Two regimes, decided by whether the family's weights are **bounded below**:

| Weight behaviour | Families | \(\epsilon(\beta)\) | \(\lambda_{\mathrm{ess}}\) | Truncation |
|---|---|---|---|---|
| \(w(\eta)\equiv \mathrm{wt}/\phi\) bounded below | Gaussian (identity), Gamma (log link) | \(\ge\epsilon_{\min}>0\) **uniformly in \(\beta\)** | \(\le \frac{1}{1+\epsilon_{\min}}<1\) | **not needed**: \(\delta_1=0\), \(E_{\lambda^\star}=\mathcal X\) |
| \(w(\eta)\to0\) as \(|\eta|\to\infty\) | binomial (logit, probit, cloglog), Poisson (log), Gamma (inverse) | \(>0\) pointwise, \(\inf_\beta\epsilon=0\) | \(=1\) when \(\operatorname{supp}\pi=\mathcal X\) | **necessary**: \(\delta_1>0\) |

In the first regime the rate is a **global** attained eigenvalue,
\(\sup_\beta T(\beta)\le 1/(1+\epsilon_{\min})<1\) with
\(\epsilon_{\min}=\min_j\lambda_{\min}\bigl(P_b^{-1/2}Z_j^\top Z_j P_b^{-1/2}\bigr)/\phi\),
and Proposition R2 degenerates to the familiar unrestricted LMM statement (no
truncation, \(\pi(C_d^c)=0\)). In the second regime \(\lambda_{\mathrm{ess}}=1\), so
\(\delta_1>0\) is **not** an artefact of the proof: no single \(\lambda^\star<1\)
covers the full support, and the restriction to \(C_d\) is essential. \(\square\)

### E.3 Lemma R2a‴ (propriety of the flat-limit posterior: log-concavity + full rank + estimability)

This supplies **(H1)** at \(\Lambda_\gamma=0\), which is the one ingredient
Corollary R2a″ (b3) cannot get from the matrix algebra.

**Hypotheses.**

| | Condition | Package gate |
|---|---|---|
| **(P1)** | \(\ell_j\) **concave** in \(\beta_j\), with weights \(w_{j,i}(\eta)>0\) at finite \(\eta\) | **(H2)** (GLM families with log-concave likelihood) |
| **(P2)** | \(Z_j\) **full column rank** for each \(j\) | `groupef.rank` |
| **(P3)** | stacked \(\bigl[H_1^\top\cdots H_J^\top\bigr]^\top\) has **rank \(q\)** \(\iff P_{11}^{\mathrm{RE}}\succ0\) | **(H3a)** `popef.rank_ok` |
| **(P4)** | **finite group-wise MLE** \(\hat\beta_j\) for each \(j\) | **(H3b)** `groupef.estimable` |

**Statement.** Under (P1)–(P4), the flat-\(\gamma\) posterior
\[
\pi_0(\gamma,\beta\mid y)
\;\propto\;
\exp\bigl(g(\gamma,\beta)\bigr),
\qquad
g(\gamma,\beta)
:=
\sum_{j=1}^{J}\ell_j(\beta_j)
-\tfrac12\sum_{j=1}^{J}(\beta_j-H_j\gamma)^\top P_b\,(\beta_j-H_j\gamma),
\]
is **proper**: \(\int_{\mathbb R^q}\int_{\mathbb R^{Jp_{\mathrm{re}}}}e^{g}\,d\beta\,d\gamma<\infty\).
Hence **(H1)** holds at \(\Lambda_\gamma=0\), \(\{\pi_{\Lambda_\gamma}\}_{\Lambda_\gamma\preceq\Lambda_0}\)
is tight, and Corollary R2a″ (b3) gives a **single \(\lambda^\star<1\) valid uniformly
as the population prior goes flat**.

*Proof.*

**1. \(g\) is concave and bounded above.** Each \(\ell_j\) is concave by (P1) and the
quadratic term is concave, so \(g\) is concave on
\(\mathbb R^q\times\mathbb R^{Jp_{\mathrm{re}}}\). By (P4) each
\(\ell_j\le\ell_j(\hat\beta_j)<\infty\) and the quadratic term is \(\le0\), so
\(g\le\sum_j\ell_j(\hat\beta_j)<\infty\).

**2. Each \(\ell_j\) has compact superlevel sets.** By (P2) the map
\(\beta_j\mapsto\eta_j=Z_j\beta_j\) is injective, and with \(w_{j,i}>0\) (P1) the
log-likelihood is **strictly** concave in \(\beta_j\); by (P4) its maximum is attained
at a finite \(\hat\beta_j\). A strictly concave function with an attained maximum
decays **at least linearly** along every ray: for each unit \(v\),
\[
\ell_{j}^{\infty}(v)
:=\lim_{t\to\infty}\frac{\ell_j(\hat\beta_j+tv)-\ell_j(\hat\beta_j)}{t}
\;<\;0,
\qquad
c_j:=-\max_{\|v\|=1}\ell_j^{\infty}(v)>0
\]
(the max is attained since \(\ell_j^\infty\) is concave, positively homogeneous, u.s.c.
on the compact sphere). This is exactly "finite, unique `glm()` solution with positive
IRLS weights" — i.e. no separation and no flat direction.

**3. Every direction recedes.** Take a direction \((d_\gamma,d_\beta)\ne0\).

- *Case A: \(d_{\beta,j}\ne H_j d_\gamma\) for some \(j\).* Then along the ray the
  quadratic term \(-\tfrac12\sum_j\|\beta_j+td_{\beta,j}-H_j(\gamma+td_\gamma)\|^2_{P_b}\)
  decreases **quadratically** in \(t\) (since \(P_b\succ0\)), while
  \(\sum_j\ell_j\) is bounded above by Step 1. So \(g\to-\infty\) quadratically.
- *Case B: \(d_{\beta,j}=H_j d_\gamma\) for all \(j\).* The quadratic term is
  **constant** along the ray — this is the flat direction created by the improper
  \(\gamma\) prior, and it is the only place propriety can fail. If \(d_\gamma=0\)
  then \(d_\beta=0\), excluded; so \(d_\gamma\ne0\), and by **(P3)** there is a group
  \(j_0\) with \(H_{j_0}d_\gamma\ne0\). By Step 2,
  \(\ell_{j_0}(\beta_{j_0}+tH_{j_0}d_\gamma)\to-\infty\) at least linearly with slope
  \(\le -c_{j_0}\|H_{j_0}d_\gamma\|\), while every other \(\ell_k\) is bounded above.
  So \(g\to-\infty\) at least **linearly**.

**4. Integrability.** \(g\) is concave, bounded above, and its recession function is
strictly negative in every direction; by compactness of the unit sphere there is
\(c>0\) and \(x_0=(\hat\gamma,\hat\beta)\) with
\(g(x)\le g(x_0)-c\,\|x-x_0\|\) for \(\|x-x_0\|\) large. Hence
\[
\int e^{g}
\;\le\;
\mathrm{const}\cdot e^{g(x_0)}\int e^{-c\|x-x_0\|}dx
\;<\;\infty . \qquad\square
\]

**Sharpness.** (P4) is **sufficient, not necessary**. The exact requirement in Case B
is that the **pooled profile** \(\gamma\mapsto\sum_j\ell_j(H_j\gamma)\) have compact
superlevel sets — i.e. a finite MLE for the GLM with pooled design rows
\(z_{j,i}^\top H_j\). Groups can borrow strength: a single separated group need not
destroy propriety if the others pin down \(\gamma\). The package's per-group
`groupef.estimable` gate is therefore **conservative** — it certifies propriety, but
failing it does not prove impropriety. *(Specializing to binomial with a flat prior
recovers the classical statement that posterior propriety under a flat prior
corresponds to existence of the MLE, i.e. absence of complete / quasi-complete
separation.)*

**Why this matters here.** It closes the loop between the two roles of estimability:

| Role of **(H3b)** | Where |
|---|---|
| Keeps \(1-\lambda^\star\) usably large (quantitative margin) | (iv), quantitative refinement |
| **Guarantees (H1) in the weak-prior limit**, hence a \(\lambda^\star\) uniform in \(\Lambda_\gamma\) | **Lemma R2a‴** + Cor. R2a″ (b3) |

For a **fixed proper** \(\Lambda_\gamma\succ0\) neither role is needed for mere
existence of \(\lambda^\star\) (step (iv)); both become essential as
\(\Lambda_\gamma\downarrow0\).

**Remark.** Steps (i)–(v), Corollary R2a″ and Lemma R2a‴ are complete for Lemma R2a′. The **existence** half of
claim (3) rests on **(H1)** plus rank/positivity; **(H3b)** and **(H2)** control the
**size** of \(1-\lambda^\star\). An earlier draft of this note attributed existence
itself to **(H3b)**; that was too strong — with a proper posterior, \(T<1\) pointwise
already forces \(F(t)\to0\).

### E.4 Degenerate example: separation degrades \(\lambda^\star\) (and kills it only if (H1) fails)

When **(H3b) fails** (no finite group-wise MLE), \(T(\beta)\le 1\) still holds at every
finite \(\beta\), and under a **proper** prior claim (3) still holds — but
\(1-\lambda^\star\) collapses, so the certified set is useless. If the population prior
is taken **flat** (\(\Lambda_\gamma=0\)) so that **(H1) fails**, claim (3) fails outright
— this is precisely the case excluded by hypothesis **(P4)** of Lemma R2a‴, and here
the pooled profile also has no finite maximizer, so \(\pi_0\) really is improper.
The example below is univariate only for readability; the mechanism
(\(\epsilon_j(\beta_j)\downarrow0\) along a recession direction of the group likelihood)
is the general one.

**Model (\(J=1\), one population mean, one group intercept).**

\[
\begin{aligned}
&y_i \mid \beta \;\stackrel{\text{ind}}{\sim}\; \mathrm{Binomial}(n_i,\, p_i),
\qquad \logit(p_i) = \beta,
\quad i=1,\ldots,n, \\
&\beta \mid \gamma \sim N(\gamma,\, \lambda_b^{-1}),
\qquad
\gamma \sim N(\mu_0,\, \lambda_\gamma^{-1})
\qquad\text{(or flat limit \(\lambda_\gamma\to 0\)).}
\end{aligned}
\]

Hyper-design: \(H_1=1\), \(q=1\), \(P_b=\lambda_b\), \(P_{11}^{\mathrm{RE}}=\lambda_b\)
(so **(i)** \(P_{11}^{\mathrm{RE}}\succ 0\) holds — rank of \(W\) is not the issue).

**Data: all successes.** \(y_i = n_i\) for every \(i\). Then \(\ell_j(\beta_j)\) has
**no finite maximum** (\(\beta_j\to\infty\)); **(H3b) fails**
(`groupef.estimable = FALSE`). The package blocks production sampling.

**Rate on the separation tail.** Fisher weights
\(w_i(\beta)=n_i p_i(\beta)\bigl(1-p_i(\beta)\bigr)\) with
\(p_i=\mathrm{logit}^{-1}(\beta)\). The group Gram matrix is
\(\mathcal{G}(\beta)=Z^\top W(\beta) Z = \bigl(\sum_i w_i(\beta)\bigr)\) (scalar
intercept). Hence \(W(\beta):=\sum_i w_i(\beta)\downarrow 0\) as \(\beta\to\infty\).
For \(J=1\), \(q=1\), Lemma R2a′ (iii) reduces to
\[
T(\beta)
=
\lambda_{\max}\bigl(A^{(0)}(\beta)\bigr)
=
\frac{\lambda_b}{\lambda_b + W(\beta)}.
\]
So \(T(\beta)\uparrow 1\) monotonically along the separation tail:
\[
\beta\to\infty
\;\Longrightarrow\;
W(\beta)\to 0
\;\Longrightarrow\;
T(\beta)\to 1^{-}.
\]
*(Same mechanism as Chapter-C03 logit §2.4: \(W_j\ll\lambda_b\) drives
\(\lambda^\star\uparrow 1\).)*

**Proper prior: claim (3) survives, but \(\lambda^\star\) is pushed to the ceiling.**
The binomial likelihood is bounded by \(1\), so with proper Gaussian priors on
\(\beta\mid\gamma\) and \(\gamma\) the posterior is proper and **(H1)** holds. Step 4
therefore still gives \(F(t)\to0\), and \(\lambda^\star<1\) exists. What breaks is the
**magnitude**: the posterior is pushed far up the separation tail (only the prior
resists), so for \(\delta_1=0.05\) the required \(\lambda^\star\) sits extremely close
to \(1\),
\[
1-\lambda^\star \;\approx\; \frac{W(\beta_{0.95})}{\lambda_b},
\qquad
W(\beta)=\textstyle\sum_i n_i p_i(\beta)(1-p_i(\beta)) \;\to\;0 ,
\]
with \(\beta_{0.95}\) the \(0.95\)-quantile of \(\pi(\beta\mid y)\), which grows with
the prior scale. The set \(E_{\lambda^\star}\) is then nearly all of the state space,
the \(\gamma\)-chain rate is \(\approx1\), and the Proposition R2 bound
\((1-\varepsilon_d)^n+\delta\) is valid but practically vacuous.

**Flat \(\gamma\) prior: (H1) fails and claim (3) fails.** In the limit
\(\lambda_\gamma\to0\) with all successes, there is no proper posterior on
\((\gamma,\beta)\): mass escapes along \(\beta\sim\gamma\to\infty\). Step 4 is exactly
what fails — \(\pi\) is not a probability measure, continuity from above does not
apply, and there is no \(\lambda^\star<1\) with a small spectral tail. This is the
sense in which the **flat-\(\gamma\) rate matrix \(A^{(0)}\) must be paired with a
proper posterior**: \(A^{(0)}\) is used as the **worst-case** rate (part (v)), not as
the actual prior.

**Weak \(\gamma\) prior does not help the rate.** At \(\Lambda_\gamma=0\), \(T\) is
already maximal; stronger \(\gamma\) priors only **decrease** \(\lambda_{\max}(A)\)
pointwise (part (v)), i.e. they help mixing but shrink \(E_{\lambda^\star}\) as a
\(\gamma\)-set — the tension is the same Pareto trade-off as §7.5.

**Takeaway.** Estimability (**H3b**) is not what makes \(\lambda^\star\) **exist**; it
is what keeps \(1-\lambda^\star\) **usably large** (and, with weak priors, what keeps
the posterior proper so **(H1)** holds at all). The package's
`groupef.estimable` gate is therefore a **quality / conditioning** gate for the
certificate, not a logical prerequisite for its existence.

---


## Appendix F — Proof of Proposition R2 (assembly)

This appendix combines Appendices A (refresh measure), B/C (minorization constant
and certified set), D (exhaustion, restricted-kernel stationarity), and E (spectral
cutoff) into the proof of Proposition R2 as stated in §8, then lifts the result from
the $\gamma$-chain to the joint $(\gamma,\beta)$ chain.

### F.1 Why $Q$ is fixed before $C$


Rosenthal's minorization uses a **single** refresh measure valid for every state in
the small set. We therefore specify \(Q\) **before** choosing any set. (The reverse
order — pick \(C\), then take \(Q\propto\inf_{x\in C}p(x,\cdot)\) — yields a measure
that depends on \(C\), is supported on \(C\), and is not Gaussian, so the family of
sets cannot be enlarged coherently and \(Q\) cannot be sampled.) With \(Q\) fixed, the
profile \(\varepsilon(\cdot)\) is a function of the model alone, and the **constant is
then chosen before the set**: pick \(\varepsilon_d\), and let \(C\) be the region where
it holds. Remark R2b-0 compares this with the alternative of fixing a set first — an
HPD region — and minimising \(\varepsilon\) over it.

**The \(\gamma\)-chain.** In two-block Gibbs (\(\beta\) first, then \(\gamma\)) the
\(\gamma\)-marginal is itself Markov, with kernel

\[
q(\gamma' \mid \gamma)
\;=\;
\int
\phi_q\bigl(\gamma';\, m(\beta),\, P_{11}^{-1}\bigr)\,
\pi(\beta \mid \gamma, y)\, d\beta,
\qquad
m(\beta)
=
P_{11}^{-1}\Bigl(\Lambda_\gamma\mu_0 + \sum_{j=1}^{J} H_j^\top P_b\,\beta_j\Bigr),
\]

with \(H_j\) the hyper design (\(\mathcal W_j\) of §7.1) and

\[
P_{11}
\;=\;
\underbrace{\Lambda_\gamma}_{\text{population prior}}
+
\underbrace{\sum_{j=1}^{J} H_j^\top P_b H_j}_{P_{11}^{\mathrm{RE}}},
\qquad P_b:=\Psi^{-1} .
\]

**(Q1) Precision.** The operative requirement is on the **precision**, not the
covariance:

\[
\boxed{\
\Lambda_Q \;:=\; \Sigma_Q^{-1}\;\succeq\;\Lambda_q(\gamma'\mid\gamma)
\;:=\;-\nabla^2_{\gamma'}\log q(\gamma'\mid\gamma)
\quad\text{for all }\gamma,\gamma',
\qquad\text{achieved by }\ \Sigma_Q:=P_{11}^{-1}=\Sigma^\star .
\ }
\]

**Why precision and not covariance.** The minorization needs
\(g(\gamma')=\log q(\gamma'\mid\gamma)-\log q_Q(\gamma')\) bounded below, and
\(\nabla^2 g=\Lambda_Q-\Lambda_q\). So what must dominate is the **curvature of the
log-density**: the refresh measure has to be at least as **concentrated** as the
transition kernel, or its tails die first and the ratio \(q/q_Q\) plunges to \(0\).
For Gaussians, \(\Sigma_Q\preceq\Sigma\) and \(\Lambda_Q\succeq\Lambda_q\) are the
same statement (inversion is antitone), which is why the covariance form is the one
usually quoted. **They part company for GLMMs**, where \(q(\cdot\mid\gamma)\) is a
mixture with no single \(\Sigma\) but a perfectly good pointwise \(\Lambda_q(\gamma')\).
Everything below is therefore written in precisions.

Three facts justify the choice \(\Lambda_Q=P_{11}\).

1. **It is exact, for every family.** The likelihood never involves \(\gamma\) —
   \(\gamma\) enters only through the Gaussian RE prior \(\beta_j\mid\gamma\sim
   N(H_j\gamma,\Psi)\) — so the Block-2 conditional is **exactly**
   \(\gamma\mid\beta \sim N(m(\beta), P_{11}^{-1})\), with no Laplace step. This is
   \(\Sigma^\star=P_\gamma^{-1}\) of `MINORIZATION_GAUSSIAN_REFRESH.md` §1.1 and the
   same \(P_{11}\) that normalizes the rate matrix \(A=P_{11}^{-1/2}SP_{11}^{-1/2}\).
2. **Precision dominance is an identity, not an assumption.** \(q(\cdot\mid\gamma)\)
   is a mixture of Gaussians with **common** covariance \(\Sigma^\star\) and random
   mean \(m(\beta)\). Differentiating twice through the mixture (Remark R2b-1) gives,
   for **every** GLM family and every \(\gamma,\gamma'\),
   \[
   \boxed{\ \Lambda_Q-\Lambda_q(\gamma'\mid\gamma)
   \;=\;\Lambda_Q\,\mathrm{Cov}\bigl(m(\beta)\,\big|\,\gamma,\gamma',y\bigr)\,\Lambda_Q
   \;\succeq\;0\ }
   \]
   where the covariance is over the **bridge** law \(\pi(\beta\mid\gamma,\gamma',y)\)
   of the RE block given both endpoints of the step. Dominance therefore holds
   automatically; it is **strict** exactly when that bridge law does not concentrate
   \(m(\beta)\) on a point, which is **(H3a)** plus nondegenerate \(V_j\).
   *(Under Gaussian closure this reduces to the familiar total-variance form
   \(\Sigma=\Sigma^\star+\mathrm{Cov}(m(\beta)\mid\gamma,y)\succeq\Sigma^\star\).)*
3. **It is a universal, and tight, upper bound.** The identity gives
   \(\Lambda_q\preceq P_{11}\) for every model, state and step, with the gap
   vanishing as the bridge covariance degenerates — so \(P_{11}\) cannot be lowered
   *in general*. It need not be optimal for a **given** model: since the leading
   constant is \(\sqrt{\det(\Lambda_q\Lambda_Q^{-1})}\), *decreasing* in
   \(\Lambda_Q\), any valid \(\Lambda_Q\succeq\sup_{\gamma,\gamma'}\Lambda_q\) smaller
   than \(P_{11}\) gives a **better** \(\varepsilon\). That is a sharpening, not a
   correction; \(P_{11}\) is the choice that needs no model-specific work.

**Why the \(\beta\)-conditional draw cannot be minorized directly.** It is tempting
to minorize the exact Gaussian step \(k(\gamma'\mid\beta)=\phi_q(\gamma';
m(\beta),\Sigma^\star)\) against \(q_Q=\phi_q(\cdot;\gamma^\star,\Sigma^\star)\). Here the
precision gap is **exactly zero** — \(\beta\) is held fixed, so the bridge covariance
in (2) vanishes — and with equal precisions the log-ratio is **affine** in \(\gamma'\),
\[
\log\frac{k(\gamma'\mid\beta)}{q_Q(\gamma')}
=
\bigl(m(\beta)-\gamma^\star\bigr)^\top P_{11}\gamma' + \text{const},
\]
so its infimum over \(\gamma'\in\mathbb R^q\) is \(-\infty\) and no \(\varepsilon>0\)
exists. **Integrating \(\beta\) out is essential**: it is the only source of the
strict precision gap, and hence of a log-ratio that is strictly convex with a finite
minimum.

**(Q2) Mean.**

\[
\boxed{\ \gamma^\star \;=\; M(\gamma^\star),
\qquad
M(\gamma) := E[\gamma'\mid\gamma]
=
P_{11}^{-1}\Bigl(\Lambda_\gamma\mu_0 + \sum_{j=1}^{J} H_j^\top P_b\, b_j(\gamma)\Bigr),
\quad
b_j(\gamma) := E[\beta_j\mid\gamma,y].\ }
\]

Existence and uniqueness of this fixed point is **Lemma 1A.2** for a proper prior and
**Lemma 1B.2** in the flat limit; the mechanism is **Lemma R2d** (contraction) and
**Cor. R2d-5** (variational form). Why the fixed point rather than any other centre is
the **Centering remark** in §7.3.3. By Cor. R2d-5 this \(\gamma^\star\) is precisely the **mode
of the \(\gamma\)-marginal posterior** \(\pi(\gamma\mid y)\), which is what makes it
computable in the flat-prior limit where the contraction argument degrades.


### F.2 Ingredients


The proof is seven short steps. \(Q\) is fixed first (**Section 1**), the profile
\(\varepsilon(\cdot)\) is then established (**Section 2**), and the certified set is
read off as one of its superlevel sets:

| # | Needed | Supplied by (§7.3) |
|---|---|---|
| 1 | \(\Sigma^\star=P_{11}^{-1}\succ0\) | **Lemma 1A.1** if \(\Lambda_\gamma\succ0\); **Lemma 1B.1** under **(H3a)** if \(\Lambda_\gamma=0\) |
| 2 | \(\gamma^\star\) with \(M(\gamma^\star)=\gamma^\star\), finite | **Lemma 1A.2** (strong convexity) / **Lemma 1B.2** (coercivity, needs **(H3b)**) |
| 3 | **precision** dominance \(\Lambda_q\preceq\Lambda_Q\) | **Remark R2b-1(1)** — an identity, free |
| 4 | \(\varepsilon(\cdot)>0\), unique minimiser, **continuous**, **log-concave** with peak at \(\gamma^\star\) | **Lemmas 2A.1–2A.3, 2A.3-F** / **2B.1–2B.2** — existence, no formula |
| 5 | \(\widetilde C_d=\{\varepsilon\ge\varepsilon_d\}\) convex **compact**, \(\min_{\widetilde C_d}\varepsilon=\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)>0\) attained; \(\pi_\gamma(\widetilde C_d^{\,c})<\delta_2\) | **Lemma 2A.3-G/2B.3** (compactness, needs only **(H3a)**) + **Lemma 2A.4/2B.4** (the constant, by continuity + compactness) + **Lemma R2a** (mass) |
| 6 | \(\lambda^\star\) with \(\pi_\gamma(E_{\lambda^\star}^c)<\delta_1\), \(E\) closed | **Cor. R2d-4**; \(\delta_1=0\) if \(\Lambda_\gamma\succ0\) |
| 7 | restricted kernel with \(\pi_\gamma(\cdot\mid C)\) stationary; TV bound against \(\pi_\gamma\) | **Lemma R2c** + **Remark R2c-3**; **Proposition R1** (§4) |


### F.3 Proof of Proposition R2


**Step 1 (fix \(Q\) — §7.3 Section 1).** Take \(\Sigma_Q=\Sigma^\star=P_{11}^{-1}\) by
**(Q1)**: this exists by **Lemma 1A.1** when \(\Lambda_\gamma\succ0\), and by
**Lemma 1B.1** under **(H3a)** when \(\Lambda_\gamma=0\). Take \(\gamma^\star\) by **(Q2)**:
this exists, is unique and is finite by **Lemma 1A.2** (\(\Phi\) is
\(\Lambda_\gamma\)-strongly convex) or **Lemma 1B.2** (\(\Phi_0\) is coercive, needing
**(H3b)**). Set \(Q=N(\gamma^\star,\Sigma^\star)\) with density \(q_Q\) — **Cor. 1A.3** /
**Cor. 1B.3**. Note \(Q\) depends on the model only: **not** on \(\delta\) and not on
any set.

**Step 2 (precision dominance).** By the identity of Remark R2b-1,
\(\Lambda_Q-\Lambda_q(\gamma'\mid\gamma)
=\Lambda_Q\mathrm{Cov}(m(\beta)\mid\gamma,\gamma',y)\Lambda_Q\succeq0\) for every
family: dominance holds **by identity**, with nothing to verify. So the log-ratio
\(g(\gamma'\mid\gamma)=\log q(\gamma'\mid\gamma)-\log q_Q(\gamma')\) is convex, and it
is strictly so under **(H3a)**. *(Covariance dominance \(\Sigma\succeq\Sigma^\star\) is
the same statement only when \(q\) is Gaussian; the precision form is what the argument
actually uses.)*

**Step 3 (the profile — §7.3 Section 2).** By **Lemmas 2A.1–2A.3** when
\(\Lambda_\gamma\succ0\), and **Lemmas 2B.1–2B.2** when \(\Lambda_\gamma=0\), the
minorization profile
\(\varepsilon(\gamma)=\inf_{\gamma'}q(\gamma'\mid\gamma)/q_Q(\gamma')\) — which always
lies in \([0,1]\) — is **strictly positive** for every \(\gamma\), attained at a
**unique** minimiser, and **continuous** in \(\gamma\). By **Lemma 2A.3-F** it is in
addition **log-concave**, with \(\nabla\log\varepsilon(\gamma)=P_{11}(\gamma^\star-M(\gamma))\),
so it attains its global maximum exactly at the refresh centre \(\gamma^\star\). The
two prior regimes differ only in how \(Q\) is produced; the \(\varepsilon\)-argument
itself is prior-free (Remark R2b-5). **No value of \(\varepsilon(\gamma)\) is
computed.**

**Step 4 (the level, the set, its minimum, and the budget split).** Fix
\(\delta\in(0,1)\) and \(\delta_1,\delta_2>0\) with \(\delta_1+\delta_2\le\delta\). For
\(d>0\) put

\[
\varepsilon_d\;:=\;e^{-d}\,\varepsilon(\gamma^\star)\;\in\;(0,1),
\qquad
\widetilde C_d\;:=\;\bigl\{\gamma:\varepsilon(\gamma)\ge\varepsilon_d\bigr\} .
\]

The set is **closed** by continuity of \(\varepsilon\) (Step 3), **convex** by
log-concavity (Step 3), and **compact** by **Lemma 2A.3-G** (resp. **2B.3**), which
sandwiches it between concentric ellipsoids
\(\{\|\gamma-\gamma^\star\|^2_{P_{11}^{\mathrm{RE}}}\le2d\}\subseteq\widetilde C_d
\subseteq\{\|\gamma-\gamma^\star\|^2_{S_\flat}\le2d\}\); it contains \(\gamma^\star\)
in its interior since \(\varepsilon(\gamma^\star)>\varepsilon_d\). Compactness here
consumes only **(H3a)** — it comes from nondegeneracy of the random-effect conditional
covariances, not from coercivity of the posterior potential — so unlike on the HPD
route it costs nothing extra in the flat limit.

With the set compact and \(\varepsilon\) continuous, **Weierstrass applies**, and by
**Lemma 2A.4** (resp. **2B.4**)

\[
\min_{\gamma\in\widetilde C_d}\varepsilon(\gamma)\;=\;\varepsilon_d
\quad\text{is attained on }\partial\widetilde C_d=\{\Psi=d\},
\]

so \(\varepsilon_d\) is not merely valid on \(\widetilde C_d\) but **sharp** there.
Strict positivity is \(\varepsilon(\gamma^\star)>0\) (**Lemma 2A.2**, resp. **2B.1**, at
the single point \(\gamma^\star\)). The gain over the HPD route is that this minimum has
a **closed form in \(d\)**, since the set was indexed by the function being minimised.
By **Lemma R2a** (exhaustion, which uses \(\varepsilon>0\) *pointwise*) the family is
nondecreasing with union \(\mathbb R^q\), so
\(\pi_\gamma(\widetilde C_d^{\,c})\downarrow0\) and one may choose
\(d=d(\delta_2)<\infty\) with \(\pi_\gamma(\widetilde C_d^{\,c})<\delta_2\);
quantitatively
\(\pi_\gamma(\widetilde C_d^{\,c})\le\pi_\gamma(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d)\).

**Step 5 (spectral truncation and the certified set).** By **Cor. R2d-4** pick
\(\lambda^\star\in(0,1)\) with \(\pi_\gamma(E_{\lambda^\star}^c)<\delta_1\), where
\(E_{\lambda^\star}=\{\gamma:\kappa(\gamma)\le\lambda^\star\}\) is closed. (If
\(\Lambda_\gamma\succ0\) this is vacuous: \(E_{\lambda^\star}=\mathbb R^q\) and
\(\delta_1=0\).) Put

\[
\boxed{\;C \;:=\; \widetilde C_d\cap E_{\lambda^\star}\;}
\]

— **compact**, being a closed subset of the compact \(\widetilde C_d\) (Lemma 2A.3-G),
and convex whenever the spectral truncation is vacuous — with, by the union bound,
\[
\pi_\gamma(C^c)\;\le\;\pi_\gamma(\widetilde C_d^{\,c})+\pi_\gamma(E_{\lambda^\star}^c)
\;<\;\delta_2+\delta_1\;\le\;\delta .
\]
So \(\pi_\gamma(C)\ge1-\delta>0\), conditioning on \(C\) is well defined, and the rate
on \(C\) is bounded by \(\lambda^\star<1\).

**Step 6 (minorization survives the intersection).** This is where Q-first pays. For
\(\gamma\in C\subseteq\widetilde C_d\), **Lemma 2A.5(a)** (resp. **2B.5**) gives
\(q(\gamma'\mid\gamma)\ge\varepsilon_d\,q_Q(\gamma')\) for **all**
\(\gamma'\in\mathbb R^q\), so integrating over an arbitrary Borel set,
\[
q(\gamma,A)\;\ge\;\varepsilon_d\,Q(A)
\qquad\forall\,\gamma\in C,\ \forall\,A\in\mathcal B(\mathbb R^q).
\]
Shrinking the set from \(\widetilde C_d\) to \(C\) costs **nothing**: the bound
\(\varepsilon(\gamma)\ge\varepsilon_d\) is a pointwise statement about \(\gamma\), so it
is inherited by every subset, and \(Q\) never depended on any set. Contrast the C-first
route, where intersecting two sets forces the refresh measure — and hence the whole
minorization — to be rebuilt from scratch. At this stage \(C\) is a small set of the
**unrestricted** \(\gamma\)-chain and no mass factor has been paid.

**Step 7 (restrict the sampler, pay once, and conclude).** By **Lemma R2c** applied to
\(\tilde C=C\times\mathbb R^{Jp_{\mathrm{re}}}\) (**Remark R2c-3**), the two-block
Gibbs kernel of the truncated target — exact \(\beta\)-draw, then
\(\gamma'\sim N(m(\beta),\Sigma^\star)\) truncated to \(C\) — stays in \(C\) and leaves
\(\pi_\gamma(\cdot\mid C)\) invariant. Its state space is now \(C\), so the refresh
measure must be a probability measure on \(C\) as well: take \(Q_C:=Q(\cdot\mid C)\).
Applying **Lemma 2A.5(b)** with \(C\) in place of \(\widetilde C_d\) — truncation only
*raises* the \(\gamma'\)-density on \(C\), so the Step 6 bound is inherited — yields
\[
q(\gamma,A\mid C)\;\ge\;\varepsilon_d\,Q(A)\;=\;\bigl[\varepsilon_d\,Q(C)\bigr]\,Q_C(A)
\;=\;\varepsilon\,Q_C(A),
\qquad
\varepsilon:=\varepsilon_d\,Q(C)\in(0,1],
\]
for all \(\gamma\in C\) and \(A\subseteq C\), the middle equality being
\(Q(A)=Q(C)Q_C(A)\) for \(A\subseteq C\), i.e. the definition \(Q_C(A)=Q(A)/Q(C)\).
The single factor \(Q(C)\) is the entire price of confining the sampler: it is the
refresh mass that Step 6 would have deposited outside \(C\), which a chain restricted
to \(C\) cannot use. And \(Q(C)>0\) for free — Step 5 gives \(\pi_\gamma(C)\ge1-\delta>0\)
and \(\pi_\gamma\) has a strictly positive Lebesgue density, so \(C\) has positive
Lebesgue measure, while \(Q\) is a nondegenerate Gaussian with strictly positive
density. When the spectral truncation is vacuous (\(\Lambda_\gamma\succ0\)) the mass is
bounded explicitly from below, \(Q(C)=Q(\widetilde C_d)\ge\Pr(\chi^2_q\le2d)\), by the
inscribed-ellipsoid argument of the §7.3.2 preamble — which works because \(Q\) and
\(\widetilde C_d\) are concentric at \(\gamma^\star\) (**Centering remark**, item 3).
Collecting Steps 4 and 7, the certified constant is explicit in \(d\) up to one scalar,
\(\varepsilon\ge\varepsilon(\gamma^\star)\,e^{-d}\Pr(\chi^2_q\le2d)\).

Rosenthal's Doeblin proposition now applies on \(C\), giving
\(\|q_n(\gamma,\cdot\mid C)-\pi_\gamma(\cdot\mid C)\|_{TV}\le(1-\varepsilon)^n\),
and **Proposition R1** converts to the full \(\gamma\)-target:
\[
\|q_n(\gamma,\cdot \mid C) - \pi_\gamma\|_{TV}
\;\le\;
(1-\varepsilon)^{n} + \pi_\gamma(C^c)
\;<\;
(1-\varepsilon)^{n} + \delta .
\]
\(\square\)

**Remark R2-J (lifting to the joint chain).** The joint two-block chain is a
**de-initialization** of the \(\gamma\)-chain: \((\gamma_n,\beta_n)\) is obtained from
\(\gamma_{n-1}\) by drawing \(\beta_n\sim\pi(\beta\mid\gamma_{n-1},y)\) and then
\(\gamma_n\sim\pi(\gamma\mid\beta_n)\). Since the extra \(\beta\)-draw is from the
**exact** conditional of the target, the data-processing inequality for total
variation gives, for the joint law,
\[
\|\mathcal L(\gamma_n,\beta_n)-\pi\|_{TV}\;\le\;\|\mathcal L(\gamma_{n-1})-\pi_\gamma\|_{TV},
\]
so every bound proved for the \(\gamma\)-chain transfers to the joint chain with the
same constants and a one-step offset. Scoping R2 to \(\gamma\) therefore costs
nothing, and buys a state space of dimension \(q\) instead of
\(q+Jp_{\mathrm{re}}\) — which is exactly why \(\varepsilon\) is not exponentially
small in \(J\).


## Appendix G — Scope and the symmetric-case checklist

### G.1 Symmetric case (sharper constants)


When the fixture satisfies **S1–S4** of `MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md`
§4 (balanced logit/probit, symmetric \(\gamma\)-marginal about \(\mu^\star\)):

| Step | Source |
|---|---|
| \(Q=N(\gamma^\star,\Sigma^\star)\) with \(\gamma^\star\) the fixed point of \(M\) | §7.3.1A (Lemmas 1A.1–1A.3); here \(\gamma^\star=\mu^\star\), the symmetry centre, and \(\gamma^\star=E[\gamma\mid y]\) by Cor. R2d-3 — consistently, symmetry makes mode \(=\) mean (Cor. R2d-5) |
| Flat population prior \(\Lambda_\gamma=0\) | §7.3.1B (Lemmas 1B.1–1B.3): \(\Sigma^{\star(0)}=(\sum_j\mathcal W_j^\top\Psi^{-1}\mathcal W_j)^{-1}\), \(\gamma^{\star(0)}\) the mode of \(\pi_0(\gamma\mid y)\) |
| \(\varepsilon_d\), \(Q\) with **certified constants** | Symmetric note §7.4, §8 (sharpens Remark R2b-1's \((1+\kappa)^{-q/2}e^{-d}\)) |
| \(\pi_\gamma(\widetilde C_d^{\,c}) < \delta_2\) with **explicit radii** | Raise the deficiency level \(d\) (Lemma R2a exhaustion + Remark R2a-2); symmetry makes \(\Psi\) even about \(\mu^\star\), so the \(\chi^2\) envelope of R2a-2 is exact to second order |
| \(\pi(T>\lambda^\star)<\delta_1\), \(\lambda^\star\in(0,1)\) | Appendix A.1 — **proved** from (H3a)+rank+**(H1)**; margin from (H2)/(H3b) |
| \(\pi_\gamma(\kappa>\lambda^\star)<\delta_1\) for the \(\gamma\)-chain | Cor. R2d-4 (§7.3); **vacuous** when \(\Lambda_\gamma\succ0\) |
| Attained rate \(\lambda_{C_d}=\max_{C_d}T\le\lambda^\star<1\) | Appendix A.2 (a) |
| Weak-prior limit: eigenvalue attained, ceiling \(1/(1+\epsilon)\) prior-free | Appendix A.2 (b1)–(b2) |
| One \(\lambda^\star\) uniform in \(\Lambda_\gamma\) (needs tightness / proper \(\pi_0\)) | Appendix A.2 (b3) |
| \(\pi_0\) proper at \(\Lambda_\gamma=0\) (**(H1)** in the flat limit) | Appendix A.3 — (H2)+rank+(H3b); equivalently \(\Phi_0=\sum_j\mathcal D_j\) coercive, Lemma 1B.2 |
| \(Q_0\), \(\gamma^{\star(0)}\), \(\varepsilon(\cdot\mid0)\) well defined at \(\Lambda_\gamma=0\) | §7.3.1B, §7.3.2B — needs (H3a) for \(\Sigma^{\star(0)}\), (H3b) for \(\gamma^{\star(0)}\) |
| \(\varepsilon_d(\Lambda_\gamma)\to\varepsilon_d(0)\): constants do not collapse | Lemma 2B.4 |
| \(\delta_1=0\) admissible (Gaussian, Gamma-log) | Appendix A.2 (c); also whenever \(\Lambda_\gamma\succ0\), Cor. R2d-4(a) |
| \(C = \widetilde C_d \cap E_{\lambda^\star}\), \(\pi_\gamma(C^c)<\delta\) | §7.2 Steps 4–5 |
| Restricted Gibbs \(\to \pi_\gamma(\cdot \mid C)\) | Lemma R2c + Remark R2c-3 (§7.3) |
| Full \(\pi_\gamma\) via R1; joint chain by de-initialization | §4; §7.2 Step 7 and Remark R2-J |

The main proof is §7.2; the symmetric fixture matters only because it supplies
**sharper constants** (\(\varepsilon_d\), radii) than Remark R2b-1 — and because it
supplies constants at all, which outside Gaussian closure §7.3.2 does not. Note that
S1–S4 also make the centre unambiguous: symmetry forces
\(M(\mu^\star)=\mu^\star\), so **(Q2)** is satisfied by inspection and the fixed-point
iteration of Lemma R2d is unnecessary.

Removing S1–S4 is **open item 1** of the symmetric note §10; Proposition R2 is the
**general** target.


### G.2 What Proposition R2 does not claim


- **The production sampler.** R2 certifies the Gibbs sampler **of the truncated
  target** \(\pi_\gamma(\cdot\mid C)\) (Lemma R2c, Remark R2c-3), not the
  unrestricted kernel the package runs. Approach B handles escape for the full
  kernel.
- **A formula for the minorization constant.** This is the sharpest limitation, and
  the one most easily mistaken. The \(d\)-dependence is exact and the constant is
  **sharp**: \(\min_{\widetilde C_d}\varepsilon=e^{-d}\varepsilon(\gamma^\star)\),
  attained (Lemma 2A.4). But the **scalar \(\varepsilon(\gamma^\star)\) is not
  evaluated**. What is proved in general is existence and regularity: \(\varepsilon\)
  is continuous and log-concave with \(\varepsilon(\gamma)>0\) everywhere and a unique
  inner minimiser (§7.3.2A for \(\Lambda_\gamma\succ0\), §7.3.2B for
  \(\Lambda_\gamma=0\)). See Remark R2b-0. The associated tail budget is bounded —
  \(\pi_\gamma(\widetilde C_d^{\,c})\le\pi_\gamma(\|\gamma-\gamma^\star\|^2_{P_{11}}>2d)\)
  — but turning that into a value for \(d\) needs a posterior spread constant \(\rho\),
  which is supplied only under Gaussian closure (Remark R2a-2).
- **An exact shape for the certified set.** \(\widetilde C_d\) is proved **compact and
  convex** and is trapped between two concentric ellipsoids (Lemma 2A.3-G), but outside
  Gaussian closure it is not itself an ellipsoid and no exact inequality description of
  its boundary is available, nor an exact value of \(Q(\widetilde C_d)\) — only the
  inscribed bound \(\ge\Pr(\chi^2_q\le2d)\). Membership testing is cheap (evaluate the
  profile, or equivalently \(\pi(\gamma\mid y)/q_Q(\gamma)\) via the complementarity
  identity); describing the boundary is not. The outer ellipsoid is also not sharp: it
  uses a global curvature ceiling \(G_j\), so for families with unbounded curvature
  (Poisson, Gamma with log link) one must fall back on the coarser radius \(1+d/c\).
- **The closed forms are Gaussian-closure only.** The boxed \(\mathcal D\), the
  ellipsoidal certified sets, and the leading constant \((1+\kappa)^{-q/2}e^{-d}\)
  **all** require \(\Lambda_q\) constant in \(\gamma'\), i.e. LMMs. For GLMMs the
  curvature sandwich of Remark R2b-1(3) still holds but bounds \(\nabla^2g\) from the
  wrong side; a certified number would need a Cramér–Rao lower bound on the bridge
  covariance plus the strong-convexity estimate of Remark R2b-0, and how conservative
  that is has not been quantified (Remark R2b-3). Note what has **stopped** being
  closure-dependent: convexity of the certified set and concavity of \(\log\varepsilon\)
  are now theorems for every family (Lemma 2A.3-F).
- **Numerically certified constants.** Even under closure,
  \(\varepsilon\ge(1+\kappa)^{-q/2}e^{-d}Q(C)\) is explicit but not yet *evaluated*:
  it needs a numerical \(\kappa\) (Lemma R2d), a \(d\) meeting the tail budget, and
  \(Q(C)\). See §7.2 Status items 1–3.
- **Not a caveat: precision dominance.** Worth recording, since an earlier draft of
  this note listed it as one. \(\Lambda_q\preceq\Lambda_Q\) is an identity
  (Remark R2b-1), and \(\varepsilon(\gamma)>0\) needs only that \(m(\beta)\) has full
  support, which is **(H3a)**. No uniform curvature gap is assumed, and no family is
  excluded on that ground.
- **Uniform** \(\varepsilon\) independent of \(\delta\): shrinking \(\delta\) forces a
  larger radius \(d\) and hence a smaller \(\varepsilon_d\). R2 does not remove that
  tension; it makes the exchange rate explicit and polynomial,
  \(\varepsilon_d\gtrsim\varepsilon(\gamma^\star)\delta_2^{\rho}\) (Remark R2a-2).
- **Not a caveat: compactness of \(C\).** Earlier drafts listed this first as a gap,
  then as an expensive property, then briefly as unnecessary. It is a **theorem, and it
  is load-bearing**: Lemma 2A.3-G, from
  \(\nabla^2\Psi=\sum_jH_j^\top P_bV_j(\gamma)P_bH_j\succ0\), which holds because the
  random-effect conditional covariances are nondegenerate. Together with continuity of
  \(\varepsilon\) it is precisely what produces the minorization constant, by
  Weierstrass: \(\varepsilon_d=\min_{\widetilde C_d}\varepsilon\), attained
  (Lemma 2A.4). It also gives the attained rate \(\lambda_{C_d}\) (Appendix A.2), a
  bounded numerical search, and the safe-region language of §6. What is cheap about it
  is only its **price**: **(H3a)** and nothing else, and **prior-free**, so it transfers
  to \(\Lambda_\gamma=0\) verbatim — unlike coercivity of \(\bar\Phi\), which there
  needs **(H3b)** — and it comes with an explicit containment ellipsoid.
- **Closed-form** \(R_\gamma(\delta)\) for all links: Poisson/cloglog need
  link-specific tail bounds (`CONDITIONAL_TAIL_BOUNDS_ONE_DIM_GLMM.md`), which is
  also what \(\rho\) in Remark R2a-2 requires.

---


---

## References (in repo)

- `inst/MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` — symmetric GLM minorization (§4–§8); R2 proof program
- `inst/CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` §I.3–§I.4 — \(P\)-blocks, rate matrix \(A\); Appendix A.1 proof
- `inst/LOGIT_STATIC_TAIL_CERTIFICATION.md` §4.5–§4.6, §7.1 — tail mass / Approach A
- `inst/MINORIZATION_GAUSSIAN_REFRESH.md` — \(Q\), \(\varepsilon_d\)
- `inst/PREFLIGHT_model_setup.md` — full rank, estimability (H3); Appendix A.1, A.3
- `lmebayes/references/minor.pdf` — Rosenthal (1995), Proposition 2
