# MVN-calibrated certified set for population minorization

Design note for **Chapter C05** / restricted two-block Gibbs certification.
Calibrate the minorization set \(\widetilde C_d\) from a **reference**
\(\mathcal N(0,I_q)\) tail probability, instead of inverting
\(\pi_\gamma(\widetilde C_d^{\,c})\) numerically.

**Companions:** `CHAPTER_C05_IMPLEMENTATION.md`, `restricted_gibbs_minorization_v2.md`,
`GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md`.

**Status:** draft design; not shipped.

---

## 0. Notation (fixed — do not overload)

Use **Chapter C05 symbols only**. Population parameters are **\(\gamma\)** (code:
`gamma`), never \(\delta\).

| Symbol | Meaning | Role |
|---|---|---|
| **\(\delta\)** | **Tail / escape probability budget** | User input; calibrates \(r(\delta)\). Target: \(\pi_\gamma(\widetilde C_d^{\,c}\mid y)\lesssim\delta\) (exact in MVN limit). |
| **\(r(\delta)\)** | Reference radius | \(r(\delta)^2=q_{1-\delta}(\chi^2_q)\) from \(\mathcal N(0,I_q)\) (§2). |
| **\(d\)** | Deficiency level | \(d=r(\delta)^2/2\); certified set \(\widetilde C_d=\{\Psi\le d\}\). |
| **\(\widetilde C_d\)** | **Central minorization set** | Where refresh minorization is proved; \(\pi_\gamma(\widetilde C_d)\approx 1-\delta\). |
| **\(Q\)** | Refresh law | \(Q=\mathcal N(\gamma^\star,P_{11}^{-1})\) (fixed). |
| **\(\varepsilon(\gamma)\)** | Kernel–to–\(Q\) profile | \(\varepsilon(\gamma)=\inf_{\gamma'} q(\gamma'\mid\gamma)/q_Q(\gamma')\). |
| **\(\varepsilon(\gamma^\star)\)** | **Mode profile** (Stage 2) | One scalar from the model; **not** fixed by \(\delta\). |
| **\(\varepsilon_d\)** | **Kernel floor on \(\widetilde C_d\)** | \(\varepsilon_d=e^{-d}\varepsilon(\gamma^\star)\); unrestricted minorization \(q\ge\varepsilon_d Q\). |
| **\(\varepsilon\)** | **Theorem 2 multiplier** | \(\varepsilon=\varepsilon_d\,Q(\widetilde C_d)\); restricted chain only. |

**Do not use \(\varepsilon\) for the tail budget** (that caused confusion with
\(\varepsilon\), \(\varepsilon_d\), and \(\varepsilon(\gamma^\star)\)).

**Certificate triple** (Theorem 2): \((\widetilde C_d,\,\varepsilon,\,Q)\) with escape
mass budget \(\delta\).

---

## 1. Goal

Given a tail budget **\(\delta\)**, construct:

1. **Set** \(\widetilde C_d\) on which Block‑2 minorization holds.
2. **Constants** \(\varepsilon_d\) and \(\varepsilon=\varepsilon_d Q(\widetilde C_d)\)
   for the restricted certificate.

**Idea.** Fix **\(r(\delta)\)** from a reference Gaussian, define \(\widetilde C_d\)
by a posterior **ratio cut** to \(\gamma^\star\), set **\(d=r(\delta)^2/2\)**. As
\(\pi_\gamma\to\mathcal N\) (many groups), \(\pi_\gamma(\widetilde C_d^{\,c})\to\delta\).

---

## 2. Reference radius \(r(\delta)\) from \(\mathcal N(0,I_q)\)

In **Q-whitened** coordinates \(u=P_{11}^{1/2}(\gamma-\gamma^\star)\), the refresh
density is \(q_Q(u)\propto\exp(-\|u\|^2/2)\). Under \(u\sim\mathcal N(0,I_q)\),

\[
\exp\!\Bigl(-\tfrac12\|u\|^2\Bigr)
\;\le\;
\exp\!\Bigl(-\tfrac{r^2}{2}\Bigr)
\quad\Longleftrightarrow\quad
\|u\|^2 \;\ge\; r^2.
\]

Define **\(r(\delta)\)** by

\[
\boxed{
\Pr_{u\sim\mathcal N(0,I_q)}
\bigl(
\|u\|^2 \ge r(\delta)^2
\bigr)
\;=\;
\delta.
}
\]

Equivalently,

\[
\boxed{
r(\delta)^2 \;=\; q_{1-\delta}(\chi^2_q),
\qquad
d \;=\; \tfrac{r(\delta)^2}{2} \;=\; \tfrac12\,q_{1-\delta}(\chi^2_q).
}
\]

R: `r2 <- qchisq(1-delta, df=q)`; `d <- 0.5 * r2`.

---

## 3. Certified set \(\widetilde C_d\) via posterior ratio cut

Let

\[
R(\gamma)
\;:=\;
\frac{\pi_\gamma(\gamma\mid y)}{\pi_\gamma(\gamma^\star\mid y)}.
\]

**Central set** (minorization region), with **\(r=r(\delta)\)** fixed in §2:

\[
\boxed{
\widetilde C_d
\;:=\;
\Bigl\{
\gamma :
R(\gamma)
\;\ge\;
\exp\!\bigl(-\tfrac{r(\delta)^2}{2}\bigr)
\;=\;
e^{-d}
\Bigr\}.
}
\]

**Escape set:** \(\widetilde C_d^{\,c}=\{R(\gamma)<e^{-d}\}\).

Under MVN closure with \(J\to I\), \(R(\gamma)=\exp(-\|u\|^2/2)\) and
\(\widetilde C_d=\{\|u\|^2\le r(\delta)^2\}\), so
\(\pi_\gamma(\widetilde C_d^{\,c})=\delta\). More generally
\(\pi_\gamma(\widetilde C_d^{\,c}\mid y)\to\delta\) as \(\pi_\gamma\) concentrates.

This matches the standard C05 set \(\widetilde C_d=\{\Psi\le d\}\) in the Gaussian
limit (§4).

---

## 4. Connection to Chapter C05

| Object | Formula / role |
|---|---|
| \(\Psi(\gamma)\) | \(\log[\pi_\gamma(\gamma)/q_Q(\gamma)]-\log[\pi_\gamma(\gamma^\star)/q_Q(\gamma^\star)]\) |
| Standard C05 set | \(\widetilde C_d=\{\gamma:\varepsilon(\gamma)\ge\varepsilon_d\}=\{\Psi\le d\}\) |
| This note | Same \(d=r(\delta)^2/2\); ratio cut implements \(\widetilde C_d\) without inverting \(\pi(\widetilde C_d^{\,c})\) for \(d\) |
| Escape budget | \(\pi_\gamma(\widetilde C_d^{\,c})\le\delta\) (Theorem 2 error term) |

---

## 5. Minorization constants once \(\delta\) and \(d\) are fixed

From `restricted_gibbs_minorization_v2.md` §5 (Lemma 16–17).

### 5.1 Kernel floor \(\varepsilon_d\) (on the central set)

For every \(\gamma\in\widetilde C_d\), **unrestricted** one-step kernel:

\[
q(\gamma,A)\;\ge\;\varepsilon_d\,Q(A),
\qquad
\boxed{
\varepsilon_d \;=\; e^{-d}\,\varepsilon(\gamma^\star)
\;=\;
\exp\!\bigl(-\tfrac{r(\delta)^2}{2}\bigr)\,\varepsilon(\gamma^\star).
}
\]

This is the **largest** constant valid on all of \(\widetilde C_d\) (sharp on
\(\partial\widetilde C_d=\{\Psi=d\}\)).

### 5.2 Theorem 2 multiplier \(\varepsilon\) (restricted chain)

For the **restricted** sweep on \(\widetilde C_d\) with truncated refresh
\(Q_{\widetilde C_d}=Q(\cdot\mid\widetilde C_d)\):

\[
q(\gamma,A\mid\widetilde C_d)\;\ge\;\varepsilon\,Q_{\widetilde C_d}(A),
\qquad
\boxed{
\varepsilon \;=\; \varepsilon_d\,Q(\widetilde C_d).
}
\]

**Mass discount.** Corollary 15: \(Q(\widetilde C_d)\ge\Pr(\chi^2_q\le 2d)=\Pr(\chi^2_q\le r(\delta)^2)=1-\delta\). Hence

\[
\varepsilon
\;\ge\;
\varepsilon(\gamma^\star)\,
e^{-d}\,(1-\delta)
\;=\;
\varepsilon(\gamma^\star)\,
\exp\!\bigl(-\tfrac{r(\delta)^2}{2}\bigr)\,(1-\delta).
\]

### 5.3 What \(\delta\) does and does not fix

| Fixed by \(\delta\) (via \(r,d\)) | Not fixed by \(\delta\) |
|---|---|
| \(r(\delta)\), \(d=r^2/2\) | \(\varepsilon(\gamma^\star)\) (Stage 2) |
| Ratio cut \(e^{-d}\) defining \(\widetilde C_d\) | Full kernel \(q(\gamma'\mid\gamma)\) |
| Factor \(e^{-d}\) in \(\varepsilon_d\) | |
| Lower bound \(Q(\widetilde C_d)\ge 1-\delta\) | |

**Stage 2** (`CHAPTER_C05_IMPLEMENTATION.md` §3): compute \(\varepsilon(\gamma^\star)\)
once (closure: \(\det(I+\tilde J)^{-1/2}\); else numerical).

---

## 6. End-to-end recipe

1. Input: dimension \(q\), tail budget **\(\delta\)**, mode \(\gamma^\star\), \(P_{11}\).
2. `r2 <- qchisq(1-delta, df=q)`; **`d <- r2/2`**.
3. **Set:** \(\widetilde C_d=\{R(\gamma)\ge e^{-d}\}\).
4. **Mode profile:** \(\varepsilon(\gamma^\star)\) (Stage 2).
5. **Floor:** \(\varepsilon_d = exp(-d) * eps_star`.
6. **Multiplier:** \(\varepsilon = \varepsilon_d * Q(C_d)\) (lower bound `eps_star * exp(-d) * (1-delta)`).
7. **Theorem 2:** \(\|q_n(\cdot\mid C_d)-\pi_\gamma\|_{TV}\le (1-\varepsilon)^n+\delta\).

---

## 7. What this route does not solve

- **Finite-\(J\) gap:** \(\pi_\gamma(\widetilde C_d^{\,c})\) may differ from \(\delta\)
  before the MVN limit — bound with `GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md` or
  sandwich (Lemma 14).
- **\(R(\gamma)\)** off closure needs \(\bar\Phi\) / integrated likelihood.
- **LSD population envelope** (`BUILDING_LSD_MARGINAL_POPULATION_DENSITY.md`) is a
  separate **majorization** path, not this **minorization** set construction.

---

## 8. Open items

1. Finite-\(J\) error: \(\pi(\widetilde C_d^{\,c})\) vs \(\delta\).
2. Ratio cut vs \(\{\Psi\le d\}\) when \(\bar\Phi>0\).
3. MVN scaling heuristics for log-concave non-Gaussian likelihoods (`Chapter-A08` §2.13).

---

## See also

- `inst/CHAPTER_C05_IMPLEMENTATION.md`
- `inst/BUILDING_LSD_MARGINAL_POPULATION_DENSITY.md`
- `inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md`
