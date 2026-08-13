# Rate-\(A^c\) bounds by likelihood (draft)

Likelihood-specific **§3-style** bounds on the target escape probability
\(\pi(A^c \mid \gamma, y)\), parallel to
`LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §3, but using the **rate safe
region** \(A / A^c\) from `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md`
Part IV (average shrinkage / \(\lambda^\star\)).

**Status.** Draft. §1–§6 give likelihood-specific bounds on
\(\pi(A^c \mid \gamma, y)\) (tilted-Gaussian tails via certified inner boxes).
§8–§11 complete the Markov-chain argument: Foster–Lyapunov drift on
\(\gamma\), bounds on \(P^n(x, A^c)\) and \(\pi(A^c)\), coupling gap, and
assembled TV. Normal and Gamma(log) are trivial (\(A^c = \varnothing\) when
certified). Open items in §13.

**Companion notes.**

- `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` — \(P\)-blocks, \(\omega_j\),
  \(\bar\omega\), choosing \(\tau\) (§1.4, §IV.5).
- `LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` — full logit single-group
  argument (§3–§7).
- `ELLIPSOID_TV_BOUND.md` — generic TV split (Lemma 1); uses a different
  safe-region definition, not developed here.
- `MINORIZATION_GAUSSIAN_REFRESH.md` — Gaussian \(Q\) with
  \(\Sigma^\star\), \(\mu^\star\); inner \(\gamma'_\star(\gamma)\) and
  \(\delta(\gamma)\) (standalone minorization note).
- `MINORIZATION_SYMMETRIC_NON_GAUSSIAN.md` — symmetric logit/probit (and
  related) population posteriors: mean = mode = center; non-Gaussian \(q\);
  numeric two-stage \(\varepsilon_d\) and \(\delta\)-level set \(C_d\).

---

## 0. Common setup and target quantity

### 0.1 Model and update order

**Population parameter** \(\gamma\) (block updated **second** in the package;
scalar specialization \(\gamma \sim N(\mu_0, \sigma_0^2)\) with precision
\(\lambda_\gamma = \sigma_0^{-2}\)). **Random effects** \(\beta_j\) (updated
**first**), one scalar RE per group (\(p_{\mathrm{re}}=1\), \(z_{j,i}=1\),
\(H_j = x_j \equiv 1\)):

\[
\beta_j \mid \gamma \;\sim\; N(\gamma,\, \lambda_b^{-1}),
\qquad
j = 1,\ldots,J.
\]

Group \(j\) contributes likelihood terms \(\ell_j(\beta_j) = \sum_{i\in j}
\ell_{j,i}(\eta_{j,i})\) with \(\eta_{j,i} = \mathrm{offset}_{j,i} + \beta_j\)
(scalar intercept RE per group; offsets absorb fixed effects).

Write \(\beta = (\beta_1,\ldots,\beta_J)^\top\) and the **target**
\(\pi(\beta, \gamma \mid y)\). Given \(\gamma\), groups factorize:

\[
\pi(\beta \mid \gamma, y) \;=\; \prod_{j=1}^{J} \pi(\beta_j \mid \gamma, y_j),
\qquad
\pi(\beta_j \mid \gamma, y_j) \;\propto\;
\exp\bigl(\ell_j(\beta_j)\bigr)\,\phi(\beta_j; \gamma, \lambda_b^{-1}),
\]

where \(\phi(\cdot; \gamma, \lambda_b^{-1})\) is the Gaussian prior density.

### 0.2 Rate safe region and escape functional

Fix certification threshold \(\tau \in (0,1)\). Define shrinkage weights
(§IV.2 of `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md`)

\[
\omega_j(\beta) := \frac{\lambda_b}{\lambda_b + W_j(\beta)},
\qquad
W_j(\beta) := \sum_{i\in j} w_{j,i}(\eta_{j,i}(\beta)),
\qquad
\bar\omega(\beta) := \frac{1}{J}\sum_{j=1}^{J} \omega_j(\beta).
\]

With flat \(\gamma\) prior (\(\lambda_\gamma = 0\)) the Gibbs rate is
\(\lambda^\star(\beta) = \lambda_b\,\bar\omega(\beta)\). The **rate safe
set** and **unsafe complement** are

\[
A := \bigl\{\beta : \lambda^\star(\beta) < \tau\bigr\},
\qquad
A^c := \bigl\{\beta : \lambda^\star(\beta) \ge \tau\bigr\}.
\]

(With \(\lambda_\gamma > 0\), replace \(\tau\) by the equivalent threshold on
\(\bar\alpha = \bar\omega/\lambda_b\); see Part IV of the \(P\)-matrix note.)

The **one-step escape functional** for the \(\gamma\)-chain (block updated
second, after a \(\beta\)-draw) is

\[
g(\gamma) \;:=\; \pi(A^c \mid \gamma, y)
\;=\;
\int \mathbf 1\{\beta \in A^c\}\,\pi(\beta \mid \gamma, y)\,d\beta
\;=\;
E_{\pi(\cdot\mid\gamma,y)}\bigl[\mathbf 1\{\beta \in A^c\}\bigr].
\]

This is the object bounded in §1–§6 for each likelihood. §8–§11 bound
\(P^n(x, A^c)\), \(\pi(A^c)\), and the escape gap
\(|\pi(A^c) - P^n(x, A^c)|\) for the \(\gamma\)-marginal chain.

### 0.3 Shared lemmas (all GLM cases)

These are copied from `LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §3; the
only likelihood-specific input is the **score bound** \(L_j\).

**Lemma T (tangent bound).** If \(\ell_j''(\beta_j) \le 0\) everywhere and
\(|\ell_j'(\beta_j)| \le L_j\) for all \(\beta_j\), then for any \(b_0\),
\(\ell_j(\beta_j) \le \ell_j(b_0) + L_j|\beta_j - b_0|\).

**Lemma G (tilted Gaussian tail).** For \(Z \sim N(\gamma, \lambda_b^{-1})\)
with density \(\phi\), and \(c \ge 0\),

\[
\int_t^\infty e^{c z}\,\phi(z;\gamma,\lambda_b^{-1})\,dz
\;=\;
\exp\!\Big(c\gamma + \tfrac{c^2}{2\lambda_b}\Big)\,
\Phi\!\Big(-\tfrac{t - \gamma - c/\lambda_b}{\lambda_b^{-1/2}}\Big).
\]

**Lemma Z (normalizer lower bound).** With
\(Z_j(\gamma) := \int e^{\ell_j(\beta_j)}\phi(\beta_j;\gamma,\lambda_b^{-1})\,d\beta_j\),

\[
Z_j(\gamma, y_j) \;\ge\;
\exp\!\Big(\ell_j(\gamma) - L_j\,\lambda_b^{-1/2}\sqrt{2/\pi}\Big).
\]

*Proof.* Jensen plus \(E|\beta_j - \gamma| = \lambda_b^{-1/2}\sqrt{2/\pi}\)
under the prior. \(\blacksquare\)

**Template (single-group escape).** For a **one-sided or two-sided tail event**
\(E_j \subset \mathbb R\) (e.g.\ \(|\beta_j| > b_j^\star\)), if
\(A^c \subseteq \bigcup_j E_j\), then

\[
\pi(A^c \mid \gamma, y)
\;\le\;
\sum_{j=1}^{J} \pi(E_j \mid \gamma, y_j)
\;\le\;
\sum_{j=1}^{J} G_j(\gamma),
\]

where each \(G_j\) is the **tilted-Gaussian bound** of §2.2 applied to group
\(j\) (one- or two-sided tail as in §2–§5).

For a **sufficient inner box** \(B \subseteq A\),

\[
\pi(A^c \mid \gamma, y) \;\le\; \pi(B^c \mid \gamma, y)
\;\le\;
\sum_{j=1}^{J} \pi(\beta_j \notin B_j \mid \gamma, y_j),
\]

which is the route used for multi-group GLMs in §2–§5 when an explicit
\(B\) is certified inside \(A\).

### 0.4 Choosing the inner box \(B\)

**Yes — \(B\) is a modeller's choice.** The rate-unsafe set \(A^c\) is
fixed by the likelihood, priors, and threshold \(\tau\) via
\(\bar\omega(\beta) \ge \tau/\lambda_b\) (flat \(\gamma\)). The box \(B\) is
an **auxiliary** set you introduce for bounding: any measurable
\(B \subseteq A\) gives the valid inequality
\(\pi(A^c \mid \gamma, y) \le \pi(B^c \mid \gamma, y)\). You are not
required to use a box; it is the standard closed-form route when \(A^c\) has
no simple description.

**Typical construction (logit).** Pick margins \(b_j^\star > 0\) (often a
common \(b^\star\) for all \(j\)), set
\(\kappa_j := N_j\, p(b^\star)(1-p(b^\star))\), and define
\[
B \;:=\; \prod_{j=1}^{J} [-b^\star, b^\star]
\quad\text{or}\quad
B \;:=\; \prod_{j=1}^{J} [-b_j^\star, b_j^\star].
\]
**Verify** \(B \subseteq A\) by showing \(\bar\omega(\beta) < \tau/\lambda_b\)
for every \(\beta \in B\) — in practice via a **floor** on each
\(W_j(\beta) \ge \kappa_j\) inside \(B\) and the average inequality
\((1/J)\sum_j 1/(\lambda_b + \kappa_j) < c_\tau\) (§2.3). Offsets and
design enter only through which \(\beta\) values keep all \(|\eta_{j,i}|\le
b^\star\).

**Tradeoff.** Larger \(B\) (bigger \(b^\star\)) keeps more of \(A\) inside
\(B\) and is desirable if certification succeeds, but \(\pi(B^c \mid \gamma,y)\)
is smaller (tighter escape bound). Smaller \(B\) is easier to certify
\(B \subseteq A\) (stronger curvature floors) but \(\pi(B^c \mid \gamma,y)\)
is **larger** — a looser upper bound on \(\pi(A^c \mid \gamma,y)\). Different
valid choices of \(B\) yield different upper bounds; optimize \(b^\star\) (or
asymmetric per-group margins) subject to \(B \subseteq A\).

**Poisson / cloglog (lower floors).** Pick margins \(s_j > 0\) (rate floors
\(e^{\eta_{j,i}} \ge s_j\), equivalently \(\eta_{j,i} \ge \log s_j\)) and
certify
\[
\frac{1}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + n_j s_j} \;<\; c_\tau
\quad\Longrightarrow\quad
B \;:=\; \prod_{j=1}^{J} [\log s_j, \infty) \;\subseteq\; A.
\]
Tune each \(s_j\) (subject to the **average** inequality) to keep \(B\) as
large as possible for a tight \(\pi(B^c \mid \gamma,y)\) bound. **Special
case:** when groups are exchangeable (e.g.\ common \(n_j\), \(z_{j,i}=1\),
equal offsets within group), one may take \(s_j \equiv s\); then
\((1/J)\sum_j 1/(\lambda_b + n s) < c_\tau\) and
\(B = \{\beta : \beta_j \ge \log s\ \forall j\}\).

**Other shapes.** \(B\) need not be a product of intervals: any set with
\(B \subseteq A\) and tractable \(\pi(B^c \mid \gamma,y)\) works (e.g.\
ellipsoids in multivariate RE). The draft uses intervals because they match
§2.2–§3 tilted-Gaussian tails.

---

## 1. Normal (Gaussian regression)

### 1.1 Weights and rate region

\(w_{j,i} = 1/\phi\) constant, \(W_j = n_j/\phi\), \(\omega_j = \lambda_b /
(\lambda_b + n_j/\phi)\) independent of \(\beta\). Hence \(\lambda^\star\)
is **constant** (§1.4 of the \(P\)-matrix note).

### 1.2 Closed-form bound on \(\pi(A^c \mid \gamma, y)\)

**Proposition N1.** At fixed \((\phi, \lambda_b, \lambda_\gamma, \{n_j\},
\tau)\),

\[
\pi(A^c \mid \gamma, y) \;\in\; \{0, 1\},
\qquad
\pi(A^c \mid \gamma, y) = \mathbf 1\{\lambda^\star \ge \tau\},
\]

with \(\lambda^\star\) computed once from §1.4. In particular, if priors and
design are calibrated so \(\lambda^\star < \tau\) (equivalently \(A^c =
\varnothing\)), then \(g(\gamma) \equiv 0\) for every \(\gamma\).

*Proof.* Immediate from \(\beta\)-independence of \(\lambda^\star\). Normal
is the only likelihood where \(A^c\) can be eliminated by threshold/design
choice alone (§IV.6). \(\blacksquare\)

**Remark.** No tangent/tail machinery is needed; §0.3 is vacuous here.

---

## 2. Logit (binomial, canonical link)

### 2.1 Weights and rate region

\(w_{j,i} = n_{j,i}\, p_{j,i}(1-p_{j,i})\), \(W_j(\beta) = \sum_{i\in j}
w_{j,i}\). Shrinkage \(\omega_j(\beta) = \lambda_b / (\lambda_b + W_j(\beta))\).
**Two-sided** geometry: \(W_j \to 0\) as \(|\eta_{j,i}| \to \infty\) (both
tails), so \(A^c\) is typically a **thick** region outside a central box in
\(\beta\)-space (§IV.3).

### 2.2 Single group (\(J = 1\)) — complete bound

This is `LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §2–§3, restated with
rate notation. Fix \(b^\star > 0\), \(\kappa := N\, p(b^\star)(1-p(b^\star))\)
with \(N = \sum_i n_i\). The **sufficient inner safe interval**

\[
B = [-b^\star, b^\star]
\;\subseteq\;
A
\quad\text{whenever}\quad
\frac{\lambda_b}{\lambda_b + \kappa} \;<\; \tau/\lambda_b
\]

(certify via §IV.3: pick \(\kappa\) from the average inequality). On the
**tail event** \(E = \{|b| > b^\star\} \supseteq A^c \cap B^c\),

\[
\pi(E \mid \gamma, y) \;\le\; G(\gamma)
\]

with

\[
\boxed{
G(\gamma) :=
\exp\!\Big(L\,\lambda_b^{-1/2}\sqrt{2/\pi} + \tfrac{L^2}{2\lambda_b}\Big)
\left[
\Phi\!\Big(-\tfrac{b^\star - \gamma - L/\lambda_b}{\lambda_b^{-1/2}}\Big)
+
\Phi\!\Big(-\tfrac{b^\star + \gamma - L/\lambda_b}{\lambda_b^{-1/2}}\Big)
\right],
}
\]

\(L := \sum_i \max(y_i, n_i - y_i) \le N\), \(0 < y_i < n_i\) (no
separation).

*Proof.* Lemmas T, G, Z with \(L_j = L\); same proof as Proposition 1 in
`LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md`. \(\blacksquare\)

Hence \(\pi(A^c \mid \gamma, y) \le G(\gamma)\) whenever \(A^c \subseteq E\)
(e.g.\ when \(B \subseteq A\)).

### 2.3 Multiple groups (\(J \ge 1\)) — sufficient box bound

**Inner box (chosen by the modeller; §0.4).** Choose per-group margins
\(b_j^\star > 0\) such that \(|\eta_{j,i}(\beta_j)| \le b_j^\star\) for all
\(i \in j\) whenever \(|\beta_j| \le b_j^\star\) (offsets fixed). Pick the
largest uniform margin \(b^\star := \min_j b_j^\star\) and floor
\(\kappa_j := N_j\, p(b^\star)(1-p(b^\star))\), \(N_j = \sum_{i\in j} n_{j,i}\).
**Certify** your chosen box
\[
B \;:=\; \prod_{j=1}^{J} [-b^\star, b^\star] \;\subseteq\; A
\]
by verifying (sufficient condition)
\[
\frac{1}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + \kappa_j}
\;<\;
c_\tau
\]
(with \(\kappa_j\) the min curvature of \(W_j\) on \(B\); see §0.4 for the
\(\bar\omega\) logic).

**Proposition L1 (multi-group logit, sufficient escape bound).**

\[
\pi(A^c \mid \gamma, y)
\;\le\;
\pi(B^c \mid \gamma, y)
\;=\;
1 - \prod_{j=1}^{J} \Bigl(1 - \pi(E_j \mid \gamma, y_j)\Bigr)
\;\le\;
\sum_{j=1}^{J} G_j(\gamma),
\]

where \(E_j = \{|\beta_j| > b^\star\}\) and \(G_j\) is the single-group
formula in §2.2 with \((N, L, b^\star) \mapsto (N_j, L_j, b^\star)\),
\(L_j := \sum_{i\in j} \max(y_{j,i}, n_{j,i} - y_{j,i})\).

*Proof.* Factorization §0.1; union bound on \(\{\beta_j \notin [-b^\star,
b^\star]\}\). Each tail bound is Proposition 1 applied group-wise. \(\blacksquare\)

**Gap (open).** A **necessary** closed form for \(\pi(A^c \mid \gamma, y)\)
from \(\bar\omega(\beta) \ge \tau/\lambda_b\) directly — without passing
through a sufficient box — is **not** derived here. The average can fail
when many groups are mildly separated even if no single group hits
\(|b_j| > b^\star\).

---

## 3. Probit (binomial, probit link)

### 3.1 Weights and rate region

\(w_{j,i} = n_{j,i}\, \varphi(\eta_{j,i})^2 / \bigl[\Phi(\eta_{j,i})
\Phi(-\eta_{j,i})\bigr]^2\), bounded by \(2n_{j,i}/\pi\). **Mostly lower-tail**
mixing risk: \(W_j \to 0\) as \(\eta_{j,i} \to -\infty\); upper tail is
bounded. Rate \(A^c\) is again defined by \(\bar\omega(\beta) \ge \tau/\lambda_b\).

### 3.2 Score bound

**Lemma P1 (probit score).** For binomial probit with \(0 < y_i < n_i\),

\[
|\ell_j'(\beta_j)| \;\le\; L_j^{\mathrm{pr}}
\;:=\;
\sum_{i\in j} n_{j,i}\, \frac{\varphi(0)}{\Phi(0)\,\Phi(-0)}
\;=\;
\sum_{i\in j} n_{j,i}\, \frac{2}{\pi}
\;\le\;
\frac{2 N_j}{\pi},
\]

since \(|y_i - n_i \Phi(\eta)| \le n_i\) and \(\varphi(\eta)/[\Phi(\eta)
\Phi(-\eta)] \le 2/\pi\) for all \(\eta\) (Mill's ratio bound at \(\eta=0\)
is conservative).

*Proof sketch.* Optimize \(|y - n\Phi(\eta)|\) over \(y \in (0,n)\) and
bound \(\varphi/(\Phi\Phi)\). \(\blacksquare\)

### 3.3 Closed-form sufficient bound

Choose a **lower floor** \(\eta_j^{\mathrm{floor}}\) such that
\(w_{j,i} \ge w_{\min} > 0\) whenever \(\eta_{j,i} \ge \eta_j^{\mathrm{floor}}\).
With common intercept RE, take \(\beta_j \ge \beta^{\mathrm{floor}}\) for all
\(j\) on a half-space box \(B\). Certify \(B \subseteq A\) via

\[
\frac{1}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + n_j w_{\min}}
\;<\;
c_\tau.
\]

**Proposition P2.** On the one-sided tail \(E_j = \{\beta_j < \beta^{\mathrm{floor}}\}\),

\[
\pi(E_j \mid \gamma, y_j) \;\le\; G_j^{\mathrm{pr}}(\gamma)
\]

with the same tilted-Gaussian template as §2.2, replacing \((L, b^\star)\) by
\((L_j^{\mathrm{pr}}, \beta^{\mathrm{floor}})\) and a **single lower tail**
(only the first \(\Phi\)-term in \(G\)).

**Proposition P3 (multi-group sufficient).**

\[
\pi(A^c \mid \gamma, y) \;\le\; \sum_{j=1}^{J} G_j^{\mathrm{pr}}(\gamma)
\]

when \(A^c \subseteq \bigcup_j \{\beta_j < \beta^{\mathrm{floor}}\}\) is
certified from the floor construction.

*Status.* Upper-tail separation is bounded for probit, so \(A^c\) from
\(\bar\omega\) is not purely one-sided; a two-sided sufficient box
\([\beta^{\mathrm{floor}}, \beta^{\mathrm{ceil}}]^J\) can be used when both
floors are calibrated (§IV.3). Necessary bounds without a box are open.

---

## 4. Poisson (log link)

### 4.1 Weights and rate region

\(w_{j,i} = e^{\eta_{j,i}}\). **Lower-tail** rate problem: \(W_j \to 0\) as
\(\eta_{j,i} \to -\infty\), so \(\omega_j \uparrow\) and \(\bar\omega\) can
cross the threshold. Rate unsafe region: \(\bar\omega(\beta) \ge \tau/\lambda_b\).

### 4.2 Score bound on a floor strip

Fix \(\eta_j^{\mathrm{floor}} \in \mathbb R\) for each group \(j\). On
\(\{\eta_{j,i} \ge \eta_j^{\mathrm{floor}}\ \forall i \in j\}\),

\[
|y_{j,i} - e^{\eta_{j,i}}| \;\le\; \max\bigl(\max_i y_{j,i},\, e^{\eta_j^{\mathrm{floor}}}\bigr)
\;=: M_j^{\mathrm{po}},
\]

so \(|\ell_j'(\beta_j)| \le n_j M_j^{\mathrm{po}}\) on that strip (sum of
absolute Poisson scores).

### 4.3 Closed-form sufficient bound

Choose \(s_j > 0\) with \(\log s_j = \eta_j^{\mathrm{floor}}\) and certify

\[
\frac{1}{J}\sum_{j=1}^{J} \frac{1}{\lambda_b + n_j s_j} \;<\; c_\tau
\quad\Longrightarrow\quad
B = \prod_{j=1}^{J} [\log s_j, \infty) \;\subseteq\; A.
\]

(On each group, \(e^{\eta_{j,i}} \ge s_j\) for all \(i \in j\) when
\(\beta_j \ge \log s_j\) and offsets are fixed.)

**Proposition Po1.** On \(E_j = \{\beta_j < \log s_j\}\),

\[
\pi(E_j \mid \gamma, y_j) \;\le\; G_j^{\mathrm{po}}(\gamma)
\]

with tilted-Gaussian template, \(L_j = n_j M_j^{\mathrm{po}}\), lower tail
only, threshold \(t_j = \log s_j\).

**Remark (common floor).** When \(n_j \equiv n\) and the design is
exchangeable across groups (\(z_{j,i}=1\), equal offsets per group), one may
take \(s_j \equiv s\). Then the certification reduces to
\((1/J)\sum_j 1/(\lambda_b + n s) = 1/(\lambda_b + n s) < c_\tau\) and
\(B = \{\beta : \beta_j \ge \log s\ \forall j\}\), as in the single-margin
special case.

**Proposition Po2 (multi-group sufficient).**

\[
\pi(A^c \mid \gamma, y) \;\le\; \sum_{j=1}^{J} G_j^{\mathrm{po}}(\gamma)
\]

under the floor certification above.

**Remark.** Without a floor, the Poisson score is unbounded on all of
\(\mathbb R\) (Lemma T fails); the floor strip supplies a finite \(L_j\).
This is structural, not a proof gap.

---

## 5. Complementary log-log (binomial, cloglog link)

### 5.1 Weights and rate region

\(w_{j,i} = n_{j,i}\, e^{2\eta_{j,i}}(1-\mu_{j,i})/\mu_{j,i}\) with
\(\mu = 1 - e^{-\eta}\). Lower tail: \(W_j \to 0\) as \(\eta \to -\infty\)
(rate \(A^c\)); same bounding template as Poisson §4 with cloglog weights.

### 5.2 Score bound on a floor strip

On \(\eta_{j,i} \ge \eta_j^{\mathrm{floor}}\), \(\mu_{j,i} \ge
\mu_{\min,j} := 1 - e^{-\eta_j^{\mathrm{floor}}}\), and

\[
w_{j,i} \;\le\; n_{j,i}\, \frac{e^{2\eta_j^{\mathrm{floor}}}}{\mu_{\min,j}},
\qquad
|y_{j,i} - n_{j,i}\mu_{j,i}| \;\le\; n_{j,i}.
\]

Hence \(|\ell_j'(\beta_j)| \le L_j^{\mathrm{cl}} := n_j \max(1, \mu_{\min,j}^{-1}
\cdot e^{2\eta_j^{\mathrm{floor}}})\) on the strip (conservative product bound).

### 5.3 Closed-form sufficient bound

Pick \(s_j > 0\) with \(s_j \le \mu_{\min,j}\) (e.g.\ \(s_j = \mu_{\min,j}\)
when \(\eta_j^{\mathrm{floor}} = \log s_j\) on \(\{\eta_{j,i} \ge \log s_j\}\))
and certify the Poisson-style floor
\((1/J)\sum_j 1/(\lambda_b + n_j s_j) < c_\tau\). **Proposition C1** —
identical template to §4.3 with \(L_j^{\mathrm{cl}}\) and lower tail at
\(\beta_j = \log s_j\). **Special case:** common \(s_j \equiv s\) when groups
are exchangeable, as in §4.3.

---

## 6. Gamma (log link, constant dispersion)

### 6.1 Weights and rate region

With constant \(\phi\), \(w_{j,i} = 1/\phi\), \(W_j = n_j/\phi\): same as
Normal §1.

### 6.2 Closed-form bound

**Proposition G1.** Identical to **Proposition N1** with \(\lambda^\star\)
from §1.4 / §6 of the \(P\)-matrix note (Gamma log with fixed \(\phi\)).

---

## 7. Roadmap (§1–§6 \(\to\) §8–§11)

| Stage | Section | Output |
|---|---|---|
| Rate safe set + \(g(\gamma)\) | §0 | \(A\), \(A^c\), \(g(\gamma)=\pi(A^c\mid\gamma,y)\) |
| Likelihood-specific tails | §1–§6 | Certified \(B\subseteq A\); \(\bar G(\gamma)\) with \(g\le \bar G\) |
| \(\gamma\)-drift | §8 | \(\lambda_{\mathrm{drift}}\), \(C_{\mathrm{drift}}\) |
| Escape mass | §9 | \(P^n(x,A^c)\le U_n(x)\), \(\pi(A^c)\le U_\infty\) (or \(E[\bar G]\)) |
| Gap | §10 | \(\|\pi(A^c)-P^n(x,A^c)\|\) via coupling |
| Full TV | §11 | On-\(A\) + escape terms |

Consolidated likelihood table: §12.

---

## 8. Foster–Lyapunov drift on the \(\gamma\)-chain

All statements below are about the **\(\gamma\)-marginal** chain (block
updated **second**; §0.1). Sweep \(n\): draw \(\beta_{j,n} \mid \gamma_{n-1}\)
for each group, then \(\gamma_n \mid \beta_n\). Initialize \(\gamma_0 = x\).

### 8.1 Exact Gaussian \(\gamma\)-update

Given \(\beta_n = (\beta_{1,n},\ldots,\beta_{J,n})\), conjugacy gives

\[
\gamma_n \mid \beta_n \;\sim\;
N\bigl(w\,\bar\beta_n + (1-w)\mu_0,\;\sigma_\gamma^2\bigr),
\qquad
\bar\beta_n := \frac{1}{J}\sum_{j=1}^{J}\beta_{j,n},
\]

\[
w := \frac{J\lambda_b}{J\lambda_b + \lambda_\gamma}\in(0,1),
\qquad
\sigma_\gamma^2 := (J\lambda_b + \lambda_\gamma)^{-1}.
\]

(\(J=1\): \(w = \lambda_b/(\lambda_b+\lambda_\gamma)\).)

### 8.2 Brascamp–Lieb and mean-bias bounds

**Lemma BL.** For each \(j\), \(\mathrm{Var}(\beta_j \mid \gamma, y) \le
\lambda_b^{-1}\) (log-concavity; prior precision floor \(\lambda_b\)).

**Lemma Bias.** With \(\delta_j(\gamma) := E[\beta_j \mid \gamma, y] - \gamma\),

\[
|\delta_j(\gamma)| \;\le\; \frac{L_j}{\lambda_b},
\qquad
\Bigl|\,E[\bar\beta \mid \gamma, y] - \gamma\,\Bigr|
\;\le\;
\frac{\bar L}{\lambda_b},
\qquad
\bar L := \frac{1}{J}\sum_{j=1}^{J} L_j,
\]

where \(L_j\) is the likelihood score bound from §1–§6 (Lemma T, P1, Poisson
strip, etc.).

*Proof.* Integration by parts on each group, as in
`LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §4 Lemmas 5–6, with
\(\lambda_b^{-1}\) in place of RE variance. \(\blacksquare\)

### 8.3 Drift condition

Let \(V(\gamma) := (\gamma - \mu_0)^2 + 1\). Then for every
\(\gamma_{n-1} \in \mathbb R\),

\[
\boxed{
E\bigl[V(\gamma_n) \mid \gamma_{n-1}\bigr]
\;\le\;
\lambda_{\mathrm{drift}}\, V(\gamma_{n-1}) + C_{\mathrm{drift}}
}
\qquad
\lambda_{\mathrm{drift}} := \frac{1+w^2}{2} < 1,
\]

\[
C_{\mathrm{drift}}
:=
w^2\Bigl(\frac{\bar L}{\lambda_b}\Bigr)^{\!2}
\Bigl(1 + \frac{2w^2}{1-w^2}\Bigr)
+ \sigma_\gamma^2 + w^2\,\frac{\lambda_b^{-1}}{J} + 1 - \lambda_{\mathrm{drift}}.
\]

*Proof.* Compose \(\beta_n \mid \gamma_{n-1}\) (Lemmas BL–Bias) with
\(\gamma_n \mid \beta_n\) (§8.1); Young's inequality on the mean shift and
law of total variance on the variance, as in Proposition 2 of
`LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §4. \(\blacksquare\)

**Remark.** §8 uses only \(\ell_j'' \le 0\) and \(|\ell_j'| \le L_j\); the
likelihood enters through \(\bar L\). Normal/Gamma need no score bound for
drift (any finite \(L_j\) suffices if \(A^c = \varnothing\)).

---

## 9. From drift to \(P^n(x, A^c)\) and \(\pi(A^c)\)

Fix a **certified inner box** \(B \subseteq A\) (§0.4; likelihood-specific
choices in §1–§6). Write \(B = \prod_{j=1}^{J} B_j\) with intervals
\(B_j \subset \mathbb R\). Because \(A^c \subseteq B^c\),

\[
g(\gamma) = \pi(A^c \mid \gamma, y) \;\le\; \pi(B^c \mid \gamma, y).
\]

Define per-margin constants

\[
K_j^\star \;:=\; \inf_{\beta_j \notin B_j} \bigl(\beta_j - \mu_0\bigr)^2 + 1
\]

(the nearest escape of coordinate \(j\) from \(\mu_0\) across \(B_j^c\)).

### 9.1 Markov bound (drift route)

**Lemma M (pointwise escape).** For every \(\gamma\),

\[
g(\gamma)
\;\le\;
\pi(B^c \mid \gamma, y)
\;\le\;
\sum_{j=1}^{J}
\frac{
2V(\gamma) + 2\bigl(L_j/\lambda_b\bigr)^2 + \lambda_b^{-1}
}{
K_j^\star
}
\;=:\;
\Psi_B(\gamma).
\]

*Proof.* For each \(j\), \(\pi(\beta_j \notin B_j \mid \gamma, y_j) \le
E[(\beta_j-\mu_0)^2 \mid \gamma]/K_j^\star\) (Markov on \(B_j^c\)). Lemma
BL–Bias gives \(E[(\beta_j-\mu_0)^2 \mid \gamma] \le 2V(\gamma) +
2(L_j/\lambda_b)^2 + \lambda_b^{-1}\). Sum over \(j\); union bound on
\(B^c\). \(\blacksquare\)

Let \(P_\gamma^n(x,\cdot)\) be the \(n\)-step kernel of the \(\gamma\)-chain
started at \(\gamma_0 = x\). After sweep \(n\), \(\beta_n \mid \gamma_{n-1}\)
is drawn from \(\pi(\beta \mid \gamma_{n-1}, y)\), so

\[
P^n(x, A^c)
=
E\bigl[g(\gamma_{n-1}) \mid \gamma_0 = x\bigr]
\;\le\;
E\bigl[\Psi_B(\gamma_{n-1}) \mid \gamma_0 = x\bigr].
\]

**Proposition Esc-n (started at \(x\)).**

\[
P^n(x, A^c)
\;\le\;
U_n(x)
:=
\sum_{j=1}^{J}
\frac{
2\lambda_{\mathrm{drift}}^{\,n-1} V(x)
+ 2C_{\mathrm{drift}}/(1-\lambda_{\mathrm{drift}})
+ 2(L_j/\lambda_b)^2 + \lambda_b^{-1}
}{
K_j^\star
}.
\]

*Proof.* For each \(j\), iterate §8.3:
\(E[2V(\gamma_{n-1})\mid x] \le 2\lambda_{\mathrm{drift}}^{n-1}V(x) +
2C_{\mathrm{drift}}/(1-\lambda_{\mathrm{drift}})\). Substitute into Lemma M
and sum. \(\blacksquare\)

**Proposition Esc-\(\infty\) (stationary average).** Under stationarity of
the \(\gamma\)-chain,

\[
\pi(A^c)
=
E_{\pi_\gamma}\bigl[g(\gamma)\bigr]
\;\le\;
E_{\pi_\gamma}\bigl[\Psi_B(\gamma)\bigr]
\;\le\;
U_\infty
:=
\sum_{j=1}^{J}
\frac{
2C_{\mathrm{drift}}/(1-\lambda_{\mathrm{drift}})
+ 2(L_j/\lambda_b)^2 + \lambda_b^{-1}
}{
K_j^\star
}
\;=\;
\lim_{n\to\infty} U_n(x).
\]

*Proof.* Meyn–Tweedie (1993, Ch. 14): \(E_{\pi_\gamma}[V] \le
C_{\mathrm{drift}}/(1-\lambda_{\mathrm{drift}})\). \(\blacksquare\)

### 9.2 Tilted-Gaussian bound (§1–§6 route)

When §2–§5 supply \(\bar G(\gamma) \ge \pi(B^c \mid \gamma, y) \ge g(\gamma)\),

\[
\pi(A^c)
=
E_{\pi_\gamma}[g(\gamma)]
\;\le\;
E_{\pi_\gamma}\bigl[\bar G(\gamma)\bigr].
\]

Examples: \(\bar G = G\) (logit \(J=1\), §2.2), \(\bar G = \sum_j G_j\)
(multi-group logit §2.3, probit §3, Poisson §4, cloglog §5). Use whichever
of \(U_\infty\) and \(E[\bar G]\) is smaller for \(\pi(A^c)\); for
\(P^n(x,A^c)\) the drift route §9.1 is the default (tilted bounds need not
integrate the chain state).

### 9.3 Per-likelihood inputs for §9

| Likelihood | Need \(A^c = \varnothing\)? | Certified \(B\) | \(L_j\) | \(K_j^\star\) |
|---|---|---|---|---|
| **Normal** | Yes (§1): \(g \equiv 0\) | — | — | — |
| **Gamma (log)** | Yes (§6): \(g \equiv 0\) | — | — | — |
| **Logit** | No | \(\prod_j [-b^\star,b^\star]\) (§2.3) | §2.2 | \(\min\bigl((b^\star-\mu_0)^2,(b^\star+\mu_0)^2\bigr)+1\) if \(\mu_0 \in B_j\); else nearest boundary |
| **Probit** | No | \(\prod_j [\beta^{\mathrm{floor}},\infty)\) or two-sided (§3) | \(L_j^{\mathrm{pr}}\) | \((\beta^{\mathrm{floor}}-\mu_0)^2+1\) (lower tail); add upper margin if two-sided |
| **Poisson** | No | \(\prod_j [\log s_j,\infty)\) (§4) | \(n_j M_j^{\mathrm{po}}\) | \((\log s_j - \mu_0)^2+1\) when \(\mu_0 \ge \log s_j\) |
| **Cloglog** | No | \(\prod_j [\log s_j,\infty)\) (§5) | \(L_j^{\mathrm{cl}}\) | same as Poisson |

**Logit \(J=1\).** §9.1 and §2.2 can be used together:
\(P^n(x,A^c) \le U_n(x)\) and \(\pi(A^c) \le \min\{U_\infty,\, E_{\pi_\gamma}[G(\gamma)]\}\).

**Normal / Gamma.** Skip §9: \(P^n(x,A^c) = \pi(A^c) = 0\) when
\(\lambda^\star < \tau\).

---

## 10. Bounding \(\big|\pi(A^c) - P^n(x, A^c)\big|\)

The naive bound \(\|g\|_\infty \,\|P_\gamma^{n-1}(x,\cdot)-\pi_\gamma\|_{TV}\)
fails because \(g(\gamma) \to 1\) as \(\gamma \to \pm\infty\) for GLMs. Use
the same coupling on the \(\gamma\)-chain as
`LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §6.

### 10.1 Coupling and minorization

Let \((\gamma_n^x)\) start at \(\gamma_0^x = x\) and \((\gamma_n^\pi)\) start
from \(\gamma_0^\pi \sim \pi_\gamma\) (stationary copy).

**Functions (fixed for this subsection).** Write \(\beta \in \mathbb R^J\).
For each current state \(\gamma \in \mathbb R\), define the **group
posterior** \(\pi_\gamma(\beta) := \pi(\beta \mid \gamma, y)\). Define

\[
m(\beta)
\;:=\;
w\bar\beta + (1-w)\mu_0,
\qquad
\bar\beta := \tfrac{1}{J}\sum_{j=1}^{J}\beta_j,
\]

(the mean of the Gaussian \(\gamma\)-update given \(\beta\); §8.1), and the
**update kernel**

\[
\kappa(\gamma', \beta)
\;:=\;
\phi_q\!\bigl(\gamma';\, m(\beta),\,\Sigma_{\mathrm{upd}}\bigr),
\qquad
\gamma' \in \mathbb R^q,\;\beta \in \mathbb R^J,
\]

with \(\Sigma_{\mathrm{upd}}\) the **fixed** population-update covariance from the
hierarchy (scalar §8.1: \(\Sigma_{\mathrm{upd}}=\sigma_\gamma^2\); multivariate
§10.7). This factor is **always multivariate Normal** in \(\gamma'\) regardless of
the group likelihood (logit, Poisson, etc.). The one-step
**\(\gamma\)-kernel** is the expectation of \(\kappa\) over that draw:

\[
q(\gamma', \gamma)
\;:=\;
\int_{\mathbb R^J} \kappa(\gamma', \beta)\,\pi_\gamma(\beta)\,d\beta
\;=\;
E_{\beta \sim \pi_\gamma}\bigl[\kappa(\gamma', \beta)\bigr]
\;\equiv\;
q_\gamma(\gamma' \mid \gamma).
\]

(One sweep: \(\beta \mid \gamma\), then \(\gamma' \mid \beta \sim
N(m(\beta), \Sigma_{\mathrm{upd}})\).)

Fix \(d > 2C_{\mathrm{drift}}/(1-\lambda_{\mathrm{drift}})\) and the **small
set** (current state)

\[
C_d := \{\gamma : V(\gamma) \le d\}
= [\mu_0 - \sqrt{d-1},\,\mu_0 + \sqrt{d-1}].
\]

**Notation.** Write \(\varepsilon_d := \exp(\min_{\gamma\in C_d}\delta(\gamma))\) for
the **certified minorization constant** on the \(d\)-level set \(C_d\) (§10.3;
`MINORIZATION_GAUSSIAN_REFRESH.md` §5.1.2). Equivalently
\(q(\gamma',\gamma)\ge \varepsilon_d\, q_Q(\gamma')\) for all \(\gamma\in C_d\).

Choose a **refresh density** \(q_Q\) (probability measure \(Q\) on
\(\mathbb R\)) for the merged \(\gamma'\)-draw. The default in this note is
**Gaussian**

\[
Q \;:=\; N(c,\,\tau^2),
\qquad
q_Q(\gamma') \;:=\; \phi(\gamma';\, c,\, \tau^2),
\]

with fixed \((c,\tau)\) or \((c,\Sigma_Q)\) (§10.3, §10.7). Rosenthal allows any \(Q\); Gaussian \(Q\)
is chosen so that \(g(\cdot,\gamma)=\log q(\cdot,\gamma)-\log q_Q\) is
**convex** in \(\gamma'\) when \(\Sigma \succ \Sigma_Q\) on the Normal template
(§10.5–§10.7). If

\[
q(\gamma', \gamma) \;\ge\; \varepsilon_d\, q_Q(\gamma')
\qquad
\forall\,\gamma \in C_d,\;\gamma' \in \mathbb R^q
\]

(scalar: \(q=1\), \(\mathbb R^1 \cong \mathbb R\); multivariate: §10.7),

then \(q(\gamma,\cdot) \ge \varepsilon_d Q(\cdot)\) on \(C_d\) (standard
minorization). Couple the two \(\gamma\)-chains as in
`LOGIT_SINGLE_GROUP_SAFE_REGION.md` §4.6: run independently until both visit
\(C_d\), then merge with probability \(\varepsilon_d\) using \(Q\); \(T\) =
coupling time. On \(\{T \le n-1\}\), draw a single
\(\beta_n \sim \pi(\cdot \mid \gamma_{n-1}^x, y)\) for both copies.

*Remark (uniform refresh).* A uniform \(Q\) on an interval is valid but
**not** recommended here: \(\log q_Q\) is flat on the support, so
\(g(\cdot,\gamma)\) inherits the **concavity** of \(\log q\) and the inner
minimum sits at **endpoints**, not at an interior FOC. Gaussian \(Q\) with
\(\tau^2 < V\) makes \(g\) **convex** and \(\gamma'_*(\gamma)\) an interior
critical point (§10.3).

**Theorem Coupling (Rosenthal 1995).** For any \(0 < r < 1\),

\[
P(T > k)
\;\le\;
R_k
:=
(1-\varepsilon_d)^{rk}
+
\alpha^{-k}(\alpha\Lambda)^{rk}
\Bigl[1 + \frac{C_{\mathrm{drift}}}{1-\lambda_{\mathrm{drift}}} + V(x)\Bigr],
\]

\[
\alpha^{-1} := \frac{1 + 2C_{\mathrm{drift}} + \lambda_{\mathrm{drift}} d}{1+d} < 1,
\qquad
\Lambda := 1 + 2(\lambda_{\mathrm{drift}} d + C_{\mathrm{drift}}).
\]

*Existence.* For each \(\gamma \in C_d\), if \(g(\cdot,\gamma)\) is strictly
convex and coercive on \(\mathbb R\) (Normal with \(\tau^2 < V\): §10.5), the
inner minimum \(\gamma'_*(\gamma)\) exists uniquely; \(\delta(\gamma)\) is
continuous in \(\gamma\); the outer minimum over compact \(C_d\) is attained.
**Value:** §10.3–§10.4.

### 10.2 Cauchy–Schwarz gap bound

**Proposition Gap.**

\[
\big|\pi(A^c) - P^n(x, A^c)\big|
\;\le\;
\Bigl(\sqrt{P^n(x,A^c)} + \sqrt{\pi(A^c)}\Bigr)\,
\sqrt{R_{n-1}}
\;\le\;
\Bigl(\sqrt{U_n(x)} + \sqrt{U_\infty}\Bigr)\,
\sqrt{R_{n-1}}.
\]

*Proof.* Identical to Proposition 4 of
`LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` §6, with \(g(\gamma) =
\pi(A^c \mid \gamma, y)\) and §9 upper bounds on \(P^n(x,A^c)\),
\(\pi(A^c)\). \(\blacksquare\)

**Remark.** \(R_{n-1} \to 0\) as \(n \to \infty\); \(U_n(x)\), \(U_\infty
\to 0\) as margins tighten (\(b^\star \uparrow\), \(s_j \downarrow\) toward
the rate boundary subject to \(B \subseteq A\)). The gap product vanishes in
\(n\) and shrinks as the certified box grows.

### 10.3 Computing \(\varepsilon_d\)

**Refresh density (default).** Fix

\[
Q = N(c,\,\tau^2),
\qquad
q_Q(\gamma') = \phi(\gamma';\, c,\, \tau^2),
\]

(scalar \(q=1\); **multivariate** \(\gamma' \in \mathbb R^q\):
§10.7). **Coupling parameters** \((c,\tau)\) or \((c,\Sigma_Q)\) are chosen once
for all \(\gamma \in C_d\) (§10.5–§10.7). On the Normal template, require
\(\tau^2 < V\) (scalar) or \(\Sigma \succ \Sigma_Q\) (multivariate Loewner:
kernel covariance **more diffuse** than \(Q\); equivalently **smaller precision**
\(-\nabla_{\gamma'}^2 \log q \prec \Sigma_Q^{-1}\); §10.7) so that

\[
\partial_{\gamma'}^2 g(\gamma', \gamma)
= \frac{1}{\tau^2} - \frac{1}{V} \;>\; 0
\qquad\text{(scalar)},
\]

or \(\nabla_{\gamma'}^2 g = \Sigma_Q^{-1} - \Sigma^{-1} \succ 0\) (multivariate),

i.e.\ \(g(\cdot,\gamma)\) is **strictly convex** in \(\gamma'\) and any
interior critical point is the **global** inner minimizer.

**Two-stage minimization (global minorization).** For each \(\gamma \in C_d\),
define the **log gap**

\[
g(\gamma', \gamma)
\;:=\;
\log q(\gamma', \gamma) - \log q_Q(\gamma').
\]

**Inner problem (worst \(\gamma'\) at fixed \(\gamma\)).** When \(g(\cdot,\gamma)\)
is convex on \(\mathbb R\) (scalar) or \(\mathbb R^q\) (§10.7),

\[
\boxed{
\gamma'_*(\gamma)
\;=\;
\frac{V\,c - \tau^2 M(\gamma)}{V - \tau^2},
\qquad
\partial_{\gamma'} g\bigl(\gamma'_*(\gamma), \gamma\bigr) = 0,
}
\]

(scalar Normal; multivariate FOC in §10.7). Here \(M(\gamma) := E[\gamma' \mid \gamma]\)
(§10.5–§10.7; for general GLMs solve \(\nabla_{\gamma'} g = 0\) numerically). Write

\[
\delta(\gamma)
\;:=\;
g\bigl(\gamma'_*(\gamma), \gamma\bigr)
\;=\;
\log q\bigl(\gamma'_*(\gamma), \gamma\bigr) - \log q_Q\bigl(\gamma'_*(\gamma)\bigr).
\]

**Outer problem (worst \(\gamma\) on the small set).**

\[
\boxed{
\gamma_\star
\;\in\;
\arg\min_{\gamma \in C_d} \delta(\gamma),
\qquad
\log \varepsilon_d
\;=\;
\delta(\gamma_\star)
\;=\;
\min_{\gamma \in C_d,\;\gamma' \in \mathbb R^q}
\bigl[\log q(\gamma', \gamma) - \log q_Q(\gamma')\bigr].
}
\]

(\(q=1\): \(\mathbb R^1 \cong \mathbb R\).)

Equivalently (monotonicity of \(\log\)),

\[
\boxed{
\varepsilon_d
\;=\;
\frac{q\bigl(\gamma'_*(\gamma_\star),\, \gamma_\star\bigr)}{q_Q\bigl(\gamma'_*(\gamma_\star)\bigr)}
\;=\;
\exp\bigl(\delta(\gamma_\star)\bigr).
}
\]

Convexity of \(g(\cdot,\gamma)\) implies \(g(\gamma',\gamma) \ge \delta(\gamma)\)
for **every** \(\gamma' \in \mathbb R^q\), so the binding pair
\((\gamma_\star, \gamma'_*(\gamma_\star))\) certifies §10.1 on all of
\(C_d \times \mathbb R^q\). (Do **not** use \(I\) for an identity matrix in
interval notation — package \(I\) is the identity matrix.)

**Choosing \((c,\tau)\).** **Preferred:** \(\tau^2 = \sigma_\gamma^2\) (§8.1),
\(V = \mathrm{Var}(\gamma' \mid \gamma)\) (§10.5), and
\(c = E[\gamma \mid y] = E_{\pi_\gamma}[\gamma]\) (§10.7 **Choosing the mean
\(c\)**). Do **not** set \(c = M(\gamma)\) (state-dependent). Prior mean
\(\mu_0\) is **not** the default center for \(Q\) (drift uses \(\mu_0\); refresh
uses the posterior).

**Numerical evaluation.**

| Step | Action |
|---|---|
| 1 | Fix \(d\), \(C_d\), \((c,\tau)\), and likelihood inputs (§1–§6). |
| 2 | For each \(\gamma \in C_d\): solve \(\partial_{\gamma'} g = 0\) for \(\gamma'_*(\gamma)\); record \(\delta(\gamma)\). |
| 3 | Minimize \(\delta(\gamma)\) over \(\gamma \in C_d\); set \(\varepsilon_d = \exp(\delta(\gamma_\star))\). |
| 4 | Evaluate \(q(\gamma', \gamma)\) via §10.4 (\(J=1\): `integrate()`; \(J\ge1\): quadrature or MC under \(\pi(\beta\mid\gamma,y)\)). |
| 5 | Refine \((\gamma_\star, \gamma'_*(\gamma_\star))\) and \((c,\tau)\) by local optimization from grid starts. |
| 6 | Plug \(\varepsilon_d\) into \(R_{n-1}\); optionally optimize \(d\) and \(r \in (0,1)\). |

**Likelihood notes.**

- **Normal / Gamma** with conjugate group updates: if \(\pi(\beta_j \mid
  \gamma,y_j)\) is Gaussian, \(q(\gamma', \gamma)\) is a Gaussian convolution
  (§10.5: \(\gamma'_\star(\gamma) = E[\gamma' \mid \gamma]\)). Usually
  \(A^c = \varnothing\) anyway (§1, §6).
- **Logit / probit / Poisson / cloglog:** no closed form for the \(\beta\)
  integral; numeric §10.3 is required for a **quantitative** gap bound.
  Drift constants in §8 use score bounds \(L_j\); the kernel integral uses the
  **exact** \(\pi(\beta \mid \gamma,y)\), not the floor-strip approximation.

**If \(\varepsilon_d\) is too small to be useful.** The second term in \(R_k\)
(through \(\alpha^{-1}<1\)) still gives geometric decay in \(k\) without
regeneration; it is typically looser. For applied diagnostics,
`ELLIPSOID_TV_BOUND.md` §2.1 recommends estimating \(P^n(x,A^c)\) from
post-burn-in draws and treating the \(4P^n(x,A^c)\) term as the dominant
empirical piece of §11.

### 10.4 First-order conditions for \(\min q(\gamma', \gamma)\)

**Notation (§10.3 vs joint FOCs).** §10.3 uses \(\gamma'_*(\gamma)\) for the
inner minimizer of \(g(\gamma',\gamma)=\log q - \log q_Q\) and \(\gamma_\star\)
for \(\arg\min_{C_d}\delta(\gamma)\). Block 2 gives
\(\kappa(\gamma',\beta)=\phi_q(\gamma';m(\beta),\Sigma_{\mathrm{upd}})\) always;
set \(\Sigma_Q=\Sigma_{\mathrm{upd}}\). Convexity uses **precision dominance**
\(\Lambda_q := -\nabla_{\gamma'}^2\log q \prec \Sigma_Q^{-1}\) (§10.7); on the
Normal template,
\(\partial_{\gamma'} g = -\Sigma^{-1}(\gamma'-M)+\Sigma_Q^{-1}(\gamma'-c)\). Below,
\((\gamma'_\star,\gamma_\star)\) denotes a **joint** critical point when a uniform
box diagnostic is used; for §10.1 coupling, \(\gamma'_\star = \gamma'_*(\gamma_\star)\)
at the Gaussian inner minimizer.

Fix \(\gamma\) and write \(E_\gamma[\cdot] :=
E_{\beta \sim \pi_\gamma}[\cdot]\). Define the **scalar functionals**
(\(m\) and \(\kappa\) from §10.1)

\[
v(\gamma', \gamma)
\;:=\;
E_\gamma\bigl[\kappa(\gamma', \beta)\bigr]
\;=\;
q(\gamma', \gamma),
\]

\[
u(\gamma', \gamma)
\;:=\;
E_\gamma\bigl[m(\beta)\,\kappa(\gamma', \beta)\bigr].
\]

Both \(u\) and \(v\) are **real-valued functions** of \((\gamma', \gamma)\);
all randomness is integrated out. For fixed \(\gamma\), write \(v_{\gamma}(\gamma')
:= v(\gamma', \gamma)\) and \(u_{\gamma}(\gamma') := u(\gamma', \gamma)\).

**Derivative in \(\gamma'\).** Since
\(\partial_{\gamma'} \kappa(\gamma', \beta)
= \kappa(\gamma', \beta)\,(m(\beta)-\gamma')/\sigma_\gamma^2\),

\[
\partial_{\gamma'} v(\gamma', \gamma)
\;=\;
\frac{1}{\sigma_\gamma^2}\,
\bigl(u(\gamma', \gamma) - \gamma'\, v(\gamma', \gamma)\bigr).
\]

**Interior FOC in \(\gamma'\).** At \((\gamma'_\star, \gamma_\star)\) with
\(v(\gamma'_\star, \gamma_\star) > 0\),

\[
\partial_{\gamma'} v(\gamma'_\star, \gamma_\star) = 0
\quad\Longleftrightarrow\quad
u(\gamma'_\star, \gamma_\star)
\;=\;
\gamma'_\star\, v(\gamma'_\star, \gamma_\star),
\]

equivalently

\[
\boxed{
\gamma'_\star
\;=\;
\frac{u(\gamma'_\star, \gamma_\star)}{v(\gamma'_\star, \gamma_\star)}
\;=\;
\frac{
E_{\gamma_\star}\bigl[m(\beta)\,\kappa(\gamma'_\star, \beta)\bigr]
}{
E_{\gamma_\star}\bigl[\kappa(\gamma'_\star, \beta)\bigr]
}.
}
\]

Define \(\mu_\kappa(\gamma', \gamma) := u(\gamma', \gamma)/v(\gamma', \gamma)\)
when \(v > 0\). Then \(\gamma'_\star = \mu_\kappa(\gamma'_\star, \gamma_\star)\)
(fixed point of \(\gamma' \mapsto \mu_\kappa(\gamma', \gamma_\star)\)).

**Derivative in \(\gamma\).** With prior score function
\(s(\beta, \gamma) := \partial_\gamma \log \phi(\beta \mid \gamma)
= \sum_j \mathcal W_j^\top \Psi_j^{-1}(\beta_j - \mathcal W_j\gamma)\),

\[
\partial_\gamma v(\gamma', \gamma)
\;=\;
E_\gamma\bigl[s(\beta, \gamma)\,\kappa(\gamma', \beta)\bigr]
\;-\;
E_\gamma\bigl[s(\beta, \gamma)\bigr]\, v(\gamma', \gamma).
\]

Write \(w_s(\gamma', \gamma) :=
E_\gamma[s(\beta, \gamma)\,\kappa(\gamma', \beta)]\). The interior FOC
\(\partial_\gamma v(\gamma'_\star, \gamma_\star) = 0\) is

\[
w_s(\gamma'_\star, \gamma_\star)
\;=\;
E_{\gamma_\star}\bigl[s(\beta, \gamma_\star)\bigr]\,
v(\gamma'_\star, \gamma_\star),
\qquad\text{i.e.}\quad
\frac{w_s(\gamma'_\star, \gamma_\star)}{v(\gamma'_\star, \gamma_\star)}
\;=\;
E_{\gamma_\star}\bigl[s(\beta, \gamma_\star)\bigr].
\]

**Minimizing \(\varepsilon_d\).** With Gaussian \(Q\) (§10.3),
\(\varepsilon_d = \exp(\min_{C_d} \delta(\gamma))\). For a **uniform-interval
diagnostic**, \(\varepsilon_{\mathrm{box}} = \ell_Q \min_{(\gamma',\gamma)}
v(\gamma', \gamma)\) on a fixed box remains available (§10.5 below) but is not
the §10.1 refresh measure. For GLMs on a box, \(\gamma'_\star =
\mu_\kappa(\gamma'_\star, \gamma_\star)\) is typically a **non-trivial** fixed
point (§10.5 contrasts the Normal case).

### 10.5 Normal (conjugate) reduction

When each group likelihood is **Normal** (§1) and \(\pi(\beta_j \mid \gamma,y_j)\)
is Gaussian (conjugate update), the general §10.4 machinery simplifies
sharply. **Remark:** if \(\lambda^\star < \tau\) so \(A^c = \varnothing\), escape
terms in §9–§11 vanish anyway (§1); this subsection describes the kernel
\(q(\gamma', \gamma)\) only.

**Marginal kernel.** With \(m(\beta) = w\bar\beta + (1-w)\mu_0\) from §10.1,
\(\beta \mid \gamma, y\) is Gaussian, hence the one-step update satisfies

\[
\gamma' \mid \gamma \;\sim\; N\bigl(M(\gamma),\, \mathrm{Var}(\gamma' \mid \gamma)\bigr),
\qquad
M(\gamma) := E[\gamma' \mid \gamma].
\]

(\(\mathrm{Var}(\gamma' \mid \gamma)\) constant in \(\gamma\) for fixed data and
priors). **Closed form (Normal §1).** Conjugate updates give, for each group,
\(\beta_j \mid \gamma, y_j \sim N(\mu_j(\gamma),\, \tau_j^{-1})\) with data
precision \(W_j := n_j/\phi\) (§1.1) and

\[
\tau_j \;:=\; \lambda_b + W_j,
\qquad
\mathrm{Var}(\beta_j \mid \gamma, y) \;=\; \tau_j^{-1}
\]

(independent across \(j\) given \(\gamma\)). With
\(\bar\beta := J^{-1}\sum_j \beta_j\), \(\sigma_\gamma^2\) and \(w\) from §8.1,
and \(\gamma' \mid \beta \sim N(m(\beta), \sigma_\gamma^2)\),

\[
\boxed{
\mathrm{Var}(\gamma' \mid \gamma)
\;=\;
\sigma_\gamma^2
\;+\;
w^2\,\mathrm{Var}(\bar\beta \mid \gamma, y)
\;=\;
\frac{1}{J\lambda_b + \lambda_\gamma}
\;+\;
\frac{w^2}{J^2}\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}
}.
\]

Equivalently, using \(w = J\lambda_b/(J\lambda_b + \lambda_\gamma)\),

\[
\mathrm{Var}(\gamma' \mid \gamma)
\;=\;
\frac{1}{J\lambda_b + \lambda_\gamma}
\left[
1
\;+\;
\frac{J\,\lambda_b^2}{(J\lambda_b + \lambda_\gamma)\,J}
\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}
\right]
\;=\;
\frac{1}{J\lambda_b + \lambda_\gamma}
\left[
1
\;+\;
\frac{\lambda_b^2}{J\lambda_b + \lambda_\gamma}
\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}
\right].
\]

*Special cases.* \(J=1\): \(w = \lambda_b/(\lambda_b + \lambda_\gamma)\) and
\(\mathrm{Var}(\gamma' \mid \gamma)
= (\lambda_b + \lambda_\gamma)^{-1}
+ w^2/(\lambda_b + W_1)\).
Equal groups \(W_j \equiv W\):
\(\mathrm{Var}(\gamma' \mid \gamma)
= (J\lambda_b + \lambda_\gamma)^{-1}
+ w^2/\bigl(J(\lambda_b + W)\bigr)\).

Therefore

\[
v(\gamma', \gamma) = q(\gamma', \gamma)
= \phi\bigl(\gamma';\, M(\gamma),\, \mathrm{Var}(\gamma' \mid \gamma)\bigr).
\]

The §10.4 fixed point \(\gamma' = \mu_\kappa(\gamma', \gamma)\) collapses to
\(\gamma' = M(\gamma) = E[\gamma' \mid \gamma]\) (mode of \(v(\cdot,\gamma)\)).

**Slice value at \(\gamma' = M(\gamma)\).** The kernel is maximal in \(\gamma'\) at
the conditional mean:

\[
v\bigl(M(\gamma), \gamma\bigr)
= \max_{\gamma'} v(\gamma', \gamma)
= \frac{1}{\sqrt{2\pi\,\mathrm{Var}(\gamma' \mid \gamma)}},
\]

which is **constant in \(\gamma\)** on the interior plateau. (On a bounded interval
\([\gamma'_-,\gamma'_+]\), the **minimum** of \(v(\cdot,\gamma)\) is at an
endpoint; see **Minorization constant** below.)

Hence the **value**
\(\min_{\gamma \in C_d,\,\gamma' \in [\gamma'_-,\gamma'_+]} v(\gamma', \gamma)\)
is determined by \(\mathrm{Var}(\gamma' \mid \gamma)\), endpoint distance from
\(M(\gamma)\), and interval feasibility
\(M(\gamma) \in [\gamma'_-,\gamma'_+]\); **\(\gamma_\star\)** is generally
**not unique** on the interior plateau (any \(\gamma \in C_d\) with
\(M(\gamma)\) in the interval is tied). Unique \((\gamma_\star, \gamma'_\star)\)
appear when **boundaries** of \(C_d\) or \([\gamma'_-,\gamma'_+]\) bind.

**Minorization constant \(\varepsilon_d\).** From §10.3,
\(\varepsilon_d = \exp(\min_{C_d}\delta(\gamma))\).

**Peak plug-in (diagnostic, not \(\varepsilon_d\)).**
\(\varepsilon_{\mathrm{peak}} = \ell_Q \min_{\gamma \in C_d,\,\gamma' \in [\gamma'_-,\gamma'_+]}
v(\gamma', \gamma)\) with \(\ell_Q := \gamma'_+ - \gamma'_-\). Substitute
\(v(\gamma', \gamma) = \phi(\gamma';\, M(\gamma),\, \mathrm{Var}(\gamma' \mid \gamma))\)
and the closed form above. Write

\[
V_\gamma
\;:=\;
\mathrm{Var}(\gamma' \mid \gamma)
\;=\;
\frac{1}{J\lambda_b + \lambda_\gamma}
\;+\;
\frac{w^2}{J^2}\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}
\qquad
\bigl(W_j = n_j/\phi \text{ from §1.1}\bigr).
\]

*Peak evaluation (slice at \(\gamma' = M(\gamma)\)).* At the mode,
\(v(M(\gamma), \gamma) = 1/\sqrt{2\pi V_\gamma}\). The peak plug-in
\(\varepsilon_{\mathrm{peak}}/\ell_Q = v(M(\gamma), \gamma)\) gives

\[
\boxed{
\varepsilon_{\mathrm{peak}}
=
\frac{\ell_Q}{\sqrt{2\pi V_\gamma}}
}
\]

Substituting the closed form for \(V_\gamma\),

\[
\varepsilon_{\mathrm{peak}}
=
\frac{\ell_Q}{
\sqrt{2\pi\left(
\frac{1}{J\lambda_b + \lambda_\gamma}
+ \frac{w^2}{J^2}\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}
\right)}
}.
\]

\[
\boxed{
1 - \varepsilon_{\mathrm{peak}}
=
\frac{\sqrt{2\pi V_\gamma} - \ell_Q}{\sqrt{2\pi V_\gamma}}
}
\]

Write \(D_\gamma := \sqrt{2\pi V_\gamma}\). Then \(1 - \varepsilon_{\mathrm{peak}} = (D_\gamma - \ell_Q)/D_\gamma\),
and substituting \(V_\gamma\),

\[
D_\gamma
=
\sqrt{
2\pi\left(
\frac{1}{J\lambda_b + \lambda_\gamma}
+ \frac{w^2}{J^2}\sum_{j=1}^{J} \frac{1}{\lambda_b + W_j}
\right)
}.
\]

This is tight **at the mode only**; the uniform-interval
\(\varepsilon_{\mathrm{box}} \le \varepsilon_d\) below is an auxiliary bound.
**Coupling (§10.1):** merge probability \(\varepsilon_d = \exp(\delta(\gamma_\star))\)
from Gaussian \(Q\) (§10.3); non-merge \(1-\varepsilon_d\); hence
\(P(T > k) \le (1-\varepsilon_d)^{rk} + \cdots\).

*Certified minimum on a fixed interval.* If \([\gamma'_-,\gamma'_+]\) is not
tight around \(M(\gamma)\), the box minimum is at an endpoint. With
\(\Delta(\gamma) := \max\{M(\gamma)-\gamma'_-,\,\gamma'_+ - M(\gamma)\}\),

\[
\min_{\gamma' \in [\gamma'_-,\gamma'_+]} v(\gamma', \gamma)
=
\frac{1}{\sqrt{2\pi V_\gamma}}
\exp\!\Bigl(-\frac{\Delta(\gamma)^2}{2 V_\gamma}\Bigr),
\]

\[
\boxed{
\varepsilon_{\mathrm{box}}
=
\ell_Q \cdot
\min_{\gamma \in C_d}
\frac{1}{\sqrt{2\pi V_\gamma}}
\exp\!\Bigl(-\frac{\Delta(\gamma)^2}{2 V_\gamma}\Bigr)
}
\]

At any minimizer \(\gamma_\star\),
\(\varepsilon_{\mathrm{box}} \le \ell_Q / \sqrt{2\pi V_\gamma}\).

At the worst-case \(\gamma_\star \in C_d\),

\[
\boxed{
1 - \varepsilon_{\mathrm{box}}
=
\frac{
\sqrt{2\pi V_\gamma}
- \ell_Q \exp\!\bigl(-\Delta(\gamma_\star)^2 / (2 V_\gamma)\bigr)
}{
\sqrt{2\pi V_\gamma}
}
}
\]

With \(D_\gamma := \sqrt{2\pi V_\gamma}\),

\[
1 - \varepsilon_{\mathrm{box}}
=
\frac{
D_\gamma - \ell_Q \exp\!\bigl(-\Delta(\gamma_\star)^2 / (2 V_\gamma)\bigr)
}{
D_\gamma
}.
\]

**Symmetric chord and minorization region (Normal; auxiliary).** *This block
uses a **uniform** \(q_Q\) on an interval for box diagnostics only; §10.1
coupling uses **Gaussian** \(Q\) (§10.3–§10.5 **Gaussian refresh \(Q\)**).*
On Normal (§1),
\(V_\gamma \equiv V\) is **constant in \(\gamma\)** (§10.5 closed form). Fix
\(\gamma \in C_d\), half-width \(h > 0\), and the chord
\[
[\gamma'_-,\gamma'_+]
=
\bigl[M(\gamma)-h,\, M(\gamma)+h\bigr],
\qquad
\ell_Q = 2h,
\qquad
q_Q = \mathrm{Uniform}[\gamma'_-,\gamma'_+].
\]
Then \(v(\gamma',\gamma)=\phi(\gamma';\,M(\gamma),V)\) and

\[
\min_{\gamma' \in [\gamma'_-,\gamma'_+]} v(\gamma', \gamma)
=
\frac{1}{\sqrt{2\pi V}}
\exp\!\Bigl(-\frac{h^2}{2V}\Bigr),
\]

\[
\boxed{
\varepsilon_{\mathrm{box}}(\gamma,h)
=
\ell_Q \cdot \min_{\gamma'} v
=
\frac{2h}{\sqrt{2\pi V}}
\exp\!\Bigl(-\frac{h^2}{2V}\Bigr).
}
\]

**Minorization inequality (§10.1, Normal chord).** With
\(q(\gamma',\gamma)=v(\gamma',\gamma)\) and \(q_Q(\gamma')=1/\ell_Q=1/(2h)\) on
\([\gamma'_-,\gamma'_+]\), the requirement is

\[
\boxed{
q(\gamma', \gamma)
\;\ge\;
\varepsilon_{\mathrm{box}}(\gamma,h)\, q_Q(\gamma')
\qquad
\forall\,\gamma' \in \bigl[M(\gamma)-h,\, M(\gamma)+h\bigr].
}
\]

Substituting \(q=\phi(\cdot;\,M(\gamma),V)\), \(\ell_Q=2h\), and
\(\varepsilon_{\mathrm{box}}(\gamma,h)\) above,

\[
\boxed{
\phi\bigl(\gamma';\, M(\gamma), V\bigr)
\;\ge\;
\frac{1}{2h}
\cdot
\frac{2h}{\sqrt{2\pi V}}
\exp\!\Bigl(-\frac{h^2}{2V}\Bigr)
\;=\;
\frac{1}{\sqrt{2\pi V}}
\exp\!\Bigl(-\frac{h^2}{2V}\Bigr)
\qquad
\forall\,\gamma' \in \bigl[M(\gamma)-h,\, M(\gamma)+h\bigr].
}
\]

Equivalently \(|\gamma' - M(\gamma)| \le h\) (see below).

**Why this \(\varepsilon\) is valid, and where equality holds.** Write
\[
q_{\min}
:=
\min_{\gamma' \in [M(\gamma)-h,\,M(\gamma)+h]} q(\gamma',\gamma)
=
\frac{1}{\sqrt{2\pi V}}
\exp\!\Bigl(-\frac{h^2}{2V}\Bigr)
=
q\bigl(M(\gamma)+h,\,\gamma\bigr)
=
q\bigl(M(\gamma)-h,\,\gamma\bigr).
\]

*Step 1 — definition (§10.3).* Set
\(\varepsilon_{\mathrm{box}}(\gamma,h) := \ell_Q\, q_{\min} = 2h\, q_{\min}\).
Then for **every** \(\gamma'\) on the chord,
\[
q(\gamma',\gamma) \;\ge\; q_{\min}
= \frac{\varepsilon_{\mathrm{box}}(\gamma,h)}{\ell_Q}
= \varepsilon_{\mathrm{box}}(\gamma,h)\, q_Q(\gamma'),
\]
since \(q_Q \equiv 1/\ell_Q\) on the interval. So the minorization inequality
holds **by construction**.

*Step 2 — equality at the endpoints.* At \(\gamma' = M(\gamma) \pm h\),
\[
q(M(\gamma)\pm h,\,\gamma)
=
\frac{1}{\sqrt{2\pi V}}
\exp\!\Bigl(-\frac{h^2}{2V}\Bigr)
=
q_{\min}
=
\varepsilon_{\mathrm{box}}(\gamma,h)\, q_Q(M(\gamma)\pm h).
\]
So the bound is **tight** (equality) at both endpoints.

*Step 3 — strict inequality in the interior.* \(q(\cdot,\gamma)\) is unimodal
with unique maximum at \(\gamma' = M(\gamma)\):
\[
q(M(\gamma),\gamma)
=
\frac{1}{\sqrt{2\pi V}}
\;>\;
\frac{1}{\sqrt{2\pi V}}
\exp\!\Bigl(-\frac{h^2}{2V}\Bigr)
=
\varepsilon_{\mathrm{box}}(\gamma,h)\, q_Q(\gamma')
\qquad (h > 0).
\]

*Step 4 — “small enough” / largest valid constant on the chord.* Any
\(\varepsilon' > \varepsilon_{\mathrm{box}}(\gamma,h)\) would violate the
inequality at \(\gamma' = M(\gamma) \pm h\), where \(q = q_{\min}\). So
\(\varepsilon_{\mathrm{box}}\) is the **largest** \(\varepsilon\) for which
minorization holds on this chord (not an arbitrary conservative fraction of
something larger). It is still a **small** merge probability in absolute terms
(\(\varepsilon_{\mathrm{box}} \le 2/\sqrt{2\pi e} < 1\); see \(h_\star\) below).

*Step 5 — equivalent form.* Since
\(\phi(\gamma';M,V)/\phi(M\pm h;M,V) = \exp\bigl(h^2/(2V) - (\gamma'-M)^2/(2V)\bigr)\),
the inequality \(q \ge q_{\min}\) is equivalent to \(|\gamma' - M(\gamma)| \le h\).

So for **fixed \(\gamma\)**, minorization holds on the whole chord
\(\{\gamma\}\times[M(\gamma)-h,\,M(\gamma)+h]\), with **equality** at
\(\gamma' = M(\gamma)\pm h\) and **strict inequality** in the interior.

**One global interval on \(C_d\).** If a **single** \([\gamma'_-,\gamma'_+]\) is
used for all \(\gamma \in C_d\) (not re-centered at each \(M(\gamma)\)), with
\(\Delta(\gamma) := \max\{M(\gamma)-\gamma'_-,\,\gamma'_+ - M(\gamma)\}\) and
\(\varepsilon_{\mathrm{box}} = \ell_Q \min_{\gamma\in C_d,\,\gamma'\in[\gamma'_-,\gamma'_+]} v\),

\[
\boxed{
\mathcal R
:=
\bigl\{
(\gamma,\gamma') \in C_d \times [\gamma'_-,\gamma'_+]
\;:\;
\phi\bigl(\gamma';\, M(\gamma), V\bigr)
\;\ge\;
\varepsilon_{\mathrm{box}}/\ell_Q
\bigr\}.
}
\]

Minorization holds **on all of** \(C_d \times [\gamma'_-,\gamma'_+]\) iff
\(\mathcal R = C_d \times [\gamma'_-,\gamma'_+]\). In general \(\mathcal R\) is
the set where the inequality is satisfied; the **binding** pairs
\((\gamma_\star,\gamma'_\star)\) lie on \(\partial\mathcal R\) (endpoints in
\(\gamma'\), and possibly \(\gamma_\star \in \partial C_d\)). If the interval is
symmetric about \(c\) with half-width \(H\), then at each \(\gamma\),
\(\varepsilon_{\mathrm{box}}/\ell_Q = (2\pi V)^{-1/2}
\exp\bigl(-(\Delta(\gamma))^2/(2V)\bigr)\) and \(\mathcal R\) is the region
between the level curves \(| \gamma' - M(\gamma)|\) and the threshold set by
\(\Delta(\gamma)\).

**Gaussian refresh \(Q\) for §10.1 coupling (Normal).** With
\(q(\gamma',\gamma)=\phi(\gamma';\,M(\gamma),V)\), \(V_\gamma\equiv V\), and
\(q_Q(\gamma')=\phi(\gamma';\,c,\,\tau^2)\), require \(\tau^2 < V\) so
\(g(\cdot,\gamma)\) is convex. **Canonical choice:** \(\tau^2=\sigma_\gamma^2\)
(§8.1 block-2 conditional variance) and \(V=\mathrm{Var}(\gamma'\mid\gamma)\)
(§10.5); then \(V=\sigma_\gamma^2+w^2\mathrm{Var}(\bar\beta\mid\gamma,y)>\tau^2\).
Then §10.3 gives

\[
\gamma'_*(\gamma)
=
\frac{V\,c - \tau^2 M(\gamma)}{V - \tau^2},
\qquad
\boxed{
\delta(\gamma)
=
\frac{1}{2}\log\frac{\tau^2}{V}
\;-\;
\frac{\bigl(M(\gamma)-c\bigr)^2}{2\,(V+\tau^2)}
}.
\]

**Outer min.** \(\gamma_\star\) maximizes \((M(\gamma)-c)^2\) on \(C_d\)
(typically \(\gamma_\star \in \partial C_d\)). With \(c=m_\gamma\) and
\(M(\gamma) = \mu_0 + w(\gamma-\mu_0) + \text{bias}\),

\[
\boxed{
\varepsilon_d
=
\exp\bigl(\delta(\gamma_\star)\bigr)
=
\sqrt{\frac{\tau^2}{V}}\;
\exp\!\Bigl(
-\frac{\bigl(M(\gamma_\star)-c\bigr)^2}{2\,(V+\tau^2)}
\Bigr).
}
\]

*Template.* Take \(c=m_\gamma\), \(\tau^2 = \rho V\) with \(\rho \in (0,1)\). At
\(M(\gamma)=c\) (rare unless \(M(m_\gamma)=m_\gamma\)), \(\varepsilon(\gamma) = \sqrt{\rho}\).
**Preferred:** \(\tau^2=\sigma_\gamma^2\), \(c=m_\gamma\) (§10.7). Two-block:
\(V > \sigma_\gamma^2\) and \(\varepsilon(\gamma)=\sigma_\gamma/\sqrt{V}\) when
\(M(\gamma)=c\).

*Envelope of means.* On \(C_d\), \(M(\gamma) \in [M_-,M_+]\) with
\(M_\pm = \mu_0 \pm w(\sqrt{d-1}+\bar L/\lambda_b)\) (Lemma Bias). With
\(c=m_\gamma\), the outer penalty uses \(\|M(\gamma_\star)-m_\gamma\|_{(V+\tau^2)}^2\)
when \(\gamma_\star \in \partial C_d\) binds.

*GLM / multivariate.* See §10.7 for vector \(\gamma' \in \mathbb R^q\); GLMs use
the same \(g\) with numeric \(\nabla_{\gamma'} g = 0\).

*Uniform interval bounds (auxiliary).* The chord / \(\varepsilon_{\mathrm{box}}\)
calculations below with \(\mathrm{Uniform}[\gamma'_-,\gamma'_+]\) remain useful
**diagnostics** for box-shaped minorization but are **not** the refresh measure
used in §10.1 coupling.

**Optimal half-width (fixed \(\gamma\), Normal).** Maximizing
\(\varepsilon_{\mathrm{box}}(\gamma,h)\) over \(h > 0\) gives

\[
\boxed{
h_\star = \sqrt{V},
\qquad
\varepsilon_{\max}
=
\frac{2}{\sqrt{2\pi e}}
\approx 0.757,
}
\]

independent of \(V\) at the interior optimum. Require \(\varepsilon_{\mathrm{box}} \le 1\)
(valid merge probability); here \(\varepsilon_{\max} < 1\).

*Equal groups* (\(W_j \equiv W\)):
\(
V_\gamma = (J\lambda_b + \lambda_\gamma)^{-1} + w^2/\bigl(J(\lambda_b + W)\bigr)
\),
so replace the sum by \(J/(\lambda_b + W)\) in every display above.
*Single group* (\(J=1\)):
\(
V_\gamma = (\lambda_b + \lambda_\gamma)^{-1} + w^2/(\lambda_b + W_1)
\).

**Contrast (GLM).** For logit/probit/Poisson/cloglog, \(v(\gamma', \gamma)\) is
not Gaussian in \(\gamma'\); the interior \(\gamma'\)-FOC remains
\(\gamma'_\star = \mu_\kappa(\gamma'_\star, \gamma_\star)\) with
\(\mu_\kappa \neq E[\gamma' \mid \gamma]\) in general, and \(v\) typically
varies with \(\gamma\) on \(C_d\), so both coordinates of the minimizer are
informative.

### 10.6 Logit: half-width FOC with \(q\), \(u\), \(v\) plugged in

**Logit data and likelihood (§2).** Group \(j\) has binomial counts
\(y_{j,i}\), \(n_{j,i}\) with \(0 < y_{j,i} < n_{j,i}\), canonical link, scalar
intercept RE \(\eta_{j,i} = \mathrm{offset}_{j,i} + \beta_j\). Write

\[
p_{j,i}(\beta)
:=
\operatorname{expit}(\eta_{j,i}),
\qquad
w_{j,i}(\beta)
:=
n_{j,i}\, p_{j,i}(\beta)\bigl(1-p_{j,i}(\beta)\bigr),
\qquad
W_j(\beta)
:=
\sum_{i\in j} w_{j,i}(\beta),
\]

\[
\ell_j(\beta_j)
:=
\sum_{i\in j}
\Bigl(
y_{j,i}\,\eta_{j,i}
-
n_{j,i}\,\log\bigl(1+e^{\eta_{j,i}}\bigr)
\Bigr),
\qquad
\ell_j'(\beta_j)
=
\sum_{i\in j}
\bigl(y_{j,i} - n_{j,i}\, p_{j,i}(\beta)\bigr),
\qquad
\ell_j''(\beta_j)
=
-\,W_j(\beta).
\]

The **group posterior** (§0.1) is
\(\pi(\beta_j \mid \gamma, y_j) \propto
\exp(\ell_j(\beta_j))\,\phi(\beta_j;\gamma,\lambda_b^{-1})\).
The weights \(w_{j,i}=n_{j,i}p_{j,i}(1-p_{j,i})\) and \(W_j\) enter **only through
\(\pi_\gamma(\beta)\)** (curvature / mass of \(\beta \mid \gamma,y\)), not through
\(\kappa(\gamma',\beta)\) itself.

Fix \(\gamma \in C_d\) and \(M(\gamma) := E[\gamma' \mid \gamma]\). From §10.1–§10.4,

\[
\kappa(\gamma', \beta)
=
\phi\!\bigl(\gamma';\, m(\beta),\, \sigma_\gamma^2\bigr),
\qquad
m(\beta) = w\bar\beta + (1-w)\mu_0,
\]

\[
q(\gamma', \gamma)
=
v(\gamma', \gamma)
=
E_\gamma\bigl[\kappa(\gamma', \beta)\bigr],
\qquad
u(\gamma', \gamma)
=
E_\gamma\bigl[m(\beta)\,\kappa(\gamma', \beta)\bigr],
\]

with \(E_\gamma[\cdot]\) over \(\beta \sim \pi(\beta \mid \gamma, y)\). Then

\[
q_{\gamma'}(\gamma', \gamma)
=
\frac{1}{\sigma_\gamma^2}
\Bigl[
u(\gamma', \gamma) - \gamma'\, q(\gamma', \gamma)
\Bigr].
\]

**Explicit \(\beta\)-integrals (logit).** With
\(Z_j(\gamma) := \int e^{\ell_j(\beta_j)}\phi(\beta_j;\gamma,\lambda_b^{-1})\,d\beta_j\),

\[
\pi_\gamma(\beta)
=
\prod_{j=1}^{J}
\frac{
e^{\ell_j(\beta_j)}\,\phi(\beta_j;\gamma,\lambda_b^{-1})
}{
Z_j(\gamma)
},
\]

\[
q(\gamma', \gamma)
=
\frac{1}{\prod_j Z_j(\gamma)}
\int
\Bigl(\prod_j e^{\ell_j(\beta_j)}\phi(\beta_j;\gamma,\lambda_b^{-1})\Bigr)
\,\phi\bigl(\gamma';\, m(\beta),\, \sigma_\gamma^2\bigr)\, d\beta,
\]

\[
u(\gamma', \gamma)
=
\frac{1}{\prod_j Z_j(\gamma)}
\int
m(\beta)\,
\Bigl(\prod_j e^{\ell_j(\beta_j)}\phi(\beta_j;\gamma,\lambda_b^{-1})\Bigr)
\,\phi\bigl(\gamma';\, m(\beta),\, \sigma_\gamma^2\bigr)\, d\beta.
\]

Every \(p_{j,i}(\beta)\) and \(W_j(\beta)\) sits inside \(\ell_j(\beta_j)\) (and
\(Z_j\)) in the integrand. **Single group** (\(J=1\)): drop \(j\), write
\(y_i,n_i,p_i(\beta),W(\beta)=\sum_i n_i p_i(\beta)(1-p_i(\beta))\).

**Symmetric chord parametrization.** Candidate interval
\([\gamma'_-, \gamma'_+] = [M(\gamma)-h,\, M(\gamma)+h]\), half-width \(h>0\):

\[
\varepsilon(\gamma, h)
=
2h\, q\bigl(M(\gamma)+h,\, \gamma\bigr).
\]

**Equal endpoint values:**

\[
\boxed{
\frac{
\int
e^{\ell(\beta)}\,\phi(\beta;\gamma,\lambda_b^{-1})\,
\phi\bigl(M(\gamma)-h;\, m(\beta),\, \sigma_\gamma^2\bigr)\, d\beta
}{
\int e^{\ell(\beta)}\,\phi(\beta;\gamma,\lambda_b^{-1})\, d\beta
}
=
\frac{
\int
e^{\ell(\beta)}\,\phi(\beta;\gamma,\lambda_b^{-1})\,
\phi\bigl(M(\gamma)+h;\, m(\beta),\, \sigma_\gamma^2\bigr)\, d\beta
}{
\int e^{\ell(\beta)}\,\phi(\beta;\gamma,\lambda_b^{-1})\, d\beta
}
}
\]

(multi-group: \(\ell(\beta)=\sum_j\ell_j(\beta_j)\), product prior factors as
above; \(J=1\) scalar \(\beta\).)

**Half-width FOC** \(q(M+h,\gamma)+h\,q_{\gamma'}(M+h,\gamma)=0\):

\[
\boxed{
q(M+h,\gamma)
+
\frac{h}{\sigma_\gamma^2}
\Bigl[
u(M+h,\gamma) - \bigl(M(\gamma)+h\bigr)\, q(M+h,\gamma)
\Bigr]
= 0,
}
\qquad
\boxed{
\mu_\kappa\bigl(M(\gamma)+h,\, \gamma\bigr)
=
M(\gamma) + h - \frac{\sigma_\gamma^2}{h},
}
\]

with \(\mu_\kappa=u/v\) and \(q,u\) the logit integrals above (evaluate numerically
via likelihood-subgradient normalization, §10.3).

**Large-\(n\) / Laplace remark.** At \(\hat\beta_j(\gamma)\) with
\(W_j = W_j(\hat\beta_j(\gamma)) = \sum_i n_{j,i}\, p_{j,i}(1-p_{j,i})\) from §2.1,
\(\pi_\gamma\) is approximately
\(N(\hat\beta_j(\gamma), [\lambda_b+W_j]^{-1})\). Then §10.5 applies with
\[
V_\gamma
=
\frac{1}{J\lambda_b+\lambda_\gamma}
+
\frac{w^2}{J^2}\sum_{j=1}^{J}\frac{1}{\lambda_b+W_j(\hat\beta_j(\gamma))},
\]
and \(h \approx \sqrt{V_\gamma}\) in the \(h\)-FOC (§10.5–§10.6). Here
\(n_{j,i}\) and \(p(1-p)\) enter **explicitly** through each \(W_j\).

### 10.7 Multivariate Gaussian refresh (\(\gamma' \in \mathbb R^q\))

This section extends §10.1–§10.3 to **vector** population states
\(\gamma, \gamma' \in \mathbb R^q\) (e.g.\ full random-effect mean vector when
\(p_{\mathrm{re}} > 1\); see §13 Open Item 3 for rate-matrix / \(A\) geometry).
The merge draw uses a **fixed** multivariate Gaussian refresh density

\[
Q = N(c,\,\Sigma_Q),
\qquad
q_Q(\gamma') = \phi_q(\gamma';\, c,\, \Sigma_Q),
\]

with \(c \in \mathbb R^q\) and \(\Sigma_Q \succ 0\). The one-step kernel density
\(q(\gamma', \gamma)\) is the marginal in \(\gamma'\) after integrating the
block update (§10.1); on the **Normal template**,

\[
q(\gamma', \gamma) = \phi_q\bigl(\gamma';\, M(\gamma),\, \Sigma(\gamma)\bigr),
\qquad
M(\gamma) := E[\gamma' \mid \gamma],
\]

with \(\Sigma(\gamma) = \mathrm{Cov}(\gamma' \mid \gamma)\) (constant in \(\gamma\)
under the same conjugate Normal closure as §10.5 when applicable).

**Log gap and convexity.** Define \(g(\gamma', \gamma) = \log q(\gamma', \gamma)
- \log q_Q(\gamma')\). For **Gaussian** \(q(\cdot,\gamma)\),

\[
\nabla_{\gamma'} g(\gamma', \gamma)
=
-\,\Sigma(\gamma)^{-1}\bigl(\gamma' - M(\gamma)\bigr)
+
\Sigma_Q^{-1}\bigl(\gamma' - c\bigr),
\]

\[
\nabla_{\gamma'}^2 g
=
\Sigma_Q^{-1} - \Sigma(\gamma)^{-1}.
\]

**Precision form (Normal and non-Gaussian).** Write the kernel **precision**
(when defined)

\[
\Lambda_q(\gamma', \gamma)
\;:=\;
-\,\nabla_{\gamma'}^2 \log q(\gamma', \gamma).
\]

For Gaussian \(q\), \(\Lambda_q = \Sigma^{-1}\). The refresh uses
\(\Lambda_Q := \Sigma_Q^{-1}\). Then

\[
\nabla_{\gamma'}^2 g
=
\Lambda_Q - \Lambda_q.
\]

Require **Loewner dominance** (transition kernel **less sharp** than \(Q\) in
every direction):

\[
\boxed{
\Lambda_q(\gamma', \gamma)
\;\prec\;
\Lambda_Q
\;=\;
\Sigma_Q^{-1}
\qquad\text{(equivalently }\;
\Sigma(\gamma) \succ \Sigma_Q \text{ when } q \text{ is Normal).}
}
\]

Then \(\nabla_{\gamma'}^2 g \succ 0\) and \(g(\cdot,\gamma)\) is **strictly convex**
on \(\mathbb R^q\). (Scalar \(q=1\): \(\Lambda_q = 1/V\), \(\Lambda_Q = 1/\tau^2\),
\(V > \tau^2\).)

**Why this holds in the two-block sampler.** Set \(\Sigma_Q := \Sigma_{\mathrm{upd}}\)
from \(\gamma' \mid \beta\) (canonical choice below). The conditional
\(\kappa(\gamma',\beta)\) has Hessian \(-\Sigma_{\mathrm{upd}}^{-1}\) **independent
of \(\beta\)** and the likelihood. The **marginal** kernel \(q(\gamma',\gamma)\)
integrates \(\beta\): even when \(\pi(\beta\mid\gamma,y)\) is non-Gaussian (GLM),
\[
\Sigma
=
\mathrm{Cov}(\gamma'\mid\gamma)
=
\Sigma_{\mathrm{upd}}
+
\mathrm{Cov}\bigl(m(\beta)\mid\gamma,y\bigr)
\]
(§10.5 Normal; same identity for the second moment under any \(\pi_\gamma\)), so
the marginal is **more diffuse** and \(\Lambda_q \prec \Sigma_{\mathrm{upd}}^{-1}\)
at least in the **Normal closure** (exact \(\Lambda_q = \Sigma^{-1}\)). For
**GLMs**, verify \(\Lambda_q \prec \Lambda_Q\) on \(C_d \times \mathbb R^q\) from
\(\pi(\beta\mid\gamma,y)\) (numeric check of \(\nabla_{\gamma'}^2 g\), or a
curvature bound); the design target is that integrating \(\beta\) never **increases**
precision of the \(\gamma'\)-marginal above the block-2 value
\(\Sigma_{\mathrm{upd}}^{-1}\).

**Inner problem.** The unique minimizer satisfies \(\nabla_{\gamma'} g = 0\):

\[
\boxed{
\gamma'_*(\gamma)
=
\bigl(\Sigma_Q^{-1} - \Sigma^{-1}\bigr)^{-1}
\Bigl(\Sigma_Q^{-1} c - \Sigma^{-1} M(\gamma)\Bigr),
}
\]

(when \(\Sigma\) is constant in \(\gamma\), write \(\Sigma\) for
\(\Sigma(\gamma)\)). Equivalent forms use
\((\Sigma^{-1} - \Sigma_Q^{-1})^{-1}(\Sigma^{-1} M(\gamma) - \Sigma_Q^{-1} c)\)
with the matching sign convention from the FOC above.

**Inner value (Normal template).**

\[
\boxed{
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma_Q}{\det \Sigma}
\;-\;
\frac{1}{2}
\bigl(M(\gamma)-c\bigr)^\top
\bigl(\Sigma + \Sigma_Q\bigr)^{-1}
\bigl(M(\gamma)-c\bigr).
}
\]

(Check \(q=1\): recovers §10.5 \(\delta(\gamma)\). When \(M(\gamma)=c\),
\(\varepsilon(\gamma) = \sqrt{\det(\Sigma_Q \Sigma^{-1})}\).)

**Affine \(M(\gamma)\), \(\delta\)-level set \(C_d\), closed \(\varepsilon_d\)
(Normal conjugate).** With vector \(\gamma\), conjugate Normal Block~1 updates and
Block~2 GLS (§9.2 of `MINORIZATION_GAUSSIAN_REFRESH.md`; `notation.md`) give
\(M(\gamma) = \mathcal{L}\gamma + b(y)\) with \(\mathcal{L}\in\mathbb R^{q\times q}\).
Set \(c = \mu^\star := E[\gamma \mid y]\) (Hypothesis 2). Then

\[
M(\gamma) - \mu^\star
=
\mathcal{L}\bigl(\gamma - \mu^\star\bigr) + b_0(y),
\qquad
b_0(y) := \mathcal{L}\,\mu^\star + b(y) - \mu^\star,
\]

and \(\delta(\gamma)\) is the determinant term minus
\(\tfrac12 (\gamma-\mu_\delta)^\top\Gamma_\gamma(\gamma-\mu_\delta)\) (§5.1.1 of
the minorization note). Define the **coupling small set** from the **same**
quadratic (§5.1.2):

\[
\boxed{
C_d
=
\bigl\{
\gamma :
\bigl(\gamma - \mu_\delta\bigr)^\top \Gamma_\gamma \bigl(\gamma - \mu_\delta\bigr)
\le d
\bigr\},
\qquad
\Gamma_\gamma
=
\mathcal{L}^\top (\Sigma + \Sigma_Q)^{-1} \mathcal{L}.
}
\]

Then \(\min_{C_d}\delta(\gamma) = \tfrac12\log(\det\Sigma_Q/\det\Sigma) - d/2\) and

\[
\boxed{
\varepsilon_d
=
\sqrt{\det(\Sigma_Q \Sigma^{-1})}
\;
\exp\!\Bigl(-\frac{d}{2}\Bigr)
}
\]

(any \(\gamma_\star \in \partial C_d\)). See §6 of
`MINORIZATION_GAUSSIAN_REFRESH.md` (scalar intercept §9.1.2; general
\(\varepsilon_d = \exp(\delta(\mu^\star)-d/2)\) on the \(\delta\)-level set).

**Outer problem.** As in §10.3,

\[
\gamma_\star \in \arg\min_{\gamma \in C_d} \delta(\gamma),
\qquad
\boxed{
\varepsilon_d
=
\exp\bigl(\delta(\gamma_\star)\bigr)
=
\frac{
q\bigl(\gamma'_*(\gamma_\star),\, \gamma_\star\bigr)
}{
q_Q\bigl(\gamma'_*(\gamma_\star)\bigr)
}.
}
\]

Convexity implies \(g(\gamma', \gamma) \ge \delta(\gamma)\) for all
\(\gamma' \in \mathbb R^q\), certifying §10.1 minorization on
\(C_d \times \mathbb R^q\).

**Choosing \((c,\Sigma_Q)\).** *Canonical two-block Gibbs choice (recommended).*
Separate the **second-block conditional** from the **marginal one-step kernel**
(§10.1, §8.1). Given \(\beta\), the population draw is

\[
\gamma' \mid \beta \;\sim\; N\bigl(m(\beta),\,\Sigma_{\mathrm{upd}}\bigr),
\qquad
\kappa(\gamma',\beta) = \phi_q(\gamma';\, m(\beta),\, \Sigma_{\mathrm{upd}}),
\]

with \(\Sigma_{\mathrm{upd}}\) the **fixed** update covariance from the
hierarchy (scalar §8.1: \(\Sigma_{\mathrm{upd}} = \sigma_\gamma^2\);
multivariate: the population-prior / block-2 covariance). After integrating
\(\beta \sim \pi(\beta\mid\gamma,y)\),

\[
q(\gamma',\gamma) = \phi_q\bigl(\gamma';\, M(\gamma),\, \Sigma\bigr),
\qquad
\Sigma := \mathrm{Cov}(\gamma' \mid \gamma).
\]

**Law of total variance** gives the split used for refresh:

\[
\boxed{
\Sigma
\;=\;
\Sigma_{\mathrm{upd}}
\;+\;
\mathrm{Cov}\bigl(m(\beta) \mid \gamma, y\bigr),
}
\]

(scalar Normal §10.5: \(V = \sigma_\gamma^2 + w^2\,\mathrm{Var}(\bar\beta \mid
\gamma,y)\).) Take the refresh covariance **equal to the conditional** (same as
the Gibbs population step given \(\beta\)), not the marginal:

\[
\boxed{
\Sigma_Q \;:=\; \Sigma_{\mathrm{upd}}.
}
\]

**Choosing the mean \(c\).** With \(\Sigma_Q\) fixed, only the Mahalanobis term in
\(\delta(\gamma)\) depends on \(c\):

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma_Q}{\det \Sigma}
\;-\;
\frac{1}{2}
\bigl\|M(\gamma)-c\bigr\|_{(\Sigma+\Sigma_Q)^{-1}}^{2}.
\]

For **fixed** \(c\) (required by §10.1), \(\varepsilon_d = \exp(\min_{C_d}\delta)\)
is maximized by minimizing
\(\max_{\gamma \in C_d} \|M(\gamma)-c\|_{(\Sigma+\Sigma_Q)^{-1}}^2\) — a
**minimax center** of the set \(\{M(\gamma): \gamma \in C_d\}\), not by matching
\(c\) to each fiber’s kernel mean separately.

*Per-\(\gamma\) optimum (not usable for coupling).* At a **fixed** \(\gamma\),
the best center would be \(c = M(\gamma) = E[\gamma' \mid \gamma]\) (posterior /
predictive mean of the **next** population state given **current** \(\gamma\)),
which zeros the penalty at that \(\gamma\). Rosenthal’s merge step requires **one**
\(Q\) for all \(\gamma \in C_d\), so \(c\) cannot track \(M(\gamma)\).

*Full Gaussian case (recommended).* Define the **posterior mean** of the
population parameter (marginal on \(\gamma\), not the prior hypermean \(\mu_0\)):

\[
\boxed{
m_\gamma \;:=\; E[\gamma \mid y] \;=\; E_{\pi_\gamma}[\gamma]
\;=\;
\int \gamma\,\pi_\gamma(\gamma)\,d\gamma.
}
\]

Set

\[
\boxed{
c \;:=\; m_\gamma.
}
\]

In the conjugate Normal template, \(m_\gamma\) is the mean of the \(\gamma\)-marginal
and generally **differs from \(\mu_0\)** once data enter (only when
\(\pi_\gamma\) is centered at \(\mu_0\) do they coincide). The refresh draw
\(Q\) is centered where the **posterior** places population mass, not at the
prior mode.

*Why not \(c=\mu_0\)?* \(\mu_0\) enters the **prior** and the drift small set
\(C_d\) (§8–§10.1), but minorization is certified on \(C_d \times \mathbb R\)
with kernel means \(M(\gamma)\) tied to **updated** \(\beta \mid \gamma,y\).
Centering \(Q\) at \(m_\gamma\) aligns the merge with the stationary /
posterior location of the \(\gamma\)-chain.

*Why not \(c=M(\gamma)\)?* That is optimal **separately at each** \(\gamma\)
(zeros the Mahalanobis term at that fiber) but cannot be used as a fixed \(c\).
At \(\gamma = m_\gamma\), the kernel mean is \(M(m_\gamma)\), which need not
equal \(m_\gamma\) (Lemma Bias); still use \(c=m_\gamma\), not \(c=M(m_\gamma)\),
unless they are numerically equal.

*Refinement.* To maximize the **certified** \(\varepsilon_d = \exp(\min_{C_d}\delta)\)
over \(c\), solve the minimax problem
\(\min_c \max_{\gamma\in C_d} \|M(\gamma)-c\|_{(\Sigma+\Sigma_Q)^{-1}}^2\)
(convex program in \(c\)); \(c=m_\gamma\) is the default anchor and often near
the minimax when \(m_\gamma \in C_d\). **GLM:** compute \(m_\gamma\) from
\(\pi_\gamma\) (MC or quadrature); same rule for \(M(\gamma)\) in §10.4.

Then \(\Sigma - \Sigma_Q = \mathrm{Cov}(m(\beta)\mid\gamma,y)\succcurlyeq 0\), so
\(\Sigma \succ \Sigma_Q\) whenever the \(\beta\)-draw spreads the update mean
(strictly: \(w\neq 0\) and nontrivial \(\mathrm{Cov}(\bar\beta\mid\gamma,y)\) in
the scalar case). Hence \(\nabla_{\gamma'}^2 g = \Sigma^{-1}-\Sigma_Q^{-1}\succ
0\) (§10.7 convexity) **without** tuning a free \(\rho\).

At \(M(\gamma)=c\),

\[
\delta(\gamma)
=
\frac{1}{2}\log \frac{\det \Sigma_Q}{\det \Sigma}
=
-\frac{1}{2}\log \frac{\det \Sigma}{\det \Sigma_Q},
\qquad
\varepsilon(\gamma)
=
\sqrt{\det\bigl(\Sigma_Q \Sigma^{-1}\bigr)}.
\]

(scalar: \(\varepsilon(\gamma)=\sigma_\gamma/\sqrt{V}\).) The outer penalty term
\(\|M(\gamma)-c\|_{(\Sigma+\Sigma_Q)^{-1}}^2\) still depends on \(\gamma\);
\(\gamma_\star\) on \(C_d\) balances that against the determinant ratio.

*Alternative template.* \(\Sigma_Q = \rho\,\Sigma\) with \(\rho \in (0,1)\) if a
single scale knob is desired; at \(M(\gamma)=c\), \(\varepsilon(\gamma) =
\rho^{q/2}\). Refine \(\rho\) and \(c\) by maximizing \(\varepsilon_d\) after the
outer min (same two-stage recipe as §10.3).

**Small set \(C_d\).** For coupling with \(c = \mu^\star\) and affine \(M(\gamma)\),
use the **penalty ellipsoid** matched to \(\delta(\gamma)\):

\[
C_d
=
\bigl\{
\gamma :
(\gamma-\mu_\delta)^\top \Gamma_\gamma (\gamma-\mu_\delta) \le d
\bigr\},
\qquad
\Gamma_\gamma = \mathcal{L}^\top (\Sigma+\Sigma_Q)^{-1}\mathcal{L},
\]

with \(\mu_\delta\), \(\mathcal{L}\) from `MINORIZATION_GAUSSIAN_REFRESH.md` §5.1.1
and §9.1. Then \(\varepsilon_d = \exp(\delta(\mu^\star)-d/2)\) (§6 there).
For drift, intersect with the Foster–Lyapunov sublevel set from §8 if needed.
More generally, \(C_d = \{\gamma : V(\gamma) \le d\}\) for a vector Lyapunov
\(V\) once extended to \(\mathbb R^q\).

**GLM / non-Gaussian \(q\).** Keep \(\Sigma_Q = \Sigma_{\mathrm{upd}}\) and
\(Q = N(c,\Sigma_Q)\). The block-2 conditional Hessian is always
\(-\Sigma_{\mathrm{upd}}^{-1}\); certify
\(\Lambda_q(\gamma',\gamma) \prec \Sigma_{\mathrm{upd}}^{-1}\) (equivalently
\(\nabla_{\gamma'}^2 g \succ 0\)) on the relevant range. Then solve
\(\nabla_{\gamma'} g = 0\) for \(\gamma'_*(\gamma)\) and minimize
\(\delta(\gamma)\) over \(\gamma \in C_d\) (§10.3 table; §10.4 for \(\partial_{\gamma'}
v\)).

**Contrast (uniform ellipsoid refresh).** A fixed ellipsoid
\(\{\gamma' : (\gamma'-c)^\top \Sigma_Q^{-1}(\gamma'-c) \le r^2\}\) with
\(\mathrm{Uniform}\) \(Q\) makes \(\log q_Q\) **flat** on the support; \(g\) is
then concave in the \(\log\)-concave Normal direction — same endpoint issue as
the scalar uniform interval. **Gaussian** \(Q\) on all of \(\mathbb R^q\) is the
default for interior FOC-based \(\gamma'_*(\gamma)\).

---

## 11. Assembled TV bound

Write \(\pi_A(\cdot) := \pi(\cdot \mid \beta \in A)\) and
\(P_A^n(x,\cdot)\) for the target and \(n\)-step kernel restricted to \(A\)
(the on-\(A\) comparison term from `vignette("Chapter-C05")` / `ELLIPSOID_TV_BOUND.md`).

**Theorem TV (rate safe set).**

\[
\big\|P^n(x,\cdot) - \pi(\cdot)\big\|_{TV}
\;\le\;
\big\|P_A^n(x,\cdot) - \pi_A(\cdot)\big\|_{TV}
\;+\;
4\,P^n(x, A^c)
\;+\;
2\,\big|\pi(A^c) - P^n(x, A^c)\big|.
\]

Substitute §9 for \(P^n(x,A^c)\) and \(\pi(A^c)\), and §10 for the gap. For
**Normal / Gamma** with certified \(A^c = \varnothing\), the last two terms
vanish and TV reduces to the on-\(A\) term only.

**Per-likelihood summary.**

| Likelihood | \(P^n(x,A^c)\), \(\pi(A^c)\) | Gap | On-\(A\) term |
|---|---|---|---|
| Normal, Gamma | \(0\) | \(0\) | Heuristic / ellipsoid (§1, §6) |
| Logit | §9.1 + §2.2–§2.3; §9.3 | §10 | `Chapter-C05` |
| Probit | §9.1 + §3; §9.3 | §10 | same |
| Poisson, cloglog | §9.1 + §4–§5; §9.3 | §10 | same |

---

## 12. Summary table (conditional escape \(\pi(A^c \mid \gamma,y)\))

| Likelihood | \(W_j\) depends on \(\beta\)? | Closed \(\pi(A^c \mid \gamma,y)\) bound | Method |
|---|---|---|---|
| Normal | No | **Exact** \(0/1\) | §1 |
| Gamma (log) | No | **Exact** \(0/1\) | §6 |
| Logit, \(J=1\) | Yes | **Tail** \(G(\gamma)\); **drift** §9.1 | §2.2, §9 |
| Logit, \(J\ge1\) | Yes | **Sufficient** \(\sum_j G_j(\gamma)\); **drift** §9.1 | §2.3, §9 |
| Probit | Yes | **Sufficient** box + §9.1 | §3, §9 |
| Poisson | Yes | **Sufficient** floors \(s_j\) + §9.1 | §4, §9 |
| Cloglog | Yes | **Sufficient** floors \(s_j\) + §9.1 | §5, §9 |

---

## 13. Open items

1. **Multi-group logit:** Necessary closed form for \(\pi(A^c \mid \gamma,y)\)
   from \(\bar\omega \ge \tau/\lambda_b\) without a sufficient box.
2. **Probit / Poisson / cloglog:** Sharpen \(L_j\) constants; two-sided
   necessary regions for \(\bar\omega\).
3. **Non-scalar RE** (\(p_{\mathrm{re}} > 1\)): replace \(\omega_j\) by
   eigenvalues of the group coupling block; rate safe sets use ellipsoids (§12);
   **minorization / coupling:** §10.7 (Gaussian \(Q = N(c,\Sigma_Q)\) on
   \(\mathbb R^q\), two-stage \(\gamma'_*(\gamma)\) then \(\gamma_\star\));
   §8–§11 drift on vector \(\gamma\) still to be aligned with `Chapter-C05`.
4. **Weak \(\gamma\) prior:** bounds above are conditional on \(\gamma\);
   marginal \(\pi(A^c \mid y)\) requires outer integration (§9.2 with
   \(E_{\pi_\gamma}[\bar G(\gamma)]\)).
5. **Minorization constant \(\varepsilon_d\)** (§10.3–§10.7): Gaussian refresh
   \(Q = N(c,\tau^2)\) or \(N(c,\Sigma_Q)\), two-stage \(\gamma'_*(\gamma)\) then
   \(\gamma_\star\); Normal closed \(\delta(\gamma)\) (scalar §10.5, vector §10.7);
   GLM uses numeric \(\nabla_{\gamma'} g = 0\) (§10.4); optimize \((c,\Sigma_Q)\),
   \(d\), and \(r\).
6. **On-\(A\) TV term** (§11): still the `Chapter-C05` heuristic for GLMs;
   not likelihood-closed in this note.

---

## References

- `CHAPTER_C03_P_AND_A_MATRICES_BY_LIKELIHOOD.md` — rate matrix, \(\omega_j\),
  \(A/A^c\), §IV.5 certification strategies.
- `LOGIT_SINGLE_GROUP_SAFE_TIGHT_ARGUMENT.md` — template proofs for tilted
  tails (§2.2), drift (§8), coupling (§10).
- `LOGIT_SINGLE_GROUP_SAFE_REGION.md` — §4.6 minorization on
  \(C_d \times [\gamma'_-,\gamma'_+]\) (same notation as §10.3).
- `ELLIPSOID_TV_BOUND.md` — generic TV split (Lemma 1); on-\(A\) term in §11.
