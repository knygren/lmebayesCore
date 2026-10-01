# Foundational Lemmas and a General Theorem on Total Variation Distance Between Multivariate Normal Densities

**Companion to:** Nygren, "On the total variation distance between multivariate normal densities with applications to two-block Gibbs samplers"

These results build toward a single general theorem comparing total variation distances across pairs of multivariate normal densities, using only the sorted eigenvalues of each pair's precision-ratio matrix. Lemma 1 shows that any such comparison reduces, without loss of generality, to a diagonal reference matrix. Lemma 2 shows that, at this diagonal level, the true joint total variation distance for one pair is at least as large as for another pair whenever their (diagonal) eigenvalue structures agree except at a single coordinate, where the first pair's value is at least as large. Lemma 3 generalizes this, by repeated application of Lemma 2, to the case of full coordinatewise eigenvalue domination between two diagonal matrices. Theorem 1 then removes the diagonality requirement entirely, using Lemmas 1 and 3 together to show that plain sorted-eigenvalue domination is sufficient for arbitrary (non-commuting, non-simultaneously-diagonalizable) covariance matrices. A final section examines how far these ideas extend beyond the Gaussian case to independent coordinates with general symmetric log-concave marginals: a proven analog (Lemma 2′) requires an exactly shared reference density, and the natural relaxations — pointwise curvature-ratio domination, ratio domination combined with gap domination, and even that combination under the further restriction that pair 2 be exactly Gaussian — are each disproven with explicit counterexamples.

---

## Lemma 1 (Reduction to Diagonal Form; TV Distance Depends Only on the Eigenvalues)

**Lemma 1.** Let $\Sigma_1,\Sigma_2$ be any two positive definite $n\times n$ matrices, and let $M := \Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2} = C\Lambda C^T$ be an eigendecomposition, with $C$ orthogonal and $\Lambda = \mathrm{diag}(k_1,\dots,k_n)$. Then there exist positive definite matrices $\Sigma_3,\Sigma_4$ such that:

(a) $\Sigma_4^{1/2}\Sigma_3^{-1}\Sigma_4^{1/2}$ is diagonal, with eigenvalues exactly $k_1,\dots,k_n$ (the same eigenvalues as $M$, in the same order); and

(b) $\big\|N(0,\Sigma_1) - N(0,\Sigma_2)\big\|_{TV} = \big\|N(0,\Sigma_3) - N(0,\Sigma_4)\big\|_{TV}$.

In particular, if the eigenvalues of $M$ are all $\ge 1$ (the hypothesis of Lemma 2 of the main paper), the same holds for $\Sigma_4^{1/2}\Sigma_3^{-1}\Sigma_4^{1/2}$.

**Proof.** Take $\Sigma_3 := \Lambda^{-1}$ and $\Sigma_4 := I_n$.

*(a)* $\Sigma_4^{1/2}\Sigma_3^{-1}\Sigma_4^{1/2} = I^{1/2}\,\Lambda\,I^{1/2} = \Lambda$, already diagonal with eigenvalues $k_1,\dots,k_n$.

*(b)* Define the invertible linear map $T(x) := C^T\Sigma_2^{-1/2}x$. Total variation distance is invariant under any invertible (indeed, any bijective measurable) transformation of the sample space: since $T$ maps Borel sets to Borel sets bijectively,
$$\big\|T_\#P - T_\#Q\big\|_{TV} = \sup_A\big|P(T^{-1}A)-Q(T^{-1}A)\big| = \sup_B\big|P(B)-Q(B)\big| = \|P-Q\|_{TV}.$$

Push both mean-zero normal measures through $T$. Using orthogonality of $C$ ($C^TC=I$):
$$T\Sigma_2T^T = C^T\Sigma_2^{-1/2}\Sigma_2\Sigma_2^{-1/2}C = C^TC = I = \Sigma_4, \qquad\text{so } T_\#N(0,\Sigma_2)=N(0,\Sigma_4).$$
Since $\Sigma_2^{-1/2}\Sigma_1\Sigma_2^{-1/2} = M^{-1} = C\Lambda^{-1}C^T$:
$$T\Sigma_1T^T = C^T\Sigma_2^{-1/2}\Sigma_1\Sigma_2^{-1/2}C = C^T\big(C\Lambda^{-1}C^T\big)C = \Lambda^{-1} = \Sigma_3, \qquad\text{so } T_\#N(0,\Sigma_1)=N(0,\Sigma_3).$$
Combining,
$$\big\|N(0,\Sigma_1)-N(0,\Sigma_2)\big\|_{TV} = \big\|T_\#N(0,\Sigma_1)-T_\#N(0,\Sigma_2)\big\|_{TV} = \big\|N(0,\Sigma_3)-N(0,\Sigma_4)\big\|_{TV}. \qquad\blacksquare$$

**Remark.** Lemma 1 shows that the total variation distance between two mean-zero multivariate normal densities depends on $(\Sigma_1,\Sigma_2)$ *only* through the eigenvalues of $\Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2}$ — the eigenbasis $C$ and the overall scale of $\Sigma_2$ are TV-irrelevant, since both can be removed by an invertible linear change of variables without changing the distance. This is precisely the invariance principle used implicitly in the opening equalities of the proof of Lemma 2 in the main paper (the passage through $f(\cdot\mid \Sigma_2^{-1/2}\mu, [\Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2}]^{-1})$ and the rotation by $C^{-1}$); stating it here as a standalone lemma makes it directly citable wherever a reduction to diagonal form is needed.

---

## Lemma 2 (Comparison of TV Distances for Diagonal Pairs Differing in One Eigenvalue)

**Lemma 2.** Let $\Sigma_1,\Sigma_2,\Sigma_3,\Sigma_4$ be positive definite $n\times n$ matrices such that
$$M := \Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2} = \mathrm{diag}(k_1,\dots,k_n), \qquad M' := \Sigma_4^{1/2}\Sigma_3^{-1}\Sigma_4^{1/2} = \mathrm{diag}(k_1,\dots,k_{j-1},k_j',k_{j+1},\dots,k_n)$$
are both diagonal and agree in every coordinate except $j$, where $k_j' \ge k_j$. Then
$$\big\|N(0,\Sigma_1)-N(0,\Sigma_2)\big\|_{TV} \;\le\; \big\|N(0,\Sigma_3)-N(0,\Sigma_4)\big\|_{TV}.$$
The inequality is strict except in a measure-zero degenerate case (identified in the proof).

No sign or magnitude restriction is needed on the shared coordinates $k_i$, $i\ne j$.

**Proof.** By Lemma 1, $\|N(0,\Sigma_1)-N(0,\Sigma_2)\|_{TV} = \|N(0,M^{-1})-N(0,I)\|_{TV}$ and $\|N(0,\Sigma_3)-N(0,\Sigma_4)\|_{TV} = \|N(0,M'^{-1})-N(0,I)\|_{TV}$ (here the eigendecomposition step of Lemma 1 is trivial, since $M,M'$ are already diagonal). It therefore suffices to prove the inequality for diagonal $\Lambda=M$, $\Lambda'=M'$ against the common reference $I$.

For $\kappa\ge k_j$, let $\Lambda_\kappa$ denote the diagonal matrix agreeing with $\Lambda$ (equivalently $\Lambda'$) except that coordinate $j$ is set to $\kappa$, and define
$$\phi(\kappa) := \big\|N(0,\Lambda_\kappa^{-1}) - N(0,I)\big\|_{TV}, \qquad\text{so that } \phi(k_j) = \|N(0,\Sigma_1)-N(0,\Sigma_2)\|_{TV},\ \ \phi(k_j') = \|N(0,\Sigma_3)-N(0,\Sigma_4)\|_{TV}.$$

Write $\mu_\kappa := N(0,\Lambda_\kappa^{-1})$, $\nu:=N(0,I)$, and factor $\mu_\kappa(x) = h(x_{-j})\,\varphi_\kappa(x_j)$, where $\varphi_\kappa$ is the $N(0,1/\kappa)$ density and $h(x_{-j}) = \prod_{i\ne j} N(0,1/k_i)$-density, independent of $\kappa$.

The crossing region $A_\kappa = \{x:\mu_\kappa(x)\ge\nu(x)\}$ is $\{x: \sum_i(1-k_i)x_i^2 + \log k_i \ge 0\}$ (with $k_j$ replaced by $\kappa$ in the $j$-th term). Isolating coordinate $j$ and writing $S(x_{-j}) := \sum_{i\ne j}\big[(1-k_i)x_i^2+\log k_i\big]$ (independent of $\kappa$), the condition becomes, for $\kappa>1$,
$$x_j^2 \;\le\; \tau(x_{-j},\kappa) := \frac{S(x_{-j})+\log\kappa}{\kappa-1},$$
an interval $[-\sqrt\tau,\sqrt\tau]$ in $x_j$ (empty when $\tau<0$). Hence
$$\phi(\kappa) = \int h(x_{-j})\Big[2\Phi\big(\sqrt{\kappa\tau}\big)-1\Big]\mathbb 1(\tau\ge0)\,dx_{-j}.$$

Since $A_\kappa$ is the TV-optimal region (density equality holds exactly on its boundary), the envelope theorem permits differentiating in $\kappa$ while holding $A_\kappa$ fixed:
$$\phi'(\kappa) = \mathbb E_{\mu_\kappa}\!\left[\Big(\tfrac{1}{2\kappa}-\tfrac{X_j^2}{2}\Big)\mathbb 1_{A_\kappa}\right],$$
using $\partial\varphi_\kappa/\partial\kappa(x_j) = \varphi_\kappa(x_j)\big[\tfrac1{2\kappa}-\tfrac{x_j^2}2\big]$.

Conditional on $x_{-j}$ (with $\tau=\tau(x_{-j},\kappa)\ge0$), $X_j\sim N(0,1/\kappa)$, and the standard truncated-normal second-moment identity gives
$$\mathbb E\big[X_j^2\mathbb 1(X_j^2\le\tau)\big] = \tfrac1\kappa\big(2\Phi(\sqrt{\kappa\tau})-1\big) - \tfrac{2}{\sqrt\kappa}\sqrt\tau\,\varphi(\sqrt{\kappa\tau}),$$
while $\mathbb P(X_j^2\le\tau) = 2\Phi(\sqrt{\kappa\tau})-1$. Substituting, the $\big(2\Phi(\sqrt{\kappa\tau})-1\big)$ terms cancel *exactly*, leaving
$$\frac{1}{2\kappa}\mathbb P(A_\kappa\mid x_{-j}) - \frac12\,\mathbb E\big[X_j^2\mathbb 1_{A_\kappa}\mid x_{-j}\big] \;=\; \sqrt{\tau/\kappa}\,\varphi\big(\sqrt{\kappa\tau}\big).$$

Integrating over $x_{-j}$:
$$\phi'(\kappa) = \int h(x_{-j})\,\sqrt{\tau(x_{-j},\kappa)/\kappa}\;\varphi\big(\sqrt{\kappa\,\tau(x_{-j},\kappa)}\big)\,\mathbb 1\big(\tau(x_{-j},\kappa)\ge0\big)\,dx_{-j} \;\ge\; 0,$$
with equality only if $\tau(x_{-j},\kappa)=0$ for almost every $x_{-j}$ — the measure-zero degenerate case referred to in the statement. Since $\phi$ is differentiable on $[k_j,k_j']$ with $\phi'\ge0$ throughout, the fundamental theorem of calculus gives
$$\big\|N(0,\Sigma_3)-N(0,\Sigma_4)\big\|_{TV} - \big\|N(0,\Sigma_1)-N(0,\Sigma_2)\big\|_{TV} \;=\; \phi(k_j') - \phi(k_j) \;=\; \int_{k_j}^{k_j'}\phi'(\kappa)\,d\kappa \;\ge\; 0,$$
which is the claimed inequality, strict unless $\phi'\equiv0$ on $[k_j,k_j']$. $\blacksquare$

**Numerical verification.** The derivative formula above was checked against direct central finite-difference differentiation of $\phi$ (computed via exact 2D numerical integration of the joint density difference) at $k_2=9$ fixed and $k_1\in\{2,4,6.6,8,8.9\}$: agreement to $10^{-6}$ at every point tested. Direct computation of the two endpoint distances for $(k_j,k_j')=(6.6,8.9)$ with $k_2=9$ fixed gives $\phi(6.6)=0.642293 \le \phi(8.9)=0.674243$, confirming the inequality on this case. The inequality was further confirmed across a fine grid $k_1\in[1.5,9]$ (14 points, $k_2=9$ fixed): $\phi$ increasing throughout, with no reversals. A further check with an eigenvalue below 1 elsewhere ($k_2=0.5$, outside the paper's usual $\ge1$ hypothesis) also confirmed the inequality, consistent with the proof placing no sign restriction on $k_i$, $i\ne j$.

**Remark.** This is a statement about the *true* total variation distance, obtained by direct differentiation and integration of the exact joint density difference — it is not derived from, and does not always agree in direction with, the triangle-inequality upper bound of Lemma 2 in the main paper (that bound, built from a chain of intermediate matrices, can in fact decrease under exactly the kind of single-coordinate eigenvalue increase for which this Lemma 2 guarantees the true distance is ordered correctly). The two objects — the true distance and the main paper's upper bound on it — should not be assumed to move together under coordinatewise perturbations.

---

## Lemma 3 (Comparison of TV Distances Under Full Diagonal Eigenvalue Domination)

**Lemma 3.** Let $\Sigma_1,\Sigma_2,\Sigma_3,\Sigma_4$ be positive definite $n\times n$ matrices such that
$$M := \Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2} = \mathrm{diag}(k_1,\dots,k_n), \qquad M' := \Sigma_4^{1/2}\Sigma_3^{-1}\Sigma_4^{1/2} = \mathrm{diag}(k_1',\dots,k_n')$$
are both diagonal, with $k_i \le k_i'$ for every $i=1,\dots,n$ (no ordering among the $k_i$ themselves, and no sign restriction, is required). Then
$$\big\|N(0,\Sigma_1)-N(0,\Sigma_2)\big\|_{TV} \;\le\; \big\|N(0,\Sigma_3)-N(0,\Sigma_4)\big\|_{TV}.$$

**Proof.** Define a chain of diagonal matrices $\Lambda^{(0)},\Lambda^{(1)},\dots,\Lambda^{(n)}$ interpolating between $M$ and $M'$ coordinate by coordinate:
$$\Lambda^{(0)} := M = \mathrm{diag}(k_1,k_2,\dots,k_n), \qquad \Lambda^{(m)} := \mathrm{diag}(k_1',\dots,k_m',\,k_{m+1},\dots,k_n) \ \ (m=1,\dots,n),$$
so that $\Lambda^{(n)} = M'$, and $\Lambda^{(m)}$ differs from $\Lambda^{(m-1)}$ in exactly coordinate $m$, where $\Lambda^{(m)}$'s value ($k_m'$) is $\ge$ $\Lambda^{(m-1)}$'s value ($k_m$).

By Lemma 2 applied at each step (with $\Sigma_2=\Sigma_4=I$ and the two diagonal matrices being $\Lambda^{(m-1)},\Lambda^{(m)}$),
$$\big\|N(0,(\Lambda^{(m-1)})^{-1})-N(0,I)\big\|_{TV} \;\le\; \big\|N(0,(\Lambda^{(m)})^{-1})-N(0,I)\big\|_{TV}, \qquad m=1,\dots,n.$$
Chaining these $n$ inequalities by transitivity,
$$\big\|N(0,(\Lambda^{(0)})^{-1})-N(0,I)\big\|_{TV} \;\le\; \big\|N(0,(\Lambda^{(n)})^{-1})-N(0,I)\big\|_{TV},$$
i.e. $\|N(0,M^{-1})-N(0,I)\|_{TV} \le \|N(0,M'^{-1})-N(0,I)\|_{TV}$. By Lemma 1, $\|N(0,\Sigma_1)-N(0,\Sigma_2)\|_{TV} = \|N(0,M^{-1})-N(0,I)\|_{TV}$ and $\|N(0,\Sigma_3)-N(0,\Sigma_4)\|_{TV} = \|N(0,M'^{-1})-N(0,I)\|_{TV}$, giving the claim. $\blacksquare$

**Numerical verification.** Testing on the exact eigenvalue pair used to disprove the analogous claim for the main paper's Lemma 2 upper bound (see the companion supplement, §6): $k^f=(2,8)$, $k^g=(8,9)$, which satisfies $k_i^f\le k_i^g$ for both coordinates. Direct 2D numerical integration of the true joint densities gives
$$\|N(0,(\Lambda^f)^{-1})-N(0,I)\|_{TV} = 0.495983, \qquad \big(\text{intermediate, coordinate 1 raised}\big) = 0.650123, \qquad \|N(0,(\Lambda^g)^{-1})-N(0,I)\|_{TV} = 0.663003,$$
confirming both individual chain steps and the overall inequality $0.495983 \le 0.663003$.

**Remark (relation to the main paper's Lemma 2 and its disproof of the analogous bound comparison).** It is essential to distinguish Lemma 3 here from a superficially similar but false claim disproven elsewhere in the companion supplement (§6): there, it is shown that the main paper's Lemma 2 upper bound $\sum_i d_i$ — an artificial quantity built from a *different* telescoping chain (in which intermediate matrices are constant across a trailing *block* of coordinates, not diagonal with independently-set entries) — does **not** respect plain eigenvalue domination; an explicit example gives $k^f=(6.6,9)\le k^g=(8.9,9)$ with $\sum_i d_i^f > \sum_i d_i^g$. Lemma 3 does not contradict this. The two results concern different objects: Lemma 3 is a statement about the **true** total variation distance (verified above, on the very same eigenvalues that disprove the bound-level claim, to move in the *correct* direction), obtained by chaining the exact two-pair comparison of Lemma 2 here; the disproven claim in §6 concerns only the main paper's *upper bound* on that distance, which can be loose enough to reverse direction even while the true distance behaves exactly as Lemma 3 guarantees.

---

## Theorem 1 (General Comparison of TV Distances Under Sorted Eigenvalue Domination — No Diagonality or Commuting Assumption)

Lemma 3 required $M$ and $M'$ to already be diagonal — in effect, that both pairs be expressed in a shared, fixed coordinate system. This assumption can be removed entirely. The following theorem applies to **arbitrary** positive definite $\Sigma_1,\Sigma_2,\Sigma_3,\Sigma_4$, with no requirement that any of the four matrices commute or share an eigenbasis, either within a pair or between the two pairs.

**Theorem 1.** Let $\Sigma_1,\Sigma_2,\Sigma_3,\Sigma_4$ be any positive definite $n\times n$ matrices. Let $1 \le k_1\le\cdots\le k_n$ be the sorted eigenvalues of $M:=\Sigma_2^{1/2}\Sigma_1^{-1}\Sigma_2^{1/2}$, and $1\le k_1'\le\cdots\le k_n'$ the sorted eigenvalues of $M':=\Sigma_4^{1/2}\Sigma_3^{-1}\Sigma_4^{1/2}$. If
$$k_i \le k_i' \qquad\text{for every } i=1,\dots,n,$$
then
$$\big\|N(0,\Sigma_1)-N(0,\Sigma_2)\big\|_{TV} \;\le\; \big\|N(0,\Sigma_3)-N(0,\Sigma_4)\big\|_{TV}.$$

**Proof.** Apply Lemma 1 to the pair $(\Sigma_1,\Sigma_2)$, choosing the eigendecomposition $M=C\Lambda C^T$ with $\Lambda=\mathrm{diag}(k_1,\dots,k_n)$ in ascending sorted order. This gives $\Sigma_1^\dagger := \Lambda^{-1}$, $\Sigma_2^\dagger := I$ with
$$\big\|N(0,\Sigma_1)-N(0,\Sigma_2)\big\|_{TV} = \big\|N(0,\Sigma_1^\dagger)-N(0,\Sigma_2^\dagger)\big\|_{TV}.$$
Independently, apply Lemma 1 to the pair $(\Sigma_3,\Sigma_4)$, choosing the eigendecomposition $M'=C'\Lambda'C'^T$ with $\Lambda'=\mathrm{diag}(k_1',\dots,k_n')$, also in ascending sorted order. This gives $\Sigma_3^\dagger := \Lambda'^{-1}$, $\Sigma_4^\dagger := I$ with
$$\big\|N(0,\Sigma_3)-N(0,\Sigma_4)\big\|_{TV} = \big\|N(0,\Sigma_3^\dagger)-N(0,\Sigma_4^\dagger)\big\|_{TV}.$$

Crucially, the two applications of Lemma 1 are entirely independent: $C$ and $C'$ need not be related, and $\Sigma_1,\Sigma_2$ need not commute with $\Sigma_3,\Sigma_4$ in any sense. Nonetheless, both reductions land in the same diagonal-versus-identity form, and by construction $\Lambda,\Lambda'$ are diagonal matrices whose $i$-th entries are the two matrices' respective sorted eigenvalues, so that $k_i\le k_i'$ holds entrywise, position by position, by hypothesis.

Lemma 3 now applies directly to $\Sigma_1^\dagger,\Sigma_2^\dagger,\Sigma_3^\dagger,\Sigma_4^\dagger$ (which are diagonal, i.e. satisfy Lemma 3's hypothesis exactly), giving
$$\big\|N(0,\Sigma_1^\dagger)-N(0,\Sigma_2^\dagger)\big\|_{TV} \;\le\; \big\|N(0,\Sigma_3^\dagger)-N(0,\Sigma_4^\dagger)\big\|_{TV}.$$
Combining the two Lemma 1 equalities with this inequality gives the claim. $\blacksquare$

**Why no shared structure is needed.** The proof works because Lemma 1 diagonalizes each pair *on its own terms*, independently discarding that pair's own eigenbasis and overall scale (§ Remark following Lemma 1). The comparison in Lemma 3 is then applied not to the original matrices, nor to any shared rotation of them, but purely to the two resulting lists of sorted numbers $k_1,\dots,k_n$ and $k_1',\dots,k_n'$. Since Lemma 3's proof (chaining Lemma 2) only ever manipulates diagonal matrices entry-by-entry, it is entirely indifferent to what coordinate system produced those entries — so the two independent reductions can be combined even though $C\ne C'$ in general.

**Numerical verification.** Tested with $n=2$ and $\Sigma_1,\Sigma_2$ built from one random rotation and $\Sigma_3,\Sigma_4$ built from an unrelated second rotation (so that the two pairs share no eigenbasis whatsoever), with sorted eigenvalues $(k_1,k_2)=(3,6)$ for the first pair and $(k_1',k_2')=(4,7)$ for the second (satisfying $k_i\le k_i'$). Direct 2D numerical integration of the true joint densities gives
$$\|N(0,\Sigma_1)-N(0,\Sigma_2)\|_{TV} = 0.478956 \;\le\; \|N(0,\Sigma_3)-N(0,\Sigma_4)\|_{TV} = 0.551542,$$
confirming the theorem in a case where Lemma 3 alone (which requires literal diagonality) does not apply directly, and only the reduction via Lemma 1 makes the comparison possible.

**Remark.** Theorem 1 is the fully general resolution, at the level of the *true* total variation distance, of the question first posed for the main paper's Lemma 2 upper bound and disproven there in that form (companion supplement, §6): plain sorted eigenvalue domination — the simplest, most easily checked condition on two pairs of covariance matrices — is **sufficient** to order the true distances, unconditionally, with no requirement of commuting matrices, a shared eigenbasis, or the stronger ratio-domination condition used elsewhere (companion supplement, §3) to bound the main paper's telescoped sum. The main paper's Lemma 2 bound remains useful as a computable ceiling (companion supplement, §3–5), but Theorem 1 shows that the *true* distances themselves obey the simpler, more intuitive ordering that one would naturally have hoped for from the outset.

---

## Theorem 2 (Eigenvalue Monotonicity Along the Floor-to-True Homotopy)

This applies Theorem 1 along a specific path of matrices, arising when the true and floor likelihoods are both Gaussian and interpolated by a homotopy parameter $t\in[0,1]$ — the setting used throughout the restricted-sampler comparison (Applied Context, below).

**Setup.** $\gamma\in\mathbb R^d$. $L,\kappa_{\mathrm{true}}\in\mathbb R^{d\times d}$ symmetric positive definite, with $L\preceq\kappa_{\mathrm{true}}$ in the Loewner order. Define $\kappa_t:=(1-t)L+t\,\kappa_{\mathrm{true}}$, $V_t:=(\Psi^{-1}+\kappa_t)^{-1}$, and $M_t:=P_{12}V_tP_{21}P_{11}^{-1}$, where $P_{11}:=\Lambda_\gamma+H^\top\Psi^{-1}H$ and $P_{12}:=H^\top\Psi^{-1}=P_{21}^\top$ (the precision of $\gamma\mid\beta$; see Applied Context below for the full derivation of these quantities from the joint negative log-density). By Lemma 1, $M_t$'s eigenvalues are real and positive, via its similarity to the symmetric matrix $\tilde M_t:=P_{11}^{-1/2}P_{12}V_tP_{21}P_{11}^{-1/2}$.

**Theorem 2.** *Every eigenvalue of $M_t$ is monotonically decreasing in $t$, for every dimension $d$, with no hypothesis beyond $L\preceq\kappa_{\mathrm{true}}$.*

**Proof.**
1. $\kappa_t'=\kappa_{\mathrm{true}}-L\succeq0$, so $\kappa_t$ is Loewner-increasing in $t$ — immediate from linearity.
2. Matrix inversion reverses the Loewner order on positive definite matrices ($A\preceq B\Rightarrow B^{-1}\preceq A^{-1}$, standard), so $V_t=(\Psi^{-1}+\kappa_t)^{-1}$ is Loewner-**decreasing** in $t$.
3. Congruence preserves the Loewner order: for any matrix $C$, $A\succeq B\Rightarrow CAC^\top\succeq CBC^\top$, immediate from comparing quadratic forms $x^\top CAC^\top x=(C^\top x)^\top A(C^\top x)$. Applying this with $C=P_{11}^{-1/2}P_{12}$ shows $\tilde M_t$ is Loewner-decreasing in $t$.
4. The Courant–Fischer min-max characterization of the eigenvalues of a symmetric matrix is monotone in the Loewner order, coordinate by coordinate — so every eigenvalue of $\tilde M_t$, and hence of $M_t$ (same eigenvalues, by similarity), decreases individually as $t$ increases. $\blacksquare$

**Consequence via Theorem 1 — the shared mean, possibly moving with $t$, is essential to state explicitly.** Define $\gamma^\star(t)$ as the mean of $\pi_t$ **under the unrestricted ($r=\infty$) model** — this is the only place $\gamma^\star(t)$ is defined, via the fully-Gaussian case where it is given directly by the standard conjugate mean formula. $\pi_t$ and $Q_t^{(l)}$ share this mean at each $t$; the sampler is understood to be started at $\gamma^\star(t)$. When the model is exchange-symmetric, $\gamma^\star(t)\equiv0$ for every $t$ (No-Drift, established elsewhere); in general, $\gamma^\star(t)$ moves with $t$, and comparing $\pi_t$ to $Q_t^{(l)}$ correctly requires tracking it. This causes no difficulty here: **TV distance is translation-invariant** — $\|N(\mu,\Sigma_1)-N(\mu,\Sigma_2)\|_{TV}=\|N(0,\Sigma_1)-N(0,\Sigma_2)\|_{TV}$ for any $\mu$ — so Theorem 1 applies unchanged to $\pi_t=N(\gamma^\star(t),\Sigma_t)$ and $Q_t^{(l)}=N(\gamma^\star(t),\Sigma_t^{(l)})$ at each $t$, *provided only that the two objects share the same mean at that $t$*, not that the mean is constant across $t$. With this understood:
$$\phi(t):=\|\pi_t-Q_t^{(l)}\|_{TV} \text{ is strictly decreasing in } t \text{ on } [0,1],$$
for arbitrary dimension $d$ and arbitrary $l$, in a single argument — no chain rule, no covariance bookkeeping, and no assumption that $\gamma^\star(t)$ is fixed.

**Why this matters, not just as a technicality.** Silently writing $\pi_t=N(0,\Sigma_t)$ throughout — assuming the mean stays at $0$ for every $t$ — is exactly the assumption that fails without exchange symmetry, and doing so is what produced a spurious sign reversal in an unrelated non-Gaussian example examined earlier: the sampler stayed anchored at a fixed point while the true target's mean drifted away from it, and the resulting discrepancy was a location mismatch, not a genuine spread comparison. The theorem above is stated to make explicit that the *only* thing required is a shared, correctly-tracked mean at each $t$ — never a fixed one.

### Corollary (Sufficiently Large Ball Suffices)

**Setup.** Suppose $\beta$ is additionally restricted to a ball of radius $r$ (as in the restricted-sampler construction below), so that $\phi(t;r)$ denotes the resulting, generally non-Gaussian, TV distance at ball radius $r$, with $\phi(t;\infty)$ the unrestricted (Theorem 2) value. **The sampler at every finite $r$ is started at the same $\gamma^\star(t)$ defined above — the $r=\infty$ limiting mean, fixed as a function of $t$ once and for all — not at any $r$-dependent "true" mean of the ball-restricted model.** This is a deliberate choice, not merely a simplification: the actual mean of the ball-restricted $\pi_t$ generally differs from $\gamma^\star(t)$ at finite $r$ (truncation shifts it, by the same position-dependent mechanism noted elsewhere), and tracking that shifted quantity would reintroduce exactly the complication this corollary is built to avoid. Anchoring at the $r=\infty$ value throughout — a single, well-defined, $r$-independent function of $t$ — is what makes the convergence in step 1 below unambiguous: as $r\to\infty$, both the restricted $\pi_t$ and the restricted $Q_t^{(l)}$ (both still started at this same fixed $\gamma^\star(t)$) converge to their unrestricted counterparts, with no ambiguity about which reference point is being tracked.

**Corollary.** *There exists $\tilde r^\star<\infty$ such that for every $r>\tilde r^\star$, $\phi'(t;r)<0$ for every $t\in[0,1]$ simultaneously.*

**Proof.**
1. As $r\to\infty$, the ball-restricted model (sampler started at $\gamma^\star(t)$) converges to the unrestricted Gaussian model, so $\phi(t;r)\to\phi(t;\infty)$.
2. By Theorem 2's consequence above, $\phi'(t;\infty)<0$ for every $t\in[0,1]$ — strictly, unconditionally, in full generality.
3. Define $M(r):=\max_{t\in[0,1]}\phi'(t;r)$. Assuming $\phi'(t;r)$ is continuous in $r$, $M$ is continuous in $r$ (maximum of a continuous function over the compact set $[0,1]$, standard). Then $M(\infty)=\max_{t\in[0,1]}\phi'(t;\infty)<0$, since it is the maximum, over a compact set, of a function that is strictly negative at every point of that set.
4. By continuity of $M$ at $r=\infty$, there exists $\tilde r^\star$ such that $M(r)<0$ for every $r>\tilde r^\star$ — i.e. $\phi'(t;r)<0$ for every $t\in[0,1]$ simultaneously. $\blacksquare$

**Why this is the right proof strategy here.** The corollary does not require Theorem 2's Loewner-order argument to survive ball restriction — it only needs step 2's endpoint fact, which holds exactly and unconditionally by Theorem 2 itself. The compactness step then transports that single endpoint fact to a neighborhood of finite $r$, without any need to re-derive eigenvalue monotonicity under truncation (where $V$ becomes $\gamma$-dependent and the clean Loewner-order chain no longer directly applies), and without needing to characterize how the ball-restricted model's own mean shifts at finite $r$ — anchoring at $\gamma^\star(t)$ throughout sidesteps that question entirely. This mirrors the identical strategy used for the non-Gaussian case (companion roadmap document, restricted-sampler comparison), where the same compactness argument transports the $L=0$ endpoint fact to a finite ball threshold — here, the endpoint fact is a fully proven theorem rather than a numerically-supported one, so the corollary is unconditional wherever Theorem 2 applies.

### Closed form for $\phi'(t;\infty)$, multivariate, accounting for a moving mean

This makes explicit, in closed form, what Theorem 2's eigenvalue argument proves implicitly. To support generalization beyond the Gaussian case, every object below is defined first in terms of the fixed tilt $\log g(\beta):=\ell_{\mathrm{true}}(\beta)-\ell_{\mathrm{floor}}(\beta)$ and expectations of it — exactly the objects used in the companion roadmap's non-Gaussian treatment — with the Gaussian-specific algebra (Theorem 2's $\delta,M_t,V_t$) entering only afterward, as the closed-form *evaluation* of these general definitions, not as their starting point.

**General definitions (model- and dimension-independent).** Let
$$
\bar g_t(\gamma) := \mathbb E_{\pi_t(\beta\mid\gamma,y)}\big[\log g(\beta)\big]
$$
be the conditional expectation of the tilt given $\gamma$. The general tilting-derivative identity $\partial_t\pi_t(\gamma)=\pi_t(\gamma)\big[\bar g_t(\gamma)-\bar{\bar g}_t\big]$ (with $\bar{\bar g}_t:=\mathbb E_{\pi_t}[\bar g_t(\Gamma)]$) gives, on integrating over $A_t^c$:
$$
J_\pi(t) := \mathrm{Cov}_{\pi_t}\big(\bar g_t(\Gamma),\ \mathbf 1_{A_t^c}(\Gamma)\big).
$$
For the sampler side, reusing the $l$-term construction from the companion roadmap's non-Gaussian treatment directly — $\tilde h_m^c(\beta)$ the probability the remaining $l-m$ Gibbs steps land in $A_t^c$, starting from $\beta$ at step $m$ — gives
$$
J_Q(t) := \sum_{m=1}^{l} \mathbb E_{\gamma_{m-1}}\Big[\mathrm{Cov}_{\pi_t(\beta\mid\gamma_{m-1},y)}\big(\tilde h_m^c(\beta),\ \log g(\beta)\big)\Big].
$$
Both $\bar g_t$ and $\tilde h_m^c$ are defined for arbitrary $\gamma$-dimension $d$ and arbitrary log-concave $\ell_{\mathrm{true}}$; nothing here assumes Gaussianity yet. With this, $\phi'(t;\infty)=J_Q(t)-J_\pi(t)$ up to the mean-movement contribution addressed next — the general statement of this program's central identity.

**Verified**: for $d=1$ (the running scalar example), $J_\pi(t)$ computed directly from this general definition matches the value obtained via the Gaussian-specific quadratic-form route below to 6+ significant figures ($0.054401$ both ways) — confirming the two routes are the same object, general definition versus closed-form evaluation.

**Where the fixed tilt enters explicitly once we specialize to Gaussian.** When $\ell_{\mathrm{true}},\ell_{\mathrm{floor}}$ are both Gaussian, $\log g(\beta)=-\tfrac12(\beta-\beta_0)^\top\delta(\beta-\beta_0)$ with $\delta:=\kappa_{\mathrm{true}}-L=-(\log g)''(\beta)$, constant. Evaluating $\bar g_t(\gamma)$ under this quadratic tilt (via the standard Gaussian quadratic-form expectation identity) and carrying the resulting derivatives through gives:
$$
V_t' = -V_t\,\delta\,V_t \qquad\text{(verified: matches numerical differentiation of }V_t\text{ to }6\times10^{-12}\text{)},
$$
$$
M_t' = -P_{12}V_t\,\delta\,V_tP_{21}P_{11}^{-1},
$$
so that every derivative used from here on is built from $\delta$ alone — the curvature of $\log g$, and the only place the true-versus-floor comparison enters this evaluation.

**Setup.** $\pi_t=N(\gamma^\star(t),\Sigma_t)$, $Q_t^{(l)}=N(\gamma^\star(t),\Sigma_t^{(l)})$. The standard multivariate Gaussian score identity gives, for either density (writing $\Sigma$ generically for $\Sigma_t$ or $\Sigma_t^{(l)}$):
$$
\partial_t\log N(\gamma;\gamma^\star(t),\Sigma) =
\underbrace{(\gamma-\gamma^\star(t))^\top\Sigma^{-1}(\gamma^\star)'(t)}_{\text{mean-movement}}
$$
$$
+\ \underbrace{\tfrac12(\gamma-\gamma^\star(t))^\top\Sigma^{-1}\Sigma'\Sigma^{-1}(\gamma-\gamma^\star(t))}_{\text{spread-change}}
$$
$$
-\ \tfrac12\mathrm{tr}(\Sigma^{-1}\Sigma').
$$
This is exactly $\bar g_t(\gamma)-\bar{\bar g}_t$ evaluated for the Gaussian tilt — the spread-change piece is $\bar g_t(\gamma)$ itself (up to the additive constant that a covariance discards), and the mean-movement piece is the new term arising because $\gamma^\star(t)$ may move. Since $\mathbb E_{\pi_t}[\partial_t\log\pi_t]=0$ identically and constants drop out of a covariance, integrating $\partial_t\pi_t=\pi_t\cdot\partial_t\log\pi_t$ over $A_t^c$ gives:
$$
\phi'(t;\infty) = \big[J_Q(t)-J_\pi(t)\big] + \big[N_Q(t)-N_\pi(t)\big],
$$
$$
J_\pi(t)=\mathrm{Cov}_{\pi_t}\Big(\tfrac12(\Gamma-\gamma^\star)^\top\Sigma_t^{-1}\Sigma_t'\Sigma_t^{-1}(\Gamma-\gamma^\star),\ \mathbf 1_{A_t^c}(\Gamma)\Big),
$$
$$
N_\pi(t):=\mathrm{Cov}_{\pi_t}\big((\Gamma-\gamma^\star)^\top\Sigma_t^{-1}(\gamma^\star)',\ \mathbf 1_{A_t^c}(\Gamma)\big),
$$
with $J_Q,N_Q$ defined identically under $Q_t^{(l)}$ using $\Sigma_t^{(l)},(\Sigma_t^{(l)})'$ in place of $\Sigma_t,\Sigma_t'$ — this $J_\pi(t)$ is the Gaussian closed-form evaluation of the general $J_\pi(t)$ defined above from $\bar g_t$, confirmed to agree numerically.

**$J_\pi(t)$ written using $\log g(\beta)$ directly — comparable in structure to $J_Q(t)$.** Since $\mathbf 1_{A_t^c}(\Gamma)$ depends only on $\gamma$, the law of total covariance gives, for the full joint $\pi_t(\gamma,\beta\mid y)$:
$$
\mathrm{Cov}_{\pi_t}\big(\bar g_t(\Gamma),\mathbf 1_{A_t^c}(\Gamma)\big) = \mathrm{Cov}_{\pi_t(\gamma,\beta)}\big(\log g(\beta),\ \mathbf 1_{A_t^c}(\Gamma)\big),
$$
exactly, with no approximation (verified directly: the joint-form covariance, computed by direct integration, matches $J_\pi(t_0)=0.05441$ to the precision of the numerical crossing point used). This is deliberately the same structure as $J_Q(t)$'s own definition — a covariance against $\log g(\beta)$ directly, under a joint density — with the only difference being *which* joint density and *which* paired weight: $\pi_t(\gamma,\beta)$ paired with $\mathbf 1_{A_t^c}(\Gamma)$ here, versus $\pi_t(\beta\mid\gamma_{m-1},y)$ paired with $\tilde h_m^c(\beta)$ there.

**The mean-movement terms vanish exactly — a proof, not an assumption.** Because $\pi_t$ and $Q_t^{(l)}$ share the same mean $\gamma^\star(t)$, their crossing set $A_t^c=\{\gamma:Q_t^{(l)}(\gamma)>\pi_t(\gamma)\}$ depends on $\gamma$ only through the quadratic form $(\gamma-\gamma^\star(t))^\top\big[(\Sigma_t^{(l)})^{-1}-\Sigma_t^{-1}\big](\gamma-\gamma^\star(t))$ — an **ellipsoid centered exactly at $\gamma^\star(t)$**, symmetric under the reflection $(\gamma-\gamma^\star(t))\mapsto-(\gamma-\gamma^\star(t))$. Both $\pi_t$ and $Q_t^{(l)}$ are themselves symmetric under this same reflection (Gaussians centered at $\gamma^\star(t)$). The mean-movement weight $(\gamma-\gamma^\star)^\top\Sigma^{-1}(\gamma^\star)'$ is linear, hence **odd** under this reflection, while $\mathbf 1_{A_t^c}$ is **even** — a covariance between an odd and an even function, under a measure symmetric about the same center, is exactly zero:
$$
N_\pi(t) = N_Q(t) = 0, \qquad \text{regardless of the value of } (\gamma^\star)'(t).
$$
This is the precise, rigorous reason a moving mean does not affect $\phi'(t;\infty)$ — not because the mean happens to be fixed, but because the geometry of comparing two same-mean Gaussians is symmetric about that shared mean no matter how it moves. (Verified directly: for an artificial, nontrivial $\gamma^\star(t)$ with $(\gamma^\star)'(t_0)\ne0$, both $N_\pi(t_0)$ and $N_Q(t_0)$ computed by direct integration vanish to numerical precision.)

**The resulting closed form.**
$$
\boxed{\phi'(t;\infty) = J_Q(t) - J_\pi(t), \qquad J_\pi(t)=\mathrm{Cov}_{\pi_t(\gamma,\beta)}\big(\log g(\beta),\ \mathbf 1_{A_t^c}(\Gamma)\big),}
$$
with $J_Q(t)$ defined analogously, given in full just below.

**Restated in full, as the general form this entire derivation evaluates:**
$$
\boxed{J_Q(t) := \sum_{m=1}^{l} \mathbb E_{\gamma_{m-1}}\Big[\mathrm{Cov}_{\pi_t(\beta\mid\gamma_{m-1},y)}\big(\tilde h_m^c(\beta),\ \log g(\beta)\big)\Big].}
$$

**Each object inside this, explicit for the multivariate Gaussian case** (verified: the conditional mean below matches the mode of the true joint negative log-density directly, to numerical precision):
$$
\boxed{
\begin{aligned}
\gamma_{m-1} &\sim N\big(\gamma^\star(t),\ \Sigma_t^{(m-1)}\big), \\[4pt]
\pi_t(\beta\mid\gamma_{m-1},y) &= N\Big(\beta_0(t) + V_tP_{12}^\top(\gamma_{m-1}-\gamma^\star(t)),\ V_t\Big), \qquad \beta_0(t):=H\gamma^\star(t), \\[4pt]
\log g(\beta) &= -\tfrac12(\beta-\beta_0(t))^\top\delta\,(\beta-\beta_0(t)), \\[4pt]
\gamma_l \mid \beta &\sim N\Big(\gamma^\star(t) + M_t^{\,l-m}P_{12}^\top P_{11}^{-1}(\beta-\beta_0(t)),\ \ \Xi_m\Big), \\[4pt]
\Xi_m &:= M_t^{\,l-m}P_{11}^{-1}\big(M_t^{\,l-m}\big)^\top + \Sigma_t^{(l-m)}, \\[4pt]
\tilde h_m^c(\beta) &= \mathbb P\big(\gamma_l\in A_t^c \mid \beta\big) = \int_{A_t^c} N\big(\gamma;\ \gamma^\star(t)+M_t^{\,l-m}P_{12}^\top P_{11}^{-1}(\beta-\beta_0(t)),\ \Xi_m\big)\,d\gamma.
\end{aligned}
}
$$
Every object above is fully explicit **except** $\tilde h_m^c(\beta)$ itself: for $d=1$ it reduces to a difference of two normal CDFs (used throughout the companion roadmap's scalar treatment), but for $d>1$, $A_t^c$ is a genuine ellipsoid (not an interval), so $\tilde h_m^c(\beta)$ is a non-central, ellipsoid-restricted Gaussian probability — expressible via a non-central chi-square CDF in the basis that diagonalizes $\Xi_m$ against the ellipsoid's own defining quadratic form, but not through an elementary closed form the way the scalar CDF-difference is. This is stated as an integral rather than forced into a false closed form.

**Correction flag:** the $M_t$ used in the $\gamma_l\mid\beta$ mean and $\Xi_m$ above has a transpose-ordering bug identified later in this document (§ "Multivariate extension," corrected using Nygren's own $\hat M=P_{11}^{-1}P_{12}P_{22}^{-1}P_{21}$) — this box should be re-derived with $\hat M$ before being relied upon.

Everything derived from here on — $M_t,V_t,C_m(t),W_m(\gamma)$ — is the Gaussian closed-form content of $\tilde h_m^c(\beta)$ and $\log g(\beta)$ evaluated explicitly; the definition itself is model-independent and requires no Gaussianity, which is precisely what lets it serve, unchanged, as the general target for the non-Gaussian case in the companion roadmap.

### Fully evaluated, scalar case ($d=1$): closed forms for $J_\pi(t)$ and each $J_Q^{(m)}(t)$

**$J_\pi(t)$**, via the standard truncated-normal second-moment identity applied to $\pi_t=N(\gamma^\star,\Sigma_t)$:
$$
\boxed{J_\pi(t) = \delta\,s_t^2\,\Sigma_t\,z_0\,\varphi(z_0), \qquad z_0 := \frac{c}{\sqrt{\Sigma_t}}, \qquad s_t = V_tP_{12}.}
$$

**Each term $J_Q^{(m)}(t)$**, derived from the second-order Stein identity $\mathrm{Cov}(f(x),x^2)=V^2\mathbb E[f'']+2\mu V\mathbb E[f']$ applied twice — once for the inner ($\beta$) layer, once for the outer ($\gamma_{m-1}$) layer (verified against direct nested numerical integration, ground truth matching to 10 decimal places after correcting a sign error caught in that verification):
$$
\boxed{J_Q^{(m)}(t) = \frac{\delta B_m^2\,c}{\bar S_m^3}\,\varphi\!\left(\frac{c}{\bar S_m}\right)\Big[V_t^2 + 2V_t\,\sigma_{m-1}^2\Big],}
$$
$$
B_m := a_t^{\,l-m}\frac{P_{12}}{P_{11}}, \qquad S_m^2 := \frac{a_t^{2(l-m)}}{P_{11}}+\Sigma_t^{(l-m)}, \qquad \tilde S_m^2 := B_m^2V_t+S_m^2,
$$
$$
\sigma_{m-1}^2 := s_t^2\,\Sigma_t^{(m-1)}, \qquad \bar S_m^2 := \tilde S_m^2 + B_m^2\sigma_{m-1}^2.
$$
At $m=1$ ($\sigma_0^2=0$, since $\gamma_0=\gamma^\star$ exactly, a point mass): this reduces to $J_Q^{(1)}(t)=\delta(P_{12}/P_{11})^2\,cV_t^2/(\Sigma_t^{(1)})^{3/2}\cdot\varphi(z_1)$, matching the $l=1$-specific formula derived earlier exactly.

**Implications for the convergence rate.** $B_m^2=a_t^{2(l-m)}(P_{12}/P_{11})^2$ governs each summand's size and decays geometrically in $l-m$ — the distance back from the final step — so the sum's shape does not spread out as $l$ grows: terms near $m=l$ stay $O(1)$, while terms far from $m=l$ are suppressed by $a_t^{2(l-m)}$ regardless of how large $l$ is. Consequently:
- the number of terms that matter is bounded independent of $l$ — only $O(1/\log(1/a_t^2))$ terms near the end carry essentially the whole sum, for any sufficiently large $l$;
- the gap $J_\pi(t)-\sum_{m=1}^lJ_Q^{(m)}(t)$ inherits this same $a_t^2$ geometric rate — it is exactly the missing contribution from steps beyond $l$, each smaller than the last by a factor of $a_t^2$ through the same $B_m^2$ mechanism;
- this matches the rate already established for $\Sigma_t-\Sigma_t^{(l)}=\Sigma_t\,a_t^{2l}$ exactly (Remark 13), since $J_Q^{(l)}(t)\to J_\pi(t)$ is ultimately driven by the same underlying convergence $\Sigma_t^{(l)}\to\Sigma_t$.

### Multivariate extension — general, exact, via Nygren's own $A$ and $C$ (corrected)

**An earlier pass at this section was wrong, due to a genuine bug in the covariance recursion — corrected here after re-checking directly against the source paper (Nygren, Claim 3/Remark 13).** I had used $M_t\Sigma^{(m-1)}M_t^\top$ with $M_t:=P_{12}P_{22}^{-1}P_{21}P_{11}^{-1}$ ($P_{11}^{-1}$ on the right), plus an inhomogeneous term with a spurious extra factor. The correct recursion (Nygren's own) uses
$$
\hat M := P_{11}^{-1}P_{12}P_{22}^{-1}P_{21} \qquad (P_{11}^{-1}\text{ on the left}; \ \hat M = M_t^\top),
$$
$$
\Sigma_t^{(l)} = \big[P_{11}^{-1}+\hat MP_{11}^{-1}\big] + \hat M\,\Sigma_t^{(l-1)}\,\hat M^\top, \qquad \Sigma_t^{(0)}:=0,
$$
verified against Nygren's closed form $\Sigma_t^{(l)}=P_{11}^{-1/2}\big[\sum_{i=1}^{2l}A^{i-1}\big]P_{11}^{-1/2}$ exactly ($d=3$, $l=1,2,3$, matching to machine precision). With this corrected recursion, **there is no special case and no separate diagonalization needed**:

$$
\boxed{\Sigma_t \text{ and } \Sigma_t^{(l)} \text{ are simultaneously diagonalized, for every } l, \text{ by the single transform } u:=C^\top P_{11}^{1/2}\gamma,}
$$

where $C$ diagonalizes $A=\tilde M_t=P_{11}^{-1/2}P_{12}P_{22}^{-1}P_{21}P_{11}^{-1/2}$ itself — the **same** $C$ used throughout Theorem 2, with no dependence on $l$. Explicitly (verified numerically to match Nygren's closed form exactly, $d=3$):
$$
(\Sigma_t)_{u,kk} = \frac{1}{1-a_k}, \qquad (\Sigma_t^{(l)})_{u,kk} = \frac{1-a_k^{2l}}{1-a_k}, \qquad \kappa_k^{(l)} := \frac{(\Sigma_t)_{u,kk}}{(\Sigma_t^{(l)})_{u,kk}} = \frac{1}{1-a_k^{2l}},
$$
**for every eigen-direction $k$, every $l\ge1$, with no commuting assumption on $L,\kappa_{\mathrm{true}}$ anywhere.** Lemma 1's $M,C$ for the specific pair $(\Sigma_t,\Sigma_t^{(l)})$ are, with the correct recursion, exactly $A$ and $C$ themselves — not a separate object requiring its own computation, as an earlier (buggy) pass concluded.

**Consequence for the $J_\pi,J_Q$ machinery.** Since $\Sigma_t,\Sigma_t^{(l)}$ (and, by the same recursion, every intermediate $\Sigma_t^{(m)}$) share the *same* eigenbasis $u$, the sum-over-eigen-directions closed forms hold exactly, with no coupled-ellipsoid obstruction:
$$
J_\pi(t) = \sum_{k=1}^d \delta_k\,s_{t,k}^2\,\Sigma_{t,k}\,z_{0,k}\,\varphi(z_{0,k}), \qquad
J_Q^{(m)}(t) = \sum_{k=1}^d \frac{\delta_kB_{m,k}^2c_k}{\bar S_{m,k}^3}\,\varphi\!\left(\frac{c_k}{\bar S_{m,k}}\right)\Big[V_{t,k}^2+2V_{t,k}\sigma_{m-1,k}^2\Big],
$$
each term built from eigen-direction $k$'s own scalar quantities, in the shared $u$-basis — genuinely separable, not merely approximately so. (An intermediate attempt at fixing this, using a generalized eigendecomposition of the pair $\Xi_0=\Psi^{-1}+L,\ \Xi_1=\Psi^{-1}+\kappa_{\mathrm{true}}$, was unnecessary — that construction is valid but redundant once the recursion bug above is fixed, since $A$ itself already provides the shared eigenbasis directly, for every $l$, with no extra machinery.)

**What still needs correcting: the $\gamma_l\mid\beta$ propagation box above uses the same buggy $M_t$.** In the "Each object inside this" box earlier in this section, the $\gamma_l\mid\beta$ mean and $\Xi_m$ should be re-derived using $\hat M$ in place of $M_t$ throughout before being relied upon (the $C_m(t),W_m(\gamma)$ unrolling that previously appeared in this section has been removed here, superseded by the clean per-eigen-direction result above). Given every $\Sigma_t^{(m)}$ now shares the single eigenbasis $u$, the non-central chi-square obstruction flagged for $\tilde h_m^c(\beta)$ in $d>1$ likely simplifies substantially — plausible, given this correction, but not yet verified and should be checked directly rather than assumed.

---

A natural question is how much of the above survives beyond the Gaussian case, when each coordinate is independent with a general symmetric log-concave marginal density. The natural analog of a precision-matrix eigenvalue is the **curvature function** $c(x) := -(\log f)''(x) \ge 0$ (constant for Gaussians; generally varying with $x$ otherwise).

### What generalizes unconditionally

For independent-coordinate product densities, the log-likelihood ratio between any two such densities is additive across coordinates, and the total variation distance's crossing region remains the set where this sum is nonnegative. The envelope theorem — differentiating a family of densities while holding the (TV-optimal) crossing region fixed — requires no Gaussian structure and applies to any such family.

### Lemma 2′ (proven): matching reference exactly

**Claim.** Let $f_1,f_2,f_3,f_4$ be symmetric, twice-differentiable, log-concave densities on $\mathbb R$ with curvature functions $c_1,c_2,c_3,c_4$. If $f_2=f_4$ (exactly, everywhere) and $c_2(x)\le c_1(x)\le c_3(x)$ for every $x$, then $\|f_1-f_2\|_{TV}\le\|f_3-f_4\|_{TV}$.

**Proof idea.** Build the interpolating curvature $c_t(x):=(1-t)c_1(x)+tc_3(x)$ and recover $f_t$ via $\log f_t(x)=-\int_0^x\!\int_0^s c_t(u)\,du\,ds + A(t)$. Since $c_t(x)$ is a convex combination of $c_1(x),c_3(x)$, both $\ge c_2(x)$ by hypothesis, $c_t(x)\ge c_2(x)$ for every $t$, so $h_t(x):=\log(f_t(x)/f_2(x))$ is concave and even, giving a crossing region that is always an interval around the mode. Writing $D(x):=\int_0^x\!\int_0^s[c_3(u)-c_1(u)]\,du\,ds\ge0$ (even, increasing in $|x|$), the envelope theorem gives $\phi'(t) = -\mathrm{Cov}_{f_t}(D(X),\mathbb 1_{A_t}(X))$, and since $D(X)$ increases in $|X|$ while $\mathbb 1_{A_t}(X)$ decreases in $|X|$, Chebyshev's covariance inequality gives $\phi'(t)\ge0$ for every $t$. Integrating from $0$ to $1$ gives the claim.

This requires $f_2=f_4$ **exactly**, mirroring how Theorem 1's reduction to a shared reference (via Lemma 1's whitening transform) is what let the Gaussian case handle $\Sigma_2\ne\Sigma_4$ — a reduction with no known analog for general log-concave shapes, since rotating or rescaling a non-Gaussian coordinate does not preserve the curvature structure in a controllable way.

### Disproof: pointwise curvature-ratio ordering is not sufficient

The natural candidate condition, matching Lemma 2's convention that pair 2 ($f_3,f_4$) is hypothesized to hold the larger quantity: $c_2(x)/c_1(x) \ge c_4(x)/c_3(x)$ for every $x$ — equivalently, writing $\rho_1(x):=c_1(x)/c_2(x)$ (pair 1's own ratio) and $\rho_2(x):=c_3(x)/c_4(x)$ (pair 2's own ratio), the condition $\rho_2(x)\ge\rho_1(x)$ for every $x$, hoping this gives $\|f_1-f_2\|_{TV}\le\|f_3-f_4\|_{TV}$. This is the pointwise analog of the ratio-domination condition that succeeds in the matrix/Gaussian setting (Lemma 2–3, §3 of the companion supplement), and is the natural thing to hope generalizes, since ratio (not difference, and not raw curvature value) is the correct Gaussian invariant.

**Initial tests were consistent with this condition holding.** Several deliberately adversarial constructions — different shapes for $c_2$ versus $c_4$, ratio functions concentrating their advantage near or away from the mode, and one test with the ratio margin as small as $\varepsilon=10^{-6}$ at a single point — all satisfied $\rho_2\ge\rho_1$ everywhere while also satisfying $\|f_1-f_2\|_{TV}\le\|f_3-f_4\|_{TV}$, suggesting (initially) that pointwise ratio ordering might hold in general.

**A clean, robust counterexample.** Take
$$c_2(x) = 1+\tanh(x^2), \qquad c_1(x) = \rho_1\,(1+\tanh(x^2)) \quad(\rho_1 \text{ constant, bounded/saturating shape}),$$
$$c_4(x) = 1+x^4, \qquad c_3(x) = 2\,(1+x^4) \quad(\rho_2(x)\equiv 2, \text{ constant, "wild" quartic-growth shape}).$$
Sweeping pair 1's ratio $\rho_1$ up toward pair 2's ratio $\rho_2=2$ from below (so that $\rho_2\ge\rho_1$ holds throughout, as required):

| $\rho_1$ | $\|f_1-f_2\|_{TV}$ | $\|f_3-f_4\|_{TV}$ | $\rho_2\ge\rho_1$ holds? | Hoped ordering $\|f_1-f_2\|\le\|f_3-f_4\|$? |
|---|---|---|---|---|
| $1.500$ | $0.089299$ | $0.132463$ | Yes | Yes |
| $1.800$ | $0.129109$ | $0.132463$ | Yes | Yes |
| $1.900$ | $0.140851$ | $0.132463$ | Yes | **No** |
| $1.950$ | $0.146481$ | $0.132463$ | Yes | **No** |
| $1.990$ | $0.150876$ | $0.132463$ | Yes | **No** |
| $1.999$ | $0.151852$ | $0.132463$ | Yes | **No** |
| $1.9999$ | $0.151949$ | $0.132463$ | Yes | **No** |

The transition occurs between $\rho_1=1.8$ (holds) and $\rho_1=1.9$ (fails), and the violation is substantial and stable, not a numerical artifact of an infinitesimal margin: at $\rho_1=1.9$, a full $10\%$ gap remains between $\rho_2=2$ and $\rho_1$ (pair 2 still comfortably holds the larger ratio, as required), yet $\|f_1-f_2\|_{TV}$ already exceeds $\|f_3-f_4\|_{TV}$ by roughly $6\%$, growing to over $14\%$ as $\rho_1\to2^-$. So **pair 1 — the pair with the smaller ratio — ends up with the larger true distance**, exactly reversing the hoped conclusion. This is confirmed stable across truncation ranges $L=100,150,250$.

**Why it fails.** Pointwise curvature-*ratio* ordering exactly matches Lemma 1's invariant when $f_2=f_4$ share the same shape (constant curvature, or more generally identical functional form) — in that regime the ratio alone determines total variation distance, and ordering it suffices. But once $f_2$ and $f_4$ have genuinely different shapes (here, a bounded/saturating one versus quartic tail decay), the *shape itself* — not merely the ratio at each point — independently affects how much separation a given ratio produces: at a matched ratio of $2$, the bounded shape produces $TV=0.15196$ while the wild shape produces only $TV=0.13246$ (verified directly). When pair 2's ratio advantage over pair 1 is comfortably large ($\rho_2/\rho_1 \gtrsim 1.1$ in this example), the ratio effect dominates and the naive ordering happens to hold; but as pair 1's ratio approaches pair 2's, the shape difference becomes the deciding factor, and the pair with the shape that produces more separation per unit ratio (here, the bounded shape) wins out even while holding the smaller ratio. No pointwise scalar comparison at each $x$ captures enough about the *global* relationship between two differently-shaped log-concave densities to determine the ordering of their true total variation distances from two different references.

**Note: this counterexample does not additionally satisfy gap domination.** One might hope the failure above is an artifact of testing ratio domination in isolation, and that also requiring $c_3(x)-c_4(x)\ge c_1(x)-c_2(x)$ everywhere (gap domination, tested on its own earlier and disproven separately in the companion supplement) would rule it out. It does not need to — because this specific counterexample already fails the gap condition at exactly the values where the TV violation occurs. With $c_1=\rho_1(1+\tanh(x^2))$, $c_2=1+\tanh(x^2)$, $c_3=2(1+x^4)$, $c_4=1+x^4$, the gap difference $\big[c_3(x)-c_4(x)\big]-\big[c_1(x)-c_2(x)\big] = (1+x^4) - (\rho_1-1)(1+\tanh(x^2))$ dips **negative** near $x\approx0.6$ for every $\rho_1\ge1.9$ (e.g., a minimum of $-0.151$ at $x\approx0.64$ when $\rho_1=1.95$): $1+x^4$ grows very slowly right off the origin, while $\tanh(x^2)$ rises quickly, so pair 1's gap can briefly exceed pair 2's even while pair 2's *ratio* stays dominant throughout. This counterexample therefore isolates the failure of ratio domination alone; it says nothing about the stronger combined hypothesis, addressed next.

### A stronger disproof: ratio domination and gap domination together are still not sufficient

A natural strengthening is to require **both** conditions to hold simultaneously — pair 2 having the larger ratio *and* the larger gap, everywhere — hoping the combination is restrictive enough to force the correct ordering. It is not.

**A combined counterexample.** Replace pair 2's reference with $c_4(x)=1+x^2$ (in place of the earlier $1+x^4$), keeping the ratio fixed at the same value:
$$c_2(x) = 1+\tanh(x^2), \qquad c_1(x) = \rho_1\,(1+\tanh(x^2)), \qquad c_4(x) = 1+x^2, \qquad c_3(x) = 2\,(1+x^2).$$
At $\rho_1=1.95$:

- **Ratio domination holds:** $\rho_2=2 \ge \rho_1=1.95$.
- **Gap domination holds everywhere on $\mathbb R$:** $\big[c_3(x)-c_4(x)\big]-\big[c_1(x)-c_2(x)\big] = (1+x^2) - 0.95\,(1+\tanh(x^2))$ attains its minimum value $0.050$ exactly at $x=0$ and is strictly positive for every $x\ne0$ (the gap only widens further out, since $c_3-c_4=1+x^2\to\infty$ while $c_1-c_2\to 2(\rho_1-1)$, a bounded constant). This can be seen directly from the small-$x$ expansion $\tanh(x^2)\approx x^2$, giving $\big[c_3-c_4\big]-\big[c_1-c_2\big] \approx (1+x^2)(2-\rho_1) > 0$ for every $\rho_1<2$, with no dip away from the origin — confirmed by a fine numerical grid over $x\in[0,10]$.
- **Yet the true distances are still ordered the wrong way:**
$$\|f_1-f_2\|_{TV} = 0.146481 \;>\; \|f_3-f_4\|_{TV} = 0.143588.$$

Pair 2 — holding both the strictly larger ratio and the strictly larger gap, everywhere on the real line, with clean positive margins in both — still produces the **smaller** true distance. Imposing both natural pointwise conditions simultaneously is therefore still not sufficient.

### An even stronger disproof: restricting pair 2 to be Gaussian does not help either

A further natural attempt to salvage the result is to restrict pair 2 itself to be Gaussian — hoping that removing pair 2's freedom to have an unusual shape (and hence pinning $\|f_3-f_4\|_{TV}$ to the classical, well-understood value $\psi_1(\rho_2)$ of Lemma 1) is enough to force the ordering, since Gaussian intuitively feels like a "typical" or "maximal" shape. It is not: pair 1 alone, using a smooth, everywhere-log-concave shape with **no curvature singularity**, still achieves more separation per unit ratio than any Gaussian can.

**The key building block.** Consider the curvature family $c(x;\kappa) = \dfrac{\kappa}{1+x^2}$ (finite and positive everywhere, including at $x=0$, unlike a raw power-law curvature which would blow up at the origin). Integrating twice gives the closed form
$$f(x;\kappa) = \frac{1}{Z(\kappa)}\exp\!\Big(-\kappa\big[x\arctan(x) - \tfrac12\ln(1+x^2)\big]\Big),$$
verified by finite differences to satisfy $(\log f)''(x) = -\kappa/(1+x^2)$ exactly at every point tested. This density is symmetric, unimodal, strictly log-concave everywhere (curvature is strictly positive for all $x$), and its tails are exponential-type — heavier than Gaussian's — since the curvature *decreases* toward $0$ as $|x|\to\infty$ rather than staying constant. At matched ratio $\kappa_1/\kappa_2$, this shape produces **more** total variation distance than a Gaussian pair at the same ratio (e.g., at ratio $2$: $0.196551$ for this shape versus $0.166064$ for Gaussian; at ratio $1.5$: $0.117458$ versus $0.097776$; at ratio $3$: $0.300601$ versus $0.259359$ — confirmed at every ratio tested).

**The counterexample.** Let pair 1 use this family with $c_2(x) = \dfrac{1}{1+x^2}$, $c_1(x) = \dfrac{1.95}{1+x^2}$ (so $\rho_1(x)\equiv1.95$ exactly, a genuine scale family), and force pair 2 to be **exactly Gaussian**: $c_4(x)=2$, $c_3(x)=4$ (so $\rho_2=2$; the overall level $2$, rather than $1$, is chosen only to make the *absolute* gap large enough — by Lemma 1, this rescaling leaves $\|f_3-f_4\|_{TV}$ unchanged, since it depends only on the ratio $\rho_2$).

- **Ratio domination holds:** $\rho_2 = 2 \ge \rho_1 = 1.95$.
- **Gap domination holds everywhere:** $c_1(x)-c_2(x) = \dfrac{0.95}{1+x^2}$, maximized at $x=0$ with value $0.95$ and strictly decreasing toward $0$ as $|x|\to\infty$; since $c_3-c_4=2 \ge 0.95$, the Gaussian pair's constant gap dominates pair 1's gap at every $x$, with comfortable margin even at the worst point $x=0$.
- **The true distances are nonetheless ordered the wrong way:**
$$\|f_1-f_2\|_{TV} = 0.189756 \;>\; \|f_3-f_4\|_{TV} = 0.166064.$$

So even after eliminating pair 2's freedom to have an unusual shape entirely — forcing it to be the "textbook" Gaussian case where Lemma 1 gives an exact, classical formula for its distance — pair 1's shape alone is enough to produce a larger true distance from a strictly smaller ratio, with gap domination holding simultaneously and with clean positive margins throughout. This rules out the natural remaining hope that Gaussian shapes are somehow extremal (maximizing separation for a given ratio) among log-concave densities: they are not, and smooth, heavier-tailed alternatives can always produce strictly more separation at any fixed ratio, which is precisely why no combination of pointwise conditions on $c_1,\dots,c_4$ tested so far is sufficient once $f_2\ne f_4$ is allowed to vary in shape.

**Standing conclusion.** Pointwise curvature-ratio domination, gap domination, their conjunction, and even their conjunction under the further restriction that pair 2 be Gaussian, are all disproven by explicit, robust, verified counterexamples using genuinely log-concave densities over the entire real line. Lemma 2′ (requiring $f_2=f_4$ exactly) remains the only proven route to comparing true total variation distances for non-Gaussian symmetric log-concave coordinates; extending it to $f_2\ne f_4$ appears to require either a shape-specific reduction analogous to Lemma 1's whitening transform, or a fundamentally different proof technique, neither of which is currently available.

---

## Applied Context: Gradient and Hessian Structure for GLMM Posteriors and Samplers

The disproofs above were motivated by a concrete applied question: comparing the total variation distance between a GLMM's restricted posterior and its restricted sampler, against the analogous distance for a hypothetical model in which the likelihood is replaced by a Gaussian using a lower bound on the data precision. This section records the exact gradient/Hessian structure underlying that comparison, derived directly from the joint negative log-density, and the resulting ratio formula that any future sufficient condition would need to control.

### Setup

Let $-\log\pi(\gamma,\beta\mid y) = \tfrac12(\gamma-\mu_0)^\top\Lambda_\gamma(\gamma-\mu_0) + \tfrac12(\beta-H\gamma)^\top\Psi^{-1}(\beta-H\gamma) - \log L(\beta) + \text{const}$ (design/random-effects terms $H,\Psi$ block-stacked over groups). Define
$$P_{11} := \Lambda_\gamma + H^\top\Psi^{-1}H, \qquad P_{12} := H^\top\Psi^{-1} \ (=P_{21}^\top),$$
the precision of $\gamma\mid\beta$ (exact, since the likelihood does not depend on $\gamma$ directly).

### Posterior: gradient and Hessian

Differentiating $\log\pi(\gamma\mid y)$ under the integral defining the marginal likelihood (a score-function identity) gives
$$\nabla_\gamma\log\pi(\gamma\mid y) = -\Lambda_\gamma(\gamma-\mu_0) + P_{12}\big(\mathbb E[\beta\mid\gamma,y]-H\gamma\big).$$
Setting this to zero shows the mode $\gamma^\star$ is a **linear function of $\mathbb E[\beta\mid\gamma^\star,y]$**:
$$\gamma^\star = P_{11}^{-1}\Big[\Lambda_\gamma\mu_0 + P_{12}\,\mathbb E[\beta\mid\gamma^\star,y]\Big].$$
Differentiating again, using $\dfrac{d}{d\gamma}\mathbb E[\beta\mid\gamma,y] = \mathrm{Cov}(\beta\mid\gamma,y)\,\Psi^{-1}H$:
$$\boxed{c_{\mathrm{post}}(\gamma) := -\nabla^2_\gamma\log\pi(\gamma\mid y) = P_{11} - P_{12}\,V(\gamma)\,P_{21}, \qquad V(\gamma):=\mathrm{Cov}(\beta\mid\gamma,y).}$$

### Sampler kernel: gradient and Hessian

$\gamma\mid\beta,y \sim N\big(m(\beta),P_{11}^{-1}\big)$ with $m(\beta)=P_{11}^{-1}[\Lambda_\gamma\mu_0+P_{12}\beta]$, linear in $\beta$ directly. The marginal one-step kernel is $q(\gamma'\mid\gamma)=\mathbb E_{\beta\sim\pi(\beta\mid\gamma,y)}\big[N(\gamma';m(\beta),P_{11}^{-1})\big]$. Running the identical differentiation, with the **bridge law** $\pi(\beta\mid\gamma,\gamma',y)\propto N(\gamma';m(\beta),P_{11}^{-1})\,\pi(\beta\mid\gamma,y)$ playing the conditioning role:
$$\nabla_{\gamma'}\log q(\gamma'\mid\gamma) = -P_{11}\big(\gamma'-\mathbb E[m(\beta)\mid\gamma,\gamma',y]\big),$$
the same linear-in-conditional-mean structure, and differentiating once more (using linearity of $m(\beta)$ in $\beta$ to push the covariance through, $\mathrm{Cov}(m(\beta)\mid\cdot)=P_{11}^{-1}P_{12}\,U(\gamma,\gamma')\,P_{21}P_{11}^{-1}$):
$$\boxed{c_{\mathrm{sampler}}(\gamma,\gamma') := -\nabla^2_{\gamma'}\log q(\gamma'\mid\gamma) = P_{11} - P_{12}\,U(\gamma,\gamma')\,P_{21}, \qquad U(\gamma,\gamma'):=\mathrm{Cov}(\beta\mid\gamma,\gamma',y).}$$

**Both curvatures share the identical form $P_{11}-P_{12}(\cdot)P_{21}$**, differing only in which distribution the covariance is taken under: $V(\gamma)$ conditions on $\gamma$ alone; $U(\gamma,\gamma')$ conditions on the bridge (both endpoints).

### The ratio structure

Writing $c_{\mathrm{post}}^{\mathrm{hyp}}(\gamma) = P_{11}-P_{12}V^{\mathrm{hyp}}(\gamma)P_{21}$ and $c_{\mathrm{sampler}}^{\mathrm{hyp}}(\gamma,\gamma')=P_{11}-P_{12}U^{\mathrm{hyp}}(\gamma,\gamma')P_{21}$ for the hypothetical (Gaussian-likelihood-at-the-floor) model — using the *same* $P_{11},P_{12}$, since those depend only on the prior and random-effects layer, not the likelihood — the ratios needed for the Lemma-2-extension machinery are
$$\rho_{\mathrm{true}}(\gamma,\gamma') = \frac{c_{\mathrm{sampler}}(\gamma,\gamma')}{c_{\mathrm{post}}(\gamma')} = \frac{P_{11}-P_{12}U(\gamma,\gamma')P_{21}}{P_{11}-P_{12}V(\gamma')P_{21}}, \qquad \rho_{\mathrm{hyp}}(\gamma,\gamma') = \frac{P_{11}-P_{12}U^{\mathrm{hyp}}(\gamma,\gamma')P_{21}}{P_{11}-P_{12}V^{\mathrm{hyp}}(\gamma')P_{21}}.$$
Since $P_{11},P_{12}$ cancel out of the comparison entirely, **the entire question reduces to how $U$ compares to $V$ under the true model, versus how $U^{\mathrm{hyp}}$ compares to $V^{\mathrm{hyp}}$ under the hypothetical model** — a statement purely about the relative gap between bridge-conditioning and single-point-conditioning, in each model separately.

### What is known, and what is not

**Known, and provable without further work.** Since the hypothetical model replaces the likelihood with a strictly less informative Gaussian (using the lower-bound precision $L$ rather than the true, larger likelihood curvature), the resulting posterior uncertainty about $\beta$ is larger under the hypothetical model at every conditioning level:
$$V^{\mathrm{hyp}}(\gamma) \;\succeq\; V(\gamma), \qquad U^{\mathrm{hyp}}(\gamma,\gamma') \;\succeq\; U(\gamma,\gamma') \qquad \text{for every } \gamma,\gamma'.$$
By congruence with the same $P_{12},P_{21}$, this gives $c_{\mathrm{post}}^{\mathrm{hyp}}(\gamma)\preceq c_{\mathrm{post}}(\gamma)$ and $c_{\mathrm{sampler}}^{\mathrm{hyp}}(\gamma,\gamma')\preceq c_{\mathrm{sampler}}(\gamma,\gamma')$ pointwise — **exactly the standing "positive perturbation from a floor" hypothesis** ($\tilde c\ge0$) tested throughout the disproofs above, now derived as a genuine structural fact of the model rather than an assumption.

**Not known, and this is the crux of the difficulty.** The ratio comparison $\rho_{\mathrm{true}}$ vs. $\rho_{\mathrm{hyp}}$ depends on the *relative* gap between bridge- and single-point-conditioning in each model — i.e., on quantities like $V(\gamma')-U(\gamma,\gamma')$ (how much extra information the bridge reveals, under the true model) versus the same difference under the hypothetical model. The pointwise domination facts above say nothing about these gaps, only about the absolute levels $V,U$ versus $V^{\mathrm{hyp}},U^{\mathrm{hyp}}$. Moreover, since both models are restricted to $\widetilde B(\delta_2)$, none of $V,U,V^{\mathrm{hyp}},U^{\mathrm{hyp}}$ have closed forms even in the Gaussian-hypothetical case — truncating a Gaussian to a ball shrinks its covariance by an amount that itself depends on $\gamma,\gamma'$ through where the (untruncated) conditional mass sits relative to the ball, exactly the mechanism flagged earlier as producing the localized, floor-touching-at-the-boundary perturbation shape that broke every pointwise domination condition tested. **This is precisely the missing ingredient**: knowing $c_{\mathrm{true}}\succeq c_{\mathrm{hyp}}$ pointwise (now proven) is not sufficient — as every disproof above demonstrates — without additional control on the *shape* of the gap $c_{\mathrm{true}}-c_{\mathrm{hyp}}$, and the restricted-sampler setting is exactly the case where that shape is not available in closed form.
