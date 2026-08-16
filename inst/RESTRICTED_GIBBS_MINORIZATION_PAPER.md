# Restricted Two-Block Gibbs: Minorization on \(C_d\) and Distance to \(\pi\)

## 1. Introduction

### 1.1 Problem and setting

### 1.2 Statement of the main results (informal)

### 1.3 Relation to Rosenthal (1995)

## 2. The Conditional Bound

### 2.1 Two-block Gibbs sampler and the restricted kernel

### 2.2 Total variation and \(L^1\) distance

### 2.3 Truncation: \(\pi(\cdot\mid C_d)\) vs. \(\pi\)

### 2.4 Theorem 1 — geometric convergence of the restricted chain

## 3. The Hierarchical GLMM Model

### 3.1 Model and hypotheses (H1)–(H3)

### 3.2 Propriety of the flat-prior posterior

### 3.3 The mean map and the marginal \(\gamma\)-chain

### 3.4 Theorem 2 — existence of a certified safe set (statement)

**Theorem 2 (existence of a certified safe set).** Under (H1), (H2), (H3a)
— and (H3b) as well if \(\Lambda_\gamma = 0\) — for every \(\delta > 0\)
there exist

- a compact, convex set \(C = \widetilde C_d \subseteq \mathbb R^q\), the
  superlevel set \(\{\gamma : \varepsilon(\gamma) \ge \varepsilon_d\}\) of the
  minorization profile, with \(d = d(\delta)\) chosen so that
  \(\pi_\gamma(C^c) < \delta\);
- a restricted \(\gamma\)-chain kernel \(q(\cdot \mid \cdot\,; C)\) — the
  two-block sweep confined to \(C\);
- a constant \(\varepsilon = \varepsilon_d\, Q(C) \in (0,1]\) and the truncated
  Gaussian \(Q_C = Q(\cdot \mid C)\), where \(Q = N(\gamma^\star, \Sigma^\star)\)
  is the refresh measure of §4,

such that \(\pi_\gamma(\cdot \mid C)\) is stationary for \(q(\cdot \mid \cdot\,; C)\)
and \(q(\gamma, A \mid C) \ge \varepsilon\, Q_C(A)\) for all \(\gamma \in C\),
\(A \subseteq C\) — hence, by Theorem 1, for \(\gamma \in C\),

\[
\|q_n(\gamma,\cdot \mid C) - \pi_\gamma\|_{TV} \;\le\; (1-\varepsilon)^n + \delta.
\]

## 4. The Refresh Measure \(Q\)

### 4.1 Proper population prior

### 4.2 Flat-prior limit

### 4.3 Summary: existence and nondegeneracy of \(Q\)

## 5. The Minorization Constant and the Certified Set

### 5.1 The minorization profile \(\varepsilon(\gamma)\)

### 5.2 Convexity, coercivity, and compactness of \(\widetilde C_d\)

### 5.3 The attained constant \(\varepsilon_d\)

### 5.4 Flat-prior limit

### 5.5 Closed form under Gaussian closure

## 6. Proof of Theorem 2

### 6.1 Assembly of §4–§5

### 6.2 Lifting to the joint \((\gamma,\beta)\) chain

## 7. Scope and Extensions

### 7.1 Symmetric case: sharper constants

### 7.2 What is and is not certified

### 7.3 Open problems

## Appendix — Proofs

### A.1 Proof of §3.2 (propriety of the flat-prior posterior)

### A.2 Proofs for §4

### A.3 Proofs for §5

### A.4 Proof of Theorem 2

## References
