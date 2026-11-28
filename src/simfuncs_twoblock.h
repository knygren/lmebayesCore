// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-

/**
 * @file simfuncs_twoblock.h
 * @brief Indexed two-block Gibbs sweep (groups / population), staged in
 *        lmebayesCore for migration to glmbayesCore.
 *
 * @namespace glmbayes::sim::twoblock
 *
 * @section Layout
 *   Names are resolved once in R (.two_block_spec / .two_block_init). Every
 *   step below works on plain matrices, vectors, and 1-based integer index
 *   vectors held in a TwoBlockIndexSpec; indices are shifted only at the
 *   point of use. Each step has an R twin with the same inputs, output
 *   shape, and loop order:
 *
 *     C++                    R
 *     groups_prior_mean      .two_block_groups_prior_mean
 *     groups_draw            .two_block_groups_draw
 *     groups_place           .two_block_groups_place
 *     population_response    .two_block_population_response
 *     population_draw        .two_block_population_draw
 *     population_write       .two_block_population_write
 *     two_block_sweep        .two_block_sweep_r
 *
 *   State layout, identical in both languages:
 *     b            numeric array  J x p_re x n (column-major)
 *     fixef        list of p_re numeric matrices, n x p_k (stored column order)
 *     iters        numeric matrix n x p_re
 *     iters_ranef  numeric vector n
 *
 *   In this header, `i` (chain), `k` (component), and `m` (sweep) are
 *   0-based.
 *
 * @section Implementation
 *   - two_block_index_spec.cpp
 *   - two_block_groups_indexed.cpp
 *   - two_block_population_indexed.cpp
 *   - two_block_population_chol.cpp
 *   - two_block_sweep.cpp
 *
 * @section Migration
 *   Merge these declarations into glmbayesCore's simfuncs_twoblock.h (same
 *   namespace). See inst/DEV_STAGING.md.
 */

#ifndef LMEBAYES_SIMFUNCS_TWOBLOCK_H
#define LMEBAYES_SIMFUNCS_TWOBLOCK_H

#include "RcppArmadillo.h"

#include <string>
#include <vector>

namespace glmbayes {
namespace sim {
namespace twoblock {

/// Plain-number sweep inputs, converted once from the R `spec` list built by
/// `.two_block_spec()`. Index vectors stay 1-based.
struct TwoBlockIndexSpec {
  explicit TwoBlockIndexSpec(const Rcpp::List& spec);

  // Groups block: y ~ D b_j, one row block per group.
  Rcpp::NumericVector y;
  Rcpp::NumericMatrix D;
  Rcpp::NumericVector offset1;
  Rcpp::NumericVector wt1;
  Rcpp::NumericVector disp1;
  Rcpp::NumericMatrix P_D;
  Rcpp::List obs_rows;
  Rcpp::IntegerVector d_col;

  // Population block, one entry per random-effect component k.
  std::vector<Rcpp::NumericMatrix> W;
  std::vector<Rcpp::IntegerVector> w_for_group;
  std::vector<Rcpp::IntegerVector> b_for_wrow;
  std::vector<Rcpp::IntegerVector> w_col;
  std::vector<Rcpp::NumericVector> mu2;
  std::vector<Rcpp::NumericMatrix> P2;
  Rcpp::NumericVector disp2;

  Rcpp::Function f2;
  Rcpp::Function f3;
  Rcpp::Function f2_gauss;
  Rcpp::Function f3_gauss;
  std::string family;
  std::string link;
  int gridtype;
  int n_envopt;

  // Ridge prior for the groups block, from .two_block_ridge_attach().
  // When `ridge` is set, group j is drawn with prior precision ridge_P[[j]]
  // and prior mean ridge_A[, , j] %*% mu_D[, j] + ridge_c[, j].
  // ridge_A is an ncol(D) x ncol(D) x J array, ridge_c is ncol(D) x J.
  bool ridge;
  Rcpp::NumericVector ridge_A;
  Rcpp::NumericMatrix ridge_c;
  Rcpp::List ridge_P;

  int l1;    // ncol(D)
  int J;     // number of groups
  int p_re;  // number of components
};

/// Groups-block draw for one chain: `coef` and `coef_mode` are J x ncol(D);
/// `iters` is the mean envelope candidate count over groups.
struct GroupsDraw {
  Rcpp::NumericMatrix coef;
  Rcpp::NumericMatrix coef_mode;
  double iters;
};

/// Population-block draw for one component, in W_k column order.
struct PopulationDraw {
  Rcpp::NumericVector coef;
  Rcpp::NumericVector coef_mode;
};

/// Optional hooks into two_block_sweep. The default does nothing, so a
/// null observer and a default observer give the same run.
class SweepObserver {
 public:
  virtual ~SweepObserver() {}
  /// Return true and fill `coef` (J x ncol(D)) and `iters` to replace the
  /// groups draw of chain `i` in sweep `m`.
  virtual bool groups_source(int m, int i, Rcpp::NumericMatrix& coef,
                             double& iters) { return false; }
  /// Return true and fill `coef` (length p_k) to replace the population
  /// draw of component `k`, chain `i`, sweep `m`.
  virtual bool population_source(int m, int i, int k,
                                 Rcpp::NumericVector& coef) { return false; }
  virtual void groups_seen(int m, int i, const Rcpp::NumericMatrix& mu_D,
                           const Rcpp::NumericMatrix& coef, double iters) {}
  virtual void population_seen(int m, int i, int k,
                               const Rcpp::NumericVector& y_k,
                               const Rcpp::NumericVector& coef) {}
};

/// Groups prior mean for chain `i`: ncol(D) x J. Row d_col[k], column j is
/// sum_c W_k[w_for_group_k[j], c] * fixef_k[i, w_col_k[c]], summed over c in
/// increasing order.
Rcpp::NumericMatrix groups_prior_mean(const TwoBlockIndexSpec& spec,
                                      const Rcpp::List& fixef, int i);

/// Ridge prior mean, ncol(D) x J: column j is
/// ridge_A[, , j] %*% mu_D[, j] + ridge_c[, j], summed over the column of
/// ridge_A in increasing order. Requires `spec.ridge`.
Rcpp::NumericMatrix groups_ridge_mean(const TwoBlockIndexSpec& spec,
                                      const Rcpp::NumericMatrix& mu_D);

/// One groups draw for every group, given the prior mean `mu_D`. With
/// `spec.ridge`, the draw uses the ridge prior built from `mu_D`.
GroupsDraw groups_draw(const TwoBlockIndexSpec& spec,
                       const Rcpp::NumericMatrix& mu_D);

/// One group's rows of the groups block and its prior precision. These do
/// not change across sweeps, so they are sliced once per two_block_sweep.
struct GroupData {
  Rcpp::NumericVector y;
  Rcpp::NumericMatrix x;
  Rcpp::NumericVector offset;
  Rcpp::NumericVector wt;
  Rcpp::NumericMatrix P;
  double disp;
};

/// GroupData for every group, from obs_rows, offset1, wt1, disp1 and P_D
/// (ridge_P[[j]] when `spec.ridge`).
std::vector<GroupData> groups_data(const TwoBlockIndexSpec& spec);

/// Groups prior means for all chains: an ncol(D) x J x n array (ridge mean
/// when `spec.ridge`). Element (a, j, i) is at a + l1 * (j + J * i).
/// `mu_D`, when given, receives groups_prior_mean for each chain (before the
/// ridge mean), as passed to SweepObserver::groups_seen.
Rcpp::NumericVector groups_prior_mean_all(
    const TwoBlockIndexSpec& spec, const Rcpp::List& fixef, int n,
    std::vector<Rcpp::NumericMatrix>* mu_D = nullptr);

/// Group `j` drawn for the chains in `chains`: `coef` and `coef_mode` are
/// chains.size() x ncol(D), `iters` has one entry per chain. `timing` is
/// seconds in mode, setup, eval, finalize, simulate (rNormalGLMGroupChains).
/// `setup_timing` and `mode_timing` split setup and mode finding into their
/// parts, `setup_compare` counts the
/// R / C++ grid setup comparisons and `mode_compare` is the optim / vmmin
/// comparison (all named; empty when the installed glmbayesCore does not
/// return them). `mode_hessian` is ncol(D)^2 x chains.size() (column c is
/// chain c's mode Hessian A1, column-major), `warm_compare` the warm /
/// cold start comparison, `hessian_compare` the exact / finite-difference
/// Hessian comparison, `newton_compare` the Newton / vmmin comparison,
/// `sim_compare` the serial / parallel simulation comparison and
/// `eval_compare` the legacy / batch envelope evaluation comparison.
struct GroupChainsDraw {
  Rcpp::NumericMatrix coef;
  Rcpp::NumericMatrix coef_mode;
  Rcpp::NumericMatrix mode_hessian;
  Rcpp::NumericVector iters;
  Rcpp::NumericVector timing;
  Rcpp::NumericVector setup_timing;
  Rcpp::NumericVector mode_timing;
  Rcpp::NumericVector setup_compare;
  Rcpp::NumericVector mode_compare;
  Rcpp::NumericVector warm_compare;
  Rcpp::NumericVector hessian_compare;
  Rcpp::NumericVector newton_compare;
  Rcpp::NumericVector sim_compare;
  Rcpp::NumericVector eval_compare;
};

/// One draw of group `j` per chain in `chains`, prior means taken from
/// `mu_all` (groups_prior_mean_all). `check` compares each chain's mode and
/// envelope with rNormalGLM and stops on a difference (validation only).
/// `mode_method`: 0 optim, 1 vmmin with C++ f2/f3, 2 optim compared against
/// vmmin, 3 Newton (C++ f2/f3 + exact Hessian), 4 vmmin compared against 3,
/// 5 Newton (fam f2_f3_f4), 6 vmmin compared against 5. `setup_method`: 0 R grid setup, 1 C++, 2 both compared exactly
/// (rNormalGLMGroupChains). `start` is empty or ncol(D) x chains.size()
/// starting modes (a column with an NA starts cold); `warm_method`: 0 cold,
/// 1 warm start, 2 both compared (cold kept). `hess_method`: 0 finite
/// differences, 1 exact (finite differences if not finite or not PD),
/// 2 both compared (finite differences kept).
GroupChainsDraw group_chains_draw(const TwoBlockIndexSpec& spec,
                                  const GroupData& gd,
                                  const Rcpp::NumericVector& mu_all, int j,
                                  const std::vector<int>& chains,
                                  bool check = false, int mode_method = 0,
                                  int setup_method = 0,
                                  const Rcpp::NumericMatrix& start =
                                      Rcpp::NumericMatrix(0, 0),
                                  int warm_method = 0, int hess_method = 0);

/// Columns of a groups draw in stored `b` order: J x p_re.
Rcpp::NumericMatrix groups_place(const TwoBlockIndexSpec& spec,
                                 const Rcpp::NumericMatrix& coef);

/// Population response for component `k`, in the row order of W_k.
Rcpp::NumericVector population_response(const TwoBlockIndexSpec& spec,
                                        const Rcpp::NumericMatrix& b_i, int k);

/// One population draw for component `k`.
PopulationDraw population_draw(const TwoBlockIndexSpec& spec,
                               const Rcpp::NumericVector& y_k, int k);

/// Block-2 draw path from getOption("lmebayesCore.pop_draw"): "reg" calls
/// population_draw; "chol" (or unset) calls population_draw_chol;
/// "compare" runs both from the same R RNG state, keeps the population_draw
/// result and fills PopulationCompare. Uses the R API.
enum PopMethod {
  POP_REG = 0,
  POP_CHOL = 1,
  POP_COMPARE = 2
};
int population_draw_method();

/// Posterior of component `k` with W_k, P2_k, mu2_k and disp2_k fixed:
/// precision Q = W' W / disp2 + P2 = R' R (R upper triangular). Rinv is
/// R^{-1}, Wt is W' / disp2 and Pmu is P2 mu2, so the mean for response y
/// is Rinv Rinv' (Wt y + Pmu).
struct PopulationFactor {
  int q = 0;
  int n_w = 0;
  arma::mat Rinv;
  arma::mat Wt;
  arma::vec Pmu;
};
PopulationFactor population_factor(const TwoBlockIndexSpec& spec, int k);

/// The population_draw posterior from a PopulationFactor: coef_mode is the
/// mean, coef is mean + Rinv z with z = q norm_rand() draws from R's RNG,
/// the same count and order as glmbayesCore rNormalReg.
PopulationDraw population_draw_chol(const PopulationFactor& f,
                                    const Rcpp::NumericVector& y_k);

/// Counts and maxima accumulated with POP_COMPARE: calls, flagged
/// (max_mean or max_draw above 1e-10), max_mean and max_draw
/// (|chol - reg| / (1 + |reg|) for coef_mode and coef), and seconds in each
/// path.
struct PopulationCompare {
  double calls = 0.0;
  double flagged = 0.0;
  double max_mean = 0.0;
  double max_draw = 0.0;
  double t_reg = 0.0;
  double t_chol = 0.0;
};

/// population_draw and population_draw_chol from the same R RNG state
/// (.Random.seed saved and restored between them); returns the
/// population_draw result and leaves the RNG where that draw left it.
PopulationDraw population_draw_compare(PopulationCompare& pc,
                                       const TwoBlockIndexSpec& spec,
                                       const PopulationFactor& f,
                                       const Rcpp::NumericVector& y_k, int k);

/// fixef_k[i, w_col_k] <- coef, in place.
void population_write(const TwoBlockIndexSpec& spec, Rcpp::NumericMatrix& fixef_k,
                      int i, int k, const Rcpp::NumericVector& coef);

/// Progress prefix "[stage] sweep m/M phase: "; an empty `phase` gives the
/// per-sweep prefix used for a single chain. `m` is 1-based.
std::string sweep_progbar_prefix(const std::string& stage_label, int m,
                                 int inner_sweeps, const std::string& phase);

/// `inner_sweeps` sweeps over all chains. Inputs are cloned. Returns a list
/// with `b`, `fixef`, `iters`, `iters_ranef` in the layout above, plus the
/// per-sweep cross-chain summary of the stacked population coefficients
/// (components in order, columns in stored order; p_tot columns in all):
/// `sweep_mean` (p_tot x M) and `sweep_cov` (p_tot x p_tot x M, denominator
/// n - 1, NA when n < 2).
/// `use_mode` keeps each block's posterior mode instead of its draw (used
/// only for exact R / C++ comparison).
/// `rwt` in [0, 1) is the ridge weight: the groups draws use observation
/// weights `wt1 * (1 - rwt)`, i.e. the group likelihood raised to 1 - rwt.
/// The matching ridge term (rwt times the quadratic expansion at the
/// reference point) is in the spec's ridge prior when `spec.ridge` is set.
/// `rwt = 0` leaves the weights unchanged.
/// `group_major` (default) runs the groups block as an outer loop over groups
/// and an inner draw over all chains (group_chains_draw); false keeps the
/// chain-major loop (groups_draw per chain), which matches the R twin's loop
/// order. Both give the same conditional draws; the random number stream is
/// used in a different order. `group_check` is passed to group_chains_draw
/// (validation only, doubles the groups cost).
/// The list also has `timing`: seconds over the call in the block-1 parts
/// mode, setup, eval, finalize, simulate (group-major only; 0 otherwise),
/// then block1 and block2 wall time. `mode_method` is passed to
/// group_chains_draw (group-major only), and so is `setup_method`. The list
/// also has `setup_timing`, `mode_timing` and `setup_compare` (summed over
/// calls) and
/// `mode_compare` (zero unless mode_method is 2): counts and times summed
/// over calls, max_* the maximum.
/// `warm_method` (group-major only): 0 cold starts; 1 starts each chain's
/// group mode search from b_old + A1_old^{-1} P (mu_new - mu_old), the
/// previous sweep's mode, Hessian and prior mean for that chain and group in
/// this call (the first sweep starts cold); 2 runs both and keeps the cold
/// result. The list also has `warm_compare`, combined like mode_compare.
/// `hess_method` (group-major only) is passed to group_chains_draw; the list
/// also has `hessian_compare`, `newton_compare`, `sim_compare` and
/// `eval_compare`, combined
/// the same way.
/// The population draws follow population_draw_method(); the factors are
/// built once per call, and `pop_compare` holds the PopulationCompare fields
/// (all zero unless the method is POP_COMPARE).
Rcpp::List two_block_sweep(const TwoBlockIndexSpec& spec,
                           const Rcpp::NumericVector& b,
                           const Rcpp::List& fixef,
                           const Rcpp::NumericMatrix& iters,
                           const Rcpp::NumericVector& iters_ranef,
                           int inner_sweeps, bool use_mode, bool progbar,
                           const std::string& stage_label,
                           double rwt = 0.0,
                           SweepObserver* observer = nullptr,
                           bool group_major = true,
                           bool group_check = false,
                           int mode_method = 0,
                           int setup_method = 0,
                           int warm_method = 0,
                           int hess_method = 0);

}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes

#endif
