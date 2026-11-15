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

/// One groups draw for every group, given the prior mean `mu_D`.
GroupsDraw groups_draw(const TwoBlockIndexSpec& spec,
                       const Rcpp::NumericMatrix& mu_D);

/// Columns of a groups draw in stored `b` order: J x p_re.
Rcpp::NumericMatrix groups_place(const TwoBlockIndexSpec& spec,
                                 const Rcpp::NumericMatrix& coef);

/// Population response for component `k`, in the row order of W_k.
Rcpp::NumericVector population_response(const TwoBlockIndexSpec& spec,
                                        const Rcpp::NumericMatrix& b_i, int k);

/// One population draw for component `k`.
PopulationDraw population_draw(const TwoBlockIndexSpec& spec,
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
Rcpp::List two_block_sweep(const TwoBlockIndexSpec& spec,
                           const Rcpp::NumericVector& b,
                           const Rcpp::List& fixef,
                           const Rcpp::NumericMatrix& iters,
                           const Rcpp::NumericVector& iters_ranef,
                           int inner_sweeps, bool use_mode, bool progbar,
                           const std::string& stage_label,
                           SweepObserver* observer = nullptr);

}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes

#endif
