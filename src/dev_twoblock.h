// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-

/**
 * @file dev_twoblock.h
 * @brief Development-only diagnostics for the indexed two-block sweep.
 *
 * @namespace lmebayes::dev::twoblock
 *
 *   Replay (tape), input/draw logging (record), and groups prior-mean
 *   tracking, all attached to glmbayes::sim::twoblock::two_block_sweep
 *   through a SweepObserver. The production loop is unchanged; nothing here
 *   migrates to glmbayesCore (see inst/DEV_STAGING.md).
 *
 * @section Implementation
 *   - dev_two_block_replay.cpp
 *   - dev_exports.cpp
 */

#ifndef LMEBAYES_DEV_TWOBLOCK_H
#define LMEBAYES_DEV_TWOBLOCK_H

#include "RcppArmadillo.h"
#include "simfuncs_twoblock.h"

namespace lmebayes {
namespace dev {
namespace twoblock {

/// Zero array with dim c(d1, d2, d3) or c(d1, d2, d3, d4).
Rcpp::NumericVector zero_array(int d1, int d2, int d3, int d4 = -1);

/// Empty log for `n` chains and `M` sweeps: `mu_D` (ncol(D) x J x n x M),
/// `block1` (J x ncol(D) x n x M), `iters1` (n x M), and lists `y`
/// (nrow(W_k) x n x M) and `block2` (ncol(W_k) x n x M) over components.
/// A filled log is also a valid tape.
Rcpp::List log_new(const glmbayes::sim::twoblock::TwoBlockIndexSpec& spec,
                   int n, int M);

/// Replays a tape (when given), records a log (when `record`), and always
/// accumulates the groups prior-mean sums `prior_sum` / `prior_sumsq`
/// (ncol(D) x J x M, summed over chains).
class ReplayObserver : public glmbayes::sim::twoblock::SweepObserver {
 public:
  ReplayObserver(const glmbayes::sim::twoblock::TwoBlockIndexSpec& spec,
                 int n, int M, Rcpp::Nullable<Rcpp::List> tape, bool record);

  bool groups_source(int m, int i, Rcpp::NumericMatrix& coef,
                     double& iters) override;
  bool population_source(int m, int i, int k,
                         Rcpp::NumericVector& coef) override;
  void groups_seen(int m, int i, const Rcpp::NumericMatrix& mu_D,
                   const Rcpp::NumericMatrix& coef, double iters) override;
  void population_seen(int m, int i, int k, const Rcpp::NumericVector& y_k,
                       const Rcpp::NumericVector& coef) override;

  /// The log, or NULL when not recording.
  Rcpp::RObject log() const;

  Rcpp::NumericVector prior_sum;
  Rcpp::NumericVector prior_sumsq;

 private:
  int l1_, J_, n_;
  bool use_tape_, record_;
  Rcpp::NumericVector tape_b1_;
  Rcpp::NumericMatrix tape_it1_;
  Rcpp::List tape_b2_;
  Rcpp::List log_;
  Rcpp::NumericVector log_mu_, log_b1_;
  Rcpp::NumericMatrix log_it1_;
  Rcpp::List log_y_, log_b2_;
};

}  // namespace twoblock
}  // namespace dev
}  // namespace lmebayes

#endif
