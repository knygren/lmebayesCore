// DEV ONLY -- does not migrate to glmbayesCore.
// Tape replay, input/draw log, and groups prior-mean sums for the indexed
// two-block sweep, attached to two_block_sweep() as a SweepObserver.
// R twin: R/dev_two_block_replay.R.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "dev_twoblock.h"

using namespace Rcpp;
using glmbayes::sim::twoblock::TwoBlockIndexSpec;

namespace lmebayes {
namespace dev {
namespace twoblock {

NumericVector zero_array(int d1, int d2, int d3, int d4) {
  const int len = d1 * d2 * d3 * (d4 < 0 ? 1 : d4);
  NumericVector out(len, 0.0);
  if (d4 < 0) {
    out.attr("dim") = IntegerVector::create(d1, d2, d3);
  } else {
    out.attr("dim") = IntegerVector::create(d1, d2, d3, d4);
  }
  return out;
}

List log_new(const TwoBlockIndexSpec& spec, int n, int M) {
  List y(spec.p_re);
  List block2(spec.p_re);
  for (int k = 0; k < spec.p_re; ++k) {
    y[k] = zero_array(spec.W[k].nrow(), n, M);
    block2[k] = zero_array(spec.W[k].ncol(), n, M);
  }
  return List::create(
      Named("mu_D") = zero_array(spec.l1, spec.J, n, M),
      Named("block1") = zero_array(spec.J, spec.l1, n, M),
      Named("iters1") = NumericMatrix(n, M),
      Named("y") = y,
      Named("block2") = block2);
}

ReplayObserver::ReplayObserver(const TwoBlockIndexSpec& spec, int n, int M,
                               Nullable<List> tape, bool record)
    : prior_sum(zero_array(spec.l1, spec.J, M)),
      prior_sumsq(zero_array(spec.l1, spec.J, M)),
      l1_(spec.l1),
      J_(spec.J),
      n_(n),
      use_tape_(tape.isNotNull()),
      record_(record) {
  if (use_tape_) {
    const List t(tape);
    tape_b1_ = t["block1"];
    tape_it1_ = as<NumericMatrix>(t["iters1"]);
    tape_b2_ = t["block2"];
  }
  if (record_) {
    log_ = log_new(spec, n, M);
    log_mu_ = log_["mu_D"];
    log_b1_ = log_["block1"];
    log_it1_ = as<NumericMatrix>(log_["iters1"]);
    log_y_ = log_["y"];
    log_b2_ = log_["block2"];
  }
}

bool ReplayObserver::groups_source(int m, int i, NumericMatrix& coef,
                                   double& iters) {
  if (!use_tape_) {
    return false;
  }
  coef = NumericMatrix(J_, l1_);
  for (int c = 0; c < l1_; ++c) {
    for (int j = 0; j < J_; ++j) {
      coef(j, c) = tape_b1_[j + J_ * (c + l1_ * (i + n_ * m))];
    }
  }
  iters = tape_it1_(i, m);
  return true;
}

bool ReplayObserver::population_source(int m, int i, int k, NumericVector& coef) {
  if (!use_tape_) {
    return false;
  }
  const NumericVector t2 = tape_b2_[k];
  const IntegerVector dim = t2.attr("dim");
  const int p_k = dim[0];
  coef = NumericVector(p_k);
  for (int c = 0; c < p_k; ++c) {
    coef[c] = t2[c + p_k * (i + n_ * m)];
  }
  return true;
}

void ReplayObserver::groups_seen(int m, int i, const NumericMatrix& mu_D,
                                 const NumericMatrix& coef, double iters) {
  for (int j = 0; j < J_; ++j) {
    for (int a = 0; a < l1_; ++a) {
      const int idx = a + l1_ * (j + J_ * m);
      prior_sum[idx] = prior_sum[idx] + mu_D(a, j);
      prior_sumsq[idx] = prior_sumsq[idx] + mu_D(a, j) * mu_D(a, j);
    }
  }
  if (!record_) {
    return;
  }
  for (int j = 0; j < J_; ++j) {
    for (int a = 0; a < l1_; ++a) {
      log_mu_[a + l1_ * (j + J_ * (i + n_ * m))] = mu_D(a, j);
    }
  }
  for (int c = 0; c < l1_; ++c) {
    for (int j = 0; j < J_; ++j) {
      log_b1_[j + J_ * (c + l1_ * (i + n_ * m))] = coef(j, c);
    }
  }
  log_it1_(i, m) = iters;
}

void ReplayObserver::population_seen(int m, int i, int k, const NumericVector& y_k,
                                     const NumericVector& coef) {
  if (!record_) {
    return;
  }
  NumericVector ly = log_y_[k];
  NumericVector lb2 = log_b2_[k];
  const int r_k = y_k.size();
  const int p_k = coef.size();
  for (int r = 0; r < r_k; ++r) {
    ly[r + r_k * (i + n_ * m)] = y_k[r];
  }
  for (int c = 0; c < p_k; ++c) {
    lb2[c + p_k * (i + n_ * m)] = coef[c];
  }
}

RObject ReplayObserver::log() const {
  return record_ ? RObject(log_) : RObject(R_NilValue);
}

}  // namespace twoblock
}  // namespace dev
}  // namespace lmebayes
