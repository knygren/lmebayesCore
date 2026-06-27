// two_block_block1.cpp
// Block 1 batch helpers for run_sweep_outer_chains_v6 (R-oracle semantics).
// Migrated incrementally from the R driver in R/two_block_batch_gibbs.R.

#include "simfuncs.h"

namespace glmbayes {
namespace sim {

/// Mean envelope candidate count across groups from one Block 1 draw
/// (port of .two_block_block1_iters_mean in two_block_batch_gibbs.R).
double two_block_block1_iters_mean(const Rcpp::List& block_out) {
  if (!block_out.containsElementNamed("block_results")) {
    return 1.0;
  }
  Rcpp::List br = Rcpp::as<Rcpp::List>(block_out["block_results"]);
  if (br.size() == 0) {
    return 1.0;
  }
  double sum = 0.0;
  for (int b = 0; b < br.size(); ++b) {
    Rcpp::List out_b = Rcpp::as<Rcpp::List>(br[b]);
    if (!out_b.containsElementNamed("iters") || Rf_isNull(out_b["iters"])) {
      sum += 1.0;
      continue;
    }
    SEXP it = out_b["iters"];
    if (Rf_isMatrix(it)) {
      Rcpp::NumericMatrix im(it);
      sum += im(0, 0);
    } else {
      Rcpp::NumericVector iv(it);
      sum += iv[0];
    }
  }
  return sum / static_cast<double>(br.size());
}

} // namespace sim
} // namespace glmbayes
