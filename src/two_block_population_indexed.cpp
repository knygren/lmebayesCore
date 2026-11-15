// Population block of the indexed two-block sweep: response, draw, and
// write-back into the stored fixef matrix. R twins:
// R/two_block_population_indexed.R.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "rng_utils.h"
#include "simfuncs.h"
#include "simfuncs_twoblock.h"

using namespace Rcpp;

namespace glmbayes {
namespace sim {
namespace twoblock {

NumericVector population_response(const TwoBlockIndexSpec& spec,
                                  const NumericMatrix& b_i, int k) {
  const IntegerVector& rows = spec.b_for_wrow[k];
  NumericVector y_k(rows.size());
  for (int r = 0; r < rows.size(); ++r) {
    y_k[r] = b_i(rows[r] - 1, k);
  }
  return y_k;
}

PopulationDraw population_draw(const TwoBlockIndexSpec& spec,
                               const NumericVector& y_k, int k) {
  const NumericVector& mu = spec.mu2[k];
  const int n_w = y_k.size();
  List out;
  {
    glmbayes::rng::RSyncScope rng_sync;
    out = glmbayes::sim::rNormalReg(
        1, y_k, spec.W[k], mu, spec.P2[k], NumericVector(n_w, 0.0),
        NumericVector(n_w, 1.0), spec.disp2[k], spec.f2_gauss, spec.f3_gauss,
        mu, "gaussian", "identity", spec.gridtype);
  }
  const NumericMatrix draw = out["coefficients"];
  const NumericVector mode = out["coef.mode"];
  PopulationDraw d;
  d.coef = NumericVector(draw.ncol());
  for (int c = 0; c < draw.ncol(); ++c) {
    d.coef[c] = draw(0, c);
  }
  d.coef_mode = NumericVector(mode.begin(), mode.end());
  return d;
}

void population_write(const TwoBlockIndexSpec& spec, NumericMatrix& fixef_k,
                      int i, int k, const NumericVector& coef) {
  const IntegerVector& wc = spec.w_col[k];
  for (int c = 0; c < wc.size(); ++c) {
    fixef_k(i, wc[c] - 1) = coef[c];
  }
}

}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes
