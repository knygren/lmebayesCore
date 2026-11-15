// Groups block of the indexed two-block sweep: prior mean, draw, and
// placement into stored `b` order. R twins: R/two_block_groups_indexed.R.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "rng_utils.h"
#include "simfuncs_groups.h"
#include "simfuncs_twoblock.h"

using namespace Rcpp;

namespace glmbayes {
namespace sim {
namespace twoblock {

NumericMatrix groups_prior_mean(const TwoBlockIndexSpec& spec,
                                const List& fixef, int i) {
  const int J = spec.J;
  NumericMatrix mu_D(spec.l1, J);
  for (int k = 0; k < spec.p_re; ++k) {
    const NumericMatrix& Wk = spec.W[k];
    const IntegerVector& wg = spec.w_for_group[k];
    const IntegerVector& wc = spec.w_col[k];
    const NumericMatrix fk = fixef[k];
    NumericVector s(J);
    for (int c = 0; c < Wk.ncol(); ++c) {
      const double g = fk(i, wc[c] - 1);
      for (int j = 0; j < J; ++j) {
        s[j] = s[j] + Wk(wg[j] - 1, c) * g;
      }
    }
    for (int j = 0; j < J; ++j) {
      mu_D(spec.d_col[k] - 1, j) = s[j];
    }
  }
  return mu_D;
}

GroupsDraw groups_draw(const TwoBlockIndexSpec& spec, const NumericMatrix& mu_D) {
  List out;
  {
    glmbayes::rng::RSyncScope rng_sync;
    out = glmbayes::sim::group::rNormalGLMGroups(
        1, spec.y, spec.D, spec.offset1, spec.wt1, spec.disp1, mu_D,
        List::create(spec.P_D), true, spec.obs_rows, spec.f2, spec.f3,
        spec.family, spec.link, spec.gridtype, spec.n_envopt, false, false,
        false);
  }

  const List gr = out["group_results"];
  long double s = 0.0;
  for (int b = 0; b < gr.size(); ++b) {
    const List g = gr[b];
    const NumericVector it = g["iters"];
    s += it[0];
  }
  GroupsDraw d;
  d.coef = as<NumericMatrix>(out["coefficients"]);
  d.coef_mode = as<NumericMatrix>(out["coef.mode"]);
  d.iters = static_cast<double>(s / gr.size());
  return d;
}

NumericMatrix groups_place(const TwoBlockIndexSpec& spec, const NumericMatrix& coef) {
  const int J = coef.nrow();
  NumericMatrix out(J, spec.p_re);
  for (int k = 0; k < spec.p_re; ++k) {
    for (int j = 0; j < J; ++j) {
      out(j, k) = coef(j, spec.d_col[k] - 1);
    }
  }
  return out;
}

}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes
