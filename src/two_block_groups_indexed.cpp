// Groups block of the indexed two-block sweep: prior mean, draw, and
// placement into stored `b` order. R twins: R/two_block_groups_indexed.R.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <glmbayesCore/glmbayesCore_capi.h>

#include "rng_utils.h"
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

NumericMatrix groups_ridge_mean(const TwoBlockIndexSpec& spec,
                                const NumericMatrix& mu_D) {
  const int l1 = spec.l1;
  const int J = spec.J;
  NumericMatrix out(l1, J);
  for (int j = 0; j < J; ++j) {
    const R_xlen_t base = static_cast<R_xlen_t>(l1) * l1 * j;
    for (int a = 0; a < l1; ++a) {
      double s = 0.0;
      for (int b = 0; b < l1; ++b) {
        s = s + spec.ridge_A[base + a + static_cast<R_xlen_t>(l1) * b] * mu_D(b, j);
      }
      out(a, j) = s + spec.ridge_c(a, j);
    }
  }
  return out;
}

GroupsDraw groups_draw(const TwoBlockIndexSpec& spec, const NumericMatrix& mu_D) {
  const RObject n1 = wrap(1);
  const RObject prior_by_block = wrap(true);
  const RObject family = wrap(spec.family);
  const RObject link = wrap(spec.link);
  const RObject gridtype = wrap(spec.gridtype);
  const RObject n_envopt = wrap(spec.n_envopt);
  const RObject off = wrap(false);
  const List P_blocks = spec.ridge ? spec.ridge_P : List::create(spec.P_D);
  const NumericMatrix mu = spec.ridge ? groups_ridge_mean(spec, mu_D) : mu_D;
  List out;
  {
    glmbayes::rng::RSyncScope rng_sync;
    out = glmbayesCore_rNormalGLMGroups(
        n1, spec.y, spec.D, spec.offset1, spec.wt1, spec.disp1, mu,
        P_blocks, prior_by_block, spec.obs_rows, spec.f2, spec.f3,
        family, link, gridtype, n_envopt, off, off, off);
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
