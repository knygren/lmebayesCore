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

std::vector<GroupData> groups_data(const TwoBlockIndexSpec& spec) {
  const int l1 = spec.l1;
  const int J = spec.J;
  const R_xlen_t l2 = spec.y.size();
  if (spec.obs_rows.size() != J) {
    Rcpp::stop("length(obs_rows) must equal the number of groups");
  }
  if (spec.disp1.size() != 1 && spec.disp1.size() != J) {
    Rcpp::stop("length(disp1) must be 1 or the number of groups");
  }
  if (spec.offset1.size() != 1 && spec.offset1.size() != l2) {
    Rcpp::stop("length(offset1) must be 1 or length(y)");
  }
  if (spec.wt1.size() != 1 && spec.wt1.size() != l2) {
    Rcpp::stop("length(wt1) must be 1 or length(y)");
  }
  std::vector<GroupData> out(J);
  for (int j = 0; j < J; ++j) {
    const IntegerVector rows = spec.obs_rows[j];
    const int nj = rows.size();
    GroupData& g = out[j];
    g.y = NumericVector(nj);
    g.x = NumericMatrix(nj, l1);
    g.offset = NumericVector(nj);
    g.wt = NumericVector(nj);
    for (int r = 0; r < nj; ++r) {
      const int row = rows[r] - 1;
      g.y[r] = spec.y[row];
      g.offset[r] = spec.offset1.size() == 1 ? spec.offset1[0] : spec.offset1[row];
      g.wt[r] = spec.wt1.size() == 1 ? spec.wt1[0] : spec.wt1[row];
      for (int a = 0; a < l1; ++a) {
        g.x(r, a) = spec.D(row, a);
      }
    }
    g.P = spec.ridge ? as<NumericMatrix>(spec.ridge_P[j]) : spec.P_D;
    g.disp = spec.disp1.size() == 1 ? spec.disp1[0] : spec.disp1[j];
  }
  return out;
}

NumericVector groups_prior_mean_all(const TwoBlockIndexSpec& spec,
                                    const List& fixef, int n,
                                    std::vector<NumericMatrix>* mu_D) {
  const int l1 = spec.l1;
  const int J = spec.J;
  NumericVector out(static_cast<R_xlen_t>(l1) * J * n);
  out.attr("dim") = IntegerVector::create(l1, J, n);
  if (mu_D) mu_D->assign(n, NumericMatrix());
  for (int i = 0; i < n; ++i) {
    const NumericMatrix mi = groups_prior_mean(spec, fixef, i);
    const NumericMatrix m = spec.ridge ? groups_ridge_mean(spec, mi) : mi;
    const R_xlen_t base = static_cast<R_xlen_t>(l1) * J * i;
    for (int j = 0; j < J; ++j) {
      for (int a = 0; a < l1; ++a) {
        out[base + a + static_cast<R_xlen_t>(l1) * j] = m(a, j);
      }
    }
    if (mu_D) (*mu_D)[i] = mi;
  }
  return out;
}

GroupChainsDraw group_chains_draw(const TwoBlockIndexSpec& spec,
                                  const GroupData& gd,
                                  const NumericVector& mu_all, int j,
                                  const std::vector<int>& chains, bool check) {
  const int l1 = spec.l1;
  const int J = spec.J;
  const int nc = static_cast<int>(chains.size());
  NumericMatrix mu(l1, nc);
  for (int c = 0; c < nc; ++c) {
    const R_xlen_t base =
        static_cast<R_xlen_t>(l1) * (j + static_cast<R_xlen_t>(J) * chains[c]);
    for (int a = 0; a < l1; ++a) {
      mu(a, c) = mu_all[base + a];
    }
  }
  const RObject disp = wrap(gd.disp);
  const RObject family = wrap(spec.family);
  const RObject link = wrap(spec.link);
  const RObject gridtype = wrap(spec.gridtype);
  const RObject n_envopt = wrap(spec.n_envopt);
  const RObject check_r = wrap(check);
  const RObject off = wrap(false);
  List out;
  {
    glmbayes::rng::RSyncScope rng_sync;
    out = glmbayesCore_rNormalGLMGroupChains(
        gd.y, gd.x, gd.offset, gd.wt, disp, mu, gd.P, spec.f2, spec.f3,
        family, link, gridtype, n_envopt, off, check_r, off);
  }
  GroupChainsDraw d;
  d.coef = as<NumericMatrix>(out["coefficients"]);
  d.coef_mode = as<NumericMatrix>(out["coef.mode"]);
  d.iters = as<NumericVector>(out["iters"]);
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
