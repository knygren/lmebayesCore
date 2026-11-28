// Population block draw from a factor of the posterior precision built once
// per two_block_sweep call, and the comparison against population_draw
// (glmbayesCore rNormalReg, which refits the same posterior per draw).

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "rng_utils.h"
#include "simfuncs_twoblock.h"

#include <chrono>
#include <cmath>
#include <string>

using namespace Rcpp;

namespace glmbayes {
namespace sim {
namespace twoblock {

namespace {

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
      .count();
}

double max_rel_diff(const NumericVector& a, const NumericVector& ref) {
  if (a.size() != ref.size()) return R_PosInf;
  double m = 0.0;
  for (R_xlen_t i = 0; i < a.size(); ++i) {
    const double d = std::fabs(a[i] - ref[i]) / (1.0 + std::fabs(ref[i]));
    if (!(d <= m)) m = d;
  }
  return m;
}

}  // namespace

int population_draw_method() {
  SEXP opt = Rf_GetOption1(Rf_install("lmebayesCore.pop_draw"));
  if (Rf_isNull(opt)) return POP_CHOL;
  if (!Rf_isString(opt) || Rf_length(opt) != 1) {
    stop("options(lmebayesCore.pop_draw) must be \"reg\", \"chol\" or "
         "\"compare\".");
  }
  const std::string s = CHAR(STRING_ELT(opt, 0));
  if (s == "reg") return POP_REG;
  if (s == "chol") return POP_CHOL;
  if (s == "compare") return POP_COMPARE;
  stop("options(lmebayesCore.pop_draw) must be \"reg\", \"chol\" or "
       "\"compare\".");
  return POP_REG;
}

PopulationFactor population_factor(const TwoBlockIndexSpec& spec, int k) {
  const NumericMatrix& W = spec.W[k];
  const NumericMatrix& P = spec.P2[k];
  const NumericVector& mu = spec.mu2[k];
  const int n_w = W.nrow();
  const int q = W.ncol();
  if (P.nrow() != q || P.ncol() != q || mu.size() != q) {
    stop("population_factor: P2 and mu2 for component %d do not match W.",
         k + 1);
  }
  const double disp = spec.disp2[k];
  if (!(disp > 0.0)) {
    stop("population_factor: disp2 for component %d must be positive.", k + 1);
  }
  const arma::mat Wm(const_cast<double*>(W.begin()), n_w, q, false);
  const arma::mat Pm(const_cast<double*>(P.begin()), q, q, false);
  const arma::vec mum(const_cast<double*>(mu.begin()), q, false);

  PopulationFactor f;
  f.q = q;
  f.n_w = n_w;
  f.Wt = arma::trans(Wm) / disp;
  f.Pmu = Pm * mum;
  arma::mat R;
  if (!arma::chol(R, f.Wt * Wm + Pm)) {
    stop("population_factor: the posterior precision for component %d is "
         "not positive definite.", k + 1);
  }
  f.Rinv = arma::inv(arma::trimatu(R));
  return f;
}

PopulationDraw population_draw_chol(const PopulationFactor& f,
                                    const NumericVector& y_k) {
  if (y_k.size() != f.n_w) {
    stop("population_draw_chol: length(y_k) must equal nrow(W).");
  }
  const arma::vec y(const_cast<double*>(y_k.begin()), f.n_w, false);
  const arma::vec m = f.Rinv * (arma::trans(f.Rinv) * (f.Wt * y + f.Pmu));
  arma::vec z(f.q);
  {
    glmbayes::rng::RSyncScope rng_sync;
    for (int a = 0; a < f.q; ++a) z[a] = R::norm_rand();
  }
  const arma::vec b = m + f.Rinv * z;
  PopulationDraw d;
  d.coef = NumericVector(b.begin(), b.end());
  d.coef_mode = NumericVector(m.begin(), m.end());
  return d;
}

PopulationDraw population_draw_compare(PopulationCompare& pc,
                                       const TwoBlockIndexSpec& spec,
                                       const PopulationFactor& f,
                                       const NumericVector& y_k, int k) {
  RNGScope scope;
  Environment genv = Environment::global_env();
  PutRNGstate();
  if (!genv.exists(".Random.seed")) {
    stop("population_draw_compare: .Random.seed is not set.");
  }
  const IntegerVector saved = clone(as<IntegerVector>(genv.get(".Random.seed")));

  auto t0 = std::chrono::steady_clock::now();
  const PopulationDraw d_reg = population_draw(spec, y_k, k);
  pc.t_reg += seconds_since(t0);

  PutRNGstate();
  const IntegerVector after = clone(as<IntegerVector>(genv.get(".Random.seed")));
  genv.assign(".Random.seed", saved);
  GetRNGstate();

  t0 = std::chrono::steady_clock::now();
  const PopulationDraw d_chol = population_draw_chol(f, y_k);
  pc.t_chol += seconds_since(t0);

  genv.assign(".Random.seed", after);
  GetRNGstate();

  const double dm = max_rel_diff(d_chol.coef_mode, d_reg.coef_mode);
  const double dd = max_rel_diff(d_chol.coef, d_reg.coef);
  pc.calls += 1.0;
  if (!(dm <= 1e-10) || !(dd <= 1e-10)) pc.flagged += 1.0;
  if (!(dm <= pc.max_mean)) pc.max_mean = dm;
  if (!(dd <= pc.max_draw)) pc.max_draw = dd;
  return d_reg;
}

}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes
