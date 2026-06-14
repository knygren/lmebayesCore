// twoBlockGibbsStaged.cpp
// Staged two-block Gibbs driver: pilot replicate chains, Hotelling chi-squared
// test vs fixef_start, then main replicate chains.
// Reuses two_block_rNormal_reg_v2_cpp_export for both sampling stages.
// Pilot eigenvalue upper bounds and TV calibration run in the R wrapper (Phase 2b).

#include "RcppArmadillo.h"
#include "simfuncs.h"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace glmbayes {
namespace sim {

namespace {

using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::RObject;

// Stack fixef_draws (List of n x q_k matrices) into n x p_dim matrix.
NumericMatrix fixef_draws_to_matrix(const List& fixef_draws, int n_chains) {
  const int p_re = fixef_draws.size();
  int p_dim = 0;
  std::vector<int> qk(static_cast<size_t>(p_re));
  for (int k = 0; k < p_re; ++k) {
    NumericMatrix M = Rcpp::as<NumericMatrix>(fixef_draws[k]);
    if (M.nrow() != n_chains) {
      Rcpp::stop("fixef_draws[[%d]] has %d rows; expected %d.", k + 1, M.nrow(), n_chains);
    }
    qk[static_cast<size_t>(k)] = M.ncol();
    p_dim += M.ncol();
  }
  NumericMatrix X(n_chains, p_dim);
  int col0 = 0;
  for (int k = 0; k < p_re; ++k) {
    NumericMatrix M = Rcpp::as<NumericMatrix>(fixef_draws[k]);
    const int qk_k = qk[static_cast<size_t>(k)];
    for (int c = 0; c < qk_k; ++c) {
      for (int r = 0; r < n_chains; ++r) {
        X(r, col0 + c) = M(r, c);
      }
    }
    col0 += qk_k;
  }
  return X;
}

NumericVector fixef_start_to_vector(const List& fixef_start) {
  const int p_re = fixef_start.size();
  int p_dim = 0;
  for (int k = 0; k < p_re; ++k) {
    p_dim += Rcpp::as<NumericVector>(fixef_start[k]).size();
  }
  NumericVector out(p_dim);
  int j = 0;
  for (int k = 0; k < p_re; ++k) {
    NumericVector v = Rcpp::as<NumericVector>(fixef_start[k]);
    for (int i = 0; i < v.size(); ++i) {
      out[j++] = v[i];
    }
  }
  return out;
}

List colmeans_fixef_draws(const List& fixef_draws, const List& fixef_start) {
  const int p_re = fixef_draws.size();
  List out(p_re);
  for (int k = 0; k < p_re; ++k) {
    NumericMatrix M = Rcpp::as<NumericMatrix>(fixef_draws[k]);
    const int qk = M.ncol();
    NumericVector mu(qk);
    for (int c = 0; c < qk; ++c) {
      double s = 0.0;
      for (int r = 0; r < M.nrow(); ++r) {
        s += M(r, c);
      }
      mu[c] = s / static_cast<double>(M.nrow());
    }
    NumericVector start_k = Rcpp::as<NumericVector>(fixef_start[k]);
    RObject start_nm = start_k.attr("names");
    if (!start_nm.isNULL()) {
      mu.attr("names") = start_nm;
    }
    out[k] = mu;
  }
  RObject list_nm = fixef_start.attr("names");
  if (!list_nm.isNULL()) {
    out.attr("names") = list_nm;
  }
  return out;
}

List pilot_chisq_test(
    const List& pilot_fixef_draws,
    const List& fixef_mode,
    int n_pilot
) {
  NumericMatrix X = fixef_draws_to_matrix(pilot_fixef_draws, n_pilot);
  const int p_dim = X.ncol();
  if (p_dim < 1) {
    Rcpp::stop("pilot chi-squared test requires at least one hyper-parameter.");
  }
  NumericVector mode_vec = fixef_start_to_vector(fixef_mode);

  arma::mat Xa(const_cast<double*>(X.begin()), X.nrow(), X.ncol(), false);
  arma::rowvec mu = arma::mean(Xa, 0);
  arma::vec mode_a(mode_vec.begin(), mode_vec.size(), false);
  arma::vec d = mu.t() - mode_a;

  arma::mat Xc = Xa.each_row() - mu;
  arma::mat S;
  if (n_pilot > 1) {
    S = (Xc.t() * Xc) / static_cast<double>(n_pilot - 1);
  } else {
    S = arma::zeros<arma::mat>(p_dim, p_dim);
  }
  double ridge = 1e-8 * arma::trace(S) / static_cast<double>(p_dim);
  if (!std::isfinite(ridge) || ridge <= 0.0) {
    ridge = 1e-8;
  }
  arma::mat Sinv = arma::inv(S + ridge * arma::eye<arma::mat>(p_dim, p_dim));
  const double Q = static_cast<double>(
    n_pilot * arma::as_scalar(d.t() * Sinv * d)
  );
  const double p_val = R::pchisq(
    Q, static_cast<double>(p_dim), FALSE, FALSE
  );

  return List::create(
    Named("Q") = Q,
    Named("df") = p_dim,
    Named("p_value") = p_val,
    Named("n_pilot") = n_pilot
  );
}

} // anonymous namespace

List two_block_rNormal_reg_staged_cpp_export(
    int n_main,
    int m_convergence_main,
    int n_pilot,
    int m_convergence_pilot,
    const NumericVector& y,
    const NumericMatrix& x,
    const RObject& block,
    const List& x_hyper,
    const List& prior_list_block1,
    const RObject& dispersion_block1,
    const RObject& ddef_block1,
    const List& pfamily_list,
    const List& fixef_start,
    const Rcpp::CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar_main,
    bool progbar_pilot
) {
  if (n_main < 1) {
    Rcpp::stop("'n_main' must be at least 1.");
  }
  if (m_convergence_main < 1) {
    Rcpp::stop("'m_convergence_main' must be at least 1.");
  }
  if (n_pilot < 0) {
    Rcpp::stop("'n_pilot' must be non-negative.");
  }
  if (n_pilot > 0 && m_convergence_pilot < 1) {
    Rcpp::stop("'m_convergence_pilot' must be at least 1 when n_pilot > 0.");
  }

  const bool run_pilot = n_pilot > 0;
  List fixef_main_start = Rcpp::clone(fixef_start);
  int m_convergence_used = m_convergence_main;

  List pilot_out;
  List pilot_chisq;

  if (run_pilot) {
    pilot_out = two_block_rNormal_reg_v2_cpp_export(
      n_pilot, m_convergence_pilot, y, x, block, x_hyper,
      prior_list_block1, dispersion_block1, ddef_block1,
      pfamily_list, fixef_start, group_levels,
      family, link, f2, f3, f2_gauss, f3_gauss,
      offset, wt, Gridtype, n_envopt,
      use_parallel, use_opencl, verbose, progbar_pilot
    );

    List pilot_fixef = pilot_out["fixef_draws"];
    pilot_chisq = pilot_chisq_test(pilot_fixef, fixef_start, n_pilot);
    fixef_main_start = colmeans_fixef_draws(pilot_fixef, fixef_start);
  }

  List main_out = two_block_rNormal_reg_v2_cpp_export(
    n_main, m_convergence_used, y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    pfamily_list, fixef_main_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt, Gridtype, n_envopt,
    use_parallel, use_opencl, verbose, progbar_main
  );

  List res = Rcpp::clone(main_out);
  res["fixef_main_start"] = fixef_main_start;
  res["m_convergence_used"] = m_convergence_used;
  res["n_main"] = n_main;
  res["n_pilot"] = n_pilot;
  res["m_convergence_main"] = m_convergence_main;
  if (run_pilot) {
    res["pilot"] = pilot_out;
    res["pilot_chisq"] = pilot_chisq;
    res["m_convergence_pilot"] = m_convergence_pilot;
  }
  return res;
}

} // namespace sim
} // namespace glmbayes
