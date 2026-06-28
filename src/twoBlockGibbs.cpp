// twoBlockGibbs.cpp
// C++ port of the two_block_rNormal_reg() Gibbs loop (port-only).
//
// Mirrors the R loop in R/two_block_rNormal_reg.R verbatim:
//   for (i in 1:n) {
//     fixef <- fixef_start                          # "replicate" sampling
//     for (m in 1:m_convergence) {
//       mu_all  <- .two_block_mu_all(fixef, x_hyper, re_names, group_levels)
//       Block 1 <- block_rNormalReg() / block_rNormalGLM()   (n = 1)
//       Block 2 <- multi_rNormal_reg() list-x branch (rNormal_reg gaussian
//                  per RE component, n = 1)
//     }
//   }
// Block 1 goes through the existing block_rNormalReg_cpp_export /
// block_rNormalGLM_cpp_export (per-iteration prior normalization unchanged).
// Block 2 replicates rNormal_reg()'s gaussian-branch prep per component per
// iteration (Sigma -> P via inv_sympd == chol2inv(chol()), PD checks), then
// calls the rNormalReg core.  f2/f3 R closures and the optim call inside
// rNormalGLM are untouched.

#include "RcppArmadillo.h"
#include "simfuncs.h"
#include "progress_utils.h"

#include <iomanip>

#include <string>
#include <vector>

namespace glmbayes {
namespace sim {

namespace {

using Rcpp::CharacterVector;
using Rcpp::Function;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

// Mirror R's is.null(): a present-but-NULL element counts as absent.
inline bool has_non_null(const List& pl, const char* name) {
  return pl.containsElementNamed(name) && !Rf_isNull(pl[name]);
}

inline bool is_symmetric_mat(const NumericMatrix& M, double tol = 1e-8) {
  if (M.nrow() != M.ncol()) return false;
  for (int i = 0; i < M.nrow(); ++i) {
    for (int j = i + 1; j < M.ncol(); ++j) {
      if (std::fabs(M(i, j) - M(j, i)) > tol) return false;
    }
  }
  return true;
}

// R .check_symmetric_pd() / rNormal_reg PD check: eigen, tol 1e-6.
inline void check_pd(const NumericMatrix& M, const char* label) {
  arma::mat M2(const_cast<double*>(M.begin()), M.nrow(), M.ncol(), false);
  arma::vec ev = arma::eig_sym(M2);
  const double tol = 1e-6;
  if (ev.min() < -tol * std::fabs(ev.max())) {
    Rcpp::stop("'%s' is not positive definite.", label);
  }
}

// Port of .two_block_block1_prior_list(): list(mu, dispersion, ddef, P/Sigma).
// dispersion/ddef are forwarded as-is (possibly NULL); the block exports
// treat present-but-NULL as absent, matching R semantics.
List block1_prior_list(
    const NumericMatrix& mu_all,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1
) {
  List out = List::create(
    Rcpp::Named("mu") = mu_all,
    Rcpp::Named("dispersion") = dispersion_block1,
    Rcpp::Named("ddef") = ddef_block1
  );
  if (has_non_null(prior_list_block1, "P")) {
    out["P"] = prior_list_block1["P"];
  }
  if (has_non_null(prior_list_block1, "Sigma")) {
    out["Sigma"] = prior_list_block1["Sigma"];
  }
  if (!has_non_null(out, "P") && !has_non_null(out, "Sigma")) {
    Rcpp::stop("prior_list_block1 must contain 'P' or 'Sigma'.");
  }
  return out;
}

// Port of .validate_normal_prior_list() + rNormal_reg() gaussian-branch prior
// prep for one Block 2 component: mu, P (from chol2inv(chol(Sigma)) when P is
// absent; inv_sympd is the same LAPACK path), required dispersion.
struct Block2Prior {
  NumericVector mu;
  NumericMatrix P;
  double dispersion;
};

Block2Prior block2_prior_prep(const List& pl, int j1 /*1-based*/, int p) {
  if (!has_non_null(pl, "mu")) {
    Rcpp::stop("prior_list[[%d]] must contain 'mu'.", j1);
  }
  if (!has_non_null(pl, "Sigma") && !has_non_null(pl, "P")) {
    Rcpp::stop("prior_list[[%d]] must contain 'Sigma' or 'P'.", j1);
  }

  NumericVector mu = Rcpp::as<NumericVector>(pl["mu"]);
  if (mu.size() != p) {
    Rcpp::stop("prior_list[[%d]]$mu must have length ncol(x) = %d.", j1, p);
  }

  NumericMatrix P;
  if (has_non_null(pl, "Sigma")) {
    NumericMatrix S = Rcpp::as<NumericMatrix>(pl["Sigma"]);
    if (S.nrow() != p || S.ncol() != p) {
      Rcpp::stop("prior_list[[%d]]$Sigma must be %d x %d.", j1, p, p);
    }
    if (!is_symmetric_mat(S)) {
      Rcpp::stop("prior_list[[%d]]$Sigma must be symmetric.", j1);
    }
    check_pd(S, "Sigma");
    if (!has_non_null(pl, "P")) {
      // rNormal_reg(): P <- 0.5 * (chol2inv(chol(Sigma)) + t(...))
      arma::mat S2(const_cast<double*>(S.begin()), p, p, false);
      arma::mat Pinv = arma::inv_sympd(S2);
      P = NumericMatrix(Rcpp::wrap(0.5 * (Pinv + Pinv.t())));
    }
  }
  if (has_non_null(pl, "P")) {
    P = Rcpp::as<NumericMatrix>(pl["P"]);
    if (P.nrow() != p || P.ncol() != p) {
      Rcpp::stop("prior_list[[%d]]$P must be %d x %d.", j1, p, p);
    }
  }
  if (!is_symmetric_mat(P)) {
    Rcpp::stop("matrix P must be symmetric");
  }
  check_pd(P, "P");

  // rNormal_reg(): gaussian requires an explicit dispersion (ddef rules).
  bool ddef;
  if (pl.containsElementNamed("ddef")) {
    SEXP dd = pl["ddef"];
    ddef = Rf_isLogical(dd) && Rf_length(dd) >= 1 &&
           LOGICAL(dd)[0] == TRUE;
  } else {
    ddef = !has_non_null(pl, "dispersion");
  }
  if (ddef || !has_non_null(pl, "dispersion")) {
    Rcpp::stop(
      "For gaussian() models, dNormal() requires an explicit dispersion "
      "(prior_list[[%d]]). Omitted or NULL dispersion is not allowed", j1
    );
  }
  double dispersion = Rcpp::as<NumericVector>(pl["dispersion"])[0];

  Block2Prior out;
  out.mu = mu;
  out.P = P;
  out.dispersion = dispersion;
  return out;
}

} // anonymous namespace

List two_block_rNormal_reg_cpp_export(
    int n,
    int m_convergence,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    const List& x_hyper,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const List& prior_list_block2,
    const List& fixef_start,
    const CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Function& f2,
    const Function& f3,
    const Function& f2_gauss,
    const Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (m_convergence < 1) {
    Rcpp::stop("'m_convergence' must be at least 1.");
  }
  const int p_re = x.ncol();
  const int J = group_levels.size();
  if (x_hyper.size() != p_re) {
    Rcpp::stop("length(x_hyper) must equal ncol(x) = %d.", p_re);
  }
  if (prior_list_block2.size() != p_re) {
    Rcpp::stop("length(prior_list_block2) must equal ncol(x) = %d.", p_re);
  }
  if (fixef_start.size() != p_re) {
    Rcpp::stop("length(fixef_start) must equal ncol(x) = %d.", p_re);
  }

  const bool is_gaussian = (family == "gaussian");

  MuAllBuilder mu_builder(x_hyper, group_levels);

  // fixef state (component order fixed by the R wrapper)
  std::vector<NumericVector> fixef_start_v(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_start_v[j] = Rcpp::as<NumericVector>(fixef_start[j]);
  }
  std::vector<NumericVector> fixef = fixef_start_v;

  // Storage: fixef draws per component (n x q_k); b draws as J x p_re x n.
  List fixef_draws(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_draws[j] = NumericMatrix(n, fixef_start_v[j].size());
  }
  NumericVector b_arr(Rcpp::Dimension(J, p_re, n));

  NumericMatrix mu_all;
  NumericMatrix b_i;
  CharacterVector group_ids;
  bool have_ids = false;

  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar) {
      glmbayes::progress::progress_bar(
        static_cast<double>(i + 1), static_cast<double>(n)
      );
    }

    fixef = fixef_start_v;

    for (int m = 0; m < m_convergence; ++m) {

      mu_all = mu_builder.build(fixef);
      List pl1 = block1_prior_list(
        mu_all, prior_list_block1, dispersion_block1, ddef_block1
      );

      List block_i;
      if (is_gaussian) {
        // R: block_rNormalReg(n = 1, ...) with default Gridtype = 2L.
        block_i = block_rNormalReg_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3, 2
        );
      } else {
        block_i = block_rNormalGLM_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3,
          family, link, Gridtype, n_envopt,
          use_parallel, use_opencl, verbose
        );
      }

      b_i = Rcpp::as<NumericMatrix>(block_i["coefficients"]);
      if (b_i.nrow() != J || b_i.ncol() != p_re) {
        Rcpp::stop(
          "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
          b_i.nrow(), b_i.ncol(), J, p_re
        );
      }
      if (!have_ids) {
        List bi_info = block_i["block_info"];
        group_ids = Rcpp::as<CharacterVector>(bi_info["ids"]);
        have_ids = true;
      }

      // Block 2: multi_rNormal_reg() list-x branch, n = 1 per component.
      for (int j = 0; j < p_re; ++j) {
        const NumericMatrix& X_j = mu_builder.X[j];
        if (X_j.nrow() != b_i.nrow()) {
          Rcpp::stop(
            "nrow(x_hyper[[%d]]) (%d) must equal nrow(y) (%d).",
            j + 1, X_j.nrow(), b_i.nrow()
          );
        }
        Block2Prior pr = block2_prior_prep(
          List(prior_list_block2[j]), j + 1, X_j.ncol()
        );
        NumericVector y_j = b_i(Rcpp::_, j);
        NumericVector offset_j(X_j.nrow(), 0.0);
        NumericVector wt_j(X_j.nrow(), 1.0);

        List out_j = rNormalReg(
          1, y_j, X_j, pr.mu, pr.P, offset_j, wt_j,
          pr.dispersion, f2_gauss, f3_gauss, pr.mu,
          "gaussian", "identity", 2
        );
        NumericMatrix coef_j = Rcpp::as<NumericMatrix>(out_j["coefficients"]);
        fixef[j] = NumericVector(coef_j(0, Rcpp::_));
      }
    }

    for (int j = 0; j < p_re; ++j) {
      NumericMatrix fd = fixef_draws[j];
      const NumericVector& fj = fixef[j];
      if (fj.size() != fd.ncol()) {
        Rcpp::stop(
          "Block 2 draw for component %d has length %d; expected %d.",
          j + 1, fj.size(), fd.ncol()
        );
      }
      fd(i, Rcpp::_) = fj;
    }
    for (int j = 0; j < p_re; ++j) {
      for (int g = 0; g < J; ++g) {
        b_arr[g + J * (j + p_re * i)] = b_i(g, j);
      }
    }
  }

  if (progbar) {
    glmbayes::progress::progress_bar_finish();
  }

  List fixef_last(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_last[j] = fixef[j];
  }

  return List::create(
    Rcpp::Named("fixef_draws") = fixef_draws,
    Rcpp::Named("b_draws") = b_arr,
    Rcpp::Named("fixef_last") = fixef_last,
    Rcpp::Named("b_last") = b_i,
    Rcpp::Named("mu_all_last") = mu_all,
    Rcpp::Named("group_ids") = group_ids
  );
}

// ===========================================================================
// v2: pfamily-based Block 2 contract (development track).
//
// Differences vs the v1 driver above:
//   * Block 2 priors arrive as glmbayesCore pfamily objects (dNormal or
//     dIndependent_Normal_Gamma); the type string in pf$pfamily selects the
//     Block 2 update (native analogue of rglmb()'s simfun dispatch).
//   * Block 2 priors are parsed once before the loop (no RNG is consumed in
//     the prep, so the dNormal draw stream is identical to v1).
//   * A tau2 working vector tracks the Block 2 dispersion per RE component;
//     for ING components it will be updated by the joint (gamma_k, tau2_k)
//     draw and fed back into the Block 1 prior precision (next milestone).
//     The values at each stored draw are returned as dispersion_fixef_draws.
// ===========================================================================

namespace {

struct Block2PriorV2 {
  NumericVector mu;
  NumericMatrix P;
  double dispersion;       // dNormal: fixed tau2_k; ING: initial tau2_k
  bool is_ing;
  double shape;
  double rate;
  double max_disp_perc;
  SEXP disp_lower;         // NumericVector or R_NilValue (rglmb semantics)
  SEXP disp_upper;
};

// Shared mu / Sigma-or-P parsing (same checks and LAPACK path as
// block2_prior_prep above; duplicated so the v1 code stays byte-identical
// until v2 is promoted).
void block2_mu_P_prep_v2(
    const List& pl, int j1 /*1-based*/, int p,
    NumericVector& mu_out, NumericMatrix& P_out
) {
  if (!has_non_null(pl, "mu")) {
    Rcpp::stop("pfamily_list[[%d]]$prior_list must contain 'mu'.", j1);
  }
  if (!has_non_null(pl, "Sigma") && !has_non_null(pl, "P")) {
    Rcpp::stop(
      "pfamily_list[[%d]]$prior_list must contain 'Sigma' or 'P'.", j1
    );
  }

  NumericVector mu = Rcpp::as<NumericVector>(pl["mu"]);
  if (mu.size() != p) {
    Rcpp::stop(
      "pfamily_list[[%d]]$prior_list$mu must have length ncol(x) = %d.",
      j1, p
    );
  }

  NumericMatrix P;
  if (has_non_null(pl, "Sigma")) {
    NumericMatrix S = Rcpp::as<NumericMatrix>(pl["Sigma"]);
    if (S.nrow() != p || S.ncol() != p) {
      Rcpp::stop(
        "pfamily_list[[%d]]$prior_list$Sigma must be %d x %d.", j1, p, p
      );
    }
    if (!is_symmetric_mat(S)) {
      Rcpp::stop(
        "pfamily_list[[%d]]$prior_list$Sigma must be symmetric.", j1
      );
    }
    check_pd(S, "Sigma");
    if (!has_non_null(pl, "P")) {
      arma::mat S2(const_cast<double*>(S.begin()), p, p, false);
      arma::mat Pinv = arma::inv_sympd(S2);
      P = NumericMatrix(Rcpp::wrap(0.5 * (Pinv + Pinv.t())));
    }
  }
  if (has_non_null(pl, "P")) {
    P = Rcpp::as<NumericMatrix>(pl["P"]);
    if (P.nrow() != p || P.ncol() != p) {
      Rcpp::stop(
        "pfamily_list[[%d]]$prior_list$P must be %d x %d.", j1, p, p
      );
    }
  }
  if (!is_symmetric_mat(P)) {
    Rcpp::stop("matrix P must be symmetric");
  }
  check_pd(P, "P");

  mu_out = mu;
  P_out = P;
}

Block2PriorV2 block2_prior_prep_v2(const List& pf, int j1 /*1-based*/, int p) {
  if (!has_non_null(pf, "pfamily") || !has_non_null(pf, "prior_list")) {
    Rcpp::stop(
      "pfamily_list[[%d]] must be a pfamily object (e.g. dNormal() or "
      "dIndependent_Normal_Gamma()) with 'pfamily' and 'prior_list' fields.",
      j1
    );
  }
  const std::string ptype = Rcpp::as<std::string>(pf["pfamily"]);
  List pl(pf["prior_list"]);

  Block2PriorV2 out;
  out.is_ing = false;
  out.shape = NA_REAL;
  out.rate = NA_REAL;
  out.max_disp_perc = 0.99;
  out.disp_lower = R_NilValue;
  out.disp_upper = R_NilValue;

  block2_mu_P_prep_v2(pl, j1, p, out.mu, out.P);

  if (ptype == "dNormal") {
    // Same dispersion requirement as the v1 gaussian Block 2 path.
    bool ddef;
    if (pl.containsElementNamed("ddef")) {
      SEXP dd = pl["ddef"];
      ddef = Rf_isLogical(dd) && Rf_length(dd) >= 1 && LOGICAL(dd)[0] == TRUE;
    } else {
      ddef = !has_non_null(pl, "dispersion");
    }
    if (ddef || !has_non_null(pl, "dispersion")) {
      Rcpp::stop(
        "For gaussian() Block 2, dNormal() requires an explicit dispersion "
        "(pfamily_list[[%d]]).", j1
      );
    }
    out.dispersion = Rcpp::as<NumericVector>(pl["dispersion"])[0];
  } else if (ptype == "dIndependent_Normal_Gamma") {
    out.is_ing = true;
    if (!has_non_null(pl, "shape") || !has_non_null(pl, "rate")) {
      Rcpp::stop(
        "pfamily_list[[%d]]: dIndependent_Normal_Gamma requires 'shape' "
        "and 'rate'.", j1
      );
    }
    out.shape = Rcpp::as<NumericVector>(pl["shape"])[0];
    out.rate = Rcpp::as<NumericVector>(pl["rate"])[0];
    if (!(out.shape > 0.0) || !(out.rate > 0.0)) {
      Rcpp::stop(
        "pfamily_list[[%d]]: 'shape' and 'rate' must be positive.", j1
      );
    }
    if (has_non_null(pl, "max_disp_perc")) {
      out.max_disp_perc = Rcpp::as<NumericVector>(pl["max_disp_perc"])[0];
    }
    if (has_non_null(pl, "disp_lower")) {
      out.disp_lower = pl["disp_lower"];
    }
    if (has_non_null(pl, "disp_upper")) {
      out.disp_upper = pl["disp_upper"];
    }
    // Initial tau2_k for the first Block 1 sweep of each replicate chain:
    // the conservative disp_lower plug-in (consistent with the lmebayes
    // calibration); refined by the joint Block 2 draw within the sweep.
    out.dispersion = has_non_null(pl, "disp_lower")
      ? Rcpp::as<NumericVector>(pl["disp_lower"])[0]
      : NA_REAL;
  } else {
    Rcpp::stop(
      "pfamily_list[[%d]]: unsupported pfamily '%s' (allowed: dNormal, "
      "dIndependent_Normal_Gamma).", j1, ptype.c_str()
    );
  }

  return out;
}

void v3_reseed_r(int seed_val) {
  Rcpp::Environment base("package:base");
  Rcpp::Function set_seed = base["set.seed"];
  set_seed(seed_val);
}

void two_block_driver_banner_v3(
    int n,
    int m_convergence,
    int seed_offset,
    bool have_seed
) {
  Rcpp::Rcout << "--- glmbayesCore two_block_rNormal_reg_v3 (C++ driver)"
              << ": n=" << n
              << " m_convergence=" << m_convergence
              << " seed_offset=" << seed_offset
              << (have_seed ? " seeded" : " unseeded")
              << " ---\n";
}

void two_block_driver_banner_v4(
    int n,
    int m_convergence,
    int seed_offset,
    bool have_seed,
    int p_dim
) {
  Rcpp::Rcout << "--- glmbayesCore two_block_rNormal_reg_v4 (C++ driver)"
              << ": n=" << n
              << " m_convergence=" << m_convergence
              << " seed_offset=" << seed_offset
              << (have_seed ? " seeded" : " unseeded")
              << "; fixef_temp " << p_dim << "-col load/store ON"
              << ", tau2_temp load/store ON ---\n";
}


void two_block_driver_banner_v5(
    int n,
    int m_convergence,
    int seed_offset,
    bool have_seed,
    int p_dim
) {
  Rcpp::Rcout << "--- glmbayesCore two_block_rNormal_reg_v5 (C++ driver)"
              << ": n=" << n
              << " m_convergence=" << m_convergence
              << " seed_offset=" << seed_offset
              << (have_seed ? " seeded" : " unseeded")
              << "; fixef_temp " << p_dim << "-col load/store ON"
              << ", tau2_temp load/store ON"
              << "; loop order: sweep-outer (m then i) ---\n";
}

int fixef_col_offset(
    const std::vector<int>& q_k,
    int j
) {
  int col0 = 0;
  for (int k = 0; k < j; ++k) {
    col0 += q_k[static_cast<size_t>(k)];
  }
  return col0;
}

void two_block_assign_fixef_component(
    std::vector<Rcpp::NumericVector>& fixef,
    int j,
    const Rcpp::NumericVector& gamma_new,
    const Rcpp::NumericVector& template_v
) {
  if (gamma_new.size() != template_v.size()) {
    Rcpp::stop(
      "Block 2 draw for component %d has length %d; expected %d.",
      j + 1, gamma_new.size(), template_v.size()
    );
  }
  fixef[static_cast<size_t>(j)] = Rcpp::clone(gamma_new);
  if (template_v.hasAttribute("names")) {
    fixef[static_cast<size_t>(j)].attr("names") = template_v.attr("names");
  }
}

void two_block_pack_fixef_row(
    Rcpp::NumericMatrix& fixef_temp,
    int i,
    int p_re,
    const std::vector<int>& q_k,
    const std::vector<Rcpp::NumericVector>& fixef_vec
) {
  for (int j = 0; j < p_re; ++j) {
    const int col0 = fixef_col_offset(q_k, j);
    const Rcpp::NumericVector& fj = fixef_vec[static_cast<size_t>(j)];
    if (fj.size() != q_k[static_cast<size_t>(j)]) {
      Rcpp::stop(
        "fixef[[%d]] has length %d; expected %d for fixef_temp pack.",
        j + 1, fj.size(), q_k[static_cast<size_t>(j)]
      );
    }
    for (int c = 0; c < q_k[static_cast<size_t>(j)]; ++c) {
      fixef_temp(i, col0 + c) = fj[c];
    }
  }
}

void two_block_unpack_fixef_row(
    const Rcpp::NumericMatrix& fixef_temp,
    int i,
    int p_re,
    const std::vector<int>& q_k,
    std::vector<Rcpp::NumericVector>& fixef_vec,
    const std::vector<Rcpp::NumericVector>& templates
) {
  for (int j = 0; j < p_re; ++j) {
    const int col0 = fixef_col_offset(q_k, j);
    fixef_vec[static_cast<size_t>(j)] =
      Rcpp::clone(templates[static_cast<size_t>(j)]);
    for (int c = 0; c < q_k[static_cast<size_t>(j)]; ++c) {
      fixef_vec[static_cast<size_t>(j)][c] = fixef_temp(i, col0 + c);
    }
  }
}

void two_block_print_fixef_line(
    const char* prefix,
    const std::vector<Rcpp::NumericVector>& fixef_vec,
    int p_re,
    const Rcpp::CharacterVector& re_names
) {
  Rcpp::Rcout << prefix;
  Rcpp::Rcout << std::fixed << std::setprecision(4);
  for (int j = 0; j < p_re; ++j) {
    const Rcpp::NumericVector& fj = fixef_vec[static_cast<size_t>(j)];
    Rcpp::CharacterVector pnames = fj.attr("names");
    const std::string re_lab =
      (re_names.size() == p_re && !Rcpp::CharacterVector::is_na(re_names[j]))
        ? Rcpp::as<std::string>(re_names[j])
        : ("RE" + std::to_string(j + 1));
    for (int c = 0; c < fj.size(); ++c) {
      const std::string p_lab =
        (pnames.size() == fj.size() && !Rcpp::CharacterVector::is_na(pnames[c]))
          ? Rcpp::as<std::string>(pnames[c])
          : ("p" + std::to_string(c + 1));
      Rcpp::Rcout << "  " << re_lab << "::" << p_lab << "=" << fj[c];
    }
  }
  Rcpp::Rcout << "\n";
}

std::vector<Rcpp::NumericVector> two_block_mean_fixef_components(
    const Rcpp::NumericMatrix& fixef_temp,
    int n,
    int p_re,
    const std::vector<int>& q_k,
    const std::vector<Rcpp::NumericVector>& templates
) {
  const double inv_n = 1.0 / static_cast<double>(n);
  std::vector<Rcpp::NumericVector> avg(static_cast<size_t>(p_re));
  for (int j = 0; j < p_re; ++j) {
    avg[static_cast<size_t>(j)] =
      Rcpp::NumericVector(q_k[static_cast<size_t>(j)], 0.0);
    if (templates[static_cast<size_t>(j)].hasAttribute("names")) {
      avg[static_cast<size_t>(j)].attr("names") =
        templates[static_cast<size_t>(j)].attr("names");
    }
    const int col0 = fixef_col_offset(q_k, j);
    for (int i = 0; i < n; ++i) {
      for (int c = 0; c < q_k[static_cast<size_t>(j)]; ++c) {
        avg[static_cast<size_t>(j)][c] += fixef_temp(i, col0 + c);
      }
    }
    for (int c = 0; c < q_k[static_cast<size_t>(j)]; ++c) {
      avg[static_cast<size_t>(j)][c] *= inv_n;
    }
  }
  return avg;
}

void two_block_print_sweep_chain_means(
    const std::string& stage_label,
    int sweep1,
    int n,
    int J,
    int p_re,
    const std::vector<int>& q_k,
    const Rcpp::NumericMatrix& fixef_temp,
    const Rcpp::NumericVector& b_work,
    const std::vector<Rcpp::NumericVector>& templates,
    const Rcpp::CharacterVector& re_names
) {
  if (!stage_label.empty()) {
    Rcpp::Rcout << "[" << stage_label << "] ";
  }
  Rcpp::Rcout << "sweep " << sweep1 << " chain means (C++): n=" << n << "\n";

  const std::vector<Rcpp::NumericVector> fe_mean =
    two_block_mean_fixef_components(fixef_temp, n, p_re, q_k, templates);
  two_block_print_fixef_line("  fixef:", fe_mean, p_re, re_names);

  Rcpp::Rcout << "  ranef:";
  Rcpp::Rcout << std::fixed << std::setprecision(4);
  const double inv_nJ = 1.0 / static_cast<double>(n) / static_cast<double>(J);
  for (int j = 0; j < p_re; ++j) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      for (int g = 0; g < J; ++g) {
        sum += b_work[g + J * (j + p_re * i)];
      }
    }
    const std::string re_lab =
      (re_names.size() == p_re && !Rcpp::CharacterVector::is_na(re_names[j]))
        ? Rcpp::as<std::string>(re_names[j])
        : ("RE" + std::to_string(j + 1));
    Rcpp::Rcout << "  " << re_lab << "=" << sum * inv_nJ;
  }
  Rcpp::Rcout << "\n\n";
}

std::vector<Rcpp::NumericVector> two_block_sd_fixef_components(
    const Rcpp::NumericMatrix& fixef_temp,
    int n,
    int p_re,
    const std::vector<int>& q_k,
    const std::vector<Rcpp::NumericVector>& templates
) {
  std::vector<Rcpp::NumericVector> sd(static_cast<size_t>(p_re));
  for (int j = 0; j < p_re; ++j) {
    sd[static_cast<size_t>(j)] =
      Rcpp::NumericVector(q_k[static_cast<size_t>(j)], NA_REAL);
    if (templates[static_cast<size_t>(j)].hasAttribute("names")) {
      sd[static_cast<size_t>(j)].attr("names") =
        templates[static_cast<size_t>(j)].attr("names");
    }
    const int col0 = fixef_col_offset(q_k, j);
    const int qj = q_k[static_cast<size_t>(j)];
    if (n <= 1) {
      continue;
    }
    const double inv_nm1 = 1.0 / static_cast<double>(n - 1);
    for (int c = 0; c < qj; ++c) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
        sum += fixef_temp(i, col0 + c);
      }
      const double mean = sum / static_cast<double>(n);
      double ss = 0.0;
      for (int i = 0; i < n; ++i) {
        const double d = fixef_temp(i, col0 + c) - mean;
        ss += d * d;
      }
      sd[static_cast<size_t>(j)][c] = std::sqrt(ss * inv_nm1);
    }
  }
  return sd;
}

Rcpp::List two_block_snapshot_fixef_stats_cpp(
    const Rcpp::NumericMatrix& fixef_temp,
    int n,
    int p_re,
    const std::vector<int>& q_k,
    const std::vector<Rcpp::NumericVector>& templates,
    const Rcpp::CharacterVector& re_names
) {
  const std::vector<Rcpp::NumericVector> means =
    two_block_mean_fixef_components(fixef_temp, n, p_re, q_k, templates);
  const std::vector<Rcpp::NumericVector> sds =
    two_block_sd_fixef_components(fixef_temp, n, p_re, q_k, templates);
  Rcpp::List snap(p_re);
  for (int j = 0; j < p_re; ++j) {
    snap[j] = Rcpp::List::create(
      Rcpp::Named("mean") = means[static_cast<size_t>(j)],
      Rcpp::Named("sd") = sds[static_cast<size_t>(j)]
    );
  }
  snap.attr("names") = re_names;
  return snap;
}

void two_block_accumulate_fixef_sum(
    Rcpp::NumericVector& fixef_sum,
    int p_re,
    const std::vector<int>& q_k,
    const std::vector<Rcpp::NumericVector>& fixef_vec
) {
  for (int j = 0; j < p_re; ++j) {
    const int col0 = fixef_col_offset(q_k, j);
    const Rcpp::NumericVector& fj = fixef_vec[static_cast<size_t>(j)];
    for (int c = 0; c < q_k[static_cast<size_t>(j)]; ++c) {
      fixef_sum[col0 + c] += fj[c];
    }
  }
}

void two_block_print_b_coefficients(
    const char* prefix,
    const Rcpp::NumericMatrix& b_vals,
    int J,
    int p_re,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& re_names
) {
  Rcpp::Rcout << prefix;
  Rcpp::Rcout << std::fixed << std::setprecision(4);
  for (int g = 0; g < J; ++g) {
    const std::string grp_lab =
      (group_levels.size() == J && !Rcpp::CharacterVector::is_na(group_levels[g]))
        ? Rcpp::as<std::string>(group_levels[g])
        : ("g" + std::to_string(g + 1));
    for (int j = 0; j < p_re; ++j) {
      const std::string re_lab =
        (re_names.size() == p_re && !Rcpp::CharacterVector::is_na(re_names[j]))
          ? Rcpp::as<std::string>(re_names[j])
          : ("RE" + std::to_string(j + 1));
      Rcpp::Rcout << "  " << grp_lab << "::" << re_lab << "=" << b_vals(g, j);
    }
  }
  Rcpp::Rcout << "\n";
}

void two_block_print_b_delta(
    const char* prefix,
    const Rcpp::NumericMatrix& b_avg,
    const Rcpp::NumericMatrix& b_mode,
    int J,
    int p_re,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& re_names
) {
  Rcpp::Rcout << prefix;
  Rcpp::Rcout << std::fixed << std::setprecision(4);
  for (int g = 0; g < J; ++g) {
    const std::string grp_lab =
      (group_levels.size() == J && !Rcpp::CharacterVector::is_na(group_levels[g]))
        ? Rcpp::as<std::string>(group_levels[g])
        : ("g" + std::to_string(g + 1));
    for (int j = 0; j < p_re; ++j) {
      const std::string re_lab =
        (re_names.size() == p_re && !Rcpp::CharacterVector::is_na(re_names[j]))
          ? Rcpp::as<std::string>(re_names[j])
          : ("RE" + std::to_string(j + 1));
      Rcpp::Rcout << "  " << grp_lab << "::" << re_lab << "="
                  << (b_avg(g, j) - b_mode(g, j));
    }
  }
  Rcpp::Rcout << "\n";
}

void two_block_print_sweep_diagnostic(
    const std::string& stage_label,
    int m,
    int m_convergence,
    int n,
    int p_re,
    int p_dim,
    int J,
    const Rcpp::NumericVector& fixef_sum,
    const Rcpp::NumericMatrix& b_sum,
    const std::vector<Rcpp::NumericVector>& fixef_mode_v,
    const Rcpp::NumericMatrix& b_mode,
    bool have_b_mode,
    const Rcpp::List& fixef_mode_list,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels
) {
  const double inv_n = 1.0 / static_cast<double>(n);

  Rcpp::Rcout << "[" << stage_label << " sweep " << (m + 1)
              << " / " << m_convergence << " AFTER, n=" << n << "]\n";

  std::vector<Rcpp::NumericVector> fixef_avg(static_cast<size_t>(p_re));
  int col0 = 0;
  for (int j = 0; j < p_re; ++j) {
    const Rcpp::NumericVector& fj_mode = fixef_mode_v[static_cast<size_t>(j)];
    fixef_avg[static_cast<size_t>(j)] = Rcpp::clone(fj_mode);
    for (int c = 0; c < fj_mode.size(); ++c) {
      fixef_avg[static_cast<size_t>(j)][c] = fixef_sum[col0 + c] * inv_n;
      col0 += 1;
    }
  }

  two_block_print_fixef_line(
    "  fixef avg:", fixef_avg, p_re, re_names
  );
  two_block_print_fixef_line(
    "  fixef mode:", fixef_mode_v, p_re, re_names
  );

  Rcpp::Rcout << std::fixed << std::setprecision(4);
  Rcpp::Rcout << "  fixef delta (avg - mode):";
  col0 = 0;
  for (int j = 0; j < p_re; ++j) {
    const Rcpp::NumericVector& fj_mode = fixef_mode_v[static_cast<size_t>(j)];
    for (int c = 0; c < fj_mode.size(); ++c) {
      const double avg = fixef_sum[col0 + c] * inv_n;
      Rcpp::Rcout << "  " << (avg - fj_mode[c]);
      col0 += 1;
    }
  }
  Rcpp::Rcout << "\n";

  Rcpp::NumericMatrix b_avg(J, p_re);
  for (int g = 0; g < J; ++g) {
    for (int j = 0; j < p_re; ++j) {
      b_avg(g, j) = b_sum(g, j) * inv_n;
    }
  }
  two_block_print_b_coefficients(
    "  b avg (mean over chains, each group x RE):",
    b_avg, J, p_re, group_levels, re_names
  );

  if (have_b_mode && b_mode.nrow() == J && b_mode.ncol() == p_re) {
    two_block_print_b_coefficients(
      "  b mode (ICM reference, each group x RE):",
      b_mode, J, p_re, group_levels, re_names
    );
    two_block_print_b_delta(
      "  b delta (avg - mode):",
      b_avg, b_mode, J, p_re, group_levels, re_names
    );
  }

  (void)fixef_mode_list;
  (void)p_dim;
}

// Prefix for sweep-outer progress bars (matches R .two_block_progbar_prefix).
std::string two_block_progbar_prefix(
    const std::string& stage_label,
    int sweep,
    int inner_sweeps,
    const char* phase
) {
  const char* phase_label =
    (std::string(phase) == "Block1") ? "RE" : "fixef";
  if (!stage_label.empty()) {
    return "[" + stage_label + "] sweep " + std::to_string(sweep) + "/" +
           std::to_string(inner_sweeps) + " " + phase_label + ": ";
  }
  return "sweep " + std::to_string(sweep) + "/" +
         std::to_string(inner_sweeps) + " " + phase_label + ": ";
}

std::string two_block_sweep_only_prefix(
    const std::string& stage_label,
    int sweep,
    int inner_sweeps
) {
  if (!stage_label.empty()) {
    return "[" + stage_label + "] sweep " + std::to_string(sweep) + "/" +
           std::to_string(inner_sweeps) + ": ";
  }
  return "sweep " + std::to_string(sweep) + "/" +
         std::to_string(inner_sweeps) + ": ";
}

// R batch Block 2 uses rglmb()$coef.mode (conditional mean), not coefficients.
NumericVector two_block_coef_mode_from_rNormalReg(const Rcpp::List& out_j) {
  SEXP cm = out_j["coef.mode"];
  if (Rf_isNull(cm)) {
    Rcpp::stop("rNormalReg result missing coef.mode.");
  }
  if (Rcpp::is<Rcpp::NumericMatrix>(cm)) {
    Rcpp::NumericMatrix M(cm);
    if (M.nrow() == 1) {
      return Rcpp::NumericVector(M(0, Rcpp::_));
    }
    if (M.ncol() == 1) {
      Rcpp::NumericVector v(M.nrow());
      for (int r = 0; r < M.nrow(); ++r) {
        v[r] = M(r, 0);
      }
      return v;
    }
  }
  return Rcpp::as<Rcpp::NumericVector>(cm);
}

namespace {

bool charvec_contains_name(
    const Rcpp::CharacterVector& names,
    const std::string& target
) {
  for (int i = 0; i < names.size(); ++i) {
    if (Rcpp::CharacterVector::is_na(names[i])) {
      continue;
    }
    if (Rcpp::as<std::string>(names[i]) == target) {
      return true;
    }
  }
  return false;
}

double charvec_named_value(
    const Rcpp::NumericVector& values,
    const Rcpp::CharacterVector& names,
    const std::string& target
) {
  for (int i = 0; i < names.size(); ++i) {
    if (Rcpp::CharacterVector::is_na(names[i])) {
      continue;
    }
    if (Rcpp::as<std::string>(names[i]) == target) {
      return values[i];
    }
  }
  Rcpp::stop("Group level \"%s\" missing from b.", target.c_str());
  return NA_REAL;
}

Rcpp::CharacterVector matrix_rownames(const Rcpp::NumericMatrix& X_k) {
  if (!X_k.hasAttribute("dimnames")) {
    return Rcpp::CharacterVector();
  }
  Rcpp::List dn(X_k.attr("dimnames"));
  if (Rf_isNull(dn[0])) {
    return Rcpp::CharacterVector();
  }
  return Rcpp::as<Rcpp::CharacterVector>(dn[0]);
}

void two_block_update_fixef_row(
    Rcpp::NumericVector& row,
    const Rcpp::NumericVector& coef_k
) {
  if (coef_k.hasAttribute("names") && !Rf_isNull(coef_k.attr("names"))) {
    Rcpp::CharacterVector row_names = row.attr("names");
    Rcpp::CharacterVector coef_names = coef_k.attr("names");
    for (int i = 0; i < coef_names.size(); ++i) {
      const std::string nm = Rcpp::as<std::string>(coef_names[i]);
      for (int c = 0; c < row_names.size(); ++c) {
        if (Rcpp::CharacterVector::is_na(row_names[c])) {
          continue;
        }
        if (Rcpp::as<std::string>(row_names[c]) == nm) {
          row[c] = coef_k[i];
          break;
        }
      }
    }
    return;
  }
  const int ncopy = std::min(row.size(), coef_k.size());
  for (int c = 0; c < ncopy; ++c) {
    row[c] = coef_k[c];
  }
}

} // namespace

} // anonymous namespace (Block~2 exports below need glmbayes::sim linkage)

// Port of R two_block_align_b_to_xhyper().
Rcpp::NumericVector two_block_align_b_to_xhyper_cpp(
    Rcpp::NumericVector b_vec,
    Rcpp::NumericMatrix X_k,
    Rcpp::CharacterVector group_levels
) {
  const Rcpp::CharacterVector rn = matrix_rownames(X_k);
  if (rn.size() == 0 || Rf_isNull(rn)) {
    if (b_vec.size() != X_k.nrow()) {
      Rcpp::stop(
        "length(b) (%d) must equal nrow(X_hyper) (%d) when X_hyper has no rownames.",
        b_vec.size(), X_k.nrow()
      );
    }
    return b_vec;
  }

  if (b_vec.hasAttribute("names") && !Rf_isNull(b_vec.attr("names"))) {
    Rcpp::CharacterVector bnames = b_vec.attr("names");
    std::string miss;
    for (int r = 0; r < rn.size(); ++r) {
      if (Rcpp::CharacterVector::is_na(rn[r])) {
        continue;
      }
      const std::string lev = Rcpp::as<std::string>(rn[r]);
      if (!charvec_contains_name(bnames, lev)) {
        if (!miss.empty()) {
          miss += ", ";
        }
        miss += lev;
      }
    }
    if (!miss.empty()) {
      Rcpp::stop("Group level(s) missing from b: %s", miss.c_str());
    }
    Rcpp::NumericVector out(rn.size());
    for (int r = 0; r < rn.size(); ++r) {
      if (Rcpp::CharacterVector::is_na(rn[r])) {
        out[r] = NA_REAL;
        continue;
      }
      out[r] = charvec_named_value(
        b_vec, bnames, Rcpp::as<std::string>(rn[r])
      );
    }
    return out;
  }

  if (b_vec.size() != group_levels.size() || b_vec.size() != X_k.nrow()) {
    Rcpp::stop(
      "b and X_hyper row counts disagree (b: %d, X_hyper: %d, group_levels: %d).",
      b_vec.size(), X_k.nrow(), group_levels.size()
    );
  }

  std::string miss;
  for (int r = 0; r < rn.size(); ++r) {
    if (Rcpp::CharacterVector::is_na(rn[r])) {
      continue;
    }
    const std::string lev = Rcpp::as<std::string>(rn[r]);
    if (!charvec_contains_name(group_levels, lev)) {
      if (!miss.empty()) {
        miss += ", ";
      }
      miss += lev;
    }
  }
  if (!miss.empty()) {
    Rcpp::stop(
      "X_hyper rownames do not match group_levels; missing in groups: %s",
      miss.c_str()
    );
  }

  Rcpp::NumericVector out(rn.size());
  for (int r = 0; r < rn.size(); ++r) {
    if (Rcpp::CharacterVector::is_na(rn[r])) {
      out[r] = NA_REAL;
      continue;
    }
    const std::string lev = Rcpp::as<std::string>(rn[r]);
    out[r] = charvec_named_value(b_vec, group_levels, lev);
  }
  return out;
}

// Lazy stats::gaussian() for rglmb Block~2 calls (matches R batch driver).
Rcpp::List two_block_gaussian_family_r() {
  Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function gaussian = stats["gaussian"];
  return Rcpp::as<Rcpp::List>(gaussian());
}

// .two_block_rglmb_iter_count() — envelope candidates from one rglmb draw.
int two_block_rglmb_iter_count_r(const Rcpp::List& fit_k) {
  if (!fit_k.containsElementNamed("iters") || Rf_isNull(fit_k["iters"])) {
    return 1;
  }
  SEXP it = fit_k["iters"];
  if (Rcpp::is<Rcpp::NumericMatrix>(it)) {
    Rcpp::NumericMatrix M(it);
    return static_cast<int>(M(0, 0));
  }
  Rcpp::NumericVector v(it);
  return static_cast<int>(v[0]);
}

// Block~2 via rglmb(..., pfamily = pf) — same path as .two_block_block2_one_chain.
Rcpp::NumericVector two_block_block2_rglmb_gamma(
    const Rcpp::NumericMatrix& b_i,
    int col_j,
    const Rcpp::NumericMatrix& X_j,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::List& pfamily_j,
    bool is_ing,
    double& tau2_j,
    double& iters_add
) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("glmbayesCore");
  Rcpp::Function rglmb = pkg["rglmb"];

  const int J = group_levels.size();
  Rcpp::NumericVector b_col(J);
  for (int g = 0; g < J; ++g) {
    b_col[g] = b_i(g, col_j);
  }
  b_col.attr("names") = group_levels;

  Rcpp::NumericVector y_j = two_block_align_b_to_xhyper_cpp(
    b_col, X_j, group_levels
  );

  Rcpp::List fit_k = rglmb(
    Rcpp::Named("n") = 1,
    Rcpp::Named("y") = y_j,
    Rcpp::Named("x") = X_j,
    Rcpp::Named("family") = two_block_gaussian_family_r(),
    Rcpp::Named("pfamily") = pfamily_j,
    Rcpp::Named("verbose") = false,
    Rcpp::Named("use_parallel") = false
  );

  if (is_ing) {
    Rcpp::NumericVector disp_j =
      Rcpp::as<Rcpp::NumericVector>(fit_k["dispersion"]);
    tau2_j = disp_j[0];
    iters_add = static_cast<double>(two_block_rglmb_iter_count_r(fit_k));
  } else {
    iters_add = 1.0;
  }

  return two_block_coef_mode_from_rNormalReg(fit_k);
}

// Port of R two_block_block2_one_chain() for one replicate chain.
Rcpp::List two_block_block2_one_chain_cpp_export(
    const Rcpp::NumericMatrix& b_i,
    const Rcpp::List& fixef_rows,
    const Rcpp::NumericVector& tau2_i,
    const Rcpp::NumericVector& iters_i,
    const Rcpp::List& x_hyper,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::List& pfamily_list,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names
) {
  const int p_re = re_names.size();
  Rcpp::List fixef_out = Rcpp::clone(fixef_rows);
  Rcpp::NumericVector tau2_out = Rcpp::clone(tau2_i);
  Rcpp::NumericVector iters_out = Rcpp::clone(iters_i);

  for (int k = 0; k < p_re; ++k) {
    const std::string re_k = Rcpp::as<std::string>(re_names[k]);
    Rcpp::NumericMatrix X_k = Rcpp::as<Rcpp::NumericMatrix>(x_hyper[re_k]);
    Rcpp::List pfamily_k = Rcpp::as<Rcpp::List>(pfamily_list[re_k]);
    const std::string ptype_k = Rcpp::as<std::string>(ptypes[k]);
    const bool is_ing = (ptype_k == "dIndependent_Normal_Gamma");

    double tau2_k = tau2_out[k];
    double it_add = 0.0;
    Rcpp::NumericVector gamma_k = two_block_block2_rglmb_gamma(
      b_i,
      k,
      X_k,
      group_levels,
      pfamily_k,
      is_ing,
      tau2_k,
      it_add
    );

    Rcpp::NumericVector row_k = Rcpp::as<Rcpp::NumericVector>(fixef_out[k]);
    two_block_update_fixef_row(row_k, gamma_k);
    fixef_out[k] = row_k;
    if (is_ing) {
      tau2_out[k] = tau2_k;
      iters_out[k] = iters_out[k] + it_add;
    } else {
      iters_out[k] = iters_out[k] + 1.0;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("fixef") = fixef_out,
    Rcpp::Named("tau2") = tau2_out,
    Rcpp::Named("iters") = iters_out
  );
}

namespace {

void two_block_reorder_b_to_group_levels(
    Rcpp::NumericMatrix& b_i,
    const Rcpp::CharacterVector& block_ids,
    const Rcpp::CharacterVector& group_levels
) {
  b_i = glmbayes::sim::two_block_reorder_b_to_group_levels(
    b_i, block_ids, group_levels
  );
}

void store_b_chain(
    Rcpp::NumericVector& b_store,
    int i,
    int J,
    int p_re,
    const Rcpp::NumericMatrix& b_i
) {
  for (int j = 0; j < p_re; ++j) {
    for (int g = 0; g < J; ++g) {
      b_store[g + J * (j + p_re * i)] = b_i(g, j);
    }
  }
}

void load_b_chain(
    const Rcpp::NumericVector& b_store,
    int i,
    int J,
    int p_re,
    Rcpp::NumericMatrix& b_i
) {
  b_i = Rcpp::NumericMatrix(J, p_re);
  for (int j = 0; j < p_re; ++j) {
    for (int g = 0; g < J; ++g) {
      b_i(g, j) = b_store[g + J * (j + p_re * i)];
    }
  }
}

void pack_two_block_chain_draw(
    int i,
    int p_re,
    int J,
    const std::vector<Rcpp::NumericVector>& fixef,
    const std::vector<double>& tau2,
    const Rcpp::NumericMatrix& b_i,
    Rcpp::List& fixef_draws,
    Rcpp::NumericVector& b_arr,
    Rcpp::NumericMatrix& disp_draws
) {
  for (int j = 0; j < p_re; ++j) {
    Rcpp::NumericMatrix fd = fixef_draws[j];
    const Rcpp::NumericVector& fj = fixef[static_cast<size_t>(j)];
    if (fj.size() != fd.ncol()) {
      Rcpp::stop(
        "Block 2 draw for component %d has length %d; expected %d.",
        j + 1, fj.size(), fd.ncol()
      );
    }
    fd(i, Rcpp::_) = fj;
  }
  for (int j = 0; j < p_re; ++j) {
    for (int g = 0; g < J; ++g) {
      b_arr[g + J * (j + p_re * i)] = b_i(g, j);
    }
    disp_draws(i, j) = tau2[static_cast<size_t>(j)];
  }
}

} // anonymous namespace

List two_block_rNormal_reg_v2_cpp_export(
    int n,
    int m_convergence,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    const List& x_hyper,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const List& pfamily_list,
    const List& fixef_start,
    const CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Function& f2,
    const Function& f3,
    const Function& f2_gauss,
    const Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
) {
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (m_convergence < 1) {
    Rcpp::stop("'m_convergence' must be at least 1.");
  }
  const int p_re = x.ncol();
  const int J = group_levels.size();
  if (x_hyper.size() != p_re) {
    Rcpp::stop("length(x_hyper) must equal ncol(x) = %d.", p_re);
  }
  if (pfamily_list.size() != p_re) {
    Rcpp::stop("length(pfamily_list) must equal ncol(x) = %d.", p_re);
  }
  if (fixef_start.size() != p_re) {
    Rcpp::stop("length(fixef_start) must equal ncol(x) = %d.", p_re);
  }

  const bool is_gaussian = (family == "gaussian");

  MuAllBuilder mu_builder(x_hyper, group_levels);

  // Block 2 priors parsed once (constant across draws; no RNG consumed).
  std::vector<Block2PriorV2> pr2(p_re);
  bool any_ing = false;
  for (int j = 0; j < p_re; ++j) {
    const NumericMatrix& X_j = mu_builder.X[j];
    pr2[j] = block2_prior_prep_v2(List(pfamily_list[j]), j + 1, X_j.ncol());
    if (pr2[j].is_ing) any_ing = true;
  }

  // fixef state (component order fixed by the R wrapper)
  std::vector<NumericVector> fixef_start_v(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_start_v[j] = Rcpp::as<NumericVector>(fixef_start[j]);
  }
  std::vector<NumericVector> fixef = fixef_start_v;

  // tau2 working vector: dNormal components stay fixed; ING components are
  // re-drawn jointly with gamma_k each sweep (and reset per stored draw).
  std::vector<double> tau2_start(p_re);
  for (int j = 0; j < p_re; ++j) {
    tau2_start[j] = pr2[j].dispersion;
  }
  std::vector<double> tau2 = tau2_start;

  // Base Block 1 prior precision (p_re x p_re), needed only when ING
  // components feed the current tau2 back into Block 1 each sweep.
  // ING rows/cols are overridden with diag(1/tau2_j); dNormal rows keep
  // their fixed entries from prior_list_block1.
  NumericMatrix base_P1;
  if (any_ing) {
    if (has_non_null(prior_list_block1, "P")) {
      base_P1 = Rcpp::as<NumericMatrix>(prior_list_block1["P"]);
    } else if (has_non_null(prior_list_block1, "Sigma")) {
      NumericMatrix S1 = Rcpp::as<NumericMatrix>(prior_list_block1["Sigma"]);
      arma::mat S1a(const_cast<double*>(S1.begin()),
                    S1.nrow(), S1.ncol(), false);
      arma::mat P1inv = arma::inv_sympd(S1a);
      base_P1 = NumericMatrix(Rcpp::wrap(0.5 * (P1inv + P1inv.t())));
    } else {
      Rcpp::stop("prior_list_block1 must contain 'P' or 'Sigma'.");
    }
    if (base_P1.nrow() != p_re || base_P1.ncol() != p_re) {
      Rcpp::stop(
        "prior_list_block1 P/Sigma must be %d x %d.", p_re, p_re
      );
    }
  }

  // Storage: fixef draws per component (n x q_k); b draws as J x p_re x n;
  // tau2 values at each stored draw (n x p_re); total Block 2 candidates
  // generated per stored draw (n x p_re; summed over the m_convergence
  // inner sweeps -- ING components count envelope candidates until
  // acceptance, dNormal components count 1 per conjugate draw).
  List fixef_draws(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_draws[j] = NumericMatrix(n, fixef_start_v[j].size());
  }
  NumericVector b_arr(Rcpp::Dimension(J, p_re, n));
  NumericMatrix disp_draws(n, p_re);
  NumericMatrix iters_draws(n, p_re);

  NumericMatrix mu_all;
  NumericMatrix b_i;
  CharacterVector group_ids;
  bool have_ids = false;

  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar) {
      glmbayes::progress::progress_bar(
        static_cast<double>(i + 1), static_cast<double>(n)
      );
    }

    fixef = fixef_start_v;
    tau2 = tau2_start;

    for (int m = 0; m < m_convergence; ++m) {

      mu_all = mu_builder.build(fixef);
      List pl1;
      if (any_ing) {
        // Refresh the Block 1 prior precision from the current tau2: ING
        // rows/cols become diag(1/tau2_j); dNormal rows keep base entries.
        NumericMatrix P1 = Rcpp::clone(base_P1);
        for (int j = 0; j < p_re; ++j) {
          if (!pr2[j].is_ing) continue;
          for (int c = 0; c < p_re; ++c) {
            P1(j, c) = 0.0;
            P1(c, j) = 0.0;
          }
          P1(j, j) = 1.0 / tau2[j];
        }
        List pl1_base = List::create(
          Rcpp::Named("P") = P1
        );
        pl1 = block1_prior_list(
          mu_all, pl1_base, dispersion_block1, ddef_block1
        );
      } else {
        pl1 = block1_prior_list(
          mu_all, prior_list_block1, dispersion_block1, ddef_block1
        );
      }

      List block_i;
      if (is_gaussian) {
        block_i = block_rNormalReg_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3, 2
        );
      } else {
        block_i = block_rNormalGLM_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3,
          family, link, Gridtype, n_envopt,
          use_parallel, use_opencl, verbose
        );
      }

      b_i = Rcpp::as<NumericMatrix>(block_i["coefficients"]);
      if (b_i.nrow() != J || b_i.ncol() != p_re) {
        Rcpp::stop(
          "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
          b_i.nrow(), b_i.ncol(), J, p_re
        );
      }
      if (!have_ids) {
        List bi_info = block_i["block_info"];
        group_ids = Rcpp::as<CharacterVector>(bi_info["ids"]);
        have_ids = true;
      }

      // Block 2: per-component dispatch on the pfamily type.
      for (int j = 0; j < p_re; ++j) {
        const NumericMatrix& X_j = mu_builder.X[j];
        if (X_j.nrow() != b_i.nrow()) {
          Rcpp::stop(
            "nrow(x_hyper[[%d]]) (%d) must equal nrow(y) (%d).",
            j + 1, X_j.nrow(), b_i.nrow()
          );
        }
        NumericVector y_j = b_i(Rcpp::_, j);
        NumericVector offset_j(X_j.nrow(), 0.0);
        NumericVector wt_j(X_j.nrow(), 1.0);

        if (!pr2[j].is_ing) {
          // dNormal: conjugate draw at fixed dispersion (identical to v1).
          List out_j = rNormalReg(
            1, y_j, X_j, pr2[j].mu, pr2[j].P, offset_j, wt_j,
            pr2[j].dispersion, f2_gauss, f3_gauss, pr2[j].mu,
            "gaussian", "identity", 2
          );
          NumericMatrix coef_j =
            Rcpp::as<NumericMatrix>(out_j["coefficients"]);
          fixef[j] = NumericVector(coef_j(0, Rcpp::_));
          iters_draws(i, j) += 1.0;
        } else {
          // ING: joint (gamma_k, tau2_k) draw via the likelihood-subgradient
          // envelope sampler (same path as rglmb with a
          // dIndependent_Normal_Gamma pfamily; rIndepNormalGammaReg, n = 1).
          List out_j = rIndepNormalGammaReg(
            1, y_j, X_j, pr2[j].mu, pr2[j].P, offset_j, wt_j,
            pr2[j].shape, pr2[j].rate, pr2[j].max_disp_perc,
            Rcpp::Nullable<NumericVector>(pr2[j].disp_lower),
            Rcpp::Nullable<NumericVector>(pr2[j].disp_upper),
            2 /*Gridtype*/, 1 /*n_envopt*/,
            false /*use_parallel: n = 1 is serial anyway*/,
            false /*use_opencl*/, false /*verbose*/, false /*progbar*/
          );
          // "out" is p x n (coefficients in columns, original scale).
          NumericMatrix coef_j = Rcpp::as<NumericMatrix>(out_j["out"]);
          NumericVector gamma_j(coef_j.nrow());
          for (int c = 0; c < coef_j.nrow(); ++c) {
            gamma_j[c] = coef_j(c, 0);
          }
          fixef[j] = gamma_j;
          NumericVector disp_j = Rcpp::as<NumericVector>(out_j["disp_out"]);
          tau2[j] = disp_j[0];
          // Candidates generated by the envelope accept-reject sampler for
          // this accepted joint draw (>= 1; reading it consumes no RNG).
          NumericVector iters_j =
            Rcpp::as<NumericVector>(out_j["iters_out"]);
          iters_draws(i, j) += iters_j[0];
        }
      }
    }

    for (int j = 0; j < p_re; ++j) {
      NumericMatrix fd = fixef_draws[j];
      const NumericVector& fj = fixef[j];
      if (fj.size() != fd.ncol()) {
        Rcpp::stop(
          "Block 2 draw for component %d has length %d; expected %d.",
          j + 1, fj.size(), fd.ncol()
        );
      }
      fd(i, Rcpp::_) = fj;
    }
    for (int j = 0; j < p_re; ++j) {
      for (int g = 0; g < J; ++g) {
        b_arr[g + J * (j + p_re * i)] = b_i(g, j);
      }
      disp_draws(i, j) = tau2[j];
    }
  }

  if (progbar) {
    glmbayes::progress::progress_bar_finish();
  }

  List fixef_last(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_last[j] = fixef[j];
  }

  return List::create(
    Rcpp::Named("fixef_draws") = fixef_draws,
    Rcpp::Named("b_draws") = b_arr,
    Rcpp::Named("fixef_last") = fixef_last,
    Rcpp::Named("b_last") = b_i,
    Rcpp::Named("mu_all_last") = mu_all,
    Rcpp::Named("group_ids") = group_ids,
    Rcpp::Named("dispersion_fixef_draws") = disp_draws,
    Rcpp::Named("iters_fixef_draws") = iters_draws,
    Rcpp::Named("any_ing") = any_ing
  );
}

List two_block_rNormal_reg_v3_cpp_export(
    int n,
    int m_convergence,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    const List& x_hyper,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const List& pfamily_list,
    const List& fixef_start,
    const CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Function& f2,
    const Function& f3,
    const Function& f2_gauss,
    const Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    Rcpp::Nullable<int> seed,
    int seed_offset,
    bool progbar
) {
  // Two-block Gibbs sampler (v3): independent short chains in C++.
  //
  // Indexing:
  //   p_re  = ncol(x) = number of random-effect *columns* in Z (e.g. intercept,
  //           slope RE terms), not the number of hyperparameters.
  //   J     = number of grouping-factor *levels* (rows of b_i).
  //   n     = number of independent chains (stored draws / replicates).
  //   m     = inner Gibbs sweeps per chain (m_convergence); only the state
  //           after the last sweep is stored for each chain.
  //
  // Loop hierarchy (outer -> inner):
  //   for i in 0..n-1     independent chains (sequential draws from R RNG)
  //     for m in 0..m_convergence-1   one full two-block sweep
  //       Block 1 once    joint draw of all group-level b (J x p_re matrix)
  //       for j in 0..p_re-1   Block 2 per RE column (hyperparameter gamma_k)
  //
  // Block 1 (level-1 / observation model): given hyperparameters fixef, draw
  //   random effects b_{g,k} for every group g and RE column k from the
  //   conditional of y | b (Gaussian reg or GLM blocked by the factor).
  //   block_rNormalReg / block_rNormalGLM handle grouping internally; there
  //   is no explicit loop over factor levels here.
  //
  // Block 2 (level-2 / hyperparameter model): given Block-1 draws b[,k] as
  //   pseudo-response (length J), draw hyperparameters gamma_k via regression
  //   on x_hyper[[k]] (one pfamily component per RE column j).

  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (m_convergence < 1) {
    Rcpp::stop("'m_convergence' must be at least 1.");
  }
  const int p_re = x.ncol();
  const int J = group_levels.size();
  if (x_hyper.size() != p_re) {
    Rcpp::stop("length(x_hyper) must equal ncol(x) = %d.", p_re);
  }
  if (pfamily_list.size() != p_re) {
    Rcpp::stop("length(pfamily_list) must equal ncol(x) = %d.", p_re);
  }
  if (fixef_start.size() != p_re) {
    Rcpp::stop("length(fixef_start) must equal ncol(x) = %d.", p_re);
  }

  const bool is_gaussian = (family == "gaussian");

  MuAllBuilder mu_builder(x_hyper, group_levels);

  // Setup (no RNG): parse Block 2 pfamily priors, one entry per RE column j.
  std::vector<Block2PriorV2> pr2(p_re);
  bool any_ing = false;
  for (int j = 0; j < p_re; ++j) {
    const NumericMatrix& X_j = mu_builder.X[j];
    pr2[j] = block2_prior_prep_v2(List(pfamily_list[j]), j + 1, X_j.ncol());
    if (pr2[j].is_ing) any_ing = true;
  }

  // Deep snapshot: each chain resets from this copy (same as a fresh n = 1
  // v2 call with the current R RNG stream).
  std::vector<NumericVector> fixef_start_v(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_start_v[j] =
      Rcpp::clone(Rcpp::as<NumericVector>(fixef_start[j]));
  }
  std::vector<NumericVector> fixef(p_re);

  // Working state for Block 2 hyperparameters (gamma_k) and ING dispersions.
  // fixef[j] = hyperparameter vector for RE column j; tau2[j] = Block 2
  // dispersion for that component (fixed for dNormal, updated for ING).
  std::vector<double> tau2_start(p_re);
  for (int j = 0; j < p_re; ++j) {
    tau2_start[j] = pr2[j].dispersion;
  }
  std::vector<double> tau2 = tau2_start;

  // Block 1 prior template for ING: p_re x p_re precision; ING rows refreshed
  // each sweep from current tau2 (see inner loop over j below).
  NumericMatrix base_P1;
  if (any_ing) {
    if (has_non_null(prior_list_block1, "P")) {
      base_P1 = Rcpp::as<NumericMatrix>(prior_list_block1["P"]);
    } else if (has_non_null(prior_list_block1, "Sigma")) {
      NumericMatrix S1 = Rcpp::as<NumericMatrix>(prior_list_block1["Sigma"]);
      arma::mat S1a(const_cast<double*>(S1.begin()),
                    S1.nrow(), S1.ncol(), false);
      arma::mat P1inv = arma::inv_sympd(S1a);
      base_P1 = NumericMatrix(Rcpp::wrap(0.5 * (P1inv + P1inv.t())));
    } else {
      Rcpp::stop("prior_list_block1 must contain 'P' or 'Sigma'.");
    }
    if (base_P1.nrow() != p_re || base_P1.ncol() != p_re) {
      Rcpp::stop(
        "prior_list_block1 P/Sigma must be %d x %d.", p_re, p_re
      );
    }
  }

  // Output buffers:
  //   fixef_draws[[j]]     n x q_k  Block 2 hyperparameters per chain
  //   b_arr                J x p_re x n  Block 1 random effects per chain
  //   disp_draws, iters_draws   n x p_re  Block 2 tau^2 and sampler counts
  List fixef_draws(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_draws[j] = NumericMatrix(n, fixef_start_v[j].size());
  }
  NumericVector b_arr(Rcpp::Dimension(J, p_re, n));
  NumericMatrix disp_draws(n, p_re);
  NumericMatrix iters_draws(n, p_re);

  NumericMatrix mu_all;
  NumericMatrix b_i;
  CharacterVector group_ids;
  bool have_ids = false;

  const bool have_seed = seed.isNotNull();
  two_block_driver_banner_v3(n, m_convergence, seed_offset, have_seed);
  // Match run_short_chains_v2: R chain bar when n > 1, inner bar when n == 1.
  const bool chain_progbar = progbar && n > 1;
  const bool inner_progbar = progbar && n == 1;

  // ---- Outer loop: n independent short chains (stored draws) ----
  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();

    if (have_seed) {
      v3_reseed_r(Rcpp::as<int>(seed.get()) + seed_offset + i + 1);
    }

    if (chain_progbar || inner_progbar) {
      glmbayes::progress::progress_bar(
        static_cast<double>(i + 1), static_cast<double>(n)
      );
    }

    // Reset working state to fixef_start at the beginning of each chain
    // (same as a fresh two_block_rNormal_reg_v2(n = 1) call).
    for (int j = 0; j < p_re; ++j) {
      fixef[j] = Rcpp::clone(fixef_start_v[j]);
    }
    tau2 = tau2_start;

    // ---- Inner loop: m_convergence full two-block Gibbs sweeps ----
    // Each sweep: Block 1 (all b) then Block 2 (each RE column j).  Only the
    // final sweep's state is retained and packed after this loop completes.
    for (int m = 0; m < m_convergence; ++m) {

      // Build group-level means mu_{k,g} from current fixef (Block 2 -> Block 1).
      mu_all = mu_builder.build(fixef);
      List pl1;
      if (any_ing) {
        // Refresh Block 1 prior precision from ING tau2_k (loop over RE cols,
        // not factor levels).
        NumericMatrix P1 = Rcpp::clone(base_P1);
        for (int j = 0; j < p_re; ++j) {
          if (!pr2[j].is_ing) continue;
          for (int c = 0; c < p_re; ++c) {
            P1(j, c) = 0.0;
            P1(c, j) = 0.0;
          }
          P1(j, j) = 1.0 / tau2[j];
        }
        List pl1_base = List::create(
          Rcpp::Named("P") = P1
        );
        pl1 = block1_prior_list(
          mu_all, pl1_base, dispersion_block1, ddef_block1
        );
      } else {
        pl1 = block1_prior_list(
          mu_all, prior_list_block1, dispersion_block1, ddef_block1
        );
      }

      // ---- Block 1: one joint draw of b (J groups x p_re RE columns) ----
      // Conditions on observation-level y, design Z, and current hyperprior pl1.
      // Grouping factor levels are handled inside the block_* export.
      List block_i;
      if (is_gaussian) {
        block_i = block_rNormalReg_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3, 2
        );
      } else {
        block_i = block_rNormalGLM_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3,
          family, link, Gridtype, n_envopt,
          use_parallel, use_opencl, verbose
        );
      }

      b_i = Rcpp::as<NumericMatrix>(block_i["coefficients"]);
      if (b_i.nrow() != J || b_i.ncol() != p_re) {
        Rcpp::stop(
          "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
          b_i.nrow(), b_i.ncol(), J, p_re
        );
      }
      if (!have_ids) {
        List bi_info = block_i["block_info"];
        group_ids = Rcpp::as<CharacterVector>(bi_info["ids"]);
        have_ids = true;
      }

      // ---- Block 2: one hyperparameter draw per RE column j ----
      // Pseudo-response y_j = b_i[, j] (length J, one value per group level).
      // x_hyper[[j]] is the level-2 design for RE column j.
      for (int j = 0; j < p_re; ++j) {
        const NumericMatrix& X_j = mu_builder.X[j];
        if (X_j.nrow() != b_i.nrow()) {
          Rcpp::stop(
            "nrow(x_hyper[[%d]]) (%d) must equal nrow(y) (%d).",
            j + 1, X_j.nrow(), b_i.nrow()
          );
        }
        NumericVector y_j = b_i(Rcpp::_, j);
        NumericVector offset_j(X_j.nrow(), 0.0);
        NumericVector wt_j(X_j.nrow(), 1.0);

        if (!pr2[j].is_ing) {
          // dNormal Block 2: conjugate draw of gamma_k at fixed tau2[j].
          List out_j = rNormalReg(
            1, y_j, X_j, pr2[j].mu, pr2[j].P, offset_j, wt_j,
            pr2[j].dispersion, f2_gauss, f3_gauss, pr2[j].mu,
            "gaussian", "identity", 2
          );
          NumericMatrix coef_j =
            Rcpp::as<NumericMatrix>(out_j["coefficients"]);
          fixef[j] = NumericVector(coef_j(0, Rcpp::_));
          iters_draws(i, j) += 1.0;
        } else {
          // ING Block 2: joint draw of (gamma_k, tau2_k); tau2[j] feeds next
          // sweep's Block 1 prior row via the P1 refresh above.
          List out_j = rIndepNormalGammaReg(
            1, y_j, X_j, pr2[j].mu, pr2[j].P, offset_j, wt_j,
            pr2[j].shape, pr2[j].rate, pr2[j].max_disp_perc,
            Rcpp::Nullable<NumericVector>(pr2[j].disp_lower),
            Rcpp::Nullable<NumericVector>(pr2[j].disp_upper),
            2 /*Gridtype*/, 1 /*n_envopt*/,
            false /*use_parallel: n = 1 is serial anyway*/,
            false /*use_opencl*/, false /*verbose*/, false /*progbar*/
          );
          NumericMatrix coef_j = Rcpp::as<NumericMatrix>(out_j["out"]);
          NumericVector gamma_j(coef_j.nrow());
          for (int c = 0; c < coef_j.nrow(); ++c) {
            gamma_j[c] = coef_j(c, 0);
          }
          fixef[j] = gamma_j;
          NumericVector disp_j = Rcpp::as<NumericVector>(out_j["disp_out"]);
          tau2[j] = disp_j[0];
          NumericVector iters_j =
            Rcpp::as<NumericVector>(out_j["iters_out"]);
          iters_draws(i, j) += iters_j[0];
        }
      }
    }

    // Pack chain i: store Block 2 hyperparameters and Block 1 b after the
    // final inner sweep (m = m_convergence - 1).
    for (int j = 0; j < p_re; ++j) {
      NumericMatrix fd = fixef_draws[j];
      const NumericVector& fj = fixef[j];
      if (fj.size() != fd.ncol()) {
        Rcpp::stop(
          "Block 2 draw for component %d has length %d; expected %d.",
          j + 1, fj.size(), fd.ncol()
        );
      }
      fd(i, Rcpp::_) = fj;
    }
    for (int j = 0; j < p_re; ++j) {
      for (int g = 0; g < J; ++g) {
        // b_arr[g, j, i]: random effect for group level g, RE column j, chain i.
        b_arr[g + J * (j + p_re * i)] = b_i(g, j);
      }
      disp_draws(i, j) = tau2[j];
    }
  }

  if (chain_progbar || inner_progbar) {
    glmbayes::progress::progress_bar_finish();
  }

  List fixef_last(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_last[j] = fixef[j];
  }

  return List::create(
    Rcpp::Named("fixef_draws") = fixef_draws,
    Rcpp::Named("b_draws") = b_arr,
    Rcpp::Named("fixef_last") = fixef_last,
    Rcpp::Named("b_last") = b_i,
    Rcpp::Named("mu_all_last") = mu_all,
    Rcpp::Named("group_ids") = group_ids,
    Rcpp::Named("dispersion_fixef_draws") = disp_draws,
    Rcpp::Named("iters_fixef_draws") = iters_draws,
    Rcpp::Named("any_ing") = any_ing
  );
}


List two_block_rNormal_reg_v4_cpp_export(
    int n,
    int m_convergence,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    const List& x_hyper,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const List& pfamily_list,
    const List& fixef_start,
    const CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Function& f2,
    const Function& f3,
    const Function& f2_gauss,
    const Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    Rcpp::Nullable<int> seed,
    int seed_offset,
    bool progbar,
    std::string stage_label,
    bool diag_sweeps,
    SEXP fixef_mode,
    SEXP b_mode
) {
  // Two-block Gibbs sampler (v4): independent short chains in C++.
  //
  // Indexing:
  //   p_re  = ncol(x) = number of random-effect *columns* in Z (e.g. intercept,
  //           slope RE terms), not the number of hyperparameters.
  //   J     = number of grouping-factor *levels* (rows of b_i).
  //   n     = number of independent chains (stored draws / replicates).
  //   m     = inner Gibbs sweeps per chain (m_convergence); only the state
  //           after the last sweep is stored for each chain.
  //
  // Loop hierarchy (outer -> inner):
  //   for i in 0..n-1     independent chains (sequential draws from R RNG)
  //     for m in 0..m_convergence-1   one full two-block sweep
  //       Block 1 once    joint draw of all group-level b (J x p_re matrix)
  //       for j in 0..p_re-1   Block 2 per RE column (hyperparameter gamma_k)
  //
  // Block 1 (level-1 / observation model): given hyperparameters fixef, draw
  //   random effects b_{g,k} for every group g and RE column k from the
  //   conditional of y | b (Gaussian reg or GLM blocked by the factor).
  //   block_rNormalReg / block_rNormalGLM handle grouping internally; there
  //   is no explicit loop over factor levels here.
  //
  // Block 2 (level-2 / hyperparameter model): given Block-1 draws b[,k] as
  //   pseudo-response (length J), draw hyperparameters gamma_k via regression
  //   on x_hyper[[k]] (one pfamily component per RE column j).
  
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (m_convergence < 1) {
    Rcpp::stop("'m_convergence' must be at least 1.");
  }
  const int p_re = x.ncol();
  const int J = group_levels.size();
  if (x_hyper.size() != p_re) {
    Rcpp::stop("length(x_hyper) must equal ncol(x) = %d.", p_re);
  }
  if (pfamily_list.size() != p_re) {
    Rcpp::stop("length(pfamily_list) must equal ncol(x) = %d.", p_re);
  }
  if (fixef_start.size() != p_re) {
    Rcpp::stop("length(fixef_start) must equal ncol(x) = %d.", p_re);
  }
  
  const bool is_gaussian = (family == "gaussian");
  
  MuAllBuilder mu_builder(x_hyper, group_levels);
  
  // Setup (no RNG): parse Block 2 pfamily priors, one entry per RE column j.
  std::vector<Block2PriorV2> pr2(p_re);
  bool any_ing = false;
  for (int j = 0; j < p_re; ++j) {
    const NumericMatrix& X_j = mu_builder.X[j];
    pr2[j] = block2_prior_prep_v2(List(pfamily_list[j]), j + 1, X_j.ncol());
    if (pr2[j].is_ing) any_ing = true;
  }
  
  // Deep snapshot: each chain resets from this copy (same as a fresh n = 1
  // v2 call with the current R RNG stream).
  std::vector<NumericVector> fixef_start_v(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_start_v[j] =
      Rcpp::clone(Rcpp::as<NumericVector>(fixef_start[j]));
  }

  std::vector<NumericVector> fixef_mode_v(p_re);
  if (diag_sweeps && !Rf_isNull(fixef_mode)) {
    List fixef_mode_in = Rcpp::as<List>(fixef_mode);
    if (fixef_mode_in.size() != p_re) {
      Rcpp::stop("length(fixef_mode) must equal ncol(x) = %d.", p_re);
    }
    for (int j = 0; j < p_re; ++j) {
      fixef_mode_v[j] =
        Rcpp::clone(Rcpp::as<NumericVector>(fixef_mode_in[j]));
    }
  } else {
    fixef_mode_v = fixef_start_v;
  }

  List fixef_mode_list_for_print;
  if (diag_sweeps && !Rf_isNull(fixef_mode)) {
    fixef_mode_list_for_print = Rcpp::as<List>(fixef_mode);
  } else {
    fixef_mode_list_for_print = fixef_start;
  }

  NumericMatrix b_mode_mat;
  bool have_b_mode = false;
  if (diag_sweeps && !Rf_isNull(b_mode)) {
    b_mode_mat = Rcpp::as<NumericMatrix>(b_mode);
    if (b_mode_mat.nrow() == J && b_mode_mat.ncol() == p_re) {
      have_b_mode = true;
    }
  }

  CharacterVector re_names;
  if (x.hasAttribute("dimnames")) {
    List dn = x.attr("dimnames");
    if (!Rf_isNull(dn[1])) {
      re_names = dn[1];
    }
  }
  if (re_names.size() != p_re) {
    re_names = CharacterVector(p_re);
    for (int j = 0; j < p_re; ++j) {
      re_names[j] = "RE" + std::to_string(j + 1);
    }
  }

  std::vector<NumericVector> fixef(p_re);
  
  // Working state for Block 2 hyperparameters (gamma_k) and ING dispersions.
  // fixef[j] = hyperparameter vector for RE column j; tau2[j] = Block 2
  // dispersion for that component (fixed for dNormal, updated for ING).
  std::vector<double> tau2_start(p_re);
  for (int j = 0; j < p_re; ++j) {
    tau2_start[j] = pr2[j].dispersion;
  }
  std::vector<double> tau2 = tau2_start;
  
  // Block 1 prior template for ING: p_re x p_re precision; ING rows refreshed
  // each sweep from current tau2 (see inner loop over j below).
  NumericMatrix base_P1;
  if (any_ing) {
    if (has_non_null(prior_list_block1, "P")) {
      base_P1 = Rcpp::as<NumericMatrix>(prior_list_block1["P"]);
    } else if (has_non_null(prior_list_block1, "Sigma")) {
      NumericMatrix S1 = Rcpp::as<NumericMatrix>(prior_list_block1["Sigma"]);
      arma::mat S1a(const_cast<double*>(S1.begin()),
                    S1.nrow(), S1.ncol(), false);
      arma::mat P1inv = arma::inv_sympd(S1a);
      base_P1 = NumericMatrix(Rcpp::wrap(0.5 * (P1inv + P1inv.t())));
    } else {
      Rcpp::stop("prior_list_block1 must contain 'P' or 'Sigma'.");
    }
    if (base_P1.nrow() != p_re || base_P1.ncol() != p_re) {
      Rcpp::stop(
        "prior_list_block1 P/Sigma must be %d x %d.", p_re, p_re
      );
    }
  }
  
  // Output buffers:
  //   fixef_draws[[j]]     n x q_k  Block 2 hyperparameters per chain
  //   b_arr                J x p_re x n  Block 1 random effects per chain
  //   disp_draws, iters_draws   n x p_re  Block 2 tau^2 and sampler counts
  List fixef_draws(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_draws[j] = NumericMatrix(n, fixef_start_v[j].size());
  }
  NumericVector b_arr(Rcpp::Dimension(J, p_re, n));
  NumericMatrix disp_draws(n, p_re);
  NumericMatrix iters_draws(n, p_re);
  
  NumericMatrix mu_all;
  NumericMatrix b_i;
  CharacterVector group_ids;
  bool have_ids = false;
  
  const bool have_seed = seed.isNotNull();
  // Match run_short_chains_v2: R chain bar when n > 1, inner bar when n == 1.
  const bool chain_progbar = progbar && n > 1;
  const bool inner_progbar = progbar && n == 1;
  
  
  int p_dim = 0;
  std::vector<int> q_k(p_re);
  for (int j = 0; j < p_re; ++j) {
    q_k[j] = fixef_start_v[j].size();
    p_dim += q_k[j];
  }

  NumericVector fixef_sum;
  NumericMatrix b_sum;
  if (diag_sweeps) {
    fixef_sum = NumericVector(p_dim, 0.0);
    b_sum = NumericMatrix(J, p_re);
    std::fill(b_sum.begin(), b_sum.end(), 0.0);
    if (stage_label.empty()) {
      stage_label = "v4";
    }
  }
  
  NumericMatrix fixef_temp=NumericMatrix(n,p_dim);
  NumericMatrix tau2_temp = NumericMatrix(n, p_re);
  
  
  for (int i = 0; i < n; ++i) {
    two_block_pack_fixef_row(fixef_temp, i, p_re, q_k, fixef_start_v);
    for (int j = 0; j < p_re; ++j) {
      tau2_temp(i, j) = tau2_start[j];
    }
  }

  if (diag_sweeps) {
    Rcpp::Rcout << "[" << stage_label << " at entry, n=" << n
                << ", fixef_start loaded into fixef_temp]\n";
    two_block_print_fixef_line(
      "  fixef_start (ICM mode):",
      fixef_mode_v, p_re, re_names
    );
  }

  two_block_driver_banner_v4(n, m_convergence, seed_offset, have_seed, p_dim);
  
  // ---- Outer loop: n independent short chains (stored draws) ----
  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();

    if (chain_progbar || inner_progbar) {
      glmbayes::progress::progress_bar(
        static_cast<double>(i + 1), static_cast<double>(n)
      );
    }

    // ---- Inner loop: m_convergence full two-block Gibbs sweeps ----
    // Each sweep: Block 1 (all b) then Block 2 (each RE column j).  Only the
    // final sweep's state is retained and packed after this loop completes.
    for (int m = 0; m < m_convergence; ++m) {

      if (diag_sweeps && i == 0) {
        std::fill(fixef_sum.begin(), fixef_sum.end(), 0.0);
        std::fill(b_sum.begin(), b_sum.end(), 0.0);
        Rcpp::Rcout << "--- glmbayesCore v4 sweep m=" << (m + 1)
                    << " / " << m_convergence << " ---\n";
        Rcpp::Rcout << "[" << stage_label << " sweep " << (m + 1)
                    << " / " << m_convergence << " BEFORE load (chain 1), n="
                    << n << "]\n";
      }

      if (m == 0) {
        if (have_seed) {
          v3_reseed_r(Rcpp::as<int>(seed.get()) + seed_offset + i + 1);
        }
        for (int j = 0; j < p_re; ++j) {
          fixef[j] = Rcpp::clone(fixef_start_v[j]);
        }
        tau2 = tau2_start;
      }

      two_block_unpack_fixef_row(
        fixef_temp, i, p_re, q_k, fixef, fixef_start_v
      );
      for (int j = 0; j < p_re; ++j) {
        tau2[j] = tau2_temp(i, j);
      }

      if (diag_sweeps && i == 0) {
        two_block_print_fixef_line(
          "  fixef after load (chain 1):",
          fixef, p_re, re_names
        );
      }
      
      // Build group-level means mu_{k,g} from current fixef (Block 2 -> Block 1).
      mu_all = mu_builder.build(fixef);
      List pl1;
      if (any_ing) {
        // Refresh Block 1 prior precision from ING tau2_k (loop over RE cols,
        // not factor levels).
        NumericMatrix P1 = Rcpp::clone(base_P1);
        for (int j = 0; j < p_re; ++j) {
          if (!pr2[j].is_ing) continue;
          for (int c = 0; c < p_re; ++c) {
            P1(j, c) = 0.0;
            P1(c, j) = 0.0;
          }
          P1(j, j) = 1.0 / tau2[j];
        }
        List pl1_base = List::create(
          Rcpp::Named("P") = P1
        );
        pl1 = block1_prior_list(
          mu_all, pl1_base, dispersion_block1, ddef_block1
        );
      } else {
        pl1 = block1_prior_list(
          mu_all, prior_list_block1, dispersion_block1, ddef_block1
        );
      }
      
      // ---- Block 1: one joint draw of b (J groups x p_re RE columns) ----
      // Conditions on observation-level y, design Z, and current hyperprior pl1.
      // Grouping factor levels are handled inside the block_* export.
      List block_i;
      if (is_gaussian) {
        block_i = block_rNormalReg_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3, 2
        );
      } else {
        block_i = block_rNormalGLM_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3,
          family, link, Gridtype, n_envopt,
          use_parallel, use_opencl, verbose
        );
      }
      
      b_i = Rcpp::as<NumericMatrix>(block_i["coefficients"]);
      if (b_i.nrow() != J || b_i.ncol() != p_re) {
        Rcpp::stop(
          "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
          b_i.nrow(), b_i.ncol(), J, p_re
        );
      }
      if (!have_ids) {
        List bi_info = block_i["block_info"];
        group_ids = Rcpp::as<CharacterVector>(bi_info["ids"]);
        have_ids = true;
      }
      
      // ---- Block 2: one hyperparameter draw per RE column j ----
      // Pseudo-response y_j = b_i[, j] (length J, one value per group level).
      // x_hyper[[j]] is the level-2 design for RE column j.
      for (int j = 0; j < p_re; ++j) {
        const NumericMatrix& X_j = mu_builder.X[j];
        if (X_j.nrow() != b_i.nrow()) {
          Rcpp::stop(
            "nrow(x_hyper[[%d]]) (%d) must equal nrow(y) (%d).",
            j + 1, X_j.nrow(), b_i.nrow()
          );
        }
        NumericVector y_j = b_i(Rcpp::_, j);
        NumericVector offset_j(X_j.nrow(), 0.0);
        NumericVector wt_j(X_j.nrow(), 1.0);
        
        if (!pr2[j].is_ing) {
          // dNormal Block 2: conjugate draw of gamma_k at fixed tau2[j].
          List out_j = rNormalReg(
            1, y_j, X_j, pr2[j].mu, pr2[j].P, offset_j, wt_j,
            pr2[j].dispersion, f2_gauss, f3_gauss, pr2[j].mu,
            "gaussian", "identity", 2
          );
          NumericMatrix coef_j =
            Rcpp::as<NumericMatrix>(out_j["coefficients"]);
          two_block_assign_fixef_component(
            fixef, j, NumericVector(coef_j(0, Rcpp::_)), fixef_start_v[j]
          );
          iters_draws(i, j) += 1.0;
        } else {
          // ING Block 2: joint draw of (gamma_k, tau2_k); tau2[j] feeds next
          // sweep's Block 1 prior row via the P1 refresh above.
          List out_j = rIndepNormalGammaReg(
            1, y_j, X_j, pr2[j].mu, pr2[j].P, offset_j, wt_j,
            pr2[j].shape, pr2[j].rate, pr2[j].max_disp_perc,
            Rcpp::Nullable<NumericVector>(pr2[j].disp_lower),
            Rcpp::Nullable<NumericVector>(pr2[j].disp_upper),
            2 /*Gridtype*/, 1 /*n_envopt*/,
            false /*use_parallel: n = 1 is serial anyway*/,
            false /*use_opencl*/, false /*verbose*/, false /*progbar*/
          );
          NumericMatrix coef_j = Rcpp::as<NumericMatrix>(out_j["out"]);
          NumericVector gamma_j(coef_j.nrow());
          for (int c = 0; c < coef_j.nrow(); ++c) {
            gamma_j[c] = coef_j(c, 0);
          }
          two_block_assign_fixef_component(
            fixef, j, gamma_j, fixef_start_v[j]
          );
          NumericVector disp_j = Rcpp::as<NumericVector>(out_j["disp_out"]);
          tau2[j] = disp_j[0];
          NumericVector iters_j =
            Rcpp::as<NumericVector>(out_j["iters_out"]);
          iters_draws(i, j) += iters_j[0];
        }
      }
      
      
      two_block_pack_fixef_row(fixef_temp, i, p_re, q_k, fixef);
      for (int j = 0; j < p_re; ++j) {
        tau2_temp(i, j) = tau2[j];
      }

      if (diag_sweeps) {
        two_block_accumulate_fixef_sum(fixef_sum, p_re, q_k, fixef);
        for (int g = 0; g < J; ++g) {
          for (int j = 0; j < p_re; ++j) {
            b_sum(g, j) += b_i(g, j);
          }
        }
        if (i == n - 1) {
          two_block_print_sweep_diagnostic(
            stage_label, m, m_convergence, n, p_re, p_dim, J,
            fixef_sum, b_sum, fixef_mode_v, b_mode_mat, have_b_mode,
            fixef_mode_list_for_print, re_names, group_levels
          );
        }
      }

      if (m == m_convergence - 1) {
        pack_two_block_chain_draw(
          i, p_re, J, fixef, tau2, b_i,
          fixef_draws, b_arr, disp_draws
        );
      }
    } // Enf of sweep loop
  }  /// End of chain loop
  
  if (chain_progbar || inner_progbar) {
    glmbayes::progress::progress_bar_finish();
  }
  
  List fixef_last(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_last[j] = fixef[j];
  }
  
  return List::create(
    Rcpp::Named("fixef_draws") = fixef_draws,
    Rcpp::Named("b_draws") = b_arr,
    Rcpp::Named("fixef_last") = fixef_last,
    Rcpp::Named("b_last") = b_i,
    Rcpp::Named("mu_all_last") = mu_all,
    Rcpp::Named("group_ids") = group_ids,
    Rcpp::Named("dispersion_fixef_draws") = disp_draws,
    Rcpp::Named("iters_fixef_draws") = iters_draws,
    Rcpp::Named("any_ing") = any_ing
  );
}


List two_block_rNormal_reg_v5_cpp_export(
    int n,
    int m_convergence,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    const List& x_hyper,
    const List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const List& pfamily_list,
    const List& fixef_start,
    const CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Function& f2,
    const Function& f3,
    const Function& f2_gauss,
    const Function& f3_gauss,
    const NumericVector& offset,
    const NumericVector& wt,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    Rcpp::Nullable<int> seed,
    int seed_offset,
    bool progbar,
    std::string stage_label,
    bool diag_sweeps,
    SEXP fixef_mode,
    SEXP b_mode
) {
  // Two-block Gibbs sampler (v5): independent short chains in C++.
  //
  // Indexing:
  //   p_re  = ncol(x) = number of random-effect *columns* in Z (e.g. intercept,
  //           slope RE terms), not the number of hyperparameters.
  //   J     = number of grouping-factor *levels* (rows of b_i).
  //   n     = number of independent chains (stored draws / replicates).
  //   m     = inner Gibbs sweeps per chain (m_convergence); only the state
  //           after the last sweep is stored for each chain.
  //
  // Loop hierarchy (outer -> inner):
  //   for m in 0..m_convergence-1   one full two-block sweep across all chains
  //     for i in 0..n-1     independent chains (sequential draws from R RNG)
  //       Block 1 once    joint draw of all group-level b (J x p_re matrix)
  //       for j in 0..p_re-1   Block 2 per RE column (hyperparameter gamma_k)
  //
  // Block 1 (level-1 / observation model): given hyperparameters fixef, draw
  //   random effects b_{g,k} for every group g and RE column k from the
  //   conditional of y | b (Gaussian reg or GLM blocked by the factor).
  //   block_rNormalReg / block_rNormalGLM handle grouping internally; there
  //   is no explicit loop over factor levels here.
  //
  // Block 2 (level-2 / hyperparameter model): given Block-1 draws b[,k] as
  //   pseudo-response (length J), draw hyperparameters gamma_k via regression
  //   on x_hyper[[k]] (one pfamily component per RE column j).
  
  if (n < 1) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (m_convergence < 1) {
    Rcpp::stop("'m_convergence' must be at least 1.");
  }

  const int p_re = x.ncol();
  const int J = group_levels.size();
  if (x_hyper.size() != p_re) {
    Rcpp::stop("length(x_hyper) must equal ncol(x) = %d.", p_re);
  }
  if (pfamily_list.size() != p_re) {
    Rcpp::stop("length(pfamily_list) must equal ncol(x) = %d.", p_re);
  }
  if (fixef_start.size() != p_re) {
    Rcpp::stop("length(fixef_start) must equal ncol(x) = %d.", p_re);
  }
  
  const bool is_gaussian = (family == "gaussian");
  
  MuAllBuilder mu_builder(x_hyper, group_levels);
  
  // Setup (no RNG): parse Block 2 pfamily priors, one entry per RE column j.
  std::vector<Block2PriorV2> pr2(p_re);
  bool any_ing = false;
  for (int j = 0; j < p_re; ++j) {
    const NumericMatrix& X_j = mu_builder.X[j];
    pr2[j] = block2_prior_prep_v2(List(pfamily_list[j]), j + 1, X_j.ncol());
    if (pr2[j].is_ing) any_ing = true;
  }
  
  // Deep snapshot: each chain resets from this copy (same as a fresh n = 1
  // v2 call with the current R RNG stream).
  std::vector<NumericVector> fixef_start_v(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_start_v[j] =
      Rcpp::clone(Rcpp::as<NumericVector>(fixef_start[j]));
  }

  std::vector<NumericVector> fixef_mode_v(p_re);
  if (diag_sweeps && !Rf_isNull(fixef_mode)) {
    List fixef_mode_in = Rcpp::as<List>(fixef_mode);
    if (fixef_mode_in.size() != p_re) {
      Rcpp::stop("length(fixef_mode) must equal ncol(x) = %d.", p_re);
    }
    for (int j = 0; j < p_re; ++j) {
      fixef_mode_v[j] =
        Rcpp::clone(Rcpp::as<NumericVector>(fixef_mode_in[j]));
    }
  } else {
    fixef_mode_v = fixef_start_v;
  }

  List fixef_mode_list_for_print;
  if (diag_sweeps && !Rf_isNull(fixef_mode)) {
    fixef_mode_list_for_print = Rcpp::as<List>(fixef_mode);
  } else {
    fixef_mode_list_for_print = fixef_start;
  }

  NumericMatrix b_mode_mat;
  bool have_b_mode = false;
  if (diag_sweeps && !Rf_isNull(b_mode)) {
    b_mode_mat = Rcpp::as<NumericMatrix>(b_mode);
    if (b_mode_mat.nrow() == J && b_mode_mat.ncol() == p_re) {
      have_b_mode = true;
    }
  }

  CharacterVector re_names;
  if (x.hasAttribute("dimnames")) {
    List dn = x.attr("dimnames");
    if (!Rf_isNull(dn[1])) {
      re_names = dn[1];
    }
  }
  if (re_names.size() != p_re) {
    re_names = CharacterVector(p_re);
    for (int j = 0; j < p_re; ++j) {
      re_names[j] = "RE" + std::to_string(j + 1);
    }
  }

  std::vector<NumericVector> fixef(p_re);
  
  // Working state for Block 2 hyperparameters (gamma_k) and ING dispersions.
  // fixef[j] = hyperparameter vector for RE column j; tau2[j] = Block 2
  // dispersion for that component (fixed for dNormal, updated for ING).
  std::vector<double> tau2_start(p_re);
  for (int j = 0; j < p_re; ++j) {
    tau2_start[j] = pr2[j].dispersion;
  }
  std::vector<double> tau2 = tau2_start;
  
  // Block 1 prior template for ING: p_re x p_re precision; ING rows refreshed
  // each sweep from current tau2 (see inner loop over j below).
  NumericMatrix base_P1;
  if (any_ing) {
    if (has_non_null(prior_list_block1, "P")) {
      base_P1 = Rcpp::as<NumericMatrix>(prior_list_block1["P"]);
    } else if (has_non_null(prior_list_block1, "Sigma")) {
      NumericMatrix S1 = Rcpp::as<NumericMatrix>(prior_list_block1["Sigma"]);
      arma::mat S1a(const_cast<double*>(S1.begin()),
                    S1.nrow(), S1.ncol(), false);
      arma::mat P1inv = arma::inv_sympd(S1a);
      base_P1 = NumericMatrix(Rcpp::wrap(0.5 * (P1inv + P1inv.t())));
    } else {
      Rcpp::stop("prior_list_block1 must contain 'P' or 'Sigma'.");
    }
    if (base_P1.nrow() != p_re || base_P1.ncol() != p_re) {
      Rcpp::stop(
        "prior_list_block1 P/Sigma must be %d x %d.", p_re, p_re
      );
    }
  }
  
  // Output buffers:
  //   fixef_draws[[j]]     n x q_k  Block 2 hyperparameters per chain
  //   b_arr                J x p_re x n  Block 1 random effects per chain
  //   disp_draws, iters_draws   n x p_re  Block 2 tau^2 and sampler counts
  List fixef_draws(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_draws[j] = NumericMatrix(n, fixef_start_v[j].size());
  }
  NumericVector b_arr(Rcpp::Dimension(J, p_re, n));
  NumericMatrix disp_draws(n, p_re);
  NumericMatrix iters_draws(n, p_re);
  NumericVector iters_ranef(n);
  std::fill(iters_ranef.begin(), iters_ranef.end(), 0.0);
  List sweep_stats(m_convergence);
  
  NumericMatrix mu_all;
  NumericMatrix b_i;
  CharacterVector group_ids;
  bool have_ids = false;
  
  const bool have_seed = seed.isNotNull();
  
  
  int p_dim = 0;
  std::vector<int> q_k(p_re);
  for (int j = 0; j < p_re; ++j) {
    q_k[j] = fixef_start_v[j].size();
    p_dim += q_k[j];
  }

  NumericVector fixef_sum;
  NumericMatrix b_sum;
  if (diag_sweeps) {
    fixef_sum = NumericVector(p_dim, 0.0);
    b_sum = NumericMatrix(J, p_re);
    std::fill(b_sum.begin(), b_sum.end(), 0.0);
    if (stage_label.empty()) {
      stage_label = "v5";
    }
  }
  
  NumericMatrix fixef_temp=NumericMatrix(n,p_dim);
  NumericMatrix tau2_temp = NumericMatrix(n, p_re);
  
  
  for (int i = 0; i < n; ++i) {
    two_block_pack_fixef_row(fixef_temp, i, p_re, q_k, fixef_start_v);
    for (int j = 0; j < p_re; ++j) {
      tau2_temp(i, j) = tau2_start[j];
    }
  }

  // Live sweep diagnostics disabled; use progbar only (lmerb / v6 style).
  // if (diag_sweeps) { ... entry fixef_start print ... }

  two_block_driver_banner_v5(n, m_convergence, seed_offset, have_seed, p_dim);

  const bool chain_progbar = progbar && n > 1;
  const bool sweep_progbar = progbar && n <= 1;
  NumericVector b_work(Rcpp::Dimension(J, p_re, n));

  // ---- Outer loop: m_convergence full two-block Gibbs sweeps ----
  for (int m = 0; m < m_convergence; ++m) {
    Rcpp::checkUserInterrupt();

    const int sweep1 = m + 1;
    const std::string prefix_re =
      two_block_progbar_prefix(stage_label, sweep1, m_convergence, "Block1");
    const std::string prefix_fe =
      two_block_progbar_prefix(stage_label, sweep1, m_convergence, "Block2");
    const std::string prefix_sweep =
      two_block_sweep_only_prefix(stage_label, sweep1, m_convergence);

    // ---- Block 1 (RE): all chains ----
    for (int i = 0; i < n; ++i) {
      Rcpp::checkUserInterrupt();

      if (chain_progbar) {
        glmbayes::progress::progress_bar(
          static_cast<double>(i + 1), static_cast<double>(n), prefix_re
        );
      }

      two_block_unpack_fixef_row(
        fixef_temp, i, p_re, q_k, fixef, fixef_start_v
      );
      for (int j = 0; j < p_re; ++j) {
        tau2[j] = tau2_temp(i, j);
      }

      mu_all = mu_builder.build(fixef);
      List pl1;
      if (any_ing) {
        NumericMatrix P1 = Rcpp::clone(base_P1);
        for (int j = 0; j < p_re; ++j) {
          if (!pr2[j].is_ing) continue;
          for (int c = 0; c < p_re; ++c) {
            P1(j, c) = 0.0;
            P1(c, j) = 0.0;
          }
          P1(j, j) = 1.0 / tau2[j];
        }
        List pl1_base = List::create(Rcpp::Named("P") = P1);
        pl1 = block1_prior_list(
          mu_all, pl1_base, dispersion_block1, ddef_block1
        );
      } else {
        pl1 = block1_prior_list(
          mu_all, prior_list_block1, dispersion_block1, ddef_block1
        );
      }

      List block_i;
      if (is_gaussian) {
        block_i = block_rNormalReg_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3, 2
        );
      } else {
        block_i = block_rNormalGLM_cpp_export(
          1, y, x, block, pl1, R_NilValue, offset, wt, f2, f3,
          family, link, Gridtype, n_envopt,
          use_parallel, use_opencl, verbose
        );
      }

      b_i = Rcpp::as<NumericMatrix>(block_i["coefficients"]);
      if (b_i.nrow() != J || b_i.ncol() != p_re) {
        Rcpp::stop(
          "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
          b_i.nrow(), b_i.ncol(), J, p_re
        );
      }
      {
        List bi_info = block_i["block_info"];
        CharacterVector block_ids =
          Rcpp::as<CharacterVector>(bi_info["ids"]);
        if (!have_ids) {
          group_ids = block_ids;
          have_ids = true;
        }
        two_block_reorder_b_to_group_levels(b_i, block_ids, group_levels);
      }
      iters_ranef[i] += two_block_block1_iters_mean(block_i);
      store_b_chain(b_work, i, J, p_re, b_i);
    }

    if (chain_progbar) {
      glmbayes::progress::progress_bar_finish(false);
    }

    // ---- Block 2 (fixef): all chains ----
    for (int i = 0; i < n; ++i) {
      Rcpp::checkUserInterrupt();

      if (chain_progbar) {
        glmbayes::progress::progress_bar(
          static_cast<double>(i + 1), static_cast<double>(n), prefix_fe
        );
      }

      two_block_unpack_fixef_row(
        fixef_temp, i, p_re, q_k, fixef, fixef_start_v
      );
      for (int j = 0; j < p_re; ++j) {
        tau2[j] = tau2_temp(i, j);
      }
      load_b_chain(b_work, i, J, p_re, b_i);

      for (int j = 0; j < p_re; ++j) {
        const NumericMatrix& X_j = mu_builder.X[j];
        if (X_j.nrow() != b_i.nrow()) {
          Rcpp::stop(
            "nrow(x_hyper[[%d]]) (%d) must equal nrow(y) (%d).",
            j + 1, X_j.nrow(), b_i.nrow()
          );
        }
        double it_j = 0.0;
        NumericVector gamma_j = two_block_block2_rglmb_gamma(
          b_i,
          j,
          X_j,
          group_levels,
          List(pfamily_list[j]),
          pr2[j].is_ing,
          tau2[j],
          it_j
        );
        two_block_assign_fixef_component(
          fixef, j, gamma_j, fixef_start_v[j]
        );
        iters_draws(i, j) += it_j;
      }

      two_block_pack_fixef_row(fixef_temp, i, p_re, q_k, fixef);
      for (int j = 0; j < p_re; ++j) {
        tau2_temp(i, j) = tau2[j];
      }

      if (m == m_convergence - 1) {
        pack_two_block_chain_draw(
          i, p_re, J, fixef, tau2, b_i,
          fixef_draws, b_arr, disp_draws
        );
      }
    }

    if (chain_progbar) {
      glmbayes::progress::progress_bar_finish(m == m_convergence - 1);
    }

    if (sweep_progbar) {
      glmbayes::progress::progress_bar(
        static_cast<double>(sweep1), static_cast<double>(m_convergence),
        prefix_sweep
      );
      glmbayes::progress::progress_bar_finish(m == m_convergence - 1);
    }

    sweep_stats[m] = two_block_snapshot_fixef_stats_cpp(
      fixef_temp, n, p_re, q_k, fixef_start_v, re_names
    );

    if (m == 0) {
      two_block_print_sweep_chain_means(
        stage_label, sweep1, n, J, p_re, q_k,
        fixef_temp, b_work, fixef_start_v, re_names
      );
    }

    // Live per-sweep diagnostics disabled (defer to future sweep_history).
    // if (diag_sweeps) { two_block_print_sweep_diagnostic(...); }
  } // End of sweep loop
  
  List fixef_last(p_re);
  for (int j = 0; j < p_re; ++j) {
    fixef_last[j] = fixef[j];
  }
  
  return List::create(
    Rcpp::Named("fixef_draws") = fixef_draws,
    Rcpp::Named("b_draws") = b_arr,
    Rcpp::Named("fixef_last") = fixef_last,
    Rcpp::Named("b_last") = b_i,
    Rcpp::Named("mu_all_last") = mu_all,
    Rcpp::Named("group_ids") = group_ids,
    Rcpp::Named("dispersion_fixef_draws") = disp_draws,
    Rcpp::Named("iters_fixef_draws") = iters_draws,
    Rcpp::Named("iters_ranef_draws") = iters_ranef,
    Rcpp::Named("sweep_stats") = sweep_stats,
    Rcpp::Named("any_ing") = any_ing
  );
}



} // namespace sim
} // namespace glmbayes
