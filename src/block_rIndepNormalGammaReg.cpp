// block_rIndepNormalGammaReg.cpp
// Block envelope ING sampler: BlockEnvelopeCentering, BlockEnvelopeBuild,
// BlockEnvelopeSim, and rIndepNormalGammaRegBlock orchestrator.
// New code only — does not modify existing sampler implementations.

#include "RcppArmadillo.h"
#include "Envelopefuncs.h"
#include "simfuncs.h"
#include <cmath>
#include <string>

namespace glmbayes {
namespace env {

namespace {

using Rcpp::CharacterVector;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::Nullable;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::RObject;

static inline bool has_non_null(const List& pl, const char* name) {
  return pl.containsElementNamed(name) && !Rf_isNull(pl[name]);
}

static inline bool is_symmetric_mat(const NumericMatrix& M, double tol = 1e-10) {
  if (M.nrow() != M.ncol()) return false;
  for (int i = 0; i < M.nrow(); ++i) {
    for (int j = i + 1; j < M.ncol(); ++j) {
      if (std::fabs(M(i, j) - M(j, i)) > tol) return false;
    }
  }
  return true;
}

static void check_P_pd(const NumericMatrix& P, const char* what) {
  arma::mat P_arma(const_cast<double*>(P.begin()), P.nrow(), P.ncol(), false);
  try {
    arma::mat L = arma::chol(P_arma);
    (void)L;
  } catch (...) {
    Rcpp::stop("%s must be positive definite.", what);
  }
}

List prior_list_to_P_Sigma(List pl) {
  if (!has_non_null(pl, "mu")) {
    Rcpp::stop("prior_list must contain 'mu'.");
  }
  NumericVector mu = pl["mu"];
  List out = List::create(Rcpp::Named("mu") = mu);

  if (has_non_null(pl, "P")) {
    NumericMatrix P = pl["P"];
    if (!is_symmetric_mat(P)) {
      Rcpp::stop("prior precision matrix P must be symmetric.");
    }
    check_P_pd(P, "P");
    NumericMatrix Sigma = Rcpp::wrap(arma::inv(arma::mat(
      const_cast<double*>(P.begin()), P.nrow(), P.ncol(), false
    )));
    out["P"] = P;
    out["Sigma"] = Sigma;
    return out;
  }
  if (has_non_null(pl, "Sigma")) {
    NumericMatrix Sigma = pl["Sigma"];
    if (!is_symmetric_mat(Sigma)) {
      Rcpp::stop("prior covariance Sigma must be symmetric.");
    }
    arma::mat S(const_cast<double*>(Sigma.begin()), Sigma.nrow(), Sigma.ncol(), false);
    arma::mat P = arma::inv_sympd(S);
    out["Sigma"] = Sigma;
    out["P"] = NumericMatrix(Rcpp::wrap(0.5 * (P + P.t())));
    return out;
  }
  Rcpp::stop("prior_list must contain 'P' or 'Sigma'.");
}

List base_prior_block(List pl, int l1) {
  List ps = prior_list_to_P_Sigma(pl);
  NumericVector mu = ps["mu"];
  if (mu.size() != l1) {
    Rcpp::stop("length(mu) must equal ncol(x) (%d).", l1);
  }
  NumericMatrix P = ps["P"];
  if (P.nrow() != l1 || P.ncol() != l1) {
    Rcpp::stop("dim(P) or dim(Sigma) must be %d x %d.", l1, l1);
  }
  return List::create(
    Rcpp::Named("mu") = mu,
    Rcpp::Named("Sigma") = ps["Sigma"],
    Rcpp::Named("P") = P
  );
}

List rep_list_blocks(const List& one, int k) {
  List out(k);
  for (int j = 0; j < k; ++j) {
    out[j] = Rcpp::clone(one);
  }
  return out;
}

List normalize_prior_for_block_ing(
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    const List& block_info,
    int l1
) {
  const int k = block_info["k"];

  if (!Rf_isNull(prior_lists_sexp)) {
    List prior_lists(prior_lists_sexp);
    if (!Rcpp::is<List>(prior_lists)) {
      Rcpp::stop("'prior_lists' must be a list.");
    }
    if (prior_lists.size() == 1) {
      List one = base_prior_block(List(prior_lists[0]), l1);
      return rep_list_blocks(one, k);
    }
    if (prior_lists.size() != static_cast<R_xlen_t>(k)) {
      Rcpp::stop("'prior_lists' must have length 1 or k = %d.", k);
    }
    List out(k);
    for (int j = 0; j < k; ++j) {
      out[j] = base_prior_block(List(prior_lists[j]), l1);
    }
    return out;
  }

  if (Rf_isNull(prior_list_sexp)) {
    Rcpp::stop("Provide 'prior_list' or 'prior_lists'.");
  }

  List prior_list(prior_list_sexp);

  if (has_non_null(prior_list, "dispersion")) {
    Rcpp::stop(
      "BlockEnvelopeCentering requires unknown measurement dispersion; "
      "prior_list must not contain a fixed 'dispersion'."
    );
  }

  SEXP mu_sexp = has_non_null(prior_list, "mu") ? (SEXP)prior_list["mu"] : R_NilValue;
  if (Rcpp::is<NumericMatrix>(mu_sexp)) {
    NumericMatrix mu_mat(mu_sexp);
    if (mu_mat.nrow() != l1) {
      Rcpp::stop("nrow(prior_list$mu) must equal ncol(x) (%d).", l1);
    }
    if (mu_mat.ncol() == 1) {
      List pl_copy = Rcpp::clone(prior_list);
      pl_copy["mu"] = mu_mat(Rcpp::_, 0);
      List one = base_prior_block(pl_copy, l1);
      return rep_list_blocks(one, k);
    }
    if (mu_mat.ncol() != k) {
      Rcpp::stop("ncol(prior_list$mu) must equal number of blocks k = %d.", k);
    }

    List P_list;
    List Sigma_list;
    bool have_P_list = false;
    bool have_S_list = false;

    if (has_non_null(prior_list, "P")) {
      SEXP P_sexp = prior_list["P"];
      if (Rcpp::is<List>(P_sexp)) {
        P_list = List(P_sexp);
        have_P_list = P_list.size() == 1 || P_list.size() == static_cast<R_xlen_t>(k);
      } else if (Rcpp::is<NumericMatrix>(P_sexp)) {
        P_list = List::create(P_sexp);
        have_P_list = true;
      }
    }
    if (has_non_null(prior_list, "Sigma")) {
      SEXP S_sexp = prior_list["Sigma"];
      if (Rcpp::is<List>(S_sexp)) {
        Sigma_list = List(S_sexp);
        have_S_list = Sigma_list.size() == 1 || Sigma_list.size() == static_cast<R_xlen_t>(k);
      } else if (Rcpp::is<NumericMatrix>(S_sexp)) {
        Sigma_list = List::create(S_sexp);
        have_S_list = true;
      }
    }

    List out(k);
    for (int j = 0; j < k; ++j) {
      List pl_j = List::create(Rcpp::Named("mu") = mu_mat(Rcpp::_, j));
      if (have_P_list) {
        pl_j["P"] = P_list[std::min(j, static_cast<int>(P_list.size()) - 1)];
      } else if (have_S_list) {
        pl_j["Sigma"] = Sigma_list[std::min(j, static_cast<int>(Sigma_list.size()) - 1)];
      } else {
        Rcpp::stop("prior_list must contain 'P' or 'Sigma'.");
      }
      out[j] = base_prior_block(pl_j, l1);
    }
    return out;
  }

  List one = base_prior_block(prior_list, l1);
  return rep_list_blocks(one, k);
}

inline int check_index_1based(int idx, int n, const char* what) {
  if (idx < 1 || idx > n) {
    Rcpp::stop("%s index %d out of range [1, %d]", what, idx, n);
  }
  return idx - 1;
}

NumericVector slice_numeric(const NumericVector& v, const IntegerVector& rows) {
  const int m = rows.size();
  NumericVector out(m);
  for (int i = 0; i < m; ++i) {
    out[i] = v[check_index_1based(rows[i], v.size(), "row")];
  }
  return out;
}

NumericMatrix slice_matrix_rows(const NumericMatrix& x, const IntegerVector& rows) {
  const int m = rows.size();
  const int l1 = x.ncol();
  NumericMatrix out(m, l1);
  for (int i = 0; i < m; ++i) {
    const int r = check_index_1based(rows[i], x.nrow(), "row");
    for (int j = 0; j < l1; ++j) {
      out(i, j) = x(r, j);
    }
  }
  return out;
}

double sum_wt(const NumericVector& wt) {
  double s = 0.0;
  for (int i = 0; i < wt.size(); ++i) s += wt[i];
  return s;
}

struct BlockWlsInit {
  int rank;
  double rss_wls;
  double n_w_j;
  double dispersion_wls;
  bool identifiable;
};

BlockWlsInit block_wls_init(
    const NumericVector& y_j,
    const NumericMatrix& x_j,
    const NumericVector& offset_j,
    const NumericVector& wt_j
) {
  BlockWlsInit out;
  out.n_w_j = sum_wt(wt_j);

  const int n_j = y_j.size();
  NumericVector ystar(n_j);
  for (int i = 0; i < n_j; ++i) {
    ystar[i] = y_j[i] - offset_j[i];
  }

  Rcpp::Function lm_wfit("lm.wfit");
  List fit = lm_wfit(
    Rcpp::_["x"] = x_j,
    Rcpp::_["y"] = ystar,
    Rcpp::_["w"] = wt_j
  );

  NumericVector res = fit["residuals"];
  out.rss_wls = 0.0;
  for (int i = 0; i < res.size(); ++i) {
    out.rss_wls += res[i] * res[i];
  }
  out.rank = Rcpp::as<int>(fit["rank"]);
  out.identifiable = out.rank >= 1;
  const double denom = std::max(out.n_w_j - static_cast<double>(out.rank), 1.0);
  out.dispersion_wls = out.rss_wls / denom;
  return out;
}

List block_envelope_centering_one(
    const NumericVector& y_j,
    const NumericMatrix& x_j,
    const NumericVector& mu_j,
    const NumericMatrix& P_j,
    const NumericVector& offset_j,
    const NumericVector& wt_j,
    double dispersion2,
    int rank_j,
    bool identifiable
) {
  const int l1 = x_j.ncol();
  const int l2_j = x_j.nrow();
  const double n_w_j = sum_wt(wt_j);

  if (!identifiable || rank_j < 1) {
    return List::create(
      Rcpp::Named("RSS_post") = 0.0,
      Rcpp::Named("b_post_mean") = mu_j,
      Rcpp::Named("rank") = rank_j,
      Rcpp::Named("n_w_j") = n_w_j,
      Rcpp::Named("identifiable") = false
    );
  }

  const arma::mat X = Rcpp::as<arma::mat>(x_j);
  const arma::vec Y = Rcpp::as<arma::vec>(y_j);
  const arma::vec off = Rcpp::as<arma::vec>(offset_j);
  const arma::vec wv = Rcpp::as<arma::vec>(wt_j);
  const arma::vec mu_vec = Rcpp::as<arma::vec>(mu_j);
  const arma::mat P_arma = Rcpp::as<arma::mat>(P_j);

  arma::mat Xw(l2_j, l1);
  arma::vec yw(l2_j);
  for (int i = 0; i < l2_j; ++i) {
    const double sw = std::sqrt(wv[i]);
    Xw.row(i) = sw * X.row(i);
    yw[i] = (Y[i] - off[i]) * sw;
  }

  const arma::mat RA = arma::chol(P_arma);
  const arma::vec z_bot = RA * mu_vec;
  const arma::mat XtWX = X.t() * (arma::diagmat(wv) * X);

  const double s = 1.0 / std::sqrt(dispersion2);
  arma::mat W(l2_j + l1, l1);
  arma::vec z(l2_j + l1);
  W.rows(0, l2_j - 1) = s * Xw;
  W.rows(l2_j, l2_j + l1 - 1) = RA;
  z.rows(0, l2_j - 1) = s * yw;
  z.rows(l2_j, l2_j + l1 - 1) = z_bot;

  const arma::mat WtW = W.t() * W;
  const arma::mat IR = arma::inv(arma::trimatu(arma::chol(WtW)));
  const arma::mat Sigma = IR * arma::trans(IR);
  const arma::vec b2_fast = Sigma * (W.t() * z);

  const arma::vec r_fast = Y - X * b2_fast - off;
  const double rss_at_mean = arma::dot(wv, r_fast % r_fast);
  const double trace_term = arma::trace(XtWX * Sigma);
  const double RSS_post = rss_at_mean + trace_term;

  return List::create(
    Rcpp::Named("RSS_post") = RSS_post,
    Rcpp::Named("b_post_mean") = NumericVector(b2_fast.begin(), b2_fast.end()),
    Rcpp::Named("rank") = rank_j,
    Rcpp::Named("n_w_j") = n_w_j,
    Rcpp::Named("identifiable") = true
  );
}

Nullable<double> resolve_nullable_bound(
    Nullable<double> arg_val,
    const List& prior_list,
    const char* name
) {
  if (arg_val.isNotNull()) {
    return arg_val;
  }
  if (has_non_null(prior_list, name)) {
    NumericVector v = prior_list[name];
    if (v.size() >= 1) {
      return Nullable<double>(Rcpp::wrap(static_cast<double>(v[0])));
    }
  }
  return R_NilValue;
}

}  // anonymous namespace

List BlockEnvelopeCentering(
    NumericVector y,
    NumericMatrix x,
    SEXP block,
    SEXP prior_list_sexp,
    SEXP prior_lists_sexp,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Nullable<double> disp_lower,
    Nullable<double> disp_upper,
    int p_re,
    int n_rss_iter,
    bool verbose
) {
  const int l2 = y.size();
  const int l1 = x.ncol();
  if (x.nrow() != l2) {
    Rcpp::stop("nrow(x) must equal length(y).");
  }
  if (n_rss_iter < 1) {
    Rcpp::stop("'n_rss_iter' must be at least 1.");
  }

  if (offset.size() == 1) offset = Rcpp::rep(offset[0], l2);
  if (wt.size() == 1) wt = Rcpp::rep(wt[0], l2);
  if (offset.size() != l2) {
    Rcpp::stop("length(offset) must be 1 or length(y).");
  }
  if (wt.size() != l2) {
    Rcpp::stop("length(wt) must be 1 or length(y).");
  }

  List block_info = glmbayes::sim::normalize_block_cpp(block, l2);
  const int k = block_info["k"];
  CharacterVector ids = block_info["ids"];
  List row_blocks = block_info["rows"];

  List prior_block = normalize_prior_for_block_ing(
    prior_list_sexp, prior_lists_sexp, block_info, l1
  );

  List prior_list = Rf_isNull(prior_list_sexp) ? List() : List(prior_list_sexp);
  if (max_disp_perc <= 0.0 || max_disp_perc >= 1.0) {
    if (has_non_null(prior_list, "max_disp_perc")) {
      max_disp_perc = Rcpp::as<double>(prior_list["max_disp_perc"]);
    } else {
      max_disp_perc = 0.99;
    }
  }

  disp_lower = resolve_nullable_bound(disp_lower, prior_list, "disp_lower");
  disp_upper = resolve_nullable_bound(disp_upper, prior_list, "disp_upper");

  if (p_re < 1) {
    p_re = l1;
  }

  const double n_w = sum_wt(wt);
  const double n_prior_implied = 2.0 * shape - 1.0 - static_cast<double>(p_re);
  if (n_prior_implied > n_w) {
    Rcpp::stop(
      "dIndependent_Normal_Gamma prior implies n_prior = %g effective prior "
      "observations, but the data supply only n_w = sum(weights) = %g. The "
      "dispersion envelope requires n_prior <= n_w (prior weight pwt <= 0.5); "
      "weaken the prior (smaller shape) or supply more data.",
      n_prior_implied, n_w
    );
  }

  std::vector<BlockWlsInit> wls_init(k);
  double rss_wls_sum = 0.0;
  double denom_sum = 0.0;
  for (int j = 0; j < k; ++j) {
    IntegerVector rows = row_blocks[j];
    List pb = prior_block[j];
    NumericVector y_j = slice_numeric(y, rows);
    NumericMatrix x_j = slice_matrix_rows(x, rows);
    NumericVector offset_j = slice_numeric(offset, rows);
    NumericVector wt_j = slice_numeric(wt, rows);

    wls_init[j] = block_wls_init(y_j, x_j, offset_j, wt_j);
    rss_wls_sum += wls_init[j].rss_wls;
    denom_sum += std::max(wls_init[j].n_w_j - static_cast<double>(wls_init[j].rank), 1.0);

    if (verbose) {
      Rcpp::Rcout << "BlockEnvelopeCentering block " << (j + 1)
                  << " rank=" << wls_init[j].rank
                  << " n_obs=" << rows.size() << std::endl;
    }
  }

  double dispersion2 = rss_wls_sum / denom_sum;
  double RSS_post_pooled = NA_REAL;

  List blocks(k);
  for (int iter = 0; iter < n_rss_iter; ++iter) {
    RSS_post_pooled = 0.0;
    for (int j = 0; j < k; ++j) {
      IntegerVector rows = row_blocks[j];
      List pb = prior_block[j];
      NumericVector mu_j = pb["mu"];
      NumericMatrix P_j = pb["P"];
      NumericVector y_j = slice_numeric(y, rows);
      NumericMatrix x_j = slice_matrix_rows(x, rows);
      NumericVector offset_j = slice_numeric(offset, rows);
      NumericVector wt_j = slice_numeric(wt, rows);

      List one = block_envelope_centering_one(
        y_j, x_j, mu_j, P_j, offset_j, wt_j,
        dispersion2,
        wls_init[j].rank,
        wls_init[j].identifiable
      );

      const double RSS_j = Rcpp::as<double>(one["RSS_post"]);
      RSS_post_pooled += RSS_j;

      if (iter == n_rss_iter - 1) {
        blocks[j] = List::create(
          Rcpp::Named("id") = Rcpp::wrap(Rcpp::as<std::string>(ids[j])),
          Rcpp::Named("n_obs") = rows.size(),
          Rcpp::Named("n_w_j") = one["n_w_j"],
          Rcpp::Named("rank") = one["rank"],
          Rcpp::Named("identifiable") = one["identifiable"],
          Rcpp::Named("RSS_post") = RSS_j,
          Rcpp::Named("mu") = mu_j,
          Rcpp::Named("b_post_mean") = one["b_post_mean"],
          Rcpp::Named("dispersion_wls") = wls_init[j].dispersion_wls
        );
      }
    }

    const double shape2 = shape + n_w / 2.0;
    const double rate2 = rate + RSS_post_pooled / 2.0;
    dispersion2 = rate2 / (shape2 - 1.0);
  }

  const double shape2 = shape + n_w / 2.0;
  const double rate3 = rate + RSS_post_pooled / 2.0;

  RObject disp_lower_out = disp_lower.isNotNull()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_lower.get())))
    : RObject(R_NilValue);
  RObject disp_upper_out = disp_upper.isNotNull()
    ? RObject(Rcpp::wrap(Rcpp::as<double>(disp_upper.get())))
    : RObject(R_NilValue);

  return List::create(
    Rcpp::Named("dispersion") = dispersion2,
    Rcpp::Named("RSS_post") = RSS_post_pooled,
    Rcpp::Named("n_w") = n_w,
    Rcpp::Named("p_re") = p_re,
    Rcpp::Named("l1") = l1,
    Rcpp::Named("l2") = l2,
    Rcpp::Named("k") = k,
    Rcpp::Named("shape2") = shape2,
    Rcpp::Named("rate3") = rate3,
    Rcpp::Named("max_disp_perc") = max_disp_perc,
    Rcpp::Named("disp_lower") = disp_lower_out,
    Rcpp::Named("disp_upper") = disp_upper_out,
    Rcpp::Named("block_info") = block_info,
    Rcpp::Named("blocks") = blocks,
    Rcpp::Named("prior_lists") = prior_block
  );
}

}  // namespace env
}  // namespace glmbayes
