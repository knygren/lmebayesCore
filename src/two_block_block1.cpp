// two_block_block1.cpp
// Block 1 batch helpers for run_sweep_outer_chains_v6 (R-oracle semantics).
// Migrated incrementally from the R driver in R/two_block_batch_gibbs.R.

#include "simfuncs.h"
#include "progress_utils.h"

#include <string>
#include <vector>

namespace glmbayes {
namespace sim {

namespace {

void set_matrix_dimnames(
    Rcpp::NumericMatrix& M,
    const Rcpp::CharacterVector& rownames,
    SEXP colnames_src
) {
  if (rownames.size() == 0 && Rf_isNull(colnames_src)) {
    return;
  }
  Rcpp::List dn(2);
  if (rownames.size() > 0) {
    dn[0] = rownames;
  } else {
    dn[0] = R_NilValue;
  }
  if (Rf_isNull(colnames_src)) {
    dn[1] = R_NilValue;
  } else {
    dn[1] = colnames_src;
  }
  M.attr("dimnames") = dn;
}

std::vector<Rcpp::NumericVector> fixef_vec_from_list(
    const Rcpp::List& fixef,
    const Rcpp::CharacterVector& re_names
) {
  const int p_re = re_names.size();
  std::vector<Rcpp::NumericVector> out(static_cast<size_t>(p_re));
  for (int i = 0; i < p_re; ++i) {
    const std::string nm = Rcpp::as<std::string>(re_names[i]);
    if (!fixef.containsElementNamed(nm.c_str())) {
      Rcpp::stop("fixef must contain element '%s'.", nm.c_str());
    }
    out[static_cast<size_t>(i)] =
      Rcpp::as<Rcpp::NumericVector>(fixef[nm.c_str()]);
  }
  return out;
}

/// Forward optional list element (mirror R \code{$} on absent name → \code{NULL}).
SEXP optional_list_elt(const Rcpp::List& pl, const char* name) {
  return pl.containsElementNamed(name) ? pl[name] : R_NilValue;
}

/// Port of \code{.two_block_batch_fixef_chain(batch, i)}.
Rcpp::List fixef_list_from_batch_chain(
    const Rcpp::List& batch_fixef,
    int chain_i,
    const Rcpp::CharacterVector& re_names
) {
  if (chain_i < 1) {
    Rcpp::stop("'chain_i' must be at least 1.");
  }
  const int row_i = chain_i - 1;
  const int p_re = re_names.size();
  Rcpp::List out(p_re);
  out.attr("names") = re_names;
  for (int k = 0; k < p_re; ++k) {
    const std::string nm = Rcpp::as<std::string>(re_names[k]);
    if (!batch_fixef.containsElementNamed(nm.c_str())) {
      Rcpp::stop("batch_fixef must contain element '%s'.", nm.c_str());
    }
    Rcpp::NumericMatrix mat =
      Rcpp::as<Rcpp::NumericMatrix>(batch_fixef[nm.c_str()]);
    if (row_i >= mat.nrow()) {
      Rcpp::stop(
        "chain index %d exceeds nrow(batch_fixef[['%s']]) (%d).",
        chain_i, nm.c_str(), mat.nrow()
      );
    }
    Rcpp::NumericVector row = mat(row_i, Rcpp::_);
    SEXP col_dn = R_NilValue;
    SEXP mat_dn = mat.attr("dimnames");
    if (!Rf_isNull(mat_dn)) {
      Rcpp::List mat_dnl(mat_dn);
      if (mat_dnl.size() >= 2 && !Rf_isNull(mat_dnl[1])) {
        col_dn = mat_dnl[1];
      }
    }
    if (!Rf_isNull(col_dn)) {
      row.attr("names") = col_dn;
    }
    out[k] = row;
  }
  return out;
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

/// Mirror \code{block_rNormalGLM()} / \code{block_rNormalReg()} dimnames on
/// \code{coefficients} (\code{colnames(x)}, \code{rownames} from
/// \code{block_info$ids}).
void set_block_draw_coefficient_dimnames(
    Rcpp::NumericMatrix& coef_draw,
    const Rcpp::NumericMatrix& x,
    const Rcpp::List& block_info
) {
  Rcpp::List dn(2);
  SEXP old_dn = coef_draw.attr("dimnames");
  if (!Rf_isNull(old_dn)) {
    Rcpp::List old_dnl(old_dn);
    if (old_dnl.size() >= 1 && !Rf_isNull(old_dnl[0])) {
      dn[0] = old_dnl[0];
    }
    if (old_dnl.size() >= 2 && !Rf_isNull(old_dnl[1])) {
      dn[1] = old_dnl[1];
    }
  }

  SEXP x_dn = x.attr("dimnames");
  if (!Rf_isNull(x_dn)) {
    Rcpp::List x_dnl(x_dn);
    if (x_dnl.size() >= 2 && !Rf_isNull(x_dnl[1])) {
      dn[1] = x_dnl[1];
    }
  }
  if (block_info.containsElementNamed("ids")) {
    SEXP ids = block_info["ids"];
    if (!Rf_isNull(ids)) {
      dn[0] = ids;
    }
  }
  if (!Rf_isNull(dn[0]) || !Rf_isNull(dn[1])) {
    coef_draw.attr("dimnames") = dn;
  }
}

/// Row ids for reorder (default \code{block_ids = rownames(b_draw)} in R).
SEXP matrix_rownames_sexp(const Rcpp::NumericMatrix& M) {
  SEXP dn = M.attr("dimnames");
  if (Rf_isNull(dn)) {
    return R_NilValue;
  }
  Rcpp::List dnl(dn);
  if (dnl.size() < 1 || Rf_isNull(dnl[0])) {
    return R_NilValue;
  }
  return dnl[0];
}

} // namespace

MuAllBuilder::MuAllBuilder(
    const Rcpp::List& x_hyper,
    const Rcpp::CharacterVector& group_levels
) {
  p_re = x_hyper.size();
  J = group_levels.size();
  X.reserve(p_re);
  row_idx.resize(static_cast<size_t>(p_re));
  for (int i = 0; i < p_re; ++i) {
    Rcpp::NumericMatrix X_k = Rcpp::as<Rcpp::NumericMatrix>(x_hyper[i]);
    X.push_back(X_k);

    Rcpp::CharacterVector rn;
    SEXP dn = X_k.attr("dimnames");
    if (!Rf_isNull(dn)) {
      Rcpp::List dnl(dn);
      if (dnl.size() >= 1 && !Rf_isNull(dnl[0])) {
        rn = Rcpp::CharacterVector(dnl[0]);
      }
    }

    std::vector<int>& idx = row_idx[static_cast<size_t>(i)];
    idx.resize(static_cast<size_t>(J));
    if (rn.size() == 0) {
      if (X_k.nrow() != J) {
        Rcpp::stop(
          "nrow(X_hyper[[%d]]) (%d) must equal length(group_levels) (%d).",
          i + 1, X_k.nrow(), J
        );
      }
      for (int j = 0; j < J; ++j) {
        idx[static_cast<size_t>(j)] = j;
      }
    } else {
      for (int j = 0; j < J; ++j) {
        const std::string lev = Rcpp::as<std::string>(group_levels[j]);
        int found = -1;
        for (int r = 0; r < rn.size(); ++r) {
          if (!Rcpp::CharacterVector::is_na(rn[r]) &&
              Rcpp::as<std::string>(rn[r]) == lev) {
            found = r;
            break;
          }
        }
        if (found < 0) {
          Rcpp::stop(
            "group level(s) not found in rownames(X_hyper[[%d]]): %s",
            i + 1, lev.c_str()
          );
        }
        idx[static_cast<size_t>(j)] = found;
      }
    }
  }
}

Rcpp::NumericMatrix MuAllBuilder::build(
    const std::vector<Rcpp::NumericVector>& fixef
) const {
  Rcpp::NumericMatrix mu_all(p_re, J);
  for (int i = 0; i < p_re; ++i) {
    const Rcpp::NumericMatrix& X_k = X[static_cast<size_t>(i)];
    const Rcpp::NumericVector& gamma_k = fixef[static_cast<size_t>(i)];
    if (gamma_k.size() != X_k.ncol()) {
      Rcpp::stop(
        "length(fixef[[%d]]) (%d) must equal ncol(X_hyper[[%d]]) (%d).",
        i + 1, gamma_k.size(), i + 1, X_k.ncol()
      );
    }
    const std::vector<int>& idx = row_idx[static_cast<size_t>(i)];
    for (int j = 0; j < J; ++j) {
      const int r = idx[static_cast<size_t>(j)];
      double s = 0.0;
      for (int c = 0; c < X_k.ncol(); ++c) {
        s += X_k(r, c) * gamma_k[c];
      }
      mu_all(i, j) = s;
    }
  }
  return mu_all;
}

Rcpp::NumericMatrix two_block_build_mu_all(
    const Rcpp::List& x_hyper,
    const Rcpp::List& fixef,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels
) {
  const std::vector<Rcpp::NumericVector> fixef_v =
    fixef_vec_from_list(fixef, re_names);
  MuAllBuilder builder(x_hyper, group_levels);
  Rcpp::NumericMatrix mu_all = builder.build(fixef_v);
  set_matrix_dimnames(mu_all, re_names, group_levels);
  return mu_all;
}

/// Refresh Block~1 prior precision for ING components from chain \code{tau2}
/// (port of \code{.two_block_block1_prior_with_tau2}; v5 \code{any_ing} semantics).
Rcpp::List two_block_block1_prior_with_tau2(
    const Rcpp::List& base_prior,
    const Rcpp::NumericVector& tau2_vec,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::NumericMatrix& mu_all
) {
  const int p_re = re_names.size();
  if (tau2_vec.size() != p_re) {
    Rcpp::stop("length(tau2_vec) must equal length(re_names).");
  }
  if (ptypes.size() != p_re) {
    Rcpp::stop("length(ptypes) must equal length(re_names).");
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("mu") = mu_all,
    Rcpp::Named("dispersion") = optional_list_elt(base_prior, "dispersion"),
    Rcpp::Named("ddef") = optional_list_elt(base_prior, "ddef")
  );

  bool any_ing = false;
  for (int k = 0; k < p_re; ++k) {
    if (Rcpp::as<std::string>(ptypes[k]) == "dIndependent_Normal_Gamma") {
      any_ing = true;
      break;
    }
  }

  if (!any_ing) {
    out["P"] = base_prior["P"];
    return out;
  }

  if (!base_prior.containsElementNamed("P")) {
    Rcpp::stop("base_prior must contain 'P' when refreshing ING precision.");
  }

  Rcpp::NumericMatrix P1 =
    Rcpp::clone(Rcpp::as<Rcpp::NumericMatrix>(base_prior["P"]));
  if (P1.nrow() != p_re || P1.ncol() != p_re) {
    Rcpp::stop("dim(base_prior$P) must be p_re x p_re.");
  }

  for (int k = 0; k < p_re; ++k) {
    if (Rcpp::as<std::string>(ptypes[k]) != "dIndependent_Normal_Gamma") {
      continue;
    }
    for (int c = 0; c < p_re; ++c) {
      P1(k, c) = 0.0;
      P1(c, k) = 0.0;
    }
    P1(k, k) = 1.0 / tau2_vec[k];
  }
  out["P"] = P1;
  return out;
}

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

/// Reorder Block~1 coefficient rows to \code{group_levels}
/// (port of R \code{match(group_levels, rownames(b_draw))} in draw_one_chain).
/// When \code{block_ids} is \code{NULL}, returns \code{b_draw} unchanged.
Rcpp::NumericMatrix two_block_reorder_b_to_group_levels(
    Rcpp::NumericMatrix b_draw,
    SEXP block_ids_sexp,
    const Rcpp::CharacterVector& group_levels
) {
  if (Rf_isNull(block_ids_sexp)) {
    return b_draw;
  }
  Rcpp::CharacterVector block_ids(block_ids_sexp);
  if (block_ids.size() == 0) {
    return b_draw;
  }

  const int J = group_levels.size();
  if (b_draw.nrow() != J || block_ids.size() != J) {
    Rcpp::stop("Block 1 group dimension mismatch during reorder.");
  }

  SEXP col_dn = R_NilValue;
  SEXP b_dn = b_draw.attr("dimnames");
  if (!Rf_isNull(b_dn)) {
    Rcpp::List b_dnl(b_dn);
    if (b_dnl.size() >= 2 && !Rf_isNull(b_dnl[1])) {
      col_dn = b_dnl[1];
    }
  }

  bool aligned = true;
  for (int g = 0; g < J; ++g) {
    if (Rcpp::CharacterVector::is_na(block_ids[g]) ||
        Rcpp::CharacterVector::is_na(group_levels[g]) ||
        Rcpp::as<std::string>(block_ids[g]) !=
          Rcpp::as<std::string>(group_levels[g])) {
      aligned = false;
      break;
    }
  }
  if (aligned) {
    Rcpp::NumericMatrix out = Rcpp::clone(b_draw);
    set_matrix_dimnames(out, group_levels, col_dn);
    return out;
  }

  Rcpp::NumericMatrix out(J, b_draw.ncol());
  for (int g = 0; g < J; ++g) {
    const std::string lev = Rcpp::as<std::string>(group_levels[g]);
    int src = -1;
    for (int r = 0; r < J; ++r) {
      if (!Rcpp::CharacterVector::is_na(block_ids[r]) &&
          Rcpp::as<std::string>(block_ids[r]) == lev) {
        src = r;
        break;
      }
    }
    if (src < 0) {
      Rcpp::stop("Block 1 group ids do not match group_levels.");
    }
    out(g, Rcpp::_) = b_draw(src, Rcpp::_);
  }
  set_matrix_dimnames(out, group_levels, col_dn);
  return out;
}

/// Block~1 prep + draw for one replicate chain.
/// Mirrors \code{.two_block_block1_prep_one_chain} then
/// \code{.two_block_block1_draw_one_chain} with all piecewise C++ flags
/// \code{TRUE} (\code{use_cpp_mu_all}, \code{use_cpp_prior_tau2},
/// \code{use_cpp_reorder}, \code{use_cpp_iters}).
Rcpp::List two_block_block1_one_chain_impl(
    int chain_i,
    const Rcpp::List& batch_fixef,
    const Rcpp::NumericVector& tau2_i,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& Z,
    SEXP groups,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    const Rcpp::List& x_hyper,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::List& block1_prior,
    bool is_gaussian,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    const std::string& family,
    const std::string& link,
    int Gridtype,
    int n_envopt
) {
  // ---- prep (.two_block_block1_prep_one_chain; piecewise C++ flags TRUE) ----
  Rcpp::List fixef_i = fixef_list_from_batch_chain(
    batch_fixef, chain_i, re_names
  );
  Rcpp::NumericMatrix mu_all = two_block_build_mu_all(
    x_hyper, fixef_i, re_names, group_levels
  );
  Rcpp::List prior_list = two_block_block1_prior_with_tau2(
    block1_prior, tau2_i, ptypes, re_names, mu_all
  );

  // ---- draw (.two_block_block1_draw_one_chain; piecewise C++ flags TRUE) ----
  Rcpp::List block_out;
  if (is_gaussian) {
    block_out = block_rNormalReg_cpp_export(
      1, y, Z, groups, prior_list, R_NilValue, offset, wt,
      f2_gauss, f3_gauss, Gridtype
    );
  } else {
    block_out = block_rNormalGLM_cpp_export(
      1, y, Z, groups, prior_list, R_NilValue, offset, wt,
      f2, f3, family, link, Gridtype, n_envopt,
      false, false, false
    );
  }

  Rcpp::List block_info = Rcpp::as<Rcpp::List>(block_out["block_info"]);
  Rcpp::NumericMatrix b_draw =
    Rcpp::as<Rcpp::NumericMatrix>(block_out["coefficients"]);
  const int J = group_levels.size();
  const int p_re = re_names.size();
  if (b_draw.nrow() != J || b_draw.ncol() != p_re) {
    Rcpp::stop(
      "Block 1 returned a %d x %d coefficient matrix; expected %d x %d.",
      b_draw.nrow(), b_draw.ncol(), J, p_re
    );
  }

  set_block_draw_coefficient_dimnames(b_draw, Z, block_info);
  b_draw = two_block_reorder_b_to_group_levels(
    b_draw, matrix_rownames_sexp(b_draw), group_levels
  );
  const double iters_mean = two_block_block1_iters_mean(block_out);

  return Rcpp::List::create(
    Rcpp::Named("mu_all") = mu_all,
    Rcpp::Named("prior_list") = prior_list,
    Rcpp::Named("b") = b_draw,
    Rcpp::Named("iters_mean") = iters_mean
  );
}

/// Block~1 prep + draw for all replicate chains (v6 batch driver).
Rcpp::List two_block_block1_all_chains_impl(
    Rcpp::NumericVector b_store,
    Rcpp::NumericVector iters_ranef,
    const Rcpp::List& batch_fixef,
    const Rcpp::NumericMatrix& batch_tau2,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& Z,
    SEXP groups,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    const Rcpp::List& x_hyper,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::List& block1_prior,
    bool is_gaussian,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    const std::string& family,
    const std::string& link,
    int Gridtype,
    int n_envopt,
    bool progbar,
    const std::string& progbar_prefix,
    bool progbar_finish_newline
) {
  Rcpp::IntegerVector b_dim = b_store.attr("dim");
  if (b_dim.size() != 3) {
    Rcpp::stop("'b_store' must be a 3-dimensional array (J x p_re x n).");
  }
  const int J = b_dim[0];
  const int p_re = b_dim[1];
  const int n = b_dim[2];
  if (batch_tau2.nrow() != n || batch_tau2.ncol() != p_re) {
    Rcpp::stop("dim(batch_tau2) must be c(n, p_re).");
  }
  if (iters_ranef.size() != n) {
    Rcpp::stop("length(iters_ranef) must equal n.");
  }
  if (re_names.size() != p_re) {
    Rcpp::stop("length(re_names) must equal p_re.");
  }

  Rcpp::NumericVector b_out = Rcpp::clone(b_store);
  Rcpp::NumericVector iters_out = Rcpp::clone(iters_ranef);
  const bool show_bar = progbar && n > 1;

  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt();
    if (show_bar) {
      glmbayes::progress::progress_bar(
        static_cast<double>(i + 1),
        static_cast<double>(n),
        progbar_prefix
      );
    }

    Rcpp::NumericVector tau2_i = batch_tau2(i, Rcpp::_);
    Rcpp::List one = two_block_block1_one_chain_impl(
      i + 1, batch_fixef, tau2_i, y, Z, groups, offset, wt,
      x_hyper, re_names, group_levels, ptypes, block1_prior,
      is_gaussian, f2, f3, f2_gauss, f3_gauss, family, link,
      Gridtype, n_envopt
    );
    Rcpp::NumericMatrix b_i = Rcpp::as<Rcpp::NumericMatrix>(one["b"]);
    store_b_chain(b_out, i, J, p_re, b_i);
    iters_out[i] += Rcpp::as<double>(one["iters_mean"]);
  }

  if (show_bar) {
    glmbayes::progress::progress_bar_finish(progbar_finish_newline);
  }

  return Rcpp::List::create(
    Rcpp::Named("b") = b_out,
    Rcpp::Named("iters_ranef") = iters_out
  );
}

} // namespace sim
} // namespace glmbayes
