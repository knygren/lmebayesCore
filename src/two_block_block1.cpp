// two_block_block1.cpp
// Block 1 batch helpers for run_sweep_outer_chains_v6 (R-oracle semantics).
// Migrated incrementally from the R driver in R/two_block_batch_gibbs.R.

#include "simfuncs.h"

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

} // namespace sim
} // namespace glmbayes
