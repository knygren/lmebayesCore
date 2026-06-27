// two_block_block1.cpp
// Block 1 batch helpers for run_sweep_outer_chains_v6 (R-oracle semantics).
// Migrated incrementally from the R driver in R/two_block_batch_gibbs.R.

#include "simfuncs.h"

namespace glmbayes {
namespace sim {

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
  dn[0] = rownames;
  dn[1] = Rf_isNull(colnames_src) ? R_NilValue : colnames_src;
  M.attr("dimnames") = dn;
}

} // namespace

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
