// TwoBlockIndexSpec: the R `spec` list from .two_block_spec(), converted
// once per call so the sweep steps do not look names up in a List.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "simfuncs_twoblock.h"

using namespace Rcpp;

namespace glmbayes {
namespace sim {
namespace twoblock {

namespace {

template <typename T>
std::vector<T> list_elements(const List& x) {
  std::vector<T> out;
  out.reserve(x.size());
  for (int k = 0; k < x.size(); ++k) {
    out.push_back(as<T>(x[k]));
  }
  return out;
}

}  // namespace

TwoBlockIndexSpec::TwoBlockIndexSpec(const List& spec)
    : y(as<NumericVector>(spec["y"])),
      D(as<NumericMatrix>(spec["D"])),
      offset1(as<NumericVector>(spec["offset1"])),
      wt1(as<NumericVector>(spec["wt1"])),
      disp1(as<NumericVector>(spec["disp1"])),
      P_D(as<NumericMatrix>(spec["P_D"])),
      obs_rows(as<List>(spec["obs_rows"])),
      d_col(as<IntegerVector>(spec["d_col"])),
      W(list_elements<NumericMatrix>(as<List>(spec["W"]))),
      w_for_group(list_elements<IntegerVector>(as<List>(spec["w_for_group"]))),
      b_for_wrow(list_elements<IntegerVector>(as<List>(spec["b_for_wrow"]))),
      w_col(list_elements<IntegerVector>(as<List>(spec["w_col"]))),
      mu2(list_elements<NumericVector>(as<List>(spec["mu2"]))),
      P2(list_elements<NumericMatrix>(as<List>(spec["P2"]))),
      disp2(as<NumericVector>(spec["disp2"])),
      f2(as<Function>(spec["f2"])),
      f3(as<Function>(spec["f3"])),
      f2_gauss(as<Function>(spec["f2_gauss"])),
      f3_gauss(as<Function>(spec["f3_gauss"])),
      family(as<std::string>(spec["family"])),
      link(as<std::string>(spec["link"])),
      gridtype(as<int>(spec["gridtype"])),
      n_envopt(as<int>(spec["n_envopt"])),
      ridge(false),
      l1(D.ncol()),
      J(obs_rows.size()),
      p_re(d_col.size()) {
  const std::size_t p = static_cast<std::size_t>(p_re);
  if (W.size() != p || w_for_group.size() != p || b_for_wrow.size() != p ||
      w_col.size() != p || mu2.size() != p || P2.size() != p ||
      disp2.size() != p_re) {
    stop("spec: every population entry must have one element per component.");
  }
  if (spec.containsElementNamed("ridge") && as<bool>(spec["ridge"])) {
    ridge = true;
    ridge_A = as<NumericVector>(spec["ridge_A"]);
    ridge_c = as<NumericMatrix>(spec["ridge_c"]);
    ridge_P = as<List>(spec["ridge_P"]);
    if (ridge_A.size() != static_cast<R_xlen_t>(l1) * l1 * J ||
        ridge_c.nrow() != l1 || ridge_c.ncol() != J || ridge_P.size() != J) {
      stop("spec: ridge_A, ridge_c and ridge_P must have one entry per group.");
    }
  }
}

}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes
