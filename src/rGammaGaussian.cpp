#include "RcppArmadillo.h"
using namespace Rcpp;


namespace glmbayes {
namespace sim {

Rcpp::List rGammaGaussian(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector beta,
    Rcpp::NumericVector wt,
    Rcpp::NumericVector alpha,
    double shape,
    double rate,
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    bool verbose
) {
  // --------------------------------------------------------------
  // Placeholder implementation
  // This shell compiles and returns a minimal list with the
  // expected structure. You will fill in the conjugate sampler later.
  // --------------------------------------------------------------
  
  Rcpp::List out;
  out["dispersion"] = Rcpp::NumericVector(n, NA_REAL);
  out["iters"]      = Rcpp::IntegerVector(n, 0);
  
  return out;
}

} // namespace sim
} // namespace glmbayes