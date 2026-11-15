// DEV ONLY -- does not migrate to glmbayesCore.
// Exports for the two-block replay / record / prior-tracking diagnostics.
// R wrappers: R/dev_rcpp_wrappers.R.

#include "RcppArmadillo.h"
#include "dev_twoblock.h"

using glmbayes::sim::twoblock::TwoBlockIndexSpec;

// [[Rcpp::export]]
Rcpp::List dev_two_block_log_new_cpp_export(const Rcpp::List& spec, int n, int M) {
  return lmebayes::dev::twoblock::log_new(TwoBlockIndexSpec(spec), n, M);
}

// [[Rcpp::export]]
Rcpp::List dev_two_block_sweep_cpp_export(const Rcpp::List& spec,
                                          const Rcpp::NumericVector& b,
                                          const Rcpp::List& fixef,
                                          const Rcpp::NumericMatrix& iters,
                                          const Rcpp::NumericVector& iters_ranef,
                                          int inner_sweeps, bool use_mode,
                                          bool progbar,
                                          const std::string& stage_label,
                                          Rcpp::Nullable<Rcpp::List> tape,
                                          bool record) {
  const TwoBlockIndexSpec s(spec);
  const Rcpp::IntegerVector dim = b.attr("dim");
  if (dim.size() != 3) {
    Rcpp::stop("b must be a J x p_re x n array.");
  }
  lmebayes::dev::twoblock::ReplayObserver obs(s, dim[2], inner_sweeps, tape,
                                              record);
  const Rcpp::List st = glmbayes::sim::twoblock::two_block_sweep(
      s, b, fixef, iters, iters_ranef, inner_sweeps, use_mode, progbar,
      stage_label, 0.0, &obs);
  return Rcpp::List::create(
      Rcpp::Named("b") = st["b"],
      Rcpp::Named("fixef") = st["fixef"],
      Rcpp::Named("iters") = st["iters"],
      Rcpp::Named("iters_ranef") = st["iters_ranef"],
      Rcpp::Named("sweep_mean") = st["sweep_mean"],
      Rcpp::Named("sweep_cov") = st["sweep_cov"],
      Rcpp::Named("log") = obs.log(),
      Rcpp::Named("prior_sum") = obs.prior_sum,
      Rcpp::Named("prior_sumsq") = obs.prior_sumsq);
}
