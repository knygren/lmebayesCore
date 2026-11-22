#include "RcppArmadillo.h"
#include "simfuncs_twoblock.h"

using namespace glmbayes::sim::twoblock;

// =============================================================================
// Two-block indexed: groups / population steps and the sweep loop
// Callers: .rGLMM_sweep_cpp (via .two_block_sweep_cpp); the per-step exports
//          back the R/C++ parity checks (via .two_block_groups_*_cpp,
//          .two_block_population_*_cpp).
// Index arguments `i` (chain) and `k` (component) are 1-based, as in R.
// Migration: move into glmbayesCore export_wrappers.cpp, two-block tier.
// =============================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix two_block_groups_prior_mean_cpp_export(const Rcpp::List& spec,
                                                           const Rcpp::List& fixef,
                                                           int i) {
  return groups_prior_mean(TwoBlockIndexSpec(spec), fixef, i - 1);
}

// [[Rcpp::export]]
Rcpp::List two_block_groups_draw_cpp_export(const Rcpp::List& spec,
                                            const Rcpp::NumericMatrix& mu_D) {
  const GroupsDraw d = groups_draw(TwoBlockIndexSpec(spec), mu_D);
  return Rcpp::List::create(
      Rcpp::Named("coef") = d.coef,
      Rcpp::Named("coef_mode") = d.coef_mode,
      Rcpp::Named("iters") = d.iters);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix two_block_groups_place_cpp_export(const Rcpp::List& spec,
                                                      const Rcpp::NumericMatrix& coef) {
  return groups_place(TwoBlockIndexSpec(spec), coef);
}

// [[Rcpp::export]]
Rcpp::NumericVector two_block_population_response_cpp_export(
    const Rcpp::List& spec, const Rcpp::NumericMatrix& b_i, int k) {
  return population_response(TwoBlockIndexSpec(spec), b_i, k - 1);
}

// [[Rcpp::export]]
Rcpp::List two_block_population_draw_cpp_export(const Rcpp::List& spec,
                                                const Rcpp::NumericVector& y_k,
                                                int k) {
  const PopulationDraw d = population_draw(TwoBlockIndexSpec(spec), y_k, k - 1);
  return Rcpp::List::create(
      Rcpp::Named("coef") = d.coef,
      Rcpp::Named("coef_mode") = d.coef_mode);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix two_block_population_write_cpp_export(
    const Rcpp::List& spec, const Rcpp::NumericMatrix& fixef_k, int i, int k,
    const Rcpp::NumericVector& coef) {
  Rcpp::NumericMatrix out = Rcpp::clone(fixef_k);
  population_write(TwoBlockIndexSpec(spec), out, i - 1, k - 1, coef);
  return out;
}

// [[Rcpp::export]]
Rcpp::List two_block_sweep_cpp_export(const Rcpp::List& spec,
                                      const Rcpp::NumericVector& b,
                                      const Rcpp::List& fixef,
                                      const Rcpp::NumericMatrix& iters,
                                      const Rcpp::NumericVector& iters_ranef,
                                      int inner_sweeps, bool use_mode,
                                      bool progbar,
                                      const std::string& stage_label,
                                      double rwt, bool group_major,
                                      bool group_check) {
  return two_block_sweep(TwoBlockIndexSpec(spec), b, fixef, iters, iters_ranef,
                         inner_sweeps, use_mode, progbar, stage_label, rwt,
                         nullptr, group_major, group_check);
}
