// lmebayesCore-only Rcpp exports. glmbayesCore's own exports are reached
// from R via glmbayesCore:::<name>_cpp_export (R/rcpp_wrappers*.R).

#include "RcppArmadillo.h"
#include "simfuncs.h"

using namespace glmbayes::sim;

// [[Rcpp::export]]
Rcpp::List two_block_rNormal_reg_cpp_export(
    int n,
    int m_convergence,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    const Rcpp::List& x_hyper,
    const Rcpp::List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const Rcpp::List& prior_list_block2,
    const Rcpp::List& fixef_start,
    const Rcpp::CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    int Gridtype = 2,
    int n_envopt = 1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false,
    bool progbar = true
) {
  return glmbayes::sim::two_block_rNormal_reg_cpp_export(
    n, m_convergence, y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    prior_list_block2, fixef_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt, Gridtype, n_envopt,
    use_parallel, use_opencl, verbose, progbar
  );
}

// [[Rcpp::export]]
Rcpp::List two_block_rNormal_reg_v5_cpp_export(
    int n,
    int m_convergence,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    SEXP block,
    const Rcpp::List& x_hyper,
    const Rcpp::List& prior_list_block1,
    SEXP dispersion_block1,
    SEXP ddef_block1,
    const Rcpp::List& pfamily_list,
    const Rcpp::List& fixef_start,
    const Rcpp::CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    int Gridtype = 2,
    int n_envopt = 1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false,
    Rcpp::Nullable<int> seed = R_NilValue,
    int seed_offset = 0,
    bool progbar = true,
    std::string stage_label = "",
    bool diag_sweeps = false,
    SEXP fixef_mode = R_NilValue,
    SEXP b_mode = R_NilValue
) {
  return glmbayes::sim::two_block_rNormal_reg_v5_cpp_export(
    n, m_convergence, y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    pfamily_list, fixef_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt, Gridtype, n_envopt,
    use_parallel, use_opencl, verbose, seed, seed_offset, progbar,
    stage_label, diag_sweeps, fixef_mode, b_mode
  );
}

// [[Rcpp::export]]
Rcpp::List two_block_rNormal_reg_staged_cpp_export(
    int n_main,
    int m_convergence_main,
    int n_pilot,
    int m_convergence_pilot,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::RObject& block,
    const Rcpp::List& x_hyper,
    const Rcpp::List& prior_list_block1,
    const Rcpp::RObject& dispersion_block1,
    const Rcpp::RObject& ddef_block1,
    const Rcpp::List& pfamily_list,
    const Rcpp::List& fixef_start,
    const Rcpp::CharacterVector& group_levels,
    const std::string& family,
    const std::string& link,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& wt,
    int Gridtype = 2,
    int n_envopt = 1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false,
    bool progbar_main = true,
    bool progbar_pilot = false
) {
  return glmbayes::sim::two_block_rNormal_reg_staged_cpp_export(
    n_main, m_convergence_main, n_pilot, m_convergence_pilot,
    y, x, block, x_hyper,
    prior_list_block1, dispersion_block1, ddef_block1,
    pfamily_list, fixef_start, group_levels,
    family, link, f2, f3, f2_gauss, f3_gauss,
    offset, wt, Gridtype, n_envopt,
    use_parallel, use_opencl, verbose, progbar_main, progbar_pilot
  );
}

// [[Rcpp::export]]
Rcpp::List two_block_block1_one_chain_v2_cpp_export(
    const Rcpp::List& fixef_i,
    const Rcpp::NumericVector& tau2_i,
    const Rcpp::List& design,
    const Rcpp::List& block1_prior,
    SEXP family,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss
) {
  return glmbayes::sim::two_block_block1_one_chain_v2_impl(
    fixef_i, tau2_i, design, block1_prior, family, ptypes,
    re_names, group_levels, f2, f3, f2_gauss, f3_gauss
  );
}

// [[Rcpp::export]]
Rcpp::List two_block_block1_all_chains_v2_internal_cpp_export(
    const Rcpp::List& fixef,
    int chain_i,
    const Rcpp::NumericMatrix& tau2,
    Rcpp::NumericVector b_store,
    Rcpp::NumericVector iters_ranef,
    const Rcpp::List& design,
    const Rcpp::List& block1_prior,
    SEXP family,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    bool use_cpp_tau2_row,
    bool use_cpp_b_slice,
    bool use_cpp_iters_ranef_add
) {
  return glmbayes::sim::two_block_block1_all_chains_v2_internal_impl(
    fixef, chain_i, tau2, b_store, iters_ranef, design, block1_prior, family,
    ptypes, re_names, group_levels, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add
  );
}

// [[Rcpp::export]]
Rcpp::List two_block_block1_all_chains_v2_internal_loop_cpp_export(
    int n,
    const Rcpp::List& fixef,
    const Rcpp::NumericMatrix& tau2,
    const Rcpp::NumericMatrix& b_in_master,
    const Rcpp::NumericVector& iters_ranef_in,
    const Rcpp::List& design,
    const Rcpp::List& block1_prior,
    SEXP family,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    bool use_cpp_tau2_row,
    bool use_cpp_b_slice,
    bool use_cpp_iters_ranef_add,
    bool show_bar,
    const std::string& progbar_prefix,
    bool progbar_finish_newline
) {
  return glmbayes::sim::two_block_block1_all_chains_v2_internal_loop_impl(
    n, fixef, tau2, b_in_master, iters_ranef_in, design, block1_prior, family, ptypes,
    re_names, group_levels, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add,
    show_bar, progbar_prefix, progbar_finish_newline
  );
}

// [[Rcpp::export]]
void two_block_block1_all_chains_cpp_export(
    int n,
    const Rcpp::List& fixef,
    const Rcpp::NumericMatrix& tau2,
    Rcpp::NumericVector& b,
    Rcpp::NumericVector& iters_ranef,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels,
    const Rcpp::List& design,
    const Rcpp::List& block1_prior,
    SEXP family,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::Function& f2,
    const Rcpp::Function& f3,
    const Rcpp::Function& f2_gauss,
    const Rcpp::Function& f3_gauss,
    bool use_cpp_tau2_row,
    bool use_cpp_b_slice,
    bool use_cpp_iters_ranef_add,
    bool show_bar,
    const std::string& progbar_prefix,
    bool progbar_finish_newline
) {
  glmbayes::sim::two_block_block1_all_chains_impl(
    n, fixef, tau2, b, iters_ranef, re_names, group_levels, design,
    block1_prior, family, ptypes, f2, f3, f2_gauss, f3_gauss,
    use_cpp_tau2_row, use_cpp_b_slice, use_cpp_iters_ranef_add,
    show_bar, progbar_prefix, progbar_finish_newline
  );
}

