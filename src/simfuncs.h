// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file simfuncs.h
 * @brief Simulation and posterior sampling routines for glmbayes.
 *
 * @namespace glmbayes::sim
 * @brief Core Normal, Normal–Gamma, and independent Normal–Gamma samplers.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - rNormalGLM.cpp
 *     - rNormalGLMBlocks.cpp
 *     - rNormalRegBlocks.cpp
 *     - block_utils.cpp
 *     - rIndepNormalGammaReg.cpp
 *     - rNormalGammaReg.cpp
 *     - rGammaGaussian.cpp
 *     - rGammaGamma.cpp
 *     - glmb_Standardize_Model.cpp
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - export_wrappers.cpp 
 *
 * @section Responsibilities
 *   Provides simulation kernels for:
 *     - Normal GLM posterior draws (standardized and unstandardized)
 *     - Normal–Gamma regression models
 *     - Independent Normal–Gamma regression models (standard and parallel variants)
 *     - Gamma–Gaussian and Gamma–Gamma models
 *
 *   All routines:
 *     - assume validated inputs from R wrappers,
 *     - use Rcpp::List for structured return objects,
 *   
 *   Some routines:
 *     - rely on envelope objects and f2/f3 functions for accept–reject sampling,
 *     - support optional parallelization and OpenCL where applicable.
 */

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef GLMBAYES_SIM_H
#define GLMBAYES_SIM_H


// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include <vector>

using namespace Rcpp;


namespace glmbayes {

namespace sim {

Rcpp::List  rNormalGLM_std(int n,
                               NumericVector y,
                               NumericMatrix x,
                               NumericMatrix mu,
                               NumericMatrix P,
                               NumericVector alpha,
                               NumericVector wt,
                               Function f2,
                               Rcpp::List  Envelope,
                               Rcpp::CharacterVector   family,
                               Rcpp::CharacterVector   link, 
                               int progbar=1,
                               bool verbose = false                                 
                                 );

Rcpp::List rNormalGLM(
    int n,
    NumericVector y,
    NumericMatrix x,
    NumericVector mu,
    NumericMatrix P,
    NumericVector offset,
    NumericVector wt,
    double dispersion,
    Function f2,
    Function f3,
    NumericVector start,
    std::string family = "binomial",
    std::string link = "logit",
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
);

/// C++ counterpart to R \c rNormalGLM_reg_block(): loop over \c row_blocks,
/// call \c rNormalGLM() per block. Not exported to R until Phase 2.
Rcpp::List rNormalGLMBlocks(
    int n,
    NumericVector y,
    NumericMatrix x,
    NumericVector offset,
    NumericVector wt,
    NumericVector dispersion,
    NumericMatrix mu,
    List P_blocks,
    bool prior_by_block,
    List row_blocks,
    Function f2,
    Function f3,
    std::string family,
    std::string link,
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false
);

Rcpp::List rNormalReg(int n,NumericVector y,NumericMatrix x,
                         NumericVector mu,NumericMatrix P,
                         NumericVector offset,NumericVector wt,
                         double dispersion,
                         Function f2,Function f3,
                         NumericVector start,
                         std::string family="gaussian",
                         std::string link="identity",
                         int Gridtype=2      
);

/// C++ block loop for Gaussian regression: calls rNormalReg() per block.
Rcpp::List rNormalRegBlocks(
    int n,
    NumericVector y,
    NumericMatrix x,
    NumericVector offset,
    NumericVector wt,
    NumericVector dispersion,
    NumericMatrix mu,
    List P_blocks,
    bool prior_by_block,
    List row_blocks,
    Function f2,
    Function f3,
    int Gridtype = 2
);

Rcpp::List normalize_block_cpp(SEXP block, int l2);

Rcpp::List block_rNormalReg_cpp_export(
    int n,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    const NumericVector& offset,
    const NumericVector& wt,
    const Function& f2,
    const Function& f3,
    int Gridtype
);

Rcpp::List block_rNormalGLM_cpp_export(
    int n,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    const NumericVector& offset,
    const NumericVector& wt,
    const Function& f2,
    const Function& f3,
    const std::string& family,
    const std::string& link,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose
);

Rcpp::List block_rNormalGLM_cpp_export(
    int n,
    const NumericVector& y,
    const NumericMatrix& x,
    SEXP block,
    SEXP prior_list,
    SEXP prior_lists,
    const NumericVector& offset,
    const NumericVector& wt,
    const Function& f2,
    const Function& f3,
    const std::string& family,
    const std::string& link,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose
);

/// C++ port of the two_block_rNormal_reg() Gibbs loop (twoBlockGibbs.cpp).
/// Per inner step: mu_all -> Block 1 via block_rNormal{Reg,GLM}_cpp_export ->
/// Block 2 via rNormalReg() per RE component (port-only; algorithm unchanged).
Rcpp::List two_block_rNormal_reg_cpp_export(
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
);

/// v2 of the two-block Gibbs driver (twoBlockGibbs.cpp): Block 2 priors are
/// pfamily objects (dNormal / dIndependent_Normal_Gamma); dispatch on the
/// pfamily type string. Returns dispersion_fixef_draws in addition to the
/// v1 fields. Development track; v1 above is the frozen regression baseline.
Rcpp::List two_block_rNormal_reg_v2_cpp_export(
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
);

/// v5 short-chain driver (twoBlockGibbs.cpp): sweep-outer loop order
/// (m then i) with per-sweep chain progress bars.
Rcpp::List two_block_rNormal_reg_v5_cpp_export(
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
);

/// Staged v2 driver (twoBlockGibbsStaged.cpp): pilot replicate chains,
/// Hotelling chi-squared vs fixef_start, then main replicate chains.
/// Pilot eigenvalue upper bounds are computed in the R wrapper (Phase 2b).
Rcpp::List two_block_rNormal_reg_staged_cpp_export(
    int n_main,
    int m_convergence_main,
    int n_pilot,
    int m_convergence_pilot,
    const NumericVector& y,
    const NumericMatrix& x,
    const Rcpp::RObject& block,
    const List& x_hyper,
    const List& prior_list_block1,
    const Rcpp::RObject& dispersion_block1,
    const Rcpp::RObject& ddef_block1,
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
    bool progbar_main,
    bool progbar_pilot
);

/// Cached \code{X_hyper} row lookups for \code{build_mu_all} (two_block_block1.cpp).
struct MuAllBuilder {
  int p_re;
  int J;
  std::vector<Rcpp::NumericMatrix> X;
  std::vector<std::vector<int>> row_idx;

  MuAllBuilder(
      const Rcpp::List& x_hyper,
      const Rcpp::CharacterVector& group_levels
  );
  Rcpp::NumericMatrix build(
      const std::vector<Rcpp::NumericVector>& fixef
  ) const;
};

/// Block~1: \code{mu_all} matrix (port of \code{build_mu_all()}).
Rcpp::NumericMatrix two_block_build_mu_all(
    const Rcpp::List& x_hyper,
    const Rcpp::List& fixef,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::CharacterVector& group_levels
);

/// Block~1: ING \code{P} refresh from \code{tau2} (two_block_block1.cpp).
Rcpp::List two_block_block1_prior_with_tau2(
    const Rcpp::List& base_prior,
    const Rcpp::NumericVector& tau2_vec,
    const Rcpp::CharacterVector& ptypes,
    const Rcpp::CharacterVector& re_names,
    const Rcpp::NumericMatrix& mu_all
);

/// Block~1: mean envelope iters across groups (two_block_block1.cpp).
double two_block_block1_iters_mean(const Rcpp::List& block_out);

/// All-chains step A: \code{batch$tau2[chain_i, ]} (two_block_block1.cpp).
Rcpp::NumericVector batch_tau2_chain_row(
    const Rcpp::NumericMatrix& batch_tau2,
    int chain_i
);

/// All-chains step C: \code{batch$b[, , chain_i] <- b_draw} (two_block_block1.cpp).
void batch_b_assign_slice(
    Rcpp::NumericVector& b_store,
    int chain_i,
    const Rcpp::NumericMatrix& b_draw
);

/// All-chains step D: \code{batch$iters_ranef[chain_i] += iters_mean}
/// (two_block_block1.cpp).
void batch_iters_ranef_add(
    Rcpp::NumericVector& iters_ranef,
    int chain_i,
    double iters_mean
);

/// Block~1: reorder coefficient rows to group_levels (two_block_block1.cpp).
Rcpp::NumericMatrix two_block_reorder_b_to_group_levels(
    Rcpp::NumericMatrix b_draw,
    SEXP block_ids,
    const Rcpp::CharacterVector& group_levels
);

/// Block~1 prep + draw for one chain (two_block_block1.cpp).
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
);

/// Block~1 prep + draw for all chains (two_block_block1.cpp).
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
);

/// Block~2 align (twoBlockGibbs.cpp): port of R two_block_align_b_to_xhyper().
Rcpp::NumericVector two_block_align_b_to_xhyper_cpp(
    Rcpp::NumericVector b_vec,
    Rcpp::NumericMatrix X_k,
    Rcpp::CharacterVector group_levels
);

/// Block~2 one-chain update (twoBlockGibbs.cpp): native align + rglmb per RE.
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
);

Rcpp::List rNormalGLM_optim_poisson_log(
    const NumericVector& parin,
    const NumericVector& y,
    const NumericMatrix& x,
    const NumericVector& mu1,
    const NumericMatrix& P,
    const NumericVector& alpha,
    const NumericVector& wt2
);


Rcpp::List  rIndepNormalGammaReg_std(int n,NumericVector y,NumericMatrix x,
                                          NumericMatrix mu, /// This is typically standardized to be a zero vector
                                          NumericMatrix P, /// Part of prior precision shifted to the likelihood
                                          NumericVector alpha,NumericVector wt,
                                          Function f2,Rcpp::List  Envelope,
                                          Rcpp::List  gamma_list,
                                          Rcpp::List  UB_list,
                                          Rcpp::CharacterVector   family,Rcpp::CharacterVector   link,
                                          bool progbar=true,
                                          bool verbose=false);

Rcpp::List rIndepNormalGammaReg_std_parallel(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,   // typically standardized to be a zero vector
    Rcpp::NumericMatrix P,    // part of prior precision shifted to the likelihood
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    Rcpp::Function f2,        // kept for signature parity
    Rcpp::List Envelope,
    Rcpp::List gamma_list,
    Rcpp::List UB_list,
    Rcpp::CharacterVector family,
    Rcpp::CharacterVector link,
    bool progbar = true,
    bool verbose = false
);


Rcpp::List rNormalGammaReg(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
    double shape,
    double rate,
    Rcpp::Nullable<double> max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose
);


Rcpp::List rIndepNormalGammaReg(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    int Gridtype,
    int n_envopt,
    bool use_parallel,
    bool use_opencl,
    bool verbose,
    bool progbar
);

Rcpp::List rGammaGaussian(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector beta,
    Rcpp::NumericVector wt,
    Rcpp::NumericVector alpha,
    double shape,
    double rate,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose = false
);


Rcpp::List rGammaGamma(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector beta,
    Rcpp::NumericVector wt,
    Rcpp::NumericVector alpha,
    double shape,
    double rate,
    double max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose = false
);

Rcpp::List glmb_Standardize_Model(
    NumericVector y, 
    NumericMatrix x,   // Original design matrix (to be adjusted)
    NumericMatrix P,   // Prior Precision Matrix (to be adjusted)
    NumericMatrix bstar, // Posterior Mode from optimization (to be adjusted)
    NumericMatrix A1  // Precision for Log-Posterior at posterior mode (to be adjusted)
);


}
}

#endif