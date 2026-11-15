// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-

/**
 * @file simfuncs_groups.h
 * @brief Staging shim for glmbayesCore's glmbayes::sim::group kernels.
 *
 * @namespace glmbayes::sim::group
 *
 *   glmbayesCore names the grouped Normal-GLM kernel
 *   glmbayes::sim::group::rNormalGLMGroups; lmebayesCore still has it as
 *   glmbayes::sim::rNormalGLMBlocks with the same arguments and return
 *   value. Staged code calls the glmbayesCore name through this forwarder.
 *
 * @section Migration
 *   Drop this file; glmbayesCore's own simfuncs_groups.h declares the real
 *   function.
 */

#ifndef LMEBAYES_SIMFUNCS_GROUPS_SHIM_H
#define LMEBAYES_SIMFUNCS_GROUPS_SHIM_H

#include "RcppArmadillo.h"
#include "simfuncs.h"

#include <string>

namespace glmbayes {
namespace sim {
namespace group {

inline Rcpp::List rNormalGLMGroups(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
    Rcpp::NumericVector dispersion,
    Rcpp::NumericMatrix mu,
    Rcpp::List P_blocks,
    bool prior_by_block,
    Rcpp::List row_blocks,
    Rcpp::Function f2,
    Rcpp::Function f3,
    std::string family,
    std::string link,
    int Gridtype = 2,
    int n_envopt = -1,
    bool use_parallel = true,
    bool use_opencl = false,
    bool verbose = false) {
  return glmbayes::sim::rNormalGLMBlocks(
      n, y, x, offset, wt, dispersion, mu, P_blocks, prior_by_block,
      row_blocks, f2, f3, family, link, Gridtype, n_envopt, use_parallel,
      use_opencl, verbose);
}

}  // namespace group
}  // namespace sim
}  // namespace glmbayes

#endif
