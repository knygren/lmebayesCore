// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file rng_utils.h
 * @brief Thread‑safe uniform sampling and truncated‑distribution utilities
 *        used throughout glmbayes.
 *
 * @namespace glmbayes::rng
 * @brief Random‑number generators for uniform draws and truncated Normal
 *        and truncated inverse‑Gamma distributions, including numerically
 *        stable “central truncation” variants.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - rng_utils.cpp
 *     - (optionally) supporting truncated‑distribution source files
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - RcppParallel‑based samplers (Normal, Normal–Gamma, and GLM samplers)
 *     - Envelope construction routines requiring truncated inverse‑Gamma draws
 *     - Diagnostic and simulation utilities across the glmbayes backend
 *
 * @section Responsibilities
 *   Provides:
 *     - a thread‑safe uniform RNG (`runif_safe`) suitable for parallel contexts,
 *     - truncated inverse‑Gamma samplers (`rinvgamma_ct_safe`, `rinvgamma_ct`)
 *       for dispersion parameters restricted to \([d_{\min}, d_{\max}]\),
 *     - truncated Normal sampler (`rnorm_ct`) for draws restricted to
 *       \([l_{\mathrm{lower}}, l_{\mathrm{upper}}]\),
 *     - log‑density evaluation for truncated inverse‑Gamma distributions
 *       (`log_p_inv_gamma_ct_safe`).
 *
 *   The suffix `_ct` denotes “central truncation,” meaning all functions
 *   operate on a parameter restricted to a finite interval and are designed
 *   to remain numerically stable even when the bounds are tight or the shape
 *   and rate parameters induce heavy tails.
 */

#ifndef RNG_UTILS_H
#define RNG_UTILS_H

#include "RcppArmadillo.h"

namespace glmbayes{

namespace rng {

/// Put one around every R-RNG kernel call made inside a longer C++ loop.
/// The R callbacks f2/f3 reach Rcpp exports whose own RNGScope reloads
/// .Random.seed on entry, and Rcpp scopes do not nest. Without a scope that
/// writes R's RNG state back after each call, the next call restarts from
/// the stream the previous one already used.
class RSyncScope {
 public:
  RSyncScope() {}

 private:
  RSyncScope(const RSyncScope&);
  RSyncScope& operator=(const RSyncScope&);
  Rcpp::RNGScope scope_;
};

// Thread-safe uniform RNG [0, 1)
double runif_safe();

double rinvgamma_ct_safe(double shape,
                       double rate,
                       double disp_upper,
                       double disp_lower);



double  rnorm_ct(double lgrt,double lglt,double mu,double sigma);

double rinvgamma_ct(double shape,double rate,double disp_upper,double disp_lower);

double log_p_inv_gamma_ct_safe(double disp_lower,
                               double disp_upper,
                               double shape,
                               double rate);

}
}
#endif
