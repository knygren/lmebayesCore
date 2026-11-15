// Typed C++ entry points into glmbayesCore, implemented in
// glmbayesCore_calls.cpp on top of <glmbayesCore/glmbayesCore_capi.h>
// (same layout as glmbayesCore's kernel_loader.cpp over opencltools).

#ifndef LMEBAYESCORE_GLMBAYESCORE_CALLS_H
#define LMEBAYESCORE_GLMBAYESCORE_CALLS_H

#include "RcppArmadillo.h"

#include <string>

namespace lmecore {

void progress_bar(double x, double N, const std::string& prefix = "");

void progress_bar_finish(bool newline = true);

Rcpp::List group_rNormalReg(
    int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x,
    SEXP group, SEXP prior_list, SEXP prior_lists,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2, const Rcpp::Function& f3, int Gridtype);

Rcpp::List group_rNormalGLM(
    int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x,
    SEXP group, SEXP prior_list, SEXP prior_lists,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2, const Rcpp::Function& f3,
    const std::string& family, const std::string& link, int Gridtype,
    int n_envopt, bool use_parallel, bool use_opencl, bool verbose);

Rcpp::List rNormalReg(
    int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu, const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& wt,
    double dispersion, const Rcpp::Function& f2, const Rcpp::Function& f3,
    const Rcpp::NumericVector& start, const std::string& family,
    const std::string& link, int Gridtype);

Rcpp::List rIndepNormalGammaReg(
    int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu, const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& wt,
    double shape, double rate, double max_disp_perc,
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    int Gridtype, int n_envopt, bool use_parallel, bool use_opencl,
    bool verbose, bool progbar);

}  // namespace lmecore

#endif
