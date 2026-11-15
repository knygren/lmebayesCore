#include "glmbayesCore_calls.h"

#include <glmbayesCore/glmbayesCore_capi.h>

namespace lmecore {

using Rcpp::RObject;
using Rcpp::wrap;

void progress_bar(double x, double N, const std::string& prefix) {
  const RObject x_ = wrap(x);
  const RObject N_ = wrap(N);
  const RObject prefix_ = wrap(prefix);
  glmbayesCore_progress_bar(x_, N_, prefix_);
}

void progress_bar_finish(bool newline) {
  const RObject newline_ = wrap(newline);
  glmbayesCore_progress_bar_finish(newline_);
}

Rcpp::List group_rNormalReg(
    int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x,
    SEXP group, SEXP prior_list, SEXP prior_lists,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2, const Rcpp::Function& f3, int Gridtype) {
  const RObject n_ = wrap(n);
  const RObject Gridtype_ = wrap(Gridtype);
  return glmbayesCore_group_rNormalReg(
      n_, y, x, group, prior_list, prior_lists, offset, wt, f2, f3,
      Gridtype_);
}

Rcpp::List group_rNormalGLM(
    int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x,
    SEXP group, SEXP prior_list, SEXP prior_lists,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& wt,
    const Rcpp::Function& f2, const Rcpp::Function& f3,
    const std::string& family, const std::string& link, int Gridtype,
    int n_envopt, bool use_parallel, bool use_opencl, bool verbose) {
  const RObject n_ = wrap(n);
  const RObject family_ = wrap(family);
  const RObject link_ = wrap(link);
  const RObject Gridtype_ = wrap(Gridtype);
  const RObject n_envopt_ = wrap(n_envopt);
  const RObject use_parallel_ = wrap(use_parallel);
  const RObject use_opencl_ = wrap(use_opencl);
  const RObject verbose_ = wrap(verbose);
  return glmbayesCore_group_rNormalGLM(
      n_, y, x, group, prior_list, prior_lists, offset, wt, f2, f3, family_,
      link_, Gridtype_, n_envopt_, use_parallel_, use_opencl_, verbose_);
}

Rcpp::List rNormalReg(
    int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu, const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& wt,
    double dispersion, const Rcpp::Function& f2, const Rcpp::Function& f3,
    const Rcpp::NumericVector& start, const std::string& family,
    const std::string& link, int Gridtype) {
  const RObject n_ = wrap(n);
  const RObject dispersion_ = wrap(dispersion);
  const RObject family_ = wrap(family);
  const RObject link_ = wrap(link);
  const RObject Gridtype_ = wrap(Gridtype);
  return glmbayesCore_rNormalReg(
      n_, y, x, mu, P, offset, wt, dispersion_, f2, f3, start, family_,
      link_, Gridtype_);
}

Rcpp::List rIndepNormalGammaReg(
    int n, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& mu, const Rcpp::NumericMatrix& P,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& wt,
    double shape, double rate, double max_disp_perc,
    Rcpp::Nullable<Rcpp::NumericVector> disp_lower,
    Rcpp::Nullable<Rcpp::NumericVector> disp_upper,
    int Gridtype, int n_envopt, bool use_parallel, bool use_opencl,
    bool verbose, bool progbar) {
  const RObject n_ = wrap(n);
  const RObject shape_ = wrap(shape);
  const RObject rate_ = wrap(rate);
  const RObject max_disp_perc_ = wrap(max_disp_perc);
  const RObject Gridtype_ = wrap(Gridtype);
  const RObject n_envopt_ = wrap(n_envopt);
  const RObject use_parallel_ = wrap(use_parallel);
  const RObject use_opencl_ = wrap(use_opencl);
  const RObject verbose_ = wrap(verbose);
  const RObject progbar_ = wrap(progbar);
  return glmbayesCore_rIndepNormalGammaReg(
      n_, y, x, mu, P, offset, wt, shape_, rate_, max_disp_perc_,
      disp_lower, disp_upper, Gridtype_, n_envopt_, use_parallel_,
      use_opencl_, verbose_, progbar_);
}

}  // namespace lmecore
