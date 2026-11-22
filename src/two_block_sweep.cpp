// Indexed two-block sweep loop. R twin: .two_block_sweep_r() in
// R/two_block_sweep_indexed.R; same loop order and progress output when
// group_major is false. The default group-major groups block loops over
// groups outside and chains inside, with progress counted in groups.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "glmbayesCore_calls.h"
#include "simfuncs_twoblock.h"

using namespace Rcpp;

namespace glmbayes {
namespace sim {
namespace twoblock {

std::string sweep_progbar_prefix(const std::string& stage_label, int m,
                                 int inner_sweeps, const std::string& phase) {
  std::string out = stage_label.empty() ? std::string("")
                                        : std::string("[") + stage_label + "] ";
  out += "sweep " + std::to_string(m) + "/" + std::to_string(inner_sweeps);
  if (!phase.empty()) {
    out += " " + phase;
  }
  return out + ": ";
}

namespace {

// Cross-chain mean and covariance (denominator n - 1, NA when n < 2) of the
// stacked population coefficients after sweep `mi`; components in order,
// columns in stored order. R twin: colMeans / stats::cov in .two_block_sweep_r.
void sweep_summary(const List& fixef, const std::vector<int>& col_off, int n,
                   int p_tot, int mi, std::vector<double>& stacked,
                   NumericMatrix& sweep_mean, NumericVector& sweep_cov) {
  for (int k = 0; k < static_cast<int>(col_off.size()); ++k) {
    const NumericMatrix fk = fixef[k];
    for (int c = 0; c < fk.ncol(); ++c) {
      for (int i = 0; i < n; ++i) {
        stacked[i + static_cast<std::size_t>(n) * (col_off[k] + c)] = fk(i, c);
      }
    }
  }
  std::vector<double> mean(p_tot);
  for (int a = 0; a < p_tot; ++a) {
    long double s = 0.0;
    for (int i = 0; i < n; ++i) {
      s += stacked[i + static_cast<std::size_t>(n) * a];
    }
    mean[a] = static_cast<double>(s / n);
    sweep_mean(a, mi) = mean[a];
  }
  const R_xlen_t base = static_cast<R_xlen_t>(p_tot) * p_tot * mi;
  for (int a = 0; a < p_tot; ++a) {
    for (int c = a; c < p_tot; ++c) {
      double v = NA_REAL;
      if (n > 1) {
        long double s = 0.0;
        for (int i = 0; i < n; ++i) {
          s += (stacked[i + static_cast<std::size_t>(n) * a] - mean[a]) *
               (stacked[i + static_cast<std::size_t>(n) * c] - mean[c]);
        }
        v = static_cast<double>(s / (n - 1));
      }
      sweep_cov[base + a + static_cast<R_xlen_t>(p_tot) * c] = v;
      sweep_cov[base + c + static_cast<R_xlen_t>(p_tot) * a] = v;
    }
  }
}

}  // namespace

List two_block_sweep(const TwoBlockIndexSpec& spec_in, const NumericVector& b_in,
                     const List& fixef_in, const NumericMatrix& iters_in,
                     const NumericVector& iters_ranef_in, int inner_sweeps,
                     bool use_mode, bool progbar, const std::string& stage_label,
                     double rwt, SweepObserver* observer, bool group_major,
                     bool group_check) {
  if (!(rwt >= 0.0 && rwt < 1.0)) {
    stop("rwt must be in [0, 1); the ridge term needed at rwt = 1 is not implemented yet.");
  }
  // spec_in.wt1 shares memory with the R object spec$wt1, and assigning a
  // same-length sugar expression to an Rcpp vector writes into its existing
  // memory. Scale a clone and rebind to it instead.
  TwoBlockIndexSpec spec = spec_in;
  NumericVector wt_scaled = clone(spec_in.wt1);
  for (R_xlen_t r = 0; r < wt_scaled.size(); ++r) {
    wt_scaled[r] *= (1.0 - rwt);
  }
  spec.wt1 = wt_scaled;

  const IntegerVector dim = b_in.attr("dim");
  if (dim.size() != 3) {
    stop("b must be a J x p_re x n array.");
  }
  const int J = dim[0];
  const int p_re = dim[1];
  const int n = dim[2];
  if (J != spec.J || p_re != spec.p_re) {
    stop("b does not match spec (J x p_re).");
  }
  SweepObserver no_observer;
  SweepObserver& obs = observer ? *observer : no_observer;

  NumericVector b = clone(b_in);
  List fixef(p_re);
  for (int k = 0; k < p_re; ++k) {
    fixef[k] = clone(as<NumericMatrix>(fixef_in[k]));
  }
  NumericMatrix iters = clone(iters_in);
  NumericVector iters_ranef = clone(iters_ranef_in);

  std::vector<int> col_off(p_re);
  int p_tot = 0;
  for (int k = 0; k < p_re; ++k) {
    col_off[k] = p_tot;
    p_tot += as<NumericMatrix>(fixef[k]).ncol();
  }
  NumericMatrix sweep_mean(p_tot, inner_sweeps);
  NumericVector sweep_cov(static_cast<R_xlen_t>(p_tot) * p_tot * inner_sweeps);
  sweep_cov.attr("dim") = IntegerVector::create(p_tot, p_tot, inner_sweeps);
  std::vector<double> stacked(static_cast<std::size_t>(n) * p_tot);

  std::vector<GroupData> gdata;
  if (group_major) gdata = groups_data(spec);

  for (int m = 1; m <= inner_sweeps; ++m) {
    checkUserInterrupt();
    const int mi = m - 1;
    const std::string prefix_groups =
        sweep_progbar_prefix(stage_label, m, inner_sweeps, "RE");
    const std::string prefix_population =
        sweep_progbar_prefix(stage_label, m, inner_sweeps, "fixef");

    if (group_major) {
      std::vector<NumericMatrix> mu_D;
      const NumericVector mu_all = groups_prior_mean_all(spec, fixef, n, &mu_D);

      // Replayed chains keep their coef and skip the group loop.
      std::vector<NumericMatrix> coef(n);
      std::vector<double> it1(n, 0.0);
      std::vector<int> active;
      active.reserve(n);
      for (int i = 0; i < n; ++i) {
        if (!obs.groups_source(mi, i, coef[i], it1[i])) {
          coef[i] = NumericMatrix(J, spec.l1);
          active.push_back(i);
        }
      }

      std::vector<double> iters_sum(n, 0.0);
      if (!active.empty()) {
        for (int j = 0; j < J; ++j) {
          const GroupChainsDraw d1 =
              group_chains_draw(spec, gdata[j], mu_all, j, active, group_check);
          const NumericMatrix& cj = use_mode ? d1.coef_mode : d1.coef;
          for (int c = 0; c < static_cast<int>(active.size()); ++c) {
            const int i = active[c];
            for (int a = 0; a < spec.l1; ++a) {
              coef[i](j, a) = cj(c, a);
            }
            iters_sum[i] += d1.iters[c];
          }
          if (progbar) {
            lmecore::progress_bar(j + 1.0, static_cast<double>(J),
                                  prefix_groups);
          }
        }
        for (int i : active) {
          it1[i] = iters_sum[i] / J;
        }
      }
      if (progbar) {
        lmecore::progress_bar_finish(false);
      }

      for (int i = 0; i < n; ++i) {
        obs.groups_seen(mi, i, mu_D[i], coef[i], it1[i]);
        const NumericMatrix b_i = groups_place(spec, coef[i]);
        for (int k = 0; k < p_re; ++k) {
          for (int j = 0; j < J; ++j) {
            b[j + J * (k + p_re * i)] = b_i(j, k);
          }
        }
        iters_ranef[i] = iters_ranef[i] + it1[i];
      }
    } else {
      for (int i = 0; i < n; ++i) {
        const NumericMatrix mu_D = groups_prior_mean(spec, fixef, i);
        NumericMatrix coef;
        double it1 = 0.0;
        if (!obs.groups_source(mi, i, coef, it1)) {
          const GroupsDraw d1 = groups_draw(spec, mu_D);
          coef = use_mode ? d1.coef_mode : d1.coef;
          it1 = d1.iters;
        }
        obs.groups_seen(mi, i, mu_D, coef, it1);
        const NumericMatrix b_i = groups_place(spec, coef);
        for (int k = 0; k < p_re; ++k) {
          for (int j = 0; j < J; ++j) {
            b[j + J * (k + p_re * i)] = b_i(j, k);
          }
        }
        iters_ranef[i] = iters_ranef[i] + it1;
        if (progbar) {
          lmecore::progress_bar(i + 1.0, static_cast<double>(n),
                                           prefix_groups);
        }
      }
      if (progbar) {
        lmecore::progress_bar_finish(false);
      }
    }

    for (int i = 0; i < n; ++i) {
      NumericMatrix b_i(J, p_re);
      for (int k = 0; k < p_re; ++k) {
        for (int j = 0; j < J; ++j) {
          b_i(j, k) = b[j + J * (k + p_re * i)];
        }
      }
      for (int k = 0; k < p_re; ++k) {
        const NumericVector y_k = population_response(spec, b_i, k);
        NumericVector coef;
        if (!obs.population_source(mi, i, k, coef)) {
          const PopulationDraw d2 = population_draw(spec, y_k, k);
          coef = use_mode ? d2.coef_mode : d2.coef;
        }
        obs.population_seen(mi, i, k, y_k, coef);
        NumericMatrix fk = fixef[k];
        population_write(spec, fk, i, k, coef);
        iters(i, k) = iters(i, k) + 1.0;
      }
      if (progbar) {
        lmecore::progress_bar(i + 1.0, static_cast<double>(n),
                                         prefix_population);
      }
    }
    if (progbar) {
      lmecore::progress_bar_finish(m == inner_sweeps);
    }
    sweep_summary(fixef, col_off, n, p_tot, mi, stacked, sweep_mean, sweep_cov);
    if (progbar && n <= 1) {
      lmecore::progress_bar(
          static_cast<double>(m), static_cast<double>(inner_sweeps),
          sweep_progbar_prefix(stage_label, m, inner_sweeps, ""));
      lmecore::progress_bar_finish(m == inner_sweeps);
    }
  }

  return List::create(
      Named("b") = b,
      Named("fixef") = fixef,
      Named("iters") = iters,
      Named("iters_ranef") = iters_ranef,
      Named("sweep_mean") = sweep_mean,
      Named("sweep_cov") = sweep_cov);
}

}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes
