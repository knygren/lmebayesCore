// Indexed two-block sweep loop. R twin: .two_block_sweep_r() in
// R/two_block_sweep_indexed.R; same loop order and progress output when
// group_major is false. The default group-major groups block loops over
// groups outside and chains inside, with progress counted in groups.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "glmbayesCore_calls.h"
#include "simfuncs_twoblock.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>

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

// lmecore::progress_bar rewrites the whole line and flushes per character,
// so redraw only when the shown percentage or bar length changes (40 dots,
// as in glmbayesCore's progress_bar), and always at x = N.
class ThrottledBar {
 public:
  ThrottledBar(bool on, double N, const std::string& prefix)
      : on_(on && N > 0.0), N_(N), prefix_(prefix) {}
  void update(double x) {
    if (!on_) return;
    const double f = x / N_;
    const long key = std::lround(100.0 * f) * 64 + std::lround(40.0 * f);
    if (key == last_ && x < N_) return;
    last_ = key;
    lmecore::progress_bar(x, N_, prefix_);
  }

 private:
  bool on_;
  double N_;
  std::string prefix_;
  long last_ = -1;
};

// Fold one named rNormalGLMGroupChains vector (setup_timing, mode_compare)
// into `acc` by name: max_* entries take the maximum, all others are summed.
void named_add(NumericVector& acc, const NumericVector& x) {
  if (x.size() == 0) return;
  if (acc.size() == 0) {
    acc = clone(x);
    return;
  }
  const CharacterVector nm = x.names();
  for (R_xlen_t i = 0; i < x.size() && i < acc.size(); ++i) {
    const std::string field = as<std::string>(nm[i]);
    if (field.compare(0, 4, "max_") == 0) {
      if (!(x[i] <= acc[i])) acc[i] = x[i];
    } else {
      acc[i] += x[i];
    }
  }
}

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
                     bool group_check, int mode_method, int setup_method,
                     int warm_method, int hess_method) {
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

  const int pop_method = population_draw_method();
  std::vector<PopulationFactor> pop_f;
  if (pop_method != POP_REG) {
    pop_f.reserve(p_re);
    for (int k = 0; k < p_re; ++k) pop_f.push_back(population_factor(spec, k));
  }
  PopulationCompare pop_cmp;

  // Seconds over the whole call: block-1 parts from rNormalGLMGroupChains
  // (group-major only), then block-1 and block-2 wall time.
  double t_parts[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double t_block1 = 0.0, t_block2 = 0.0;
  // rNormalGLMGroupChains setup_timing and mode_compare, combined over calls
  // (see header).
  NumericVector setup_parts;
  NumericVector mode_parts;
  NumericVector setup_cmp;
  NumericVector mode_cmp;
  NumericVector warm_cmp;
  NumericVector hess_cmp;
  NumericVector newton_cmp;
  NumericVector sim_cmp;
  NumericVector eval_cmp;

  // Previous mode, prior mean and mode Hessian per (chain i, group j), at
  // index i * J + j; filled only when warm starts are on.
  const bool warm = group_major && warm_method != 0;
  const int l1 = spec.l1;
  const std::size_t n_slots = warm ? static_cast<std::size_t>(n) * J : 0;
  std::vector<char> warm_has(n_slots, 0);
  std::vector<double> warm_b(n_slots * l1);
  std::vector<double> warm_mu(n_slots * l1);
  std::vector<double> warm_A1(n_slots * l1 * l1);

  for (int m = 1; m <= inner_sweeps; ++m) {
    checkUserInterrupt();
    const int mi = m - 1;
    const std::string prefix_groups =
        sweep_progbar_prefix(stage_label, m, inner_sweeps, "RE");
    const std::string prefix_population =
        sweep_progbar_prefix(stage_label, m, inner_sweeps, "fixef");

    auto t_start = std::chrono::steady_clock::now();
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
        const int na = static_cast<int>(active.size());
        ThrottledBar bar(progbar, static_cast<double>(J), prefix_groups);
        for (int j = 0; j < J; ++j) {
          NumericMatrix start(0, 0);
          if (warm) {
            start = NumericMatrix(l1, na);
            std::fill(start.begin(), start.end(), NA_REAL);
            const arma::mat Pj(const_cast<double*>(gdata[j].P.begin()), l1, l1,
                               false);
            for (int c = 0; c < na; ++c) {
              const std::size_t s = static_cast<std::size_t>(active[c]) * J + j;
              if (!warm_has[s]) continue;
              const double* mu_new =
                  &mu_all[static_cast<R_xlen_t>(l1) *
                          (j + static_cast<R_xlen_t>(J) * active[c])];
              arma::vec dmu(l1);
              for (int a = 0; a < l1; ++a) {
                dmu[a] = mu_new[a] - warm_mu[s * l1 + a];
              }
              const arma::mat A1(&warm_A1[s * l1 * l1], l1, l1, false);
              arma::vec step;
              if (!arma::solve(step, A1, Pj * dmu,
                               arma::solve_opts::no_approx)) {
                continue;
              }
              for (int a = 0; a < l1; ++a) {
                start(a, c) = warm_b[s * l1 + a] + step[a];
              }
            }
          }
          const GroupChainsDraw d1 =
              group_chains_draw(spec, gdata[j], mu_all, j, active, group_check,
                                mode_method, setup_method, start, warm_method,
                                hess_method);
          const NumericMatrix& cj = use_mode ? d1.coef_mode : d1.coef;
          for (int p = 0; p < 5 && p < d1.timing.size(); ++p) {
            t_parts[p] += d1.timing[p];
          }
          named_add(setup_parts, d1.setup_timing);
          named_add(mode_parts, d1.mode_timing);
          named_add(setup_cmp, d1.setup_compare);
          named_add(mode_cmp, d1.mode_compare);
          named_add(warm_cmp, d1.warm_compare);
          named_add(hess_cmp, d1.hessian_compare);
          named_add(newton_cmp, d1.newton_compare);
          named_add(sim_cmp, d1.sim_compare);
          named_add(eval_cmp, d1.eval_compare);
          if (warm && d1.mode_hessian.nrow() == l1 * l1 &&
              d1.mode_hessian.ncol() == na) {
            for (int c = 0; c < na; ++c) {
              const std::size_t s = static_cast<std::size_t>(active[c]) * J + j;
              const double* mu_new =
                  &mu_all[static_cast<R_xlen_t>(l1) *
                          (j + static_cast<R_xlen_t>(J) * active[c])];
              for (int a = 0; a < l1; ++a) {
                warm_b[s * l1 + a] = d1.coef_mode(c, a);
                warm_mu[s * l1 + a] = mu_new[a];
              }
              for (int k = 0; k < l1 * l1; ++k) {
                warm_A1[s * l1 * l1 + k] = d1.mode_hessian(k, c);
              }
              warm_has[s] = 1;
            }
          }
          for (int c = 0; c < static_cast<int>(active.size()); ++c) {
            const int i = active[c];
            for (int a = 0; a < spec.l1; ++a) {
              coef[i](j, a) = cj(c, a);
            }
            iters_sum[i] += d1.iters[c];
          }
          bar.update(j + 1.0);
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
      ThrottledBar bar(progbar, static_cast<double>(n), prefix_groups);
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
        bar.update(i + 1.0);
      }
      if (progbar) {
        lmecore::progress_bar_finish(false);
      }
    }
    t_block1 += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t_start).count();

    t_start = std::chrono::steady_clock::now();
    ThrottledBar pop_bar(progbar, static_cast<double>(n), prefix_population);
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
          const PopulationDraw d2 =
              pop_method == POP_CHOL
                  ? population_draw_chol(pop_f[k], y_k)
                  : pop_method == POP_COMPARE
                        ? population_draw_compare(pop_cmp, spec, pop_f[k],
                                                  y_k, k)
                        : population_draw(spec, y_k, k);
          coef = use_mode ? d2.coef_mode : d2.coef;
        }
        obs.population_seen(mi, i, k, y_k, coef);
        NumericMatrix fk = fixef[k];
        population_write(spec, fk, i, k, coef);
        iters(i, k) = iters(i, k) + 1.0;
      }
      pop_bar.update(i + 1.0);
    }
    if (progbar) {
      lmecore::progress_bar_finish(m == inner_sweeps);
    }
    t_block2 += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t_start).count();
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
      Named("sweep_cov") = sweep_cov,
      Named("timing") = NumericVector::create(
          Named("mode") = t_parts[0], Named("setup") = t_parts[1],
          Named("eval") = t_parts[2], Named("finalize") = t_parts[3],
          Named("simulate") = t_parts[4], Named("block1") = t_block1,
          Named("block2") = t_block2),
      Named("setup_timing") = setup_parts,
      Named("mode_timing") = mode_parts,
      Named("setup_compare") = setup_cmp,
      Named("mode_compare") = mode_cmp,
      Named("warm_compare") = warm_cmp,
      Named("hessian_compare") = hess_cmp,
      Named("newton_compare") = newton_cmp,
      Named("sim_compare") = sim_cmp,
      Named("eval_compare") = eval_cmp,
      Named("pop_compare") = NumericVector::create(
          Named("calls") = pop_cmp.calls, Named("flagged") = pop_cmp.flagged,
          Named("max_mean") = pop_cmp.max_mean,
          Named("max_draw") = pop_cmp.max_draw,
          Named("t_reg") = pop_cmp.t_reg, Named("t_chol") = pop_cmp.t_chol));
}

}  // namespace twoblock
}  // namespace sim
}  // namespace glmbayes
