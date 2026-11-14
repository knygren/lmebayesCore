// C++ port of the portable two-block sweep in R/rGLMM_sweep_indexed.R.
//
// Each port_* function has an R twin (.two_block_port_*) with the same
// inputs, the same output shape, and the same loop order. All index vectors
// in `spec` are 1-based, as in R; they are shifted only at the point of use.
// The draws come from the same samplers the R steps call through their
// exports: rNormalGLMBlocks (Block 1) and rNormalReg (Block 2).
//
// State layout, identical to R:
//   b            numeric array  J x p_re x n (column-major)
//   fixef        list of p_re numeric matrices, n x p_k
//   iters        numeric matrix n x p_re
//   iters_ranef  numeric vector n

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "progress_utils.h"
#include "simfuncs.h"

using namespace Rcpp;

namespace {

// R: .two_block_port_block1_prior_mean(spec, fixef, i); `i` is 0-based here.
NumericMatrix port_block1_prior_mean(const List& spec, const List& fixef, int i) {
  const NumericMatrix D = spec["D"];
  const List obs_rows = spec["obs_rows"];
  const IntegerVector d_col = spec["d_col"];
  const List W = spec["W"];
  const List w_for_group = spec["w_for_group"];
  const List w_col = spec["w_col"];
  const int l1 = D.ncol();
  const int J = obs_rows.size();

  NumericMatrix mu_D(l1, J);
  for (int k = 0; k < d_col.size(); ++k) {
    const NumericMatrix Wk = W[k];
    const IntegerVector wg = w_for_group[k];
    const IntegerVector wc = w_col[k];
    const NumericMatrix fk = fixef[k];
    NumericVector s(J);
    for (int c = 0; c < Wk.ncol(); ++c) {
      const double g = fk(i, wc[c] - 1);
      for (int j = 0; j < J; ++j) {
        s[j] = s[j] + Wk(wg[j] - 1, c) * g;
      }
    }
    for (int j = 0; j < J; ++j) {
      mu_D(d_col[k] - 1, j) = s[j];
    }
  }
  return mu_D;
}

// R: .two_block_port_block1_draw(spec, mu_D)
List port_block1_draw(const List& spec, const NumericMatrix& mu_D) {
  const NumericVector y = spec["y"];
  const NumericMatrix D = spec["D"];
  const NumericVector offset1 = spec["offset1"];
  const NumericVector wt1 = spec["wt1"];
  const NumericVector disp1 = spec["disp1"];
  const NumericMatrix P_D = spec["P_D"];
  const List obs_rows = spec["obs_rows"];
  const Function f2 = spec["f2"];
  const Function f3 = spec["f3"];
  const std::string family = as<std::string>(spec["family"]);
  const std::string link = as<std::string>(spec["link"]);
  const int gridtype = as<int>(spec["gridtype"]);
  const int n_envopt = as<int>(spec["n_envopt"]);

  // The R callbacks f2/f3 enter their own (non-nesting) RNGScope and reload
  // .Random.seed; this scope writes the RNG state back after every draw, as
  // the R step's .Call does, so later calls do not reuse the same stream.
  RNGScope rng_scope;
  List out = glmbayes::sim::rNormalGLMBlocks(
    1, y, D, offset1, wt1, disp1, mu_D, List::create(P_D), true, obs_rows,
    f2, f3, family, link, gridtype, n_envopt, false, false, false
  );

  const List gr = out["group_results"];
  long double s = 0.0;
  for (int b = 0; b < gr.size(); ++b) {
    const List g = gr[b];
    const NumericVector it = g["iters"];
    s += it[0];
  }
  return List::create(
    Named("coef") = out["coefficients"],
    Named("coef_mode") = out["coef.mode"],
    Named("iters") = static_cast<double>(s / gr.size())
  );
}

// R: .two_block_port_block1_place(spec, coef)
NumericMatrix port_block1_place(const List& spec, const NumericMatrix& coef) {
  const IntegerVector d_col = spec["d_col"];
  const int J = coef.nrow();
  const int p_re = d_col.size();
  NumericMatrix out(J, p_re);
  for (int k = 0; k < p_re; ++k) {
    for (int j = 0; j < J; ++j) {
      out(j, k) = coef(j, d_col[k] - 1);
    }
  }
  return out;
}

// R: .two_block_port_block2_response(spec, b_i, k); `k` is 0-based here.
NumericVector port_block2_response(const List& spec, const NumericMatrix& b_i, int k) {
  const List b_for_wrow = spec["b_for_wrow"];
  const IntegerVector rows = b_for_wrow[k];
  NumericVector y_k(rows.size());
  for (int r = 0; r < rows.size(); ++r) {
    y_k[r] = b_i(rows[r] - 1, k);
  }
  return y_k;
}

// R: .two_block_port_block2_draw(spec, y_k, k); `k` is 0-based here.
List port_block2_draw(const List& spec, const NumericVector& y_k, int k) {
  const List W = spec["W"];
  const List mu2 = spec["mu2"];
  const List P2 = spec["P2"];
  const NumericVector disp2 = spec["disp2"];
  const NumericMatrix Wk = W[k];
  const NumericVector mu = mu2[k];
  const NumericMatrix P = P2[k];
  const Function f2g = spec["f2_gauss"];
  const Function f3g = spec["f3_gauss"];
  const int gridtype = as<int>(spec["gridtype"]);
  const int n_w = y_k.size();

  // See port_block1_draw: one RNGScope per draw, as for the R step's .Call.
  RNGScope rng_scope;
  List out = glmbayes::sim::rNormalReg(
    1, y_k, Wk, mu, P, NumericVector(n_w, 0.0), NumericVector(n_w, 1.0),
    disp2[k], f2g, f3g, mu, "gaussian", "identity", gridtype
  );
  const NumericMatrix draw = out["coefficients"];
  const NumericVector mode = out["coef.mode"];
  NumericVector coef(draw.ncol());
  for (int c = 0; c < draw.ncol(); ++c) {
    coef[c] = draw(0, c);
  }
  return List::create(
    Named("coef") = coef,
    Named("coef_mode") = NumericVector(mode.begin(), mode.end())
  );
}

// R: .two_block_port_block2_write(spec, fixef_k, i, k, coef); `i`, `k` 0-based.
void port_block2_write(const List& spec, NumericMatrix& fixef_k, int i, int k,
                       const NumericVector& coef) {
  const List w_col = spec["w_col"];
  const IntegerVector wc = w_col[k];
  for (int c = 0; c < wc.size(); ++c) {
    fixef_k(i, wc[c] - 1) = coef[c];
  }
}

// R: .two_block_progbar_prefix(stage_label, m, inner_sweeps, phase), with
// phase "" giving the per-sweep prefix used when there is a single chain.
std::string port_progbar_prefix(const std::string& stage_label, int m,
                                int inner_sweeps, const std::string& phase) {
  std::string out = stage_label.empty() ? std::string("")
                                        : std::string("[") + stage_label + "] ";
  out += "sweep " + std::to_string(m) + "/" + std::to_string(inner_sweeps);
  if (!phase.empty()) {
    out += " " + phase;
  }
  return out + ": ";
}

NumericVector port_array(int d1, int d2, int d3, int d4 = -1) {
  const int len = d1 * d2 * d3 * (d4 < 0 ? 1 : d4);
  NumericVector out(len, 0.0);
  if (d4 < 0) {
    out.attr("dim") = IntegerVector::create(d1, d2, d3);
  } else {
    out.attr("dim") = IntegerVector::create(d1, d2, d3, d4);
  }
  return out;
}

// R: .two_block_port_log_new(spec, n, M)
List port_log_new(const List& spec, int n, int M) {
  const NumericMatrix D = spec["D"];
  const List obs_rows = spec["obs_rows"];
  const List W = spec["W"];
  const int l1 = D.ncol();
  const int J = obs_rows.size();
  List y(W.size());
  List block2(W.size());
  for (int k = 0; k < W.size(); ++k) {
    const NumericMatrix Wk = W[k];
    y[k] = port_array(Wk.nrow(), n, M);
    block2[k] = port_array(Wk.ncol(), n, M);
  }
  return List::create(
    Named("mu_D") = port_array(l1, J, n, M),
    Named("block1") = port_array(J, l1, n, M),
    Named("iters1") = NumericMatrix(n, M),
    Named("y") = y,
    Named("block2") = block2
  );
}

// R: .two_block_port_sweep_r(spec, state, inner_sweeps, use_mode, ...)
// `m` is 1-based in the loop, as in R; array offsets use m - 1.
List port_sweep(const List& spec, const NumericVector& b_in, const List& fixef_in,
                const NumericMatrix& iters_in, const NumericVector& iters_ranef_in,
                int inner_sweeps, bool use_mode, bool progbar,
                const std::string& stage_label, Nullable<List> tape_in,
                bool record) {
  const IntegerVector dim = b_in.attr("dim");
  if (dim.size() != 3) {
    stop("b must be a J x p_re x n array.");
  }
  const int J = dim[0];
  const int p_re = dim[1];
  const int n = dim[2];
  const NumericMatrix D = spec["D"];
  const int l1 = D.ncol();

  NumericVector b = clone(b_in);
  List fixef(p_re);
  for (int k = 0; k < p_re; ++k) {
    fixef[k] = clone(as<NumericMatrix>(fixef_in[k]));
  }
  NumericMatrix iters = clone(iters_in);
  NumericVector iters_ranef = clone(iters_ranef_in);
  const bool use_tape = tape_in.isNotNull();
  List tape = use_tape ? List(tape_in) : List();
  NumericVector tape_b1;
  NumericMatrix tape_it1;
  List tape_b2;
  if (use_tape) {
    tape_b1 = tape["block1"];
    tape_it1 = as<NumericMatrix>(tape["iters1"]);
    tape_b2 = tape["block2"];
  }
  List log = record ? port_log_new(spec, n, inner_sweeps) : List();
  const int J1 = as<List>(spec["obs_rows"]).size();
  NumericVector prior_sum = port_array(l1, J1, inner_sweeps);
  NumericVector prior_sumsq = port_array(l1, J1, inner_sweeps);
  NumericVector log_mu, log_b1;
  NumericMatrix log_it1;
  List log_y, log_b2;
  if (record) {
    log_mu = log["mu_D"];
    log_b1 = log["block1"];
    log_it1 = as<NumericMatrix>(log["iters1"]);
    log_y = log["y"];
    log_b2 = log["block2"];
  }

  for (int m = 1; m <= inner_sweeps; ++m) {
    checkUserInterrupt();
    const int mi = m - 1;
    const std::string prefix_b1 = port_progbar_prefix(stage_label, m, inner_sweeps, "RE");
    const std::string prefix_b2 = port_progbar_prefix(stage_label, m, inner_sweeps, "fixef");
    for (int i = 0; i < n; ++i) {
      const NumericMatrix mu_D = port_block1_prior_mean(spec, fixef, i);
      for (int j = 0; j < J1; ++j) {
        for (int a = 0; a < l1; ++a) {
          const int idx = a + l1 * (j + J1 * mi);
          prior_sum[idx] = prior_sum[idx] + mu_D(a, j);
          prior_sumsq[idx] = prior_sumsq[idx] + mu_D(a, j) * mu_D(a, j);
        }
      }
      NumericMatrix coef;
      double it1 = 0.0;
      if (!use_tape) {
        const List d1 = port_block1_draw(spec, mu_D);
        coef = use_mode ? as<NumericMatrix>(d1["coef_mode"])
                        : as<NumericMatrix>(d1["coef"]);
        it1 = as<double>(d1["iters"]);
      } else {
        coef = NumericMatrix(J, l1);
        for (int c = 0; c < l1; ++c) {
          for (int j = 0; j < J; ++j) {
            coef(j, c) = tape_b1[j + J * (c + l1 * (i + n * mi))];
          }
        }
        it1 = tape_it1(i, mi);
      }
      if (record) {
        for (int j = 0; j < J; ++j) {
          for (int a = 0; a < l1; ++a) {
            log_mu[a + l1 * (j + J * (i + n * mi))] = mu_D(a, j);
          }
        }
        for (int c = 0; c < l1; ++c) {
          for (int j = 0; j < J; ++j) {
            log_b1[j + J * (c + l1 * (i + n * mi))] = coef(j, c);
          }
        }
        log_it1(i, mi) = it1;
      }
      const NumericMatrix b_i = port_block1_place(spec, coef);
      for (int k = 0; k < p_re; ++k) {
        for (int j = 0; j < J; ++j) {
          b[j + J * (k + p_re * i)] = b_i(j, k);
        }
      }
      iters_ranef[i] = iters_ranef[i] + it1;
      if (progbar) {
        glmbayes::progress::progress_bar(i + 1.0, static_cast<double>(n), prefix_b1);
      }
    }
    if (progbar) {
      glmbayes::progress::progress_bar_finish(false);
    }
    for (int i = 0; i < n; ++i) {
      NumericMatrix b_i(J, p_re);
      for (int k = 0; k < p_re; ++k) {
        for (int j = 0; j < J; ++j) {
          b_i(j, k) = b[j + J * (k + p_re * i)];
        }
      }
      for (int k = 0; k < p_re; ++k) {
        const NumericVector y_k = port_block2_response(spec, b_i, k);
        const List w_col = spec["w_col"];
        const int p_k = as<IntegerVector>(w_col[k]).size();
        NumericVector coef;
        if (!use_tape) {
          const List d2 = port_block2_draw(spec, y_k, k);
          coef = use_mode ? as<NumericVector>(d2["coef_mode"])
                          : as<NumericVector>(d2["coef"]);
        } else {
          const NumericVector t2 = tape_b2[k];
          coef = NumericVector(p_k);
          for (int c = 0; c < p_k; ++c) {
            coef[c] = t2[c + p_k * (i + n * mi)];
          }
        }
        if (record) {
          NumericVector ly = log_y[k];
          NumericVector lb2 = log_b2[k];
          const int r_k = y_k.size();
          for (int r = 0; r < r_k; ++r) {
            ly[r + r_k * (i + n * mi)] = y_k[r];
          }
          for (int c = 0; c < p_k; ++c) {
            lb2[c + p_k * (i + n * mi)] = coef[c];
          }
        }
        NumericMatrix fk = fixef[k];
        port_block2_write(spec, fk, i, k, coef);
        iters(i, k) = iters(i, k) + 1.0;
      }
      if (progbar) {
        glmbayes::progress::progress_bar(i + 1.0, static_cast<double>(n), prefix_b2);
      }
    }
    if (progbar) {
      glmbayes::progress::progress_bar_finish(m == inner_sweeps);
    }
    if (progbar && n <= 1) {
      glmbayes::progress::progress_bar(static_cast<double>(m),
                                       static_cast<double>(inner_sweeps),
                                       port_progbar_prefix(stage_label, m,
                                                           inner_sweeps, ""));
      glmbayes::progress::progress_bar_finish(m == inner_sweeps);
    }
  }

  return List::create(
    Named("b") = b,
    Named("fixef") = fixef,
    Named("iters") = iters,
    Named("iters_ranef") = iters_ranef,
    Named("log") = record ? RObject(log) : RObject(R_NilValue),
    Named("prior_sum") = prior_sum,
    Named("prior_sumsq") = prior_sumsq
  );
}

}  // namespace

// Exports for step-by-step comparison with the R twins. `i` and `k` are
// 1-based, as in R.

// [[Rcpp::export]]
NumericMatrix two_block_port_block1_prior_mean_cpp_export(const List& spec,
                                                          const List& fixef,
                                                          int i) {
  return port_block1_prior_mean(spec, fixef, i - 1);
}

// [[Rcpp::export]]
List two_block_port_block1_draw_cpp_export(const List& spec,
                                           const NumericMatrix& mu_D) {
  return port_block1_draw(spec, mu_D);
}

// [[Rcpp::export]]
NumericMatrix two_block_port_block1_place_cpp_export(const List& spec,
                                                     const NumericMatrix& coef) {
  return port_block1_place(spec, coef);
}

// [[Rcpp::export]]
NumericVector two_block_port_block2_response_cpp_export(const List& spec,
                                                        const NumericMatrix& b_i,
                                                        int k) {
  return port_block2_response(spec, b_i, k - 1);
}

// [[Rcpp::export]]
List two_block_port_block2_draw_cpp_export(const List& spec,
                                           const NumericVector& y_k, int k) {
  return port_block2_draw(spec, y_k, k - 1);
}

// [[Rcpp::export]]
NumericMatrix two_block_port_block2_write_cpp_export(const List& spec,
                                                     const NumericMatrix& fixef_k,
                                                     int i, int k,
                                                     const NumericVector& coef) {
  NumericMatrix out = clone(fixef_k);
  port_block2_write(spec, out, i - 1, k - 1, coef);
  return out;
}

// [[Rcpp::export]]
List two_block_port_sweep_cpp_export(const List& spec, const NumericVector& b,
                                     const List& fixef, const NumericMatrix& iters,
                                     const NumericVector& iters_ranef,
                                     int inner_sweeps, bool use_mode,
                                     bool progbar, const std::string& stage_label,
                                     Nullable<List> tape = R_NilValue,
                                     bool record = false) {
  return port_sweep(spec, b, fixef, iters, iters_ranef, inner_sweeps, use_mode,
                    progbar, stage_label, tape, record);
}
