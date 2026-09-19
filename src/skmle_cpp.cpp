#include <RcppArmadillo.h>
#include <algorithm>
#include <memory>
#include <nloptrAPI.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Box-Cox transform g(x) = (1 + s*x)^(1/s) (s > 0), exp(x) (s = 0).
// Outside the feasibility region (base = 1 + s*x <= 0) the hazard is
// mathematically 0 but the NLL depends on log(t_val). We therefore floor
// the *value* at 2.22e-16 (machine eps): log(2.22e-16) ~= -36 becomes a
// soft infeasibility penalty that SLSQP can see and move away from,
// rather than a hard excluded region. Callers still guard log() with
// `t_val > 0` so a future change to this floor cannot produce -Inf.
static constexpr double kBoundaryValueFloor = 2.220446e-16;

inline double trans_fun(double x, double s) {
  if (s == 0.0) return std::exp(x);
  const double base = s * x + 1.0;
  if (base <= 0.0) return kBoundaryValueFloor;
  return std::pow(base, 1.0 / s);
}

// The derivative functions below are not logged; they feed into gradient
// and sandwich-variance sums. In the infeasibility region:
//   * g'(x)   -> 0     for 0 < s <= 1
//   * g'/g    -> +Inf  (all s > 0)
//   * g''/g   -> +/-Inf
// Returning the 0 limit is safe and consistent with the ineqmat hard
// constraints that actually gate feasibility at the data points; the
// previous 2.22e-16 floor silently understated derivatives that truly
// diverge, biasing SLSQP's gradient search for s > 1.
inline double trans_fun_d(double x, double s) {
  if (s == 0.0) return std::exp(x);
  const double base = s * x + 1.0;
  if (base <= 0.0) return 0.0;
  return std::pow(base, 1.0 / s - 1.0);
}

inline double trans_fun_d1o1(double x, double s) {
  if (s == 0.0) return 1.0;
  const double base = s * x + 1.0;
  if (base <= 0.0) return 0.0;
  return 1.0 / base;
}

inline double trans_fun_d12o1(double x, double s) {
  if (s == 0.0) return std::exp(x);
  const double base = s * x + 1.0;
  if (base <= 0.0) return 0.0;
  return std::pow(base, 1.0 / s - 2.0);
}

// RAII wrapper for nlopt_opt so the optimizer is destroyed on any
// unwind path (including Rcpp interrupt exceptions raised inside
// `Rcpp::checkUserInterrupt()`).
struct NloptOptDeleter {
  void operator()(nlopt_opt opt) const noexcept {
    if (opt != nullptr) nlopt_destroy(opt);
  }
};
using NloptOptPtr =
    std::unique_ptr<std::remove_pointer<nlopt_opt>::type, NloptOptDeleter>;

// Data structure to pass to nlopt
struct skmle_data {
  int n;
  int p;      // number of covariates
  int gammap; // number of spline basis functions
  double s;
  double h;
  const mat *covariates;
  const mat *bsmat;
  const vec *X;
  const vec *obs_times;
  const vec *delta;
  const vec *kerval;

  // Quadrature points, on [-1, 1]; mapped to [0, tau] where tau is the end of
  // follow-up.  Fixing tau at 1 would silently integrate the cumulative hazard
  // over the wrong interval whenever the data are not on the unit scale.
  const vec *lq_x;
  const vec *lq_w;
  double tau;

  // Matrices pre-computed for quadrature points
  const mat *bsmat_tt_all;
  const mat *kerval_tt_all; // n_quad x n matrix

  // Matrix for inequality constraints
  const mat *ineqmat;
};

// Objective function for NLOPT
double nll_obj(unsigned n_vars, const double *x, double *grad,
               void *my_func_data) {
  Rcpp::checkUserInterrupt();

  skmle_data *data = (skmle_data *)my_func_data;

  vec beta(const_cast<double *>(x), data->p, false);
  vec gamma(const_cast<double *>(x + data->p), data->gammap, false);

  // 1. loglik1_1 and loglik1_1_d
  vec alphaX = (*data->bsmat) * gamma;
  vec inner1 = alphaX + (*data->covariates) * beta;

  double res1 = 0.0;
  vec d_beta = zeros<vec>(data->p);
  vec d_gamma = zeros<vec>(data->gammap);
  int n_obs = data->X->n_elem;

  for (int i = 0; i < n_obs; ++i) {
    // Guard on nonzero, not positive.  A weight function supplied by a caller
    // is not required to be nonnegative, and a `> 0` test drops those rows
    // silently rather than failing.  The objective, its gradient, and the
    // bread and meat in calc_A/calc_B must all sum over the SAME rows or they
    // describe different estimators, so the four guards move together.  The
    // kernels shipped here are nonnegative, for which the two forms agree.
    if ((*data->delta)[i] == 1.0 && (*data->kerval)[i] != 0) {
      double t_val = trans_fun(inner1[i], data->s);
      if (t_val > 0)
        res1 += std::log(t_val) * (*data->kerval)[i];

      if (grad) {
        double d_val = trans_fun_d1o1(inner1[i], data->s) * (*data->delta)[i] *
                       (*data->kerval)[i];
        d_beta += d_val * trans((*data->covariates).row(i)); // column vector
        d_gamma += d_val * trans((*data->bsmat).row(i));     // column vector
      }
    }
  }

  // 2. loglik2_inner_1 and loglik2_inner_1_d
  double res2 = 0.0;
  vec d2_beta = zeros<vec>(data->p);
  vec d2_gamma = zeros<vec>(data->gammap);

  int n_quad = data->lq_x->n_elem;

  for (int q = 0; q < n_quad; ++q) {
    double tt = 0.5 * data->tau * ((*data->lq_x)[q] + 1.0);
    double weight = 0.5 * data->tau * (*data->lq_w)[q];

    // alpha_tt is 1x1 vector (scalar) because tt is scalar
    double alpha_tt = dot(trans((*data->bsmat_tt_all).row(q)), gamma);

    double q_res2 = 0.0;
    vec q_d2_beta = zeros<vec>(data->p);
    vec q_d2_gamma = zeros<vec>(data->gammap);

    for (int i = 0; i < n_obs; ++i) {
      if (tt < (*data->X)[i]) {
        double k_tt = (*data->kerval_tt_all)(q, i);
        if (k_tt > 0) {
          double in_val = alpha_tt + dot((*data->covariates).row(i), beta);
          double t_val = trans_fun(in_val, data->s);
          q_res2 += t_val * k_tt; // sum over i

          if (grad) {
            double td_val = trans_fun_d(in_val, data->s) * k_tt;
            q_d2_beta += td_val * trans((*data->covariates).row(i));
            q_d2_gamma += td_val * trans((*data->bsmat_tt_all).row(q));
          }
        }
      }
    }

    res2 += weight * q_res2;
    if (grad) {
      d2_beta += weight * q_d2_beta;
      d2_gamma += weight * q_d2_gamma;
    }
  }

  if (grad) {
    for (int j = 0; j < data->p; ++j) {
      grad[j] = -(d_beta[j] - d2_beta[j]) / data->n;
    }
    for (int j = 0; j < data->gammap; ++j) {
      grad[data->p + j] = -(d_gamma[j] - d2_gamma[j]) / data->n;
    }
  }

  return -(res1 - res2) / data->n;
}

// Inequality constraints for s != 0
void ineq_constraints(unsigned m, double *result, unsigned n_vars,
                      const double *x, double *grad, void *my_func_data) {
  skmle_data *data = (skmle_data *)my_func_data;

  vec beta(const_cast<double *>(x), data->p, false);
  vec gamma(const_cast<double *>(x + data->p), data->gammap, false);

  // ineqmat is (num constraints) x (p + gammap)
  // Constraint: -ineqmat * x - 1 / s <= 0
  int num_constraints = data->ineqmat->n_rows;

  for (int i = 0; i < num_constraints; ++i) {
    double val = 0.0;
    for (int j = 0; j < data->p; ++j)
      val += (*data->ineqmat)(i, j) * x[j];
    for (int j = 0; j < data->gammap; ++j)
      val += (*data->ineqmat)(i, data->p + j) * x[data->p + j];

    result[i] = -val - 1.0 / data->s;

    if (grad) {
      for (unsigned int k = 0; k < n_vars; ++k) {
        grad[i * n_vars + k] = -(*data->ineqmat)(i, k);
      }
    }
  }
}

// One fit, shared by the R-facing export and the cross-validation loop below.
//
// Factored out rather than copied.  The CV loop needs exactly this optimiser
// set-up, and a second copy of it is a second place for the objective, the
// constraint tolerance or the algorithm to drift away from the one the user
// gets from skmle().
struct FitResult {
  int status;
  double minimum;
  std::vector<double> solution;
};

static FitResult fit_core(int n, int p, int gammap, double s, double h,
                          double tau, const arma::mat &covariates,
                          const arma::mat &bsmat, const arma::vec &X,
                          const arma::vec &obs_times, const arma::vec &delta,
                          const arma::vec &kerval, const arma::vec &lq_x,
                          const arma::vec &lq_w, const arma::mat &bsmat_tt_all,
                          const arma::mat &kerval_tt_all,
                          const arma::mat &ineqmat, int maxeval,
                          double xtol_rel) {

  skmle_data data = {n,
                     p,
                     gammap,
                     s,
                     h,
                     &covariates,
                     &bsmat,
                     &X,
                     &obs_times,
                     &delta,
                     &kerval,
                     &lq_x,
                     &lq_w,
                     tau,
                     &bsmat_tt_all,
                     &kerval_tt_all,
                     &ineqmat};

  int n_vars = p + gammap;
  std::vector<double> x(n_vars, 0.0);

  NloptOptPtr opt(nlopt_create(NLOPT_LD_SLSQP, n_vars));
  if (s != 0.0 && ineqmat.n_rows > 0) {
    std::vector<double> tol(ineqmat.n_rows, 1e-8);
    nlopt_add_inequality_mconstraint(opt.get(), ineqmat.n_rows,
                                     ineq_constraints, &data, tol.data());
  }

  nlopt_set_min_objective(opt.get(), nll_obj, &data);
  nlopt_set_xtol_rel(opt.get(), xtol_rel);
  nlopt_set_maxeval(opt.get(), maxeval);

  double minf = 0.0;
  nlopt_result res = nlopt_optimize(opt.get(), x.data(), &minf);

  return FitResult{static_cast<int>(res), minf, x};
}

// [[Rcpp::export]]
List skmle_cpp_fit(int n, int p, int gammap, double s, double h, double tau,
                   const arma::mat &covariates, const arma::mat &bsmat,
                   const arma::vec &X, const arma::vec &obs_times,
                   const arma::vec &delta, const arma::vec &kerval,
                   const arma::vec &lq_x, const arma::vec &lq_w,
                   const arma::mat &bsmat_tt_all,
                   const arma::mat &kerval_tt_all, const arma::mat &ineqmat,
                   int maxeval, double xtol_rel) {

  FitResult fr =
      fit_core(n, p, gammap, s, h, tau, covariates, bsmat, X, obs_times, delta,
               kerval, lq_x, lq_w, bsmat_tt_all, kerval_tt_all, ineqmat,
               maxeval, xtol_rel);

  return List::create(Named("status") = fr.status,
                      Named("minimum") = fr.minimum,
                      Named("solution") = fr.solution);
}

// [[Rcpp::export]]
arma::mat calc_A(const arma::vec &beta, const arma::vec &gamma, double s,
                 double h, bool one_sided, const arma::mat &covariates,
                 const arma::mat &bsmat,
                 const arma::vec &X, const arma::vec &obs_times,
                 const arma::vec &delta, const arma::vec &kerval,
                 const arma::mat &bsmat_XX, int n_subj) {
  int n = X.n_elem;
  int p = covariates.n_cols;
  arma::mat A_est = arma::zeros<arma::mat>(p, p);

  arma::vec alpha_XX = bsmat_XX * gamma;
  arma::vec cov_beta = covariates * beta;

  for (int i = 0; i < n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();

    // Nonzero, not positive: see the note in nll_obj.
    if (delta[i] == 1.0 && kerval[i] != 0) {
      double S0_sum = 0.0;
      arma::vec S1_sum = arma::zeros<arma::vec>(p);
      arma::mat S2_sum = arma::zeros<arma::mat>(p, p);

      for (int k = 0; k < n; ++k) {
        if (X[i] <= X[k]) {
          double dist_XX = X[i] - obs_times[k];
          if (!one_sided || dist_XX > 0) {
            double kerval_XX =
                std::max((1 - std::pow(dist_XX / h, 2)) * 0.75, 0.0) / h;
            if (kerval_XX > 0) {
              double inner = alpha_XX[k] + cov_beta[k];
              double term = trans_fun_d12o1(inner, s) * kerval_XX;
              S0_sum += term;
              S1_sum += term * trans(covariates.row(k));
              S2_sum += term * (trans(covariates.row(k)) * covariates.row(k));
            }
          }
        }
      }

      if (S0_sum > 0) {
        arma::vec S1 = S1_sum / S0_sum;
        arma::mat S2 = S2_sum / S0_sum;
        double inner_i = alpha_XX[i] + cov_beta[i];
        double t_val = trans_fun_d1o1(inner_i, s);
        arma::mat outerprod = (S2 - S1 * trans(S1)) * std::pow(t_val, 2);
        A_est += outerprod * kerval[i];
      }
    }
  }
  return A_est / n_subj;
}

// [[Rcpp::export]]
arma::mat calc_B(const arma::vec &beta, const arma::vec &gamma, double s,
                 double h, double tau, bool one_sided, const arma::mat &covariates,
                 const arma::mat &bsmat,
                 const arma::vec &X, const arma::vec &obs_times,
                 const arma::vec &delta, const arma::vec &kerval,
                 const arma::vec &id, const arma::mat &bsmat_XX,
                 const arma::vec &lq_x, const arma::vec &lq_w,
                 const arma::mat &bsmat_tt_all, const arma::mat &kerval_tt_all,
                 int n_subj) {
  int n = X.n_elem;
  int p = covariates.n_cols;

  arma::vec alpha_XX = bsmat_XX * gamma;
  arma::vec cov_beta = covariates * beta;

  // use vectors instead of maps keyed by id
  std::vector<arma::vec> id_to_bb1(n_subj, arma::zeros<arma::vec>(p));
  std::vector<int> id_counts(n_subj, 0);
  for (int i = 0; i < n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();

    int idx = static_cast<int>(std::round(id[i])) - 1;
    if (idx < 0 || idx >= n_subj)
      continue;
    id_counts[idx] += 1;

    // Nonzero WEIGHT, positive LAG.  The two tests are unrelated: the lag test
    // is the one-sided restriction and stays `> 0`; the weight test only asks
    // whether the row contributes at all.  See the note in nll_obj.
    if (delta[i] == 1.0 && kerval[i] != 0 && (X[i] - obs_times[i]) > 0) {
      double S0_sum = 0.0;
      arma::vec S1_sum = arma::zeros<arma::vec>(p);
      for (int k = 0; k < n; ++k) {
        if (X[i] <= X[k]) {
          double dist_XX = X[i] - obs_times[k];
          if (!one_sided || dist_XX > 0) {
            double kerval_XX =
                std::max((1 - std::pow(dist_XX / h, 2)) * 0.75, 0.0) / h;
            if (kerval_XX > 0) {
              double inner = alpha_XX[k] + cov_beta[k];
              double term = trans_fun_d12o1(inner, s) * kerval_XX;
              S0_sum += term;
              S1_sum += term * trans(covariates.row(k));
            }
          }
        }
      }
      if (S0_sum > 0) {
        arma::vec S1 = S1_sum / S0_sum;
        double inner_i = alpha_XX[i] + cov_beta[i];
        double t_val = trans_fun_d1o1(inner_i, s);
        id_to_bb1[idx] += (S1 - trans(covariates.row(i))) * t_val * kerval[i];
      }
    }
  }

  int n_quad = lq_x.n_elem;
  std::vector<arma::vec> id_to_bb2(n_subj, arma::zeros<arma::vec>(p));

  for (int q = 0; q < n_quad; ++q) {
    Rcpp::checkUserInterrupt();

    double tt = 0.5 * tau * (lq_x[q] + 1.0);
    double weight = 0.5 * tau * lq_w[q];
    double alpha_tt = dot(trans(bsmat_tt_all.row(q)), gamma);

    // The S0/S1 aggregates over `k` depend only on `q` (via tt and alpha_tt),
    // not on the reference index `i` or covariate index `j`. Hoist the
    // `k`-loop out of the i/j loops and precompute once per quadrature node.
    double S0_tt_q = 0.0;
    arma::vec S1_tt_q_vec = arma::zeros<arma::vec>(p);
    for (int k = 0; k < n; ++k) {
      if (tt <= X[k]) {
        double dist_k = tt - obs_times[k];
        if (!one_sided || dist_k > 0) {
          double k_tt_k =
              std::max((1 - std::pow(dist_k / h, 2)) * 0.75, 0.0) / h;
          if (k_tt_k > 0) {
            double in_val = alpha_tt + dot(covariates.row(k), beta);
            double td = trans_fun_d12o1(in_val, s) * k_tt_k;
            S0_tt_q += td;
            S1_tt_q_vec += td * trans(covariates.row(k));
          }
        }
      }
    }
    if (S0_tt_q <= 0) continue;

    arma::vec S1_tt_q = S1_tt_q_vec / S0_tt_q;

    // per-subject r1 contribution for this quadrature node
    std::vector<arma::vec> id_to_r1_vec(n_subj, arma::zeros<arma::vec>(p));
    for (int i = 0; i < n; ++i) {
      int idx = static_cast<int>(std::round(id[i])) - 1;
      if (idx < 0 || idx >= n_subj) continue;
      if (tt <= X[i]) {
        double k_tt = kerval_tt_all(q, i);
        if (k_tt > 0) {
          double in_val = alpha_tt + dot(covariates.row(i), beta);
          double factor = trans_fun_d(in_val, s) * k_tt;
          id_to_r1_vec[idx] +=
              factor * (S1_tt_q - trans(covariates.row(i)));
        }
      }
    }

    for (int subj = 0; subj < n_subj; ++subj) {
      int count = id_counts[subj];
      if (count > 0) {
        id_to_bb2[subj] +=
            weight * id_to_r1_vec[subj] / static_cast<double>(count);
      }
    }
  }

  arma::mat B_est = arma::zeros<arma::mat>(p, p);
  for (int subj = 0; subj < n_subj; ++subj) {
    arma::vec bb_diff = id_to_bb1[subj] - id_to_bb2[subj];
    B_est += bb_diff * trans(bb_diff);
  }

  return B_est / n_subj;
}

// [[Rcpp::export]]
double skmle_eval_nll_cpp(int n, int p, int gammap, double s, double h,
                          double tau,
                          const arma::vec &beta, const arma::vec &gamma,
                          const arma::mat &covariates, const arma::mat &bsmat,
                          const arma::vec &X, const arma::vec &obs_times,
                          const arma::vec &delta, const arma::vec &kerval,
                          const arma::vec &lq_x, const arma::vec &lq_w,
                          const arma::mat &bsmat_tt_all,
                          const arma::mat &kerval_tt_all) {

  arma::mat empty_ineq;
  skmle_data data = {n,
                     p,
                     gammap,
                     s,
                     h,
                     &covariates,
                     &bsmat,
                     &X,
                     &obs_times,
                     &delta,
                     &kerval,
                     &lq_x,
                     &lq_w,
                     tau,
                     &bsmat_tt_all,
                     &kerval_tt_all,
                     &empty_ineq};

  int n_vars = p + gammap;
  std::vector<double> x(n_vars, 0.0);
  for (int i = 0; i < p; ++i)
    x[i] = beta[i];
  for (int j = 0; j < gammap; ++j)
    x[p + j] = gamma[j];

  return nll_obj(n_vars, x.data(), nullptr, &data);
}


// ---------------------------------------------------------------------------
// Cross-validation
// ---------------------------------------------------------------------------
//
// The whole grid runs here rather than in R for one reason: the held-out score
// evaluates the model's transformation, and evaluating it in R meant a second
// implementation of trans_fun() living in R/utils.R, with nothing comparing the
// two and no test that could, because trans_fun is not exported.  Scoring here
// means the criterion and the objective call the same function.
//
// The weight is an R function and is evaluated in R, once per candidate
// bandwidth, vectorised over every lag at once.  It is NOT reimplemented here.
// A previous version described the weight to C++ as polynomial coefficients so
// the loop could rebuild it, which meant the package held two descriptions of
// every kernel and needed a test to keep them agreeing.  One R call per
// candidate costs nothing next to K fits and leaves exactly one description.
// It also means any R function works, including weights that are not
// polynomial, not symmetric and not nonnegative.

// Held-out negative log-likelihood per held-out subject.
//
// No kernel and no bandwidth appear here.  The pieces that depend on the data
// alone -- the basis at the quadrature nodes, the quadrature weights, the row
// each node reads its carried-forward covariate from -- are built once in R by
// locf_score() and passed in; only beta and gamma change per candidate.
static double locf_nll_cpp(const arma::mat &node_bs, const arma::vec &node_w,
                           const arma::uvec &node_row,
                           const arma::uvec &node_subj, const arma::mat &ev_bs,
                           const arma::uvec &ev_row, const arma::vec &ev_delta,
                           int n_subj, const arma::mat &Z,
                           const arma::vec &beta, const arma::vec &gamma,
                           double s, const arma::uvec &keep) {

  arma::vec cumhaz = arma::zeros<arma::vec>(n_subj);
  if (node_bs.n_rows > 0) {
    arma::vec lin = node_bs * gamma + Z.rows(node_row) * beta;
    for (arma::uword i = 0; i < lin.n_elem; ++i) {
      cumhaz[node_subj[i]] += node_w[i] * trans_fun(lin[i], s);
    }
  }

  arma::vec evlin = ev_bs * gamma + Z.rows(ev_row) * beta;

  double acc = 0.0;
  for (arma::uword j = 0; j < keep.n_elem; ++j) {
    const arma::uword m = keep[j];
    acc += ev_delta[m] * std::log(trans_fun(evlin[m], s)) - cumhaz[m];
  }
  return -acc / static_cast<double>(keep.n_elem);
}

// [[Rcpp::export]]
arma::mat skmle_cv_cpp(int p, int gammap, double s, double tau,
                       const arma::vec &h_grid, int K,
                       const arma::uvec &fold_id_subj, const arma::uvec &id_vec,
                       const arma::mat &covariates, const arma::mat &bsmat,
                       const arma::vec &X, const arma::vec &obs_times,
                       const arma::vec &delta, const arma::vec &lq_x,
                       const arma::vec &lq_w, const arma::mat &bsmat_tt_all,
                       const arma::vec &tts, const arma::mat &node_bs,
                       const arma::vec &node_w, const arma::uvec &node_row,
                       const arma::uvec &node_subj, const arma::mat &ev_bs,
                       const arma::uvec &ev_row, const arma::vec &ev_delta,
                       int n_subj, int maxeval, double xtol_rel, bool quiet,
                       Rcpp::Function weight, bool one_sided) {

  const int n_h = h_grid.n_elem;
  const int n_row = X.n_elem;
  const int n_quad = tts.n_elem;

  // Zbs is the constraint design, needed only when s != 0.
  arma::mat Zbs;
  if (s != 0.0) Zbs = arma::join_rows(covariates, bsmat);

  // Held-out subject codes per fold, and the row-level training mask.
  std::vector<arma::uvec> keep_by_fold(K);
  for (int k = 0; k < K; ++k) {
    keep_by_fold[k] = arma::find(fold_id_subj == static_cast<arma::uword>(k + 1));
  }

  arma::mat out(2, n_h);

  for (int hi = 0; hi < n_h; ++hi) {
    Rcpp::checkUserInterrupt();
    const double h = h_grid[hi];
    if (!quiet) Rcpp::Rcout << "Evaluating bandwidth h = " << h << "\n";

    // Two R calls, both vectorised over every lag: the row lags, and the
    // quadrature-node-by-row lag matrix.  The one-sided restriction is applied
    // here rather than inside the weight, because it is the risk-set rule of
    // the model and not a property of the kernel.
    const arma::vec row_lag = X - obs_times;
    arma::vec kerval = Rcpp::as<arma::vec>(weight(Rcpp::wrap(row_lag / h))) / h;
    if (kerval.n_elem != static_cast<arma::uword>(n_row)) {
      Rcpp::stop("'weight' returned %d values for %d lags; it must return one "
                 "value per element of its argument",
                 static_cast<int>(kerval.n_elem), n_row);
    }
    if (one_sided) kerval %= arma::conv_to<arma::vec>::from(row_lag > 0);

    arma::mat node_lag(n_quad, n_row);
    for (int q = 0; q < n_quad; ++q) {
      for (int i = 0; i < n_row; ++i) node_lag(q, i) = tts[q] - obs_times[i];
    }
    arma::mat kerval_tt =
        Rcpp::as<arma::mat>(weight(Rcpp::wrap(node_lag / h))) / h;
    if (kerval_tt.n_rows != static_cast<arma::uword>(n_quad) ||
        kerval_tt.n_cols != static_cast<arma::uword>(n_row)) {
      Rcpp::stop("'weight' did not preserve the shape of its argument: given a "
                 "%dx%d matrix it returned %dx%d",
                 n_quad, n_row, static_cast<int>(kerval_tt.n_rows),
                 static_cast<int>(kerval_tt.n_cols));
    }
    if (one_sided) kerval_tt %= arma::conv_to<arma::mat>::from(node_lag > 0);

    // Same support the fit uses, so the constraint set matches the objective.
    arma::uvec feasible;
    if (s != 0.0) {
      feasible = arma::zeros<arma::uvec>(n_row);
      for (int i = 0; i < n_row; ++i) {
        const double d = X[i] - obs_times[i];
        feasible[i] = (std::fabs(d) <= h && (!one_sided || d > 0)) ? 1u : 0u;
      }
    }

    arma::vec fold_losses(K);

    for (int k = 0; k < K; ++k) {
      arma::uvec tr(n_row);
      for (int i = 0; i < n_row; ++i) {
        tr[i] = (fold_id_subj[id_vec[i]] != static_cast<arma::uword>(k + 1)) ? 1u : 0u;
      }
      const arma::uvec tr_idx = arma::find(tr);

      arma::mat ineqmat(0, p + gammap);
      if (s != 0.0) {
        const arma::uvec ci = arma::find(tr % feasible);
        if (!ci.is_empty()) ineqmat = Zbs.rows(ci);
      }

      // The objective divides by the number of training SUBJECTS, not rows.
      const arma::uvec tr_subj = arma::unique(id_vec.elem(tr_idx));

      FitResult fr = fit_core(
          static_cast<int>(tr_subj.n_elem), p, gammap, s, h, tau,
          covariates.rows(tr_idx), bsmat.rows(tr_idx), X.elem(tr_idx),
          obs_times.elem(tr_idx), delta.elem(tr_idx), kerval.elem(tr_idx), lq_x,
          lq_w, bsmat_tt_all, kerval_tt.cols(tr_idx), ineqmat, maxeval,
          xtol_rel);

      // A training fit that failed cannot produce a trustworthy held-out
      // score: let this bandwidth lose the grid outright.
      if (fr.status < 0) {
        fold_losses[k] = arma::datum::inf;
        continue;
      }

      const arma::vec sol(fr.solution);
      fold_losses[k] = locf_nll_cpp(
          node_bs, node_w, node_row, node_subj, ev_bs, ev_row, ev_delta,
          n_subj, covariates, sol.head(p), sol.tail(gammap), s,
          keep_by_fold[k]);
    }

    out(0, hi) = arma::mean(fold_losses);
    // stats::sd() is the n-1 denominator; arma::stddev defaults to the same.
    out(1, hi) = arma::stddev(fold_losses) / std::sqrt(static_cast<double>(K));
  }

  return out;
}

// Thin wrapper so the test suite can check the held-out score directly,
// against integrate() on the same carried-forward path.  The scorer is not
// reachable from R otherwise, and a criterion nothing can call is a criterion
// nothing can check.
//
// [[Rcpp::export]]
double locf_nll_cpp_r(const arma::mat &node_bs, const arma::vec &node_w,
                      const arma::uvec &node_row, const arma::uvec &node_subj,
                      const arma::mat &ev_bs, const arma::uvec &ev_row,
                      const arma::vec &ev_delta, int n_subj,
                      const arma::mat &Z, const arma::vec &beta,
                      const arma::vec &gamma, double s,
                      const arma::uvec &keep) {
  return locf_nll_cpp(node_bs, node_w, node_row, node_subj, ev_bs, ev_row,
                      ev_delta, n_subj, Z, beta, gamma, s, keep);
}
