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

// [[Rcpp::export]]
List skmle_cpp_fit(int n, int p, int gammap, double s, double h, double tau,
                   const arma::mat &covariates, const arma::mat &bsmat,
                   const arma::vec &X, const arma::vec &obs_times,
                   const arma::vec &delta, const arma::vec &kerval,
                   const arma::vec &lq_x, const arma::vec &lq_w,
                   const arma::mat &bsmat_tt_all,
                   const arma::mat &kerval_tt_all, const arma::mat &ineqmat,
                   int maxeval, double xtol_rel) {

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

  double minf;
  nlopt_result res = nlopt_optimize(opt.get(), x.data(), &minf);

  return List::create(Named("status") = static_cast<int>(res),
                      Named("minimum") = minf, Named("solution") = x);
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

