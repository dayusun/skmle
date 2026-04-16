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

  // Quadrature points
  const vec *lq_x;
  const vec *lq_w;

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
    if ((*data->delta)[i] == 1.0 && (*data->kerval)[i] > 0) {
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
    double tt = 0.5 * (*data->lq_x)[q] + 0.5;
    double weight = 0.5 * (*data->lq_w)[q];

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
List skmle_cpp_fit(int n, int p, int gammap, double s, double h,
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
                 double h, const arma::mat &covariates, const arma::mat &bsmat,
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

    if (delta[i] == 1.0 && kerval[i] > 0) {
      double S0_sum = 0.0;
      arma::vec S1_sum = arma::zeros<arma::vec>(p);
      arma::mat S2_sum = arma::zeros<arma::mat>(p, p);

      for (int k = 0; k < n; ++k) {
        if (X[i] <= X[k]) {
          double dist_XX = X[i] - obs_times[k];
          if (dist_XX > 0) {
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
                 double h, const arma::mat &covariates, const arma::mat &bsmat,
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

    if (delta[i] == 1.0 && kerval[i] > 0 && (X[i] - obs_times[i]) > 0) {
      double S0_sum = 0.0;
      arma::vec S1_sum = arma::zeros<arma::vec>(p);
      for (int k = 0; k < n; ++k) {
        if (X[i] <= X[k]) {
          double dist_XX = X[i] - obs_times[k];
          if (dist_XX > 0) {
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

    double tt = 0.5 * lq_x[q] + 0.5;
    double weight = 0.5 * lq_w[q];
    double alpha_tt = dot(trans(bsmat_tt_all.row(q)), gamma);

    // The S0/S1 aggregates over `k` depend only on `q` (via tt and alpha_tt),
    // not on the reference index `i` or covariate index `j`. Hoist the
    // `k`-loop out of the i/j loops and precompute once per quadrature node.
    double S0_tt_q = 0.0;
    arma::vec S1_tt_q_vec = arma::zeros<arma::vec>(p);
    for (int k = 0; k < n; ++k) {
      if (tt <= X[k]) {
        double dist_k = tt - obs_times[k];
        if (dist_k > 0) {
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

// Utility to generate kerfun
inline double calc_kerfun(double dist, double h) {
  if (dist > 0) {
    double res = (1.0 - std::pow(dist / h, 2.0)) * 0.75;
    if (res > 0)
      return res / h;
  }
  return 0.0;
}

// [[Rcpp::export]]
arma::vec skmle_cv_cpp(int n, int p, int gammap, double s,
                       const arma::vec &h_grid, int K, const arma::vec &fold_id,
                       const arma::vec &id_vec, const arma::mat &covariates,
                       const arma::mat &bsmat, const arma::vec &X,
                       const arma::vec &obs_times, const arma::vec &delta,
                       const arma::vec &lq_x, const arma::vec &lq_w,
                       const arma::mat &bsmat_tt_all, int maxeval,
                       double xtol_rel, bool quiet) {

  int n_h = h_grid.n_elem;
  arma::vec cv_losses = arma::zeros<arma::vec>(n_h);
  int n_obs = X.n_elem;
  int n_quad = lq_x.n_elem;

  // Quadrature points
  arma::vec tts(n_quad);
  for (int q = 0; q < n_quad; ++q) {
    tts[q] = 0.5 * lq_x[q] + 0.5;
  }

  for (int hi = 0; hi < n_h; ++hi) {
    Rcpp::checkUserInterrupt();

    double h = h_grid[hi];
    if (!quiet) {
      Rcout << "Evaluating bandwidth h = " << h << "\n";
    }
    double loss_sum = 0.0;

    for (int k = 1; k <= K; ++k) {
      Rcpp::checkUserInterrupt();

      // Find indices for test and train
      std::vector<int> test_subj;
      std::vector<int> train_subj;

      // fold_id is 1-indexed, length equals number of unique subjects
      // but id_vec links each observation to a subject (1-indexed)
      // We assume fold_id aligns with unique subject IDs 1..n

      std::vector<int> test_idx_v;
      std::vector<int> train_idx_v;

      for (int i = 0; i < n_obs; ++i) {
        int subj_id = static_cast<int>(std::round(id_vec[i]));
        int fold_for_subj = static_cast<int>(std::round(fold_id[subj_id - 1]));
        if (fold_for_subj == k) {
          test_idx_v.push_back(i);
        } else {
          train_idx_v.push_back(i);
        }
      }

      int n_train_obs = train_idx_v.size();
      int n_test_obs = test_idx_v.size();

      // Count unique subjects in train fold
      std::vector<int> unique_train_subjs;
      for (int idx : train_idx_v) {
        int subj_id = static_cast<int>(std::round(id_vec[idx]));
        if (std::find(unique_train_subjs.begin(), unique_train_subjs.end(),
                      subj_id) == unique_train_subjs.end()) {
          unique_train_subjs.push_back(subj_id);
        }
      }
      int n_train_subjs = unique_train_subjs.size();

      // Arrays for train data
      arma::mat Z_train(n_train_obs, p);
      arma::mat bsmat_train(n_train_obs, gammap);
      arma::vec X_train(n_train_obs);
      arma::vec obs_train(n_train_obs);
      arma::vec delta_train(n_train_obs);
      arma::vec kerval_train(n_train_obs);
      arma::mat kerval_tt_train(n_quad, n_train_obs);

      for (int i = 0; i < n_train_obs; ++i) {
        int orig_idx = train_idx_v[i];
        Z_train.row(i) = covariates.row(orig_idx);
        bsmat_train.row(i) = bsmat.row(orig_idx);
        X_train[i] = X[orig_idx];
        obs_train[i] = obs_times[orig_idx];
        delta_train[i] = delta[orig_idx];

        double dist = X_train[i] - obs_train[i];
        kerval_train[i] = calc_kerfun(dist, h);

        for (int q = 0; q < n_quad; ++q) {
          double dist_tt = tts[q] - obs_train[i];
          kerval_tt_train(q, i) = calc_kerfun(dist_tt, h);
        }
      }

      // Inequality constraints for train
      std::vector<int> ineq_rows;
      if (s != 0.0) {
        for (int i = 0; i < n_train_obs; ++i) {
          double dist = X_train[i] - obs_train[i];
          if (dist > 0 && dist <= h) {
            ineq_rows.push_back(i);
          }
        }
      }

      arma::mat ineqmat_train(ineq_rows.size(), p + gammap);
      for (size_t r = 0; r < ineq_rows.size(); ++r) {
        int tr_idx = ineq_rows[r];
        for (int j = 0; j < p; ++j)
          ineqmat_train(r, j) = Z_train(tr_idx, j);
        for (int j = 0; j < gammap; ++j)
          ineqmat_train(r, p + j) = bsmat_train(tr_idx, j);
      }

      // Fit model on train fold
      skmle_data train_data = {n_train_subjs,
                               p,
                               gammap,
                               s,
                               h,
                               &Z_train,
                               &bsmat_train,
                               &X_train,
                               &obs_train,
                               &delta_train,
                               &kerval_train,
                               &lq_x,
                               &lq_w,
                               &bsmat_tt_all,
                               &kerval_tt_train,
                               &ineqmat_train};

      int n_vars = p + gammap;
      std::vector<double> x_est(n_vars, 0.0);

      NloptOptPtr opt(nlopt_create(NLOPT_LD_SLSQP, n_vars));
      if (s != 0.0 && ineqmat_train.n_rows > 0) {
        std::vector<double> tol(ineqmat_train.n_rows, 1e-8);
        nlopt_add_inequality_mconstraint(opt.get(), ineqmat_train.n_rows,
                                         ineq_constraints, &train_data,
                                         tol.data());
      }
      nlopt_set_min_objective(opt.get(), nll_obj, &train_data);
      nlopt_set_xtol_rel(opt.get(), xtol_rel);
      nlopt_set_maxeval(opt.get(), maxeval);

      double minf;
      nlopt_result train_res =
          nlopt_optimize(opt.get(), x_est.data(), &minf);

      // If the training fit failed we cannot trust the test-fold NLL:
      // let this bandwidth lose the grid outright.
      if (train_res < 0) {
        loss_sum = R_PosInf;
        break;
      }

      // Extract beta and gamma
      arma::vec beta_est(x_est.data(), p);
      arma::vec gamma_est(x_est.data() + p, gammap);

      // Arrays for test data
      // Count unique subjects in test fold
      std::vector<int> unique_test_subjs;
      for (int idx : test_idx_v) {
        int subj_id = static_cast<int>(std::round(id_vec[idx]));
        if (std::find(unique_test_subjs.begin(), unique_test_subjs.end(),
                      subj_id) == unique_test_subjs.end()) {
          unique_test_subjs.push_back(subj_id);
        }
      }
      int n_test_subjs = unique_test_subjs.size();

      arma::mat Z_test(n_test_obs, p);
      arma::mat bsmat_test(n_test_obs, gammap);
      arma::vec X_test(n_test_obs);
      arma::vec obs_test(n_test_obs);
      arma::vec delta_test(n_test_obs);
      arma::vec kerval_test(n_test_obs);
      arma::mat kerval_tt_test(n_quad, n_test_obs);

      for (int i = 0; i < n_test_obs; ++i) {
        int orig_idx = test_idx_v[i];
        Z_test.row(i) = covariates.row(orig_idx);
        bsmat_test.row(i) = bsmat.row(orig_idx);
        X_test[i] = X[orig_idx];
        obs_test[i] = obs_times[orig_idx];
        delta_test[i] = delta[orig_idx];

        double dist = X_test[i] - obs_test[i];
        kerval_test[i] = calc_kerfun(dist, h);

        for (int q = 0; q < n_quad; ++q) {
          double dist_tt = tts[q] - obs_test[i];
          kerval_tt_test(q, i) = calc_kerfun(dist_tt, h);
        }
      }

      arma::mat empty_ineq_test;
      skmle_data test_data = {n_test_subjs,
                              p,
                              gammap,
                              s,
                              h,
                              &Z_test,
                              &bsmat_test,
                              &X_test,
                              &obs_test,
                              &delta_test,
                              &kerval_test,
                              &lq_x,
                              &lq_w,
                              &bsmat_tt_all,
                              &kerval_tt_test,
                              &empty_ineq_test};

      double nll_val = nll_obj(n_vars, x_est.data(), nullptr, &test_data);
      loss_sum += nll_val;
    }

    cv_losses[hi] = loss_sum / static_cast<double>(K);
  }

  return cv_losses;
}
