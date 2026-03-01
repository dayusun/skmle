#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <algorithm> // for std::count and other utilities

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Utility functions
inline double trans_fun(double x, double s) {
  if (s == 0.0) return exp(x);
  if (x < -1.0 / s) return 2.220446e-16; // Machine epsilon approx
  return std::pow(s * x + 1.0, 1.0 / s);
}

inline double trans_fun_d(double x, double s) {
  if (s == 0.0) return exp(x);
  if (x < -1.0 / s) return 2.220446e-16;
  return std::pow(s * x + 1.0, 1.0 / s - 1.0);
}

inline double trans_fun_d1o1(double x, double s) {
  if (s == 0.0) return 1.0;
  if (x < -1.0 / s) return 2.220446e-16;
  return std::pow(s * x + 1.0, -1.0);
}

inline double trans_fun_d12o1(double x, double s) {
  if (s == 0.0) return exp(x);
  if (x < -1.0 / s) return 2.220446e-16;
  return std::pow(s * x + 1.0, 1.0 / s - 2.0);
}

// Data structure to pass to nlopt
struct skmle_data {
  int n;
  int p; // number of covariates
  int gammap; // number of spline basis functions
  double s;
  double h;
  const mat* covariates;
  const mat* bsmat;
  const vec* X;
  const vec* obs_times;
  const vec* delta;
  const vec* kerval;
  
  // Quadrature points
  const vec* lq_x;
  const vec* lq_w;
  
  // Matrices pre-computed for quadrature points
  const mat* bsmat_tt_all;
  const mat* kerval_tt_all; // n_quad x n matrix
  
  // Matrix for inequality constraints
  const mat* ineqmat;
};


// Objective function for NLOPT
double nll_obj(unsigned n_vars, const double* x, double* grad, void* my_func_data) {
  skmle_data* data = (skmle_data*)my_func_data;
  
  vec beta(const_cast<double*>(x), data->p, false);
  vec gamma(const_cast<double*>(x + data->p), data->gammap, false);
  
  // 1. loglik1_1 and loglik1_1_d
  vec alphaX = (*data->bsmat) * gamma;
  vec inner1 = alphaX + (*data->covariates) * beta;
  
  double res1 = 0.0;
  vec d_beta = zeros<vec>(data->p);
  vec d_gamma = zeros<vec>(data->gammap);
  int n_obs = data->X->n_elem;
  
  for (int i = 0; i < n_obs; ++i) {
    if ( (*data->delta)[i] == 1.0 && (*data->kerval)[i] > 0 ) {
      double t_val = trans_fun(inner1[i], data->s);
      if (t_val > 0) res1 += std::log(t_val) * (*data->kerval)[i];
      
      if (grad) {
        double d_val = trans_fun_d1o1(inner1[i], data->s) * (*data->delta)[i] * (*data->kerval)[i];
        d_beta += d_val * trans((*data->covariates).row(i));  // column vector
        d_gamma += d_val * trans((*data->bsmat).row(i));      // column vector
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
    if(grad) {
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
void ineq_constraints(unsigned m, double *result, unsigned n_vars, const double* x, double* grad, void* my_func_data) {
  skmle_data* data = (skmle_data*)my_func_data;
  
  vec beta(const_cast<double*>(x), data->p, false);
  vec gamma(const_cast<double*>(x + data->p), data->gammap, false);
  
  // ineqmat is (num constraints) x (p + gammap)
  // Constraint: -ineqmat * x - 1 / s <= 0
  int num_constraints = data->ineqmat->n_rows;
  
  for (int i = 0; i < num_constraints; ++i) {
    double val = 0.0;
    for (int j = 0; j < data->p; ++j) val += (*data->ineqmat)(i, j) * x[j];
    for (int j = 0; j < data->gammap; ++j) val += (*data->ineqmat)(i, data->p + j) * x[data->p + j];
    
    result[i] = -val - 1.0 / data->s;
    
    if (grad) {
      for (int k = 0; k < n_vars; ++k) {
        grad[i * n_vars + k] = -(*data->ineqmat)(i, k);
      }
    }
  }
}

// [[Rcpp::export]]
List skmle_cpp_fit(int n, int p, int gammap, double s, double h,
                   const arma::mat& covariates,
                   const arma::mat& bsmat,
                   const arma::vec& X,
                   const arma::vec& obs_times,
                   const arma::vec& delta,
                   const arma::vec& kerval,
                   const arma::vec& lq_x,
                   const arma::vec& lq_w,
                   const arma::mat& bsmat_tt_all,
                   const arma::mat& kerval_tt_all,
                   const arma::mat& ineqmat,
                   int maxeval,
                   double xtol_rel) {
  
  skmle_data data = {n, p, gammap, s, h, 
                     &covariates, &bsmat, &X, &obs_times, &delta, &kerval,
                     &lq_x, &lq_w, &bsmat_tt_all, &kerval_tt_all, &ineqmat};
  
  int n_vars = p + gammap;
  std::vector<double> x(n_vars, 0.0);
  
  nlopt_opt opt;
  if (s == 0.0) {
    opt = nlopt_create(NLOPT_LD_SLSQP, n_vars);
  } else {
    opt = nlopt_create(NLOPT_LD_SLSQP, n_vars);
    if (ineqmat.n_rows > 0) {
      std::vector<double> tol(ineqmat.n_rows, 1e-8);
      nlopt_add_inequality_mconstraint(opt, ineqmat.n_rows, ineq_constraints, &data, tol.data());
    }
  }
  
  nlopt_set_min_objective(opt, nll_obj, &data);
  nlopt_set_xtol_rel(opt, xtol_rel);
  nlopt_set_maxeval(opt, maxeval);
  
  double minf;
  nlopt_result res = nlopt_optimize(opt, x.data(), &minf);
  
  nlopt_destroy(opt);
  return List::create(Named("status") = static_cast<int>(res),
                      Named("minimum") = minf,
                      Named("solution") = x);
}

// [[Rcpp::export]]
arma::mat calc_A(const arma::vec& beta, const arma::vec& gamma, double s, double h, const arma::mat& covariates, const arma::mat& bsmat, const arma::vec& X, const arma::vec& obs_times, const arma::vec& delta, const arma::vec& kerval, const arma::mat& bsmat_XX, int n_subj) {
  int n = X.n_elem;
  int p = covariates.n_cols;
  arma::mat A_est = arma::zeros<arma::mat>(p, p);

  arma::vec alpha_XX = bsmat_XX * gamma;
  arma::vec cov_beta = covariates * beta;

  for (int i = 0; i < n; ++i) {
    if (delta[i] == 1.0 && kerval[i] > 0) {
      double S0_sum = 0.0;
      arma::vec S1_sum = arma::zeros<arma::vec>(p);
      arma::mat S2_sum = arma::zeros<arma::mat>(p, p);

      for (int k = 0; k < n; ++k) {
        if (X[i] <= X[k]) {
          double dist_XX = X[i] - obs_times[k];
          if (dist_XX > 0) {
            double kerval_XX = std::max((1 - std::pow(dist_XX / h, 2)) * 0.75, 0.0) / h;
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
arma::mat calc_B(const arma::vec& beta, const arma::vec& gamma, double s, double h, const arma::mat& covariates, const arma::mat& bsmat, const arma::vec& X, const arma::vec& obs_times, const arma::vec& delta, const arma::vec& kerval, const arma::vec& id, const arma::mat& bsmat_XX, const arma::vec& lq_x, const arma::vec& lq_w, const arma::mat& bsmat_tt_all, const arma::mat& kerval_tt_all, int n_subj) {
  int n = X.n_elem;
  int p = covariates.n_cols;
  
  arma::vec alpha_XX = bsmat_XX * gamma;
  arma::vec cov_beta = covariates * beta;
  
  // use vectors instead of maps keyed by id
  std::vector<arma::vec> id_to_bb1(n_subj, arma::zeros<arma::vec>(p));
  std::vector<int> id_counts(n_subj, 0);
  for (int i = 0; i < n; ++i) {
    int idx = static_cast<int>(std::round(id[i])) - 1;
    if (idx < 0 || idx >= n_subj) continue;
    id_counts[idx] += 1;

    if (delta[i] == 1.0 && kerval[i] > 0 && (X[i] - obs_times[i]) > 0) {
      double S0_sum = 0.0;
      arma::vec S1_sum = arma::zeros<arma::vec>(p);
      for (int k = 0; k < n; ++k) {
        if (X[i] <= X[k]) {
          double dist_XX = X[i] - obs_times[k];
          if (dist_XX > 0) {
            double kerval_XX = std::max((1 - std::pow(dist_XX / h, 2)) * 0.75, 0.0) / h;
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
    double tt = 0.5 * lq_x[q] + 0.5;
    double weight = 0.5 * lq_w[q];
    double alpha_tt = dot(trans(bsmat_tt_all.row(q)), gamma);
    
    // vector for temporary r1 values per subject
    std::vector<double> id_to_r1(n_subj, 0.0);
    
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < n; ++i) {
        int idx = static_cast<int>(std::round(id[i])) - 1;
        if (idx < 0 || idx >= n_subj) continue;
        if (tt <= X[i]) {
          double k_tt = kerval_tt_all(q, i);
          if (k_tt > 0) {
            double S0_tt_sum = 0.0;
            double S1_tt_sum = 0.0;
            for (int k = 0; k < n; ++k) {
              if (tt <= X[k]) {
                double dist_k = tt - obs_times[k];
                if (dist_k > 0) {
                   double k_tt_k = std::max((1 - std::pow(dist_k / h, 2)) * 0.75, 0.0) / h;
                   if (k_tt_k > 0) {
                     double in_val = alpha_tt + dot(covariates.row(k), beta);
                     double td = trans_fun_d12o1(in_val, s) * k_tt_k;
                     S0_tt_sum += td;
                     S1_tt_sum += td * covariates(k, j);
                   }
                }
              }
            }
            if (S0_tt_sum > 0) {
               double S1_tt = S1_tt_sum / S0_tt_sum;
               double in_val = alpha_tt + dot(covariates.row(i), beta);
               id_to_r1[idx] += trans_fun_d(in_val, s) * k_tt * (S1_tt - covariates(i, j));
            }
          }
        }
      }
      for (int subj = 0; subj < n_subj; ++subj) {
         int count = id_counts[subj];
         if (count > 0) {
           id_to_bb2[subj](j) += weight * id_to_r1[subj] / static_cast<double>(count);
         }
      }
      std::fill(id_to_r1.begin(), id_to_r1.end(), 0.0);
    }
  }
  
  arma::mat B_est = arma::zeros<arma::mat>(p, p);
  for (int subj = 0; subj < n_subj; ++subj) {
    arma::vec bb_diff = id_to_bb1[subj] - id_to_bb2[subj];
    B_est += bb_diff * trans(bb_diff);
  }
  
  return B_est / n_subj;
}
