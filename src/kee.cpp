#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec kee_cox_estequ(const arma::vec &beta, const arma::mat &covariates,
                         const arma::vec &X, const arma::vec &obs_times,
                         const arma::vec &delta, const arma::vec &kerval,
                         double h, int n_subj) {

  int n = X.n_elem;
  int p = covariates.n_cols;

  arma::vec part1 = arma::zeros<arma::vec>(p);
  arma::vec part2 = arma::zeros<arma::vec>(p);

  arma::vec expbeta = arma::exp(covariates * beta);

  for (int i = 0; i < n; ++i) {
    if (delta[i] == 1.0 && kerval[i] > 0) {
      part1 += kerval[i] * trans(covariates.row(i));

      double S0_XX_sum = 0.0;
      arma::vec S1_XX_sum = arma::zeros<arma::vec>(p);

      // Compute S0 and S1
      for (int k = 0; k < n; ++k) {
        if (X[i] <= X[k]) {
          double dist_XX = X[i] - obs_times[k];
          if (dist_XX > 0) {
            double kerval_XX =
                std::max((1.0 - std::pow(dist_XX / h, 2.0)) * 0.75, 0.0) / h;
            if (kerval_XX > 0) {
              double term = kerval_XX * expbeta[k];
              S0_XX_sum += term;
              S1_XX_sum += term * trans(covariates.row(k));
            }
          }
        }
      }

      if (S0_XX_sum > 0) {
        arma::vec Zbar_XX = S1_XX_sum / S0_XX_sum;
        part2 += kerval[i] * Zbar_XX;
      }
    }
  }

  return (part1 - part2) / static_cast<double>(n_subj);
}

// [[Rcpp::export]]
List kee_cox_var(const arma::vec &beta, const arma::mat &covariates,
                 const arma::vec &X, const arma::vec &obs_times,
                 const arma::vec &delta, const arma::vec &kerval, double h,
                 const arma::vec &id, int n_subj) {

  int n = X.n_elem;
  int p = covariates.n_cols;

  arma::mat W = arma::zeros<arma::mat>(p, p);
  // convert id values (1..n_subj) to zero-based indices; n_subj should
  // equal the number of unique subjects and ids must be coded accordingly.
  std::vector<arma::vec> id_to_sigma(n_subj, arma::zeros<arma::vec>(p));
  arma::vec expbeta = arma::exp(covariates * beta);

  for (int i = 0; i < n; ++i) {
    if (delta[i] == 1.0 && kerval[i] > 0) {
      double S0_XX_sum = 0.0;
      arma::vec S1_XX_sum = arma::zeros<arma::vec>(p);
      arma::mat S2_XX_sum = arma::zeros<arma::mat>(p, p);

      for (int k = 0; k < n; ++k) {
        if (X[i] <= X[k]) {
          double dist_XX = X[i] - obs_times[k];
          if (dist_XX > 0) {
            double kerval_XX =
                std::max((1.0 - std::pow(dist_XX / h, 2.0)) * 0.75, 0.0) / h;
            if (kerval_XX > 0) {
              double term = kerval_XX * expbeta[k];
              S0_XX_sum += term;
              S1_XX_sum += term * trans(covariates.row(k));
              S2_XX_sum +=
                  term * (trans(covariates.row(k)) * covariates.row(k));
            }
          }
        }
      }

      if (S0_XX_sum > 0) {
        arma::vec Zbar_XX = S1_XX_sum / S0_XX_sum;
        arma::mat S2_o_S0 = S2_XX_sum / S0_XX_sum;

        W += (S2_o_S0 - Zbar_XX * trans(Zbar_XX)) * kerval[i];
        int idx = static_cast<int>(std::round(id[i])) - 1;
        if (idx >= 0 && idx < n_subj) {
          id_to_sigma[idx] += kerval[i] * (trans(covariates.row(i)) - Zbar_XX);
        }
      }
    }
  }

  W = W / static_cast<double>(n_subj);

  arma::mat Sigma = arma::zeros<arma::mat>(p, p);
  for (int j = 0; j < n_subj; ++j) {
    Sigma += id_to_sigma[j] * trans(id_to_sigma[j]);
  }
  Sigma = Sigma / std::pow(static_cast<double>(n_subj), 2.0);

  return List::create(Named("W") = W, Named("Sigma") = Sigma);
}

// [[Rcpp::export]]
List kee_additive_est(const arma::mat &covariates, const arma::vec &X,
                      const arma::vec &obs_times, const arma::vec &delta,
                      const arma::vec &kerval, double h, const arma::vec &id,
                      const arma::vec &lq_x, const arma::vec &lq_w,
                      int n_subj) {

  int n = X.n_elem;
  int p = covariates.n_cols;
  int n_quad = lq_x.n_elem;

  arma::mat A = arma::zeros<arma::mat>(p, p);
  arma::vec B = arma::zeros<arma::vec>(p);

  // To compute A: needs double integral approximated by lq sum
  // use vector indexed by subject code (1..n_subj)
  std::vector<arma::mat> id_to_A(n_subj, arma::zeros<arma::mat>(p, p));

  for (int q = 0; q < n_quad; ++q) {
    double tt = 0.5 * lq_x[q] + 0.5;
    double weight = 0.5 * lq_w[q];

    double S0_tt_sum = 0.0;
    arma::vec S1_tt_sum = arma::zeros<arma::vec>(p);

    for (int k = 0; k < n; ++k) {
      if (tt <= X[k]) {
        double dist_tt = tt - obs_times[k];
        if (dist_tt > 0) {
          double kerval_tt =
              std::max((1.0 - std::pow(dist_tt / h, 2.0)) * 0.75, 0.0) / h;
          if (kerval_tt > 0) {
            S0_tt_sum += kerval_tt;
            S1_tt_sum += kerval_tt * trans(covariates.row(k));
          }
        }
      }
    }

    if (S0_tt_sum > 0) {
      arma::vec Zbar_tt = S1_tt_sum / S0_tt_sum;

      for (int i = 0; i < n; ++i) {
        if (tt <= X[i]) {
          double dist_tt = tt - obs_times[i];
          if (dist_tt > 0) {
            double kerval_tt =
                std::max((1.0 - std::pow(dist_tt / h, 2.0)) * 0.75, 0.0) / h;
            if (kerval_tt > 0) {
              arma::vec cov_m_Zbar = trans(covariates.row(i)) - Zbar_tt;
              int idx = static_cast<int>(std::round(id[i])) - 1;
              if (idx >= 0 && idx < n_subj) {
                // cov_m_Zbar is col vec (p x 1), covariates.row(i) is row vec
                // (1 x p) -> (p x p)
                id_to_A[idx] +=
                    weight * kerval_tt * (cov_m_Zbar * covariates.row(i));
              }
            }
          }
        }
      }
    }
  }

  for (int j = 0; j < n_subj; ++j) {
    A += id_to_A[j];
  }
  A = A / static_cast<double>(n_subj);

  // Compute B
  std::vector<arma::vec> id_to_B(n_subj, arma::zeros<arma::vec>(p));

  for (int i = 0; i < n; ++i) {
    if (kerval[i] > 0) {
      double delta_i = delta[i]; // B uses delta according to kee formulation
      if (delta_i == 1.0) {
        double S0_XX_sum = 0.0;
        arma::vec S1_XX_sum = arma::zeros<arma::vec>(p);

        for (int k = 0; k < n; ++k) {
          if (X[i] <= X[k]) {
            double dist_XX = X[i] - obs_times[k];
            if (dist_XX > 0) {
              double kerval_XX =
                  std::max((1.0 - std::pow(dist_XX / h, 2.0)) * 0.75, 0.0) / h;
              if (kerval_XX > 0) {
                S0_XX_sum += kerval_XX;
                S1_XX_sum += kerval_XX * trans(covariates.row(k));
              }
            }
          }
        }

        if (S0_XX_sum > 0) {
          arma::vec Zbar_XX = S1_XX_sum / S0_XX_sum;
          int idx = static_cast<int>(std::round(id[i])) - 1;
          if (idx >= 0 && idx < n_subj) {
            id_to_B[idx] += kerval[i] * (trans(covariates.row(i)) - Zbar_XX);
          }
        }
      }
    }
  }

  for (int j = 0; j < n_subj; ++j) {
    B += id_to_B[j];
  }
  B = B / static_cast<double>(n_subj);

  // Point estimate beta
  // beta = A^(-1) B
  arma::vec beta = arma::solve(A, B);

  // Compute Sigma
  arma::mat Sigma = arma::zeros<arma::mat>(p, p);
  for (int j = 0; j < n_subj; ++j) {
    arma::vec b_indi = id_to_B[j];
    arma::mat a_indi = id_to_A[j];
    arma::vec diff = b_indi - a_indi * beta;
    Sigma += diff * trans(diff);
  }
  Sigma = Sigma / std::pow(static_cast<double>(n_subj), 2.0);

  return List::create(Named("est") = beta, Named("A_est") = A,
                      Named("B_est") = B, Named("Sigma_est") = Sigma);
}
