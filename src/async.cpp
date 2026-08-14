#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Kernel-weighted estimating equations for asynchronous longitudinal data,
// following Cao, Zeng and Fine (2015, JRSS-B 77, 755-776).
//
// Both estimators in that paper -- time-invariant coefficients (their Section
// 2) and time-dependent coefficients (their Section 3) -- reduce to the same
// core.  Write the estimating function over contributing (response, covariate)
// pairs m = (j, k):
//
//   U(beta) = sum_m w_m X(S_k) [ Y(T_j) - mu{beta' X(S_k)} ].
//
// The link is evaluated at the covariate row alone, never at the response
// occasion, so every pair sharing a covariate row can be collapsed BEFORE beta
// is touched:
//
//   omega_k = sum_{m: k_m = k} w_m,      c_k = sum_{m: k_m = k} w_m Y(T_{j_m}),
//   U(beta) = sum_k [ c_k - omega_k mu_k ] X_k,
//   H(beta) = sum_k omega_k mu'_k X_k X_k'.
//
// So the pair design never enters the Newton loop: it is summarised once into
// two vectors of length n_x, and each iteration costs O(n_x p^2) rather than
// O(n_pair p^2).  Since k determines the subject, the subject-level scores that
// form the sandwich meat collapse the same way.
//
// The two estimators differ only in how (omega, c) are built:
//
//   time-invariant   w_m = W_h(T_j - S_k)                 -- needs the pairs
//   time-dependent   w_m = W_h1(T_j - t) W_h2(S_k - t)    -- factorises, so
//                    omega_k = vx_k sum_{j in i} wy_j and
//                    c_k     = vx_k sum_{j in i} wy_j Y_j, no pairs at all.

using namespace Rcpp;
using namespace arma;

namespace {

// link codes: 0 identity, 1 log, 2 logistic.
inline double mu_of(double eta, int link) {
  if (link == 1) return std::exp(eta);
  if (link == 2) return 1.0 / (1.0 + std::exp(-eta));
  return eta;
}

inline double dmu_of(double eta, int link) {
  if (link == 1) return std::exp(eta);
  if (link == 2) {
    const double m = 1.0 / (1.0 + std::exp(-eta));
    return m * (1.0 - m);
  }
  return 1.0;
}

} // namespace

//' Enumerate in-window (response, covariate) pairs within subject
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List async_pairs_cpp(const IntegerVector &y_sub, const NumericVector &y_time,
                     const NumericVector &x_time, const IntegerVector &lo,
                     const IntegerVector &hi, double h, double supp_lo,
                     double supp_hi) {
  const R_xlen_t ny = y_time.size();
  std::vector<int> pj, pk;
  std::vector<double> pz;

  for (R_xlen_t j = 0; j < ny; ++j) {
    if (j % 1000 == 0) Rcpp::checkUserInterrupt();
    const int i = y_sub[j];
    if (i == NA_INTEGER || i < 0 || i >= lo.size()) continue;
    // x_time is sorted by subject then by time, so subject i owns [lo, hi).
    for (int k = lo[i]; k < hi[i]; ++k) {
      const double z = (y_time[j] - x_time[k]) / h;
      if (z < supp_lo || z > supp_hi) continue;
      pj.push_back(static_cast<int>(j));
      pk.push_back(k);
      pz.push_back(z);
    }
  }

  return List::create(_["j"] = wrap(pj), _["k"] = wrap(pk), _["z"] = wrap(pz));
}

//' Collapse weighted pairs onto covariate rows (time-invariant coefficients)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List async_pair_accum_cpp(const IntegerVector &pj, const IntegerVector &pk,
                          const arma::vec &w, const arma::vec &yv, int n_x) {
  vec omega(n_x, fill::zeros), cc(n_x, fill::zeros);
  const R_xlen_t np = pj.size();

  for (R_xlen_t m = 0; m < np; ++m) {
    if (m % 100000 == 0) Rcpp::checkUserInterrupt();
    const int k = pk[m];
    // Weights are signed for order-p ladders, so no `> 0` guard here: dropping
    // negative pairs would silently change the estimating equation.
    omega[k] += w[m];
    cc[k] += w[m] * yv[pj[m]];
  }

  return List::create(_["omega"] = omega, _["c"] = cc);
}

//' Build the collapsed weights at one target time (time-dependent coefficients)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List async_td_accum_cpp(const IntegerVector &y_sub, const arma::vec &wy,
                        const arma::vec &yv, const IntegerVector &x_sub,
                        const arma::vec &vx, int n) {
  vec Sw(n, fill::zeros), Sy(n, fill::zeros);

  const R_xlen_t ny = y_sub.size();
  for (R_xlen_t j = 0; j < ny; ++j) {
    if (wy[j] == 0.0) continue;
    const int i = y_sub[j];
    if (i == NA_INTEGER || i < 0 || i >= n) continue;
    Sw[i] += wy[j];
    Sy[i] += wy[j] * yv[j];
  }

  const R_xlen_t nx = x_sub.size();
  vec omega(nx, fill::zeros), cc(nx, fill::zeros);
  for (R_xlen_t k = 0; k < nx; ++k) {
    if (vx[k] == 0.0) continue;
    const int i = x_sub[k];
    if (i == NA_INTEGER || i < 0 || i >= n) continue;
    omega[k] = Sw[i] * vx[k];
    cc[k] = Sy[i] * vx[k];
  }

  return List::create(_["omega"] = omega, _["c"] = cc);
}

//' Solve the collapsed kernel-weighted estimating equation
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List async_solve_cpp(const arma::mat &Xs, const IntegerVector &x_sub,
                     const arma::vec &omega, const arma::vec &cc, int n,
                     int link, int maxit, double tol, const arma::vec &start) {
  const uword p = Xs.n_cols;
  const uword nx = Xs.n_rows;
  vec beta = start;
  int conv = 0;
  int iter = 0;

  if (link == 0) {
    // Identity link: one linear solve, no iteration.
    mat A(p, p, fill::zeros);
    vec b(p, fill::zeros);
    for (uword k = 0; k < nx; ++k) {
      if (omega[k] == 0.0 && cc[k] == 0.0) continue;
      const rowvec xr = Xs.row(k);
      A += omega[k] * (xr.t() * xr);
      b += cc[k] * xr.t();
    }
    if (!solve(beta, A, b, solve_opts::no_approx)) {
      beta.fill(NA_REAL);
      conv = 1;
    }
  } else {
    conv = 2; // hit maxit unless the step test fires
    for (iter = 1; iter <= maxit; ++iter) {
      Rcpp::checkUserInterrupt();
      vec U(p, fill::zeros);
      mat H(p, p, fill::zeros);
      for (uword k = 0; k < nx; ++k) {
        if (omega[k] == 0.0 && cc[k] == 0.0) continue;
        const rowvec xr = Xs.row(k);
        const double eta = dot(xr, beta);
        U += (cc[k] - omega[k] * mu_of(eta, link)) * xr.t();
        H += (omega[k] * dmu_of(eta, link)) * (xr.t() * xr);
      }
      vec step;
      if (!solve(step, H, U, solve_opts::no_approx) || !step.is_finite()) {
        beta.fill(NA_REAL);
        conv = 1;
        break;
      }
      beta += step;
      if (max(abs(step)) < tol) {
        conv = 0;
        break;
      }
    }
  }

  // One final pass at the solution for the bread and the subject-level scores.
  // Doing it here rather than reusing the last Newton iterate matters: that
  // iterate was formed at the previous beta.
  mat bread(p, p, fill::zeros);
  mat U(n, p, fill::zeros);
  if (beta.is_finite()) {
    for (uword k = 0; k < nx; ++k) {
      if (omega[k] == 0.0 && cc[k] == 0.0) continue;
      const rowvec xr = Xs.row(k);
      const double eta = dot(xr, beta);
      const int i = x_sub[k];
      if (i != NA_INTEGER && i >= 0 && i < n) {
        U.row(i) += (cc[k] - omega[k] * mu_of(eta, link)) * xr;
      }
      bread += (omega[k] * dmu_of(eta, link)) * (xr.t() * xr);
    }
  }

  return List::create(_["beta"] = beta, _["bread"] = bread, _["U"] = U,
                      _["convergence"] = conv, _["iter"] = iter);
}
