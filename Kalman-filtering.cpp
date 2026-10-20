#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
mat kalmanFiltering(mat Yobs, vec z0, mat A, mat C,
                    double sigmaQ, double sigmaR) {
  const uword p = z0.n_elem;
  const uword tTotal = Yobs.n_cols;
  vec mu = z0;
  mat Sigma = zeros<mat>(p, p);
  mat Zmat = zeros<mat>(p, tTotal);
  
  for (uword t = 0; t < tTotal; ++t) {
    // Prediction.
    vec mu_f = A * mu;
    mat Sigma_f = A * Sigma * A.t();
    Sigma_f.diag() += sigmaQ * sigmaQ;
    
    // Measurement update.
    mat S = C * Sigma_f * C.t();
    S.diag() += sigmaR * sigmaR;
    mat Kgain = solve(S, C * Sigma_f);
    Kgain = Kgain.t();
    mu = mu_f + Kgain * (Yobs.col(t) - C * mu_f);
    Sigma = Sigma_f - Kgain * C * Sigma_f;
    Zmat.col(t) = mu;
  }
  return Zmat;
}