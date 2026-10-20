#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

// Provided helper: complete, no student implementation required.
// Supports n >= p and a well-conditioned, invertible p-by-p B only.
// This is not a general-purpose GSVD algorithm for arbitrary or singular B.
// [[Rcpp::export]]
Rcpp::List computeGSVD(const arma::mat& X, const arma::mat& B) {
    if (X.n_cols == 0 || X.n_rows < X.n_cols ||
        B.n_rows != X.n_cols || B.n_cols != X.n_cols ||
        !X.is_finite() || !B.is_finite()) {
        Rcpp::stop("Require finite X (n-by-p, n >= p) and B (p-by-p).");
    }
    const double reciprocal_condition = arma::rcond(B);
    if (!std::isfinite(reciprocal_condition) || reciprocal_condition < 1e-10) {
        Rcpp::stop("This helper requires a nonsingular, well-conditioned B.");
    }

    // Form (X B^{-1})^T by solving a system, without an explicit inverse.
    arma::mat transformed_X_t;
    if (!arma::solve(transformed_X_t, B.t(), X.t(), arma::solve_opts::no_approx)) {
        Rcpp::stop("Failed to construct the transformed design.");
    }
    arma::mat U, V;
    arma::vec s;
    if (!arma::svd_econ(U, s, V, transformed_X_t.t())) {
        Rcpp::stop("SVD failed.");
    }

    arma::vec scale(s.n_elem);
    for (arma::uword i = 0; i < s.n_elem; ++i) {
        scale(i) = std::hypot(s(i), 1.0);
    }
    arma::vec alpha = s / scale;
    arma::vec beta = 1.0 / scale;
    arma::mat scaled_V = V;
    scaled_V.each_row() /= scale.t();
    arma::mat Z;
    if (!arma::solve(Z, B, scaled_V, arma::solve_opts::no_approx)) {
        Rcpp::stop("Failed to construct Z.");
    }

    // X Z = U diag(alpha), B Z = V diag(beta), alpha^2 + beta^2 = 1.
    // Return alpha and beta as R vectors, rather than p-by-1 matrices.
    return Rcpp::List::create(
        Rcpp::Named("U") = U,
        Rcpp::Named("U2") = V,
        Rcpp::Named("Z") = Z,
        Rcpp::Named("alpha") = Rcpp::NumericVector(alpha.begin(), alpha.end()),
        Rcpp::Named("beta") = Rcpp::NumericVector(beta.begin(), beta.end()));
}

// [[Rcpp::export]]
Rcpp::List tikhonovPath(const arma::mat& U, const arma::mat& Z,
                        const arma::vec& alpha, const arma::vec& beta,
                        const arma::vec& y, const arma::vec& lambdas) {
    // U contains the first p columns of U1 in the supplied GSVD.
    // beta contains GSVD diagonal entries, not regression coefficients.
    // TODO: Reuse U.t() * y and implement both GSVD sum formulas.
    int  n = U.n_rows;
  int  p = alpha.n_elem;
  int M = lambdas.n_elem;
  
  // Output matrices
  arma::mat coefficients(p, M);
  arma::mat fitted(n, M);
  
  // Reuse U.t() * y for all lambda values
  const arma::vec Uty = U.t() * y;
  
  // Precompute squared GSVD diagonal entries
  const arma::vec alpha2 = arma::square(alpha);
  const arma::vec beta2 = arma::square(beta);
  
  for (arma::uword k = 0; k < M; ++k) {
    const double lambda = lambdas(k);
    
    // Denominator: alpha_i^2 + lambda * beta_i^2
    const arma::vec denom = alpha2 + lambda * beta2;
    
    // Coefficient formula:
    // beta_hat = sum_i [alpha_i * (u_i' y) / denom_i] * z_i
    const arma::vec coef_weights = (alpha % Uty) / denom;
    coefficients.col(k) = Z * coef_weights;
    
    // Fitted-value formula:
    // y_hat = sum_i [alpha_i^2 * (u_i' y) / denom_i] * u_i
    const arma::vec fitted_weights = (alpha2 % Uty) / denom;
    fitted.col(k) = U * fitted_weights;
  }
  // Return a list with "coefficients" (p-by-M) and "fitted" (n-by-M).
  // Do not replace the formulas with a direct normal-equation solve.
    return Rcpp::List::create(
      Rcpp::Named("coefficients") = coefficients,
      Rcpp::Named("fitted") = fitted
    );
   
    
}
