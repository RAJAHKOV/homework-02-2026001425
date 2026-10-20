#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// All methods use sum of squared errors + lambda * squared coefficient norm.
// Return one LOOCV sum of squared prediction errors per positive lambda.
// [[Rcpp::export]]
arma::vec ridgeCVDirect(const arma::mat& X, const arma::vec& y,
                        const arma::vec& lambdas) {
  int n = X.n_rows;
  int p = X.n_cols;
  int M = lambdas.n_elem;
  
  arma::vec cv(M, arma::fill::zeros);
  
  for (int k = 0; k < M; ++k) {
    double lambda = lambdas(k);
    double cv_sum = 0.0;
    
    for (int i = 0; i < n; ++i) {
      // Delete observation i
      arma::uvec idx(n - 1);
      int pos = 0;
      for (int j = 0; j < n; ++j) {
        if (j != i) idx(pos++) = j;
      }
      
      arma::mat Xi = X.rows(idx);
      arma::vec yi = y.elem(idx);
      
      // Ridge fit on the remaining observations
      arma::mat A = Xi.t() * Xi
      + lambda * arma::eye(p, p);
      arma::vec b = Xi.t() * yi;
      arma::vec beta = arma::solve(A, b);
      
      double pred = arma::as_scalar(X.row(i) * beta);
      double err = y(i) - pred;
      cv_sum += err * err;
    }
    
    cv(k) = cv_sum;
  }
  
  return cv;
}

// [[Rcpp::export]]
arma::vec ridgeCVLeverage(const arma::mat& X, const arma::vec& y,
                          const arma::vec& lambdas) {
  int n=X.n_rows;
  int p=X.n_cols;
  int m=lambdas.n_elem;
  
  arma::vec cv(m,arma::fill::zeros);
  
  for (int i=0;i<m;i++)
  {
    double l=lambdas(i);
    arma::mat XtX=X.t()*X;
    arma::mat A=XtX+l*arma::eye(p,p);
    arma::vec beta=solve(A,X.t()*y);
    arma::vec yhat=X*beta;
    
    arma::mat AinvXt=solve(A,X.t());
    arma::mat S=X*AinvXt;
    arma::vec s_ii=S.diag();
    double total=0.0;
    for (int j=0;j<n;j++)
    {
      double r=y(j)-yhat(j);
      total+=pow(r/(1.0-s_ii(j)),2.0);
    }
    cv(i)=total;
  }
  return cv;
}

// [[Rcpp::export]]
arma::vec ridgeCVSVD(const arma::mat& X, const arma::vec& y,
                     const arma::vec& lambdas) {
  int n = X.n_rows;
  int p = X.n_cols;
  int M = lambdas.n_elem;
  
  arma::vec cv(M, arma::fill::zeros);
  
  // Thin SVD: X = U diag(s) V'
  arma::mat U, V;
  arma::vec s;
  
  bool ok = arma::svd_econ(U, s, V, X, "both", "std");
  if (!ok) {
    Rcpp::stop("SVD failed.");
  }
  
  arma::vec qty = U.t() * y;
  arma::vec s2 = arma::square(s);
  
  // Squared row entries of U, reused for every lambda
  arma::mat U2 = arma::square(U);
  
  for (int k = 0; k < M; ++k) {
    double lambda = lambdas(k);
    
    // Eigenvalue shrinkage factors
    arma::vec d = s2 / (s2 + lambda);
    
    // Fitted values: U diag(d) U' y
    arma::vec fitted = U * (d % qty);
    arma::vec residual = y - fitted;
    
    // Leverage: diag(U diag(d) U')
    arma::vec leverage = U2 * d;
    
    arma::vec loo_resid = residual / (1.0 - leverage);
    cv(k) = arma::dot(loo_resid, loo_resid);
  }
  
  return cv;
}