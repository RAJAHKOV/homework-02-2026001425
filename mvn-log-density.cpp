#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
double mvnLogDensity(const arma::vec& x, const arma::vec& mu,
                     const arma::mat& Sigma) {
  // TODO: Use a lower Cholesky factor and a triangular solve.
  int p=x.n_elem;
  
  mat L = arma::chol(Sigma,"lower");
  
  arma::vec r=x-mu;
  
  arma::vec v=arma::solve(L,r);
  
  double sumlogDiag = sum(log(L.diag()));
    
  double v_2=dot(v,v);
  
  const double c=log(2.0*M_PI);
  
  double logpdf=-0.5*p*c-sumlogDiag-0.5*v_2;
  
  return logpdf;
}