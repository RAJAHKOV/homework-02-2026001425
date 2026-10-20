#include <RcppArmadillo.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat kalmanSmoothing(arma::mat Yobs, arma::vec z0,
                          arma::mat A, arma::mat C,
                          double sigmaQ, double sigmaR) {
  // TODO: Store forward predictions and filtering results, then smooth backward.
  
  int T=Yobs.n_cols;
  //four-dimensional conditions
  int pz=z0.n_elem; 
  //two-dimensional observations
  int py=Yobs.n_rows;
  
  //noise covariance matrix
  arma::mat Q=sigmaQ*sigmaQ*arma::eye(pz,pz);
  arma::mat R=sigmaR*sigmaR*arma::eye(py,py);
  
  //restore prediction mean & variance and filtering mean&variance
  arma::mat z_pred(pz,T);
  arma::cube P_pred(pz,pz,T);
  arma::mat z_filt(pz,T);
  arma::cube P_filt(pz,pz,T);
  
  //Initial:z=0
  arma::vec zt_filt=z0;
  arma::mat Pt_filt=arma::zeros(pz,pz);
  
  for (int i=0;i<T;i++)
  {
    //Predict (i+1)|i
    arma::vec zt_pred=A*zt_filt;
    arma::mat Pt_pred=A*Pt_filt*A.t()+Q;
    
    z_pred.col(i)=zt_pred;
    P_pred.slice(i)=Pt_pred;
    
    //measurement
    arma::vec yt=Yobs.col(i);
    arma::mat S=C*Pt_pred*C.t()+R;
    arma::mat K=Pt_pred*C.t()*arma::inv(S);//By using the inverse of matrix
    
    zt_filt=zt_pred+K*(yt-C*zt_pred);
    Pt_filt=(arma::eye(pz,pz)-K*C)*Pt_pred;
    
    z_filt.col(i)=zt_filt;
    P_filt.slice(i)=Pt_filt;
  }
  
  arma::mat z_smooth=z_filt;
  arma::cube P_smooth=P_filt;
  
  for (int i=T-2;i>=0;i--)
  {
    arma::mat Jt=P_filt.slice(i)*A.t()*arma::inv(P_pred.slice(i+1));
    
    z_smooth.col(i)=z_filt.col(i)+Jt*(z_smooth.col(i+1)-z_pred.col(i+1));
    P_smooth.slice(i)=P_filt.slice(i)+Jt*(P_smooth.slice(i+1)-P_pred.slice(i+1))*Jt.t();
  }
  
  return z_smooth;
}
