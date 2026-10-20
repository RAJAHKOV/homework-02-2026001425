#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Only adjacent indices outside the image are needed for a 3-by-3 kernel.
int mirrorIndex(int index, int length) {
  if (index < 0) return -index - 1;
  if (index >= length) return 2 * length - index - 1;
  return index;
}

// [[Rcpp::export]]
arma::mat imageConv(const arma::mat& img, const arma::mat& kernel) {
  // TODO: Apply the kernel, reflect boundaries, and clamp pixels to [0, 1].
  int H = img.n_rows;
  int W = img.n_cols;
  arma::mat out(H, W, arma::fill::zeros);
  
  // Traverse every element of (i,j)
  for(int i = 0; i < H; i++){
    for(int j = 0; j < W; j++){
      double sum_val = 0.0;
      // a,b: -1,0,1
      for(int a = -1; a <= 1; a++){
        for(int b = -1; b <= 1; b++){
          // the coordinates of the image are i+a, j+b .
          int im_row = mirrorIndex(i + a, H);
          int im_col = mirrorIndex(j + b, W);
          // the coordinates of arrays in cpp starts from 0,so kernel corresponds to 0,1,2(a+1,b+1)
          double A = img(im_row, im_col);
          double K = kernel(a+1, b+1);
          sum_val += A * K;
        }
      }
      // truncate the mat sum_val into [0, 1]
      if(sum_val < 0) sum_val = 0.0;
      if(sum_val > 1) sum_val = 1.0;
      out(i,j) = sum_val;
    }
  }
  return arma::mat(out);
}
