#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double totCorC(arma::mat x){
  x = x.t();
  arma::colvec u = sum(x, 1);
  arma::mat d = x.each_col() - u/x.n_cols;
  arma::mat d2 = pow(d, 2.0);
  arma::colvec v = sum(d2, 1);
  arma::mat z = d.each_col() / pow(v, 0.5);
  arma::rowvec a = pow(sum(z, 0), 2.0);
  arma::colvec b = sum(d2.each_col() / v, 1);
  double(r) = sum(a) - sum(b);
  return(r);
}
