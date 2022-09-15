#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix set_constraint(NumericMatrix theta, 
                             NumericVector constraint) {
  int T = theta.ncol();
  int I = theta.nrow();
  NumericMatrix theta_new(I, T);
  for (int t = 0; t < T; t++) {
    if (theta(constraint[t], t) < 0) {
      theta_new(_, t) = -theta(_, t);
    } else {
      theta_new(_, t) = theta(_, t);
    }
  }
  return(theta_new);
}