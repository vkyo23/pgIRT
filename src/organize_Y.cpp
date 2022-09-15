#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix organize_Y(NumericMatrix Y1, 
                         NumericMatrix Y2, 
                         NumericMatrix Y3) {
  int I = Y1.nrow();
  int J = Y1.ncol();
  NumericMatrix Y(I, J);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      if (Y1(i, j) == 1) {
        Y(i, j) = 1;
      } else if (Y2(i, j) == 1) {
        Y(i, j) = 2;
      } else if (Y3(i, j) == 1) {
        Y(i, j) = 3;
      }
    }
  }
  return(Y);
}
