#include <Rcpp.h>
using namespace Rcpp;

// Get predicted Y values for bootstrap (binary)
//[[Rcpp::export]]
IntegerMatrix Calc_predY_bin(NumericVector theta, 
                             NumericVector alpha, 
                             NumericVector beta) {
  int I = theta.length(); int J = alpha.length();
  
  IntegerMatrix Y(I, J);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      double exp1 = exp(alpha[j] + beta[j] * theta[i]);
      double denom = 1 + exp1;
      double p1 = exp1 / denom;
      Y(i, j) = R::rbinom(1, p1);
    } 
  }
  return(Y);
}

// Get predicted Y values for bootstrap (binary, dynamic)
//[[Rcpp::export]]
IntegerMatrix Calc_predY_bin_dyn(NumericMatrix theta, 
                                 NumericVector alpha, 
                                 NumericVector beta,
                                 NumericVector bill_session) {
  int I = theta.nrow(); int J = alpha.length();
  
  IntegerMatrix Y(I, J);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      double exp1 = exp(alpha[j] + beta[j] * theta(i, bill_session[j]));
      double denom = 1 + exp1;
      double p1 = exp1 / denom;
      Y(i, j) = R::rbinom(1, p1);
    } 
  }
  return(Y);
}

// Get predicted Y values for bootstrap (multinomial)
//[[Rcpp::export]]
List Calc_predY_mlt(NumericVector theta, 
                    NumericMatrix alpha, 
                    NumericMatrix beta) {
  int I = theta.length(); int J = alpha.nrow();
  
  IntegerMatrix Y1(I, J);
  IntegerMatrix Y2(I, J);
  IntegerMatrix Y3(I, J);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      double exp1 = exp(alpha(j, 0) + beta(j, 0) * theta[i]);
      double exp2 = exp(alpha(j, 1) + beta(j, 1) * theta[i]);
      if (NumericVector::is_na(alpha(j, 1))) {
        exp2= 0;
      }
      double denom = 1 + exp1 + exp2;
      double p1 = exp1 / denom;
      double p2 = exp2/ denom;
      double p3 = 1 / denom;
      NumericVector p = {p1, p2, p3};
      IntegerVector y(3);
      rmultinom(1, p.begin(), 3, y.begin());
      Y1(i, j) = y[0];
      Y2(i, j) = y[1];
      Y3(i, j) = y[2];
    } 
  }
  List L = List::create(Named("Y1") = Y1,
                        Named("Y2") = Y2,
                        Named("Y3") = Y3);
  return(L);
}

// Get predicted Y values for bootstrap (multinomial, dynamic)
//[[Rcpp::export]]
List Calc_predY_mlt_dyn(NumericMatrix theta, 
                        NumericMatrix alpha, 
                        NumericMatrix beta,
                        NumericVector bill_session) {
  int I = theta.nrow(); int J = alpha.nrow();
  
  IntegerMatrix Y1(I, J);
  IntegerMatrix Y2(I, J);
  IntegerMatrix Y3(I, J);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      double exp1 = exp(alpha(j, 0) + beta(j, 0) * theta(i, bill_session[j]));
      double exp2 = exp(alpha(j, 1) + beta(j, 1) * theta(i, bill_session[j]));
      if (NumericVector::is_na(alpha(j, 1))) {
        exp2= 0;
      }
      double denom = 1 + exp1 + exp2;
      double p1 = exp1 / denom;
      double p2 = exp2/ denom;
      double p3 = 1 / denom;
      NumericVector p = {p1, p2, p3};
      IntegerVector y(3);
      rmultinom(1, p.begin(), 3, y.begin());
      Y1(i, j) = y[0];
      Y2(i, j) = y[1];
      Y3(i, j) = y[2];
    } 
  }
  List L = List::create(Named("Y1") = Y1,
                        Named("Y2") = Y2,
                        Named("Y3") = Y3);
  return(L);
}