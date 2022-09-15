#include <Rcpp.h>
using namespace Rcpp;

// get I * J matrix of omega (binary)
//[[Rcpp::export]]
NumericMatrix get_Eomega_bin(NumericVector theta, 
                             NumericVector alpha, 
                             NumericVector beta) {
  int I = theta.length(); int J = alpha.length();
  
  NumericMatrix omega(I, J);
  
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      double psi1 = alpha[j] + beta[j] * theta[i];
      
      omega(i, j) = std::tanh(psi1 / 2) / (2 * psi1);
      
      if (NumericVector::is_na(omega(i, j))) {
        omega(i, j) = 0.25;
      }
    }
  }
  return(omega);
}

// get I * J matrix of omega (binary, dynamic)
//[[Rcpp::export]]
NumericMatrix get_Eomega_bin_dyn(NumericMatrix theta, 
                                 NumericVector alpha, 
                                 NumericVector beta, 
                                 NumericVector bill_session) {
  int I = theta.nrow(); int J = alpha.length();
  
  NumericMatrix omega(I, J);
  
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      double psi1 = alpha[j] + beta[j] * theta(i, bill_session[j]);
      
      omega(i, j) = std::tanh(psi1 / 2) / (2 * psi1);
      
      if (NumericVector::is_na(omega(i, j))) {
        omega(i, j) = 0.25;
      }
    }
  }
  return(omega);
}

// get I * J matrix of omega (multinomial)
//[[Rcpp::export]]
List get_Eomega_mlt(NumericVector theta, 
                    NumericMatrix alpha, 
                    NumericMatrix beta) {
  int I = theta.length(); int J = alpha.nrow();
  
  NumericVector alpha1 = alpha(_, 0);
  NumericVector alpha2 = alpha(_, 1);
  NumericVector beta1 = beta(_, 0);
  NumericVector beta2 = beta(_, 1);
  
  NumericMatrix omega1(I, J); NumericMatrix omega2(I, J);
  
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      double psi1 = alpha1[j] + beta1[j] * theta[i];
      double psi2 = alpha2[j] + beta2[j] * theta[i];
      
      omega1(i, j) = std::tanh(psi1 / 2) / (2 * psi1);
      omega2(i, j) = std::tanh(psi2 / 2) / (2 * psi2);
      
      if (NumericVector::is_na(omega1(i, j))) {
        omega1(i, j) = 0.25;
      }
      if (NumericVector::is_na(omega2(i, j))) {
        omega2(i, j) = 0.25;
      }
    }
  }
  List L = List::create(Named("omega1") = omega1, Named("omega2") = omega2);
  return(L);
}

// get I * J matrix of omega (multinomial, dynamic)
//[[Rcpp::export]]
List get_Eomega_mlt_dyn(NumericMatrix theta, 
                        NumericMatrix alpha, 
                        NumericMatrix beta, 
                        NumericVector bill_session) {
  int I = theta.nrow(); int J = alpha.nrow();
  
  NumericVector alpha1 = alpha(_, 0);
  NumericVector alpha2 = alpha(_, 1);
  NumericVector beta1 = beta(_, 0);
  NumericVector beta2 = beta(_, 1);
  
  NumericMatrix omega1(I, J); NumericMatrix omega2(I, J);
  
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      double psi1 = alpha1[j] + beta1[j] * theta(i, bill_session[j]);
      double psi2 = alpha2[j] + beta2[j] * theta(i, bill_session[j]);
      
      omega1(i, j) = std::tanh(psi1 / 2) / (2 * psi1);
      omega2(i, j) = std::tanh(psi2 / 2) / (2 * psi2);
      
      if (NumericVector::is_na(omega1(i, j))) {
        omega1(i, j) = 0.25;
      } 
      if (NumericVector::is_na(omega2(i, j))) {
        omega2(i, j) = 0.25;
      }
    }
  }
  List L = List::create(Named("omega1") = omega1, Named("omega2") = omega2);
  return(L);
}