#include <Rcpp.h>
using namespace Rcpp;

//update alpha (binary)
//[[Rcpp::export]]
NumericVector update_alpha_bin(NumericMatrix Y, NumericMatrix omega,
                               NumericVector beta, NumericVector theta,
                               double a0, double A0) {
  
  // Rows and Cols
  int I = Y.nrow(); int J = Y.ncol();
  
  // Define a store for sampling
  NumericVector sample(J);
  
  // Calculate posteriors
  for (int j = 0; j < J; j++) {
    // For summation
    double sig_part = 0; 
    double mu_part = 0;
    double collapse = 0; 
    
    for (int i = 0; i < I; i++) {
      double omega_ij = omega(i, j);
      double k_ij = Y(i, j) - 0.5;
      if (NumericVector::is_na(Y(i, j))) {
        sig_part += 0;
        mu_part += 0;
        collapse += 0;
      } else {
        sig_part += omega_ij;
        mu_part += omega_ij * theta[i];
        collapse += k_ij;
      }
    }
    
    double sig = 1 / (sig_part + (1 / A0));
    double mu = sig * (collapse - mu_part * beta[j] + (a0 / A0));
    
    sample(j) = mu;
  }
  
  return(sample);
}

//update alpha (binary, dynamic)
//[[Rcpp::export]]
NumericVector update_alpha_bin_dyn(NumericMatrix Y, NumericMatrix omega, 
                                   NumericVector beta, NumericMatrix theta,
                                   double a0, double A0,
                                   NumericVector bill_session,
                                   NumericVector matched_bill) {
  
  // Rows and Cols
  int I = Y.nrow(); int J = Y.ncol();
  
  // Define a store for sampling
  NumericVector sample(J);
  
  // Calculate posteriors
  for (int j = 0; j < J; j++) {
    // For summation
    double sig_part = 0; 
    double mu_part = 0;
    double collapse = 0; 
    
    for (int i = 0; i < I; i++) {
      double omega_ij = omega(i, j);
      double k_ij = Y(i, j) - 0.5;
      if (NumericVector::is_na(Y(i, j))) {
        sig_part += 0;
        mu_part += 0;
        collapse += 0;
      } else {
        sig_part += omega_ij;
        mu_part += omega_ij * theta(i, bill_session[j]);
        collapse += k_ij;
      }
    }
    
    double sig = 1 / (sig_part + (1 / A0));
    double mu = sig * (collapse - mu_part * beta[j] + (a0 / A0));
    sample(j) = mu;
    
    if (sum(matched_bill) != 0) {
      if (!NumericVector::is_na(matched_bill[j])) {
        sample(j) = sample(matched_bill[j]);
      }
    }
  }
  
  return(sample);
}

//update alpha (multinomial)
//[[Rcpp::export]]
NumericMatrix update_alpha_mlt(NumericMatrix Y1, NumericMatrix Y2, 
                               List Omega, 
                               NumericMatrix beta, NumericVector theta,
                               NumericVector a0, NumericVector A0) {
  
  // Rows and Cols
  int I = Y1.nrow(); int J = Y1.ncol();
  
  NumericMatrix omega1 = Omega[0];
  NumericMatrix omega2 = Omega[1];
  NumericVector beta1 = beta(_, 0);
  NumericVector beta2 = beta(_, 1);
  
  double a0_1 = a0[0];
  double a0_2 = a0[1];
  double A0_1 = A0[0];
  double A0_2 = A0[1];
  
  // Define a store for sampling
  NumericMatrix sample(J, 2);
  
  // Calculate posteriors
  for (int j = 0; j < J; j++) {
    // For summation
    double sig_part1 = 0; double mu_part1 = 0;
    double collapse1 = 0; double collapse2 = 0;
    double sig_part2 = 0; double mu_part2 = 0;
    
    for (int i = 0; i < I; i++) {
      double omega1_ij = omega1(i, j);
      double omega2_ij = omega2(i, j);
      double k1_ij = Y1(i, j) - 0.5;
      double k2_ij = Y2(i, j) - 0.5;
      if (NumericVector::is_na(Y1(i, j))) {
        sig_part1 += 0;
        mu_part1 += 0;
        collapse1 += 0;
      } else {
        sig_part1 += omega1_ij;
        mu_part1 += omega1_ij * theta[i];
        collapse1 += k1_ij;
      }
      if (NumericVector::is_na(Y2(i, j))) {
        sig_part2 += 0;
        mu_part2 += 0;
        collapse2 += 0;
      } else {
        sig_part2 += omega2_ij;
        mu_part2 += omega2_ij * theta[i];
        collapse2 += k2_ij;
      }
    }
    
    double sig1 = 1 / (sig_part1 + (1 / std::sqrt(A0_1)));
    double mu1 = sig1 * (collapse1 - mu_part1 * beta1[j] + (a0_1 / A0_1));
    
    double sig2 = 1 / (sig_part2 + (1 / std::sqrt(A0_2)));
    double mu2 = sig2 * (collapse2 - mu_part2 * beta2[j] + (a0_2 / A0_2));
    
    sample(j, 0) = mu1;
    sample(j, 1) = mu2;
  }
  
  return(sample);
}

//update alpha (multinomial, dynamic)
//[[Rcpp::export]]
NumericMatrix update_alpha_mlt_dyn(NumericMatrix Y1, NumericMatrix Y2, 
                                   List Omega, 
                                   NumericMatrix beta, NumericMatrix theta,
                                   NumericVector a0, NumericVector A0,
                                   NumericVector bill_session,
                                   NumericVector matched_bill) {
  
  // Rows and Cols
  int I = Y1.nrow(); int J = Y1.ncol();
  
  NumericMatrix omega1 = Omega[0];
  NumericMatrix omega2 = Omega[1];
  NumericVector beta1 = beta(_, 0);
  NumericVector beta2 = beta(_, 1);
  
  double a0_1 = a0[0];
  double a0_2 = a0[1];
  double A0_1 = A0[0];
  double A0_2 = A0[1];
  
  // Define a store for sampling
  NumericMatrix sample(J, 2);
  
  // Calculate posteriors
  for (int j = 0; j < J; j++) {
    // For summation
    double sig_part1 = 0; double mu_part1 = 0;
    double collapse1 = 0; double collapse2 = 0;
    double sig_part2 = 0; double mu_part2 = 0;
    
    for (int i = 0; i < I; i++) {
      double omega1_ij = omega1(i, j);
      double omega2_ij = omega2(i, j);
      double k1_ij = Y1(i, j) - 0.5;
      double k2_ij = Y2(i, j) - 0.5;
      if (NumericVector::is_na(Y1(i, j))) {
        sig_part1 += 0;
        mu_part1 += 0;
        collapse1 += 0;
      } else {
        sig_part1 += omega1_ij;
        mu_part1 += omega1_ij * theta(i, bill_session[j]);
        collapse1 += k1_ij;
      }
      if (NumericVector::is_na(Y2(i, j))) {
        sig_part2 += 0;
        mu_part2 += 0;
        collapse2 += 0;
      } else {
        sig_part2 += omega2_ij;
        mu_part2 += omega2_ij * theta(i, bill_session[j]);
        collapse2 += k2_ij;
      }
    }
    
    double sig1 = 1 / (sig_part1 + (1 / std::sqrt(A0_1)));
    double mu1 = sig1 * (collapse1 - mu_part1 * beta1[j] + (a0_1 / A0_1));
    
    double sig2 = 1 / (sig_part2 + (1 / std::sqrt(A0_2)));
    double mu2 = sig2 * (collapse2 - mu_part2 * beta2[j] + (a0_2 / A0_2));
    
    sample(j, 0) = mu1;
    sample(j, 1) = mu2;
    
    if (!NumericVector::is_na(matched_bill[j])) {
      sample(j, 0) = sample(matched_bill[j], 0);
      if (!NumericVector::is_na(sample(j, 1)) & !NumericVector::is_na(sample(matched_bill[j], 1))) {
        sample(j, 1) = sample(matched_bill[j], 1);
      }
    }
  }
  
  return(sample);
}