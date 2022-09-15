#include <Rcpp.h>
using namespace Rcpp;

// update beta (binary)
//[[Rcpp::export]]
NumericVector update_beta_bin(NumericMatrix Y, 
                              NumericMatrix omega, 
                              NumericVector alpha, 
                              NumericVector theta,
                              double b0, 
                              double B0) {
  
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
        sig_part += omega_ij * std::pow(theta[i], 2.0);
        mu_part +=  omega_ij * theta[i];
        collapse += k_ij * theta[i];
      }
    }
    
    double sig = 1 / (sig_part + (1 / B0));
    double mu = sig * (collapse - mu_part * alpha[j] + (b0 / B0));
    
    sample(j) = mu;
  }
  
  return(sample);
}

// update beta (binary, dynamic)
//[[Rcpp::export]]
NumericVector update_beta_bin_dyn(NumericMatrix Y, 
                                  NumericMatrix omega, 
                                  NumericVector alpha, 
                                  NumericMatrix theta,
                                  double b0, 
                                  double B0, 
                                  NumericVector bill_session) {
  
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
        sig_part += omega_ij * std::pow(theta(i, bill_session[j]), 2.0);
        mu_part +=  omega_ij * theta(i, bill_session[j]);
        collapse += k_ij * theta(i, bill_session[j]);
      }
    }
    
    double sig =  1 / (sig_part + (1 / B0));
    double mu = sig * (collapse - mu_part * alpha[j] + (b0 / B0));
    
    sample(j) = mu;
  }
  
  return(sample);
}

// update beta (multinomial)
//[[Rcpp::export]]
NumericMatrix update_beta_mlt(NumericMatrix Y1, 
                              NumericMatrix Y2,
                              List Omega, 
                              NumericMatrix alpha, 
                              NumericVector theta,
                              NumericVector b0, 
                              NumericVector B0) {
  
  // Rows and Cols
  int I = Y1.nrow(); int J = Y1.ncol();
  
  NumericMatrix omega1 = Omega[0];
  NumericMatrix omega2 = Omega[1];
  NumericVector alpha1 = alpha(_, 0);
  NumericVector alpha2 = alpha(_, 1);
  
  double b0_1 = b0[0];
  double b0_2 = b0[1];
  double B0_1 = B0[0];
  double B0_2 = B0[1];
  
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
        sig_part1 += omega1_ij * std::pow(theta[i], 2.0);
        mu_part1 +=  omega1_ij * theta[i];
        collapse1 += k1_ij * theta[i];
      }
      if (NumericVector::is_na(alpha2[j]) | NumericVector::is_na(Y2(i, j))) {
        sig_part2 += 0;
        mu_part2 += 0;
        collapse2 += 0;
      } else {
        sig_part2 += omega2_ij * std::pow(theta[i], 2.0);
        mu_part2 += omega2_ij * theta[i];
        collapse2 += k2_ij * theta[i];
      }
    }
    
    double sig1 = 1 / (sig_part1 + (1 / std::sqrt(B0_1)));
    double mu1 = sig1 * (collapse1 - mu_part1 * alpha1[j] + (b0_1 / B0_1));
    double sig2 = 1 / (sig_part2 + (1 / std::sqrt(B0_2)));
    double mu2 = sig2 * (collapse2 - mu_part2 * alpha2[j] + (b0_2 / B0_2));
    
    sample(j, 0) = mu1; sample(j, 1) = mu2;
  }
  
  return(sample);
}

// update beta (multinomial, dynamic)
//[[Rcpp::export]]
NumericMatrix update_beta_mlt_dyn(NumericMatrix Y1, 
                                  NumericMatrix Y2,
                                  List Omega, 
                                  NumericMatrix alpha, 
                                  NumericMatrix theta,
                                  NumericVector b0, 
                                  NumericVector B0, 
                                  NumericVector bill_session) {
  
  // Rows and Cols
  int I = Y1.nrow(); int J = Y1.ncol();
  
  NumericMatrix omega1 = Omega[0];
  NumericMatrix omega2 = Omega[1];
  NumericVector alpha1 = alpha(_, 0);
  NumericVector alpha2 = alpha(_, 1);
  
  double b0_1 = b0[0];
  double b0_2 = b0[1];
  double B0_1 = B0[0];
  double B0_2 = B0[1];
  
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
        sig_part1 += omega1_ij * std::pow(theta(i, bill_session[j]), 2.0);
        mu_part1 +=  omega1_ij * theta(i, bill_session[j]);
        collapse1 += k1_ij * theta(i, bill_session[j]);
      }
      if (NumericVector::is_na(alpha2[j]) | NumericVector::is_na(Y2(i, j))) {
        sig_part2 += 0;
        mu_part2 += 0;
        collapse2 += 0;
      } else {
        sig_part2 += omega2_ij * std::pow(theta(i, bill_session[j]), 2.0);
        mu_part2 += omega2_ij * theta(i, bill_session[j]);
        collapse2 += k2_ij * theta(i, bill_session[j]);
      }
    }
    
    double sig1 = 1 / (sig_part1 + (1 / std::sqrt(B0_1)));
    double mu1 = sig1 * (collapse1 - mu_part1 * alpha1[j] + (b0_1 / std::sqrt(B0_1)));
    double sig2 = 1 / (sig_part2 + (1 / std::sqrt(B0_2)));
    double mu2 = sig2 * (collapse2 - mu_part2 * alpha2[j] + (b0_2 / std::sqrt(B0_2)));
    
    sample(j, 0) = mu1; sample(j, 1) = mu2;
  }
  
  return(sample);
}