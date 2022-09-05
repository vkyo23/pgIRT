//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// update theta (binary)
//[[Rcpp::export]]
NumericVector update_theta_bin(NumericMatrix Y, NumericMatrix omega, NumericVector alpha, 
                               NumericVector beta, int constraint, bool is_const) {
  int I = Y.nrow();
  int J = Y.ncol();
  
  NumericVector theta(I);
  for (int i = 0; i < I; i++) {
    double left_part = 0;
    double right_part = 0;
    for (int j = 0; j < J; j++) {
      if (NumericVector::is_na(Y(i, j))) {
        left_part += 0;
        right_part += 0;
      } else {
        double omega_ij = omega(i, j);
        double k_ij = Y(i, j) - 0.5;
        left_part += std::pow(beta[j], 2.0) * omega_ij;
        right_part += beta[j] * k_ij - alpha[j] * beta[j] * omega_ij;
      }
    }
    left_part = left_part + 1;
    theta[i] = (1 / left_part) * right_part;
  }
  
  if (is_const) {
    if (theta[constraint - 1] < 0) {
      theta = -theta;
    }
  }
  
  return(theta);
}

// update theta (binary, dynamic)
//[[Rcpp::export]]
NumericMatrix update_theta_bin_dyn(arma::mat Y, arma::mat omega, arma::vec alpha, arma::vec beta,
                                   NumericVector theta0, NumericVector Delta0, double Delta, 
                                   IntegerVector constraint, bool is_const, 
                                   arma::mat session_individual, arma::vec bill_session) {
  
  
  // Rows
  int I = Y.n_rows; 
  
  arma::vec unq_t = unique(session_individual.col(0));
  arma::vec ind_id = unique(session_individual.col(1));
  int Time = unq_t.size();
  
  arma::mat draw(I, Time);
  
  for (int i = 0; i < I; i++) {
    arma::vec X_store(Time);
    arma::mat A(Time, Time);
    arma::vec B(Time);
    //find attending session
    arma::uvec which = find(session_individual.col(1) == ind_id[i]);
    arma::vec times_i = session_individual.elem(which);
    int T = times_i.n_rows;
    
    for (int t = 0; t < T; t++) {
      //get time index
      int time = times_i[t];
      arma::uvec bill_t_index = find(bill_session == time); 
      
      //subsetting
      arma::mat Y_t = Y.cols(bill_t_index);
      arma::mat omega_t = omega.cols(bill_t_index);
      
      arma::vec alpha_t = alpha(bill_t_index);
      arma::vec beta_t = beta(bill_t_index);
      
      
      int J = bill_t_index.size();
      
      double a_part = 0; 
      double b_part = 0; 
      
      for (int j = 0; j < J; j++) {
        if (NumericVector::is_na(Y_t(i, j))) {
          a_part += 0;
          b_part += 0;
        } else {
          a_part += std::pow(beta_t[j], 2.0) * omega_t(i, j);
          b_part += beta_t[j] * (Y_t(i, j) - 0.5) - beta_t[j] * omega_t(i, j) * alpha_t[j];
        }
      }
      
      if (t == 0) {
        // first attending
        double a = a_part + 1 / Delta0[i] + 1 / Delta;
        double b = b_part + theta0[i] / Delta0[i];
        
        A(t, t) = a;
        A(t, t+1) = -1/Delta;
        A(t+1, t) = -1/Delta;
        B(t) = b;
      } else if (t != (T - 1)) {
        double a = a_part + 2 / Delta + 1;
        double b = b_part;
        
        A(t, t) = a;
        A(t, t+1) = -1/Delta;
        A(t+1, t) = -1/Delta;
        B(t) = b;
        
      } else if (t == (T - 1)) {
        double a = a_part + 1 / Delta + 1;
        double b = b_part;
        
        A(t, t) = a;
        B(t) = b;
      }
    }
    arma::uvec diag_w = find(A.diag() != 0);
    arma::mat AA = A.submat(diag_w, diag_w);
    arma::vec BB = B.rows(diag_w);
    
    arma::vec X = inv(AA) * BB;
    arma::uvec t_i = as<arma::uvec>(wrap(times_i));
    X_store.rows(t_i) = X;
    
    draw.row(i) = X_store.t();
  }
  
  
  NumericMatrix draw_2 = as<NumericMatrix>(wrap(draw));
  
  
  for (int t = 0; t < Time; t++) {
    if (is_const) {
      if (draw_2(constraint[t] - 1, t) < 0) {
        draw_2(_, t) = -draw_2(_, t);
      }
    }
    draw_2(_, t) = ifelse(draw_2(_, t) == 0, NA_REAL, draw_2(_, t));
  }
  
  return(draw_2);
}

//update theta (multinomial)
//[[Rcpp::export]]
NumericVector update_theta_mlt(NumericMatrix Y1, NumericMatrix Y2, List Omega, 
                               NumericMatrix alpha, NumericMatrix beta, int constraint,
                               bool is_const,
                               NumericVector max_cat, NumericVector num_cat) {
  int I = Y1.nrow();
  int J = Y1.ncol();
  
  NumericVector alpha1 = alpha(_, 0);
  NumericVector alpha2 = alpha(_, 1);
  NumericVector beta1 = beta(_, 0);
  NumericVector beta2 = beta(_, 1);
  NumericMatrix omega1 = Omega[0];
  NumericMatrix omega2 = Omega[1];
  
  NumericVector theta(I);
  for (int i = 0; i < I; i++) {
    // For K_j = 3
    double a_part1 = 0; double a_part2 = 0;
    double b_part1 = 0; double b_part2 = 0;
    
    // For K_j = 2
    double a_partA = 0; double b_partA = 0;
    double a_partB = 0; double b_partB = 0;
    
    for (int j = 0; j < J; j++) {
      if ((max_cat[j] == 3) & (num_cat[j] == 3)) {
        if (NumericVector::is_na(Y1(i, j))) {
          a_part1 += 0;
          b_part1 += 0;
        } else {
          a_part1 += std::pow(beta1[j], 2.0) * omega1(i, j);
          b_part1 += beta1[j] * (Y1(i, j) - 0.5) - beta1[j] * omega1(i, j) * alpha1[j];
        }
        if (NumericVector::is_na(Y2(i, j))) {
          a_part2 += 0;
          b_part2 += 0;
        } else {
          a_part2 += std::pow(beta2[j], 2.0) * omega2(i, j);
          b_part2 += beta2[j] * (Y2(i, j) - 0.5) - beta2[j] * omega2(i, j) * alpha2[j];
        }
      } else if ((max_cat[j] == 3) & (num_cat[j] == 2)) {
        if (NumericVector::is_na(Y1(i, j))) {
          a_partA += 0;
          b_partA += 0;
        } else {
          a_partA += std::pow(beta1[j], 2.0) * omega1(i, j);
          b_partA += beta1[j] * (Y1(i, j) - 0.5) - beta1[j] * omega1(i, j) * alpha1[j];
        }
      } else if (max_cat[j] == 2) {
        if (NumericVector::is_na(Y1(i, j))) {
          a_partB += 0;
          b_partB += 0;
        } else {
          a_partB += std::pow(beta1[j], 2.0) * omega1(i, j);
          b_partB += beta1[j] * (Y1(i, j) - 0.5) - beta1[j] * omega1(i, j) * alpha1[j];
        }
      }
    }
    double a_part = a_part1 + a_part2 + a_partA + a_partB; 
    double b_part = b_part1 + b_part2 + b_partA + b_partB;
    
    theta(i) = (1 / (1 + a_part)) * b_part;
  }
  
  if (is_const) {
    if (theta[constraint - 1] < 0) {
      theta = -theta;
    }
  }
  
  return(theta);
}

// update theta (multinomial, dynamic)
//[[Rcpp::export]]
NumericMatrix update_theta_mlt_dyn(arma::mat Y1, arma::mat Y2, List Omega, NumericMatrix alpha, NumericMatrix beta,
                                   NumericVector theta0, NumericVector Delta0, double Delta, 
                                   IntegerVector constraint, bool is_const,
                                   arma::mat session_individual, arma::vec bill_session, 
                                   arma::vec max_cat, arma::vec num_cat) {
  
  
  // Rows
  int I = Y1.n_rows; 
  
  arma::mat omega1 = Omega[0];
  arma::mat omega2 = Omega[1];
  arma::vec alpha1 = alpha(_, 0);
  arma::vec alpha2 = alpha(_, 1); 
  arma::vec beta1 = beta(_, 0);
  arma::vec beta2 = beta(_, 1);
  
  
  arma::vec unq_t = unique(session_individual.col(0));
  arma::vec ind_id = unique(session_individual.col(1));
  int Time = unq_t.size();
  
  arma::mat draw(I, Time);
  
  for (int i = 0; i < I; i++) {
    arma::vec X_store(Time);
    arma::mat A(Time, Time);
    arma::vec B(Time);
    //find attending session
    arma::uvec which = find(session_individual.col(1) == ind_id[i]);
    arma::vec times_i = session_individual.elem(which);
    int T = times_i.n_rows;
    
    for (int t = 0; t < T; t++) {
      //get time index
      int time = times_i[t];
      arma::uvec bill_t_index = find(bill_session == time); 
      
      //subsetting
      arma::mat Y1_t = Y1.cols(bill_t_index);
      arma::mat Y2_t = Y2.cols(bill_t_index);
      arma::mat omega1_t = omega1.cols(bill_t_index);
      arma::mat omega2_t = omega2.cols(bill_t_index);
      arma::vec alpha1_t = alpha1.rows(bill_t_index);
      arma::vec alpha2_t = alpha2.rows(bill_t_index);
      arma::vec beta1_t = beta1.rows(bill_t_index);
      arma::vec beta2_t = beta2.rows(bill_t_index);
      arma::vec max_cat_t = max_cat.rows(bill_t_index);
      arma::vec num_cat_t = num_cat.rows(bill_t_index);
      
      int J = bill_t_index.size();
      
      // For K_j = 3
      double a_part1 = 0; double a_part2 = 0;
      double b_part1 = 0; double b_part2 = 0;
      
      // For K_j = 2
      double a_partA = 0; double b_partA = 0;
      double a_partB = 0; double b_partB = 0;
      for (int j = 0; j < J; j++) {
        if ((max_cat_t[j] == 3) & (num_cat_t[j] == 3)) {
          if (NumericVector::is_na(Y1_t(i, j))) {
            a_part1 += 0;
            b_part1 += 0;
          } else {
            a_part1 += std::pow(beta1_t[j], 2.0) * omega1_t(i, j);
            b_part1 += beta1_t[j] * (Y1_t(i, j) - 0.5) - beta1_t[j] * omega1_t(i, j) * alpha1_t[j];
          }
          if (NumericVector::is_na(Y2_t(i, j))) {
            a_part2 += 0;
            b_part2 += 0;
          } else {
            a_part2 += std::pow(beta2_t[j], 2.0) * omega2_t(i, j);
            b_part2 += beta2_t[j] * (Y2_t(i, j) - 0.5) - beta2_t[j] * omega2_t(i, j) * alpha2_t[j];
          }
        } else if ((max_cat_t[j] == 3) & (num_cat_t[j] == 2)) {
          if (NumericVector::is_na(Y1_t(i, j))) {
            a_partA += 0;
            b_partA += 0;
          } else {
            a_partA += std::pow(beta1_t[j], 2.0) * omega1_t(i, j);
            b_partA += beta1_t[j] * (Y1_t(i, j) - 0.5) - beta1_t[j] * omega1_t(i, j) * alpha1_t[j];
          }
        } else if (max_cat_t[j] == 2) {
          if (NumericVector::is_na(Y1_t(i, j))) {
            a_partB += 0;
            b_partB += 0;
          } else {
            a_partB += std::pow(beta1_t[j], 2.0) * omega1_t(i, j);
            b_partB += beta1_t[j] * (Y1_t(i, j) - 0.5) - beta1_t[j] * omega1_t(i, j) * alpha1_t[j];
          }
        }
      }
      double a_part = a_part1 + a_part2 + a_partA + a_partB; 
      double b_part = b_part1 + b_part2 + b_partA + b_partB;
      
      if (t == 0) {
        // first attending
        double a = a_part + 1 / Delta0[i] + 1 / Delta;
        double b = b_part + theta0[i] / Delta0[i];
        
        A(t, t) = a;
        A(t, t+1) = -1/Delta;
        A(t+1, t) = -1/Delta;
        B(t) = b;
      } else if (t != (T - 1)) {
        double a = a_part + 2 / Delta + 1;
        double b = b_part;
        
        A(t, t) = a;
        A(t, t+1) = -1/Delta;
        A(t+1, t) = -1/Delta;
        B(t) = b;
        
      } else if (t == (T - 1)) {
        double a = a_part + 1 / Delta + 1;
        double b = b_part;
        
        A(t, t) = a;
        B(t) = b;
      }
    }
    arma::uvec diag_w = find(A.diag() != 0);
    arma::mat AA = A.submat(diag_w, diag_w);
    arma::vec BB = B.rows(diag_w);
    
    arma::vec X = inv(AA) * BB;
    arma::uvec t_i = as<arma::uvec>(wrap(times_i));
    X_store.rows(t_i) = X;
    
    draw.row(i) = X_store.t();
  }
  
  
  NumericMatrix draw_2 = as<NumericMatrix>(wrap(draw));
  
  for (int t = 0; t < Time; t++) {
    if (is_const) {
      if (draw_2(constraint[t] - 1, t) < 0) {
        draw_2(_, t) = -draw_2(_, t);
      }
    }
    draw_2(_, t) = ifelse(draw_2(_, t) == 0, NA_REAL, draw_2(_, t));
  }
  
  return(draw_2);
}