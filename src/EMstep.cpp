//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::cube get_EOmega(const arma::mat& alpha, 
                      const arma::mat& beta, 
                      const arma::vec& theta) {
  int K = alpha.n_cols;
  int I = theta.n_rows;
  int J = beta.n_rows;
  arma::cube Omega(I, J, K);
  for (int k = 0; k < K; k++) {
    for (int j = 0; j < J; j++) {
      for (int i = 0; i < I; i++) {
        double psi = alpha(j, k) + beta(j, k) * theta[i];
        if (!NumericVector::is_na(psi)) {
          Omega(i, j, k) = std::tanh(psi/2) / (2 * psi);
        } else {
          Omega(i, j, k) = 0.25;
        }
      }
    }
  }
  return Omega;
}

arma::mat update_alpha(const arma::mat& Y,
                       const arma::cube& Omega,
                       const arma::mat& beta,
                       const arma::vec& theta,
                       const arma::mat& a0,
                       const arma::mat& A0) {
  int I = Y.n_rows;
  int J = Y.n_cols;
  int K = beta.n_cols;
  arma::mat draw(J, K);
  arma::cube IYs(I, J, K);
  for (int j = 0; j < J; j++) {
    for (int k = 0; k < K; k++) {
      double sig_part = 0;
      double mu_part = 0;
      double collapse = 0;
      int kk = k + 1.0;
      for (int i = 0; i < I; i++) {
        if (!NumericVector::is_na(Y(i, j))) {
          int IY = 0.0;
          if (Y(i, j) == kk) {
            IY = 1.0;
          }
          IYs(i, j, k) = IY;
          double s_ij = IY - 0.5;
          if (k > 0) {
            if ((IYs(i, j, k-1) == 1) | (IntegerVector::is_na(IYs(i, j, k-1)))) {
              IYs(i, j, k) = NA_INTEGER;
            }
          }
          if (!IntegerVector::is_na(IYs(i, j, k))) {
            sig_part += Omega(i, j, k);
            mu_part += s_ij;
            collapse += Omega(i, j, k) * theta(i, 0);
          } 
        } else {
          IYs(i, j, k) = NA_INTEGER;
        }
      }
      sig_part = sig_part + 1 / A0(j, k);
      mu_part = mu_part - collapse * beta(j, k) + a0(j, k) / A0(j, k);
      draw(j, k) = (1 / sig_part) * mu_part;
    }
  }
  return draw;
}

arma::mat update_beta(const arma::mat& Y,
                      const arma::cube& Omega,
                      const arma::mat& alpha,
                      const arma::vec& theta,
                      const arma::mat& b0,
                      const arma::mat& B0) {
  int I = Y.n_rows;
  int J = Y.n_cols;
  int K = alpha.n_cols;
  arma::mat draw(J, K);
  arma::cube IYs(I, J, K);
  for (int j = 0; j < J; j++) {
    for (int k = 0; k < K; k++) {
      double sig_part = 0;
      double mu_part = 0;
      double collapse = 0;
      int kk = k + 1.0;
      if (!NumericVector::is_na(alpha(j, k))) {
        for (int i = 0; i < I; i++) {
          if (!NumericVector::is_na(Y(i, j))) {
            int IY = 0;
            if (Y(i, j) == kk) {
              IY = 1;
            }
            IYs(i, j, k) = IY;
            double s_ij = IY - 0.5;
            if (k > 0) {
              if ((IYs(i, j, k-1) == 1) | (IntegerVector::is_na(IYs(i, j, k-1)))) {
                IYs(i, j, k) = NA_INTEGER;
              }
            }
            if (!IntegerVector::is_na(IYs(i, j, k))) {
              sig_part += Omega(i, j, k) * std::pow(theta[i], 2.0);
              mu_part += s_ij * theta[i];
              collapse += Omega(i, j, k) * theta(i, 0);
            } 
          } else {
            IYs(i, j, k) = NA_INTEGER;
          }
        }
      } else {
        sig_part = NA_REAL;
        mu_part = NA_REAL;
        collapse = NA_REAL;
      }
      sig_part = sig_part + 1 / B0(j, k);
      mu_part = mu_part - collapse * alpha(j, k) + b0(j, k) / B0(j, k);
      draw(j, k) = (1 / sig_part) * mu_part;
    }
  }
  return draw;
}

arma::vec update_theta(const arma::mat& Y, 
                       const arma::cube& Omega, 
                       const arma::mat& alpha, 
                       const arma::mat& beta,
                       const arma::vec& max_cat,
                       const int& constraint,
                       const bool& std) {
  int I = Y.n_rows;
  int J = Y.n_cols;
  arma::vec draw(I);
  arma::cube IYs(I, J, alpha.n_cols);
  for (int i = 0; i < I; i++) {
    double sig_part = 0;
    double mu_part = 0;
    for (int j = 0; j < J; j++) {
      for (int k = 0; k < (max_cat[j]-1); k++) {
        if (!NumericVector::is_na(Y(i, j))) {
          if (!NumericVector::is_na(alpha(j, k))) {
            int IY = 0;
            int kk = k + 1.0;
            if (Y(i, j) == kk) {
              IY = 1.0;
            }
            IYs(i, j, k) = IY;
            double s_ij = IY - 0.5; 
            if (k > 0) {
              if ((IYs(i, j, k-1) == 1) | (IntegerVector::is_na(IYs(i, j, k-1)))) {
                IYs(i, j, k) = NA_INTEGER;
              }
            }
            if (!IntegerVector::is_na(IYs(i, j, k))) {
              sig_part += Omega(i, j, k) * std::pow(beta(j, k), 2.0);
              mu_part += s_ij * beta(j, k) - beta(j, k) * Omega(i, j, k) * alpha(j, k);
            } 
          }
        } else {
          IYs(i, j, k) = NA_INTEGER;
        }
      } 
    }
    sig_part = sig_part + 1.0;
    draw(i) = (1 / sig_part) * mu_part;
  }
  if (draw[constraint] < 0) {
    draw = -draw;
  }
  if (std) {
    draw = (draw - arma::mean(draw)) / arma::stddev(draw);
  }
  return draw;
}


List Mstep(const arma::mat& Y,
           const arma::cube& Omega,
           const arma::mat& alpha,
           const arma::mat& beta,
           const arma::vec& theta,
           const arma::vec& max_cat,
           const int& constraint,
           const arma::mat& a0,
           const arma::mat& A0,
           const arma::mat& b0,
           const arma::mat& B0,
           const bool& std) {
  
  arma::vec theta_new = update_theta(Y,
                                     Omega,
                                     alpha, 
                                     beta, 
                                     max_cat,
                                     constraint,
                                     std);
  arma::mat beta_new = update_beta(Y,
                                   Omega,
                                   alpha,
                                   theta_new,
                                   b0, 
                                   B0);
  arma::mat alpha_new = update_alpha(Y,
                                     Omega,
                                     beta_new,
                                     theta_new,
                                     a0,
                                     A0);
  List L = List::create(Named("alpha") = alpha_new,
                        Named("beta") = beta_new,
                        Named("theta") = theta_new);
  return L;
}

//[[Rcpp::export]]
List EMstep(const arma::mat& Y,
            const arma::mat& alpha,
            const arma::mat& beta,
            const arma::vec& theta,
            const arma::vec& max_cat,
            const int& constraint,
            const arma::mat& a0,
            const arma::mat& A0,
            const arma::mat& b0,
            const arma::mat& B0,
            const int& maxit,
            const double& tol,
            const int& verbose,
            const bool& std) {
  
  arma::mat conv(maxit, 3);
  bool converge = true;
  int iter = 0;
  List par = List::create(alpha, 
                          beta, 
                          theta);
  for (int g = 0; g < (maxit + 1); g++) {
    arma::mat alpha_old = par[0];
    arma::mat beta_old = par[1];
    arma::vec theta_old = par[2];
    
    //Estep
    arma::cube Omega = get_EOmega(alpha_old, 
                                  beta_old,
                                  theta_old);
    //Mstep
    par = Mstep(Y,
                Omega,
                alpha_old,
                beta_old,
                theta_old,
                max_cat,
                constraint,
                a0,
                A0,
                b0,
                B0,
                std);
    
    if (g != 0) {
      arma::vec tmp_alpha1 = na_omit(as<NumericVector>(wrap(alpha_old)));
      arma::vec tmp_alpha2 = na_omit(as<NumericVector>(wrap(par[0])));
      arma::vec tmp_beta1 = na_omit(as<NumericVector>(wrap(beta_old)));
      arma::vec tmp_beta2 = na_omit(as<NumericVector>(wrap(par[1])));
      conv(g-1, 0) = cor(tmp_alpha1, tmp_alpha2).min();
      conv(g-1, 1) = cor(tmp_beta1, tmp_beta2).min();
      conv(g-1, 2) = cor(as<arma::vec>(wrap(theta_old)), as<arma::vec>(wrap(par[2]))).min();
      
      bool check = (1 - conv(g-1, 0) < tol) & (1 - conv(g-1, 1) < tol) & (1 - conv(g-1, 2) < tol);
      
      if (check) {
        arma::uvec seqq = as<arma::uvec>(wrap(seq(0, g-1)));
        conv = conv.rows(seqq);
        iter = g;
        break;
      } else if (g == maxit) {
        converge = false;
        iter = maxit;
        break;
      } else if (g % verbose == 0) {
        Rcout << "Iteration " << g << ": eval = " << max(1 - conv.row(g-1)) << '\n';
      }
    }
  }
  List L = List::create(Named("alpha") = par[0],
                        Named("beta") = par[1],
                                           Named("theta") = par[2],
                                                               Named("converge") = converge,
                                                               Named("cor") = conv,
                                                               Named("iter") = iter);
  return L;
}
