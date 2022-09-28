//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::mat calc_predY(const arma::mat& alpha, 
                     const arma::mat& beta, 
                     const arma::vec& theta,
                     const arma::vec& max_cat) {
  int I = theta.n_rows;
  int J = alpha.n_rows;
  arma::mat predY(I, J);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      for (int k = 0; k < max_cat[j]-1; k++) {
        double psi = alpha(j, k) + beta(j, k) * theta[i];
        double p_k = std::exp(psi) / (1 + std::exp(psi));
        int y_ijk = R::rbinom(1, p_k);
        if (y_ijk == 1) {
          predY(i, j) = k + 1;
        }
      }
      if (predY(i, j) == 0) {
        predY(i, j) = max_cat[j];
      }
    }
  }
  return predY;
}

arma::mat calc_predY_dyn(const arma::mat& alpha, 
                         const arma::mat& beta, 
                         const arma::mat& theta,
                         const arma::vec& max_cat,
                         const arma::vec& bill_session) {
  int I = theta.n_rows;
  int J = alpha.n_rows;
  arma::mat predY(I, J);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      for (int k = 0; k < max_cat[j]-1; k++) {
        double psi = alpha(j, k) + beta(j, k) * theta(i, bill_session[j]);
        double p_k = std::exp(psi) / (1 + std::exp(psi));
        int y_ijk = R::rbinom(1, p_k);
        if (y_ijk == 1) {
          predY(i, j) = k + 1;
        }
      }
      if (predY(i, j) == 0) {
        predY(i, j) = max_cat[j];
      }
    }
  }
  return predY;
}

arma::cube calc_probY(const arma::mat& alpha, 
                      const arma::mat& beta, 
                      const arma::vec& theta,
                      const arma::vec& max_cat) {
  int I = theta.n_rows;
  int J = alpha.n_rows;
  int K = alpha.n_cols + 1;
  if (K == 2) {
    K = 1;
  }
  arma::cube probY(I, J, K);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      arma::vec p_i(max_cat[j]-1);
      for (int k = 0; k < max_cat[j]-1; k++) {
        if (!NumericVector::is_na(alpha(j, k))) {
          double psi = alpha(j, k) + beta(j, k) * theta[i];
          p_i(k) = std::exp(psi) / (1 + std::exp(psi));
          probY(i, j, k) = std::exp(psi) / (1 + std::exp(psi));
        }
        if ((K != 1) & (k == max_cat[j]-2)) {
          probY(i, j, k+1) = 1 - sum(p_i);
        }
      }
      if (max_cat[j]-1 == 0) {
        double psi = alpha(j, 0) + beta(j, 0) * theta[i];
        probY(i, j, 0) = std::exp(psi) / (1 + std::exp(psi));
      }
    }
  }
  return probY;
}

arma::cube calc_probY_dyn(const arma::mat& alpha, 
                          const arma::mat& beta, 
                          const arma::mat& theta,
                          const arma::vec& max_cat,
                          const arma::vec& bill_session) {
  int I = theta.n_rows;
  int J = alpha.n_rows;
  int K = alpha.n_cols + 1;
  if (K == 2) {
    K = 1;
  }
  arma::cube probY(I, J, K);
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < I; i++) {
      arma::vec p_i(max_cat[j]-1);
      for (int k = 0; k < max_cat[j]-1; k++) {
        if (!NumericVector::is_na(alpha(j, k))) {
          double psi = alpha(j, k) + beta(j, k) * theta(i, bill_session[j]);
          p_i(k) = std::exp(psi) / (1 + std::exp(psi));
          probY(i, j, k) = std::exp(psi) / (1 + std::exp(psi));
        }
        if ((K != 1) & (k == max_cat[j]-2)) {
          probY(i, j, k+1) = 1 - sum(p_i);
        }
      }
      if (max_cat[j]-1 == 0) {
        double psi = alpha(j, 0) + beta(j, 0) * theta(i, bill_session[j]);
        probY(i, j, 0) = std::exp(psi) / (1 + std::exp(psi));
      }
    }
  }
  return probY;
}

//[[Rcpp::export]]
List prediction(const arma::mat& alpha, 
                const arma::mat& beta, 
                const arma::mat& theta,
                const arma::vec& max_cat,
                const arma::vec& bill_session,
                const String& model,
                const String& type) {
  int I = theta.n_rows;
  int J = alpha.n_rows;
  int K = alpha.n_cols;
  if (type == "prob") {
    arma::cube out(I, J, K);
    if (model == "default") {
      out = calc_probY(alpha,
                       beta,
                       theta.col(0),
                       max_cat);
    } else {
      out = calc_probY_dyn(alpha,
                           beta,
                           theta,
                           max_cat,
                           bill_session);
    }
    List L = List::create(Named("out") = out);
    return L;
  } else {
    arma::mat out(I, J);
    if (model == "default") {
      out = calc_predY(alpha,
                       beta,
                       theta.col(0),
                       max_cat);
    } else {
      out = calc_predY_dyn(alpha,
                           beta,
                           theta,
                           max_cat,
                           bill_session);
    }
    List L = List::create(Named("out") = out);
    return L;
  }
}