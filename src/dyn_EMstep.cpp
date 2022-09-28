//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::cube dyn_get_EOmega(const arma::mat& alpha, 
                          const arma::mat& beta, 
                          const arma::mat& theta,
                          const arma::vec& bill_session) {
  int K = alpha.n_cols;
  int I = theta.n_rows;
  int J = beta.n_rows;
  arma::cube Omega(I, J, K);
  for (int k = 0; k < K; k++) {
    for (int j = 0; j < J; j++) {
      for (int i = 0; i < I; i++) {
        double psi = alpha(j, k) + beta(j, k) * theta(i, bill_session[j]);
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

arma::mat dyn_update_alpha(const arma::mat& Y,
                           const arma::cube& Omega,
                           const arma::mat& beta,
                           const arma::mat& theta,
                           const arma::mat& a0,
                           const arma::mat& A0,
                           const arma::vec& bill_session,
                           const arma::vec& matched_bill) {
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
      if (!NumericVector::is_na(beta(j, k))) {
        for (int i = 0; i < I; i++) {
          if (!NumericVector::is_na(Y(i, j))) {
            int IY = 0;
            int kk = k + 1;
            if (Y(i, j) == kk) {
              IY = 1;
            }
            IYs(i, j, k) = IY;
            if (k == 0) {
              sig_part += Omega(i, j, k);
              mu_part += IY - 0.5;
              collapse += Omega(i, j, k) * theta(i, bill_session[j]); 
            } else if ((IYs(i, j, k-1) == 1) | NumericVector::is_na(IYs(i, j, k-1))) {
              sig_part += 0;
              mu_part += 0;
              collapse += 0;
              IYs(i, j, k) = NA_REAL;
            } else {
              sig_part += Omega(i, j, k);
              mu_part += IY - 0.5;
              collapse += Omega(i, j, k) * theta(i, bill_session[j]);
            }
          } else {
            IYs(i, j, k) = NA_REAL;
          }
        }
      } else {
        sig_part = NA_REAL;
        mu_part = NA_REAL;
        collapse = NA_REAL;
      }
      draw(j, k) = (1/(1/A0(j, k) + sig_part)) * (a0(j, k)/A0(j, k) + mu_part - collapse * beta(j, k));
      if (!NumericVector::is_na(matched_bill[j])) {
        if (!NumericVector::is_na(draw(matched_bill[j], k)) & !NumericVector::is_na(beta(j, k))) {
          draw(j, k) = draw(matched_bill[j], k);
        } 
      } 
    }
  }
  return draw;
}

arma::mat dyn_update_beta(const arma::mat& Y,
                          const arma::cube& Omega,
                          const arma::mat& alpha,
                          const arma::mat& theta,
                          const arma::mat& b0,
                          const arma::mat& B0,
                          const arma::vec& bill_session) {
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
      if (!NumericVector::is_na(alpha(j, k))) {
        for (int i = 0; i < I; i++) {
          if (!NumericVector::is_na(Y(i, j))) {
            int IY = 0;
            int kk = k + 1;
            if (Y(i, j) == kk) {
              IY = 1;
            }
            IYs(i, j, k) = IY;
            if (k == 0) {
              sig_part += Omega(i, j, k) * std::pow(theta(i, bill_session[j]), 2.0);
              mu_part += (IY - 0.5) * theta(i, bill_session[j]);
              collapse += Omega(i, j, k) * theta(i, bill_session[j]); 
            } else if ((IYs(i, j, k-1) == 1) | NumericVector::is_na(IYs(i, j, k-1))) {
              sig_part += 0;
              mu_part += 0;
              collapse += 0;
              IYs(i, j, k) = NA_REAL;
            } else {
              sig_part += Omega(i, j, k) * std::pow(theta(i, bill_session[j]), 2.0);
              mu_part += (IY - 0.5) * theta(i, bill_session[j]);
              collapse += Omega(i, j, k) * theta(i, bill_session[j]); 
            }
          } else {
            IYs(i, j, k) = NA_REAL;
          }
        }
      } else {
        sig_part = NA_REAL;
        mu_part = NA_REAL;
        collapse = NA_REAL;
      }
      draw(j, k) = (1/(1/B0(j, k) + sig_part)) * (b0(j, k)/B0(j, k) + mu_part - collapse * alpha(j, k));
    }
  }
  return draw;
}

arma::mat dyn_update_theta(const arma::mat& Y, 
                           const arma::cube& Omega, 
                           const arma::mat& alpha, 
                           const arma::mat& beta,
                           const arma::vec& theta0,
                           const arma::vec& Delta0,
                           const arma::mat& Delta,
                           const arma::vec& max_cat,
                           const arma::vec& constraint,
                           const arma::mat& session_individual,
                           const arma::vec& bill_session,
                           const bool& std) {
  int I = Y.n_rows;
  arma::vec unq_t = unique(session_individual.col(0));
  arma::vec ind_id = unique(session_individual.col(1));
  int Time = unq_t.size();
  arma::mat draw(I, Time);
  arma::cube IYs(I, Y.n_cols, alpha.n_cols);
  for (int i = 0; i < I; i++) {
    arma::vec X_store(Time);
    arma::mat A(Time, Time);
    arma::vec B(Time);
    //find attending session
    arma::uvec which = find(session_individual.col(1) == ind_id[i]);
    arma::vec times_i = session_individual.elem(which);
    int T = times_i.n_rows;
    
    for (int t = 0; t < T; t++) {
      int time = times_i[t];
      arma::uvec bill_t_index = find(bill_session == time);
      //------
      arma::mat Y_t = Y.cols(bill_t_index);
      arma::mat alpha_t = alpha.rows(bill_t_index);
      arma::mat beta_t = beta.rows(bill_t_index);
      arma::vec max_cat_t = max_cat.rows(bill_t_index);
      //------
      int J = bill_t_index.size();
      //------
      double sig_part = 0;
      double mu_part = 0;
      for (int j = 0; j < J; j++) {
        //int jt = bill_t_index[j];
        for (int k = 0; k < (max_cat_t[j]-1); k++) {
          if (!NumericVector::is_na(Y_t(i, j))) {
            if (!NumericVector::is_na(alpha_t(j, k))) {
              int IY = 0;
              if (Y_t(i, j) == k + 1) {
                IY = 1;
              }
              IYs(i, bill_t_index[j], k) = IY;
              if (k == 0) {
                sig_part += Omega(i, bill_t_index[j], k) * std::pow(beta_t(j, k), 2.0);
                mu_part += (IY - 0.5) * beta_t(j, k) - Omega(i, bill_t_index[j], k) * alpha_t(j, k) * beta_t(j, k);
              } else if ((IYs(i, bill_t_index[j], k-1) == 1) | NumericVector::is_na(IYs(i, bill_t_index[j], k-1))) {
                IYs(i, bill_t_index[j], k) = NA_REAL;
              } else {
                sig_part += Omega(i, bill_t_index[j], k) * std::pow(beta_t(j, k), 2.0);
                mu_part += (IY - 0.5) * beta_t(j, k) - Omega(i, bill_t_index[j], k) * alpha_t(j, k) * beta_t(j, k);
              }
            } else {
              IYs(i, bill_t_index[j], k) = NA_REAL;
            }
          } else {
            IYs(i, bill_t_index[j], k) = NA_REAL;
          }
        }
      }
      if (t == 0) {
        double sig = sig_part + 1/Delta0[i] + 1/Delta(i, time);
        double mu = mu_part + theta0[i]/Delta0[i];
        
        A(t, t) = sig;
        A(t, t+1) = -1/Delta(i, time);
        A(t+1, t) = -1/Delta(i, time);
        B(t) = mu;
      } else if (t != (T - 1)) {
        double sig = sig_part + 2/Delta(i, time) + 1;
        double mu = mu_part;
        
        A(t, t) = sig;
        A(t, t+1) = -1/Delta(i, time);
        A(t+1, t) = -1/Delta(i, time);
        B(t) = mu;
      } else if (t == (T - 1)) {
        double sig = sig_part + 1/Delta(i, time) + 1;
        double mu = mu_part;
        A(t, t) = sig;
        B(t) = mu;
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
  
  NumericMatrix draw2 = as<NumericMatrix>(wrap(draw));
  
  for (int t = 0; t < Time; t++) {
    if (draw2(constraint[t], t) < 0) {
      draw2(_, t) = -draw2(_, t);
    }
    draw2(_, t) = ifelse(draw2(_, t) == 0, NA_REAL, draw2(_, t));
    if (std) {
      draw2(_, t) = (draw2(_, t) - mean(na_omit(draw2(_, t)))) / sd(na_omit(draw2(_, t)));
    }
  }
  draw = as<arma::mat>(wrap(draw2));
  
  return draw;
}

List dyn_Mstep(const arma::mat& Y,
               const arma::cube& Omega,
               const arma::mat& alpha,
               const arma::mat& beta,
               const arma::mat& theta,
               const arma::vec& max_cat,
               const arma::vec& constraint,
               const arma::mat& a0,
               const arma::mat& A0,
               const arma::mat& b0,
               const arma::mat& B0,
               const arma::vec& theta0,
               const arma::vec& Delta0,
               const arma::mat& Delta,
               const arma::mat& session_individual,
               const arma::vec& bill_session,
               const arma::vec& matched_bill,
               const bool& std) {
  
  arma::mat theta_new = dyn_update_theta(Y, 
                                         Omega, 
                                         alpha, 
                                         beta,
                                         theta0,
                                         Delta0,
                                         Delta,
                                         max_cat,
                                         constraint,
                                         session_individual,
                                         bill_session,
                                         std);
  arma::mat alpha_new = dyn_update_alpha(Y,
                                         Omega,
                                         beta,
                                         theta_new,
                                         a0,
                                         A0,
                                         bill_session,
                                         matched_bill);
  arma::mat beta_new = dyn_update_beta(Y,
                                       Omega,
                                       alpha_new,
                                       theta_new,
                                       b0, 
                                       B0,
                                       bill_session);
  List L = List::create(Named("alpha") = alpha_new,
                        Named("beta") = beta_new,
                        Named("theta") = theta_new);
  return L;
}

//[[Rcpp::export]]
List dyn_EMstep(const arma::mat& Y,
                const arma::mat& alpha,
                const arma::mat& beta,
                const arma::mat& theta,
                const arma::vec& max_cat,
                const arma::vec& constraint,
                const arma::mat& a0,
                const arma::mat& A0,
                const arma::mat& b0,
                const arma::mat& B0,
                const arma::vec& theta0,
                const arma::vec& Delta0,
                const arma::mat& Delta,
                const arma::mat& session_individual,
                const arma::vec& bill_session,
                const arma::vec& matched_bill, 
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
    arma::mat theta_old = par[2];
    
    //Estep
    arma::cube Omega = dyn_get_EOmega(alpha_old, 
                                      beta_old,
                                      theta_old,
                                      bill_session);
    //Mstep
    par = dyn_Mstep(Y,
                    Omega,
                    alpha,
                    beta,
                    theta,
                    max_cat,
                    constraint,
                    a0,
                    A0,
                    b0,
                    B0,
                    theta0,
                    Delta0,
                    Delta,
                    session_individual,
                    bill_session,
                    matched_bill,
                    std);
    
    if (g != 0) {
      arma::vec tmp_alpha1 = na_omit(as<NumericVector>(wrap(alpha_old)));
      arma::vec tmp_alpha2 = na_omit(as<NumericVector>(wrap(par[0])));
      arma::vec tmp_beta1 = na_omit(as<NumericVector>(wrap(beta_old)));
      arma::vec tmp_beta2 = na_omit(as<NumericVector>(wrap(par[1])));
      arma::vec tmp_theta1 = na_omit(as<NumericVector>(wrap(theta_old)));
      arma::vec tmp_theta2 = na_omit(as<NumericVector>(wrap(par[2])));
      conv(g-1, 0) = cor(tmp_alpha1, tmp_alpha2).min();
      conv(g-1, 1) = cor(tmp_beta1, tmp_beta2).min();
      conv(g-1, 2) = cor(tmp_theta1, tmp_theta2).min();
      
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

