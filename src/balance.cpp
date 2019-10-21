// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <string>
#include <algorithm>
#include <cctype>

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat construct_X(const arma::mat & X,
                      const int & n_moments=2){
  if(n_moments > 1){
    int n_cols = X.n_cols;
    int n_obs = X.n_rows;
    int start_col;
    int end_col;
    arma::mat out_mat(n_obs, n_cols*n_moments);
    out_mat(arma::span::all, arma::span(0, n_cols-1)) = X;
    for(int m = 2; m <= n_moments; m++){
      start_col = (m-1) * n_cols;
      end_col = start_col + (n_cols-1);
      out_mat(arma::span::all, arma::span(start_col, end_col)) = arma::pow(X, m);
    }
    return out_mat;
  }
  return X;
}

// -----------------------------------------------------------------------------
//
// Weight functions for various Causal Inference Estimands
//
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::vec make_wts_ATE(const arma::vec & p,
                       const arma::vec & TI){
  arma::vec wts1 = TI / p;
  arma::vec wts0 = (1 - TI) / (1 - p);
  return (wts0 + wts1);
}

// [[Rcpp::export]]
arma::vec make_wts_ATT(const arma::vec & p,
                       const arma::vec & TI){
  arma::vec wts1 = TI;
  arma::vec wts0 = (1 - TI) % p / (1 - p);
  return (wts0 + wts1);
}

// [[Rcpp::export]]
arma::vec make_wts_ATC(const arma::vec & p,
                       const arma::vec & TI){
  arma::vec wts1 = TI % (1 - p) / p;
  arma::vec wts0 = (1 - TI);
  return (wts0 + wts1);
}

// [[Rcpp::export]]
arma::vec make_wts_ATO(const arma::vec & p,
                       const arma::vec & TI){
  arma::vec wts1 = TI % (1 - p);
  arma::vec wts0 = (1 - TI) % p;
  return (wts0 + wts1);
}

// -----------------------------------------------------------------------------
//
// Useful Weighted Covariate Balance Estimators
//
// -----------------------------------------------------------------------------

//' Calculate a weighted mean conditioned upon a vector of group indicators
//'
//' @param X Vector of observed covariate values
//' @param TI Vector representing group membership, i.e., \code{TI = (0, 1, 0, ...)}, where 1 indicates group membership
//' @param w Vector of weights, these need not sum to one in each group and are normalized in the equation
//' @return Estimate of the weighted mean in the group represented by \code{TI}
//' @export
// [[Rcpp::export]]
double wtd_mean(const arma::vec & X,
                const arma::vec & TI,
                const arma::vec & w) {
  return arma::dot(X, TI % w) / arma::sum(TI % w);
}

//' Calculate an estimate of the weighted raw moment \eqn{E(X^{p})} conditioned upon a vector of group indicators
//'
//' @param X Vector of observed covariate values
//' @param TI Vector representing group membership, i.e., \code{TI = (0, 1, 0, ...)}, where 1 indicates group membership
//' @param w Vector of weights, these need not sum to one in each group and are normalized in the equation
//' @param mom Power \eqn{p} of the required moment estimate
//' @return Estimate of the weighted raw moment in the group represented by \code{TI}
//' @export
// [[Rcpp::export]]
double wtd_moment(const arma::vec & X,
                  const arma::vec & TI,
                  const arma::vec & w,
                  const int mom = 1) {
  return arma::dot(arma::pow(X, mom), TI % w) / arma::sum(TI % w);
}

//' Calculate a weighted variance conditioned upon a vector of group indicators
//'
//' @param X Vector of observed covariate values
//' @param TI Vector representing group membership, i.e., \code{TI = (0, 1, 0, ...)}, where 1 indicates group membership
//' @param w Vector of weights, these need not sum to one in each group and are normalized in the equation
//' @return Estimate of the weighted variance in the group represented by \code{TI}
//' @export
// [[Rcpp::export]]
double wtd_sd2(const arma::vec & X,
               const arma::vec & TI,
               const arma::vec & w) {
  arma::vec squares = arma::pow(X - wtd_mean(X, TI, w), 2);
  return arma::dot(squares, TI % w) / arma::sum(TI % w);
}

// [[Rcpp::export(".cov_mean_bal")]]
double cov_mean_bal(const arma::vec & X,
                    const arma::vec & TI,
                    const arma::vec & w){
  double XT = wtd_mean(X, TI, w);
  double XC = wtd_mean(X, 1-TI, w);
  double ST = wtd_sd2(X, TI, w);
  double SC = wtd_sd2(X, 1-TI, w);
  return (XT - XC) / std::sqrt((ST + SC)/2.0);
}

// [[Rcpp::export(".cov_var_bal")]]
double cov_var_bal(const arma::vec & X,
                   const arma::vec & TI,
                   const arma::vec & w){
  double ST = wtd_sd2(X, TI, w);
  double SC = wtd_sd2(X, 1-TI, w);
  return std::log(std::sqrt(ST)) - std::log(std::sqrt(SC));
}

// [[Rcpp::export]]
double est_binary_te(const arma::vec & YO,
                     const arma::vec & TI,
                     const arma::vec & wts){
  return wtd_mean(YO, TI, wts) - wtd_mean(YO, 1 - TI, wts);
}

// [[Rcpp::export]]
double std_diff_means(const arma::vec X,
                      const arma::vec TA,
                      const arma::vec wts){
  double mean_diff = wtd_mean(X, TA, wts) - wtd_mean(X, 1-TA, wts);
  double mean_var = (wtd_sd2(X, TA, wts) + wtd_sd2(X, 1 - TA, wts)) / 2;
  return mean_diff / std::sqrt(mean_var);
}

// [[Rcpp::export]]
double log_sd_ratio(const arma::vec X,
                    const arma::vec TA,
                    const arma::vec wts){
  double sd1 = std::sqrt(wtd_sd2(X, TA, wts));
  double sd0 = std::sqrt(wtd_sd2(X, 1 - TA, wts));
  return std::log(sd1) - std::log(sd0);
}



// -----------------------------------------------------------------------------
//
// Weighted Effective Sample Sizes
//
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
double group_ESS(const arma::vec & w,
                 const arma::vec & TI){
  return(std::pow(arma::sum(w % TI), 2) / arma::sum(arma::pow(w%TI,2)));
}

// [[Rcpp::export]]
arma::vec ESS(const arma::vec & w){
  return(std::pow(arma::sum(w), 2) / arma::sum(arma::pow(w,2)));
}

// -----------------------------------------------------------------------------
//
// Unused Code ---------------------------
//
// -----------------------------------------------------------------------------

// // [[Rcpp::export]]
// arma::vec adaptive_deriv(const arma::mat & X,
//                          const arma::mat & TI,
//                          const arma::vec & w,
//                          const arma::vec & lambda){
//
//   arma::vec deriv_vec(X.n_rows, arma::fill::zeros);
//
//   deriv_vec += (TI.col(0) - TI.col(1)) %  (X * X.t() * ((TI.col(0) - TI.col(1)) % w));
//   deriv_vec += (lambda(0) * TI.col(0) + lambda(1) * TI.col(1)) % w;
//   return 2 * deriv_vec;
// }
//
// // [[Rcpp::export]]
// double adaptive_loss(const arma::mat & X,
//                      const arma::mat & TI,
//                      const arma::vec & w,
//                      const arma::vec & lambda){
//   int n_exp = TI.n_cols;
//   double covbal;
//   double penalty;
//
//   if(n_exp != 2){
//     std::cout << "ERROR : Adaptive Covariate Imbalance Minimization works for 2 Groups\n";
//     return 0;
//   }
//
//   covbal  = (w.t() * ((TI.col(0) - TI.col(1)) %  (X * X.t() * ((TI.col(0) - TI.col(1)) % w)))).eval()(0,0);
//   penalty = (w.t() * ((lambda(0) * TI.col(0) + lambda(1) * TI.col(1)) % w)).eval()(0,0);
//   return covbal + penalty;
// }
