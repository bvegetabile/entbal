#' Given Propensity Scores and Binary Treatment Indicators, Create ATE Weights.
#'
#' @param ps Propensity Scores from Binary Treatment
#' @param TI Treatment Indicators where 1 if "exposed" and zero otherwise
#'
#' @return ATE Weights
#'
#' @export
make_wts_ATE <- function(ps, TI){
  return (TI / ps + (1 - TI) / (1 - ps));
}

#' Given Propensity Scores and Binary Treatment Indicators, Create ATT Weights.
#'
#' @param ps Propensity Scores from Binary Treatment
#' @param TI Treatment Indicators where 1 if "exposed" and zero otherwise
#'
#' @return ATT Weights
#' @export
make_wts_ATT <- function(ps, TI){
  return (TI + (1 - TI) * ps / (1 - ps));
}

#' Given Propensity Scores and Binary Treatment Indicators, Create ATC Weights.
#'
#' @param ps Propensity Scores from Binary Treatment
#' @param TI Treatment Indicators where 1 if "exposed" and zero otherwise
#'
#' @return ATC Weights
#'
#' @export
make_wts_ATC <- function(ps, TI){
  return (TI * (1 - ps) / ps + (1 - TI));
}

#' Given Propensity Scores and Binary Treatment Indicators, Create ATO Weights.
#'
#' @param ps Propensity Scores from Binary Treatment
#' @param TI Treatment Indicators where 1 if "exposed" and zero otherwise
#'
#' @return ATO Weights
#'
#' @export
make_wts_ATO <- function(ps, TI){
  return (TI  * (1 - ps) + (1 - TI) * ps);
}

#' Given Propensity Scores, Group Indicators, and Weights: Estimate Group Expection
#'
#' @param X A Single Vector of Covariate Information
#' @param TI Indicators where 1 in group and zero otherwise. Typically treatment/exposure indicators.
#' @param w  Weights.
#'
#' Useful for covariate balance where the groups are defined by exposure levels.
#'
#' @return Weighted Group Expectation
#'
#' @export
wtd_mean <- function(X, TI, w) {
  sum(X * TI * w) / sum(TI  * w);
}

#' Given Propensity Scores, Group Indicators, and Weights: Estimate Group Non-Central Moment
#'
#' @param X A Single Vector of Covariate Information
#' @param TI Indicators where 1 in group and zero otherwise. Typically treatment/exposure indicators.
#' @param w  Weights.
#' @param mom Power of Moment, default to 1.
#'
#' Useful for covariate balance where the groups are defined by exposure levels.
#'
#' @return Weighted Group Non-Central Moment
#'
#' @export
wtd_moment <- function(X, TI, w, mom = 1) {
  return(sum((X ** mom) * TI * w) / sum(TI * w))
}

#' Given Propensity Scores, Group Indicators, and Weights: Estimate Weighted Group Sample Variance
#'
#' @param X A Single Vector of Covariate Information
#' @param TI Indicators where 1 in group and zero otherwise. Typically treatment/exposure indicators.
#' @param w  Weights.
#'
#' Useful for covariate balance where the groups are defined by exposure levels.
#'
#' @return Weighted Group Sample Variance
#'
#' @export
wtd_sd2 <- function(X, TI, w) {
  squares <- (X - wtd_mean(X, TI, w)) ** 2
  return(sum(squares * TI * w) / sum(TI * w))
}


.cov_mean_bal <- function(X, TI, w){
  XT <- wtd_mean(X, TI, w);
  XC <- wtd_mean(X, 1-TI, w);
  ST <- wtd_sd2(X, TI, w);
  SC <- wtd_sd2(X, 1-TI, w);
  return ((XT - XC) / sqrt((ST + SC)/2.0));
}

.cov_var_bal <- function(X, TI, w){
  ST <- wtd_sd2(X, TI, w);
  SC <- wtd_sd2(X, 1-TI, w);
  return (log(sqrt(ST)) - log(sqrt(SC)));
}

#' Simple Treatment Effect Estimation
#'
#' Simple difference in weighted expectations between the two groups defined by a vector of indicators.
#'
#' @param YO A Vector of Observed Outcomes
#' @param TI Indicators where 1 in group and zero otherwise. Typically treatment/exposure indicators.
#' @param wts  Weights.
#'
#' @return Weighted Group Sample Variance
#'
#' @export
est_binary_te <- function(YO, TI, wts){
  return (wtd_mean(YO, TI, wts) - wtd_mean(YO, 1 - TI, wts));
}

#' Covariate Balance Function - Standardized Difference in Means
#'
#'
#' @param X A Vector of Observed Covariates
#' @param TA Indicators where 1 in group and zero otherwise. Typically treatment/exposure indicators.
#' @param wts  Weights.
#'
#' @return Standardized Difference in Means
#'
#' @export
std_diff_means <- function(X, TA, wts){
  mean_diff <- wtd_mean(X, TA, wts) - wtd_mean(X, 1-TA, wts);
  mean_var <- (wtd_sd2(X, TA, wts) + wtd_sd2(X, 1 - TA, wts)) / 2;
  return (mean_diff / sqrt(mean_var));
}

#' Covariate Balance Function - Ratio of Weighted Standard Deviations
#'
#'
#' @param X A Vector of Observed Covariates
#' @param TA Indicators where 1 in group and zero otherwise. Typically treatment/exposure indicators.
#' @param wts  Weights.
#'
#' @return Ratio of Weighted Standard Deviations
#'
#' @export
log_sd_ratio <- function(X, TA, wts){
  sd1 <- sqrt(wtd_sd2(X, TA, wts));
  sd0 <- sqrt(wtd_sd2(X, 1 - TA, wts));
  return (log(sd1) - log(sd0));
}

#' Covariate Balance Function - Group Effective Sample Size
#'
#' @param w weights
#' @param TI Indicators where 1 in group and zero otherwise. Typically treatment/exposure indicators.
#'
#' @return Effective Sample of Group Defined by \code{TI==1}
#'
#' @export
group_ESS <- function(w, TI){
  return(sum((w * TI)) ** 2 / sum((w ** 2) *  TI));
}

#' Covariate Balance Function - Effective Sample Size
#'
#' @param w weights
#'
#' @return Effective Sample Size from Weights.
#'
#' @export
ESS <- function(w){
  return(sum(w)**2 / sum(w ** 2));
}

# weighted functions --------------

#' Covariate Balance Function - Weighted Mean
#'
#' @param wts weights
#' @param x vector of covariates
#'
#' @return weighted mean of x
#'
#' @export
wmean <- function(wts, x){sum(x * wts) / sum(wts)}

#' Covariate Balance Function - Weighted Variance
#'
#' @param wts weights
#' @param x vector of covariates
#'
#' @return weighted variance of x
#'
#' @export
wvar <- function(wts, x){
  wm <- wmean(wts, x)
  wv <- sum(wts * (x - wm)^2) / sum(wts)
  wv
}

#' Covariate Balance Function - Weighted Correlation
#'
#' @param wts weights
#' @param x first vector of covariates
#' @param y second vector of covariates
#'
#' @return weighted correlation between x and y
#'
#' @export
wcor <- function(wts, x, y){
  wmx <- wmean(wts, x)
  wmy <- wmean(wts, y)
  wvx <- wvar(wts, x)
  wvy <- wvar(wts, y)
  topval <- sum(wts * (x - wmx) * (y - wmy)) / sum(wts)
  topval / sqrt(wvy * wvx)
}

#' Covariate Balance Function - Weighted Covariance
#'
#' @param wts weights
#' @param x first vector of covariates
#' @param y second vector of covariates
#'
#' @return weighted covariance between x and y
#'
#' @export
wcov <- function(wts, x, y){
  wmx <- wmean(wts, x)
  wmy <- wmean(wts, y)
  wvx <- wvar(wts, x)
  wvy <- wvar(wts, y)
  topval <- sum(wts * (x - wmx) * (y - wmy)) / sum(wts)
  topval
}

.lm_ps <- function(Y, x, wts){
  X <- cbind(1,x)
  W <- c(wts)
  WX <- (W * X)
  invXtWX <- solve(t(X) %*% WX)
  hatmat <- invXtWX %*% t(WX)
  betas <- hatmat %*% Y
  pseudodf <- sum(wts)^2 / sum(wts^2)

  Yhat <- X %*% betas
  resids <- Y - Yhat

  varmat <- invXtWX %*% t(X) %*% ((W^2 * c(resids)^2) * X) %*% invXtWX
  std_errs <- sqrt(diag(varmat))

  low_int <- betas - 1.96 * std_errs
  upp_int <- betas + 1.96 * std_errs

  res <- cbind(betas, std_errs, betas/std_errs,
               2 * (1 - stats::pt(abs(betas/std_errs), df = pseudodf - 1)),
               low_int, upp_int)
  colnames(res) <- c('coef', 'stderrs', 't-value', 'p-value', 'low95', 'upp95')

  return(res)
}

.ksbal <- function(X, wts) {
  xr <- range(X)
  xs <- seq(xr[1], xr[2], length.out = 100)
  b4dist <- stats::ecdf(X)
  wdist <- entbal::wtd_ecdf(X, wts)
  max(abs(b4dist(xs) - wdist(xs)))
}
