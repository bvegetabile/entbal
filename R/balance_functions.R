
make_wts_ATE <- function(p, TI){
  return (TI / p + (1 - TI) / (1 - p));
}

make_wts_ATT <- function(p, TI){
  return (TI + (1 - TI) * p / (1 - p));
}

make_wts_ATC <- function(p, TI){
  return (TI * (1 - p) / p + (1 - TI));
}

make_wts_ATO <- function(p, TI){
  return (TI  * (1 - p) + (1 - TI) * p);
}

wtd_mean <- function(X, TI, w) {
  sum(X * TI * w) / sum(TI  * w);
}

wtd_moment <- function(X, TI, w, mom = 1) {
  return(sum((X ** mom) * TI * w) / sum(TI * w))
}

wtd_sd2 <- function(X, TI, w) {
  squares <- (X - wtd_mean(X, TI, w)) ** 2
  return(sum(squares * TI * w) / sum(TI * w))
}


.cov_mean_bal <- function(X, TI, w){
  XT <- wtd_mean(X, TI, w);
  XC <- wtd_mean(X, 1-TI, w);
  ST <- wtd_sd2(X, TI, w);
  SC <- wtd_sd2(X, 1-TI, w);
  return (XT - XC) / sqrt((ST + SC)/2.0);
}

.cov_var_bal <- function(X, TI, w){
  ST <- wtd_sd2(X, TI, w);
  SC <- wtd_sd2(X, 1-TI, w);
  return (log(sqrt(ST)) - log(sqrt(SC)));
}

est_binary_te <- function(YO, TI, wts){
  return (wtd_mean(YO, TI, wts) - wtd_mean(YO, 1 - TI, wts));
}

std_diff_means <- function(X, TA, wts){
  mean_diff <- wtd_mean(X, TA, wts) - wtd_mean(X, 1-TA, wts);
  mean_var <- (wtd_sd2(X, TA, wts) + wtd_sd2(X, 1 - TA, wts)) / 2;
  return (mean_diff / sqrt(mean_var));
}


log_sd_ratio <- function(X, TA, wts){
  sd1 <- sqrt(wtd_sd2(X, TA, wts));
  sd0 <- sqrt(wtd_sd2(X, 1 - TA, wts));
  return (log(sd1) - log(sd0));
}


group_ESS <- function(w, TI){
  return(sum((w * TI) ** 2) / sum((w ** 2) *  TI));
}

ESS <- function(w){
  return(sum(w)**2 / sum(w ** 2));
}

# weighted functions --------------
wmean <- function(wts, x){sum(x * wts) / sum(wts)}

wvar <- function(wts, x){
  wm <- wmean(wts, x)
  wv <- sum(wts * (x - wm)^2) / sum(wts)
  wv
}

wcor <- function(wts, x, y){
  wmx <- wmean(wts, x)
  wmy <- wmean(wts, y)
  wvx <- wvar(wts, x)
  wvy <- wvar(wts, y)
  topval <- sum(wts * (x - wmx) * (y - wmy)) / sum(wts)
  topval / sqrt(wvy * wvx)
}

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
  W <- diag(as.vector(wts))
  invXtWX <- solve(t(X) %*% W %*% X)
  hatmat <- invXtWX %*% t(X) %*% W
  betas <- hatmat %*% Y
  pseudodf <- sum(wts)^2 / sum(wts^2)

  Yhat <- X %*% betas
  resids <- Y - Yhat

  # varmat <- sighat * invXtWX %*% t(X) %*% W %*% W %*% X %*% invXtWX
  varmat <- invXtWX %*% t(X) %*% W %*% diag(as.vector(resids)^2) %*% W %*% X %*% invXtWX
  std_errs <- sqrt(diag(varmat))

  low_int <- betas - 1.96 * std_errs
  upp_int <- betas + 1.96 * std_errs

  res <- cbind(betas, std_errs, betas/std_errs,
               2 * (1 - pt(abs(betas/std_errs), df = pseudodf - 1)),
               low_int, upp_int)
  colnames(res) <- c('coef', 'stderrs', 't-value', 'p-value', 'low95', 'upp95')

  return(res)
}

.ksbal <- function(X, wts) {
  xr <- range(X)
  xs <- seq(xr[1], xr[2], length.out = 100)
  b4dist <- ecdf(X)
  wdist <- entbal::wtd_ecdf(X, wts)
  max(abs(b4dist(xs) - wdist(xs)))
}
