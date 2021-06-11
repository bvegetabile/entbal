#' @method plot entbal_binary
#' @export
plot.entbal_binary <- function(obj,
                               which_vars = 1,
                               plot_unweighted = TRUE){

  which_vars <- ifelse(sapply(which_vars, function(x) is.character(x)),
                       sapply(which_vars, function(x) which(x == colnames(obj$X))),
                       which_vars)

  X <- obj$X[ ,which_vars]
  NX <- length(which_vars)
  TA <- obj$TA
  for(j in 1:NX){
    if(plot_unweighted == T){
      par(mfrow = c(1,2))
      .plot_bal(X = as.matrix(X)[,j],
                TA = TA,
                pwts = rep(1, length(TA)),
                cov_name = colnames(obj$X)[which_vars[j]],
                main_title = 'Unweighted Empirical CDF')
    } else{
      par(mfrow = c(1,1))
    }
    .plot_bal(X = as.matrix(X)[,j],
              TA = TA,
              pwts = obj$wts,
              cov_name = colnames(obj$X)[which_vars[j]],
              main_title = 'Weighted Empirical CDF')
  }
}

.plot_bal <- function(X, TA, pwts, cov_name = NULL, main_title = ''){
  X_range <- range(X)
  var_pts <- seq(X_range[1], X_range[2], length.out = 500)
  uniq_TA <- unique(TA)
  plot(0, xlim = X_range, ylim = c(0,1),
       xlab = paste("Covariate:", cov_name),
       ylab = 'CDF',
       main = main_title,
       pch = 19, col = rgb(0,0,0,0))
  abline(h=c(0,1), lty=3)
  abline(v=X_range, lty=3)
  for(i in 1:length(uniq_TA)){
    ta <- TA == uniq_TA[i]
    WECDF <- wtd_ecdf(X[ta], pwts[ta])
    lines(var_pts, WECDF(var_pts), lty = i+1, lwd = 2)
  }
  legend('bottomright', legend = uniq_TA, lty=1 + 1:length(uniq_TA), lwd=2)
}
