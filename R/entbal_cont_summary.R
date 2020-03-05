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

summary.entbal_cont <- function(obj,
                                show_unweighted = TRUE,
                                show_higherorder = TRUE,
                                show_parameters = FALSE,
                                n_digits = 3){

  CB_unwtd <- cbind(apply(obj$X, 2, mean),
                    apply(obj$X, 2, sd),
                    apply(obj$X, 2, function(x) cor(x, obj$TA)),
                    t(apply(obj$X, 2, function(x) summary(lm(x ~ scale(obj$TA)))$coef[2,c(1,3,4)])))
  colnames(CB_unwtd) <- c('Mean',
                          'SD',
                          'Cor.',
                          'Beta',
                          'tval',
                          'pval')

  CB_wtd <- cbind(apply(obj$X, 2, function(x) wmean(obj$wts, x)),
                  apply(obj$X, 2, function(x) sqrt(wvar(obj$wts, x))),
                  apply(obj$X, 2, function(x) wcor(obj$wts, x, A)),
                  t(apply(obj$X, 2, function(x) .lm_ps(x, obj$TA, obj$wts)[2,c(1,3,4)])),
                  apply(obj$X, 2, function(x) .ksbal(x, obj$wts)))
  colnames(CB_wtd) <- c('wtd-Mean',
                        'wtd-SD',
                        'wtd-Cor.',
                        'wtd-Beta',
                        'wtd-tval',
                        'wtd-pval',
                        'wtd-KSstat')
  baltabs <- list('before' = CB_unwtd,
                  'after' = CB_wtd)

  if(show_unweighted){
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    cat('Unweighted Summary Statistics:\n')
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    print(round(baltabs$before, digits = n_digits))
  }
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  cat('Weighted Summary Statistics:\n')
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  colnames(CB_wtd) <- gsub('wtd-', '', colnames(CB_wtd))

  print(round(CB_wtd, digits = n_digits))
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))

  orig_N <- nrow(obj$X)
  esssum <- 1/sum(obj$wts^2)

  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  cat('Original & Effective Sample Sizes:\n')
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  SS <- cbind(orig_N, esssum, esssum/orig_N)
  colnames(SS) <- c('Orig N', 'ESS', 'Ratio')
  rownames(SS) <- c('')
  print(round(SS, digits = n_digits))
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))

  if(show_parameters){
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    cat('Parameter Values from Optimization:\n')
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    pars <- obj$opt_obj$par
    check_vals <- pars %in% obj$constraints
    pars <- data.frame('Value' = round(pars, digits = n_digits),
                       'Saturated' = check_vals)
    print(pars)
    if(sum(pars[,2]) > 0 ) {
      cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
      cat('Warning: some parameters have saturated at the constraint\n')

    }
  }
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))


  invisible(baltabs)
}
