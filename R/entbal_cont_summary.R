
summary.entbal_cont <- function(obj,
                                show_unweighted = TRUE,
                                show_higherorder = TRUE,
                                show_parameters = FALSE,
                                n_digits = 3){

  if(is.null(ncol(obj$X))) {
    obj$X <- matrix(obj$X, ncol =1)
    colnames(obj$X)[1] <- c('cov')
  }

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
                  apply(obj$X, 2, function(x) wcor(obj$wts, x, obj$TA)),
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

  cat('Original & Effective Sample Sizes:\n')
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  SS <- cbind(orig_N, esssum, esssum/orig_N)
  colnames(SS) <- c('Orig N', 'ESS', 'Ratio')
  rownames(SS) <- c('')
  print(round(SS, digits = n_digits))

  if(show_parameters){
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    cat('Parameter Values from Optimization:\n')
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    pars <- obj$opt_obj$par
    check_vals <- pars %in% obj$constraints
    pars <- data.frame('Value' = round(pars, digits = n_digits),
                       'Saturated' = check_vals)
    print(pars[order(abs(pars$Value), decreasing = T),])
    if(sum(pars[,2]) > 0 ) {
      cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
      cat('Warning: some parameters have saturated at the constraint\n')
    }
  }
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))


  invisible(baltabs)
}
