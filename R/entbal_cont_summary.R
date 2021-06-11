
#' @method summary entbal_cont
#' @export
summary.entbal_cont <- function(object,
                                show_unweighted = TRUE,
                                show_higherorder = TRUE,
                                show_parameters = FALSE,
                                n_digits = 3){

  if(is.null(ncol(object$X))) {
    object$X <- matrix(object$X, ncol =1)
    colnames(object$X)[1] <- c('cov')
  }

  CB_unwtd <- cbind(apply(object$X, 2, mean),
                    apply(object$X, 2, sd),
                    apply(object$X, 2, function(x) cor(x, object$TA)),
                    t(apply(object$X, 2, function(x) summary(lm(x ~ scale(object$TA)))$coef[2,c(1,3,4)])))
  colnames(CB_unwtd) <- c('Mean',
                          'SD',
                          'Cor.',
                          'Beta',
                          'tval',
                          'pval')

  CB_wtd <- cbind(apply(object$X, 2, function(x) wmean(object$wts, x)),
                  apply(object$X, 2, function(x) sqrt(wvar(object$wts, x))),
                  apply(object$X, 2, function(x) wcor(object$wts, x, object$TA)),
                  t(apply(object$X, 2, function(x) .lm_ps(x, object$TA, object$wts)[2,c(1,3,4)])),
                  apply(object$X, 2, function(x) .ksbal(x, object$wts)))
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

  orig_N <- nrow(object$X)
  esssum <- 1/sum(object$wts^2)

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
    pars <- object$opt_obj$par
    check_vals <- pars %in% object$constraints
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
