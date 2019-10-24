extract_wts <- function(obj){
  if(!(class(obj) %in% c("entbal_multiclass", "entbal_binary"))) stop('Works with objects of class : "entbal_multiclass", "entbal_binary"')
  obj$wts
}

print_baltables <- function(X, TA, wts, show_unweighted=TRUE, n_digits = 2){
  cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  balance_table <- entbal::baltable(X, TA, n_digits=n_digits)
  weightd_table <- entbal::baltable(X, TA, wts, show_unweighted = F, n_digits=n_digits)
  colnames(weightd_table) <- colnames(balance_table)
  if(show_unweighted){
    cat('Unweighted Balance Statistics:\n')
    cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
    print(balance_table)
  }
  cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  cat('Weighted Balance Statistics:\n')
  cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  print(weightd_table)
  cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))

  t_levels <- unique(TA)
  n_uniq <- length(t_levels)
  ESS <- rep(NA, n_uniq)
  for(i in 1:n_uniq){
    TAindicators <- TA == t_levels[i]
    ESSG <- group_ESS(wts, TAindicators)
    cat(paste('TA: ', t_levels[i], ',   Original N = ', sum(TAindicators), '\n',
              paste(rep(' ', 6 + max(nchar(t_levels))), collapse = ''),
              'Weighted ESS = ', round(ESSG, n_digits), '\n', sep =''))
    ESS[i] <- ESSG
  }
  cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  invisible(list('unweighted_balance_table' = balance_table,
                 'weighted_balance_table' = weightd_table,
                 'TA_levels' = t_levels,
                 'ESS' = ESS))
}


summary.entbal_binary <- function(obj, show_unweighted = TRUE, n_digits=2){
  # cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  # balance_table <- entbal::baltable(obj$X, obj$TA, n_digits=n_digits)
  # weightd_table <- entbal::baltable(obj$X, obj$TA, obj$wts, show_unweighted = F, n_digits=n_digits)
  # colnames(weightd_table) <- colnames(balance_table)
  # if(show_unweighted){
  #   cat('Unweighted Balance Statistics:\n')
  #   cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  #   print(balance_table)
  # }
  # cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  # cat('Weighted Balance Statistics:\n')
  # cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  # print(weightd_table)
  # cat(paste(paste(rep('-', 90), collapse = ''), '\n', sep=''))
  # invisible(list('unweighted_balance_table' = balance_table,
  #                'weighted_balance_table' = weightd_table))
  #
  outtab <- print_baltables(obj$X, obj$TA, obj$wts, show_unweighted=show_unweighted, n_digits=n_digits)

}



ks_test <- function(X, TA, wts=rep(1,length(X)), n_pts=1000){
  xmin <- min(X)
  xmax <- max(X)
  int_pts <- seq(xmin, xmax, length.out = n_pts)
  t_wts <- ifelse(TA == 1, wts, 0)
  c_wts <- ifelse(TA == 0, wts, 0)
  t_fn <- wtd_ecdf(X, wts = t_wts)
  c_fn <- wtd_ecdf(X, wts = c_wts)
  return(max(abs(t_fn(int_pts) - c_fn(int_pts))))
}

ks <- function(x,z,w) {
  w[z==1] <-  w[z==1]/sum(w[z==1])
  w[z==0] <- -w[z==0]/sum(w[z==0])
  ind <- order(x)
  cumv <- abs(cumsum(w[ind]))
  cumv <- cumv[diff(x[ind]) != 0]
  return(ifelse(length(cumv) > 0, max(cumv), 0))
}

.n_uniq <- function(x){length(unique(x))}

.find_continuous <- function(X){apply(X, 2, .n_uniq) > 2}

# make_Xmat <- function(X, m = 1){
#   cont_var <- .find_continuous(X)
#   X_con <- X[, cont_var]
#   X_bin <- X[, !cont_var]
#   Xout <- cbind(X_bin, X_con)
#   if(m > 1){
#     for(i in 2:m){
#       Xout <- cbind(Xout, X_con**i)
#     }
#   }
#   Xout
# }

make_Xmat <- function(X, m = 1){
  cont_var <- .find_continuous(X)
  X_con <- X[, cont_var]
  NC <- ncol(X_con)
  X_bin <- X[, !cont_var]
  Xout <- cbind(X_bin, X_con)
  if(m > 1){
    for(i in 1:NC){
      nu <- length(unique(X_con[,i]))
      if(nu > m){
        Xout <- cbind(Xout, poly(X_con[,i], m))
      } else {
        Xout <- cbind(Xout, poly(X_con[,i], nu-1))
      }

    }
  }
  Xout
}

group_ESS <- function(w, TA){
  ESSG <- sum(w * TA)^2 / sum(w^2 * TA)
  ESSG
}



#' Calculate covariate balance statistics
#'
#' @param Xmat A matrix of covariates
#' @param TA A vector of treatment assignments (i.e, 1 if treated, 0 in control)
#' @param wts A vector of estimated (or true) weights
#' @param show_unweighted Query if the the balance table should contain unweighted estimates
#' @param n_digits Number of digits to print for the table
#' @return A table of covariate balance statistics
#' @examples
#'
baltable <- function(Xmat, TA,
                     wts = NULL,
                     show_unweighted=T,
                     n_digits = 3){
  n_obs <- length(TA)
  n_cols <- ncol(Xmat)
  exposures <- unique(TA)
  n_exp <- length(exposures)
  cov_names <- colnames(Xmat)
  b4_wts <- rep(1, n_obs)

  if(is.null(wts)){
    bal_table <- matrix(NA, nrow = n_cols, ncol = 7)
    for(d in 1:n_cols){
      Xd <- Xmat[ ,d]
      bal_table[d, 1] <- wtd_mean(Xd, TA, b4_wts)
      bal_table[d, 2] <- sqrt(wtd_sd2(Xd, TA, b4_wts))
      bal_table[d, 3] <- wtd_mean(Xd, 1-TA, b4_wts)
      bal_table[d, 4] <- sqrt(wtd_sd2(Xd, 1-TA, b4_wts))
      bal_table[d, 5] <- .cov_mean_bal(Xd, TA, b4_wts)
      if(length(unique(Xd)) > 2){
        bal_table[d, 6] <- .cov_var_bal(Xd, TA, b4_wts)
      } else {
        bal_table[d, 6] <- bal_table[d, 1] * (1 - bal_table[d, 1])
      }
      bal_table[d, 7] <- ks_test(Xd, TA, b4_wts)
    }
    colnames(bal_table) <- c('MeanGroup1',
                             'SEGroup1',
                             'MeanGroup0',
                             'SEGroup0',
                             'StdDiffMeans',
                             'LogRatioSE',
                             'MaxKS')
  } else {
    bal_table <- matrix(NA, nrow = n_cols, ncol = 14)
    for(d in 1:n_cols){
      Xd <- Xmat[ ,d]
      bal_table[d, 1] <- wtd_mean(Xd, TA, b4_wts)
      bal_table[d, 2] <- sqrt(wtd_sd2(Xd, TA, b4_wts))
      bal_table[d, 3] <- wtd_mean(Xd, 1-TA, b4_wts)
      bal_table[d, 4] <- sqrt(wtd_sd2(Xd, 1-TA, b4_wts))
      bal_table[d, 5] <- .cov_mean_bal(Xd, TA, b4_wts)
      bal_table[d, 7] <- ks_test(Xd, TA, b4_wts)

      bal_table[d, 8] <- wtd_mean(Xd, TA, wts)
      bal_table[d, 9] <- sqrt(wtd_sd2(Xd, TA, wts))
      bal_table[d, 10] <- wtd_mean(Xd, 1-TA, wts)
      bal_table[d, 11] <- sqrt(wtd_sd2(Xd, 1-TA, wts))
      bal_table[d, 12] <- .cov_mean_bal(Xd, TA, wts)
      if(length(unique(Xd)) > 2){
        bal_table[d, 6] <- .cov_var_bal(Xd, TA, b4_wts)
        bal_table[d, 13] <- .cov_var_bal(Xd, TA, wts)
      } else {
        bal_table[d, 6] <- bal_table[d, 6] * (1 - bal_table[d, 6])
        bal_table[d, 13] <- bal_table[d, 12] * (1 - bal_table[d, 12])
      }
      bal_table[d, 14] <- ks_test(Xd, TA, wts)
    }
    colnames(bal_table) <- c('MeanGroup1',
                             'SEGroup1',
                             'MeanGroup0',
                             'SEGroup0',
                             'StdDiffMeans',
                             'LogRatioSE',
                             'MaxKS',
                             'wtd-MeanGroup1',
                             'wtd-SEGroup1',
                             'wtd-MeanGroup0',
                             'wtd-SEGroup0',
                             'wtd-StdDiffMeans',
                             'wtd-LogRatioSE',
                             'wtd-MaxKS')
  }
  rownames(bal_table) <- cov_names
  if(show_unweighted == T){
    return(round(bal_table, digits = n_digits))
  } else {
    return(round(bal_table[,8:14], digits = n_digits))
  }
  return(NULL)
}
