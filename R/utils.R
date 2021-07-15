#' Returns Object Weights
#'
#' @param object object of class \code{entbal_multiclass} or \code{entbal_binary}
#'
#' @export
extract_wts <- function(object){
  if(!(class(object) %in% c("entbal_multiclass", "entbal_binary"))) stop('Works with objects of class : "entbal_multiclass", "entbal_binary"')
  object$wts
}

#' Binary Balance Table
#'
#' @param X Matrix of covariates
#' @param TA Binary indicator
#' @param wts Weights
#' @param show_unweighted Show unweighted balanced statistics? default \code{TRUE}
#' @param n_digits How many digits to round table, default 2
#'
#' @export
print_baltables <- function(X, TA, wts, show_unweighted=TRUE, n_digits = 2){
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  balance_table <- entbal::baltable(as.matrix(X), TA, n_digits=n_digits)
  weightd_table <- entbal::baltable(as.matrix(X), TA, wts, show_unweighted = F, n_digits=n_digits)

  if(is.null(ncol(weightd_table))){
    weightd_table <- matrix(weightd_table, nrow = 1)
  }
  colnames(weightd_table) <- colnames(balance_table)
  if(show_unweighted){
    cat('Unweighted Balance Statistics:\n')
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    print(balance_table)
  }
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  cat('Weighted Balance Statistics:\n')
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  print(weightd_table)
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))

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
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  invisible(list('unweighted_balance_table' = balance_table,
                 'weighted_balance_table' = weightd_table,
                 'TA_levels' = t_levels,
                 'ESS' = ESS))
}

#' @method summary entbal_binary
#' @export
summary.entbal_binary <- function(object, show_unweighted = TRUE, n_digits=2){

  cat('Reference levels for headers:\n')
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  if(is.factor(object$TA)){
    ta_lvls <- levels(object$TA)
    ta_ind <- ifelse(object$TA == object$eb_pars$which_z,1,0)
    cat(paste('Exposure 0:', ta_lvls[ta_lvls != object$eb_pars$which_z], '\n'))
    cat(paste('Exposure 1:', object$eb_pars$which_z, '\n'))
  } else {
    ta_lvls <- unique(object$TA)
    ta_ind <- ifelse(object$TA == min(ta_lvls), 0, 1)
    cat(paste('Exposure 0:', min(ta_lvls), '\n'))
    cat(paste('Exposure 1:', max(ta_lvls), '\n'))
  }

  outtab <- print_baltables(as.matrix(object$X),
                            ta_ind,
                            object$wts,
                            show_unweighted=show_unweighted,
                            n_digits=n_digits)

}

#' @method summary entbal_multiclass
#' @export
summary.entbal_multiclass <- function(object, show_unweighted = TRUE, n_digits = 2){
  estimand <- object$eb_pars$estimand
  ta_lvls <- unique(object$TA)
  NT <- length(ta_lvls)

  if(is.null(ncol(object$X))) object$X <- matrix(object$X, ncol = 1)

  outsum1 <- matrix(NA, nrow = ncol(object$X), ncol = 2 * (NT+1))
  balsum1 <- matrix(NA, nrow = ncol(object$X), ncol = 2 * NT)

  outsum2 <- matrix(NA, nrow = ncol(object$X), ncol = 2 * (NT+1))
  balsum2 <- matrix(NA, nrow = ncol(object$X), ncol = 2 * NT)

  orig_N <- rep(NA, NT)
  esssum <- rep(NA, NT)

  names(orig_N) <- paste("TA:", ta_lvls, sep = '')
  names(esssum) <- paste("TA:", ta_lvls, sep = '')
  colnames(outsum1) <- c('Mean', 'SD', paste(c('M:', 'SD:'), rep(ta_lvls,each=2),sep=''))
  colnames(outsum2) <- c('Mean', 'SD', paste(c('M:', 'SD:'), rep(ta_lvls,each=2),sep=''))
  colnames(balsum1) <- c(paste(c('M:', 'SD:'), rep(ta_lvls,each=2),sep=''))
  colnames(balsum2) <- c(paste(c('M:', 'SD:'), rep(ta_lvls,each=2),sep=''))
  rownames(outsum1) <- colnames(object$X)
  rownames(outsum2) <- colnames(object$X)
  rownames(balsum1) <- colnames(object$X)
  rownames(balsum2) <- colnames(object$X)

  if(estimand == 'ATE'){
    target_means <- apply(object$X, 2, mean)
    target_stddv <- apply(object$X, 2, sd)
  } else{
    ref_z <- object$eb_pars$which_z
    target_means <- apply(object$X[object$TA == ref_z,], 2, mean)
    target_stddv <- apply(object$X[object$TA == ref_z,], 2, sd)
  }

  outsum1[,1] <- target_means
  outsum1[,2] <- target_stddv
  outsum2[,1] <- target_means
  outsum2[,2] <- target_stddv
  for(i in 1:NT){
    ta <- object$TA == ta_lvls[i]
    nt <- sum(ta)
    group_means <- apply(object$X, 2, function(x) wtd_mean(x, ta, object$wts))
    group_stddv <- sqrt(apply(object$X, 2, function(x) wtd_sd2(x, ta, object$wts)))

    uw_group_means <- apply(object$X, 2, function(x) wtd_mean(x, ta, ta))
    uw_group_stddv <- sqrt(apply(object$X, 2, function(x) wtd_sd2(x, ta, ta)))

    outsum1[, 2*i + 1] <- uw_group_means
    outsum1[, 2*i + 2] <- uw_group_stddv
    balsum1[, 2*i - 1] <- (uw_group_means - target_means) / target_stddv
    balsum1[, 2*i] <- log(target_stddv) - log(uw_group_stddv)

    outsum2[, 2*i + 1] <- group_means
    outsum2[, 2*i + 2] <- group_stddv
    balsum2[, 2*i - 1] <- (group_means - target_means) / target_stddv
    balsum2[, 2*i] <- log(target_stddv) - log(group_stddv)

    esssum[i] <- group_ESS(object$wts, ta)
    orig_N[i] <- nt
  }

  if(show_unweighted){
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    cat('Unweighted Summary Statistics:\n')
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    print(round(outsum1, digits = n_digits))
  }
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  cat('Weighted Summary Statistics:\n')
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  print(round(outsum2, digits = n_digits))
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))

  if(show_unweighted){
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    cat('Unweighted Balance Statistics:\n')
    cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
    print(round(balsum1, digits = n_digits))
  }
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  cat('Weighted Balance Statistics:\n')
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  print(round(balsum2, digits = n_digits))
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))

  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  cat('Original & Effective Sample Sizes:\n')
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))
  SS <- cbind(orig_N, esssum, esssum/orig_N)
  colnames(SS) <- c('Orig N', 'ESS', 'Ratio')
  rownames(SS) <- paste("TA:", ta_lvls, sep = '')
  print(round(SS, digits = n_digits))
  cat(paste(paste(rep('-', 80), collapse = ''), '\n', sep=''))

  invisible(list('unweighted_summary' = outsum1,
                 'unweighted_balstats' = balsum1,
                 'weighted_summary' = outsum2,
                 'weighted_balstats' = balsum2,
                 'weighted_ess' = esssum,
                 'original_N' = orig_N))
}

#' Weighted KS Statistic
#'
#' @param X covariate vector
#' @param TA treatment indicator
#' @param wts weights
#' @param n_pts number of points to evaluate weighted ecdf
#'
#' @export
ks_test <- function(X, TA, wts=rep(1,length(X)), n_pts=100){
  xmin <- min(X)
  xmax <- max(X)
  int_pts <- seq(xmin, xmax, length.out = n_pts)
  t_wts <- ifelse(TA == 1, wts, 0)
  c_wts <- ifelse(TA == 0, wts, 0)
  t_fn <- wtd_ecdf(X, wts = t_wts)
  c_fn <- wtd_ecdf(X, wts = c_wts)
  return(max(abs(t_fn(int_pts) - c_fn(int_pts))))
}

.n_uniq <- function(x){length(unique(x))}

.find_continuous <- function(X){apply(as.matrix(X), 2, .n_uniq) > 2}

#' Make Design Matrix for Binary/Multiclass Entropy Balancing
#'
#' @param X matrix of covariates to balance
#' @param m number of moments
#'
#' @export
make_Xmat <- function(X, m = 1){
  if(is.null(ncol(X))){
    if(length(length(unique(X))) == 2){
      Xout <- X
    } else {
      nu <- length(unique(X))
      if(nu > m){
        Xout <- poly(X, m, raw = F, simple = T)
      } else {
        Xout <- poly(X, nu-1, raw = F, simple = T)
      }
    }
  } else {
    cont_var <- .find_continuous(X)
    X_con <- as.matrix(X[, cont_var])
    NC <- ncol(X_con)
    X_bin <- as.matrix(X[, !cont_var])
    Xout <- cbind(X_bin)
    if(NC > 0){
      if(m > 1){
        for(i in 1:NC){
          nu <- length(unique(X_con[,i]))
          if(nu > m){
            Xout <- cbind(Xout, poly(X_con[,i], m, raw = F, simple = T))
          } else {
            Xout <- cbind(Xout, poly(X_con[,i], nu-1, raw = F, simple = T))
          }
        }
      } else{
        Xout <- cbind(Xout, X_con)
      }
    }
  }
  Xout
}


#' Calculate covariate balance statistics
#'
#' @param Xmat A matrix of covariates
#' @param TA A vector of treatment assignments (i.e, 1 if treated, 0 in control)
#' @param wts A vector of estimated (or true) weights
#' @param show_unweighted Query if the the balance table should contain unweighted estimates
#' @param n_digits Number of digits to print for the table
#' @return A table of covariate balance statistics
#'
#' @export
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
