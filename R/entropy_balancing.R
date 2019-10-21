entbal_wts <- function(Q, C, Z){
  norm_c <- Q %*% exp( - C %*% Z )
  Q * exp( - C %*% Z ) / c(norm_c)
}

entbal_fit <- function(C, targets,
                       n_moments = 2,
                       max_iters = 1000,
                       verbose = 0,
                       optim_method = 'L-BFGS-B'){
  n_obs <- nrow(C)
  Q <- rep(1/n_obs, n_obs)
  M <- targets
  n_targets <- length(M)

  loss_func0 <- function(f){
    loss <- log(t(Q) %*% exp( - C %*% f )) + t(M) %*% f
    return(loss)
  }

  grad_func0 <- function(f){
    W <- entbal_wts(Q, C, f)
    grad <- M - t(C) %*% W
    return(grad)
  }

  f_init <- solve(t(C) %*% C) %*% M

  opt_val <- optim(par = f_init,
                   fn = loss_func0,
                   gr = grad_func0,
                   method = optim_method,
                   control = list(trace = verbose,
                                  maxit = max_iters,
                                  lmm = 25))
  return(list(optim_obj = opt_val,
              f = opt_val$par,
              wts = entbal_wts(Q, C, opt_val$par)))

}

entbal <- function(formula,
                   data = NULL,
                   R = NULL,
                   estimand = "ATE",
                   n_moments = 3,
                   max_iters = 1000,
                   verbose = FALSE,
                   lambda = NULL,
                   optim_method = 'L-BFGS-B'){

  # Cleaning up user input
  estimand <- toupper(estimand)
  formula <- formula(formula)

  # Checking if the formula has a response
  if(!attr(terms(formula, data=data), 'response')) stop('Please supply a treatment variable on the left side of the formula');

  # Dropping the intercept term
  if(attr(terms(formula, data=data), 'intercept')){
    formula <- update(terms(formula, data=data), . ~ . -1)
  }

  # Collecting the data and making a model.frame object to create the design matrix
  mf <- model.frame(formula, data = data)

  ta <- model.response(mf, 'numeric')
  designX <- model.matrix(formula, data=mf)

  n_classes <- length(unique(ta))
  n_obs <- nrow(designX)

  if(n_classes==n_obs){stop('Number of unique treatment values equals the number of observations\n  -->Continuous treatment regimes not currently supported')}
  if(n_classes==1){stop('Number of unique treatment values is one\n  -->Single treatment value is incoherent')}

  if(!estimand %in% c('ATE', 'ATT', 'ATC')){stop('Invalid estimand: Choose ATE, ATT, or ATC')}

  if(estimand == 'ATE'){
    Xmat <- make_Xmat(designX, n_moments)
    Xmat <- scale(Xmat)
    targets <- apply(Xmat, 2, mean)

    if(is.null(R)){
      XT <- Xmat[ta == 1, ]
      XC <- Xmat[ta == 0, ]
    } else {
      XT <- Xmat[ta == 1 & R == 1, ]
      XC <- Xmat[ta == 0 & R == 1, ]
    }

    wtsT <- entbal_fit(XT, targets, n_moment, max_iters, verbose, optim_method)
    wtsC <- entbal_fit(XC, targets, n_moment, max_iters, verbose, optim_method)

    conv_status_treated <- wtsT$optim_obj$convergence
    conv_status_control <- wtsC$optim_obj$convergence

    conv_status <- ifelse(conv_status_treated==0, T, F) & ifelse(conv_status_control==0, T, F)
    conv_messages <- list('treated' = wtsT$optim_obj$message,
                          'control' = wtsC$optim_obj$message)

    mf$wts <- NA
    if(is.null(R)){
      mf$wts[ta == 1] <- wtsT$wts
      mf$wts[ta == 0] <- wtsC$wts
    } else {
      mf$wts[ta == 1 & R == 1] <- wtsT$wts
      mf$wts[ta == 0 & R == 1] <- wtsC$wts
      mf$wts[R==0] <- 0
    }


  } else if (estimand == 'ATT'){
    Xmat <- make_Xmat(designX, n_moments)
    Xmat <- scale(Xmat)
    XT <- Xmat[ta == 1, ]
    targets <- apply(XT, 2, mean)

    XC <- Xmat[ta == 0, ]

    wtsC <- entbal_fit(XC, targets, n_moment, max_iters, verbose, optim_method)

    conv_status_control <- wtsC$optim_obj$convergence

    conv_status <- ifelse(conv_status_control==0, T, F)
    conv_messages <- list('control' = wtsC$optim_obj$message)

    mf$wts <- NA
    mf$wts[ta == 1] <- 1
    mf$wts[ta == 0] <- wtsC$wts

  }

  if(conv_status == F) {
    warning('Check convergence status')
  }
  res <- list('wts' = mf$wts,
              'convergence' = conv_status,
              'message' = conv_messages,
              'n_matched_moments' = n_moments,
              'X' = designX,
              'TA' = ta)
  class(res) <- c(if(n_classes > 2) "entbal_multiclass", "entbal_binary")
  res$estimand <- estimand
  res
}
