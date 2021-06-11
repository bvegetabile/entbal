entbal_mc <- function(formula,
                      data = NULL,
                      R = NULL,
                      which_Z = NULL,
                      estimand = "ATE",
                      n_moments = 3,
                      max_iters = 1000,
                      verbose = FALSE,
                      optim_method = 'L-BFGS-B',
                      bal_tol = 0.005){


  # Cleaning up user input
  estimand <- toupper(estimand)
  formula <- formula(formula)

  # Checking if the formula has a response
  if(!attr(terms(formula, data=data), 'response')) stop('Please supply a treatment variable on the left side of the formula');

  # Collecting the data and making a model.frame object to create the design matrix
  mf <- model.frame(formula, data = data)

  ta <- model.response(mf, 'numeric')
  designX <- model.matrix(formula, data=mf)


  uniq_ta <- as.vector(unique(ta))
  uniq_ta <- uniq_ta[order(uniq_ta)]
  N_uniq <- length(uniq_ta)
  NC <- ncol(designX)

  n_classes <- length(unique(ta))
  n_obs <- nrow(designX)

  if(n_classes==n_obs){stop('Number of unique treatment values equals the number of observations\n  -->Continuous treatment regimes not currently supported')}
  if(n_classes==1){stop('Number of unique treatment values is one\n  -->Single treatment value is incoherent')}

  if(!estimand %in% c('ATE', 'ATZ')){stop('Invalid estimand: Choose ATE, ATZ')}


  mf$wts <- NA
  conv_messages <- list()
  opt_obj <- list()
  all_converged <- T
  if(estimand == 'ATE'){
    Xmat <- make_Xmat(designX[,2:NC], n_moments)
    Xmat <- scale(Xmat)
    targets <- apply(Xmat, 2, mean)

    for(Z in 1:N_uniq){
      z = uniq_ta[Z]
      if(is.null(R)){
        XZ <- Xmat[ta == z, ]
      } else {
        XZ <- Xmat[ta == z & R == 1, ]
        mf$wts[R==0] <- 0
      }
      wtsZ <- entbal_fit(XZ, targets, n_moments, max_iters, verbose, optim_method, bal_tol)
      conv_status_z <- wtsZ$optim_obj$convergence
      conv_status <- ifelse(conv_status_z ==0, T, F)
      all_converged <- conv_status & all_converged
      conv_messages[[Z]] <- wtsZ$optim_obj$message
      opt_obj[[Z]] <- wtsZ$optim_obj
      if(is.null(R)){
        mf$wts[ta == z] <- wtsZ$wts
      } else {
        mf$wts[ta == z & R == 1] <- wtsZ$wts
      }
    }
  } else if (estimand == 'ATZ'){
    Xmat <- make_Xmat(designX[,2:NC], n_moments)
    Xmat <- scale(Xmat)
    targets <- apply(Xmat, 2, mean)

    Xmat <- make_Xmat(designX[,2:NC], n_moments)
    Xmat <- scale(Xmat)
    XZ <- Xmat[ta == which_Z, ]
    targets <- apply(XZ, 2, mean)

    for(Z in 1:N_uniq){
      z = uniq_ta[Z]
      if(z == which_Z){
        conv_status_z <- NA
        conv_status <- T
        all_converged <- conv_status & all_converged
        conv_messages[[Z]] <- NA
        opt_obj[[Z]] <- NA
        mf$wts[ta == z] <- 1
      } else {
        XZ <- Xmat[ta == z, ]
        wtsZ <- entbal_fit(XZ, targets, n_moments, max_iters, verbose, optim_method, bal_tol)
        conv_status_z <- wtsZ$optim_obj$convergence
        conv_status <- ifelse(conv_status_z ==0, T, F)
        all_converged <- conv_status & all_converged
        conv_messages[[Z]] <- wtsZ$optim_obj$message
        opt_obj[[Z]] <- wtsZ$optim_obj
        mf$wts[ta == z] <- wtsZ$wts
      }
    }


  } else {
    stop("You shouldn't be here")
  }

  if(all_converged == F) {
    warning('Check convergence status')
  }
  res <- list('wts' = mf$wts,
              'convergence' = conv_status,
              'message' = conv_messages,
              'n_matched_moments' = n_moments,
              'estimand' = estimand,
              'X' = designX[,2:NC],
              'TA' = ta,
              'ref_z' = which_Z,
              'TA_map' = data.frame('TA_levels' = uniq_ta,
                                    'Z_levels' = 1:N_uniq),
              'opt_obj' = opt_obj)
  class(res) <- c(if(n_classes > 2) "entbal_multiclass", "entbal_binary")
  res$estimand <- estimand
  res
}
