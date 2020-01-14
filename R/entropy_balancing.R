entbal_wts <- function(Q, C, Z){
  norm_c <- Q %*% exp( - C %*% Z )
  Q * exp( - C %*% Z ) / c(norm_c)
}

entbal_fit <- function(C, targets,
                       n_moments = 2,
                       max_iters = 1000,
                       verbose = 0,
                       optim_method = 'L-BFGS-B',
                       bal_tol = 0.0005,
                       opt_constraints = c(-100, 100)){
  n_obs <- nrow(C)
  Q <- rep(1/n_obs, n_obs)
  M <- targets
  n_targets <- length(M)

  loss_func0 <- function(f){
    # W <- entbal_wts(Q, C, f)
    # print(max(M - t(C) %*% W))
    loss <- log(t(Q) %*% exp( - C %*% f )) + t(M) %*% f
    # print(t(Q) %*% exp( - C %*% f ))
    return(loss)
  }

  grad_func0 <- function(f){
    W <- entbal_wts(Q, C, f)
    grad <- M - t(C) %*% W
    return(grad)
  }

  f_init <- solve(t(C) %*% C + diag(ncol(C))) %*% M

  if(optim_method == 'L-BFGS-B'){
    opt_val <- optim(par = f_init,
                     fn = loss_func0,
                     gr = grad_func0,
                     method = optim_method,
                     lower = opt_constraints[1],
                     upper = opt_constraints[2],
                     control = list(trace = verbose,
                                    maxit = max_iters,
                                    lmm = 5,
                                    pgtol = bal_tol))
  } else if (optim_method == 'BFGS') {
    opt_val <- optim(par = f_init,
                     fn = loss_func0,
                     gr = grad_func0,
                     method = optim_method,
                     control = list(trace = verbose,
                                    maxit = max_iters))
  } else {
    stop('Unknown optimization method: Only L-BFGS-B and BFGS supported at this point')
  }

  return(list(optim_obj = opt_val,
              f = opt_val$par,
              wts = entbal_wts(Q, C, opt_val$par)))

}

#' Compute optimal balancing weights via entropy balancing
#'
#' @param formula Typical R style formula - ex. `TA ~ X1 + X2`
#' @param data R \code{data.frame} that contains the variables listed in the formula
#' @param R Binary vector of response indictors
#' @param estimand Estimand for optimizing the weights towards.  Acceptable values of "ATE", "ATT", "ATC" where ATT takes values of 1 to be the treatment category and 0 to be the control category
#' @param n_moments The number of moments that are to be matched.  Default to 2
#' @param max_iters Maximum number of iterations for the BFGS optimization procedure
#' @param verbose Should the algorithm print progess to screen? Defaults to FALSE
#' @param optim_method Optimization Method.  Set to 'L-BFGS-B' by deault.  Specify 'BFGS' is there are optimization issues
#' @return Object that contains the weights obtained from the balancing procedure and parameters from the optimization procedure
#'
#' The object that is returned is a list that contains the following entries
#' \itemize{
#' \item{ \code{wts} - Optimal weights for the estimand of interest and matched number of moments.}
#' \item{ \code{untransformed_wts} - Latent value used to contruct the weights.  Typically unneeded.}
#' \item{ \code{minimum_overall_bal} - Resulting value of the objective function at algorithm termination.}
#' \item{ \code{convergence} - Convergence code from \code{optim} package for BFGS algorithm.}
#' \item{ \code{message} - Convergence message from \code{optim} package for BFGS algorithm.}
#' \item{ \code{n_matched_moments} - Number of moments matched in the optimization.}
#' \item{ \code{X} - Design matrix used in obtaining balance.}
#' \item{ \code{TA} - Vector of Treatment Assignments Used.}
#' }
#' @examples
#' # Example 1 - ATE Example
#' n_obs <- 500
#' X1 <- rnorm(n_obs)
#' X2 <- rnorm(n_obs)
#' p <- pnorm( 0.5 * X1 + 0.5 * X2 )
#' TA <- rbinom(n_obs, 1, p)
#' dat <- data.frame(X1 = X1, X2 = X2, TA = TA)
#' system.time(res <- entbald('TA ~ X1 + X2',
#'                           data = dat, verbose = T, optim_method = 'BFGS'))
#' summary(res)
#'
#' # Example 2 - ATT Example
#' dset <- dw99cps1
#' estwts <- entbal(TA ~ age + education + black + hispanic + married + nodegree + RE74 + RE75,
#'                       data = dset,
#'                       n_moments = 3,
#'                       verbose = 1,
#'                       estimand = "ATT")
#' summary(estwts)
#' dset$wts <- estwts$wts
#' design <- svydesign(ids=~1, weights=~wts, data = dset)
#' resp <- svyglm('RE78 ~ TA', design = design)
#'summary(resp)
#'
entbal <- function(formula,
                   data = NULL,
                   R = NULL,
                   estimand = "ATE",
                   n_moments = 3,
                   max_iters = 1000,
                   verbose = FALSE,
                   optim_method = 'L-BFGS-B',
                   bal_tol = 0.005,
                   opt_constraints = c(-100, 100)){

  # Cleaning up user input
  estimand <- toupper(estimand)
  formula <- formula(formula)

  # Checking if the formula has a response
  if(!attr(terms(formula, data=data), 'response')) stop('Please supply a treatment variable on the left side of the formula');

  # Dropping the intercept term
  # if(attr(terms(formula, data=data), 'intercept')){
  #   formula <- update(terms(formula, data=data), . ~ . -1)
  # }

  # Collecting the data and making a model.frame object to create the design matrix
  mf <- model.frame(formula, data = data)

  ta <- model.response(mf, 'numeric')
  designX <- model.matrix(formula, data=mf)

  NC <- ncol(designX)

  n_classes <- length(unique(ta))
  n_obs <- nrow(designX)

  if(n_classes==n_obs){stop('Number of unique treatment values equals the number of observations\n  -->Continuous treatment regimes not currently supported')}
  if(n_classes==1){stop('Number of unique treatment values is one\n  -->Single treatment value is incoherent')}

  if(!estimand %in% c('ATE', 'ATT', 'ATC')){stop('Invalid estimand: Choose ATE, ATT, or ATC')}

  if(estimand == 'ATE'){
    Xmat <- make_Xmat(designX[,2:NC], n_moments)
    Xmat <- scale(Xmat)
    targets <- apply(Xmat, 2, mean)

    if(is.null(R)){
      XT <- Xmat[ta == 1, ]
      XC <- Xmat[ta == 0, ]
    } else {
      XT <- Xmat[ta == 1 & R == 1, ]
      XC <- Xmat[ta == 0 & R == 1, ]
    }

    wtsT <- entbal_fit(XT, targets, n_moment, max_iters, verbose, optim_method, bal_tol, opt_constraints)
    wtsC <- entbal_fit(XC, targets, n_moment, max_iters, verbose, optim_method, bal_tol, opt_constraints)

    conv_status_treated <- wtsT$optim_obj$convergence
    conv_status_control <- wtsC$optim_obj$convergence

    conv_status <- ifelse(conv_status_treated==0, T, F) & ifelse(conv_status_control==0, T, F)
    conv_messages <- list('treated' = wtsT$optim_obj$message,
                          'control' = wtsC$optim_obj$message)

    opt_obj <- list('treated' = wtsT$optim_obj,
                    'control' = wtsC$optim_obj)


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
    Xmat <- make_Xmat(designX[,2:NC], n_moments)
    Xmat <- scale(Xmat)
    XT <- Xmat[ta == 1, ]
    targets <- apply(XT, 2, mean)
    XC <- Xmat[ta == 0, ]

    wtsC <- entbal_fit(XC, targets, n_moment, max_iters, verbose, optim_method, bal_tol, opt_constraints)

    conv_status_control <- wtsC$optim_obj$convergence

    conv_status <- ifelse(conv_status_control==0, T, F)
    conv_messages <- list('control' = wtsC$optim_obj$message)

    opt_obj <- list('control' = wtsC$optim_obj)

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
              'estimand' = estimand,
              'X' = designX[,2:NC],
              'TA' = ta,
              'opt_obj' = opt_obj)
  class(res) <- c(if(n_classes > 2) "entbal_multiclass", "entbal_binary")
  res$estimand <- estimand
  res
}
