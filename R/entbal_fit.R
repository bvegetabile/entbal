# entbal_wts <- function(Q, C, Z){
#   norm_c <- Q %*% exp( - C %*% Z )
#   Q * exp( - C %*% Z ) / c(norm_c)
# }

entbal_wts <- function(Q, C, Z){
  V <- - C %*% Z
  maxV <- max(V)
  norm_c <- Q %*% exp( V - maxV )
  Q * exp( V - maxV ) / c(norm_c)
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
    XS = - C %*% f
    maxXS = max(XS)
    loss <- maxXS + log(t(Q) %*% exp( XS - maxXS)) + t(M) %*% f
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
