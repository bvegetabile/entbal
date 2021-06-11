#' Create list of default parameters for entbal - binary
#'
#' @param estimand default: \code{ATE}, alternative is \code{ATT}
#' @param which_z if \code{estimand} is \code{ATT}, then set \code{which_z} to target value from treatment vector
#'
#' @export
ebpars_default_binary <- function(estimand = 'ATE',
                                  which_z = NULL){
  outlist <- list('exp_type' = 'binary',
                  'n_moments' = 3,
                  'max_iters' = 1000,
                  'estimand' = estimand,
                  'verbose' = T,
                  'optim_method' = 'l-bfgs-b',
                  'bal_tol' = 1e-8,
                  'opt_constraints' = c(-1000,1000))
  if(!is.null(which_z)){
    outlist$which_z <- which_z
  }
  outlist
}

#' Create list of default parameters for entbal - multi
#'
#' @param estimand default: \code{ATE}, alternative is \code{ATT}
#' @param which_z if \code{estimand} is \code{ATT}, then set \code{which_z} to target value from treatment vector
#'
#' @export
ebpars_default_multi <- function(estimand = 'ATE',
                                  which_z = NULL){
  outlist <- ebpars_default_binary(estimand, which_z)
  outlist
}

#' Create list of default parameters for entbal - continuous
#'
#' @export
ebpars_default_cont <- function(){
  outlist <- list('exp_type' = 'continuous',
                  'n_moments' = 3,
                  'max_iters' = 1000,
                  'estimand' = 'ATE',
                  'verbose' = T,
                  'optim_method' = 'l-bfgs-b',
                  'bal_tol' = 1e-8,
                  'opt_constraints' = c(-1000,1000))
  outlist
}

