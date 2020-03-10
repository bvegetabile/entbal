
.check_pars <- function(parlist, nc, no, uc){

  par_names <- names(parlist)

  if(!('exp_type' %in% par_names)) {
    if(nc == 2){
      parlist$exp_type <- 'binary'
      warning('Exposure type not set. 2 exposure levels found. Using exp_type = binary')
    } else if(nc > 2 & nc <= 5){
      parlist$exp_type <- 'multi'
      warning(sprintf('Exposure type not set. %s exposure levels found. Using exp_type = multi', nc))
    }
    else if(nc >= no/2){
      parlist$exp_type <- 'continuous'
      warning('Exposure type not set. Number of unique obs > 0.5 * n_obs. Using exp_type = continuous')
    } else {
      stop('Unable able to automatically detect exposure type.  Requires variable "exp_type": binary, multi, continuous')
    }
  } else{
    if(nc == 2 & parlist$exp_type != 'binary'){
      parlist$exp_type <- 'binary'
      warning('Exposure type mismatch. Only 2 exposure levels found. Using exp_type = binary')
    } else if((nc > 2 & nc <= 5) & (parlist$exp_type != 'multi')){
      parlist$exp_type <- 'multi'
      warning(sprintf('Exposure type mismatch. %s exposure levels found. Using exp_type = multi', nc))
    }
    else if(nc >= no/2 & parlist$exp_type != 'continuous'){
      parlist$exp_type <- 'continuous'
      warning('Exposure type mismatch. Number of unique obs > 0.5 * n_obs. Using exp_type = continuous')
    } else {
      stop('Unable able to automatically detect exposure type.  Requires variable "exp_type": binary, multi, continuous')
    }

  }

  parlist$exp_type <- tolower(parlist$exp_type)

  # # Checking that minimum require variables are included
  # min_required <- c('n_moments',
  #                   'max_iters',
  #                   'verbose',
  #                   'optim_method')
  #
  # var_check <- min_required %in% par_names
  #
  # if(sum(var_check) != 4) stop(paste('\nMissing required variables in eb_pars: ',
  #                                    paste(min_required[!var_check], collapse = ', '),
  #                                    sep = ' '))

  # Checking optimization variables
  if( ! 'n_moments' %in% par_names){
    parlist$n_moments <- 3
    warning('Number of balanced moments not set. Using n_moments = 3')
  }

  if( ! 'max_iters' %in% par_names){
    parlist$max_iters <- 1000
    warning('Max number of iterations not set. Using max_iter = 1000')
  }

  if( ! 'optim_method' %in% par_names) {
    parlist$optim_method = 'L-BFGS-B'
    warning('Optimization method not set. Using optim_method = L-BFGS-B')
  }

  if( ! 'verbose' %in% par_names) {
    parlist$verbose = F
    warning('Verbose method not set. Using verbose = FALSE')
  }

  parlist$optim_method <- toupper(parlist$optim_method)

  if(!parlist$optim_method %in% c('BFGS', 'L-BFGS-B')){
    stop('Invalid optimization method: Choose BFGS or L-BFGS-B')
  }

  if(parlist$optim_method == 'BFGS') {
    parlist$bal_tol <- 0.005
    parlist$opt_constraints <- c(-100,100)
    if('bal_tol' %in% par_names) warning('BFGS chosen... bal_tol ignored')
    if('opt_constraints' %in% par_names) warning('BFGS chosen... opt_constraints ignored')
  }

  if(parlist$optim_method == 'L-BFGS-B'){
    if( !'bal_tol' %in% par_names) {
      parlist$bal_tol <- 0.005
      warning('bal_tol not specified. set to 0.005')
    }
    if('opt_constraints' %in% par_names) {
      if(length(parlist$opt_constraints) != 2) stop('Issue in opt_constraints.  Expecting a list of two numbers.  Good default is c(-100, 100)')
    } else {
      parlist$opt_constraints <- c(-100, 100)
      warning('opt_constraints not specified. set to c(-100, 100)')
    }
  }


  # Checking entbal parameters are included

  if(parlist$exp_type == 'binary'){

    if(!('estimand' %in% par_names)) {
      parlist$estimand <- 'ATE'
      warning('No estimand chosen. Estimand set to default: ATE')
    } else {
      parlist$estimand <- toupper(parlist$estimand)
      if(!parlist$estimand %in% c('ATE', 'ATT', 'ATC')){
        stop('Invalid estimand: Choose ATE, ATT, or ATC')
      }
    }

    if(!('R' %in% par_names)) {
      parlist$R  <- NULL
    } else{
      warning('Response indicators included')
    }

  } else if(parlist$exp_type == 'multi') {

    if(!('estimand' %in% par_names)) {
      parlist$estimand <- 'ATE'
      warning('No estimand chosen. Estimand set to default: ATE')
    } else {
      parlist$estimand <- toupper(parlist$estimand)
      if(parlist$estimand %in% c('ATT', 'ATC')){
        parlist$estimand <- 'ATZ'
        warning('ATT or ATC not coherent for exp_type = "multi". Converted to ATZ.')
      }

      if(!parlist$estimand %in% c('ATE', 'ATZ')){
        stop('Invalid estimand: Choose ATE, ATZ')
      }
    }

    if((parlist$exp_type == 'multi' | parlist$exp_type == 'binary') & !(parlist$estimand == 'ATE')){
      if( !('which_z' %in% par_names)) {
        parlist$which_z <- ifelse(is.factor(uc),
                                  min(levels(uc)),
                                  min(uc))
        warning(paste('Reference exposure not chosen: eb_pars$which_z set to', parlist$which_z))
      } else{
        if(!parlist$which_z %in% uc) {
          stop(paste('Chosen exposure reference level not in unique levels. Levels found:',
                     paste(uc, collapse =', ')))
        }
      }

    }




    if(!('R' %in% par_names)) {
      parlist$R  <- NULL
    } else{
      warning('Response indicators included')
    }

  } else if(parlist$exp_type == 'continuous') {

  } else {
    stop('Invalid choice of "exp_type", choose: binary, multi, continuous')
  }

  parlist
}

# R = NULL,
# estimand = "ATE",
# n_moments = 3,
# max_iters = 1000,
# verbose = FALSE,
# optim_method = 'L-BFGS-B',
# bal_tol = 0.005,
# opt_constraints = c(-100, 100)

entbal_dev <- function(formula,
                       data = NULL,
                       eb_pars = list('exp_type' = 'binary',
                                      'n_moments' = 3,
                                      'max_iters' = 1000,
                                      'verbose' = FALSE,
                                      'optim_method' = 'L-BFGS-B'),
                       suppress_warnings = F){



  # Cleaning up user input
  formula <- formula(formula)

  # Checking if the formula has a response
  if(!attr(terms(formula, data=data), 'response')) stop('Please supply a treatment variable on the left side of the formula');

  # Collecting the data and making a model.frame object to create the design matrix
  mf <- model.frame(formula, data = data)
  ta <- model.response(mf)
  design_x <- model.matrix(formula, data=mf)

  nc <- ncol(design_x)
  uniq_ta <- unique(ta)
  n_classes <- length(uniq_ta)
  n_obs <- nrow(design_x)

  if(suppress_warnings){
    suppressWarnings(eb_pars <- .check_pars(eb_pars, n_classes, n_obs, uniq_ta))
  } else {
    eb_pars <- .check_pars(eb_pars, n_classes, n_obs, uniq_ta)
  }


  if(eb_pars$exp_type == 'binary' | eb_pars$exp_type == 'multi'){

    mf$wts <- NA
    conv_messages <- list()
    opt_obj <- list()
    all_converged <- T

    if(eb_pars$estimand == 'ATE') {

      x_mat <- make_Xmat(design_x[,2:nc], eb_pars$n_moments)
      x_mat <- scale(x_mat)
      targets <- apply(x_mat, 2, mean)
      for(k in 1:n_classes){
        z = uniq_ta[k]
        if(is.null(eb_pars$R)){
          xz <- x_mat[ta == z, ]
        } else {
          xz <- x_mat[ta == z & eb_pars$R == 1, ]
          mf$wts[eb_pars$R==0] <- 0
        }
        wts_z <- entbal_fit(C = xz,
                            targets = targets,
                            n_moments = eb_pars$n_moments,
                            max_iters = eb_pars$max_iters,
                            verbose = eb_pars$verbose,
                            optim_method = eb_pars$optim_method,
                            bal_tol = eb_pars$bal_tol,
                            opt_constraints = eb_pars$opt_constraints)

        conv_status_z <- wts_z$optim_obj$convergence
        conv_status <- ifelse(conv_status_z == 0, T, F)
        all_converged <- conv_status & all_converged
        conv_messages[[k]] <- wts_z$optim_obj$message
        opt_obj[[k]] <- wts_z$optim_obj
        if(is.null(eb_pars$R)){
          mf$wts[ta == z] <- wts_z$wts
        } else {
          mf$wts[ta == z & eb_pars$R == 1] <- wts_z$wts
        }
      }

    } else {
      x_mat <- make_Xmat(design_x[,2:nc], eb_pars$n_moments)
      x_mat <- scale(x_mat)
      xz <- x_mat[ta == eb_pars$which_z, ]
      targets <- apply(xz, 2, mean)

      for(k in 1:n_classes){
        z = uniq_ta[k]
        if(z == eb_pars$which_z){
          conv_status_z <- NA
          conv_status <- T
          all_converged <- conv_status & all_converged
          conv_messages[[k]] <- NA
          opt_obj[[k]] <- NA
          mf$wts[ta == z] <- 1
        } else {
          xz <- x_mat[ta == z, ]
          wts_z <- entbal_fit(C = xz,
                              targets = targets,
                              n_moments = eb_pars$n_moments,
                              max_iters = eb_pars$max_iters,
                              verbose = eb_pars$verbose,
                              optim_method = eb_pars$optim_method,
                              bal_tol = eb_pars$bal_tol,
                              opt_constraints = eb_pars$opt_constraints)
          conv_status_z <- wts_z$optim_obj$convergence
          conv_status <- ifelse(conv_status_z ==0, T, F)
          all_converged <- conv_status & all_converged
          conv_messages[[k]] <- wts_z$optim_obj$message
          opt_obj[[k]] <- wts_z$optim_obj
          mf$wts[ta == z] <- wts_z$wts
        }
      }
    }
    if(all_converged == F) {
      warning('Check convergence status')
    }
  } else{
    x_mat <- makeC2(XD = design_x[,2:nc],
                    A = ta,
                    n_moments = eb_pars$n_moments)

    eb_wts <- entbal_fit(C = x_mat,
                         targets = rep(0, ncol(x_mat)),
                         n_moments = 1,
                         max_iters = eb_pars$max_iters,
                         verbose = eb_pars$verbose,
                         optim_method = eb_pars$optim_method,
                         bal_tol = eb_pars$bal_tol,
                         opt_constraints = eb_pars$opt_constraints)

    opt_obj <- eb_wts$optim_obj
    conv_messages <- eb_wts$optim_obj$message
    conv_status_control <- eb_wts$optim_obj$convergence
    conv_status <- ifelse(conv_status_control==0, T, F)
    mf$wts <- eb_wts$wts

  }

  res <- list('wts' = mf$wts,
              'convergence' = conv_status,
              'message' = conv_messages,
              'opt_obj' = opt_obj,
              'eb_pars' = eb_pars,
              'X' = design_x[,2:nc],
              'TA' = ta)
  class(res) <- c(
    if(eb_pars$exp_type == 'continuous') "entbal_cont",
    if(eb_pars$exp_type == 'multi') "entbal_multiclass",
    "entbal_binary")
  res
}


