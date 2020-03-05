makeC2 <- function(XD, A, n_moments = 3){
  NC <- ncol(XD)
  NR <- nrow(XD)
  outmat <- matrix(NA, nrow = NR, ncol = 0)

  # correlation balance
  col_ids <- c()
  for(c in 1:NC){
    nuniq <- length(unique(XD[,c]))
    if(nuniq <= 1) {
      message(paste('Column: ',
                    colnames(XD)[c],
                    ', has <= 1 unique value and was not included'))
    } else if(nuniq == 2) {
      outmat <- cbind(outmat, scale(XD[,c]) * scale(A))
      col_ids <- c(col_ids, paste('corbal_', colnames(XD)[c], sep = ''))
    } else {
      for(p in 1:n_moments){
        outmat <- cbind(outmat, scale(XD[,c]**p) * scale(A))
      }
      col_ids <- c(col_ids, paste('corbal_',colnames(XD)[c], '_poly', 1:n_moments, sep =''))
    }
  }

  # marginal balance - exposure
  outmat <- cbind(outmat, poly(A, n_moments))
  col_ids <- c(col_ids, paste('margbal_exposure_poly', 1:n_moments,sep =''))

  # marginal balance - covariates
  for(c in 1:NC){
    nuniq <- length(unique(XD[,c]))
    if(nuniq == 2 ) {
      outmat <- cbind(outmat, scale(XD[,c]))
      col_ids <- c(col_ids, paste('margbal_', colnames(XD)[c], sep = ''))
    } else {
      outmat <- cbind(outmat, poly(XD[,c], n_moments))
      col_ids <- c(col_ids, paste('margbal_', colnames(XD)[c], '_poly', 1:n_moments,sep =''))
    }

  }
  colnames(outmat) <- col_ids
  outmat
}



entbal_cont <- function(formula,
                        data = NULL,
                        n_moments = 3,
                        max_iters = 1000,
                        verbose = FALSE,
                        optim_method = 'L-BFGS-B',
                        bal_tol = 1e-8,
                        opt_constraints = c(-100, 100)){

  # Cleaning up user input
  # estimand <- toupper(estimand)
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

  # if(n_classes==n_obs){stop('Number of unique treatment values equals the number of observations\n  -->Continuous treatment regimes not currently supported')}
  if(n_classes==1){stop('Number of unique treatment values is one\n  -->Single treatment value is incoherent')}
  if(n_classes==2){stop('Number of unique treatment values is two\n  -->Use binary entropy weighting')}
  if(n_classes<=5 & n_classes > 2){warning('Number of unique treatment values is between 2 and 5\n  --> Multi-valued entropy balancing may be better')}

  Xmat <- makeC2(XD = designX[,2:NC],
                 A = ta,
                 n_moments = n_moments)

  eb_wts <- entbal_fit(C = Xmat,
                       targets = rep(0, ncol(Xmat)),
                       n_moments = 1,
                       max_iters = max_iters,
                       verbose = verbose,
                       optim_method = optim_method,
                       bal_tol = bal_tol,
                       opt_constraints = opt_constraints)

  conv_status_control <- eb_wts$optim_obj$convergence
  conv_status <- ifelse(conv_status_control==0, T, F)
  mf$wts <- eb_wts$wts

  if(conv_status == F) {
    warning('Check convergence status')
  }
  res <- list('wts' = mf$wts,
              'convergence' = conv_status,
              'message' = eb_wts$optim_obj$message,
              'n_matched_moments' = n_moments,
              'X' = designX[,2:NC],
              'XD' = Xmat,
              'TA' = ta,
              'constraints' = opt_constraints,
              'opt_obj' = eb_wts$optim_obj)
  class(res) <- c(if(n_classes > 2) "entbal_cont", "entbal_binary")
  # res$estimand <- estimand
  res
}
