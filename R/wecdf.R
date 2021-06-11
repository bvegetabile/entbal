#' Weighted ECDF
#'
#' @param var_data vector of variable data
#' @param wts vector of weights
#'
#' @export
wtd_ecdf <- function (var_data, wts) {
  #-----------------------------------------------------------------------------
  # wtd_ecdf is a modification of the ecdf() function in base R.  It modifies
  # the function to be able to incorporate weights.  This is to visualize
  # balance using the empirical cumulative distribution function for continuous
  # covariates after weighting by the inverse of the propensity score (IPTW)
  #
  # Input variables
  # --- var_data : covariate values - vector of data
  # --- wts      : weights for assessing cov balance by IPTW - vector of data.
  #-----------------------------------------------------------------------------
  ord <- order(var_data)
  var_ordered <- var_data[ord]
  wts_ordered <- wts[ord]
  n <- length(var_data)
  if (n < 1)
    stop("'var_data' must have 1 or more non-missing values")
  vals <- unique(var_ordered)
  matched_vals <- match(var_ordered, vals)
  weight_list <- aggregate(wts_ordered, by=list(matched_vals), sum)
  rval <- approxfun(vals, cumsum(weight_list[,2])/sum(wts_ordered),
                    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
