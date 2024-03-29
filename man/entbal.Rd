% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy_balancing_dev.R
\name{entbal}
\alias{entbal}
\title{Compute optimal balancing weights via entropy balancing}
\usage{
entbal(
  formula,
  data = NULL,
  eb_pars = list(exp_type = "binary", n_moments = 3, max_iters = 1000, verbose = FALSE,
    bal_tol = 1e-05, opt_constraints = c(-100, 100), estimand = "ATE", optim_method =
    "L-BFGS-B"),
  suppress_warnings = F
)
}
\arguments{
\item{formula}{Typical R style formula - ex. `TA ~ X1 + X2`}

\item{data}{R \code{data.frame} that contains the variables listed in the formula}

\item{eb_pars}{R \code{list} of parameters required for entropy balancing.

The list recommends including the following variables at a minimum:
\itemize{
  \item{\code{exp_type}}: choose from \code{binary}, \code{multi}, \code{continuous}
  \item{\code{n_moments}}: number of moments to match in the entropy balancing procedure.  Recommended \code{3}
  \item{\code{max_iters}}: maximum number of iterations for the optimization routine.  Recommended \code{1000}
  \item{\code{verbose}}: logical for if the optimization should print information to the screen
  \item{\code{optim_method}}: Choose \code{L-BFGS-B} or \code{BFGS}.  Recommendation is \code{L-BFGS-B}
}
The function will attempt to set values to default values if variables are not specified and provide warnings when necessary. To ignore these warnings set the \code{suppress_warnings} variable to \code{TRUE}.

For \code{exp_type = 'binary'} or \code{exp_type = 'multi'} the following variables should also be set

\itemize{
  \item{\code{estimand}} choose from \code{ATE}, \code{ATT}, \code{ATC}, \code{ATZ}.
  \item{\code{which_z}} defines the referent variable for estimands: \code{ATT}, \code{ATC}, \code{ATZ}
}

For \code{optim_method = 'L-BFGS-B'} the following variables may also be set

\itemize{
  \item{\code{bal_tol}} tolerance for the optimization routine
  \item{\code{opt_constraints}} parameter constraints for the optimization.  Recommend \code{opt_constraints = c(-100,100)}.  If the algorithm sets some of the parameter values to the boundary, try to relax these values further.
}}

\item{suppress_warnings}{logical argument to suppress warnings of function checking \code{eb_pars}}
}
\value{
Object that contains the weights obtained from the balancing procedure and parameters from the optimization procedure

The object that is returned is a list that contains the following entries
\itemize{
\item{ \code{wts} - Optimal weights for the estimand of interest and matched number of moments.}
\item{ \code{convergence} - Convergence code from \code{optim} package}
\item{ \code{message} - Convergence message from \code{optim} package}
\item{ \code{opt_obj} - optimization object from \code{optim} package}
\item{ \code{eb_pars} - variable list the model was fit to}
\item{ \code{X} - Model frame created from data.}
\item{ \code{TA} - Vector of Exposure Assignments Used.}
}
}
\description{
Compute optimal balancing weights via entropy balancing
}
\examples{
# Binary exposure example

library(entbal)
n_obs <- 10000
X1 <- rnorm(n_obs)
X2 <- rnorm(n_obs)
A <- rbinom(n_obs, 1, plogis(X1))
A <- factor(ifelse(A == 1, 'A', 'B'))
D <- data.frame(X1, X2, A)
par_list <-  list('exp_type' = 'binary',
                  'n_moments' = 3,
                  'max_iters' = 1000,
                  'estimand' = 'ATC',
                  'verbose' = FALSE,
                  'optim_method' = 'l-bfgs-b',
                  'bal_tol' = 1e-8,
                  'opt_constraints' = c(-1000,1000),
                  'which_z' = 'A')
Q <- entbal(A ~ X1 + X2,
            data = D,
            eb_pars = par_list,
            suppress_warnings = FALSE)
out1 <- summary(Q)

# ---------------------------------------------------------------------------
# Multi-valued exposure - ATZ Example

C <- sample(1:3, n_obs, replace = TRUE)
X1 <- NA
X1[C == 1] <- rnorm(sum(C==1), mean = -0.5, sd = 3)
X1[C == 2] <- rnorm(sum(C==2), mean = 0, sd = 3)
X1[C == 3] <- rnorm(sum(C==3), mean = 0.5)
X2 <- rnorm(n_obs)
D <- data.frame(C, X1, X2)

par_list <-  list('exp_type' = 'multi',
                  'n_moments' = 3,
                  'max_iters' = 1000,
                  'estimand' = 'ATZ',
                  'verbose' = FALSE,
                  'optim_method' = 'l-bfgs-b',
                  'bal_tol' = 1e-8,
                  'opt_constraints' = c(-1000,1000),
                  'which_z' = 3)


P <- entbal(C ~ X1 + X2, data = D, eb_pars = par_list, suppress_warnings = FALSE)
out2 <- summary(P)

# ---------------------------------------------------------------------------
# Continuous exposure example

X1 <- rnorm(n_obs)
X2 <- rnorm(n_obs)
G <- rnorm(n_obs, mean = X1 - X2)
D <- data.frame(G, X1, X2)

par_list <-  list('exp_type' = 'continuous',
                  'n_moments' = 3,
                  'max_iters' = 1000,
                  'estimand' = 'ATE',
                  'verbose' = TRUE,
                  'optim_method' = 'l-bfgs-b',
                  'bal_tol' = 1e-8,
                  'opt_constraints' = c(-1000,1000))


O <- entbal(G ~ X1 + X2, data = D, eb_pars = par_list, suppress_warnings = FALSE)
out3 <- summary(O, show_parameters = TRUE)

}
