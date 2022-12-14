% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MTAR.R
\name{MTAR}
\alias{MTAR}
\title{MTAR- Estimation}
\usage{
MTAR(
  data,
  p = 1,
  regimes = 3,
  maxiter = 200,
  st,
  q = 0.1,
  constant = TRUE,
  initA,
  initB,
  initM = NULL,
  epsilon = 10^(-3),
  verbose = TRUE,
  ncores = 6
)
}
\arguments{
\item{data}{\code{array} of dependent variables of dimension \code{(m×n×T)}}

\item{p}{lag order}

\item{regimes}{number of regimes}

\item{maxiter}{Number of maximum iterations of the iterated least squares.}

\item{st}{single transition variable for all the equation of dimension \code{(T×1)}}

\item{q}{quantile for the computation of the candidate values as thresholds.}

\item{constant}{\code{TRUE} or \code{FALSE} to include or not a constant matrix.}

\item{initA}{\code{list} of initial values for the row-wise coefficient \code{(m×m)} matrices, one list for each lag containing a list of matrices for each regime.}

\item{initB}{\code{list} of initial values for the column-wise coefficient \code{(n×n)} matrices, one list for each lag containing a list of matrices for each regime.}

\item{initM}{\code{list} of initial values for the constant \code{(m×n)} matrices, a list of matrices for each regime.}

\item{epsilon}{convergence check measure}

\item{verbose}{\code{TRUE} or \code{FALSE} to show all the iteration outputs.}

\item{ncores}{Number of cores used for parallel computation. Set to \code{NULL} by default and automatically calculated.}
}
\value{
return a list containing the following:\describe{
\item{\code{A}}{a list of estimated row-wise coefficient matrices}
\item{\code{B}}{a list of estimated column-wise coefficient matrices}
\item{\code{M}}{(optional if constant = TRUE) a list of estimated constant matrices}
\item{\code{threshold}}{a vector of estimated thresholds equal to the number of regimes minus 1}
\item{\code{Sigma}}{an \code{mnxmn} covariance matrix of the error term}
\item{\code{stderr}}{a list of standard errors for the estimated coefficients}
\item{\code{fitted}}{an array of fitted values}
\item{\code{residuals}}{an array of residuals}
\item{\code{dimensions}}{a list containing the dimension of the data in each regime}
}
}
\description{
This function allows the user to estimate the coefficients of a matrix-variate threshold autoregressive model with \emph{m} regimes through iterated least squares.
}
\examples{
\donttest{
m = 3
n = 2
nsim = 1000
A = list(list(matrix(0.05, m, m), matrix(0.10, m, m), matrix(0.15, m, m)))
B = list(list(matrix(0.05, n, n), matrix(0.10, n, n), matrix(0.15, n, n)))
M = list(matrix(1, m, n), matrix(1.5, m, n), matrix(2, m, n))
simuldata = MTAR.sim(m=m, n=n, p = 1, nsim = 1000, 
regimes = 3, threshold = c(0.3, 0.7), st_type = 'trend', 
Etype = 'normal', A = A, B= B, M = M, sdEt = 0.1)

plot.ts(simuldata$data[1,1,])

initA = list(list(diag(rep(0.05,m)), diag(rep(0.05,m)), diag(rep(0.05,m))))
initB = list(list(diag(rep(0.05,n)), diag(rep(0.05,n)), diag(rep(0.05,n))))
initM = list(matrix(1, m, n), matrix(1, m, n), matrix(1, m, n))

fit.MTAR <- MTAR(simuldata$data, regimes = 3, maxiter = 30, st = simuldata$st,
initA = initA, initB = initB, initM = initM, q = 0.1,
verbose = F, ncores = 4)
# a few methods for VLSTAR
print(fit.MTAR)
summary(fit.MTAR)
plot(fit.MTAR)
predict(fit.MTAR, st.new = 1, n.ahead = 1)
coef(fit.MTAR)}

}
\references{
Liu X. and Chen E.Y. (2022), Identification and estimation of threshold matrix-variate factor models. \emph{Scandinavian Journal of Statistics}. 49: 1383-1417
}
\author{
Andrea Bucci
}
\keyword{MTAR}
