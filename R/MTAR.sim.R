#' MTAR- Simulation
#'
#' This function allows the user to simulate the data from a matrix-variate threshold autoregressive model with \emph{m} regimes.
#'
#' @name MTAR.sim
#' @rdname MTAR.sim
#' @export
#' @import matrixcalc vec
#' @import matrixNormal rmatnorm
#' @import stats rnorm
#' @param m number of rows of the matrix-variate time series
#' @param n number of columns of the matrix-variate time series
#' @param p number of lags
#' @param regimes number of regimes
#' @param Nsim number of simulated data
#' @param burnin number of data to be burned
#' @param constant \code{TRUE} or \code{FALSE} to include or not a constant matrix.
#' @param A \code{list} of values for the row-wise coefficient \code{(m×m)} matrices, one list for each lag containing a list of matrices for each regime.
#' @param B \code{list} of values for the column-wise coefficient \code{(n×n)} matrices, one list for each lag containing a list of matrices for each regime.
#' @param M \code{list} of values for the constant \code{(m×n)} matrices, a list of matrices for each regime.
#' @param threshold vector of thresholds used for the simulation equal to the number of regimes minus 1
#' @param st_type character string, specifying the type of the threshold used for the simulation. \describe{
#'  \item{\code{"trend",}}{\code{st = t/T} as a transition variable}
#'  \item{\code{"autoregressive",}}{an exogenous variable with an autoregressive dynamics}
#'}
#' @param Etype character string, specifying the type of the threshold used for the simulation. \describe{
#'  \item{\code{"identity",}}{two identity matrices are used as row-wise and column-wise covariance matrices}
#'  \item{\code{"normal",}}{the \code{mnxmn} covariance matrix is sampled from a normal distribution}
#'}
#' @param sdEt standard deviation for the sampling of the covariance matrix when \code{Etype} is equal to "normal" 
#' @return return a list containing the following:\describe{
#' \item{\code{data}}{an \code{m×n×T} array with data simulated from a MTAR model}
#' \item{\code{st}}{the transition variable used to simulate the data}
#' \item{\code{threshold}}{a vector of thresholds equal to the number of regimes minus 1}
#' \item{\code{datavec}}{an \code{mn×T} matrix with matrix data stacked in a vector}
#' }
#' @references Liu X. and Chen E.Y. (2022), Identification and estimation of threshold matrix-variate factor models. \emph{Scandinavian Journal of Statistics}. 49: 1383-1417
#' @author Andrea Bucci
#' @export
#' @keywords MTAR.sim
#' @examples
#'m = 3
#'n = 2
#'nsim = 1000
#'A = list(list(matrix(0.05, m, m), matrix(0.05, m, m), matrix(0.05, m, m)))
#'B = list(list(matrix(0.05, n, n), matrix(0.05, n, n), matrix(0.05, n, n)))
#'M = list(matrix(1, m, n), matrix(1.5, m, n), matrix(2, m, n))
#'simuldata = MTAR.sim(m=m, n=n, p = 1,Nsim = nsim, 
#' regimes = 3, threshold = c(0.3, 0.7), st_type = 'trend', 
#' Etype = 'normal', A = A, B= B, M = M, sderror = 0.1)
#' 
#' plot.ts(simuldata$data[1,1,])

MTAR.sim <- function(m = 2, n = 3, p = 1, regimes = 2, Nsim = 1000, burnin = 200, 
                     constant = TRUE, A, B, M = NULL, threshold = 0.3,
                     st_type = c('trend', 'autoregressive'), Etype=c('identity','normal'),
                     sdEt = 0.2){
  if(regimes < 2)
    stop('The number of regimes should be greater than one.')
  if(length(threshold) != (regimes-1))
    stop('The number of thresholds should be equal to the number of transitions.')
  nsim = Nsim + burnin
  if(st_type == 'trend'){
    st = rep(0, (nsim+1))
    st = (1:nsim)/nsim
    index1 = rep(0, nsim)
    index2 = rep(0, nsim)
  }else if(st_type == 'autoregressive'){
    st = rep(0, (nsim+1))
    epsilonst = rnorm((nsim+1))
    st[1] = epsilonst[1]
    for(t in 2:nsim){
      st[t] = 0.95 * st[(t-1)] + epsilonst[t]
    }
    st = st[1:(length(st)-1)]
  }
  index = list()
  for(i in 1:regimes){
    if(i == 1){
      index[[i]] =  rep(1, nsim)
      index[[i]][st >= threshold[i]] = 0  
    }else if(i > 1 & i != regimes){
      index[[i]] =  rep(0, nsim)
      index[[i]][st >= threshold[(i-1)] & st < threshold[i]] = 1
    }else{
      index[[i]] =  rep(1, nsim)
      index[[i]][st < threshold[(i-1)]] = 0  
    }
  }
  Xt = array(0, dim = c(m,n,(nsim+1)))
  Xt[,,1:p] = matrix(rnorm(m*n,mean=0,sd=1), m, n)
  for(t in (1+p):nsim){
    if(Etype == 'identity'){
      Et = matrixNormal::rmatnorm(M = matrix(0, m, n), U=sderror*diag(m), V=sderror*diag(n))
    }else if(Etype == 'normal'){
      Et = matrix(rnorm(m*n,mean=0,sd=sdEt), m, n)
    }
    Xreg = matrix(0, m, n)
    for(i in 1:regimes){
      Xval = matrix(0, m, n)
      for(k in 1:p){
        Xval = A[[k]][[i]]%*%Xt[,,(t-k)]%*%t(B[[k]][[i]]) + Xval
      }
      if(constant == TRUE){
        Xreg = (M[[i]] + Xval + Et)*index[[i]][t] + Xreg
      }else{
        Xreg = (Xval + Et)*index[[i]][t] + Xreg
      }
    }
    Xt[,,t] = Xreg
  }
  Xt = Xt[,,(burnin+1):nsim]
  xtvec = matrix(nrow = nsim, ncol = m*n)
  for(l in 1:Nsim){
    xtvec[l,] <- matrixcalc::vec(Xt[,,l])
  }
  simuldata = list(data = Xt, st = st[(burnin+1):nsim], threshold  = threshold, datavec = xtvec)
  return(simuldata)
}
