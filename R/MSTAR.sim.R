#' MsTAR- Simulation
#'
#' This function allows the user to simulate the data from a matrix-variate smooth transition autoregressive model with \emph{m} regimes.
#'
#' @name MSTAR.sim
#' @rdname MSTAR.sim
#' @export
#' @param m number of rows of the matrix-variate time series
#' @param n number of columns of the matrix-variate time series
#' @param p number of lags
#' @param regimes number of regimes
#' @param nsim number of simulated data
#' @param constant \code{TRUE} or \code{FALSE} to include or not a constant matrix.
#' @param A \code{list} of values for the row-wise coefficient \code{(m×m)} matrices, one list for each lag containing a list of matrices for each regime.
#' @param B \code{list} of values for the column-wise coefficient \code{(n×n)} matrices, one list for each lag containing a list of matrices for each regime.
#' @param M \code{list} of values for the constant \code{(m×n)} matrices, a list of matrices for each regime.
#' @param threshold vector of thresholds used for the simulation equal to the number of regimes minus 1
#' @param gamma vector of gamma used for the simulation equal to the number of regimes minus 1
#' @param st_type character string, specifying the type of the threshold used for the simulation. \describe{
#'  \item{\code{"trend",}}{\code{st = t/T} as a transition variable}
#'  \item{\code{"autoregressive",}}{an exogenous variable with an autoregressive dynamics}
#'}
#' @param Etype character string, specifying the type of the threshold used for the simulation. \describe{
#'  \item{\code{"identity",}}{two identity matrices are used as row-wise and column-wise covariance matrices}
#'  \item{\code{"normal",}}{the \code{mnxmn} covariance matrix is sampled from a normal distribution}
#'}
#' @param sdEt standard deviation for the sampling of the covariance matrix
#' @return return a list containing the following:\describe{
#' \item{\code{data}}{an \code{m×n×T} array with data simulated from a MTAR model}
#' \item{\code{st}}{the transition variable used to simulate the data}
#' \item{\code{threshold}}{a vector of thresholds equal to the number of regimes minus 1}
#' \item{\code{gamma}}{a vector of slope parameters equal to the number of regimes minus 1}
#' \item{\code{datavec}}{an \code{mn×T} matrix with matrix data stacked in a vector}
#' \item{\code{g}}{a \code{T×1} vector of transition functions}
#' }
#' @references Bucci A. (2022), A smooth transition autoregressive model for matrix-variate time series. \emph{arXiV}.
#' @author Andrea Bucci
#' @export
#' @keywords MSTAR.sim
#' @examples
#'m = 3
#'n = 2
#'nsim = 1000
#'A = list(list(matrix(0.05, m, m), matrix(0.10, m, m), matrix(0.15, m, m)))
#'B = list(list(matrix(0.05, n, n), matrix(0.10, n, n), matrix(0.15, n, n)))
#'M = list(matrix(1, m, n), matrix(1.5, m, n), matrix(2, m, n))
#'simuldata = MSTAR.sim(m=m, n=n, p = 1, nsim = 1000, 
#' regimes = 3, threshold = c(0.3, 0.7), gamma = c(10,15), st_type = 'trend', 
#' Etype = 'normal', A = A, B= B, M = M, sdEt = 0.1)
#' 
#' plot.ts(simuldata$data[1,1,])


MSTAR.sim <- function(m = 2, n = 3, p = 1, regimes = 3, nsim = 1000, constant = TRUE,
                      A, B, M = NULL, threshold = 2, gamma = 2,
                      st_type = c('trend', 'autoregressive'), Etype=c('identity','normal'), sdEt = 0.1){
  if(regimes < 2)
    stop('The number of regimes should be greater than one.')
  if(length(threshold) != (regimes-1))
    stop('The number of thresholds should be equal to the number of transitions.')
  if(length(gamma) != (regimes-1))
    stop('The number of thresholds should be equal to the number of transitions.')
  nsim = nsim+1
  c = threshold
  gamma = gamma
  if(st_type == 'trend'){
    st = rep(0, (nsim+1))
    gt = rep(0, (nsim))
    st = (1:nsim)/nsim
  }else if(st_type == 'autoregressive'){
    st = rep(0, (nsim+1))
    epsilonst = rnorm((nsim+1))
    gt = rep(0, (nsim+1))
    st[1] = epsilonst[1]
    for(t in 2:nsim){
      st[t] = 0.95 * st[(t-1)] + epsilonst[t]
      gt[t] = 1/(1 + exp(-gamma*(st[(t-1)]-c)))
    }
    st = st[1:(length(st)-1)]
    gt = gt[-1] 
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
      index[[i]][st <= threshold[(i-1)]] = 0  
    }
  }
  Xt = array(0, dim = c(m,n,(nsim+p)))
  Xt[,,p] = M[[1]]
  g <- rep(0, nsim)
  for(t in (p+1):nsim){
    if(Etype == 'identity'){
      Et = matrixNormal::rmatnorm(M = matrix(0, m, n), U=sdEt*diag(m), V=sdEt*diag(n))
    }else if(Etype == 'normal'){
      Et = matrix(rnorm(m*n,mean=0,sd=sdEt), m, n)
    }
    Xreg = matrix(0, m, n)
    for(i in 1:regimes){
      if(i == 1){
        g[t] = (1/(1 + exp(-gamma[i]*(st[(t-1)]-c[i]))))*index[[i]][t]
      }else{
        g[t] = (1/(1 + exp(-gamma[(i-1)]*(st[(t-1)]-c[(i-1)]))))*index[[i]][t]
      }
      Xval = matrix(0, m, n)
      for(k in 1:p){
        Xval = A[[k]][[i]]%*%Xt[,,(t-k)]%*%t(B[[k]][[i]]) + Xval
      }
      if(constant == TRUE){
        Xreg = (M[[i]] + Xval)*g[t]*index[[i]][t] + Xreg
      }
    }
    Xt[,,t] = Xreg + Et
  }
  Xt = Xt[,,(p+1):nsim]
  xtvec = matrix(nrow = nsim, ncol = m*n)
  for(l in 1:(nsim-p-1)){
    xtvec[l,] <- matrixcalc::vec(Xt[,,l])
  }
  simuldata = list(data = Xt, st = st[(p+1):nsim], gamma = gamma, threshold  = c, datavec = xtvec, g = g)
  return(simuldata)
}

