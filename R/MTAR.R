#' MTAR- Estimation
#'
#' This function allows the user to estimate the coefficients of a matrix-variate threshold autoregressive model with \emph{m} regimes through iterated least squares.
#'
#'
#' @param data \code{array} of dependent variables of dimension \code{(m×n×T)}
#' @param p lag order
#' @param regimes number of regimes
#' @param maxiter Number of maximum iterations of the iterated least squares.
#' @param st single transition variable for all the equation of dimension \code{(T×1)}
#' @param q quantile for the computation of the candidate values as thresholds.
#' @param constant \code{TRUE} or \code{FALSE} to include or not a constant matrix.
#' @param initA \code{list} of initial values for the row-wise coefficient \code{(m×m)} matrices, one list for each lag containing a list of matrices for each regime.
#' @param initB \code{list} of initial values for the column-wise coefficient \code{(n×n)} matrices, one list for each lag containing a list of matrices for each regime.
#' @param initM \code{list} of initial values for the constant \code{(m×n)} matrices, a list of matrices for each regime.
#' @param epsilon convergence check measure
#' @param verbose \code{TRUE} or \code{FALSE} to show all the iteration outputs.
#' @param ncores Number of cores used for parallel computation. Set to \code{NULL} by default and automatically calculated.
#' @return return a list containing the following:\describe{
#' \item{\code{A}}{a list of estimated row-wise coefficient matrices}
#' \item{\code{B}}{a list of estimated column-wise coefficient matrices}
#' \item{\code{M}}{(optional if constant = TRUE) a list of estimated constant matrices}
#' \item{\code{threshold}}{a vector of estimated thresholds equal to the number of regimes minus 1}
#' \item{\code{Sigma}}{an \code{mnxmn} covariance matrix of the error term}
#' \item{\code{stderr}}{a list of standard errors for the estimated coefficients}
#' \item{\code{dimensions}}{a list containing the dimension of the data in each regime}
#' }
#' @references Liu X. and Chen E.Y. (2022), Identification and estimation of threshold matrix-variate factor models. \emph{Scandinavian Journal of Statistics}. 49: 1383-1417
#' @author Andrea Bucci
#' @export
#' @keywords MTAR
#' @examples
#' \donttest{
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
#' 
#' initA = list(list(diag(rep(0.05,m)), diag(rep(0.05,m)), diag(rep(0.05,m))))
#' initB = list(list(diag(rep(0.05,n)), diag(rep(0.05,n)), diag(rep(0.05,n))))
#' initM = list(matrix(1, m, n), matrix(1, m, n), matrix(1, m, n))
#' 
#'fit.MTAR <- MTAR(simuldata$data, regimes = 3, maxiter = 30, st = simuldata$st,
#'initA = initA, initB = initB, initM = initM, q = 0.1,
#'verbose = F, ncores = 4)
#'# a few methods for VLSTAR
#'print(fit.MTAR)
#'summary(fit.MTAR)
#'plot(fit.MTAR)
#'predict(fit.MTAR, st.new = 1, n.ahead = 1)
#'coef(fit.MTAR)}
#'
MTAR <- function(data, p = 1, regimes = 3, maxiter = 200, st, q = 0.10,
                 constant = TRUE, initA, initB, initM = NULL,  
                 epsilon = 10^(-3), verbose = TRUE, ncores = 6){
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if(is.null(ncores)){
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      ncores <- 2L
    } else {
      ncores <- parallel::detectCores()
    }}
  m = dim(data)[1]
  n = dim(data)[2]
  Tlength = dim(data)[3]
  st = st
  st_q = quantile(st, probs = seq(0.1,0.9,q))
  message(paste('Alternative least square estimation of a MTAR\n'))
  stq = expand.grid(rep(list(st_q), (regimes-1)))
  stq1 = as.matrix(stq[rowSums(as.matrix(stq[,-ncol(stq)]) < stq[,-1]) == ncol(stq) - 1,])
  colnames(stq1) = paste('c', rep(1:(regimes-1)), sep = '')
  Alist = list()
  Blist = list()
  Mlist = list()
  datalist = list()
  errors = list()
  Treg = rep(0, regimes)
  QQ = rep(0, nrow(stq1))
  pb <- txtProgressBar(0, nrow(stq1), style = 3)
  for(i in 1:nrow(stq1)){
    setTxtProgressBar(pb, i)
    A = initA
    B = initB
    M = initM
    for(j in 1:regimes){
      if(j == 1){
        datalist[[j]] = data[,,which(st < stq1[i,j])] 
      }else if(j > 1 & j != regimes){
        datalist[[j]] = data[,,which(st >= stq1[i,(j-1)] & st < stq1[i, j])]
      }else{
        datalist[[j]] = data[,,which(st >= stq1[i,(j-1)])]
      }
      Treg[j] = dim(datalist[[j]])[3]
      iter = 0
      bdiff = matrix(1, nrow = n, ncol = n)
      while(iter < maxiter & all(bdiff > epsilon)){
        iter <- iter+1
        temp1A = matrix(0, nrow = m, ncol = m) 
        temp2A = matrix(0, nrow = m, ncol = m)
        for(k in 1:p){
          seqex = 1:p
          seqex = seqex[-k]
        for(t in (p+1):Treg[j]){
          temp3A = matrix(0, nrow = m, ncol = m)
          for(l in seqex){
            temp3A = A[[l]][[j]]%*%datalist[[j]][,,(t-l)]%*%t(B[[l]][[j]])%*%B[[k]][[j]]%*%t(data[,,(t-k)])+temp3A
          }
            if(constant == TRUE){
            temp1A = datalist[[j]][,,t] %*%B[[k]][[j]]%*%t(datalist[[j]][,,(t-k)]) - 
              M[[j]]%*%B[[k]][[j]]%*%t(datalist[[j]][,,(t-k)]) - temp3A + temp1A
          }else{
              temp1A = (datalist[[j]][,,t] %*%B[[k]][[j]]%*%t(datalist[[j]][,,(t-k)])-temp3A) + temp1A
          }
            temp2A = (datalist[[j]][,,(t-k)]%*%t(B[[k]][[j]])%*%B[[k]][[j]]%*%t(datalist[[j]][,,(t-k)])) + temp2A
          }
          A[[k]][[j]] = frob.rescale(temp1A %*% MASS::ginv(temp2A))
        }
        temp1B = matrix(0, nrow = n, ncol = n)
        temp2B = matrix(0, nrow = n, ncol = n)
        for(k in 1:p){
          seqex = 1:p
          seqex = seqex[-k]
        for(t in (p+1):Treg[j]){
          temp3B = matrix(0, nrow = n, ncol = n)
            for(l in seqex){
              temp3B = B[[l]][[j]]%*%t(datalist[[j]][,,(t-l)])%*%t(A[[l]][[j]])%*%A[[k]][[j]]%*%data[,,(t-k)]+temp3B
            }
          if(constant == TRUE){
            temp1B = t(datalist[[j]][,,t]-M[[j]]) %*%A[[k]][[j]]%*%datalist[[j]][,,(t-k)] - temp3B + temp1B
          }else{
            temp1B = t(datalist[[j]][,,t]) %*%A[[k]][[j]]%*%datalist[[j]][,,(t-k)] -temp3B + temp1B
          }
          temp2B = t(datalist[[j]][,,(t-k)])%*%t(A[[k]][[j]])%*%A[[k]][[j]]%*%datalist[[j]][,,(t-k)] + temp2B
          Bcheck = temp1B %*% MASS::ginv(temp2B)
          if(iter >1){
            bdiff = abs(B[[k]][[j]] - Bcheck) 
          }
          }
          B[[k]][[j]] = temp1B %*% MASS::ginv(temp2B)
        }
        if(constant == TRUE){
          temp1M = matrix(0, nrow = m, ncol = n)
          for(t in (p+1):Treg[j]){
            temp2M = matrix(0, nrow = m, ncol = n)
            for(k in 1:p){
              temp2M = A[[k]][[j]]%*%datalist[[j]][,,(t-k)]%*%t(B[[k]][[j]])
            }
            temp1M = datalist[[j]][,,t] - temp2M + temp1M
          }
          M[[j]] = temp1M/(Treg[j]-1)  
        }
        if(verbose==TRUE){
          print(iter)}
        }
      ####Compute sum of squared residuals####
      epsilon1 = matrix(0, nrow = Treg[j], ncol = m*n)
      epsilon1m = matrix(0, ncol = m*n, nrow = m*n)
      for(l in (p+1):Treg[j]){
        epsilon2 = matrix(0, nrow = m, ncol = n)
        for(k in 1:p){
          epsilon2 = A[[k]][[j]]%*%datalist[[j]][,,(l-k)]%*%t(B[[k]][[j]])+epsilon2
        }
        if(constant == TRUE){
          epsilon1[l,] = matrixcalc::vec(datalist[[j]][,,l] - M[[j]]- epsilon2)
        }else{
          epsilon1[l,] = matrixcalc::vec(datalist[[j]][,,l] - epsilon2)
        }
        epsilon1m = epsilon1[l,] %*% t(epsilon1[l,]) + epsilon1m
      }
      errors[[j]] = epsilon1m
    }
    Alist[[i]] = A
    Blist[[i]] = B
    if(constant == TRUE){
      Mlist[[i]] = M
    }
    errlist = lapply(errors, tr)
    QQ[i] = Reduce("+", errlist)
    if(verbose==TRUE){
      cat(paste(i, 'th', ' combination', sep = ''))}
    Sys.sleep(time = 1)
  }
  close(pb)
  if(constant == TRUE){
    theta <- list(A = Alist[[which.min(QQ)]], B = Blist[[which.min(QQ)]], 
                  M = Mlist[[which.min(QQ)]],
                  c = t(as.matrix(stq1[which.min(QQ),])))
  }else{
    theta <- list(A = Alist[[which.min(QQ)]], B = Blist[[which.min(QQ)]], 
                  c = t(as.matrix(stq1[which.min(QQ),])))
  }
  stationary_A <- rapply(theta$A, sparsevar::spectralRadius)
  stationary_B <- rapply(theta$B, sparsevar::spectralRadius)
  stat_check = (Reduce("*", stationary_A)*Reduce("*", stationary_B))<1
  index = list()
  for(i in 1:regimes){
    if(i == 1){
      index[[i]] =  which(st < theta$c[1,i])
    }else if(i > 1 & i != regimes){
      index[[i]] =  which(st >= theta$c[1,(i-1)] & st < theta$c[1,i])
    }else{
      index[[i]] =  which(st >= theta$c[1,(i-1)])
    }
  }
  fit.val = array(data = NA, dim = c(m,n,Tlength))
  res.val = array(data = NA, dim = c(m,n,Tlength))
  xvec = matrix(NA, nrow = Tlength, ncol = m*n)
  Sigma = matrix(0, nrow = m*n, ncol = m*n)
  for(i in (p+1):Tlength){
    for(j in 1:regimes){
      if(i %in% index[[j]]){
        fittmp = matrix(0, nrow = m, ncol = n)
        for(k in 1:p){
          fittmp = theta$A[[k]][[j]]%*%data[,,(i-k)]%*%t(theta$B[[k]][[j]]) 
        }
        if(constant == TRUE){
          fit.val[,,i] = theta$M[[j]] + fittmp
        }else{
          fit.val[,,i] = fittmp  
        }
      }
    }
    res.val[,,i] = data[,,i] - fit.val[,,i]
    Sigma = Sigma + matrixcalc::vec(res.val[,,i]) %*% t(matrixcalc::vec(res.val[,,i]))
    xvec[i,] = ks::vec(data[,,(i-1)])
  }
  Omega = Sigma/(Tlength-1L)
  Sigmar = list()
  Sigmar2 = list()
  for(j in 1:regimes){
    dfperm = aperm(data[,,index[[j]]], c(3,1,2))
    for(k in 1:p){
      Sigmar2[[k]] = covtosd(Sigma, c(m,n), R = 1, p = 1, A = list(theta$A[[k]][[j]], theta$B[[k]][[j]]), 
                                 AX = dfperm[-1,,])
    }
    Sigmar[[j]] = Sigmar2
  }
  # for(s in 1:pred.step){
  #   tensor[,,(T1+s)] = iterations[[maxiter]]$A %*% tensor[,,(T+s-1)] %*% t(iterations[[maxiter]]$B)
  # }
  # predictions = tensor[,,((T1+1):(T1+pred.step))]
  if(constant == TRUE){
    results = list(stationary = stat_check, A = theta$A, 
                   B = theta$B, M = theta$M, threshold = theta$c, Sigma = Omega,
                   stderr = Sigmar, dimensions = lapply(index, length)) 
  }else{
    results = list(stationary = stat_check, A = theta$A, 
                   B = theta$B,  threshold = theta$c, Sigma = Omega,
                   stderr = Sigmar, dimensions = lapply(index, length)) 
  }
  return(results)
}


