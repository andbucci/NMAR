####Trace of a matrix####
tr <- function(A) {
  n <- dim(A)[1] # get dimension of matrix
  tr <- 0 # initialize trace value
  
  # Loop over the diagonal elements of the supplied matrix and add the element to tr
  for (k in 1:n) {
    l <- A[k,k]
    tr <- tr + l
  }
  return(tr[[1]])
}

####Frobenius distance between two matrices####
Frob <- function(A, B){
  D = A-B
  d = sqrt(tr(t(D)%*%D))
  return(d)
}

####Sum of squared residuals for the MSTAR####
SSQ <- function(param, data){
  datalist = data
  Y <- datalist$Y
  m = dim(Y)[1]
  n = dim(Y)[2]
  gamma <- param[1]
  c <- param[2]
  dify <- matrix(0, ncol = m*n, nrow = m*n)
  for(z in 2:dim(Y)[3]){
    g = (1/(1 + exp(-gamma*(datalist$st[(z-1)]-c))))
    dify = matrixcalc::vec(datalist$Y[,,z] - datalist$A %*% datalist$Y[,,(z-1)]%*%t(datalist$B) - g*datalist$C%*%datalist$Y[,,(z-1)]%*%t(datalist$D))%*%
      t(matrixcalc::vec(datalist$Y[,,z] - datalist$A %*% datalist$Y[,,(z-1)]%*%t(datalist$B) - g*datalist$C%*%datalist$Y[,,(z-1)]%*%t(datalist$D)))+dify
  }
  tracemin = tr(dify)
  return(tracemin)
}


####From estimated covariance matrix to standard errors####
covtosd <- function(Sigma, dim, R, p, A, AX){
  pdim = prod(dim)
  fdim = dim**2
  ndim = sum(fdim)
  K <- length(dim)
  P <- length(R)
  Tlength = dim(AX)[1]
  Gamma <- matrix(0,sum(R)*ndim,sum(R)*ndim)
  n = 0
  for(i in c(1:p)){
    for(j in c(1:R[i])){
      for(k in c(1:(K-1))){
        r1 <- matrix(0, sum(R)*ndim,1)
        a1 <- as.vector(A[[k]])
        r1[(n*ndim + sum(fdim[0:(k-1)]) + 1): (n*ndim + sum(fdim[0:k])),] = a1
        Gamma = Gamma + r1 %*% t(r1)
        if (K==2){
          for (l in c(1:R[i])){
            if (l != j){
              r1 <- matrix(0, sum(r)*ndim, 1)
              a = as.vector(A[[k]])
              r1[(n*ndim + sum(fdim[0:(k-1)]) + 1): (n*ndim + sum(fdim[0:k])),] = a
              Gamma = Gamma + r1 %*% t(r1)
            }
          }
        }
      }
      n = n+1
    }
  }
  WT = c(); Q = list(); perm = list(); size = list()
  for (i in c(1:p)){
    for (j in c(1:R[i])){
      for (k in c(1:K)){
        if (length(Q) < K){
          perm[[k]] = c(k+1, (1:K)[-k] + 1, 1)
          size[[k]] = c(dim[k], prod(dim[-k]), Tlength)
          s = if (is.na(prod(dim[(k+1):K]))) 1 else prod(dim[(k+1):K])
          Q[[k]] = kronecker(diag(s), pm(dim[k], prod(dim[0:(k-1)])))
        }
        # AXX = AX[[i]][[j]][[k]] # AXX = rTensor::ttl(xx, A.new[[i]][[j]][-k], c(2:(K+1))[-k])
        AXX = abind::asub(AX, (2+p-i):(Tlength-i), 1, drop=FALSE)
        AXfold = array(aperm(AXX, perm[[k]]), size[[k]])
        AXI = apply(AXfold,3,function(x){kronecker(x, diag(dim[k])) %*% Q[[k]]})
        AXI.array <- array(AXI,c(dim[k]^2,pdim,Tlength))
        WT <- abind::abind(WT, AXI.array,along=1)
      }
    }
  }
  WT = aperm(WT, c(3,1,2))
  WSigma <- tensor::tensor(WT,Sigma,3,1) #t*(d1^2+d2^2+d^3)*(d1d2d3)
  EWSigmaWt <- tensor::tensor(WSigma,WT,c(3,1),c(3,1))/Tlength
  H <- tensor::tensor(WT,WT,c(3,1),c(3,1))/Tlength + Gamma #r(d1^2+d2^2+d^2)*r(d1^2+d2^2+d^3)
  Hinv <- solve(H)
  cov <- Hinv %*% EWSigmaWt %*% Hinv
  sd = list()
  for (p in c(1:P)){
    if (is.na(R[p])) stop("p != length(R)")
    if (R[p] == 0) next
    sd[[p]] <- lapply(1:R[p], function(j) {lapply(1:K, function(i) {list()})})
  }
  for (i in c(1:P)){
    for (j in c(1:R[i])){
      for (k in c(1:K)){
        left <- sum(dim^2)*sum(R[0:(i-1)]) + sum(dim^2)*(j-1) + sum((dim^2)[1:(k-1)])+1
        right <- sum(dim^2)*sum(R[0:(i-1)]) + sum(dim^2)*(j-1) + sum((dim^2)[1:k])
        sd[[i]][[j]][[k]] <- array(sqrt(diag(cov)[left:right]), c(dim[k], dim[k]))
      }
    }
  }
  return(sd)
}



# Permutation matrix pm
pm <- function(m,n){
  ## m: an array of dimensions of matrices \eqn{A_1,A_2,\cdots,A_k}
  ## n: length of time
  ## return: Permutation matrix pm
  mat <- matrix(0,m*n,m*n)
  for (i in c(1:n)){
    for (j in c(1:m)){
      mat <- mat + kronecker(em(n,m,i,j),t(em(n,m,i,j)))
    }
  }
  return(mat)
}

# Permutation matrix em
em <- function(m,n,i,j){
  ## m,n,i,j set \eqn{m \times n} zero matrix with \eqn{A_{ij} = 1}
  ## return: Permutation matrix em such that \eqn{A_{ij} = 1} and other entries equals 0.
  mat <- matrix(0,m,n)
  mat[i,j] <- 1
  return(mat)
}

####Frobenius rescale for the matrices in A####
frob.rescale <- function(A){
  m = A
  nrm = norm(m, 'f')
  Anorm = m/nrm
  return(Anorm)
}