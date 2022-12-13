MTAR.sim <- function(m = 2, n = 3, p = 1, regimes = 2, Nsim = 1000, burnin = 200, 
                     A, B, M, threshold = 0.3,
                     st_type = c('trend', 'autoregressive'), Etype=c('identity','normal'),
                     sderror = 0.2){
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
      Et = matrix(rnorm(m*n,mean=0,sd=sderror), m, n)
    }
    Xreg = matrix(0, m, n)
    for(i in 1:regimes){
      Xval = matrix(0, m, n)
      for(k in 1:p){
        Xval = A[[k]][[i]]%*%Xt[,,(t-k)]%*%t(B[[k]][[i]]) + Xval
      }
      Xreg = (M[[i]] + Xval + Et)*index[[i]][t] + Xreg
    }
    Xt[,,t] = Xreg
  }
  Xt = Xt[,,(burnin+1):nsim]
  xtvec = matrix(nrow = nsim, ncol = m*n)
  for(l in 1:Nsim){
    xtvec[l,] <- vec(Xt[,,l])
  }
  simuldata = list(data = Xt, st = st[(burnin+1):nsim], threshold  = threshold, datavec = xtvec)
  return(simuldata)
}
