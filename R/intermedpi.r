intermedpi <- function(xt,mu,V,prop){
  require(MASS)
  require(mvtnorm)
  require(matrixcalc)
  param <- list()
  m = dim(mu)[1]; p = dim(mu)[2]
  gt = matrix(0,p,1)
  for (j in 1:p){
    gt[j] = prop[j]*mvtnorm::dmvnorm(as.vector(xt),mean=as.vector(mu[,j]),sigma=V[,,j])
  }
  param$gt     = gt
  param$alphat = gt/sum(gt)
  param$at     = cbind(diag(rep(1,p-1)/prop[1:p-1]),-rep(1,p-1)/prop[p])
  param$bt     = matrix(0,m,p)
  param$Bt     = array(0,dim=c(m,m,p))
  aum          = m*(m+3)/2
  param$ct     = matrix(0,aum,p)
  param$D      = duplication.matrix(n=m)
  param$At     = array(0,dim=c(p-1,p-1,p))
  param$Ct     = array(0,dim=c(aum,aum,p))
  for (i in 1:p){
    inverV        = ginv(V[,,i])
    param$bt[,i]  = inverV %*% (xt-mu[,i])
    param$Bt[,,i] = inverV-param$bt[,i] %*% t(param$bt[,i])
    param$ct[,i]  = c(param$bt[,i],-0.5*t(param$D) %*% as.vector(param$Bt[,,i]))
    St            = (t(param$bt[,i]) %x% inverV) %*% param$D
    param$Ct[,,i] = rbind(cbind(inverV,St),cbind(t(St), 0.5*t(param$D)%*%((inverV-2*param$Bt[,,i]) %x% inverV)%*%param$D))
  }
  return(param)
}
