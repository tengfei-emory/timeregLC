multinfopi <- function(X,mu,V,prop){
  output <- list()
  m = dim(mu)[1]; p = dim(mu)[2]
  n           = dim(X)[2]
  aum         = m*(m+3)/2
  Qpi2        = matrix(0,p-1,p-1)
  Qpitheta    = matrix(0,p-1,aum*p)
  Qtheta2     = matrix(0,aum*p,aum*p)

  indInfo     = list()

  for (t in 1:n){
    indQpi2        = matrix(0,p-1,p-1)
    indQpitheta    = matrix(0,p-1,aum*p)
    indQtheta2     = matrix(0,aum*p,aum*p)
    param     = intermedpi(X[,t],mu,V,prop)
    atbar     = param$at %*% param$alphat
    st        = (param$at - atbar %x% matrix(1,1,p)) %*% diag(as.vector(param$alphat))
    ctt       = t(param$ct)
    sct       = param$ct %*% diag(as.vector(param$alphat))
    for (i in 1:p){
      Qpitheta[,(aum*(i-1)+1):(aum*i)] = Qpitheta[,(aum*(i-1)+1):(aum*i)] + st[,i] %*% t(ctt[i,])
      indQpitheta[,(aum*(i-1)+1):(aum*i)] = indQpitheta[,(aum*(i-1)+1):(aum*i)] + st[,i] %*% t(ctt[i,])
      Qtheta2[(aum*(i-1)+1):(aum*i),(aum*(i-1)+1):(aum*i)] = Qtheta2[(aum*(i-1)+1):(aum*i),(aum*(i-1)+1):(aum*i)] - param$alphat[i]*(param$Ct[,,i]-((1-param$alphat[i])*param$ct[,i] %*% t(ctt[i,])))
      indQtheta2[(aum*(i-1)+1):(aum*i),(aum*(i-1)+1):(aum*i)] = indQtheta2[(aum*(i-1)+1):(aum*i),(aum*(i-1)+1):(aum*i)] - param$alphat[i]*(param$Ct[,,i]-((1-param$alphat[i])*param$ct[,i] %*% t(ctt[i,])))
      for (j in 1:p){
        if (i != j){
          Qtheta2[(aum*(i-1)+1):(aum*i),(aum*(j-1)+1):(aum*j)] = Qtheta2[(aum*(i-1)+1):(aum*i),(aum*(j-1)+1):(aum*j)] - sct[,i] %*% t(sct[,j])
          indQtheta2[(aum*(i-1)+1):(aum*i),(aum*(j-1)+1):(aum*j)] = indQtheta2[(aum*(i-1)+1):(aum*i),(aum*(j-1)+1):(aum*j)] - sct[,i] %*% t(sct[,j])
        }
      }
    }
    Qpi2 = Qpi2 - atbar %*% t(atbar)
    indQpi2 = indQpi2 - atbar %*% t(atbar)
 
    indInfo[[t]] = - rbind(cbind(indQpi2,indQpitheta) , cbind(t(indQpitheta),indQtheta2))
  }
  
  output$infopi = - rbind(cbind(Qpi2,Qpitheta) , cbind(t(Qpitheta),Qtheta2))
  output$indInfo = indInfo
  
  return(output)
}