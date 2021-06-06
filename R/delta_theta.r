delta_theta <- function(xt,mu,V,prop,l,param){
  output = list()
  m = dim(mu)[1]; p = dim(mu)[2]
  delta_pi  = rep(0,p-1)
  aum     = m*(m+3)/2
  delta_muvar = matrix(0,aum,p)

  #param = intermedpi(xt,mu,V,prop)

  if (l != p){
    delta_pi[l] = (mvtnorm::dmvnorm(xt,mean=mu[,l],sigma=V[,,l]) * (sum(param$gt) - param$gt[l] + prop[l]*mvtnorm::dmvnorm(xt,mean=mu[,p],sigma=V[,,p])))/(sum(param$gt))^2
    for (a in setdiff(c(1:(p-1)),l)){
      delta_pi[a] = param$gt[l]*(mvtnorm::dmvnorm(xt,mean=mu[,p],sigma=V[,,p]) - mvtnorm::dmvnorm(xt,mean=mu[,a],sigma=V[,,a]))/(sum(param$gt))^2
    }
  }else if(l == p){
    for (a in setdiff(c(1:(p-1)),l)){
      delta_pi[a] = (-mvtnorm::dmvnorm(xt,mean=mu[,p],sigma=V[,,p])*(sum(param$gt) - param$gt[l]) - param$gt[l]*mvtnorm::dmvnorm(xt,mean=mu[,a],sigma=V[,,a]))/(sum(param$gt))^2
    }
  }


  delta_muvar[,l] = (sum(param$gt) - param$gt[l])*param$gt[l]*param$ct[,l]/(sum(param$gt))^2
  for (a in setdiff(c(1:(p)),l)){
    delta_muvar[,a] = -param$gt[a]*param$gt[l]*param$ct[,l]/(sum(param$gt))^2
  }

  output$delta_pi = delta_pi
  output$delta_muvar = delta_muvar
  return(output)
}
