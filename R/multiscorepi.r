multiscorepi <- function(X,mu,V,prop){
  output = list()
  m = dim(mu)[1]; p = dim(mu)[2]
  n       = dim(X)[2]
  scpi    = rep(0,p-1)
  aum     = m*(m+3)/2
  scmuvar = matrix(0,aum,p)
  sinfopi = matrix(0,aum*p+p-1,aum*p+p-1)
  indscore   =list()
  
  for (t in 1:n){
    param = intermedpi(X[,t],mu,V,prop)
    sinfa = param$at %*% param$alphat
    sinfc = param$ct %*% diag(as.vector(param$alphat))
    scpi = scpi + sinfa
    scmuvar = scmuvar + sinfc
    vecsinfc = as.vector(sinfc)
    sinfopi = sinfopi - rbind(cbind(sinfa %*% t(sinfa), sinfa %*% t(vecsinfc)), cbind(vecsinfc %*% t(sinfa), vecsinfc %*% t(vecsinfc)))
    indscore[[t]] = c(sinfa,as.vector(sinfc))
  }
  sinfopi = -sinfopi
  
  output$scpi = scpi 
  output$scmuvar = scmuvar
  output$sinfopi = sinfopi
  output$indscore = indscore
  return(output)
} 