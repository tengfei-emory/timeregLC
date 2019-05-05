JHsolver <- function(g,dg,lambda,lca){
  output <- list()
  C = length(lambda)
  n = dim(lca$z)[1]
  kappa <- cbind(rep(1,C),rbind(0,diag(C-1)))

  X = t(lca$data)
  m = dim(lca$parameters$mean)[1]

  J = 0
  H = 0

  Psi = rep(0,n)
  Psi_lambda = matrix(0,n,C)

  for (j in 1:n){
    psi = rep(0,1)
    psi_lambda = rep(0,C)
    psi_theta = rep(0,C-1+C*(m*(m+3)/2))
    param = intermedpi(X[,j],lca$parameters$mean,lca$parameters$variance$sigma,lca$parameters$pro)
    for (k in 1:C){
      psi = psi + as.vector(g(lambda%*%t(t(kappa[k,]))))*lca$z[j,k]
      psi_lambda = psi_lambda + as.vector(dg(lambda%*%t(t(kappa[k,]))))*lca$z[j,k]*kappa[k,]
      delta_theta_val = delta_theta(X[,j],lca$parameters$mean,lca$parameters$variance$sigma,lca$parameters$pro,k,param)
      delta_theta_val = c(delta_theta_val$delta_pi,as.vector(delta_theta_val$delta_muvar))
      psi_theta = psi_theta + as.vector(dg(lambda%*%t(t(kappa[k,]))))*delta_theta_val
    }

    Psi[j] = psi
    Psi_lambda[j,] = psi_lambda

    J = J + psi_lambda %o% psi_lambda
    H = H + psi_lambda %o% psi_theta
  }

  J = J
  H = H/n

  output$J = J
  output$H = H
  output$Psi = Psi
  output$Psi_lambda = Psi_lambda

  return(output)
}
