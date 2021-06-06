#' A Time-Dependent Structural Model Between Latent Classes and Competing Risks Outcomes
#'
#' @description Fits a direct binomial model (Scheike et al, 2008) with latent classes as covariates to investigate latent class effects on competing risks.
#' @param event Event object created by Event function under timereg package.
#' @param covariates Covariates matrix with n rows, where n is the sample size. Currently only continuous covariates are supported.
#' @param inference Logical indicator. Whether compute asymptotic covariance function for the time-varying latent class effect point estimates.
#' @param C Scalar indicating the number of latent classes. If specified as NULL, then C is established by the finite mixture model (McLachlan and Peel, 2000) with smallest BIC.
#' @param d Specify the cause-d failure type to consider in the competing risk model.
#' @param timepoints A vector of time points t when the model is fitted to obtain time-varying latent class effect $lambda(t)$. If not specified, then use ten evenly divided points from the first event to the last event.
#' @param control.optim Controlling parameters for DEoptim package (Mullen et al., 2011).
#' @param verbose Whether show notifications on completion of each step.
#' @return Returns point estimates for time-varying latent class effect coefficients $hat{lambda}(t)$. If inference=T, returns a list consisting $hat{lambda}(t)$ and the asymptotic covariance function of $hat{lambda}(t)$.
#' @author Teng Fei. Email: tfei@emory.edu
#' @references Scheike and Zhang (2011). Analyzing competing risk data using the r timereg package. Journal of statistical software 38. \cr
#' Scheike et al. (2008). Predicting cumulative incidence probability by direct binomial regression. Biometrika 95, 205-220. \cr
#' Mullen et al. (2011). Deoptim: An r package for global optimization by differential evolution. Journal of Statistical Software 40. \cr
#' McLachlan and Peel (2000). Finite Mixture Models. John Wiley & Sons. \cr
#' Fei, Hanfelt and Peng (under revision). Evaluating the Association between Latent Classes and  Competing Risks Outcomes with Multi-Phenotype Data.
#' @examples
#'
#' #The following example specifies all required parameters and generate a dataset.
#' lambda <- c(0.5,0.5,-1)
#' pi = c(0.3,0.35,0.35)
#' mu = matrix(c(1,1,2.5,2.5,4,4),nrow=2,ncol=3)
#' sigma1 = matrix(c(0.36,0.27,0.27,0.81),2,2)
#' sigma2 = matrix(c(0.49,0.504,0.504,0.64),2,2)
#' sigma3 = matrix(c(0.25,0.225,0.225,0.25),2,2)
#' sigma = list(sigma1,sigma2,sigma3)
#' p.cif = 0.66
#' cl=0.19
#' cu=1.09
#' dat=simulation(500,pi,mu,sigma,lambda,p.cif,cl,cu)
#'
#' #Fit the model
#' library(timereg)
#' event = Event(0,dat$ftime,dat$fstatus)
#' covariates=cbind(dat$Y.1,dat$Y.2)
#' output <- timeregLC(event,covariates,inference=T,C=3,d=1,timepoints=NULL,
#'                     control.optim=list(reltol=.00001,strategy=2,itermax=1000,trace=F),
#'                     verbose=T)
#'
#' @export


timeregLC <- function(event,covariates,inference=T,C=NULL,d=1,timepoints=NULL,
                      control.optim=list(reltol=.00001,strategy=2,itermax=1000,trace=F),
                      verbose=T){

  require(timereg)
  require(survival)
  require(DEoptim)
  require(mclust)
  require(MASS)

  #link function
  g <<- function(x) 1-exp(-exp(x))
  dg <<- function(x) exp(x-exp(x))

  #sort data by failure time
  covariates = covariates[order(event[,2]-event[,1]),]
  event = event[order(event[,2]-event[,1]),]

  #observed survival time
  ftime = event[,2]-event[,1]

  #observed failure type, subject to missing
  fstatus = event[,3]

  #censoring 1=censored
  cstatus = I(fstatus==0)
  cstatus[is.na(cstatus)] = FALSE
  cstatus = as.integer(cstatus)
  #cens = 1-as.integer(cstatus)

  #missingness indicator
  R = 1-is.na(fstatus)

  ### LCA by finite mixture model ###

  if (is.null(C)){
    lca = mclust::Mclust(covariates,modelNames='VVV',verbose=F)
    C = lca$G
  }else{
    lca = mclust::Mclust(covariates,G=C,modelNames='VVV',verbose=F)
  }

  #rank latent classes
  lca$z <- lca$z[,order(lca$parameters$mean[1,]+lca$parameters$mean[2,])]
  lca$parameters$pro <- lca$parameters$pro[order(lca$parameters$mean[1,]+lca$parameters$mean[2,])]
  lca$parameters$variance$sigma <- lca$parameters$variance$sigma[,,order(lca$parameters$mean[1,]+lca$parameters$mean[2,])]
  lca$parameters$mean <- lca$parameters$mean[,order(lca$parameters$mean[1,]+lca$parameters$mean[2,])]


  #influence function
  pi1 = lca$parameters$pro
  mean = lca$parameters$mean
  var = lca$parameters$variance$sigma
  scorevec <- multiscorepi(X=t(covariates),mu=mean,V=var,prop=pi1)
  infomat <- multinfopi(X=t(covariates),mu=mean,V=var,prop=pi1)
  inv.infomat <- solve(infomat$infopi)
  allscore <- c(scorevec$scpi,as.vector(scorevec$scmuvar))
  score <- scorevec$indscore
  invinfo <- solve(infomat$infopi/nrow(covariates))
  allphihat <- inv.infomat %*% allscore
  if (verbose==T){
    cat('LCA: Complete','\n')
  }

  ### Logistic regression model for missing failure types ###

  LR.model <- glm(R[cstatus==0] ~ covariates[cstatus==0,] + ftime[cstatus==0], family = binomial(link='logit'))
  rhat1 <- predict.glm(LR.model,type = "response")
  rhat <- rep(1,nrow(event))
  rhat[cstatus==0] <- rhat1
  Ir1 <- rhat*(1-rhat)*(1-cstatus)
  Sr <- (R - rhat)*(1-cstatus)
  Ir = 0
  for (a in 1:nrow(event)){
    Ir = Ir + Ir1[a]*(c(1,covariates[a,],dat$ftime[a]) %o% c(1,covariates[a,],dat$ftime[a]))
  }
  Ir = Ir/nrow(event)
  if (verbose==T){
    cat('Logistic Regression for missing failure type: Complete','\n')
  }

  ### Point estimates ###

  #mark missing as 0
  fstatus[is.na(fstatus)] = 0
  censored = 1-cstatus

  #inverse probability censoring weighting
  GhatX <- survival::survfit(Surv(ftime,cstatus)~1)
  Gfit <- cbind(GhatX$time, GhatX$surv)
  Gfit <- rbind(c(0, 1), Gfit)
  G <- timereg::Cpred(Gfit, ftime, strict = TRUE)[,2]

  #define time points t when \lambda(t) is solved
  if (is.null(timepoints)){
    timepoints = seq(from=min(ftime[fstatus==d]),to=max(ftime[fstatus==d]),length.out=10)
  }

  #number of latent classes C
  #create matrix \kappa
  kappa <- cbind(rep(1,C),rbind(0,diag(C-1)))

  gsum1 <<- function(lambda) as.vector(g(lambda%*%t(kappa))%*%t(lca$z))

  #initialize output point estimates
  lambda_results <- matrix(0,nrow=length(timepoints),ncol=C)

  for (i in 1:length(timepoints)){
    t = timepoints[i]
    f <- function(lambda){
      sumg <- gsum1(lambda)
      -sum(sumg*I(ftime<=t & fstatus==d)/((G+0.001)*(censored*(rhat+0.001)+1-censored)) - 0.5*sumg^2)
    }
    temp.result <- DEoptim(f,lower=rep(-5,C),upper=rep(5,C),control=control.optim)
    lambda_results[i,] <- temp.result$optim$bestmem
  }

  if (verbose==T){
    cat('Point estimates: Complete','\n')
  }

  ### Variance estimates ###

  if(inference==T){

    n = nrow(event)

    #initialize variance estimates output Sigma
    Sigma = array(0,dim=c(length(timepoints),C,C))

    #create a perturbation
    epsilon2 <- matrix(0, C, C)
    diag(epsilon2) <- 1e-10

    for (i in 1:length(timepoints)){
      t = timepoints[i]
      lambda = lambda_results[i,]

      JH = JHsolver(g,dg,lambda,lca)
      J = JH$J
      H = JH$H
      Psi = JH$Psi
      Psi_lambda = JH$Psi_lambda
      W = 0

      iota = 0
      for (k in 1:n){
        iota = iota + censored[k]*(Psi_lambda[k,] %o% c(1,covariates[k,],ftime[k]))*I(ftime[k] <= t & fstatus[k]==d)*rhat[k]*Sr[k]/(G[k]*rhat[k]^2+0.001)
      }
      iota = iota/n

      for (j in 1:n){
        xi_j = 0
        if (fstatus[j] == 0){
          denom = sum(dat$ftime >= ftime[j])
          for (k in j:n){
            xi_j = xi_j + Psi_lambda[k,]*I(ftime[k] >= ftime[j])*I(ftime[k] <= t & fstatus[k]==d)/(G[k]+0.001)
          }
          xi_j = xi_j/denom
        }

        iota_j = solve(Ir+1e-10) %*% (Sr[j]*c(1,covariates[j,],dat$ftime[j]))
        phihat = invinfo %*% score[[j]]
        w = n*solve(J+epsilon2) %*% (Psi_lambda[j,]*(I(ftime[j]<=t & fstatus[j]==d)/((G[j]+1e-3)*(censored[j]*rhat[j]+1-censored[j])) - Psi[j]) - xi_j - iota%*%iota_j - H %*% phihat)

        W = W + w %*% t(w)
      }

      W = W/n
      Sigma[i,,] = W
    }

    if (verbose==T){
      cat('Variance estimates: Complete','\n')
    }

    output <- list(lambda=lambda_results,Sigma=Sigma)
    return(output)

  }else{

    return(lambda_results)

  }

}
