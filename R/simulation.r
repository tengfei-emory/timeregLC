#' Simulate an example dataset
#'
#' @description Simulate an example dataset.
#' @param n Sample size
#' @param p.class L-dimensional marginal probability vector of latent classes. L is the number of latent classes.
#' @param mu Latent class-specific mean vectors for baseline covariates, stored in a matrix of dimension L by n. The ith column of the matrix is the mean for the ith latent class.
#' @param sigma Latent class-specific covariance matrices for baseline covariates, stored in a list. See example for detailed specifications.
#' @param lambda Underlying true latent class effect coefficients, vector of L
#' @param p.cif A scalar parameter to adjust the subdistribution survival function under Fine and Gray (1999)'s model.
#' @param cl Lower bound of censoring time.
#' @param cu Upper bound of censoring time.
#' @return Returns a data frame with event time, failure type and baseline covariates. Failure type 0 indicates censoring. Missing failure types for uncensored subjects are denoted by NAs.
#' @author Teng Fei. Email: tfei@emory.edu
#' @details Simulate latent classes from multinomial distribution, then generate baseline covariates based on finite mixture model (McLachlan and Peel, 2000) and competing risks data with 2 failure types based on Fine and Gray (1999)'s proportional subdistribution hazards model. Censoring is independently generated by uniform distribution (cl, cu). For uncensored subjects, missingness of failure type is independently generated from a logit function depending on baseline covariates and failure time.
#' @references Fine and Gray (1999). A proportional hazards model for the subdistribution of a competing risk. Journal of the American statistical association 94, 496-509.\cr
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
#' @export

simulation <- function(n,p.class,mu,sigma,lambda,p.cif,cl,cu){

  require(MASS)
  require(mvtnorm)

  C <- length(p.class)

  #Generate the true latent class 'delta'
  vector.Delta <- rmultinom(n,1,p.class)
  delta <- rep(0,n)
  for (i in 1:n){
    for (c in 1:C){
      delta[i] = delta[i] + I(vector.Delta[c,i]==1)*c
    }
  }

  #Generate covariates y
  y = matrix(ncol=dim(mu)[1],nrow=n)
  for (i in 1:n){
    y[i,] = mvrnorm(1,mu=mu[,delta[i]],Sigma=sigma[[delta[i]]])
  }

  ### Competing Risk Data ###

  Delta_star <- vector.Delta
  Delta_star[1,] = 1

  #generate the types of outcome

  P1 <- 1-(1-p.cif)^(exp(lambda%*%Delta_star))
  epsilon <- rep(0,n)
  for (i in 1:n){
    epsilon[i] <- 2 - rbinom(1,1,P1[i])
  }

  #generate the event time based on the type of outcome
  T <- rep(0,n)
  u <- runif(n)
  for (i in 1:n){
    if (epsilon[i] == 1){
      T[i] <- -log(1 - (1 - (1-u[i]*(1-(1-p.cif)^exp(lambda%*%Delta_star[,i])))^(1/exp(lambda%*%Delta_star[,i])))/p.cif)
    }
    if (epsilon[i] == 2){
      T[i] <- -log((1-u[i])^(1/exp(lambda%*%Delta_star[,i])))
    }
  }

  #generate censoring time
  c <- runif(n,cl,cu)

  #observed time
  x <- T*I(T<=c) + c*I(T>c)

  # outcome
  D <- 0*I(x == c) + epsilon*I(x < c)

  # prob of complete case
  #r <- exp(-0.35*y[,1]+0.5*y[,2]-x)/(1+exp(-0.35*y[,1]+0.5*y[,2]-x))
  r <- exp(0.25*y[,1]+0.5*y[,2])/(1+exp(0.25*y[,1]+0.5*y[,2]))
  #r <- rep(0.5,length(y[,1]))

  # complete case indicator
  R <- rep(0,n)
  for (i in 1:n){
    R[i] <- rbinom(1,1,r[i])
  }
  R = I(x<c)*R + I(x==c)

  # observed outcome

  Dobs = I(R==1)*D + I(R==0)*0

  DobsNaive = Dobs
  DobsNaive[R == 0] = NA

  # censoring

  C = 0*I(x == c) + 1*I(x < c)

  return(data.frame(ftime = x, fstatus = DobsNaive, Y = y))
}
