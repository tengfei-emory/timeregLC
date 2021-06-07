# timeregLC
A Time-Dependent Structural Model Between Latent Classes and Competing Risks Outcomes

# Installation Guide
```{r}
# install.packages("devtools")
devtools::install_github("tengfei-emory/timeregLC")
library(timeregLC)
```
Currently `timeregLC` supports R version >= 3.3.0.

# Example: analyze a simulated dataset

## Data simulation

Function `simulation` can be used to generate a dataset with baseline covariates and competing risks.

```{r}
# The following example specifies all required parameters and generate a dataset with three latent classes.

# Regression parameter of latent class effect in the structural competing risks model
lambda <- c(0.5,0.5,-1)

# Latent class proportion
pi = c(0.3,0.35,0.35)

# Mean vectors for the three classes: (1,1), (2.5,2.5) and (4,4)
mu = matrix(c(1,1,2.5,2.5,4,4),nrow=2,ncol=3)

# Covariance matrices for the three classes (as a list)
sigma1 = matrix(c(0.36,0.27,0.27,0.81),2,2)
sigma2 = matrix(c(0.49,0.504,0.504,0.64),2,2)
sigma3 = matrix(c(0.25,0.225,0.225,0.25),2,2)
sigma = list(sigma1,sigma2,sigma3)

# Parameter associated with competing risks distribution
p.cif = 0.66

# Lower bound and upper bound for uniformly distributed censoring time
cl=0.19
cu=1.09

# Main function of simulation. Here sample size is set as 500.
dat=simulation(500,pi,mu,sigma,lambda,p.cif,cl,cu)
```
Specifically, it returns a data frame of 2 latent classes with 2 baseline covariates (`Y.1` and `Y.2`), time of competing risks (`ftime`), and failure types (`fstatus`). Failure types include type `1`, `2`, censored `0`, or missing `NA`. 

## Model fitting

The analysis for the dataset `dat` can be conducted by running `timereg` function:

```{r}
library(timereg)
# create an event object 
event = Event(0,dat$ftime,dat$fstatus)

# specify the baseline covariate matrix
covariates=cbind(dat$Y.1,dat$Y.2)

# run main algorithm
fit.timeregLC <- timeregLC(event,covariates,inference=T,C=3,d=1,timepoints=NULL,
                    control.optim=list(reltol=.00001,strategy=2,itermax=1000,trace=F),
                    verbose=T)
```                   

The output list `fit.timeregLC` contains the following information:

`lambda`: time-dependent point estimates at specified time points

`Sigma`: time-dependent estimates for the asymptotic covariance estimates

Please refer to the documentation associated with the package for more details. Please report issues under this GitHub repository (tengfei-emory/timeregLC).

# References

Fei, Hanfelt and Peng (under revision). Evaluating the Association between Latent Classes and Competing Risks Outcomes with Multi-Phenotype Data.
