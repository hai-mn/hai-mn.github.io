#################################################################################
#################################################################################
# The function mahal creates a distance matrix using the Mahalanobis distance;  #
# It requires that the MASS package has been installed.                         #
# where:                                                                        #
#  z:= a vector of length(z)=n; z[i]=1 for treated and z[i]=0 for control,      #
#  X:= an n×k matrix of covariates,                                             #
#  dmat:= a distance matrix created by mahal or smahal,                         #
#  p:= a vector of length(p)=n typically containing a propensity score, and     #
#  f:= a vector of length(f)=n with a few values used in almost exact matching. #
#################################################################################
#################################################################################;

mahal <- function(z, X){ # X:= an n×k matrix of covariates
  X <- as.matrix(X)
  n <- dim(X)[1] # number of observations of X
  rownames(X) <- 1:n
  
  k <- dim(X)[2] # number of covariates of X
  m <- sum(z) # number of treated
  
  cv <- cov(X) # variance-covariance of matrix X
  out <- matrix(NA, m, n-m) # a NULL matrix with `number of treated` rows and `number of controls` columns
  
  Xt <- X[z==1,]
  Xc <- X[z==0,]
  
  ## name row and column of out matrix
  rownames(out) <- rownames(X)[z==1]
  colnames(out) <- rownames(X)[z==0]
  
  require(MASS)
  icov <- ginv(cv) # inverse matrix
  for (i in 1:m){
    out[i,] <- mahalanobis(Xc, Xt[i,], icov, inverted = T)}
  out
}

## An example: produce the Mahalanobis distance for the Welders data
#require(readxl)
#gen.tox <- read_excel("Data/table81.xlsx", sheet = "data")
#require(tidyverse)
#z=ifelse(gen.tox$Group=='Welders',1,0)
#X <- gen.tox %>%
#  mutate(R=ifelse(Race=="AA",1,0), S=ifelse(Smoker=="Y",1,0)) %>%
#  dplyr::select(Age,R,S) %>%
#  collect()
#mahal.gen.tox<-(mahal(z,X))
#print(mahal.gen.tox[,1:7])
