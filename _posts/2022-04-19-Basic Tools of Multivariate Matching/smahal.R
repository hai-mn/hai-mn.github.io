#################################################################################
#################################################################################
# The function smahal creates a distance matrix using the rank-based Mahalanobis#
# distance; It requires that the MASS package has been installed.               #
# where:                                                                        #
#  z:= a vector of length(z)=n; z[i]=1 for treated and z[i]=0 for control,      #
#  X:= an n×k matrix of covariates,                                             #
#  dmat:= a distance matrix created by mahal or smahal,                         #
#  p:= a vector of length(p)=n typically containing a propensity score, and     #
#  f:= a vector of length(f)=n with a few values used in almost exact matching. #
#################################################################################
#################################################################################;

smahal <- function(z, X){ # X:= an n×k matrix of covariates
  
  # Rank-based Mahalanobis distance
  X <- as.matrix(X)
  n <- dim(X)[1] # number of observations of X
  rownames(X) <- 1:n
  
  k <- dim(X)[2] # number of covariates of X
  m <- sum(z) # number of treated
  
  for (j in 1:k){X[,j] <- rank(X[,j])}
  cv <- cov(X)
  
  vuntied <- var(1:n)
  
  rat <- sqrt(vuntied/diag(cv))
  
  cv <- diag(rat)%*%cv%*%diag(rat)
  
  out<-matrix(NA,m,n-m)
  Xc<-X[z==0,]
  Xt<-X[z==1,]
  
  rownames(out)<-rownames(X)[z==1]
  colnames(out)<-rownames(X)[z==0]
  
  library(MASS)
  
  icov<-ginv(cv)
  
  for (i in 1:m) {
    out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)}
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
#smahal.gen.tox<-smahal(z,X)
#print(smahal.gen.tox[,1:7])
