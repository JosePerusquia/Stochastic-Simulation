############################################################
# Multivariate normal distribution                                
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Stochastic simulation
############################################################

############################################################
# Required libraries   
library(ggplot2)            # Version 3.5.2
library(GGally)             # Version 2.2.1
library(ggthemes)           # Version 5.1.0
library(scatterplot3d)      # Version 0.3-44
library(expm)               # Version 1.0-0
library(MVN)                # Version 6.3
library(here)               # Version 1.0.1
library(mvtnorm)            # Version 1.3-6
library(car)                # Version 3.1-5
############################################################

############################################################
# Source code for functions used to perform a preliminar 
# analysis on the distributional assumption of multivariate
# normality
source(here('UnivariateGraphicalTests.R'))
############################################################

############################################################
# Multivariate normal distribution 
MultNormal = function(n,mu,Sigma){
  
  # Dimensions
  p = length(mu)
  
  # Check the dimensions
  if(!all(dim(Sigma) == c(p,p))){
    stop("Dimensions of mu and Sigma are incompatible")
  }
  
  # We check that Sigma is positive definite
  eigenValues = eigen(Sigma, symmetric = TRUE)$values
  
  if(any(eigenValues <= 0)){
    stop("Sigma is not positive definite")
  }
  
  # We check that the matrix is symmetric
  if(!isSymmetric(Sigma)){
    stop("Sigma must be symmetric")
  }
  
  # Obtain the square root of Sigma, that is, the matrix
  # A such that AA^T = Sigma
  R = chol(Sigma)
  A = t(R)
  
  # Generate the random sample and transform it
  Z = matrix(rnorm(n*p), nrow=n, ncol=p)
  res = Z %*% t(A)
  res = sweep(res, 2, mu, "+")
  
  res = data.frame(x = res)
  return(res)
}
############################################################

############################################################
# Random sample of a bivariate normal distribution
mu=c(0,0)
Sigma=matrix(c(3,1,1,3),byrow = T,nrow=2)

set.seed(31415)
multnorm.sample = MultNormal(1000,mu,Sigma)

# Scatterplot
dens=dmvnorm(multnorm.sample,mean=mu,sigma=Sigma)
scatterplot3d(multnorm.sample[,1],
              multnorm.sample[,2],
              dens,pch=1,
              color="blue",cex.symbols = 0.7,
              type="p",xlab=expression(x[1]),
              ylab=expression(x[2]),zlab=expression(f(x)))

# Tests 
univariateNormalityPlots(multnorm.sample)
quadraticFormPlot(multnorm.sample)

mvn_HZ=mvn(multnorm.sample,mvn_test="hz")
mvn_HZ$multivariate_normality

mvn_RS=mvn(multnorm.sample,mvn_test="royston")
mvn_RS$multivariate_normality

mvn_MA=mvn(multnorm.sample,mvn_test="mardia")
mvn_MA$multivariate_normality
############################################################

############################################################
# Randnomness of Wishart distribution through the ellipses 
mu = c(0,0)
Sigma=diag(1,2)

# 4 random samples
X1 = as.matrix(MultNormal(100,mu,Sigma))
X2 = as.matrix(MultNormal(100,mu,Sigma))
X3 = as.matrix(MultNormal(100,mu,Sigma))
X4 = as.matrix(MultNormal(100,mu,Sigma))

# Wishart matrices
W1 = t(X1)%*%X1
W2 = t(X2)%*%X2
W3 = t(X3)%*%X3
W4 = t(X4)%*%X4

# Random ellipses
ellipse(c(0, 0), shape=W1, radius=1, col="red",lty=2,add=F)
ellipse(c(0, 0), shape=W2, radius=1, col="blue", lty=2)
ellipse(c(0, 0), shape=W3, radius=1, col="green", lty=2)
ellipse(c(0, 0), shape=W4, radius=1, col="purple", lty=2)
############################################################