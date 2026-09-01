############################################################
# Box Muller transformation for standar normal distribution
# Author: Jose A. Perusquia Cortes
# Affil:  Facultad de Ciencias-UNAM
# Module: Stochastic Simulation 
############################################################

############################################################
# Libraries
library(ggplot2)         # Version 4.0.2
library(ggthemes)        # Version 5.2.0
library(dplyr)           # Version 1.2.0
library(DescTools)       # Version 0.99.60
library(nortest)         # Version 1.0-4
library(microbenchmark)  # Version 1.5.0
############################################################

############################################################
# Function that creates standard normal random variables 
# using Box-Muller transformation
boxMuller = function(n){

  # Number of pairs required
  m = ceiling(n/2)
  
  # Generate the uniform random variables
  U = runif(m)
  V = runif(m)
  
  # Box Muller transformation
  X = sqrt(-2*log(U))*sin(2*pi*V)
  Y = sqrt(-2*log(U))*cos(2*pi*V)
    
  # Interlace pairs and retain the first n observations
  res = as.vector(rbind(X,Y))[1:n]
  
  # Results
  df = data.frame(sample = res)
  
  return(df)
}
############################################################

############################################################
# We analyse the results
set.seed(314159)
normalSample = boxMuller(1000)

# Descriptive statistics 
summary(normalSample$sample)
Kurt(normalSample$sample)
Skew(normalSample$sample)

# Histogram
ggplot(data=normalSample,aes(x=sample,y=after_stat(density)))+
  geom_histogram(bins=20,col='black',fill='darkred')+
  labs(x='',y='')+
  theme_minimal()

# Boxplot
ggplot(data=normalSample,aes(x=sample))+
  geom_boxplot(col='black',fill='darkred')+
  labs(x='',y='')+
  theme_minimal()

# qqplot
ggplot(data=normalSample,aes(sample=sample))+
  geom_qq(distribution = qnorm)+
  geom_qq_line(distribution = qnorm,col='red')+
  labs(x='',y='')+
  theme_minimal()

# Goodness of fit tests
# Does not validate the algorithm, just does not provide evidence
# against normality.
ad.test(normalSample$sample)
shapiro.test(normalSample$sample)

# Benchmark against rnorm default
microbenchmark(
  BM = boxMuller(100000),
  RDefault=rnorm(100000), 
  times=1000)
############################################################