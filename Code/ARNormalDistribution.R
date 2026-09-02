############################################################
# Acceptance-rejection for standard normal distribution
# with Cauchy instrumental distribution
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
# Acceptance-rejection algorithm
NormalAR = function(n){
  
  # Rejection constant
  M = sqrt(2*pi/exp(1))
  
  # Objects to store the proposals
  X = numeric(0)
  U = numeric(0)
  accepted = logical(0)
  
  # Number of accepted observations
  counter = 0
  
  while(counter < n){
    
    # Proposal from the instrumental distribution
    x = rcauchy(1)
    u = runif(1)
    
    # Ratio f(x)/g(x)
    h = sqrt(pi/2)*(1+x^2)*exp(-x^2/2)
    
    # Acceptance probability
    p = h/M
    
    # Acceptance/rejection
    accept = u <= p
    
    # Store proposal
    X = c(X,x)
    U = c(U,u)
    accepted = c(accepted,accept)
    
    if(accept)
      counter = counter+1
  }
  
  # Results
  df = data.frame(
    proposal = X,
    uniforms = U,
    accepted = accepted
  )
  
  return(df)
}

# Second version where we propose with batched to make it
# computationally more efficient
NormalARfast = function(n){
  
  # Rejection constant
  M = sqrt(2*pi/exp(1))
  
  # Acceptance probability
  p = 1/M
  
  # Output vector
  Y = numeric(n)
  
  # Number accepted so far
  counter = 0
  
  while(counter < n){
    
    # Number still required
    remaining = n-counter
    
    # Generate slightly more than the expected number required
    m = ceiling(1.1*remaining/p)
    
    # Vector of proposals
    X = rcauchy(m)
    U = runif(m)
    
    # Acceptance probability
    prob = sqrt(exp(1))/2*(1+X^2)*exp(-X^2/2)
    
    # Accepted proposals
    accepted = X[U <= prob]
    
    # Number we can use
    k = min(length(accepted),remaining)
    
    if(k > 0){
      Y[(counter+1):(counter+k)] = accepted[1:k]
      counter = counter+k
    }
  }
  
  return(Y)
}

############################################################

############################################################
# We analyse the results
set.seed(314159)
normalSample = NormalAR(500)

# Theoretical acceptance probability
M = sqrt(2*pi/exp(1))
1/M

# Empirical acceptance probability
mean(normalSample$accepted)

# Dataframe of accepted and rejected points
ARpoints = data.frame(
  x=normalSample$proposal,
  y=M*normalSample$uniforms*dcauchy(normalSample$proposal),
  accepted = normalSample$accepted
)

# Support for plotting the densities
sopX = seq(-6,6,by=0.01)
n_sopX = length(sopX)
gx = M*dcauchy(sopX)
fx = dnorm(sopX)

# Put it in a data frame
x = rep(sopX,2)
f = c(gx,fx)
cols = c(rep('Cauchy',n_sopX),rep('Normal',n_sopX))
df = data.frame(x = x,y = f,cols = cols)

# Generate vertical coordinates:
# conditional on X=x, Y is uniform on (0,Mg(x))
ARpoints = data.frame(
  x = normalSample$proposal,
  y = M*normalSample$uniforms*dcauchy(normalSample$proposal),
  accepted = normalSample$accepted
)

# Plot
ggplot(data = df, aes(x=x,y=y))+
  geom_line(aes(colour=cols))+
  geom_point(
    data = ARpoints,
    aes(x=x,y=y,shape=accepted),
    fill = "white",
    size = 1,
    show.legend = FALSE
  )+
  scale_shape_manual(values=c(21,24))+
  theme_minimal()+
  coord_cartesian(xlim=c(-6,6))+
  scale_color_discrete(name = "", 
                       labels = c("M*Cauchy", "Normal"))+
  labs(x='',y='')

# Accepted sample
normalSample=normalSample%>%
  filter(accepted)%>%
  rename(sample=proposal)

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
# The NormalAR is not intended to compete against the highly-
# optimized rnorm code. It's a transparent pedagogical
# implementation. 
microbenchmark(
  AR = NormalAR(1000),
  ARFast = NormalARfast(1000),
  RDefault=rnorm(1000), 
  times=1000)
############################################################