############################################################
# PRNG linear congruential algorithm
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
library(randtoolbox)     # Version 2.0.5
############################################################

############################################################
# Function that creates a sequence of length n of 
# pseudorandom numbers generated from the recursion 
# X_{n+1} = (a*X_n + b) mod M starting from Xin
PRNG = function(a,b,M,n,Xin){
  res = numeric(n+1)
  res[1] = Xin
  for(i in 1:n){
    res[i+1] = (a*res[i] + b )%%M
  }
  return(res)
}
############################################################

############################################################
# First example: X_{n+1} = (5 X_n + 1) mod 32
a   = 5
b   = 1
M   = 32
n   = 31
Xin = 7

# Generate one complete period
res = data.frame(
  iteration = 0:n,
  x = PRNG(a, b, M, n, Xin)
)

# Transform states to midpoints in (0,1)
res = res %>%
  mutate(u = (x + 0.5) / M)
N = nrow(res)

# Q-Q plot
ggplot(res, aes(sample = u)) +
  geom_qq(distribution = stats::qunif) +
  geom_qq_line(
    distribution = stats::qunif,
    colour = "darkred"
  ) +
  labs(
    x = "Teóricos",
    y = "Empíricos"
  ) +
  theme_minimal()

# Histogram
ggplot(res, aes(x = u)) +
  geom_histogram(
    bins = 10,
    colour = "black",
    fill = "darkred"
  ) +
  labs(
    x = expression(U[n]),
    y = "Frecuencia"
  ) +
  theme_minimal()

# Kolmogorov-Smirnov statistic
# Note: the conventional p-value assumes iid continuous observations
ks.test(res$u, "punif")

# Runs test
RunsTest(res$u > 0.5, exact = TRUE)

# ACF
h = floor(N / 4)

acf_obj = acf(
  res$u,
  lag.max = h,
  plot = FALSE
)

# Approximate 95% bounds under the white-noise null
acf_df = data.frame(
  lag   = as.vector(acf_obj$lag),
  acf   = as.vector(acf_obj$acf),
  upper = 1.96 / sqrt(N),
  lower = -1.96 / sqrt(N)
)

# Remove lag zero
acf_df = acf_df %>%
  filter(lag > 0)

ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0) +
  geom_hline(
    yintercept = c(-1.96 / sqrt(N), 1.96 / sqrt(N)),
    linetype = 2,
    colour = "blue"
  ) +
  labs(
    x = expression(h),
    y = "ACF"
  ) +
  theme_minimal()

# Ljung-Box test
Box.test(
  res$u,
  lag = 10,
  type = "Ljung-Box"
)

# Consecutive pairs (U_n, U_{n+1})
lag_pairs = data.frame(
  U_n  = head(res$u, -1),
  U_n1 = tail(res$u, -1)
)

ggplot(lag_pairs, aes(x = U_n, y = U_n1)) +
  geom_point(alpha = 0.6) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    x = expression(U[n]),
    y = expression(U[n+1])
  ) +
  theme_minimal()

# Serial test
serial.test(res$u)
############################################################

############################################################
# Second example: X_{n+1} = (5 X_n + 1) mod 256
a   = 5
b   = 1
M   = 256
n   = 255
Xin = 7

# Generate one complete period
res = data.frame(
  iteration = 0:n,
  x = PRNG(a, b, M, n, Xin)
)

# Transform states to midpoints in (0,1)
res = res %>%
  mutate(u = (x + 0.5) / M)
N = nrow(res)

# Q-Q plot
ggplot(res, aes(sample = u)) +
  geom_qq(distribution = stats::qunif) +
  geom_qq_line(
    distribution = stats::qunif,
    colour = "darkred"
  ) +
  labs(
    x = "Teóricos",
    y = "Empíricos"
  ) +
  theme_minimal()

# Histogram
ggplot(res, aes(x = u)) +
  geom_histogram(
    bins = 10,
    colour = "black",
    fill = "darkred"
  ) +
  labs(
    x = expression(U[n]),
    y = "Frecuencia"
  ) +
  theme_minimal()

# Kolmogorov-Smirnov statistic
# Note: the conventional p-value assumes iid continuous observations
ks.test(res$u, "punif")

# Runs test
RunsTest(res$u > 0.5, exact = TRUE)

# ACF
h = floor(N / 4)

acf_obj = acf(
  res$u,
  lag.max = h,
  plot = FALSE
)

acf_df = data.frame(
  lag   = as.vector(acf_obj$lag),
  acf   = as.vector(acf_obj$acf),
  upper = 1.96 / sqrt(N),
  lower = -1.96 / sqrt(N)
)

# Remove lag zero
acf_df = acf_df %>%
  filter(lag > 0)

ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0) +
  geom_hline(
    yintercept = c(-1.96 / sqrt(N), 1.96 / sqrt(N)),
    linetype = 2,
    colour = "blue"
  ) +
  labs(
    x = expression(h),
    y = "ACF"
  ) +
  theme_minimal()

# Ljung-Box test
Box.test(
  res$u,
  lag = 10,
  type = "Ljung-Box"
)

# Consecutive pairs (U_n, U_{n+1})
lag_pairs = data.frame(
  U_n  = head(res$u, -1),
  U_n1 = tail(res$u, -1)
)

ggplot(lag_pairs, aes(x = U_n, y = U_n1)) +
  geom_point(alpha = 0.6) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    x = expression(U[n]),
    y = expression(U[n+1])
  ) +
  theme_minimal()

# Serial test
serial.test(res$u)
############################################################

############################################################
# Mersenne twister default R algorithm

n = 256

# Generate the sample
set.seed(31415)
res = data.frame(
  iteration = 1:n,
  u = runif(n)
)
N = nrow(res)

# Q-Q plot
ggplot(res, aes(sample = u)) +
  geom_qq(distribution = stats::qunif) +
  geom_qq_line(
    distribution = stats::qunif,
    colour = "darkred"
  ) +
  labs(
    x = "Teóricos",
    y = "Empíricos"
  ) +
  theme_minimal()

# Histogram
ggplot(res, aes(x = u)) +
  geom_histogram(
    bins = 10,
    colour = "black",
    fill = "darkred"
  ) +
  labs(
    x = expression(U[n]),
    y = "Frecuencia"
  ) +
  theme_minimal()

# Kolmogorov-Smirnov statistic
# Note: the conventional p-value assumes iid continuous observations
ks.test(res$u, "punif")

# Runs test
RunsTest(res$u > 0.5, exact = TRUE)

# ACF
h = floor(N / 4)

acf_obj = acf(
  res$u,
  lag.max = h,
  plot = FALSE
)

acf_df = data.frame(
  lag   = as.vector(acf_obj$lag),
  acf   = as.vector(acf_obj$acf),
  upper = 1.96 / sqrt(N),
  lower = -1.96 / sqrt(N)
)

# Remove lag zero
acf_df = acf_df %>%
  filter(lag > 0)

ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0) +
  geom_hline(
    yintercept = c(-1.96 / sqrt(N), 1.96 / sqrt(N)),
    linetype = 2,
    colour = "blue"
  ) +
  labs(
    x = expression(h),
    y = "ACF"
  ) +
  theme_minimal()

# Ljung-Box test
Box.test(
  res$u,
  lag = 10,
  type = "Ljung-Box"
)

# Consecutive pairs (U_n, U_{n+1})
lag_pairs = data.frame(
  U_n  = head(res$u, -1),
  U_n1 = tail(res$u, -1)
)

ggplot(lag_pairs, aes(x = U_n, y = U_n1)) +
  geom_point(alpha = 0.6) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    x = expression(U[n]),
    y = expression(U[n+1])
  ) +
  theme_minimal()

# Serial test
serial.test(res$u)
############################################################