############################################################
# Buffon's needle problem
# Author: Jose A. Perusquia Cortes
# Affil:  Facultad de Ciencas-UNAM
# Module: Stochastic Simulation 
############################################################

############################################################
# Libraries
library(ggplot2)
library(ggthemes)
library(gganimate)
library(gifski)
library(dplyr)
############################################################

############################################################
# Parameters

# Number of strips
n_strips = 20

# Width of strips
w_strips = 1

# Length of the needle
l_needle = .3
############################################################

############################################################
# Grid
g_x1 = rep(0,n_strips+1)
g_x2 = rep(10,n_strips+1)

g_y1 = seq(0, n_strips, by = w_strips)
g_y2 = g_y1

grid_lines = data.frame(g_x1,g_x2,g_y1,g_y2)

ggplot(data=grid_lines)+
  geom_segment(aes(x=g_x1,xend=g_x2,y=g_y1,yend=g_y2))+
  theme_void()
############################################################

############################################################
# Needles

# Number of needles
n = 1000

# Positions on the grid
set.seed(31415)
cx = runif(n,0,10)
cy = runif(n,0,n_strips)
theta = runif(n,0,pi)

n_x1 = cx - (l_needle/2)*cos(theta)
n_y1 = cy - (l_needle/2)*sin(theta)

n_x2 = cx + (l_needle/2)*cos(theta)
n_y2 = cy + (l_needle/2)*sin(theta)

# Check if it crossed
crosses = floor(n_y1/w_strips) != floor(n_y2/w_strips)
needles = data.frame(n_x1,n_x2,n_y1,n_y2,crosses)
############################################################

############################################################
# Animated plot of the needles on the grid
needles$throw = seq_len(n)
anim=ggplot(data = grid_lines) +
  geom_segment(aes(x = g_x1,
                   xend = g_x2,
                   y = g_y1,
                   yend = g_y2)) +
  theme_void() +
  geom_segment(data = needles,
               show.legend = FALSE,
               aes(x = n_x1,
                   y = n_y1,
                   xend = n_x2,
                   yend = n_y2,
                   colour = crosses)) +
  transition_manual(throw) +
  shadow_trail(distance = 0.001)

animate(anim,fps = 20, duration = 10, 
        renderer = gifski_renderer())

# Save as GIF
anim_save("BuffonsNeedle.gif", 
          animation = last_animation())

# Complete set of needles in one plot
ggplot(data=grid_lines)+
  geom_segment(aes(x=g_x1,xend=g_x2,y=g_y1,yend=g_y2))+
  theme_void()+
  geom_segment(data=needles,show.legend = F,
               aes(x=n_x1,xend=n_x2,y=n_y1,yend=n_y2,
                   col=crosses))
############################################################

############################################################
# Approximate probability of crossing 
n_seq = seq_len(n)
p_hat = cummean(crosses)
p_theoretical = 2*l_needle/(pi*w_strips)
p_hat = data.frame(x=n_seq,y=p_hat)
p_hat$y[n]

ggplot(data=p_hat,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x=expression(n),y=expression(hat(P)))+
  geom_hline(yintercept = p_theoretical,col='darkred')

# Approximate pi
pi_hat = 2*l_needle/(w_strips*p_hat$y)
pi_hat = data.frame(x=n_seq,y=pi_hat)
pi_hat$y[n]

ggplot(data=pi_hat,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x=expression(n),y=expression(hat(pi)))+
  geom_hline(yintercept = pi,col='darkred')
############################################################
