############################################################
# Markov Chain : Rat maze
# Author: Jose A. Perusquia Cortes
# Affil:  Facultad de Ciencias - UNAM
# Module: Stochastic Simulation 
############################################################

############################################################
# Libraries
library(ggplot2)            # Version 4.0.2
library(ggthemes)           # Version 5.2.0
library(dplyr)              # Version 1.2.0
library(gganimate)          # Version 1.0.11
library(gifski)             # Version 1.32.0-2
############################################################

############################################################
# 3 x 3 Maze

# States
nodes = expand.grid(x = 1:3,y = 1:3)

# Horizontal paths
horizontal = expand.grid(x = 1:2,y = 1:3)
horizontal$xend = horizontal$x + 1
horizontal$yend = horizontal$y

# Vertical paths
vertical = expand.grid(x = 1:3,y = 1:2)
vertical$xend = vertical$x
vertical$yend = vertical$y + 1

# Combine all paths
paths = rbind(horizontal, vertical)

# Plot labyrinth
maze = ggplot() +
  # Available paths
  geom_segment(data = paths,
    aes(x = x, y = y,xend = xend, yend = yend),
    linewidth = 1)+
  # States
  geom_point(data = nodes,aes(x = x, y = y),size = 20,
    shape = 21,fill = "white",stroke = 1) +
  # Starting state
  annotate("text",x = 1, y = 1,label = "START",
           fontface = "bold",size = 4,
           col='purple4') +
  # Food state
  annotate("text",x = 3, y = 3,label = "FOOD",
           fontface = "bold",size = 4,col='orange') +
  coord_fixed(xlim = c(0.6, 3.4),ylim = c(0.6, 3.4),
              clip = "off") +
  theme_void()
maze
############################################################

############################################################
# Transition matrix
p11 = c(0,.5,0,.5,0,0,0,0,0)
p21 = c(1/3,0,1/3,0,1/3,0,0,0,0)
p31 = c(0,.5,0,0,0,.5,0,0,0)
p12 = c(1/3,0,0,0,1/3,0,1/3,0,0)
p22 = c(0,.25,0,.25,0,.25,0,.25,0)
p32 = c(0,0,1/3,0,1/3,0,0,0,1/3)
p13 = c(0,0,0,.5,0,0,0,.5,0)
p23 = c(0,0,0,0,1/3,0,1/3,0,1/3)
p33 = c(0,0,0,0,0,0,0,0,1)

P = rbind(p11,p21,p31,p12,p22,p32,p13,p23,p33)
############################################################

############################################################
# Simulate one path

# Length of the chain
m = 30

# The chain starts at the lower left position
path = numeric(m+1)
path[1] = 1

# Sequential sampling
set.seed(314159)
for(i in 2:(m+1)){
  newPos = sample(x=c(1:9),size =1, prob = P[path[i-1],])
  path[i] = newPos
}

# Plots
df_path = data.frame(m=c(0:m),Xm = path)
ggplot(data = df_path,aes(x=m,y=Xm))+
  geom_line()+
  geom_point()+
  scale_y_continuous(breaks=1:9)+
  labs(x=expression(m),y=expression(X[m]))+
  theme_minimal()

# Animation of the Markov Chain in the labyrinth
# Current position at each frame
current_position = data.frame(
  frame = 0:m,
  nodes[path,]
)

# Accumulated path at each frame
accumulated_path = lapply(0:m, function(k){
  data.frame(
    frame = k,
    nodes[path[1:(k+1)],]
  )
}) |>
  bind_rows()

anim = maze +
  # Path followed so far
  geom_path(data = accumulated_path,
            aes(x = x, y = y, group = frame),
            colour = "darkred",
    linewidth = 1) +
  # Current position of the chain
  geom_point(data = current_position,aes(x = x, y = y),
    colour = "darkred",size = 5) +
  transition_manual(frame) +
  labs(title = "Step: {current_frame}")

animation = animate(anim,nframes = m + 1,fps = 1,
                   width = 1000,height = 1000,
                   renderer = gifski_renderer())

# Save as GIF
anim_save("Maze.gif",animation = animation)
