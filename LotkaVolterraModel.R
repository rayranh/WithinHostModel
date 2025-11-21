library(deSolve)  
library(reshape2) 
library(ggplot2)
library(dplyr)


lotka_volterra <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations 

dx <- x * (alpha - beta*y) # prey 
dy <- -y*(gamma - delta*x) # predator 

    return(list(c(dx, dy)))
  }) 
}



initial_values <- c(
  x = 10 
  , y = 10 
) 


parameter_values <- c( 
  alpha = 2 # growth rate of prey 
  ,beta = 2 # rate at which prey eaten by predator 
  ,gamma = 0.5 # predator dies off at this rate 
  ,delta = 0.7)  # predator gains vitality 


time <- seq(0,100, by = 1)

out <- as.data.frame(ode(y = initial_values, 
                         times = time, 
                         func = lotka_volterra, 
                         parms = parameter_values)) 


ggplot( out, mapping = aes(time)) + 
  geom_line(aes(y = x, colour = "prey")) + 
  geom_line(aes(y = y, colour = "predator")) +  
  scale_color_manual(labels = c('prey', 'predator'), values = c("blue", "red")) + 
  labs(title = "Lotka Volterra", y = "no. individuals", x  = "time")  





