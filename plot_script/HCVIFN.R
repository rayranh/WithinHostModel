rm(list = ls())
library(deSolve) 
library(ggplot2)
library(dplyr) 
library(purrr)



HCVModel <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
    # dT <- s - dT - (1 - n)*beta*V*T 
    dI <- (1-n)*beta*V*T_const - delta*I 
    dV <- (1-e)*p*I - c*V 
    return(list(c(dI,dV)))
  }) 
}  

time_values <- seq(0,14,by=1) 

n_values <- c(0.80,0.95,0.99,1)  



initial_values <- c(  
  I = 0.05, 
  V = 1 # they are plotting relative relationships not absolute viral loads so can choose any number  
)

base_values <- c(
  beta = 3e-7, 
  delta = 0.5, 
  p = 100, 
  c = 5,
  e = 0 
)  

T_const <- (base_values["c"] * base_values["delta"])/
  (base_values["p"] * base_values["beta"])  



results <- map_dfr(n_values, function (n_val) {    
  
  parameters_values <- c(
    beta = 3e-7, 
    delta = 0.5, 
    p = 100, 
    c = 5,
    n = n_val, 
    e = 0, 
    T = T_const
  )  

  # run model
  results <- as.data.frame(
    ode(y = initial_values, 
        times = time_values, 
        func = HCVModel , 
        parms = parameters_values)
  )  
  # 
  results %>% mutate(parameter = n_val, relative_decrease = log10(V/1) ) #using n_val because n_val only takes one number at a time! 
  # 
  }) 


ggplot(data = results, aes(x = time, y = relative_decrease, colour = factor(parameter))) + geom_line() +
  coord_cartesian(ylim = c(-5,0), xlim = c(0,14)) + 
  scale_x_continuous(breaks = seq(0,14, by = 2), expand = c(0,0))  #expand gets rid of extra padding before 0 and after 14 

