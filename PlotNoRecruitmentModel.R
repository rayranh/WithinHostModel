###plotting the top 10 best likelihood scores in ggplot from csv file ###  
### Make sure same starting initial conditions ### 

#rm(list = ls())
library(deSolve)  
library(tidyr) 
library(ggplot2)
library(dplyr) 
library(readxl)  
library(purrr)

sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dB  <- -beta*Cb*B_cells - beta*Ct*B_cells  
    dCb <-  beta*Cb*B_cells + beta*Ct*B_cells  - alpha*Cb 
    dT  <- -nu_a*Cb*T_cells - nu_a*Ct*T_cells  
    dAt <-  nu_a*Cb*T_cells + nu_a*Ct*T_cells - beta*Ct*At - beta*Cb*At 
    
    Inf_At <- T_sus * (beta*Ct*At + beta*Cb*At)
    dLt <- theta*Inf_At - lambda*Lt  
    dCt <- (1-theta)*Inf_At - alpha*Ct
    
    return(list(c(dB, dCb, dT, dAt, dLt, dCt))) 
  }) 
}  

parameter_vector <- function(dat,i) { 
  param_vector <- unlist(dat[i,])  
  init <- initial_values 
  
  pars <- parameters_values
  pars[names(param_vector)] <- param_vector 
  
  # run model
  results <- as.data.frame(
    ode(y = init, 
        times = time_values, 
        func = sir_equations, 
        parms = pars))  
  
  return(results)
} 

creating_plots <- function(listofdf, i) { 
  ####################################
  ## Did you put LT in denominator? ## 
  ###################################
  
  df <- listofdf[[i]]
  
  df2 <- pivot_longer(df, cols = -time, names_to = "variable", values_to = "value")
  df_for_40000 <- df %>% mutate(
    cytolytic_scale_40000 = ((Cb+Ct)/(B_cells + Cb + Ct + T_cells + At + Lt)) * 40000
  )
  
  # All components plot
  everything <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable)) + 
    geom_line() +
    scale_color_manual(values = c("B_cells"="black","Cb"="green","T_cells"="red","At"="blue","Lt"="purple","Ct"="yellow")) +
    labs(title = "WithinHost (Simplified)", color = "Cell Type") +
    theme_minimal() + xlab("Time (Days)") + ylab("Cell Number") +
    geom_point(data = baigent1998, aes(x = time/24, y = mean.pp38),
               inherit.aes = FALSE, color = "red") 
  
  # Cytolytic pp38 scale plot
  cytolytic_plot <- ggplot(df_for_40000, aes(x = time/24, y = cytolytic_scale_40000)) + 
    geom_line(color = "#A6D854") + 
    geom_point(data = baigent1998, aes(x = time/24, y = mean.pp38),
               color = "red", linewidth = 1, inherit.aes = FALSE) +  
    labs(x = "Time (days)", y = "pp38+ cells (out of 40,000)",
         title = "Model vs observed pp38+ cells (average per bird)") +
    theme_minimal() + ylim(0,1000)
  
  print(everything)
  print(cytolytic_plot)
}


## PARAMETERS AND INITIAL VALUES ##
## IMPORTANT: for plotting we keep these in NATURAL scale (since optim_data is natural scale too)

parameters_values <- c( 
  beta   = 6.951463e-07,
  nu_a   = 4.668718e-01,
  alpha  = 0.0104,
  theta  = 0.8,
  lambda = 0.02380952,
  T_sus  = 0.005
)

initial_values <- c(
  B_cells = 2.4e6/3,
  Cb      = 1,
  T_cells = 1.5e6/3,
  At      = 0,
  Lt      = 0,
  Ct      = 0
)

time_values <- seq(0, 1080) # hours


### DATA FOR PLOT ###
baigent1998 <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", 3 ) %>% 
  mutate(mean.pp38 = as.numeric(mean.pp38))

optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/jan_9_26_SuperSimpleModel.csv") %>% 
  filter(Converged == 0) %>% 
  slice_min(Likely, n = 10) %>% 
  select(!c(Likely, Converged, X))

# creating a list of dataframes
list_of_df <- purrr::map(seq_len(nrow(optim_data)), function(i) { 
  parameter_vector(dat = optim_data, i = i)
})

pdf("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/jan_9_26_SuperSimpleModel.pdf", width = 7, height = 5)

generating_plots <- purrr::map(seq_len(nrow(optim_data)), function(i) {
  creating_plots(listofdf = list_of_df, i = i)
})

dev.off()

df <- list_of_df[[1]]
