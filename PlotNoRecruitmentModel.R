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
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
    dB <- - Pb* (beta*Cb*B_cells + beta_2*Ct*B_cells)   
    dCb <-  Pb* (beta*Cb*B_cells + beta_2*Ct*B_cells) - alpha*Cb 
    dT <- -nu_a*Cb*T_cells - nu_a*Ct*T_cells  
    dAt <- nu_a*Cb*T_cells + nu_a*Ct*T_cells - beta_2*Ct*At - beta*Cb*At 
    Inf_At <- T_sus * (beta_2*Ct*At + beta*Cb*At)
    dLt <- (1-theta)*Inf_At  
    dCt <- (theta)*Inf_At - alpha_Ct*Ct 
    dF <- -nu_f*Lt*F
    dIf <- nu_f*Lt*F
    
    return(list(c(dB, dCb, dT, dAt, dLt, dCt, dF, dIf))) 
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
  df_for_ffe <- df %>% select(time,If)  
  
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
  
  #plot for infected feather follicles  
  FFE_plot <- ggplot(data = df_for_ffe, aes(x = time/24, y = If)) + geom_line( color = "pink") + 
    geom_point(data = baigent2016, aes(x = time, y = mean_genomes), inherit.aes = FALSE, color = "red") + 
    labs(title = "WithinHost Delay", color = "Cell Type")  + xlab(label = "Time (Days)") + 
    ylab(label = "Cell Number")  +  scale_y_log10(limits = c(0.01, 10000000)) + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),   legend.text = element_text(size = 12),
          legend.title = element_text(size = 12), axis.line = element_line(color = "black")) + labs( y = "Cell Number", x = "Time (Days)", title = "Infected Feather Follicle Epithelium") 
  
  
  print(everything)
  print(cytolytic_plot) 
  print(FFE_plot)
}

## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c(
  beta      = 6.951463e-07, 
  beta_2    = 6.951463e-07,
  nu_a      = 4.668718e-01, 
  nu_f      = 0.008,
  alpha     = 0.0104, 
  alpha_Ct  = 0.0104,
  theta     = 0.8,
  T_sus     = 0.5,
  Pb        = 0.03
) 


initial_values<- c(
  B_cells = 2.4e9/3,
  Cb      = 1,
  T_cells = 7.4e8/3,
  At      = 0,
  Lt      = 0,
  Ct      = 0, 
  F       = 400000, 
  If      = 0 
)

time_values <- seq(0, 1080) # hours


### DATA FOR PLOT ###
baigent1998 <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", 3 ) %>% 
  mutate(mean.pp38 = as.numeric(mean.pp38)) 

baigent2016 <- read_xlsx("/Users/rayanhg/Downloads/baigent2016.xlsx", 2 ) %>% arrange(time)


optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Jan.31.26.NoRecruitmentModel_SuspectibleTcells_SusceptibleBcells.csv") %>% 
  filter(Converged == 0) %>% 
  slice_min(Likely, n = 10) %>% 
  select(!c(Likely, Converged,lambda))%>% 
  mutate(across(everything(), ~as.numeric(.)))


# creating a list of dataframes
list_of_df <- purrr::map(seq_len(nrow(optim_data)), function(i) { 
  parameter_vector(dat = optim_data, i = i)
})

pdf("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Jan.31.26.NoRecruitmentModel_SuspectibleTcells_SusceptibleBcells.pdf", width = 7, height = 5)

generating_plots <- purrr::map(seq_len(nrow(optim_data)), function(i) {
  creating_plots(listofdf = list_of_df, i = i)
})

dev.off()

df <- list_of_df[[1]]
