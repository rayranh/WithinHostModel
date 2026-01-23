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
    dB <- -M*B_cells - beta*Cb*B_cells - beta_2*Ct*B_cells + (g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct 
    dCb <- M*B_cells +beta*Cb*B_cells + beta_2*Ct*B_cells - alpha*Cb 
    dT <- -M*T_cells - nu_a*Cb*T_cells - nu_b*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct)))
    dAt <- M*T_cells + nu_a*Cb*T_cells + nu_b*Ct*T_cells - beta_2*Ct*At - beta*Cb*At  # beta = rate of activated T cells by Ct and Cb cells    
    dLt <- theta*(beta_2*Ct*At + beta*Cb*At)  - lambda*Lt 
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - mu*Lt5 
    dCt <- (1-theta)*(beta_2*Ct*At + beta*Cb*At) - alpha_2*Ct
    dZ <-  mu*Lt5 
    df <- -nu_f*Lt5*f
    dIf <- nu_f*Lt5*f
    return(list(c(dB, dCb, dT, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dCt, dZ, df,dIf)))
  }) 
} 


parameter_vector <- function(dat,i) { 
  param_vector <- unlist(dat[i,]) 
  pars <- parameters_values
  pars[names(param_vector)] <- abs(param_vector)
  
  # run model
  results <- as.data.frame(
    ode(y = initial_values, 
        times = time_values, 
        func = sir_equations, 
        parms = pars))  
  
  return(results)
} 

creating_plots <- function(listofdf, i) { 
  ####################################
  ## Did you put LT in denominator? ## 
  ###################################
  
  #take out a dataframe from each list and plot in ggplot/ using pivot longer because melt is old 
  df <- listofdf[[i]] # kept cytolytic infection  
  df_for_40000 <- df %>% mutate(cytolytic_scale_40000 = ((Cb+Ct)/(B_cells+Cb+Ct+T_cells+At))* 40000) # for every 40000 cells in my model how many infected 
  df_for_ffe <- df %>% select(time,If)
  
  #plotting cytolytic data 
  cytolytic_plot <- ggplot(df_for_40000, aes(x = time/24, y = cytolytic_scale_40000)) + 
    geom_line() + geom_line(data = baigent1998_summary,aes(x = time / 24, y = mean_pp38),color = "red",linewidth = 1, inherit.aes = FALSE) +  
    labs(x = "Time (days)", y = "pp38+ cells (out of 40,000)",title = "Model vs observed pp38+ cells (average per bird)") + theme_minimal() + ylim(0,200)
  
  #plot for infected feather follicles  
  FFE_plot <- ggplot(data = df_for_ffe, aes(x = time/24, y = If)) + geom_line( color = "pink") + 
    geom_point(data = baigent2016, aes(x = time, y = mean_genomes), inherit.aes = FALSE, color = "red") + 
    labs(title = "WithinHost Delay", color = "Cell Type")  + xlab(label = "Time (Days)") + 
    ylab(label = "Cell Number")  +  scale_y_log10(limits = c(0.01, 10000000)) + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),   legend.text = element_text(size = 12),
          legend.title = element_text(size = 12), axis.line = element_line(color = "black")) + labs( y = "Cell Number", x = "Time (Days)", title = "WithinHost - Blood (B and T total)") 
  
  print(cytolytic_plot) 
  print(FFE_plot)
}


# take out a named vector that is in dataframe iteratively 

## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  M = 0
  , beta =   4.214202                 #contact rate with B cells 
  , beta_2 = 62.1114629                   #contact rate with T cells 
  , nu_a = 4.668718e-01                     #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = 4.467098e-01                   #Activation rate of T cells by cytolytic T cells (hours)
  , nu_f = 6.033255e-02                      #Infection rate of follicular cells (hours)
  , mu =  0.03461908                         #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.0104                   #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.0104                 #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 929.11676                     #incoming B cells (every 15 hours)
  , g2 = 29.86632  
  , h1 = 289.4683                   #incoming T cells / determined no incoming T cells 
  , h2 = 699.0949 
  , lambda = 15.975984e-05                 #adding delay, how long latent cell 'exposed' 
)

initial_values <- c( 
  B_cells = 2.4e6/3  # from three organs 2.4e6/3 
  , Cb = 1 
  , T_cells = 1.5e6/3 # from three organs 1.5e6/3 
  , At = 0 
  , Lt = 0 
  , Lt2 = 0 
  , Lt3 = 0  
  , Lt4 = 0  
  , Lt5 = 0 
  , Ct = 0
  , Z = 0
  , f = 400000 
  , If =0 
  
) 

time_values <- seq(0, 1000) # hours  


### DATA FOR PLOT ### 
baigent2016 <- read_xlsx("/Users/rayanhg/Downloads/baigent2016.xlsx", 2 ) 
baigent1998 <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", 3 )    
optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Random_parameter_exploration_2.csv") %>% 
  filter(Converged == 0) %>% slice_min(Likely, n = 10) %>% select(!c(Likely, Converged, X)) 

baigent1998_summary <- baigent1998 %>% 
  mutate(mean.pp38 = as.numeric(mean.pp38)) %>% 
  group_by(time) %>%
  summarise(
    mean_pp38 = mean(mean.pp38, na.rm = TRUE),
    se_pp38   = sd(mean.pp38, na.rm = TRUE) / sqrt(sum(!is.na(mean.pp38))),
    n_birds   = sum(!is.na(mean.pp38)),
    .groups = "drop"
  )

#creating a list of dataframes set in motion by this function made on the fly that takes optim data and i takes in seq_len(nrow(()))
list_of_df <-  purrr::map(seq_len(nrow(optim_data)), function(i) { 
  parameter_vector(dat = optim_data, i = i)
  }) 

pdf("model_plots_csv_2.pdf", width = 7, height = 5)

generating_plots <-  purrr::map(seq_len(nrow(optim_data)), function(i) { 
  creating_plots(listofdf = list_of_df, i = i)
}) 

dev.off()



