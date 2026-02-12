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
    dB <-  -beta*Cb*B_cells - beta_2*Ct*B_cells + Pb*(g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct  
    dBr <- (1-Pb)*(g1*(Cb+Ct)/(g2+(Cb+Ct))) #probably should have alpha for cytolytic infection 
    dCb <- beta*Cb*B_cells + beta_2*Ct*B_cells - alpha*Cb 
    dT <- -nu_a*Cb*T_cells - nu_b*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct)))
    dAt <- nu_a*Cb*T_cells + nu_b*Ct*T_cells - beta_2*Ct*At - beta*Cb*At  # beta = rate of activated T cells by Ct and Cb cells   
    dLt <- theta *(beta_2*Ct*At + beta*Cb*At)  - lambda*Lt 
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - mu*Lt5
    dCt <-  (1-theta)*(beta_2*Ct*At + beta*Cb*At) - alpha_2*Ct
    dZ <- mu*Lt5
    df <- -nu_f*Lt5*f
    dIf <- nu_f*Lt5*f
    return(list(c(dB, dCb,dBr, dT, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dCt, dZ, df,dIf)))
  }) 
} 


parameter_vector <- function(dat,i) { 
  param_vector <- unlist(dat[i,])  
  init <- initial_values 
  B0 <-  2.4e9/3 
  
  pars <- parameters_values
  pars[names(param_vector)] <- param_vector 
  init["B_cells"] <- B0*pars["Pb"] 
  init["Br"] <- B0*(1 - pars["Pb"])  
  
  
  
  
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
  
  #take out a dataframe from each list and plot in ggplot/ using pivot longer because melt is old 
  df <- listofdf[[i]] # kept cytolytic infection    
  
  #turning df into long df 
  df2 <- pivot_longer(df, cols = -time, names_to = "variable", values_to = "value")
  
  df_for_40000 <- df %>% mutate(cytolytic_scale_40000 = ((Cb+Ct)/(B_cells+Cb+Ct+T_cells+At+Br+Lt+Lt2+Lt3+Lt4+Lt5))* 40000) 
    
  df_for_ffe <- df %>% select(time,If)  
  
  # #plot for all components 
  everything <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("B_cells" = "black", "T_cells" = "red", "Cb"="green", "At" = "blue", "Lt5" = "purple",
                                  "Ct" = "yellow", "Z" = "lightblue", "Br" = "orange", "f" = "magenta", "If" = "pink",
                                  "Lt2" = "chartreuse", "Lt3" = "chartreuse", "Lt4" = "chartreuse")) +
    labs(title = "WithinHost Delay", color = "Cell Type") + theme(legend.position = "right") +
    theme_minimal() + xlab(label = "Time (Days)") + ylab(label = "Cell Number") +
    geom_point(data = baigent1998, aes(x = time/24, y = mean.pp38), inherit.aes = FALSE, color = "red") 

  #plotting cytolytic data 
  cytolytic_plot <- ggplot(df_for_40000, aes(x = time/24)) + 
    geom_line(aes(y = cytolytic_scale_40000, color = "Cb+Ct")) +
    scale_color_manual(
      name = "cell type",
      values = c("Cb" = "green", "Ct" = "red")
    ) +
    geom_point(data = baigent1998, aes(x = time/24, y = mean.pp38), color = "red", size = 0.3, inherit.aes = FALSE) +
    # geom_point(data = baigent1998, aes(x = time/24, y = Bcell_no), color = "green",   size = 0.3, inherit.aes = FALSE) +
    labs(x = "Time (days)", y = "pp38+ cells (out of 40,000)",
         title = "Model vs observed pp38+ cells (average per bird)") +
    theme_minimal() + ylim(0, 1000)
  
  
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



# take out a named vector that is in dataframe iteratively 

## PARAMETERS AND INITIAL VALUES ## 

parameters_values <- c( 
  beta =  6.951463e-07                  #contact rate with B cells  
  , Pb = 0.005
  , beta_2 = 8.532282e-07                  #contact rate with T cells 
  , nu_a = 4.668718e-01                     #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = 4.467098e-01                   #Activation rate of T cells by cytolytic T cells (hours)
  , nu_f = 0.008                     #Infection rate of follicular cells (hours)
  , mu = 0.02                        #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.0104                   #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.0104                 #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 10                     #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 =10                      #incoming T cells / determined no incoming T cells 
  , h2 = 100
  , lambda =0.02380952                  #adding delay, how long latent cell 'exposed' 
)

initial_values <- c( 
  B_cells = 2.4e9/3   # from three organs 2.4e6/3 
  , Cb = 1  
  , Br =2.4e9/3 
  , T_cells =7.4e8/3# from three organs 1.5e6/3 
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


time_values <- seq(0, 1080) # hours  


### DATA FOR PLOT ### 
baigent2016 <- read_xlsx("baigent2016.xlsx", 2 ) 
baigent1998 <- read_xlsx("baigent1998.xlsx", 3 ) %>% 
  mutate(mean.pp38 = as.numeric(mean.pp38),
         Bcell_no = as.numeric(Bcell_no), 
         Tcell_no = as.numeric(Tcell_no)) %>% filter(!is.na(mean.pp38))   
optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Feb.12.26.weightedLikelihood.csv") %>% 
  filter(Converged == 0) %>% slice_min(Likely, n = 10) %>% select(c(beta, beta_2, alpha, alpha_2,nu_a,nu_b,nu_f,mu,g1,g2,h1,h2,Pb)) 


#creating a list of dataframes set in motion by this function made on the fly that takes optim data and i takes in seq_len(nrow(()))
list_of_df <-  purrr::map(seq_len(nrow(optim_data)), function(i) { 
  parameter_vector(dat = optim_data, i = i)
  }) 

 pdf("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Feb.12.26.weightedLikelihood.pdf", width = 7, height = 5)

 generating_plots <-  purrr::map(seq_len(nrow(optim_data)), function(i) {
 creating_plots(listofdf = list_of_df, i = i)
})

dev.off()

# df <- list_of_df[[1]] 
# df2 <- pivot_longer(df, cols = -time, names_to = "variable", values_to = "value")
# df_for_40000 <- df %>% mutate(cytolytic_scale_40000 = ((Cb+Ct)/(B_cells+Cb+Ct+T_cells+At+Br+Lt+Lt2+Lt3+Lt4+Lt5))* 40000) # for every 40000 cells in my model how many infected



# 
# everything <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
#   scale_color_manual(values = c("B_cells" = "black", "T_cells" = "red", "Cb"="green", "At" = "blue", "Lt5" = "purple",
#                                 "Ct" = "yellow", "Z" = "lightblue", "Br" = "orange", "f" = "magenta", "If" = "pink",
#                                 "Lt2" = "chartreuse", "Lt3" = "chartreuse", "Lt4" = "chartreuse")) +
#   labs(title = "WithinHost Delay", color = "Cell Type") + theme(legend.position = "right") +
#   theme_minimal() + xlab(label = "Time (Days)") + ylab(label = "Cell Number") +
#   geom_point(data = baigent1998, aes(x = time/24, y = mean.pp38), inherit.aes = FALSE, color = "red")
# 
#plotting cytolytic data
# cytolytic_plot <- ggplot(df_for_40000, aes(x = time/24, y = cytolytic_scale_40000)) +
#   geom_line(color = "#A6D854") + geom_point(data = baigent1998,aes(x = time / 24, y = mean.pp38),color = "red",linewidth = 1, inherit.aes = FALSE) +
#   labs(x = "Time (days)", y = "pp38+ cells (out of 40,000)",title = "Model vs observed pp38+ cells (average per bird)") + theme_minimal() + ylim(0,1000)
# 
  #plot for infected feather follicles
# FFE_plot <- ggplot(data = df_for_ffe, aes(x = time/24, y = If)) + geom_line( color = "pink") +
#   geom_point(data = baigent2016, aes(x = time, y = mean_genomes), inherit.aes = FALSE, color = "red") +
#   labs(title = "WithinHost Delay", color = "Cell Type")  + xlab(label = "Time (Days)") +
#   ylab(label = "Cell Number")  +  scale_y_log10(limits = c(0.01, 10000000)) +
#   theme(panel.grid = element_blank(), panel.background = element_blank(),   legend.text = element_text(size = 12),
#         legend.title = element_text(size = 12), axis.line = element_line(color = "black")) + labs( y = "Cell Number", x = "Time (Days)", title = "Infected Feather Follicle Epithelium")
