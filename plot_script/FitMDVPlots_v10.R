###plotting the top 10 best likelihood scores in ggplot from csv file ###  
### Make sure same starting initial conditions ### 

rm(list = ls())
library(deSolve)  
library(tidyr) 
library(ggplot2)
library(dplyr) 
library(readxl)  
library(purrr) 

setwd("~/Desktop/WithinHostModel/DataForProject/")

sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
    dB <-  -beta*Cb*B_cells - beta_2*Ct*B_cells - beta_2*Lr*B_cells + Pb*(g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct  
    dBr <- (1-Pb)*(g1*(Cb+Ct+Lr)/(g2+(Cb+Ct+Lr))) #probably should have alpha for cytolytic infection 
    dCb <- beta*Cb*B_cells + beta_2*Ct*B_cells + beta_2*Lr*B_cells - alpha*Cb 
    dT <- -nu_a*Cb*T_cells - nu_b*Ct*T_cells - nu_b*Lr*T_cells + Pt*(h1*(Cb+Ct)/(h2+(Cb+Ct))) 
    dTr <- (1 - Pt)*(h1*(Cb+Ct)/(h2+(Cb+Ct))) 
    dAt <- nu_a*Cb*T_cells + nu_b*Ct*T_cells + nu_b*Lr*T_cells - beta_2*Ct*At - beta*Cb*At - beta_2*Lr*At  # beta = rate of activated T cells by Ct and Cb cells   
    dLt <- theta *(beta_2*Ct*At + beta*Cb*At + beta_2*Lr*At)  - lambda*Lt 
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - mu*Lt5 - kappa*Lt5  
    dLr <-  kappa*Lt5 - alpha_2*Lr 
    dCt <-  (1-theta)*(beta_2*Ct*At + beta*Cb*At+ beta_2*Lr*At) - alpha_2*Ct
    dZ <- mu*Lt5
    df <- -nu_f*Lr*f
    dIf <- nu_f*Lr*f
    return(list(c(dB, dCb,dBr, dT,dTr, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dLr, dCt, dZ, df,dIf)))
  }) 
}  

# sir_equations <- function(time, variables, parameters) {
#   with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
#     dB <-  -beta*Cb*B_cells - beta_2*Ct*B_cells - beta_2*Lr*B_cells + Pb*(g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct  
#     dBr <- (1-Pb)*(g1*(Cb+Ct+Lr)/(g2+(Cb+Ct+Lr))) #probably should have alpha for cytolytic infection 
#     dCb <- beta*Cb*B_cells + beta_2*Ct*B_cells + beta_2*Lr*B_cells - alpha*Cb 
#     dT <- -nu_a*Cb*T_cells - nu_b*Ct*T_cells - nu_b*Lr*T_cells + Pt*(h1*(Cb+Ct)/(h2+(Cb+Ct))) 
#     dTr <- (1 - Pt)*(h1*(Cb+Ct)/(h2+(Cb+Ct))) 
#     dAt <- nu_a*Cb*T_cells + nu_b*Ct*T_cells + nu_b*Lr*T_cells - beta_2*Ct*At - beta*Cb*At - beta_2*Lr*At  # beta = rate of activated T cells by Ct and Cb cells   
#     dLt <- theta *(beta_2*Ct*At + beta*Cb*At + beta_2*Lr*At)  - lambda*Lt 
#     dLt2 <- lambda*Lt - lambda*Lt2 
#     dLt3 <- lambda*Lt2 - lambda*Lt3 
#     dLt4 <- lambda*Lt3 - lambda*Lt4
#     dLt5 <- lambda*Lt4 - mu*Lt5 - kappa*Lt5  
#     dLr <-  kappa*Lt5 - alpha_2*Lr 
#     dCt <-  (1-theta)*(beta_2*Ct*At + beta*Cb*At+ beta_2*Lr*At) - alpha_2*Ct
#     dZ <- mu*Lt5
#     df <- -nu_f*Lr*f
#     dIf <- nu_f*Lr*f
#     return(list(c(dB, dCb,dBr, dT,dTr, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dLr, dCt, dZ, df,dIf)))
#   }) 
# }  

parameter_vector <- function(dat,i) { 
  param_vector <- unlist(dat[i,])  
  init <- initial_values 
  B0 <-  2.4e9/3 
  T0 <- 7.4e8 / 3
  
  
  pars <- parameters_values
  pars[names(param_vector)] <- param_vector 
  init["B_cells"] <- B0*pars["Pb"] 
  init["Br"] <- B0*(1 - pars["Pb"])  
  
  
  init["T_cells"] <- T0*pars["Pt"] 
  init["Tr"] <- T0* (1 - pars["Pt"])
  
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
  
  df_for_40000 <- df %>% mutate(
                                cytolytic_scale_40000_cb = ((Cb)/(B_cells+Cb+Ct+T_cells+At+Br+Lt+Lt2+Lt3+Lt4+Lt5+Lr+Tr))* 40000, 
                                cytolytic_scale_40000_ct = ((Ct+Lr)/(B_cells+Cb+Ct+T_cells+At+Br+Lt+Lt2+Lt3+Lt4+Lt5+Lr+Tr))* 40000)  
  
  df_pfu <- df %>% mutate(pfuCells = (Cb + Ct + Lr)/(B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5 + Lr + Tr) * 10000000) 
  
  df_cyto_T <- df2 %>% filter(variable %in% c("T_cells", "At", "Lt5", "Ct", "Lr")) 
  
  df_cyto_B <- df2 %>% filter(variable %in% c("B_cells", "Cb", "At", "Lr")) 
  

  # df_PBL <- df %>% mutate(PBL_scale_10000 = (Cb+Ct+Lt+Lt2+Lt3+Lt4+Lt5)/(B_cells+T_cells+Br+Cb+Ct+At+Lt+Lt2+Lt3+Lt4+Lt5)*10000)
  
  # #plot for all components 
  everything <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("T_cells" = "red", "At" = "blue", "Lt5" = "purple", "B_cells" = "black", 
                                  "Z" = "lightblue", "Br" = "orange","Lt2" = "chartreuse", 
                                  "Lt3" = "chartreuse", "Lt4" = "chartreuse", "Lt" = "chartreuse", "Lr" = "aquamarine", 
                                  "f" = "pink","If" = "magenta", "Ct" = "brown", "Cb" = "darkorange2"), labels = c("T_cells" = "T cells")) +
    labs(title = "Within-Host Model", color = "Cell Type") + theme(legend.position = "right") + xlab(label = "Time (days post infection)") + 
    ylab(label = "Cell Number") + theme_classic() + scale_x_continuous(limits = c(0,60), breaks = seq(0,60, by = 5))
  
  
  everything_log <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("T_cells" = "red", "At" = "blue", "Lt5" = "purple", "B_cells" = "black", 
                                  "Z" = "lightblue", "Br" = "orange","Lt2" = "chartreuse", 
                                  "Lt3" = "chartreuse", "Lt4" = "chartreuse", "Lt" = "chartreuse", "Lr" = "aquamarine", 
                                  "f" = "pink","If" = "magenta", "Ct" = "brown", "Cb" = "darkorange2"), labels = c("T_cells" = "T cells")) +
    labs(title = "Within-Host Model", color = "Cell Type") + theme(legend.position = "right") + xlab(label = "Time (days post infection)") + 
    ylab(label = "Cell Number") + theme_classic() + scale_y_log10() + scale_x_continuous(limits = c(0,60), breaks = seq(0,60, by = 5)) 
  
  
  log_cyto_T <- ggplot(data = df_cyto_T, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("T_cells" = "red", "At" = "blue", "Lt5" = "purple","Ct" = "brown", "Lr" = "aquamarine"), labels = c("T_cells" = "T cells")) +
    labs(title = "Within-Host Model - log Cyto T", color = "Cell Type") + theme(legend.position = "right") + xlab(label = "Time (days post infection)") + 
    ylab(label = "Cell Number") + theme_classic() + scale_y_log10() + scale_x_continuous(limits = c(0,60), breaks = seq(0,60, by = 5))
  
  
  log_cyto_B <- ggplot(data = df_cyto_B, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("B_cells" = "black", "Cb" = "darkorange2",  "At" = "blue", "Lr" = "aquamarine"), labels = c("T_cells" = "T cells")) +
    labs(title = "Within-Host Model - log Cyto B", color = "Cell Type") + theme(legend.position = "right") + xlab(label = "Time (days post infection)") + 
    ylab(label = "Cell Number") + theme_classic() + scale_y_log10() + scale_x_continuous(limits = c(0,60), breaks = seq(0,60, by = 5))
  
  
  #plotting cytolytic data 
  cytolytic_plot <- ggplot(df_for_40000, aes(x = time/24)) +
    geom_line(aes(y = cytolytic_scale_40000_cb, color = "Cb")) + 
    geom_line(aes(y = cytolytic_scale_40000_ct, color = "Ct")) + 
    scale_color_manual(name = "cell type",values = c("Ct" = "red", "Cb" = "green")) +
    geom_point(data = baigent1998, aes(x = time/24, y = Bcell_no), color = "green",   size = 0.3, inherit.aes = FALSE) + 
    geom_point(data = baigent1998, aes(x = time/24, y = Tcell_no), color = "red",   size = 0.3, inherit.aes = FALSE) +
    labs(x = "Time (days post infection)", y = "pp38 Cell Counts per 40,000 cells", title = "Cytolytic Cells") + theme_classic()
  
  pfu_plot <- ggplot(df_pfu, aes(time/24)) + 
    geom_line(aes(y = pfuCells, color = "pfu")) + 
    scale_color_manual(name = "cell type", values = c("pfu" = "pink")) + 
    geom_point(data = powell1982_plaque, aes(x = time/24, y = pfuPer1e7)) +
    geom_errorbar(data = powell1982_plaque, aes(ymin = SE_upper, ymax = SE_lower)) +  
    geom_line(data = powell1982_plaque, aes(x = time/24, y = pfuPer1e7)) + scale_y_log10() + 
    labs(x = "Time (days post infection)", y = "PFU/10^7", title = "PFU Cells") + theme_classic()
  
  print(everything)
  print(everything_log)  
  print(log_cyto_T) 
  print(log_cyto_B)
  print(pfu_plot) 
  print(cytolytic_plot) 

  
}



# take out a named vector that is in dataframe iteratively 

## PARAMETERS AND INITIAL VALUES ## 

## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  beta =  6.951463e-07                  #contact rate with B cells   
  ,beta_2 =  6.951463e-07                  #contact rate with B cells
  , Pb = 0.01 
  , Pt = 0.01
  , nu_a = 0.1666667            #Activation rate of T cells by cytolytic B cells (hours) 
  , nu_b = 0.4467098
  , nu_f = 0.008                    #Infection rate of follicular cells (hours) 
  , kappa = 0.0001
  , mu = 1/8                       #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.0104          #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.0104          #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 10                    #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 = 10                      #incoming T cells / determined no incoming T cells 
  , h2 = 100
  , lambda = 0.02380952             # fixing delay rate to (1/(7*24))*4   
  , size_pp38 = 10  # new parameter 
  , size_pp38_Ct = 10
)


initial_values <- c( 
  B_cells = 2.4e9/3    # from three organs 2.4e6/3 
  , Cb = 1  
  , Br = 2.4e9/3       #same initial proportion here infected as in model (well mixed marbles ex) 
  , T_cells = 7.4e8/3  # from three organs 1.5e6/3  
  , Tr = 0 
  , At = 0 
  , Lt = 0 
  , Lt2 = 0 
  , Lt3 = 0  
  , Lt4 = 0  
  , Lt5 = 0  
  , Lr = 0 
  , Ct = 0
  , Z = 0
  , f = 400000 
  , If =0 
  
) 

time_values <- seq(0, 1440, by = 1)   


### DATA FOR PLOT ### 

baigent1998 <- read_xlsx("Baigent1998/baigent1998.xlsx", sheet = 3, na = "NA") %>% filter(!is.na(mean.pp38)) %>% select(time, Bcell_no,Tcell_no) 

powell1982_plaque <- read_xlsx("Powell1982/Powell1982.xlsx", sheet = 2)%>% mutate(time = dpi*24) 

optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/FitMDV_v10_HPRS16_20260720_111823.csv") %>% 
  filter(Converged == 0) %>% slice_min(Likely, n = 10) %>% arrange(Likely) %>% 
  select(c(beta, alpha, alpha_2,nu_a,nu_f,mu,Pb, size_pp38, size_pp38_Ct, g1,g2, h1,h2, beta_2,nu_b, kappa, Pt)) 


#creating a list of dataframes set in motion by this function made on the fly that takes optim data and i takes in seq_len(nrow(()))
list_of_df <-  purrr::map(seq_len(nrow(optim_data)), function(i) { 
  parameter_vector(dat = optim_data, i = i)
}) 

pdf("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/FitMDV_v10_HPRS16_20260720_111823.pdf", width = 7, height = 5)

generating_plots <-  purrr::map(seq_len(nrow(optim_data)), function(i) {
  creating_plots(listofdf = list_of_df, i = i)
})

dev.off()

