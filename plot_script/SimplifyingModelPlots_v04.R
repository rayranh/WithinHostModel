#v04 - Getting rid of Recruitment, Changed to lognormal distribution, making betas the same 

rm(list = ls())
library(deSolve)  
library(tidyr) 
library(ggplot2)
library(dplyr) 
library(readxl)  
library(purrr) 

setwd("~/Desktop/WithinHostModel/DataForProject/")

## SIR MODEL ## 
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
    dB <-  -beta*Cb*B_cells - beta*Ct*B_cells - beta*Lr*B_cells + Pb*(g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct  
    dBr <- (1-Pb)*(g1*(Cb+Ct+Lr)/(g2+(Cb+Ct+Lr))) #probably should have alpha for cytolytic infection 
    dCb <- beta*Cb*B_cells + beta*Ct*B_cells + beta*Lr*B_cells - alpha*Cb 
    dT <- -nu_a*Cb*T_cells - nu_b*Ct*T_cells - nu_b*Lr*T_cells + Pt*(h1*(Cb+Ct)/(h2+(Cb+Ct))) 
    dTr <- (1-Pt) * (h1*(Cb + Ct) / (h2 + Cb + Ct))
    dAt <- nu_a*Cb*T_cells + nu_b*Ct*T_cells + nu_b*Lr*T_cells - beta*Ct*At - beta*Cb*At - beta*Lr*At  # beta = rate of activated T cells by Ct and Cb cells   
    dLt <- theta *(beta*Ct*At + beta*Cb*At + beta*Lr*At)  - lambda*Lt 
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - mu*Lt5 - kappa*Lt5  
    dLr <-  kappa*Lt5 - alpha_2*Lr
    dCt <-  (1-theta)*(beta*Ct*At + beta*Cb*At+ beta*Lr*At) - alpha_2*Ct
    dZ <- mu*Lt5
    df <- -nu_f*Lr*f
    dIf <- nu_f*Lr*f
    return(list(c(dB, dCb,dBr, dT,dTr, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dLr, dCt, dZ, df,dIf)))
  }) 
}  


parameter_vector <- function(dat,i) { 
  param_vector <- unlist(dat[i,])  
  init <- initial_values 
  B0 <-  2.4e9/3  
  T0 <- 7.4e8/3
  
  
  pars <- parameters_values
  pars[names(param_vector)] <- param_vector 
  init["B_cells"] <- B0*pars["Pb"] 
  init["Br"] <- B0*(1 - pars["Pb"])  
  
  
  init["T_cells"] <- T0 * pars["Pt"]
  init["Tr"]       <- T0 * (1 - pars["Pt"])
  
  
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
  
  df_for_40000 <- df %>% mutate(cytolytic_scale_40000 = ((Cb+Ct+Lr)/(B_cells+Cb+Ct+T_cells+At+Br+Lt+Lt2+Lt3+Lt4+Lt5+Lr+ Tr))* 40000, 
                                cytolytic_scale_40000_cb = ((Cb)/(B_cells+Cb+Ct+T_cells+At+Br+Lt+Lt2+Lt3+Lt4+Lt5+Lr+Tr))* 40000, 
                                cytolytic_scale_40000_ct = ((Ct+Lr)/(B_cells+Cb+Ct+T_cells+At+Br+Lt+Lt2+Lt3+Lt4+Lt5+Lr+Tr))* 40000) 
  
  df_for_ffe <- df %>% dplyr::select(time,If, f) %>% mutate(IfPer10k = (If)/(If+f)* 10000 * optim_data$q_FFE[i]) 
  
  df_for_PBL <- df %>% mutate(PBL_scale_10000 = (Cb + Ct + Lt + Lt2 + Lt3 + Lt4 + Lt5 + Lr)/(B_cells+Cb+Ct+T_cells+At+Br+Lt+Lt2+Lt3+Lt4+Lt5+Lr+Tr) * 10000)
  
  df_Pfu <- df %>% mutate(Pfu_scale_1000000 = (Cb+Ct+Lr)/(B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5 + Lr+Tr) * 1000000)
  
  df_cyto_T <- df2 %>% filter(variable %in% c("T_cells", "At", "Lt5", "Ct", "Lr")) 
  
  df_cyto_B <- df2 %>% filter(variable %in% c("B_cells", "Cb", "At", "Lr")) 
  
  
  
  # #plot for all components 
  everything <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("T_cells" = "red", "At" = "blue", "Lt5" = "purple", "B_cells" = "black", 
                                  "Z" = "lightblue", "Br" = "orange","Lt2" = "chartreuse", 
                                  "Lt3" = "chartreuse", "Lt4" = "chartreuse", "Lt" = "chartreuse", "Lr" = "aquamarine", 
                                  "f" = "pink","If" = "magenta", "Ct" = "brown", "Cb" = "darkorange2", "Tr" = "goldenrod3"), labels = c("T_cells" = "T cells")) +
    labs(title = "Within-Host Model", color = "Cell Type") + theme(legend.position = "right") + xlab(label = "Time (days post infection)") + 
    ylab(label = "Cell Number") + theme_classic() + scale_x_continuous(limits = c(0,40), breaks = seq(0,40, by = 5))
  
  
  everything_log <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("T_cells" = "red", "At" = "blue", "Lt5" = "purple", "B_cells" = "black", 
                                  "Z" = "lightblue", "Br" = "orange","Lt2" = "chartreuse", 
                                  "Lt3" = "chartreuse", "Lt4" = "chartreuse", "Lt" = "chartreuse", "Lr" = "aquamarine", 
                                  "f" = "pink","If" = "magenta", "Ct" = "brown", "Cb" = "darkorange2", "Tr" = "goldenrod3"), labels = c("T_cells" = "T cells")) +
    labs(title = "Within-Host Model", color = "Cell Type") + theme(legend.position = "right") + xlab(label = "Time (days post infection)") + 
    ylab(label = "Cell Number") + theme_classic() + scale_y_log10() + scale_x_continuous(limits = c(0,40), breaks = seq(0,40, by = 5))
  
  
  log_cyto_T <- ggplot(data = df_cyto_T, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("T_cells" = "red", "At" = "blue", "Lt5" = "purple","Ct" = "brown", "Lr" = "aquamarine", "Tr" = "goldenrod3"), labels = c("T_cells" = "T cells")) +
    labs(title = "Within-Host Model - log Cyto T", color = "Cell Type") + theme(legend.position = "right") + xlab(label = "Time (days post infection)") + 
    ylab(label = "Cell Number") + theme_classic() + scale_y_log10() + scale_x_continuous(limits = c(0,40), breaks = seq(0,40, by = 5))
  
  
  log_cyto_B <- ggplot(data = df_cyto_B, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
    scale_color_manual(values = c("B_cells" = "black", "Cb" = "darkorange2",  "At" = "blue", "Lr" = "aquamarine", "Br" = "orange"), labels = c("T_cells" = "T cells")) +
    labs(title = "Within-Host Model - log Cyto B", color = "Cell Type") + theme(legend.position = "right") + xlab(label = "Time (days post infection)") + 
    ylab(label = "Cell Number") + theme_classic() + scale_y_log10() + scale_x_continuous(limits = c(0,40), breaks = seq(0,40, by = 5))
  
  #plotting cytolytic data 
  cytolytic_plot <- ggplot(df_for_40000, aes(x = time/24)) +
    geom_line(aes(y = cytolytic_scale_40000_cb, color = "Cb")) + 
    geom_line(aes(y = cytolytic_scale_40000_ct, color = "Ct")) + 
    scale_color_manual(name = "cell type",values = c("Ct" = "red", "Cb" = "green")) +
    geom_point(data = baigent1998, aes(x = time/24, y = Bcell_no), color = "green",   size = 0.3, inherit.aes = FALSE) + 
    geom_point(data = baigent1998, aes(x = time/24, y = Tcell_no), color = "red",   size = 0.3, inherit.aes = FALSE) +
    labs(x = "Time (days post infection)", y = "pp38 Cell Counts per 40,000 cells", title = "Cytolytic Cells") + theme_classic() + 
    scale_x_continuous(limits = c(0,40), breaks = seq(0,40, by = 5))
  
  
  
  #plot for infected feather follicles
  FFE_plot <- ggplot(data = df_for_ffe, aes(x = time/24, y = IfPer10k)) + geom_line( color = "pink")  +
    labs(title = "WithinHost Delay")  + xlab(label = "Time (Days)") +
    ylab(label = "MDV genomes in FFE per 10,000 cells")  +  scale_y_log10(limits = c(0.01, 1000000000000)) + 
    geom_point(data = baigent2016, aes(x = time/24, y = mean_genomes), inherit.aes = FALSE) + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),   legend.text = element_text(size = 12),
          legend.title = element_text(size = 12), axis.line = element_line(color = "black")) +
    labs( y = "Cell Number", x = "Time (days post infection)", title = "Infected Feather Follicle Epithelium") + 
    theme_classic() + scale_x_continuous(limits = c(0,40), breaks = seq(0,40, by = 5))
  
  
  PBL_plot <- ggplot(data = df_for_PBL, aes(x = time/24, y = PBL_scale_10000 )) + geom_line( color = "red")  +
    scale_y_log10(limits = c(0.01, 1000000000000)) + 
    geom_point(data = baigent2016_PBL, aes(x = time/24, y = mean), inherit.aes = FALSE) + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),   legend.text = element_text(size = 12),
          legend.title = element_text(size = 12), axis.line = element_line(color = "black")) +
    labs( y = "MDV genomes in PBL per 10,000 cells", x = "Time (days post infection)", title = "Infected PBL") + 
    theme_classic() + scale_x_continuous(limits = c(0,40), breaks = seq(0,40, by = 5))
  
  Plaque_plot <- ggplot(data = df_Pfu, aes(x = time/24, y = Pfu_scale_1000000)) + geom_line(color = "purple") + 
    scale_y_log10(limits = c(0.01, 10000000)) + geom_point(data = cortes2004_plaque, aes(x = time/24, y = genome_mean), inherit.aes = FALSE) + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),   legend.text = element_text(size = 12),
          legend.title = element_text(size = 12), axis.line = element_line(color = "black")) +
    labs( y = "PFU Splenocytes per 1000000 cells", x = "Time (days post infection)", title = "Infected Splenocyte PFU") + 
    theme_classic() + scale_x_continuous(limits = c(0,40), breaks = seq(0,40, by = 5))
  
  
  
  # diagnostic <- df %>%
  #   tidyr::pivot_longer(
  #     cols = -time,
  #     names_to = "variable",
  #     values_to = "value"
  #   ) %>%
  #   dplyr::group_by(variable) %>%
  #   dplyr::summarise(
  #     minimum       = min(value, na.rm = TRUE),
  #     maximum       = max(value, na.rm = TRUE),
  #     n_zero        = sum(value == 0, na.rm = TRUE),
  #     n_negative    = sum(value < 0, na.rm = TRUE),
  #     n_tiny        = sum(value > 0 & value < 1e-12, na.rm = TRUE),
  #     n_nonfinite   = sum(!is.finite(value)),
  #     .groups = "drop"
  #   )
  # 
  # diagnostic
  
  print(everything)
  print(everything_log) 
  print(log_cyto_T) 
  print(log_cyto_B)
  print(cytolytic_plot)  
  print(PBL_plot)
  print(FFE_plot) 
  print(Plaque_plot)
  
}



# take out a named vector that is in dataframe iteratively 

## PARAMETERS AND INITIAL VALUES ## 

## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  beta =  6.951463e-07                  #contact rate with B cells   
  # ,beta_2 =  6.951463e-07                  #contact rate with B cells
  , Pb = 0.005 
  , Pt = 0.005
  , nu_a = 0.1666667            #Activation rate of T cells by cytolytic B cells (hours) 
  , nu_b = 0.4467098
  , nu_f = 0.008                    #Infection rate of follicular cells (hours) 
  , kappa = 1e-6 
  , mu = 1/8                       #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.0104          #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.0104          #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 0                    #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 = 0                      #incoming T cells / determined no incoming T cells 
  , h2 = 100
  , lambda = 0.02380952             # fixing delay rate to (1/(7*24))*4   
  , size_pp38 = 10  # new parameter 
  , size_pp38_Ct = 10
  , q_FFE = 250
  # , q_PBL = log(115)
)


initial_values <- c( 
  B_cells = 2.4e9/3    # from three organs 2.4e6/3 
  , Cb = 1  
  , Br = 2.4e9/3       #same initial proportion here infected as in model (well mixed marbles ex) 
  , T_cells = 7.4e8/3  # from three organs 1.5e6/3  
  , Tr = 7.4e8/3
  , At = 0 
  , Lt = 0 
  , Lt2 = 0 
  , Lt3 = 0  
  , Lt4 = 0  
  , Lt5 = 0  
  , Lr = 0 
  , Ct = 0
  , Z = 0
  , f = 4e5 # (Calnek et al., 1970) range of 1.5e5 - 7e5  
  , If =0 
  
) 


time_values <- seq(0, 1080) # hours  


### DATA FOR PLOT ### 

### DATA FOR PLOT ### 
baigent2016 <- read_xlsx("Baigent2016/Unvax/baigent2016.xlsx", sheet = 2) %>% mutate(time = time*24)

baigent1998 <- read_xlsx("Baigent1998/baigent1998.xlsx", sheet = 3, na = "NA") %>% filter(!is.na(mean.pp38)) %>% dplyr::select(time, Bcell_no,Tcell_no) 

baigent2016_PBL <- read_xlsx("Baigent2016/Unvax/PBL_No_Vax_fin.xlsx") %>% filter(time > -1) %>% mutate(time = round(time*24, 0))

cortes2004_plaque <- read_xlsx("Cortes2004/RB1B_plaqueAssayPchicks.xlsx")  

optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/SimplifyingModel_v03_20260722_021844.csv") %>% 
  filter(Converged == 0) %>% slice_min(Likely, n = 10) %>% arrange(Likely) %>% 
  dplyr::select(c(beta, alpha, alpha_2,nu_a,nu_f,mu,Pb,Pt, q_FFE,size_pp38, size_pp38_Ct,nu_b, kappa)) 


#creating a list of dataframes set in motion by this function made on the fly that takes optim data and i takes in seq_len(nrow(()))
list_of_df <-  purrr::map(seq_len(nrow(optim_data)), function(i) { 
  parameter_vector(dat = optim_data, i = i)
}) 

pdf("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/.pdf", width = 7, height = 5)

generating_plots <-  purrr::map(seq_len(nrow(optim_data)), function(i) {
  creating_plots(listofdf = list_of_df, i = i)
})

dev.off()

