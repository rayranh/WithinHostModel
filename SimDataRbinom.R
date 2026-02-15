###plotting the top 10 best likelihood scores in ggplot from csv file ###  
### Make sure same starting initial conditions ### 

#rm(list = ls())
library(deSolve)  
library(tidyr) 
library(ggplot2)
library(dplyr) 
library(readxl)  
library(purrr) 
library(tictoc)

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
  
  df <- results %>% mutate(prob_Cyto_B = (Cb)/(Cb+Ct+At+Lt+Lt2+Lt3+Lt4+Lt5+B_cells+T_cells+Br), 
                             prob_Cyto_T = (Ct)/(Cb+Ct+At+Lt+Lt2+Lt3+Lt4+Lt5+B_cells+T_cells+Br), 
                             prob_Cyto = (Ct+Cb)/(Cb+Ct+At+Lt+Lt2+Lt3+Lt4+Lt5+B_cells+T_cells+Br), 
                           numInfCells = (Ct+Cb)/(Cb+Ct+At+Lt+Lt2+Lt3+Lt4+Lt5+B_cells+T_cells+Br)*40000) %>% 
    filter(time %in% obs_hourspp38) %>% select(prob_Cyto_B, prob_Cyto_T, prob_Cyto, numInfCells, time) 
  
  
  
  #counting number of birds in each timepoint 
  n_by_time <- baigent1998 %>%
    filter(time %in% obs_hourspp38) %>%
    count(time, name = "n_birds") 
  
  #attaching number of birds to each timepoint 
  cyto_df2 <- df %>% left_join(y = n_by_time, by = "time")
  
  
  return(cyto_df2)
} 



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
obs_hourspp38 <- c(72,96,120,144) 



### DATA FOR PLOT ### 
baigent2016 <- read_xlsx("baigent2016.xlsx", 2 ) 
baigent1998 <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", 3 ) %>% 
  mutate(mean.pp38 = as.numeric(mean.pp38),
         Bcell_no = as.numeric(Bcell_no), 
         Tcell_no = as.numeric(Tcell_no)) %>% filter(!is.na(mean.pp38))   
optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Feb.11.26.FittingDnbinom_mu_theta.csv") %>% 
  filter(Converged == 0) %>% slice_min(Likely, n = 10) %>% select(c(beta, beta_2, alpha, alpha_2,nu_a,nu_b,nu_f,mu,g1,g2,h1,h2,Pb)) 


 list_of_df <-  purrr::map(seq_len(nrow(optim_data)), function(i) { 
  parameter_vector(dat = optim_data, i = i)
}) 


# selecting only 1 data frame at a time 
 OneDf <- list_of_df[[1]]
 
 sim_df <- purrr::map_df(1:5000,
   function(rep_id) {
       purrr::map_df(seq_len(nrow(OneDf)), function(i) { 
         n_birds <- OneDf$n_birds[i]
         tibble(
           replicate = rep_id,
           time = OneDf$time[i],
           bird_no = seq_len(n_birds),
           pp38_sim = rnbinom(
             n = n_birds,
             size = 0.5,
             mu = OneDf$numInfCells[i] 
           )
         )
       })
     })

 
 sim_df_median <- sim_df %>%
   group_by(time) %>%
   summarise(med = median(pp38_sim),
             low = quantile(pp38_sim, 0.025), 
             high = quantile(pp38_sim, 0.975), 
             .groups = "drop") 
 
 
 
 
 obs_with_band <- baigent1998 %>%
   filter(time %in% obs_hourspp38) %>%
   left_join(sim_df_median, by = "time") %>%
   mutate(outside = mean.pp38 < low | mean.pp38 > high)
 
 n_total   <- nrow(obs_with_band)
 n_outside <- sum(obs_with_band$outside, na.rm = TRUE)
 n_inside  <- n_total - n_outside
 
 prop_outside <- n_outside / n_total
 prop_inside  <- n_inside  / n_total
 
 prop_inside
 prop_outside
 
 
 #5000 takes 60 sec to run 

 
p <- ggplot(sim_df_median, aes(x = time/24)) + geom_ribbon(aes(ymin = low, ymax = high), fill = "grey80", alpha = 0.5) + geom_line(aes(y = med, colour = "model")) +
  xlim(0,10) + geom_point(data = baigent1998, aes(x = time/24, y = mean.pp38, color = "data"), inherit.aes = FALSE) + 
  scale_color_manual(values = c("model" = "black", "data" = "red")) +  ylim(0,500) + 
  labs(title = "Simulated pp38 data", x = "time(days)", y = "No. of pp38 Infected Cells") + 
  theme(panel.grid = element_blank()) + 
  theme_classic()

print(p)

