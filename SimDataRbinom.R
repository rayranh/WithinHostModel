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
  
  df <- results %>% mutate(numInfCells = (Cb)/(Cb+Ct+At+Lt+Lt2+Lt3+Lt4+Lt5+B_cells+T_cells+Br)*40000, 
                           numInfCells_Ct = (Ct)/(Cb+Ct+At+Lt+Lt2+Lt3+Lt4+Lt5+B_cells+T_cells+Br)*40000) %>% 
    filter(time %in% obs_hourspp38) %>% select(numInfCells,numInfCells_Ct, time) 
  
  
  #counting number of birds in each timepoint 
  n_by_time <- baigent1998 %>%
    filter(time %in% obs_hourspp38) %>%
    count(time, name = "n_birds") 
  
  #attaching number of birds to each timepoint 
  cyto_df2 <- df %>% left_join(y = n_by_time, by = "time") %>% left_join(y = baigent1998, by = "time") %>% 
    select(numInfCells,numInfCells_Ct, time,n_birds) %>% group_by(time) %>% 
    mutate(size_Cb = pars["size_pp38"], size_Ct = pars["size_pp38_Ct"]) %>% 
    summarise(
      numInfCells = first(numInfCells),
      numInfCells_Ct = first(numInfCells_Ct),
      size_Cb = first(size_Cb),
      size_Ct = first(size_Ct),
      .groups = "drop"
    ) %>% 
    ungroup() 


  return(cyto_df2)
}  


parameter_vector_FFE <- function(dat,i) { 
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
  
  df_FFE <- results %>% filter(time %in% matched_time) %>% select(time, If) %>% left_join(baigent2016, by = "time")
 
  

  return(df_FFE)
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
  , size_pp38 = log(10)   # new parameter 
  , size_pp38_Ct = log(10)
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
#matching to the wack hours of baigent2016
matched_time <- c(3,6,10,17,20,26,33,38)*24 



### DATA FOR PLOT ### 
baigent2016 <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent2016.xlsx", 2) %>% mutate(time = time*24)
baigent1998 <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", 3 ) %>% 
  mutate(mean.pp38 = as.numeric(mean.pp38),
         Bcell_no = as.numeric(Bcell_no), 
         Tcell_no = as.numeric(Tcell_no)) %>% filter(!is.na(mean.pp38))   
optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/FittingForVariance_ALL.csv") %>% 
  filter(Converged == 0) %>% slice_min(Likely,n = 1) %>%  
  select(c(beta, beta_2, alpha, alpha_2,nu_a,nu_b,nu_f,mu,g1,g2,h1,h2,Pb, size_pp38, size_pp38_Ct)) 

# OneDf <- parameter_vector(dat = optim_data, 1) 

OneDf_all <- map_df(1:nrow(optim_data), ~ parameter_vector(optim_data, .x))

# OneDfFFe <- parameter_vector_FFE(dat = optim_data, 1)
OneDf_FFE_all <- map_df(1:nrow(optim_data), ~ parameter_vector_FFE(optim_data, .x)) 
 
 # ------------ Simulating Data for pp38 ------------ # 
 
 sim_df <- purrr::map_df(1:10000,
   function(rep_id) { 
     OneDf_all %>% 
       transmute(time = time, 
                 pp38_sim_Cb = rnbinom(n(), size = size_Cb , mu = numInfCells),
                 pp38_sim_Ct = rnbinom(n(), size = size_Ct , mu = numInfCells_Ct))
     })

 
 # ------------ Simulating Data for FFE ------------ # 
 
 sim_df_FFE <- purrr::map_df(1:5000,
                         function(rep_id) {
                           purrr::map_df(seq_len(nrow(OneDf_FFE_all)), function(i) { 
                             tibble(
                               replicate = rep_id,
                               time = OneDf_FFE_all$time[i],
                               FFE_sim = 10^rnorm(  
                                 n = 1, 
                                 mean = log10(OneDf_FFE_all$If[i]), 
                                 sd = 0.642 
                               )
                             )
                           })
                         }) 

 
 sim_df_median <- sim_df %>%
   group_by(time) %>%
   summarise(med_Cb = median(pp38_sim_Cb), 
             mean_Cb = mean(pp38_sim_Cb),
             low_Cb = quantile(pp38_sim_Cb, 0.025), 
             high_Cb = quantile(pp38_sim_Cb, 0.975), 
             med_Ct = median(pp38_sim_Ct), 
             mean_Ct = mean(pp38_sim_Ct),
             low_Ct = quantile(pp38_sim_Ct, 0.025), 
             high_Ct = quantile(pp38_sim_Ct, 0.975), 
             .groups = "drop") 
 
 sim_FFE_median <- sim_df_FFE %>% 
   group_by(time) %>% summarise(med = median(FFE_sim), 
                                low = quantile(FFE_sim, 0.025), 
                                high = quantile(FFE_sim, 0.975), 
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
 
 
 obs_with_band_FFE <- baigent2016 %>%
   left_join(sim_FFE_median, by = "time") %>%
   mutate(outside = mean_genomes < low | mean_genomes > high)
 
 n_total_FFE   <- nrow(obs_with_band_FFE)
 n_outside_FFE <- sum(obs_with_band_FFE$outside, na.rm = TRUE)
 n_inside_FFE  <- n_total_FFE - n_outside_FFE
 
 prop_inside_FFE <- n_inside_FFE / n_total_FFE
 prop_outside_FFE <- n_outside_FFE / n_total_FFE
 
 prop_inside_FFE
 prop_outside_FFE
 
 
 
 #5000 takes 60 sec to run 

 
p <- ggplot(sim_df_median, aes(x = time/24)) + geom_line(aes(y = mean_Cb, colour = "model B cells")) +  
  geom_line(aes(y = mean_Ct, colour = "model T cells")) + 
  geom_ribbon(aes(ymin = low_Cb, ymax = high_Cb), fill = "grey80", alpha = 0.5) + 
  geom_ribbon(aes(ymin = low_Ct, ymax = high_Ct), fill = "green", alpha = 0.5) + 
  xlim(0,10) + geom_point(data = baigent1998, aes(x = time/24, y = Bcell_no, color = "data B cells"), inherit.aes = FALSE) + 
  geom_point(data = baigent1998, aes(x = time/24, y = Tcell_no, color = "data T cells"), inherit.aes = FALSE) + 
  scale_color_manual(values = c("model B cells" = "black", "data B cells" = "red" , "data T cells" = "green", "model T cells" = "darkgreen")) + 
  labs(title = "Data Generated from Model", x = "time(days post infection)", y = "# of Infected Cells") + 
  theme(panel.grid = element_blank()) + 
  theme_classic()

p_data <- ggplot(baigent1998, aes(x = time/24, y = mean.pp38, color = "data")) + 
  geom_point() + xlim(0,10) + scale_color_manual(values = c("model" = "black", "data" = "red")) +  ylim(0,500) + 
  labs(title = "Data", x = "time(days post infection)", y = "# of Infected Cells") + 
  theme(panel.grid = element_blank()) + 
  theme_classic()

p2 <- ggplot(sim_FFE_median, aes(x = time/24))  + geom_line(aes(y = med, colour = "model")) + 
  geom_ribbon(aes(ymin = low, ymax = high), fill = "grey80", alpha = 0.5) + 
  scale_y_log10(limits = c(0.0000001, 1000000000)) +
  geom_point(data = baigent2016, aes(x = time/24, y = mean_genomes, color = "data"), inherit.aes = FALSE)+  
  # geom_errorbar(data = baigent2016, aes(ymin = lower.conf, ymax = upper.conf, color = "data")) + 
  xlim(0,40)+scale_color_manual(values = c("model" = "black", "data" = "red"))  + 
  labs(title = "Data Generated from Model ", x = "time(days post infection)", y = "Mean Number of MDV Genomes") + 
  theme(panel.grid = element_blank()) + 
  theme_classic() 

p2_data <- ggplot(baigent2016, aes(x = time/24, y = mean_genomes, color = "data")) + geom_point() + 
  geom_errorbar(aes(ymin = lower.conf, ymax = upper.conf )) + 
  scale_y_log10(limits = c(0.0000001, 1000000000)) +
  xlim(0,40) + scale_color_manual(values = c("model" = "black", "data" = "red"))  + 
  labs(title = "Data", x = "time(days post infection)", y = "Mean Number of MDV Genomes") + 
  theme(panel.grid = element_blank()) + 
  theme_classic()


print(p) 
print(p2)

