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
setwd("~/Desktop/WithinHostModel/DataForProject/") 



sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
    dB <-  -beta*Cb*B_cells - beta*Ct*B_cells + Pb*(g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct  
    dBr <- (1-Pb)*(g1*(Cb+Ct)/(g2+(Cb+Ct))) #probably should have alpha for cytolytic infection 
    dCb <- beta*Cb*B_cells + beta*Ct*B_cells - alpha*Cb 
    dT <- -nu_a*Cb*T_cells - nu_a*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct)))
    dAt <- nu_a*Cb*T_cells + nu_a*Ct*T_cells - beta*Ct*At - beta*Cb*At  # beta = rate of activated T cells by Ct and Cb cells   
    dLt <- theta *(beta*Ct*At + beta*Cb*At)  - lambda*Lt 
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - mu*Lt5
    dCt <-  (1-theta)*(beta*Ct*At + beta*Cb*At) - alpha_2*Ct
    dZ <- mu*Lt5
    df <- -nu_f*Lt5*f
    dIf <- nu_f*Lt5*f
    return(list(c(dB, dCb,dBr, dT, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dCt, dZ, df,dIf)))
  }) 
} 

#creating model predicted means 
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
  
  cyto_df2 <- df  %>%  mutate(size_Cb = pars["size_pp38"], size_Ct = pars["size_pp38_Ct"]) 


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
  

  df_FFE <- results %>% filter(time
                               %in% matched_time) %>% mutate(If_per10k = pars["q_FFE"] * (If/(If+f)) *10000) %>% 
  select(time, If_per10k) %>% left_join(baigent2016, by = "time")
 
  

  return(df_FFE)
}   


parameter_vector_PBL <- function(dat,i) {
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

 df_PBL <- results %>% filter(time %in% matched_time) %>% mutate(PBL_pred = pars["q_PBL"] * ((Cb + Ct + Lt + Lt2 + Lt3 + Lt4 + Lt5) /
                                                                               (Cb + Ct + At + Lt + Lt2 + Lt3 + Lt4 + Lt5 +
                                                                                  B_cells + T_cells + Br)) * 10000) %>% 
   select(time, PBL_pred) %>% left_join(baigent2016PBL, by = "time")

 return(df_PBL)
 }



parameters_values <- c( 
  beta =  6.951463e-07                  #contact rate with B cells  
  , Pb = 0.005
  , nu_a = 4.668718e-01                     #Activation rate of T cells by cytolytic B cells (hours)
  , nu_f = 0.008                     #Infection rate of follicular cells (hours)
  , mu = 1/8                        #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.0104                   #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.0104                 #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 0                    #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 = 0                      #incoming T cells / determined no incoming T cells 
  , h2 = 100
  , lambda =0.02380952                  #adding delay, how long latent cell 'exposed'  
  , size_pp38 = log(10)   # new parameter 
  , size_pp38_Ct = log(10) 
  , q_FFE = 115 
  , q_PBL = 115

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
  , f = 4e5  
  , If =0 
  
) 


time_values <- seq(0, 1080) # hours   
obs_hourspp38 <- c(72,96,120,144)  
#matching to the wack hours of baigent2016
matched_time <- c(3,6,10,17,20,26,33,38)*24 



### DATA FOR PLOT ### 
baigent2016 <- read_xlsx("Baigent2016/Unvax/feathers_noVax_SE.xlsx") %>% 
  mutate(mean_genomes= 10^(logged10Mean))   #switch to 3 for plot and 2 for simulation  because you fit to 2 

baigent1998 <- read_xlsx("Baigent1998/baigent1998Totalpp38.xlsx" )
baigent2016PBL <- read_xlsx("Baigent2016/Unvax/PBL_noVax_SE.xlsx")
optim_data <- read.csv("ScalingParamTestFeathers") %>% 
  filter(Converged == 0) %>% slice_min(Likely,n = 1) %>%  
  select(c(beta,  alpha, alpha_2,nu_a,nu_f,mu,Pb, size_pp38, size_pp38_Ct, q_FFE, q_PBL))  



# OneDf <- parameter_vector(dat = optim_data, 1) 

OneDf_all <- map_df(1:nrow(optim_data), ~ parameter_vector(optim_data, .x))

# OneDfFFe <- parameter_vector_FFE(dat = optim_data, 1)
OneDf_FFE_all <- map_df(1:nrow(optim_data), ~ parameter_vector_FFE(optim_data, .x))  

OneDf_PBL_all <- map_df(1:nrow(optim_data), ~ parameter_vector_PBL(optim_data, .x)) 

 
 # ------------ Simulating Data for pp38 ------------ # 
 
 sim_df <- purrr::map_df(1:5000,
   function(rep_id) { 
     OneDf_all %>% 
       transmute(time = time,  
                 replicate = rep_id,
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
                                 mean = log10(OneDf_FFE_all$If_per10k[i]), 
                                 sd = 0.352 
                               )
                             )
                           })
                         }) 

 # ------------ Simulating Data for PBL ------------ # 
 
 sim_df_PBL <- purrr::map_df(1:5000,
                             function(rep_id) {
                               purrr::map_df(seq_len(nrow(OneDf_PBL_all)), function(i) { 
                                 tibble(
                                   replicate = rep_id,
                                   time = OneDf_PBL_all$time[i],
                                   PBL_sim = 10^rnorm(  
                                     n = 1, 
                                     mean = log10(OneDf_PBL_all$PBL_pred[i]), 
                                     sd = 0.365 
                                     
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
 
 sim_PBL_median <- sim_df_PBL %>%
   group_by(time) %>% summarise(med = median(PBL_sim),
                                low = quantile(PBL_sim, 0.025),
                                high = quantile(PBL_sim, 0.975),
                                .groups = "drop") 
 

 # ------------ Simulating Data for PBL ------------ # 
 

 # obs_with_band <- baigent1998 %>%
 #   filter(time %in% obs_hourspp38) %>%
 #   left_join(sim_df_median, by = "time") %>%
 #   mutate(outside = mean.pp38 < low | mean.pp38 > high)
 # 
 # n_total   <- nrow(obs_with_band)
 # n_outside <- sum(obs_with_band$outside, na.rm = TRUE)
 # n_inside  <- n_total - n_outside
 # # 
 # prop_outside <- n_outside / n_total
 # prop_inside  <- n_inside  / n_total
 # 
 # prop_inside
 # prop_outside
 # # 
 # # 
 # obs_with_band_FFE <- baigent2016 %>%
 #   left_join(sim_FFE_median, by = "time") %>%
 #   mutate(outside = mean_genomes < low | mean_genomes > high)
 # 
 # n_total_FFE   <- nrow(obs_with_band_FFE)
 # n_outside_FFE <- sum(obs_with_band_FFE$outside, na.rm = TRUE)
 # n_inside_FFE  <- n_total_FFE - n_outside_FFE
 # 
 # prop_inside_FFE <- n_inside_FFE / n_total_FFE
 # prop_outside_FFE <- n_outside_FFE / n_total_FFE
 # 
 # prop_inside_FFE
 # prop_outside_FFE
 # 
 # 
 
 #5000 takes 60 sec to run 

 
p_cb <- ggplot(sim_df_median, aes(x = time/24)) + geom_line(aes(y = med_Cb, colour = "model B cells")) +  
  geom_ribbon(aes(ymin = low_Cb, ymax = high_Cb), fill = "grey80", alpha = 0.5) + 
  coord_cartesian(xlim = c(0,10)) + geom_point(data = baigent1998, aes(x = time/24, y = pp38BcellTotal, color = "data B cells"), inherit.aes = FALSE) + 
  scale_color_manual(values = c("model B cells" = "black", "data B cells" = "red")) + 
  labs(title = "Data Generated from Model - B cell Infection", x = "time(days post infection)", y = "# of Infected Cells") + 
  theme(panel.grid = element_blank()) + 
  theme_classic()  
p_Ct <- ggplot(sim_df_median, aes(x = time/24)) + geom_line(aes(y = med_Ct, colour = "model T cells")) +  
  geom_ribbon(aes(ymin = low_Ct, ymax = high_Ct), fill = "grey80", alpha = 0.5) + 
  xlim(0,10) +
  geom_point(data = baigent1998, aes(x = time/24, y = pp38TcellTotal, color = "data T cells"), inherit.aes = FALSE) + 
  scale_color_manual(values = c("data T cells" = "red", "model T cells" = "black")) + 
  labs(title = "Data Generated from Model - T cell Infection", x = "time(days post infection)", y = "# of Infected Cells") + 
  theme(panel.grid = element_blank()) + 
  theme_classic() 

p_Tdata <- ggplot(baigent1998, aes(x = time/24, y = pp38TcellTotal, color = "data")) + 
  geom_point() + xlim(0,10) + scale_color_manual(values = c("model" = "black", "data" = "red")) +  
  labs(title = "T cell infection Data", x = "time(days post infection)", y = "# of Infected Cells per 40,000 lymphocytes") + ylim(0,200) +
  theme(panel.grid = element_blank()) + 
  theme_classic() 

p_Bdata <- ggplot(baigent1998, aes(x = time/24, y = pp38BcellTotal, color = "data")) + 
  geom_point() + xlim(0,10) + scale_color_manual(values = c("model" = "black", "data" = "red")) +  
  labs(title = "B cell infection Data", x = "time(days post infection)", y = "# of Infected Cells per 40,000 lymphocytes") + ylim(0,600) +
  theme(panel.grid = element_blank()) + 
  theme_classic()



p_Feathers <- ggplot(sim_FFE_median, aes(x = time/24))  + geom_line(aes(y = med, colour = "model")) + 
  geom_ribbon(aes(ymin = low, ymax = high), fill = "grey80", alpha = 0.5) + 
  scale_y_log10(limits = c(0.0000001, 1000000000)) +
  geom_point(data = baigent2016, aes(x = time/24, y = mean_genomes, color = "data"), inherit.aes = FALSE)+  
  geom_errorbar(data = baigent2016, aes(ymin = mean_genomes, ymax = ConfInt, color = "data", width = 0.7)) + 
  xlim(0,40)+scale_color_manual(values = c("model" = "black", "data" = "red"))  + 
  labs(title = "Data Generated from Model - Feather Infection ", x = "time(days post infection)", y = "Mean Number of MDV Genomes") + 
  theme(panel.grid = element_blank()) + 
  theme_classic() 

p2_data <- ggplot(baigent2016, aes(x = time/24, y = mean_genomes, color = "data")) + geom_point() + 
  geom_errorbar(aes(ymin = mean_genomes, ymax = ConfInt )) + 
  scale_y_log10(limits = c(0.0000001, 1000000000)) +
  xlim(0,40) + scale_color_manual(values = c("model" = "black", "data" = "red"))  + 
  labs(title = "Data", x = "time(days post infection)", y = "Mean Number of MDV Genomes") + 
  theme(panel.grid = element_blank()) + 
  theme_classic()

p_PBL <- ggplot(sim_PBL_median, aes(x = time/24))  + geom_line(aes(y = med, colour = "model")) + 
  geom_ribbon(aes(ymin = low, ymax = high), fill = "grey80", alpha = 0.5) + 
  scale_y_log10(limits = c(0.0000001, 1000000000)) +
  geom_point(data = baigent2016PBL, aes(x = time/24, y = mean, color = "data"), inherit.aes = FALSE)+  
  geom_errorbar(data = baigent2016PBL, aes(ymin = mean, ymax = ConfInt, color = "data", width = 0.7)) + 
  xlim(0,40)+scale_color_manual(values = c("model" = "black", "data" = "red"))  + 
  labs(title = "Data Generated from Model - PBL Infection", x = "time(days post infection)", y = "Mean Number of MDV Genomes") + 
  theme(panel.grid = element_blank()) + 
  theme_classic() 

