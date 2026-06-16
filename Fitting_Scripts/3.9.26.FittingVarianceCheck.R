# This script is confirmation that 3.9.26.FittingVariance.R code in the cluster is the same model used to create the plots I used for my Comp Meeting 
rm(list = ls())
library(dplyr) 
library(tibble)
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(writexl) 
library(readxl)  
library(purrr)  
library(future) 
library(future.apply) 


## SIR MODEL ## 
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

#### LIKELIHOOD FUNCTION #### 

Likelihood_parts <- function(params){   
  pars <- parameters_values 
  init <- initial_values 
  B0 <-2.4e9/3 
  
  pars[names(params)] <- params  
  pars["beta"]   <- exp(pars["beta"])
  pars["beta_2"] <- exp(pars["beta_2"])  
  pars["alpha"]   <- exp(pars["alpha"])
  pars["alpha_2"] <- exp(pars["alpha_2"])  
  pars["nu_a"] <- exp(pars["nu_a"])
  pars["nu_b"] <- exp(pars["nu_b"]) 
  pars["nu_f"] <- exp(pars["nu_f"])
  pars["Pb"] <- plogis(pars["Pb"])  
  pars["size_pp38"] <- exp(pars["size_pp38"]) 
  pars["size_pp38_Ct"] <- exp(pars["size_pp38_Ct"])  
  
  
  
  init["B_cells"] <- B0*pars["Pb"] 
  init["Br"] <- B0* (1 - pars["Pb"])
  
  # run model
  results <- as.data.frame(
    ode(y = init, 
        times = time_values, 
        func = sir_equations, 
        parms = pars)
  ) 
  
  ### MODEL ###  
  
  # percent cytolytically infected in model
  infprob <- results %>% 
    dplyr::select(B_cells, Cb,Br, Ct, T_cells, At, Lt, Lt2, Lt3, Lt4, Lt5, time) %>% 
    mutate(
      prob_Cyto = (Cb + Ct) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5),
      inf_cells = ((Cb + Ct) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5))*40000,
      inf_cells_Cb = ((Cb) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5))*40000,  
      inf_cells_Ct = ((Ct) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5))*40000# getting rid of Lt because Lt is not even present at this time point
    ) %>% filter(time %in% obs_hourspp38) 
  
  # baigent2016 times are in days; matched_time is hours 
  ffe_dat <- baigent2016 %>%
    transmute(time = time * 24, mean_genomes)
  
  #matching model time to feathers time 
  Infected_FFE<- results %>% dplyr::select("time", "If") %>% 
    filter(time %in% matched_time) %>% arrange(time)  
  
  # join model predictions (If) to observed data by time
  ffe_join <- ffe_dat %>%
    left_join(Infected_FFE, by = "time")  # brings in If
  
  
  
  pp38_dbinom_Cb <- pp38_dat %>% filter(time %in% obs_hourspp38) %>% 
    left_join(infprob %>% select(time, inf_cells_Cb), by = "time") 
  
  pp38_dbinom_Ct <- pp38_dat %>% filter(time %in% obs_hourspp38) %>% 
    left_join(infprob %>% select(time, inf_cells_Ct), by = "time") 
  
  loglike_pp38 <- dnbinom(pp38_dbinom_Cb$Bcell_no,# from data, observed 
                          size = pars["size_pp38"], prob = pars["size_pp38"]/(pars["size_pp38"]+pp38_dbinom_Cb$inf_cells_Cb),  log = TRUE) # cytolytic cells from model 
  
  loglike_pp38_Ct <- dnbinom(pp38_dbinom_Ct$Tcell_no,# from data, observed 
                             size = pars["size_pp38_Ct"], prob = pars["size_pp38_Ct"]/(pars["size_pp38_Ct"]+pp38_dbinom_Ct$inf_cells_Ct),  log = TRUE) # cytolytic cells from model 
  
  
  sum_loglike_pp38 <- sum(loglike_pp38) + sum(loglike_pp38_Ct)
  
  
  loglike_ffe <- dnorm(
    log10(ffe_join$mean_genomes), #log10() here to match the log10 of the data and keep scale true
    mean = log10(ffe_join$If),
    sd = 0.642, # calculated from 17dpi onwards 
    log = TRUE # outer log is TRUE so it is returning ln() to keep consistent
  )  
  
  sum_loglike_ffe <- sum(loglike_ffe) 
  
  LogLHoodStor <- sum_loglike_ffe + sum_loglike_pp38
  
  list(neg_LogLHoodStor = -LogLHoodStor,
       neg_sum_loglike_ffe = -sum_loglike_ffe, 
       neg_sum_loglike_pp38 =  -sum_loglike_pp38)
  
} 



Likelihood <- function(params){  
  
  Likelihood_parts(params)$neg_LogLHoodStor
}



## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  beta =  log(6.951463e-07)                  #contact rate with B cells  
  , Pb = qlogis(0.005) 
  , beta_2 = log(8.532282e-07)                 #contact rate with T cells 
  , nu_a = log(4.668718e-01)             #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = log(4.467098e-01)             #Activation rate of T cells by cytolytic T cells (hours)
  , nu_f = log(0.008)                    #Infection rate of follicular cells (hours)
  , mu = 1/8                       #Rate of Tumor Cells (every 72 hours)
  , alpha = log(0.0104)             #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = log(0.0104)           #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 10                      #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 = 10                      #incoming T cells / determined no incoming T cells 
  , h2 = 100
  , lambda = 0.02380952             # fixing delay rate to (1/(7*24))*4   
  , size_pp38 = log(10)   # new parameter 
  , size_pp38_Ct = log(10)
)


initial_values <- c( 
  B_cells = 2.4e9/3    # from three organs 2.4e6/3 
  , Cb = 1  
  , Br = 2.4e9/3       #same initial proportion here infected as in model (well mixed marbles ex) 
  , T_cells = 7.4e8/3  # from three organs 1.5e6/3 
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
### TIME ### 

time_values <- seq(0, 1080, by = 1)   
#matching to the wack hours of baigent2016
matched_time <- c(3,6,10,17,20,26,33,38)*24 

#PP38 Time, Right now the time is only including some small amount 
obs_hourspp38 <- c(72,96,120,144) # make sure pp38 times is matching the data loaded in 


## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 
pp38_dat <- read_xlsx("Baigent1998/baigent1998.xlsx", sheet = 3, na = "NA") %>% filter(!is.na(mean.pp38)) %>% select(time, Bcell_no,Tcell_no)

baigent2016 <- read_xlsx("Baigent2016/Unvax/baigent2016.xlsx", sheet = 2 ) %>% arrange(time)


data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/3.9.26.FittingVarianceReproducingResults3.csv")  %>% 
  filter(Converged == 0) %>% arrange(Likely) 

bestrow <- data[1,]

params_check <- c(
  beta = log(bestrow$beta),
  beta_2 = log(bestrow$beta_2),
  alpha = log(bestrow$alpha),
  alpha_2 = log(bestrow$alpha_2),
  nu_a = log(bestrow$nu_a),
  nu_b = log(bestrow$nu_b),
  nu_f = log(bestrow$nu_f),
  Pb = qlogis(bestrow$Pb),
  g1 = bestrow$g1,
  g2 = bestrow$g2,
  h1 = bestrow$h1,
  h2 = bestrow$h2,
  size_pp38 = log(bestrow$size_pp38),
  size_pp38_Ct = log(bestrow$size_pp38_Ct)
)

check <- Likelihood_parts(params_check)

check


