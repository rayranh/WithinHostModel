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
  
  infprob <- results %>% 
    dplyr::select(B_cells, Cb,Br, Ct, T_cells, At, Lt, Lt2, Lt3, Lt4, Lt5, time) %>% 
    mutate(
      prob_Cyto = (Cb + Ct) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5),
      inf_cells = ((Cb + Ct) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5))*40000 # getting rid of Lt because Lt is not even present at this time point
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
  
  
  pp38_dbinom <- pp38_dat %>% filter(time %in% obs_hourspp38) %>% 
    left_join(infprob %>% select(time, inf_cells), by = "time") %>% group_by(time) %>%
    mutate(variance = sd(mean.pp38)^2, r= (inf_cells^2)/(variance-inf_cells), p = (inf_cells/variance)) 
  
  loglike_pp38 <- dnbinom(pp38_dbinom$mean.pp38,
                          size = pp38_dbinom$r, prob = pp38_dbinom$p,  log = TRUE)
  #loglike_pp38 <- dbinom(pp38_dbinom$mean.pp38, prob = pp38_dbinom$prob_Cyto, size = 40000, log = TRUE) 
  
  sum_loglike_pp38 <- sum(loglike_pp38)
  
  
  loglike_ffe <- dnorm(
    log10(ffe_join$mean_genomes),
    mean = log10(ffe_join$If),
    sd = 0.642, # calculated from 17dpi onwards 
    log = TRUE
  )  
  
  sum_loglike_ffe <- sum(loglike_ffe) 
   
  LogLHoodStor <- sum_loglike_ffe + sum_loglike_pp38 
  
  list(neg_LogLHoodStor = -LogLHoodStor,
       neg_sum_loglike_ffe = -sum_loglike_ffe, 
       neg_sum_loglike_pp38 =  -sum_loglike_pp38)
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

### TIME ### 

#matching to the wack hours of baigent2016
matched_time <- c(3,6,10,17,20,26,33,38)*24 

#PP38 Time, Right now the time is only including some small amount 
obs_hourspp38 <- c(72,96,120,144) # make sure pp38 times is matching the data loaded in 

time_values <- seq(0, 1080) # hours  


## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 
pp38_dat <- read_xlsx("baigent1998.xlsx", sheet = 3, na = "NA") %>% filter(!is.na(mean.pp38)) %>% select(time, mean.pp38)
baigent2016 <- read_xlsx("baigent2016.xlsx", 2 ) %>% arrange(time)


optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Feb.17.26.FittingDnbinom_Run2_ALL.csv") %>% 
  filter(Converged == 0) %>% slice_min(Likely, n = 10) %>% select(c(beta, beta_2, alpha, alpha_2,nu_a,nu_b,nu_f,mu,g1,g2,h1,h2,Pb)) 


list_of_df <-  purrr::map(seq_len(nrow(optim_data)), function(i) { 
  parameter_vector(dat = optim_data, i = 1)
}) 

Likelihood_10 <- list_of_df[[1]]

print(Likelihood_10)
