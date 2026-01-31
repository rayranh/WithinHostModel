#Using this code to generate random numbers for ALL parameters not keeping alpha, beta, stagnant 

rm(list = ls())
library(dplyr) 
library(tibble)
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(writexl) 
library(readxl)  
library(purrr)  




## SIR MODEL ## 
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
    dB <- - Pb* (beta*Cb*B_cells + beta_2*Ct*B_cells)   
    dCb <-  Pb* (beta*Cb*B_cells + beta_2*Ct*B_cells) - alpha*Cb 
    dT <- -nu_a*Cb*T_cells - nu_a*Ct*T_cells  
    dAt <- nu_a*Cb*T_cells + nu_a*Ct*T_cells - beta_2*Ct*At - beta*Cb*At 
    Inf_At <- T_sus * (beta_2*Ct*At + beta*Cb*At)
    dLt <- theta*Inf_At  
    dCt <- (1-theta)*Inf_At - alpha_Ct*Ct 
    dF <- -nu_f*Lt*F
    dIf <- nu_f*Lt*F
    
    return(list(c(dB, dCb, dT, dAt, dLt, dCt, dF, dIf))) 
  }) 
}  

#### LIKELIHOOD FUNCTION #### 

Likelihood <- function(params){   
  pars <- parameters_values 
  init <- initial_values 

  pars[names(params)] <- params
  pars["beta"]  <- exp(pars["beta"]) 
  pars["beta_2"]  <- exp(pars["beta_2"])
  pars["alpha"] <- exp(pars["alpha"]) 
  pars["alpha_Ct"] <- exp(pars["alpha_Ct"])
  pars["T_sus"] <- plogis(pars["T_sus"])  
  pars["Pb"] <- plogis(pars["Pb"])
  pars["nu_f"]  <- exp(pars["nu_f"]) 

  

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
    mutate(prob_Cyto = (Cb + Ct) / (B_cells + Cb + Ct + T_cells + At + Lt), 
           prob_Cyto_B = (Cb)/(Cb+Ct), prob_Cyto_T = (Ct/(Cb+Ct)))  %>% 
    filter(time %in% obs_hourspp38) 
  
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
    left_join(infprob %>% select(time, prob_Cyto, prob_Cyto_B, prob_Cyto_T), by = "time") %>% arrange(time)   
  
  LogLHoodStor <- 0 
  
  loglike_pp38 <- dbinom(pp38_dbinom$mean.pp38, prob = pp38_dbinom$prob_Cyto, size = 40000, log = TRUE) 
  
  LogLHoodStor <- LogLHoodStor + sum(loglike_pp38) 
  
  loglike_ffe <- dnorm(
    log10(ffe_join$mean_genomes),
    mean = log10(ffe_join$If),
    sd = 1.15,
    log = TRUE
  ) 
  
  LogLHoodStor <- LogLHoodStor + sum(loglike_ffe)
  
  return(-LogLHoodStor)   
  
} 



#create intervals for numbers 

parameter_intervals <-list(
  beta  = log(c(1e-08, 1e-2)), 
  beta_2  = log(c(1e-08, 1e-2)),
  alpha = log(c(0.0104, 0.041)), 
  alpha_Ct = log(c(0.0104,0.041)),
  T_sus = qlogis(c(0.0002, 0.01116)), 
  nu_f = log(c(0.005952381, 0.0104)), 
  Pb     = qlogis(c(0.0001, 0.10))
)
#taken from baigent data  


## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c(
  beta   = log(6.951463e-07), 
  beta_2   = log(6.951463e-07),
  nu_a   = 4.668718e-01, 
  nu_f   = log(0.008),
  alpha  = log(0.0104), 
  alpha_Ct = log(0.0104),
  theta  = 0.8,
  lambda = 0.005952381, #1/(24*7)  
  T_sus  = qlogis(0.005),
  Pb     = qlogis(0.03) 
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
### TIME ### 

time_values <- seq(0, 1080, by = 1)   

#PP38 Time, Right now the time is only including some small amount 
obs_hourspp38 <- c(72,96,120,144) 
#matching to the wack hours of baigent2016
matched_time <- c(3,6,10,17,20,26,33,38)*24 


## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 
pp38_dat <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", sheet = 3, na = "NA") %>% 
  filter(!is.na(mean.pp38)) %>% select(time, mean.pp38, Bcell_no, Tcell_no)
baigent2016 <- read_xlsx("/Users/rayanhg/Downloads/baigent2016.xlsx", 2 ) %>% arrange(time)


#generate random parameters 
#it can only do one number a time
rand_number_generation <- function(rng){ 
  minimum <- rng[1]
  maximum <- rng[2]
  runif(1,min = minimum, max = maximum) 
} 

# map_dbl makes it go one element at a time [[subsets]], generating random number vector put a function within a function 
making_a_set <- function() { 
  random_parameters<- map_dbl(parameter_intervals, ~rand_number_generation(.)) 
  names(random_parameters) <- names(parameter_intervals) 
  random_parameters
}  


optim_for_alpha <- function(){ 
  
  starting_parms <- making_a_set()  
  
  answeroptim <- optim(
    par = starting_parms,
    fn  = Likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 20000)
  )
  
  
  tibble::tibble(
    Likely = answeroptim$value,
    beta   = exp(answeroptim$par["beta"]), 
    beta_2   = exp(answeroptim$par["beta_2"]),
    alpha  = exp(answeroptim$par["alpha"]), 
    alpha_Ct  = exp(answeroptim$par["alpha_Ct"]),
    T_sus  = plogis(answeroptim$par["T_sus"]), 
    nu_f   = exp(answeroptim$par["nu_f"]),
    nu_a   = parameters_values["nu_a"],
    theta  = parameters_values["theta"],
    lambda = parameters_values["lambda"], 
    Pb =  plogis(answeroptim$par["Pb"]), 
    Converged = answeroptim$convergence
    )
  
}


#how many random parameter sets I want 
n_per_alpha <-10



final_df <- purrr::map_df(1:n_per_alpha, ~optim_for_alpha()) 



write.csv(final_df, file = "/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/jan_30_26_SuperSimpleModel_AddedBaigent2016.csv") 


