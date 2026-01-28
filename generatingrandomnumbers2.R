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
    dZ <- mu*Lt5
    df <- -nu_f*Lt5*f
    dIf <- nu_f*Lt5*f
    return(list(c(dB, dCb, dT, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dCt, dZ, df,dIf)))
  }) 
}  


#### LIKELIHOOD FUNCTION #### 

Likelihood <- function(params){   
  pars <- parameters_values
  
  pars[names(params)] <- abs(params)   
  
  # run model
  results <- as.data.frame(
    ode(y = initial_values, 
        times = time_values, 
        func = sir_equations, 
        parms = pars)
  ) 
  
  ### MODEL ###  
  
  # percent cytolytically infected in model
  infprob <- results %>% 
    dplyr::select(B_cells, Cb, Ct, T_cells, At, Lt, Lt2, Lt3, Lt4, Lt5, time) %>% 
    mutate(
      prob_Cyto = (Cb + Ct) / (B_cells + Cb + Ct + T_cells + At) # getting rid of Lt because Lt is not even present at this time point
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
  

  
  LogLHoodStor <- 0 
  
  pp38_dbinom <- pp38_dat %>%
    left_join(infprob %>% select(time, prob_Cyto), by = "time")
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
parameter_intervals <- list( beta_2 = c(1e-08,1e-2), 
                             beta = c(1e-08, 1e-2), 
                             alpha_2 = c(0.0104, 0.041), 
                             g1 = c(10,1000), 
                             g2 = c(1,1000), 
                             h1 = c(10,1000), 
                             h2 = c(1,1000), 
                             alpha = c(0.0104,0.041), 
                             mu = c(0.01388889, 0.05), 
                             nu_f = c(0.006, 0.01), 
                             lambda = c(0.015, 0.042))  

## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  M = 0
  , beta =  6.951463e-07                  #contact rate with B cells 
  , beta_2 = 8.532282e-07                  #contact rate with T cells 
  , nu_a = 4.668718e-01                     #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = 4.467098e-01                   #Activation rate of T cells by cytolytic T cells (hours)
  , nu_f = 0.008                     #Infection rate of follicular cells (hours)
  , mu = 0.02                        #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.0104                   #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.0104                 #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 1000                     #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 =1000                      #incoming T cells / determined no incoming T cells 
  , h2 = 100
  , lambda = 0.03                #adding delay, how long latent cell 'exposed' 
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

### TIME ### 

time_values <- seq(0, 1080, by = 1)   
#matching to the wack hours of baigent2016
matched_time <- c(3,6,10,17,20,26,33,38)*24 

#PP38 Time, Right now the time is only including some small amount 
obs_hourspp38 <- c(72,96,120,144) 


## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 
pp38_dat <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", sheet = 3, na = "NA") %>% filter(!is.na(mean.pp38)) %>% select(time, mean.pp38)
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
    control = list(maxit = 5000)
  )
  
  tibble::tibble(
    Likely   = answeroptim$value,
    beta     = answeroptim$par["beta"],
    beta_2   = answeroptim$par["beta_2"],
    alpha    = answeroptim$par["alpha"],
    alpha_2  = answeroptim$par["alpha_2"], 
    nu_f     = answeroptim$par["nu_f"], 
    lambda   = answeroptim$par["lambda"],
    mu       = answeroptim$par["mu"],
    g1       = answeroptim$par["g1"],
    g2       = answeroptim$par["g2"],
    h1       = answeroptim$par["h1"],
    h2       = answeroptim$par["h2"],
    Converged = answeroptim$convergence
  ) 
  
  
}


#how many random parameter sets I want 
n_per_alpha <-50



final_df <- purrr::map_df(1:n_per_alpha, ~optim_for_alpha())

#adding my alpha values to df 

# data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Random_parameter_exploration_1_NA.RM.csv") %>% filter(Converged != 10) %>%
#   mutate(alpha = abs(alpha), Likely = abs(Likely), beta = abs(beta), beta_2 = abs(beta_2), alpha_2 = abs(alpha_2), g1 = abs(g1), g2 = abs(g2),
#          h1 = abs(h1), h2 = abs(h2)) %>% select(-X, -Converged) %>% filter(Likely < 500)


write.csv(final_df, file = "/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Random_parameter_exploration_with_LT.csv") 

#, append = TRUE, sep = "", col.names = FALSE 

# plot(log10(data))
# ggpairs(data) 

