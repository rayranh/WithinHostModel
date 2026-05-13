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
library(future) 
library(future.apply) 


## SIR MODEL ## 
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { 
    dCb <- beta*Cb*B_cells       #assume over 50 hours that B cells do not die 
  
    return(list(c(dCb)))
  }) 
}  

#### LIKELIHOOD FUNCTION #### 

Likelihood_parts <- function(params){   
  pars <- parameters_values 
  init <- initial_values  
  
  pars[names(params)] <- params
  
  pars["beta"]   <- exp(pars["beta"])
  
  # run model
  results <- as.data.frame(
    ode(y = init, 
        times = time_values, 
        func = sir_equations, 
        parms = pars)
  ) 
  
  ### MODEL ###  
  # baigent2016 times are in days; matched_time is hours 
  ffe_dat <- BCells2015
  
  #matching model time to feathers time 
  Infected_FFE<- results %>% dplyr::select("time", "Cb") %>% 
    filter(time %in% matched_time) %>% arrange(time)  
  
  # join model predictions (If) to observed data by time
  ffe_join <- ffe_dat %>%
    left_join(Infected_FFE, by = "time") %>% 
    mutate(sd_log = (log10(uppersd) - log10(lowersd))/2)# brings in If
  
  
  loglike_ffe <- dnorm(
    log10(ffe_join$MDV_load),
    mean = log10(ffe_join$Cb),
    sd = ffe_join$sd_log, # calculated from 17dpi onwards because it is stable at that point and anything before then would drag up the estimation 
    log = TRUE
  )  
  
  sum_loglike_ffe <- sum(loglike_ffe) 
  
  list(neg_sum_loglike_ffe = -sum_loglike_ffe, 
       neg_LogLHoodStor = -sum_loglike_ffe)
  
} 



Likelihood <- function(params){  
  
  Likelihood_parts(params)$neg_LogLHoodStor
}


## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  beta =  log(6.951463e-12),
  B_cells = 1e6#contact rate with B cells  
)


initial_values <- c( 
  Cb = 1) 
### TIME ### 

time_values <- seq(0, 48, by = 1) 
matched_time <- c(4,24,48)

## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 

BCells2015 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Schermuly2015/schermuly_Bcells_fin.csv")%>% 
  mutate(MDV_load = MDV_load/100, uppersd = uppersd/100, lowersd = lowersd/100) 



# map_dbl makes it go one element at a time [[subsets]], generating random number vector put a function within a function  

starting_parms <- c(beta =log(6.951463e-6))

  
  answeroptim <- optim(
    par = starting_parms,
    fn  = Likelihood,
    method = "Brent",
    lower = log(1e-15),
    upper = log(1e-6)
  ) 
  
  parts <- Likelihood_parts(answeroptim$par)
  
optimized_param <-  tibble::tibble(
    Likely      = answeroptim$value, 
    beta     = exp(answeroptim$par[1]),
    Converged = answeroptim$convergence
  )


parameters_values["beta"] <- optimized_param$beta
# run model
results <- as.data.frame(
  ode(y = initial_values, 
      times = time_values, 
      func = sir_equations, 
      parms = parameters_values)
) 



ggplot(data = BCells2015, aes(x = time, y = MDV_load, color = "data")) + geom_line() + 
  scale_y_log10() + geom_line(data = results, aes(x = time, y =Cb, colour = "model"), inherit.aes = FALSE) + 
  scale_color_manual(values = c("data" = "black", "model" = "red"))

# write.table(
#   final_df,
#   file = "~/scratch/Feb.17.26.FittingDnbinom_mu_theta.csv",
#   sep = ",",
#   row.names = FALSE,
#   col.names = FALSE,
#   append = TRUE)
