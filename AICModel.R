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

start_time <- Sys.time()


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



Likelihood <- function(params){  
  
  Likelihood_parts(params)$neg_LogLHoodStor
}

#create intervals for numbers 

parameter_intervals <- list( beta_2 = log(c(1e-14,1e-9)), 
                             beta = log(c(1e-14, 1e-9)), 
                             alpha_2 =log(c(0.0104, 0.041)), 
                             g1 = c(10,1000), 
                             g2 = c(1,1000), 
                             h1 = c(10,1000), 
                             h2 = c(1,1000),  
                             nu_a = log(c(0.013,0.4)),             #Activation rate of T cells by cytolytic B cells (hours)
                             nu_b = log(c(0.013,0.4)),            #Activation rate of T cells by cytolytic T cells (hours)
                             alpha = log(c(0.0104,0.041)),  
                             mu = c(0.01388889, 0.05), 
                             nu_f = log(c(0.005952381, 0.0104)), 
                             Pb = qlogis(c(0.001, 0.05))) #taken from baigent pp38 data 


## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  beta =  log(6.951463e-07)                  #contact rate with B cells  
  , Pb = qlogis(0.005) 
  , beta_2 = log(8.532282e-07)                 #contact rate with T cells 
  , nu_a = log(4.668718e-01)             #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = log(4.467098e-01)             #Activation rate of T cells by cytolytic T cells (hours)
  , nu_f = log(0.008)                    #Infection rate of follicular cells (hours)
  , mu = 0.02                       #Rate of Tumor Cells (every 72 hours)
  , alpha = log(0.0104)             #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = log(0.0104)           #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 10                       #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 =10                      #incoming T cells / determined no incoming T cells 
  , h2 = 100
  , lambda = 0.02380952             # fixing delay rate to (1/(7*24))*4  
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
pp38_dat <- read_xlsx("baigent1998.xlsx", sheet = 3, na = "NA") %>% filter(!is.na(mean.pp38)) %>% select(time, mean.pp38)
baigent2016 <- read_xlsx("baigent2016.xlsx", 2 ) %>% arrange(time)
AIC_RandomStarts <- read.csv("AIC_RandomStarts.csv") 


optim_for_alpha <- function(params, i) {
  starting_parms <- unlist(params[i,])
  
  tryCatch({
    answeroptim <- optim(
      par = starting_parms,
      fn  = Likelihood,
      method = "Nelder-Mead",
      control = list(maxit = 20000)
    ) 
    
    parts <- Likelihood_parts(answeroptim$par)
    
    tibble::tibble(
      Likely   = answeroptim$value, 
      Likely_pp38 = parts$neg_sum_loglike_pp38, 
      Likely_ffe = parts$neg_sum_loglike_ffe, 
      beta     = exp(answeroptim$par["beta"]),
      beta_2   = exp(answeroptim$par["beta_2"]),
      alpha    = exp(answeroptim$par["alpha"]),
      alpha_2  = exp(answeroptim$par["alpha_2"]),
      nu_a     = exp(answeroptim$par["nu_a"]),
      nu_b     = exp(answeroptim$par["nu_b"]),
      nu_f     = exp(answeroptim$par["nu_f"]),
      mu       = answeroptim$par["mu"],
      g1       = answeroptim$par["g1"],
      g2       = answeroptim$par["g2"],
      h1       = answeroptim$par["h1"],
      h2       = answeroptim$par["h2"],
      Pb       = plogis(answeroptim$par["Pb"]),
      Converged = answeroptim$convergence,
      ErrorMsg = NA_character_
    )
    
  }, error = function(e) {
    # optional: print to console
    message("Error with starting_parms: ", paste(names(starting_parms), starting_parms, collapse = ", "))
    message("Error message: ", e$message)
    
    # return the starting params so you can inspect them later
    tibble::tibble(
      Likely   = Inf, 
      Likely_pp38 = NA_real_,
      Likely_ffe  = NA_real_,
      beta     = NA_real_,
      beta_2   = NA_real_,
      alpha    = NA_real_,
      alpha_2  = NA_real_,
      nu_a     = NA_real_,
      nu_b     = NA_real_,
      nu_f     = NA_real_,
      mu       = NA_real_,
      g1       = NA_real_,
      g2       = NA_real_,
      h1       = NA_real_,
      h2       = NA_real_,
      Pb       = NA_real_,
      Converged = NA_integer_,
      ErrorMsg = e$message
      
    )
  })
}



library(future)
library(future.apply)
library(dplyr)

# get number of CPUs Slurm gave us
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

plan(multisession, workers = workers)

results_list <- future_lapply(
  1:nrow(AIC_RandomStarts),
  function(i) optim_for_alpha(AIC_RandomStarts,i),
  future.seed = TRUE
) 

final_df <- bind_rows(results_list)


out_file <- "~/scratch/Feb.17.26.FittingDnbinom_10Parms_AIC.csv"

write.table(
  final_df,
  file = out_file,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  append = FALSE
)



end_time <- Sys.time() - start_time
print(end_time)