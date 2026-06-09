#Using this code to generate random numbers for ALL parameters not keeping alpha, beta, stagnant 

#rm(list = ls())
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
#setwd("~/Desktop/WithinHostModel/DataForProject/")

start_time <- Sys.time()

## SIR MODEL ## 
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

#### LIKELIHOOD FUNCTION #### 
Likelihood_parts <- function(params){   
  pars <- parameters_values 
  init <- initial_values 
  B0 <-2.4e9/3
  
  pars[names(params)] <- params  
  pars["beta"]   <- exp(pars["beta"])
  pars["alpha"]   <- exp(pars["alpha"])
  pars["alpha_2"] <- exp(pars["alpha_2"])  
  pars["nu_a"] <- exp(pars["nu_a"])
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
  
  #matching model time to feathers time 
  Infected_FFE<- results %>% dplyr::select("time", "If", "f") %>% 
    filter(time %in% matched_time) %>% arrange(time) %>% mutate(If_per10k = (If/(If + f))*10000) 
  
  # join model predictions (If) to observed data by time
  ffe_join <- baigent2016 %>%
    left_join(Infected_FFE, by = "time") 
  
  
  pp38_dbinom_Cb <- pp38_dat %>% filter(time %in% obs_hourspp38) %>% 
    left_join(infprob %>% select(time, inf_cells_Cb), by = "time") 
  
  pp38_dbinom_Ct <- pp38_dat %>% filter(time %in% obs_hourspp38) %>% 
    left_join(infprob %>% select(time, inf_cells_Ct), by = "time") 
  
  loglike_pp38 <- dnbinom(pp38_dbinom_Cb$pp38BcellTotal,# from data, observed 
                          size = pars["size_pp38"], prob = pars["size_pp38"]/(pars["size_pp38"]+pp38_dbinom_Cb$inf_cells_Cb),  log = TRUE) # cytolytic cells from model 
  
  loglike_pp38_Ct <- dnbinom(pp38_dbinom_Ct$pp38TcellTotal,# from data, observed 
                             size = pars["size_pp38_Ct"], prob = pars["size_pp38_Ct"]/(pars["size_pp38_Ct"]+pp38_dbinom_Ct$inf_cells_Ct),  log = TRUE) # cytolytic cells from model 
  
  
  sum_loglike_pp38 <- sum(loglike_pp38) + sum(loglike_pp38_Ct) 
  
  
  loglike_ffe <- dnorm(
    ffe_join$logged10Mean, #log10() here to match the log10 of the data and keep scale true
    mean = log10(ffe_join$If_per10k),
    sd =  0.352, # Avg SE calculated from 17dpi onwards 
    log = TRUE # outer log is TRUE so it is returning ln() to keep consistent
  )   
  
  sum_loglike_ffe <- sum(loglike_ffe)  
  
  PBLModelInfection <- results %>% 
    mutate(genomes_pbl = (Cb + Ct + Lt + Lt2 + Lt3 + Lt4 + Lt5)/(Cb+Ct+B_cells+T_cells+At+Lt+Lt2+Lt3+Lt4+Lt5+Br)*10000) %>%
    filter(time %in% matched_time) %>% select(time,Cb, Ct, Lt, Lt2, Lt3, Lt4, Lt5, genomes_pbl)
  
  PBL_join <- baigentpbl2016 %>% left_join(PBLModelInfection, by = "time")
  
  loglike_pbl <- dnorm(
    PBL_join$logged10Mean, 
    mean = log10(PBL_join$genomes_pbl), 
    sd = 0.365, 
    log = TRUE) 
  
  sumloglike_pbl <- sum(loglike_pbl) 
  
  
  ### FOR SPATZ 2007 ask about Data #### 
  
  
  
  LogLHoodStor <- sum_loglike_ffe + sum_loglike_pp38 + sumloglike_pbl
  
  list(neg_LogLHoodStor = -LogLHoodStor,
       neg_sum_loglike_ffe = -sum_loglike_ffe, 
       neg_sum_loglike_pp38 =  -sum_loglike_pp38, 
       neg_sum_loglike_pbl = -sumloglike_pbl)
} 



Likelihood <- function(params){  
  
  Likelihood_parts(params)$neg_LogLHoodStor
}

parameter_intervals <- list(beta = log(c(1e-14, 1e-5)), 
                            alpha_2 =log(c(0.0104, 0.041)), 
                            # g1 = c(10,1000),
                            # g2 = c(1,1000),
                            # h1 = c(10,1000),
                            # h2 = c(1,1000),
                            alpha = log(c(0.0034722,0.02)),  # dont say death rate 
                            # mu = c(0.01388889, 0.05), 
                            nu_f = log(c(0.005952381, 0.0104)), 
                            Pb = qlogis(c(0.001, 0.05)),
                            size_pp38 = log(c(0.05,10)),
                            size_pp38_Ct = log(c(0.05,10)))  #taken from baigent pp38 data 


## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  beta =  log(6.951463e-07)                  #contact rate with B cells  
  , Pb = qlogis(0.005) 
  , nu_a = log( 0.1666667)             #Activation rate of T cells by cytolytic B cells (hours)
  , nu_f = log(0.008)                    #Infection rate of follicular cells (hours)
  , mu = 1/8                       #Rate of Tumor Cells (every 72 hours)
  , alpha = log(0.0104)             #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = log(0.0104)           #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 0                    #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 = 0                      #incoming T cells / determined no incoming T cells 
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
  , f = 4e5 # (Calnek et al., 1970) range of 1.5e5 - 7e5  
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
pp38_dat <- read_xlsx("~/work/baigent1998Totalpp38.xlsx") 
baigent2016 <- read_xlsx("~/work/feathers_noVax_SE.xlsx") 
baigentpbl2016 <- read_xlsx( "~/work//PBL_noVax_SE.xlsx")


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

optim_for_alpha <- function() {
  starting_parms <- making_a_set()
  
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
      Likely_pbl = parts$neg_sum_loglike_pbl,
      beta     = exp(answeroptim$par["beta"]),
      alpha    = exp(answeroptim$par["alpha"]),
      alpha_2  = exp(answeroptim$par["alpha_2"]), 
      nu_a = exp(parameters_values["nu_a"]), 
      nu_f     = exp(answeroptim$par["nu_f"]),
      mu       = parameters_values["mu"],
      # g1       = answeroptim$par["g1"],
      # g2       = answeroptim$par["g2"],
      # h1       = answeroptim$par["h1"],
      # h2       = answeroptim$par["h2"],
      Pb       = plogis(answeroptim$par["Pb"]), 
      size_pp38 = exp(answeroptim$par["size_pp38"]), 
      size_pp38_Ct = exp(answeroptim$par["size_pp38_Ct"]), 
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
      Likely_pbl = NA_real_,
      beta     = NA_real_,
      alpha    = NA_real_,
      alpha_2  = NA_real_,
      nu_a     = NA_real_,
      nu_f     = NA_real_,
      mu       = NA_real_,
      # g1       = NA_real_,
      # g2       = NA_real_,
      # h1       = NA_real_,
      # h2       = NA_real_,
      Pb       = NA_real_, 
      size_pp38 = NA_real_,
      size_pp38_Ct = NA_real_,
      Converged = NA_integer_,
      ErrorMsg = e$message
      
    )
  })
}





#how many random parameter sets I want 
n_per_alpha <-10


library(future)
library(future.apply)
library(dplyr)

# get number of CPUs Slurm gave us
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
plan(multisession, workers = workers)
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")) #took out sys.tim() as set seed because array can start at same second
seed <- as.integer(Sys.time()) + task
set.seed(seed)

results_list <- future_lapply(
  1:n_per_alpha,
  function(i) optim_for_alpha(),
  future.seed = TRUE
)

final_df <- bind_rows(results_list) 


out_file <- sprintf("~/work/Jun.4.26.FittingMDVModelCluster_%03d.csv", task) #outputting separate files

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