#FitMDV_v01 = replaced old pp38 file to TOTAL organs, changed the scaling to 120,000 cells 
#FitMDV_v02 = scaled FFE to per 10k cells, adjusted SD to become SE 
#FitMDV_v04 = scaling is back to 120,000 cells, Scaled FFE to per 10k cells, NO PBL fit, Added fraction of B and T cell infections to likelihood fits, nu_f has lower rate interval random starts now 
#FitMDV_v06 = adding scale to FFE so that on average each cell can have 200 genomes  
#FitMDV_v07 = Fitting to individual Organs per 40k cells, FFE is still per 10k, added scaling param for FFE  
#FitMDV_v10 = Fitting to HPRS-16 data 
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
setwd("~/Desktop/WithinHostModel/DataForProject/")

## SIR MODEL ## 
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
    dB <-  -beta*Cb*B_cells - beta_2*Ct*B_cells - beta_2*Lr*B_cells + Pb*(g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct  
    dBr <- (1-Pb)*(g1*(Cb+Ct+Lr)/(g2+(Cb+Ct+Lr))) #probably should have alpha for cytolytic infection 
    dCb <- beta*Cb*B_cells + beta_2*Ct*B_cells + beta_2*Lr*B_cells - alpha*Cb 
    dT <- -nu_a*Cb*T_cells - nu_b*Ct*T_cells - nu_b*Lr*T_cells + Pt*(h1*(Cb+Ct)/(h2+(Cb+Ct))) 
    dTr <- (1 - Pt)*(h1*(Cb+Ct)/(h2+(Cb+Ct))) 
    dAt <- nu_a*Cb*T_cells + nu_b*Ct*T_cells + nu_b*Lr*T_cells - beta_2*Ct*At - beta*Cb*At - beta_2*Lr*At  # beta = rate of activated T cells by Ct and Cb cells   
    dLt <- theta *(beta_2*Ct*At + beta*Cb*At + beta_2*Lr*At)  - lambda*Lt 
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - mu*Lt5 - kappa*Lt5  
    dLr <-  kappa*Lt5 - alpha_2*Lr 
    dCt <-  (1-theta)*(beta_2*Ct*At + beta*Cb*At+ beta_2*Lr*At) - alpha_2*Ct
    dZ <- mu*Lt5
    df <- -nu_f*Lr*f
    dIf <- nu_f*Lr*f
    return(list(c(dB, dCb,dBr, dT,dTr, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dLr, dCt, dZ, df,dIf)))
  }) 
}  

#### LIKELIHOOD FUNCTION #### 

Likelihood_parts <- function(params){   
  pars <- parameters_values 
  init <- initial_values 
  B0 <-2.4e9/3  
  T0 <- 7.4e8/3 
  
  pars[names(params)] <- params  
  pars["beta"]   <- exp(pars["beta"])
  pars["beta_2"] <- exp(pars["beta_2"])  
  pars["alpha"]   <- exp(pars["alpha"])
  pars["alpha_2"] <- exp(pars["alpha_2"])   
  pars["kappa"] <- exp(pars["kappa"])  
  pars["nu_a"] <- exp(pars["nu_a"])
  pars["nu_b"] <- exp(pars["nu_b"]) 
  pars["nu_f"] <- exp(pars["nu_f"])
  pars["size_pp38"] <- exp(pars["size_pp38"]) 
  pars["size_pp38_Ct"] <- exp(pars["size_pp38_Ct"])
  
  
  
  init["B_cells"] <- B0*pars["Pb"] 
  init["Br"] <- B0* (1 - pars["Pb"]) 
  
  init["T_cells"] <- T0*pars["Pt"] 
  init["Tr"] <- T0* (1 - pars["Pt"])
  
  # run modelop
  results <- as.data.frame(
    ode(y = init, 
        times = time_values, 
        func = sir_equations, 
        parms = pars)
  ) 
  
  ### MODEL ###  
  
  # percent cytolytically infected in model
  infprob <- results %>% 
    dplyr::select(B_cells, Cb,Br, Ct, T_cells, At, Lt, Lt2, Lt3, Lt4, Lt5,Lr,Tr, time) %>% 
    mutate(
      prob_Cyto = (Cb + Ct + Lr) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5 + Lr + Tr),
      inf_cells = ((Cb + Ct + Lr) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5 + Lr + Tr))*40000,
      inf_cells_Cb = ((Cb) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5 + Lr + Tr))*40000,  
      inf_cells_Ct = ((Ct + Lr) / (B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5 + Lr + Tr))*40000
    ) %>% filter(time %in% obs_hourspp38) 
  
  
 
  pp38_dbinom_Cb <- pp38_dat %>% filter(time %in% obs_hourspp38) %>% 
    left_join(infprob %>% dplyr::select(time, inf_cells_Cb), by = "time") 
  
  pp38_dbinom_Ct <- pp38_dat %>% filter(time %in% obs_hourspp38) %>% 
    left_join(infprob %>% dplyr::select(time, inf_cells_Ct), by = "time") 
  
  loglike_pp38 <- dnbinom(pp38_dbinom_Cb$Bcell_no,# from data, observed 
                          size = pars["size_pp38"], prob = pars["size_pp38"]/(pars["size_pp38"]+pp38_dbinom_Cb$inf_cells_Cb),  log = TRUE) # cytolytic cells from model 
  
  loglike_pp38_Ct <- dnbinom(pp38_dbinom_Ct$Tcell_no,# from data, observed 
                             size = pars["size_pp38_Ct"], prob = pars["size_pp38_Ct"]/(pars["size_pp38_Ct"]+pp38_dbinom_Ct$inf_cells_Ct),  log = TRUE) # cytolytic cells from model 
  
  
  sum_loglike_pp38 <- sum(loglike_pp38) + sum(loglike_pp38_Ct)
  
  ### Scaling Plaque assay Data ### 
  PlaqueCells <- results %>% mutate(pfuCells = (Cb + Ct + Lr)/(B_cells + Cb + Ct + T_cells + At + Br + Lt + Lt2 + Lt3 + Lt4 + Lt5 + Lr + Tr) * 10000000) %>% 
    filter(time %in% powell1982_plaque$time)
  
  
  ### Adding Plaque Assay Data ###  
  
  loglike_Plaque1982 <- dnorm(x = powell1982_plaque$pfuPer1e7, 
                              mean = PlaqueCells$pfuCells,
                              sd = powell1982_plaque$SE, log = TRUE ) 
  
  sum_loglike_PFU <- sum(loglike_Plaque1982)
  
  LogLHoodStor <-  sum_loglike_pp38 + sum_loglike_PFU
  
  #i am making the log values negative so the smaller the log value the better
  list(neg_LogLHoodStor = -LogLHoodStor,
       neg_sum_loglike_pp38 =  -sum_loglike_pp38, 
       neg_sum_loglike_PFU = -sum_loglike_PFU) 
  
} 



Likelihood <- function(params){  
  
  Likelihood_parts(params)$neg_LogLHoodStor
}

#create intervals for numbers  


parameter_intervals <- list( beta_2 = log(c(1e-10,2e-8)), 
                             beta = log(c(1e-10, 2e-8)),              
                             alpha_2 =log(c(0.008,0.014)),          # assuming slower death rate than B cells  
                             g1 = c(50,300),
                             g2 = c(300,800),
                             h1 = c(559,600),
                             h2 = c(250,800),
                             nu_a = log(c(1e-8, 3e-8)),             
                             nu_b = log(c(1e-8, 3e-8)),             
                             alpha = log(c(0.01,0.04)),             # bursa‐dependent subpopulations of peripheral B lymphocytes in chicken blood
                             # mu = c(0.01388889, 0.05), 
                             nu_f = log(c(1e-10, 1e-8)),            # previous nu_f = log(c(0.005952381, 0.0104))
                             size_pp38 = log(c(0.05,0.5)),          # proposed fitting range based off data: 0.10 - 0.50, approx = 0.20 , 0ld = 0.05, 0.5
                             size_pp38_Ct = log(c(0.05,0.5)),       # proposed fitting range based off data: 0.10 - 0.50, approx = 0.17 , 0ld = 0.05 , 0.5 
                             kappa = log(c(1e-5, 1.3e-4)))          # Temperature-induced reactivation of Marek's disease virus-transformed T cells ex vivo   


## PARAMETERS AND INITIAL VALUES ## 
parameters_values <- c( 
  beta =  log(6.951463e-07)                  #contact rate with B cells  
  , Pb = 0.01 
  , Pt = 0.01
  , beta_2 = log(8.532282e-07)                 #contact rate with T cells 
  , nu_a = log(4.668718e-01)             #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = log(4.467098e-01)             #Activation rate of T cells by cytolytic T cells (hours)
  , nu_f = log(1e-11)                    #Infection rate of follicular cells (hours) 
  , kappa = log(1e-6)
  , mu = 1/8                       #Rate of Tumor Cells (every 72 hours)
  , alpha = log(0.0104)             #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = log(0.0104)           #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 10                     #incoming B cells (every 15 hours)
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
  , Tr = 0 
  , At = 0 
  , Lt = 0 
  , Lt2 = 0 
  , Lt3 = 0  
  , Lt4 = 0  
  , Lt5 = 0  
  , Lr = 0 
  , Ct = 0
  , Z = 0
  , f = 400000 
  , If =0 
  
) 
### TIME ### 

time_values <- seq(0, 1440, by = 1)   

#PP38 Time, Right now the time is only including some small amount 
obs_hourspp38 <- c(72,96,120,144) # make sure pp38 times is matching the data loaded in 


## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 

pp38_dat <-  read_xlsx("Baigent1998/baigent1998.xlsx", sheet = 3, na = "NA") %>% filter(!is.na(mean.pp38)) %>% 
  dplyr::select(time, Bcell_no,Tcell_no)

powell1982_plaque <- read_xlsx("Powell1982/Powell1982.xlsx", sheet = 2)%>% mutate(time = dpi*24) 

# ggplot(data = powell1982_plaque, aes(x = time/24, y = pfuPer1e7)) + 
#   geom_point() + scale_y_log10() + geom_line() + 
#   geom_errorbar(aes(ymin = SE_upper, ymax = SE_lower)) 

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


lower_bounds <- map_dbl(
  parameter_intervals,
  ~ .x[1]
)

upper_bounds <- map_dbl(
  parameter_intervals,
  ~ .x[2]
)





optim_for_alpha <- function() {
  starting_parms <- making_a_set()
  
  tryCatch({
    answeroptim <- optim(
      par = starting_parms,
      fn  = Likelihood,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = 20000)
      # method = "Nelder-Mead",
      # control = list(maxit = 20000)
    ) 
    
    parts <- Likelihood_parts(answeroptim$par)
    
    tibble::tibble(
      Likely   = answeroptim$value, 
      Likely_pp38 = parts$neg_sum_loglike_pp38, 
      Likely_PFU = parts$neg_sum_loglike_PFU,
      beta     = exp(answeroptim$par["beta"]),
      beta_2   = exp(answeroptim$par["beta_2"]),
      alpha    = exp(answeroptim$par["alpha"]),
      alpha_2  = exp(answeroptim$par["alpha_2"]), 
      Pb       = parameters_values["Pb"], 
      Pt       = parameters_values["Pt"], 
      nu_a     = exp(answeroptim$par["nu_a"]),
      nu_b     = exp(answeroptim$par["nu_b"]),
      nu_f     = exp(answeroptim$par["nu_f"]), 
      kappa    = exp(answeroptim$par["kappa"]), 
      mu       = parameters_values["mu"],
      g1       = answeroptim$par["g1"],
      g2       = answeroptim$par["g2"],
      h1       = answeroptim$par["h1"],
      h2       = answeroptim$par["h2"],
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
      Likely_PFU = NA_real_, 
      beta     = NA_real_,
      beta_2   = NA_real_,
      alpha    = NA_real_,
      alpha_2  = NA_real_, 
      kappa    = NA_real_, 
      nu_a     = NA_real_,
      nu_b     = NA_real_,
      nu_f     = NA_real_,
      mu       = NA_real_,
      g1       = NA_real_,
      g2       = NA_real_,
      h1       = NA_real_,
      h2       = NA_real_,
      Pb       = NA_real_, 
      Pt       = NA_real_, 
      size_pp38 = NA_real_, 
      size_pp38_Ct = NA_real_,
      Converged = NA_integer_,
      ErrorMsg = e$message
      
    )
  })
}



n_per_alpha <- 100  # change this to 10, 50, 100, etc.

available_workers <- max(
  1,
  parallelly::availableCores(omit = 1)
)

workers <- min(n_per_alpha, available_workers)

plan(multisession, workers = workers) 

cat("Running", n_per_alpha, "random starts using", workers, "parallel workers\n")


# run many independent model fits in parallel
results_list <- future_lapply(
  1:n_per_alpha,
  function(i) {
    optim_for_alpha()
  },
  future.seed = TRUE
)

# combine all results into one table
final_df <- bind_rows(results_list)


# save the file with a timestamp so it does not overwrite old runs
out_file <- file.path(
  "~/Desktop/WithinHostModel/CodeOutputsRandNum/",
  paste0("FitMDV_v10_HPRS16_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
)

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

cat("Saved results to:", out_file, "\n")

# optional: turn off parallel workers when finished
plan(sequential) 


end_time <- Sys.time() - start_time 

print(end_time)
