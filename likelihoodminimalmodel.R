rm(list = ls())
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(dplyr) 
library(patchwork) 
library(writexl) 
library(readxl)  
library(purrr) 

sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations 
    dB <- -M*B_cells - beta*Cb*B_cells - beta_2*Ct*B_cells + (g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct 
    dCb <- M*B_cells +beta*Cb*B_cells + beta_2*Ct*B_cells - alpha*Cb 
    dT <- -M*T_cells - nu_A*Cb*T_cells - nu_b*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct)))
    dAt <- M*T_cells + nu_A*Cb*T_cells + nu_b*Ct*T_cells - beta_2*Ct*At - beta*Cb*At   #remove mu*Lt since only the last compartment should contribute to tumors  
    dLt <- theta*(beta_2*Ct*At + beta*Cb*At)   - lambda*Lt 
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - lambda*Lt5 - mu*Lt5
    dCt <- (1-theta)*(beta_2*Ct*At + beta*Cb*At) - alpha_2*Ct
    dZ <- mu*Lt5
    df <- -nu_F*Lt5*f
    dIf <- nu_F*Lt5*f
    return(list(c(dB, dCb, dT, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dCt, dZ, df,dIf)))
  }) 
} 


Likelihood <- function(params){   
  pars <- parameters_values 
  pars["nu_f"] <- abs(params[1])  
  pars["beta"] <- abs(params[2]) 
  pars["beta_2"] <- abs(params[3]) 
  pars["nu_A"] <- abs(params[4]) 
  pars["nu_B"] <- abs(params[5]) 
  pars["mu"] <- abs(params[6])
  pars["alpha"] <- abs(params[7])
  pars["alpha2"] <- abs(params[8]) 
  pars["g1"] <- abs(params[9])
  pars["g2"] <- abs(params[10])
  pars["h1"] <- abs(params[11])
  pars["h2"] <- abs(params[12])
  pars["lambda"] <- abs(params[13])
  
  
  # 
  #creating model and replacing with new beta values 
  results <- as.data.frame(ode(y = initial_values, 
                               times = time_values, func = sir_equations, parms = pars))
  
  ### MODEL ### 

  #match model time with observed data 
  Infected_FFE<- results %>% select("time", "If") %>% filter(time %in% matched_time)   # changing to %in% because apparently better   
  #empty list 
  #result <- vector(mode = "list", length = nrow(dtotallike))   
  LogLHoodStor=0 
  
  for (i in seq_len(length(matched_time))) { 
    obs_vec_B <- Infected_FFE %>% filter(time %in% Infected_FFE$time[i]) %>% pull(If) 
    mu_B <- baigent2016$mean_genomes[i]
    loglike_B <- dnorm(log(obs_vec_B), log(mu_B), sd = 1, log = TRUE)  
    LogLHoodStor=LogLHoodStor+sum(loglike_B)
    #result[[i]] <- data.frame(time = rep(dtotallike$time[i], length(obs_vec)) , obs = obs_vec, model_val = rep( mu, length(obs_vec)) , likelihood = loglike) 
  }
  

  return(-LogLHoodStor)
}

PBL_data <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/PayneRennie1976raw.xlsx")
baigent2016 <- read_xlsx("/Users/rayanhg/Downloads/baigent2016.xlsx", 2)



#blah blah adding someing
### NOTE: 9/14/2025 used this model for presentation for collaborators to illustrate adding delay to model improves fit ### 


#get rid of macrophage and infect T cells and B cells population - could also have macrophages decay over time 
#blah blah adding someing  
# get rid of macrophages and add cB 
parameters_values <- c( 
  M = 0
  , beta = 7.108697e-7                 #contact rate with B cells 
  , beta_2 = 2.239485e-7                   #contact rate with T cells 
  , nu_A = 0.05                     #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = 0.001                    #Activation rate of T cells by cytolytic T cells (hours)
  , nu_F =35.380823e-02                      #Infection rate of follicular cells (hours)
  , mu = 0.5                        #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.55                   #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.01                  #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 10000                      #incoming B cells (every 15 hours)
  , g2 = 5    
  , h1 =100                        #incoming T cells / determined no incoming T cells 
  , h2 = 5  
  , lambda = 0.000055                #adding delay, how long latent cell 'exposed' 
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

time_values <- seq(0, 1000) # hours

sir_values_1 <- ode(
  y = initial_values, # variables = initial values 
  times = time_values,
  func = sir_equations, # ode is putting in the variables and initial values for me 
  parms = parameters_values 
) 

#added strings as factors because otherwise turns into factors for some reason 
sir_values_1 <- as.data.frame(sir_values_1, stringsAsFactors = FALSE) 
#turning it longways 
df <- melt(sir_values_1, id.vars = "time") %>% filter(!variable %in% c("Lt", "Lt2", "Lt4", "Lt3")) 



time_values <- seq(0, 1000) # hours 

#matching to the wack hours of baigent2016
matched_time <- c(3,6,10,17,20,26,33,38)*24


test_parms <- c(9.380823,7.108697e-7,2.239485e-7,0.05, 0.08,0.5,0.55,0.55,100,10,100,10,0.000055)



Likelihood(test_parms)

#### remember to transform any exp() numbers back 
optim(par = test_parms, fn = Likelihood, method = "Nelder-Mead",control = list(maxit = 5000)) 


