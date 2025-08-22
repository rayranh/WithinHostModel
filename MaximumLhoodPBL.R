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
  with(as.list(c(variables, parameters)), {
    dB <- -M*B_cells - beta*Cb*B_cells - beta_2*Ct*B_cells + (g1*(Cb+Ct)/(g2+(Cb+Ct))) - mu_o*B_cells + mu_p*Bb_cells # beta = rate of cytolytically infected B cells by Cb and Ct 
    dCb <- M*B_cells +beta*Cb*B_cells + beta_2*Ct*B_cells - alpha*Cb - mu_o*Cb+ mu_p*Cbb_cells 
    dCt <- (1-theta)*(beta_2*Ct*At + beta*Cb*At) - alpha_2*Ct - mu_t*Ct + mu_p*Ctb_cells
    dT <- -M*T_cells - nu_A*Cb*T_cells - nu_b*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct))) - mu_t*T_cells + mu_p*Tb_cells
    dAt <- M*T_cells + nu_A*Cb*T_cells + nu_b*Ct*T_cells - beta_2*Ct*At - beta*Cb*At - mu_t*At + mu_p*Atb_cells  # beta = rate of activated T cells by Ct and Cb cells 
    dLt <- theta*(beta_2*Ct*At + beta*Cb*At) - lambda*Lt  #remove mu*Lt since only the last compartment should contribute to tumors  
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - lambda*Lt5 - mu_t*Lt5 
    
    ### entering blood ###   
    dBb <- mu_o*B_cells - mu_p*Bb_cells - alpha_B*Bb_cells   
    dCbb <- mu_o*Cb - mu_p*Cbb_cells - alpha_B*Cbb_cells   # removed mu_o B cells b/c lung dynamics not being modelled 
    dTb <- mu_t*T_cells  - mu_p*Tb_cells - alpha_B*Tb_cells  
    dAtb <- mu_t*At  - mu_p*Atb_cells # no death for At cells 
    dLtb <- mu_t*Lt5 - mu_p*Ltb_cells  # population of Lt?  
    dCtb <- mu_t*Ct - mu_p*Ctb_cells    
    
    
    df <- - nu_F*Ltb_cells*f 
    dIf <- nu_F*Ltb_cells*f 
    return(list(c(dB, dCb, dT, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dBb, dCbb, dTb, dAtb, dLtb, dCtb, dCt, df,dIf)))
  }) 
}

Likelihood <- function(params){   
  pars <- parameters_values 
  pars["beta"] <- exp(params[1])  
  pars["alpha"] <-exp(params[2]) 
  pars["beta2"] <- params[3] 
  pars["mu_o"] <- params[4]
  # 
  #creating model and replacing with new beta values 
  results <- as.data.frame(ode(y = initial_values, 
                               times = time_values, func = sir_equations, parms = pars))
  
  ### MODEL ### 
  #adding up all the types of B cells 
  dtotal <-  results %>% mutate(B_total = rowSums(across(c(Bb_cells, Cbb_cells)))) %>% mutate(T_total = rowSums(across(c(Tb_cells, Atb_cells, Ltb_cells, Ctb_cells)))) %>% 
    melt(id.vars = "time") %>% filter(variable %in% c("B_total", "T_total")) 
  
  #match model time with observed data 
  dtotallike <- dtotal %>% filter(time %in% obs_hours, variable == "B_total")    
  
  #empty list 
  #result <- vector(mode = "list", length = nrow(dtotallike))   
  LogLHoodStor=0
  for (i in seq_len(nrow(dtotallike))) { 
    obs_vec <- PBL_B %>% filter(time == dtotallike$time[i] ) %>% pull(value_B)
    mu <- dtotallike$value[i]  
    loglike <- dnorm(obs_vec, mu, sd = 1, log = TRUE) 
    LogLHoodStor=LogLHoodStor+sum(loglike)
    #result[[i]] <- data.frame(time = rep(dtotallike$time[i], length(obs_vec)) , obs = obs_vec, model_val = rep( mu, length(obs_vec)) , likelihood = loglike)
  }
  
  #sum_log <- result %>% map_dbl(~sum(.x[["likelihood"]])) 
  #likelihood <- -sum(sum_log) 
  # merged_data <- data %>%
  #   left_join(subset_I, by = "time")
  # 
  # MaxLike <- merged_data %>% mutate(logLike = log(dbinom(SampledInfected, SampleSize, I))) 
  # 
  # 
  # sumLog <- -sum(MaxLike$logLike)
  print(c(-LogLHoodStor, pars["beta"], pars["alpha"], pars["beta2"], pars["mu_o"]))
  return(-LogLHoodStor)
}


PBL_data <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/PayneRennie1976raw.xlsx")




#blah blah adding someing
parameters_values <- c(
  M = 0.005
  , beta =0.01      #contact rate with B cells
  , beta_2 =0.02          #contact rate with T cells
  , nu_A = 0.01             #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = 0.01            #Activation rate of T cells by cytolytic T cells (hours)
  , nu_F =0.07              #Infection rate of follicular cells (hours)
  , mu_o =  0.0001043031           #migration of blood cells out
  , mu_p = 0.0001            #rate of circulation to lymphoid organs
  , mu_t = 0.008
  , alpha =  0.01           #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.1          #death rate of cytolytic T cells (every 48 hours)
  , alpha_B = 0.001
  , theta = 0.8            #population of activated T cells
  , g1 =  0.01                #incoming B cells (every 15 hours)
  , g2 =0.01
  , h1 = 0                     #incoming T cells / determined no incoming T cells
  , h2 = 10
  , lambda = 0.005          #adding delay, how long latent cell 'exposed'
)

initial_values <- c( 
  B_cells = 10000
  , Cb = 0 
  , T_cells =50000
  , At = 0 
  , Lt = 0 
  , Lt2 = 0 
  , Lt3 = 0  
  , Lt4 = 0  
  , Lt5 = 0  
  #blood values 
  , Bb_cells = 0
  , Cbb_cells = 0  
  , Tb_cells = 0
  , Atb_cells = 0 
  , Ltb_cells = 0 
  , Ctb_cells = 0 
  , Ct = 0
  , f = 20 
  , If =0  
  
) 

time_values <- seq(0, 1080, by = 1) # hours
obs_hours <- as.numeric((unique(PBL_data$time)) * 24) 


#data from paper 
PBL_B <- PBL_data %>% filter(variable_B == "infectBCells", variable_T == "InfectTcells") %>% mutate(time = time*24)



test_parms <- c(10.819e-4, 0.01, 0.03, 0.02)
test_parms <- c(10.819e-2, 0.01, 0.03,0.02)

Likelihood(test_parms)


optim(par = test_parms, fn = Likelihood, method = "Nelder-Mead")
