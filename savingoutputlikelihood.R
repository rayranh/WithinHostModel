library(deSolve)  
library(reshape2) 
library(ggplot2)
library(dplyr) 
library(patchwork) 
library(writexl) 
library(readxl)  
library(purrr)  
library(tibble)

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

Likelihood_info <- function(params){   
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
  # % cytolytically infected in model
  infprob <- results %>% 
    dplyr::select(B_cells, Cb, Ct, T_cells, At, Lt, Lt2, Lt3, Lt4, Lt5, time) %>% 
    mutate(
      prob_Cyto = (Cb + Ct) / (B_cells + Cb + Ct + T_cells + At + Lt + Lt2 + Lt3 + Lt4 + Lt5)
    ) 
  
  
  Infected_FFE<- results %>% dplyr::select("time", "If") %>% 
    filter(time %in% matched_time) %>% arrange(time)
  
  
  # match model time to observed pp38 (hours)
  matchedTime <- infprob %>% 
    filter(time %in% obs_hourspp38) %>% 
    mutate(
      Bprop = B_cells / (B_cells + T_cells),
      Tprop = T_cells / (B_cells + T_cells)
    )
  
  
  LogLHoodStor <- 0   
  LogLhoodbaigent2016 <- 0 
  LogLhoodbaigent1998 <- 0 
  
  # result <- data.frame() 
  for (i in seq_len(length(matched_time))) { 
    obs_vec_B <- Infected_FFE$If[i]
    mu_B <- baigent2016$mean_genomes[i]
    loglike_B <- dnorm(log(obs_vec_B), log(mu_B), sd = 1, log = TRUE)  
    LogLHoodStor=LogLHoodStor+sum(loglike_B) 
    LogLhoodbaigent2016 = LogLhoodbaigent2016 + sum(loglike_B) 
    
  } 
  
  
  
  for (i in seq_len(nrow(matchedTime))) {  
    
    p <- matchedTime$prob_Cyto[i]

    
    # ----- DIAGNOSTIC CHECK -----
    if (p < 0) { 
      p <- 0 
    }
      # DO NOT skip, DO NOT return — continue loop 

    # -----------------------------
    # observed B pp38+ cells at this time
    obs_cytoBcells <- pp38_Spleen %>% filter(time == matchedTime$time[i]) %>% pull(total)
    
    # binomial likelihood
    loglike_B_pp38 <- dbinom(x= obs_cytoBcells, size = 40000, prob = p,log  = TRUE) 
    
    LogLHoodStor <- LogLHoodStor + sum(loglike_B_pp38, na.rm = TRUE) 
    LogLhoodbaigent1998 = LogLhoodbaigent1998 + sum(loglike_B_pp38)  
    
  }  
  
  
  return(list(initial_parms = params, Likelihood_pp38 = -LogLhoodbaigent1998, 
        Likelihood_feathers = -LogLhoodbaigent2016, total_likelhood = -LogLHoodStor)) 
  
}



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
  # % cytolytically infected in model
  infprob <- results %>% 
    dplyr::select(B_cells, Cb, Ct, T_cells, At, Lt, Lt2, Lt3, Lt4, Lt5, time) %>% 
    mutate(
      prob_Cyto = (Cb + Ct) / (B_cells + Cb + Ct + T_cells + At + Lt + Lt2 + Lt3 + Lt4 + Lt5)
    ) 
  
  
  Infected_FFE<- results %>% dplyr::select("time", "If") %>% 
    filter(time %in% matched_time) %>% arrange(time)
  
  
  # match model time to observed pp38 (hours)
  matchedTime <- infprob %>% 
    filter(time %in% obs_hourspp38) %>% 
    mutate(
      Bprop = B_cells / (B_cells + T_cells),
      Tprop = T_cells / (B_cells + T_cells)
    )
  
  LogLHoodStor <- 0 
  # result <- data.frame() 
  for (i in seq_len(length(matched_time))) { 
    obs_vec_B <- Infected_FFE$If[i]
    mu_B <- baigent2016$mean_genomes[i]
    loglike_B <- dnorm(log(obs_vec_B), log(mu_B), sd = 1, log = TRUE)  
    LogLHoodStor=LogLHoodStor+sum(loglike_B) 
    
  } 
  
  
  
  for (i in seq_len(nrow(matchedTime))) {  
    
    p <- matchedTime$prob_Cyto[i]
    
    # ----- DIAGNOSTIC CHECK -----
    if (p < 0) { 
      p <- 0 
    }
    # -----------------------------
    # observed B pp38+ cells at this time
    obs_cytoBcells <- pp38_Spleen %>% filter(time == matchedTime$time[i]) %>% pull(total)
    
    # binomial likelihood
    loglike_B_pp38 <- dbinom(x= obs_cytoBcells, size = 40000, prob = p,log  = TRUE) 
    
    LogLHoodStor <- LogLHoodStor + sum(loglike_B_pp38, na.rm = TRUE) 
  } 
      
return(-LogLHoodStor)  
}

## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 
pp38_dat <- read_xlsx("baigent1998.xlsx", sheet = 3, na = "NA") %>% filter(!is.na(mean.pp38)) 
baigent2016 <- read_xlsx("/Users/rayanhg/Downloads/baigent2016.xlsx", 2 ) %>% arrange(time)


pp38_Spleen <- pp38_dat %>% mutate(Tcell_no = CD4_no+CD8_no) %>% arrange(time) %>% 
  group_by(time) %>% summarise(total = sum(mean.pp38))

#blah blah adding someing
#blah blah adding someing
parameters_values <- c( 
  M = 0
  , beta =  6.951463e-07                  #contact rate with B cells 
  , beta_2 = 8.532282e-07                  #contact rate with T cells 
  , nu_a = 4.668718e-01                     #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = 4.467098e-01                   #Activation rate of T cells by cytolytic T cells (hours)
  , nu_f = 6.033255e-02                      #Infection rate of follicular cells (hours)
  , mu = 7.165224e-01                        #Rate of Tumor Cells (every 72 hours)
  , alpha = 5.500001e-01                   #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 5.656785e-01                  #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 1000                     #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 =1000                      #incoming T cells / determined no incoming T cells 
  , h2 = 100
  , lambda = 15.975984e-05                 #adding delay, how long latent cell 'exposed' 
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


time_values <- seq(0, 1080, by = 1)   
#matching to the wack hours of baigent2016
matched_time <- c(3,6,10,17,20,26,33,38)*24 


#PP38 Time 
obs_hourspp38 <- c(72,96,120,144)

test_parms <- c( beta = 6.951463E-08                 
                 , beta_2 = 8.532282E-08
                 , alpha = 0.2300001   
                 , alpha_2 = 0.6656785  
                 , g1 = 1000
                 , g2 = 100 
                 , h1 = 1000 
                 , h2 = 100) 


### Calling Functions ### 
Likelihood(test_parms)  
report <- Likelihood_info(test_parms)  

#### remember to transform any exp() numbers back, returned -loglikelihood so smaller number is larger 
answeroptim <- optim(par = test_parms, fn = Likelihood, method = "Nelder-Mead",control = list(maxit = 5000))   


#storing answeroptim 
df_optimized <- as_tibble_row(as.list(answeroptim$par))  %>% mutate(Likelihood_total = answeroptim$value,              
                                                                    converged = answeroptim$convergence, 
                                                                    status = "Optimized Values")  

All_initial_params <- modifyList(as.list(parameters_values), as.list(report$initial_parms))
df_initial <-as_tibble_row(as.list(All_initial_params)) %>% mutate(Likelihood_pp38 = report$Likelihood_pp38, 
                                                                      Likelihood_feathers = report$Likelihood_feathers, 
                                                                      Likelihood_total = report$total_likelhood, 
                                                                      status = "Initial Values") 
df_final <- dplyr::bind_rows(df_initial, df_optimized)  

existing_files <- list.files("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputs/", pattern = "LikelihoodOutputs.csv")
# 
# write.csv(df_final, file =paste0("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputs/", Sys.Date(),"run", runid, "LikelihoodOutputs.csv" ), row.names = FALSE)
# 

