rm(list = ls())
library(dplyr) 
library(tibble)
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(writexl) 
library(readxl)  
library(purrr)  
library(tictoc) 
library(GGally)

set.seed(45)
#create intervals for numbers 
parameter_intervals <- list(beta = c(4e-8,1e-5), 
                      beta_2 = c(4e-8,1e-5), 
                      alpha_2 = c(0.0104, 0.041), 
                      g1 = c(10,1000), 
                      g2 = c(1,1000), 
                      h1 = c(10,1000), 
                      h2 = c(1,1000))  

alpha <- seq(0.01041667, 0.08333333, by = 0.001 )

#generate random parameters 

#it can only do one number a time
rand_number_generation <- function(rng, test_param){ 
  minimum <- rng[1]
  maximum <- rng[2]
  rand_num <- runif(1,min = minimum, max = maximum) 
} 

# map_dbl makes it go one element at a time [[subsets]], generating random number vector (1:8)
making_a_set <- function() { 
  random_parameters<- map_dbl(parameter_intervals, ~rand_number_generation(.)) 
  as_tibble_row(random_parameters) 
} 

#how many random parameter sets I want per alpha 
n_per_alpha <- 5 
#repeating this 5 times for each value of alpha so that its equally run  
df <- map_df(1:(length(alpha)*n_per_alpha), ~making_a_set()) 

#adding my alpha values to df 
df_final <- df %>% mutate(alpha = rep_len(alpha, length = nrow(df))) %>% arrange(alpha) 

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

## PARAMETERS AND INITIAL VALUES ## 
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

tic()
#rowwise groups by rows and c_across extracts the variables we want 
df_with_Likelihood <- df_final %>% rowwise() %>% 
  mutate(Likely = Likelihood(c(beta = beta, beta_2 = beta_2, 
                               alpha = alpha, alpha_2 = alpha_2, 
                               g1 = g1, g2 = g2, h1 = h1, h2 = h2))) %>% 
  ungroup() %>% arrange(Likely)  

toc()

ggplot(df_with_Likelihood, aes(x = alpha, y = Likely)) + geom_point()

plot(density(df_with_Likelihood$Likely))

write.csv(df_with_Likelihood, file = "/Users/rayanhg/Desktop/WithinHostModel/CodeOutputs/alpha_parameter_exploration")
