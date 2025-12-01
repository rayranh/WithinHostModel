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
  with(as.list(c(variables, parameters)), { #turning initial values and parms into vectors and then list and then applying to below equations // include death of B cells + host response death apoptosis? 
    dB <- -M*B_cells - beta*Cb*B_cells - beta_2*Ct*B_cells + (g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct 
    dCb <- M*B_cells +beta*Cb*B_cells + beta_2*Ct*B_cells - alpha*Cb 
    dT <- -M*T_cells - nu_a*Cb*T_cells - nu_b*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct)))
    dAt <- M*T_cells + nu_a*Cb*T_cells + nu_b*Ct*T_cells - beta_2*Ct*At - beta*Cb*At  # beta = rate of activated T cells by Ct and Cb cells    
    dLt <- theta*(beta_2*Ct*At + beta*Cb*At)  - lambda*Lt 
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - lambda*Lt5 - mu*Lt5
    dCt <- (1-theta)*(beta_2*Ct*At + beta*Cb*At) - alpha_2*Ct
    dZ <- mu*Lt5
    df <- -nu_f*Lt5*f
    dIf <- nu_f*Lt5*f
    return(list(c(dB, dCb, dT, dAt,dLt,dLt2, dLt3,dLt4, dLt5, dCt, dZ, df,dIf)))
  }) 
}


Likelihood <- function(params){   
  pars <- parameters_values 
  pars["nu_f"] <- abs(params[1])  
  pars["beta"] <- abs(params[2]) 
  pars["beta_2"] <- abs(params[3]) 
  pars["nu_a"] <- abs(params[4]) 
  pars["nu_b"] <- abs(params[5]) 
  pars["mu"] <- abs(params[6])
  pars["alpha"] <- abs(params[7])
  pars["alpha_2"] <- abs(params[8]) 
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
  Infected_FFE<- results %>% dplyr::select("time", "If") %>% filter(time %in% matched_time) %>% arrange(time)  # changing to %in% because apparently better  
  Model_inf_cells <- results %>% dplyr::filter(time %in% obs_hourspp38) %>% mutate(Inf_cells = (Cb+Ct)/(B_cells + T_cells + At + Cb + Ct)) %>% arrange(time)
  #empty list 
  #result <- vector(mode = "list", length = nrow(dtotallike))   
  LogLHoodStor=0 
  
  result <- data.frame() 
  for (i in seq_len(length(matched_time))) { 
    obs_vec_B <- Infected_FFE$If[i]
    mu_B <- baigent2016$mean_genomes[i]
    loglike_B <- dnorm(log(obs_vec_B), log(mu_B), sd = 1, log = TRUE)  
    LogLHoodStor=LogLHoodStor+sum(loglike_B) 
   
    # result <- rbind(result, data.frame(
    #   i = i,
    #   time = matched_time[i],
    #   my_model = log(obs_vec_B),
    #   exp.data = log(mu_B),
    #   loglike = sum(loglike_B)))
  } 
  

  for (i in seq_len(nrow(Model_inf_cells))) {  
    obs_cyto_cells <- baigent1998_cyto$total_inf[i] 
    p <- Model_inf_cells$Inf_cells[i]
    loglike_B_pp38 <- dbinom(obs_cyto_cells, 40000, prob = p, log = TRUE) #binomial doesnt like when = 1?    
    LogLHoodStor=LogLHoodStor+sum(loglike_B_pp38,na.rm= TRUE)
    
  }  
 
  # result <- print(result)
# 
#   ggplot(result, aes(x = time, y = my_model)) + geom_line() + 
#     geom_line(aes(x = time, y = exp.data, colour = "red"))
  
 return(-LogLHoodStor)
}

#feather data 
baigent2016 <- read_xlsx("/Users/rayanhg/Downloads/baigent2016.xlsx", 2)  
baigent2016 %>% arrange(time)
#pp38 data 
baigent1998 <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", 4)  
baigent1998$mean.pp38 <- as.numeric(baigent1998$mean.pp38) #turn into factor 

baigent1998_cyto <- baigent1998 %>% arrange(time) %>% group_by(time) %>% summarize(total_inf = sum(mean.pp38, na.rm = TRUE)) # summarizing all results 





parameters_values  <- c( 
  M = 0
  , beta = 7.108697e-07                 #contact rate with B cells 
  , beta_2 = 2.239485e-07                   #contact rate with T cells 
  , nu_a = 0.05                     #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = 0.08                    #Activation rate of T cells by cytolytic T cells (hours)
  , nu_f = 0.06                      #Infection rate of follicular cells (hours)
  , mu = 0.5                        #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.55                   #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.55                  #death rate of cytolytic T cells (every 48 hours)
  , theta = 0.8                     #population of activated T cells 
  , g1 = 900                      #incoming B cells (every 15 hours)
  , g2 = 100    
  , h1 =900                        #incoming T cells / determined no incoming T cells 
  , h2 = 100  
  , lambda = 0.002                #adding delay, how long latent cell 'exposed' 
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
#the hours of baigent1998 pp38 data 
obs_hourspp38 <- c(24,48,72,96,120,144,192,720)

test_parms <- c(0.06,7.108e-7, 2.239485e-07,0.05,0.08,0.5,0.55,0.55,900,100,900,100,0.002)

Likelihood(test_parms)

#### remember to transform any exp() numbers back 
answer <- optim(par = test_parms, fn = Likelihood, method = "Nelder-Mead",control = list(maxit = 5000)) 


