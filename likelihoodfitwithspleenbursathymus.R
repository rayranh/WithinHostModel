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
    ### Spleen ### 
    dB <- -M*B_cells - beta*Cb*B_cells - beta_2*Ct*B_cells + (g1*(Cb+Ct)/(g2+(Cb+Ct))) - mu_o*B_cells + mu_p*Bb_cells # beta = rate of cytolytically infected B cells by Cb and Ct 
    dCb <- M*B_cells + beta*Cb*B_cells + beta_2*Ct*B_cells - mu_o*Cb + mu_p*Cbb_cells - alpha*Cb
    dCt <- (1-theta)*(beta_2*Ct*At + beta*Cb*At) - mu_t*Ct + mu_p*Ctb_cells - alpha_2*Ct # mu_t made T cells circulate at a different rate than B cells 
    dT <- -M*T_cells - nu_A*Cb*T_cells - nu_b*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct))) - mu_t*T_cells + mu_p*Tb_cells
    dAt <- M*T_cells + nu_A*Cb*T_cells + nu_b*Ct*T_cells - beta_2*Ct*At - beta*Cb*At - mu_t*At + mu_p*Atb_cells  # beta = rate of activated T cells by Ct and Cb cells 
    dLt <- theta*(beta_2*Ct*At + beta*Cb*At) - lambda*Lt  #remove mu*Lt since only the last compartment should contribute to tumors  
    dLt2 <- lambda*Lt - lambda*Lt2 
    dLt3 <- lambda*Lt2 - lambda*Lt3 
    dLt4 <- lambda*Lt3 - lambda*Lt4
    dLt5 <- lambda*Lt4 - lambda*Lt5 - mu_t*Lt5 - epsilon*Lt5 +  mu_p*Ltb_cells   
    dZ_sp <- epsilon*Lt5 
    
    ### Bursa ### 
    dB_bu   <- -M*B_bu - beta*Cb_bu*B_bu - beta_2*Ct_bu*B_bu + (g1*(Cb_bu+Ct_bu)/(g2+(Cb_bu+Ct_bu))) - mu_o*B_bu + mu_p*Bb_cells
    dCb_bu  <-  M*B_bu + beta*Cb_bu*B_bu + beta_2*Ct_bu*B_bu - alpha*Cb_bu - mu_o*Cb_bu + mu_p*Cbb_cells
    dCt_bu  <- (1-theta)*(beta_2*Ct_bu*At_bu + beta*Cb_bu*At_bu) - alpha_2*Ct_bu - mu_t*Ct_bu + mu_p*Ctb_cells
    dT_bu   <- -M*T_bu - nu_A*Cb_bu*T_bu - nu_b*Ct_bu*T_bu + (h1*(Cb_bu+Ct_bu)/(h2+(Cb_bu+Ct_bu))) - mu_t*T_bu + mu_p*Tb_cells
    dAt_bu  <-  M*T_bu + nu_A*Cb_bu*T_bu + nu_b*Ct_bu*T_bu - beta_2*Ct_bu*At_bu - beta*Cb_bu*At_bu - mu_t*At_bu + mu_p*Atb_cells  
    dLt_bu  <-  theta*(beta_2*Ct_bu*At_bu + beta*Cb_bu*At_bu) - lambda*Lt_bu + mu_p*Ltb_cells   # no mu*Lt_bu, only terminal contributes to tumors
    dLt2_bu <-  lambda*Lt_bu - lambda*Lt2_bu 
    dLt3_bu <-  lambda*Lt2_bu - lambda*Lt3_bu 
    dLt4_bu <-  lambda*Lt3_bu - lambda*Lt4_bu
    dLt5_bu <-  lambda*Lt4_bu - lambda*Lt5_bu - mu_t*Lt5_bu - epsilon*Lt5_bu  
    dZ_bu <-  epsilon*Lt5_bu
    
    
    ### Thymus ###
    dB_th   <- -M*B_th - beta*Cb_th*B_th - beta_2*Ct_th*B_th + (g1*(Cb_th+Ct_th)/(g2+(Cb_th+Ct_th))) - mu_o*B_th + mu_p*Bb_cells
    dCb_th  <-  M*B_th + beta*Cb_th*B_th + beta_2*Ct_th*B_th - alpha*Cb_th - mu_o*Cb_th + mu_p*Cbb_cells
    dCt_th  <- (1-theta)*(beta_2*Ct_th*At_th + beta*Cb_th*At_th) - alpha_2*Ct_th - mu_t*Ct_th + mu_p*Ctb_cells
    dT_th   <- -M*T_th - nu_A*Cb_th*T_th - nu_b*Ct_th*T_th + (h1*(Cb_th+Ct_th)/(h2+(Cb_th+Ct_th))) - mu_t*T_th + mu_p*Tb_cells
    dAt_th  <-  M*T_th + nu_A*Cb_th*T_th + nu_b*Ct_th*T_th - beta_2*Ct_th*At_th - beta*Cb_th*At_th - mu_t*At_th + mu_p*Atb_cells  
    dLt_th  <-  theta*(beta_2*Ct_th*At_th + beta*Cb_th*At_th) - lambda*Lt_th + mu_p*Ltb_cells     # no mu*Lt_th, only terminal contributes to tumors
    dLt2_th <-  lambda*Lt_th - lambda*Lt2_th 
    dLt3_th <-  lambda*Lt2_th - lambda*Lt3_th 
    dLt4_th <-  lambda*Lt3_th - lambda*Lt4_th
    dLt5_th <-  lambda*Lt4_th - lambda*Lt5_th - mu_t*Lt5_th - epsilon*Lt5_th   
    dZ_th <-  epsilon*Lt5_th
    
    
    
    ### entering blood ###   
    dBb_cells <- mu_o*B_cells + mu_o*B_th + mu_o*B_bu - 3*mu_p*Bb_cells  - alpha_B*Bb_cells  
    dCbb_cells <- mu_o*Cb + mu_o*Cb_th + mu_o*Cb_bu - 3*mu_p*Cbb_cells  - alpha_B*Cbb_cells   
    dTb_cells  <- mu_t*T_cells + mu_t*T_th + mu_t*T_bu - 3*mu_p*Tb_cells - alpha_B*Tb_cells  
    dAtb_cells <- mu_t*At + mu_t*At_th + mu_t*At_bu - 3*mu_p*Atb_cells   # no death for At cells 
    dLtb_cells <- mu_t*Lt5 + mu_t*Lt5_th + mu_t*Lt5_bu - 3*mu_p*Ltb_cells - epsilon*f*Ltb_cells  
    dCtb_cells  <- mu_t*Ct + mu_t*Ct_th + mu_t*Ct_bu - 3*mu_p*Ctb_cells - alpha_B*Ctb_cells   
    
    
    ### feathers ### 
    df <- -nu_F*Ltb_cells*f 
    dIf <- nu_F*Ltb_cells*f 
    return(list(c(
      # Spleen
      dB, dCb, dCt, dT, dAt, dLt, dLt2, dLt3, dLt4, dLt5,
      
      # Bursa
      dB_bu, dCb_bu, dCt_bu, dT_bu, dAt_bu, dLt_bu, dLt2_bu, dLt3_bu, dLt4_bu, dLt5_bu,
      
      # Thymus
      dB_th, dCb_th, dCt_th, dT_th, dAt_th, dLt_th, dLt2_th, dLt3_th, dLt4_th, dLt5_th,
      
      # Blood (global + organ-specific)
      dBb_cells,dCbb_cells, dTb_cells, dAtb_cells, dLtb_cells, dCtb_cells, 
      
      # Feather & tumor
      df, dIf, dZ_sp, dZ_th, dZ_bu
    )))
    
  }) 
}

params <- c(1.347721e-15, 3.693628e-7)
Likelihood <- function(params){   
  pars <- parameters_values 
  pars["beta"] <- abs(params[1])  
  pars["beta_2"] <-  abs(params[2])
  
  
  
  # 
  #creating model and replacing with new beta values 
  results <- as.data.frame(ode(y = initial_values, 
                               times = time_values, func = sir_equations, parms = pars))
    
  long_results <- melt(results, id.vars = "time") 
  
  long_df_bu <- long_results %>% filter(variable %in% c("B_bu", "T_bu", "Ct_bu", "Cb_bu", "Lt_bu", "Z_bu", "At_bu")) %>% mutate(value = log(value + 1))
  long_df_sp <- long_results %>% filter(variable %in% c("B_cells", "T_cells", "Ct", "Cb", "Lt", "Z", "At")) %>% mutate(value = log(value+1) )
  long_df_th<- long_results %>% filter(variable %in% c("B_th", "T_th", "Ct_th", "Cb_th", "Lt_th", "Z_th", "At_th")) 
  
  model_output <- ggplot(data = long_df_th, mapping = aes(x = time, y = value, colour = variable)) + geom_line() +
    labs(title = "Model output", x = "time", y = "cells", colour = "variable" ) 
  # + scale_color_manual(values = c("At_bu" = "red") )

  
  ### MODEL ### 
  #adding up all the types of B cells and all types of T cells
  dtotal <-  results %>% mutate(B_total = rowSums(across(c(Bb_cells, Cbb_cells))/7400)) %>% mutate(T_total = rowSums(across(c(Tb_cells, Atb_cells, Ltb_cells, Ctb_cells))/7400)) %>% 
    melt(id.vars = "time") %>% filter(variable %in% c("B_total", "T_total")) 
  
  #match model time with observed PBL data 
  dtotallike_B <- dtotal %>% filter(time %in% obs_hours, variable %in% c("B_total"))   # changing to %in% because apparently better   
  dtotallike_T <- dtotal %>% filter(time %in% obs_hours, variable %in% c("T_total"))   
  
  
  #getting % CYTOLYTICALLY infected from model - SPLEEN (filtering by spleen via initial value name )
  infprob <- results %>% select(B_cells, Cb, Ct, T_cells, At,Lt, Lt2, Lt3, Lt4, Lt5, time)  %>% mutate(prob_Cb = Cb/(Cb+B_cells), prob_Ct = Ct/(T_cells+Ct+Lt+Lt2+Lt3+Lt4+Lt5+At)) 
  infprob_Bu <- results %>% select(T_bu, B_bu, Cb_bu, Ct_bu, T_bu, At_bu,Lt_bu, Lt2_bu, Lt3_bu, Lt4_bu, Lt5_bu, time)  %>% mutate(prob_Cb_bu = Cb_bu/(Cb_bu+B_bu), prob_Ct_bu = Ct_bu/(T_bu+Ct_bu+Lt_bu+Lt2_bu+Lt3_bu+Lt4_bu+Lt5_bu+At_bu)) 
  infprob_Th <- results %>% select(T_th, B_th, Cb_th, Ct_th, T_th, At_th,Lt_th, Lt2_th, Lt3_th, Lt4_th, Lt5_th, time)  %>% mutate(prob_Cb_th = Cb_th/(Cb_th+B_th), prob_Ct_th = Ct_th/(T_th+Ct_th+Lt_th+Lt2_th+Lt3_th+Lt4_th+Lt5_th+At_th)) 
  
  #match MODEL time with observed PP38 data time (in hours) 
  matchedTime <- infprob %>% filter(time %in% obs_hourspp38) %>% mutate(Bprop = B_cells/(B_cells+T_cells), Tprop = T_cells/(B_cells+T_cells)) 
  matchedTime_Bu <- infprob_Bu %>% filter(time %in% obs_hourspp38) %>% mutate(Bprop = B_bu/(B_bu+T_bu), Tprop = T_bu/(B_bu+T_bu))  
  matchedTime_Th <- infprob_Th %>% filter(time %in% obs_hourspp38) %>% mutate(Bprop = B_th/(B_th+T_th), Tprop = T_th/(B_th+T_th))  
  
  
  #empty list 
  #result <- vector(mode = "list", length = nrow(dtotallike))   
  LogLHoodStor=0 
  
  for (i in seq_len(nrow(dtotallike_B)) ) {  
    obs_vec_B <- PBL_B %>% filter(time %in% dtotallike_B$time[i]) %>% pull(value_B) 
    obs_vec_T <- PBL_B %>% filter(time %in% dtotallike_T$time[i]) %>% pull(value_T)
    mu_B <- dtotallike_B$value[i]
    if(mu_B<0){print(mu_B)}
    mu_T <- dtotallike_T$value[i]
    loglike_B <- dnorm(log(obs_vec_B), log(mu_B), sd = 1, log = TRUE)  
    loglike_T <- dnorm(log(obs_vec_T), log(mu_T), sd = 1, log = TRUE)  
    LogLHoodStor=LogLHoodStor+sum(loglike_B, loglike_T)
    #result[[i]] <- data.frame(time = rep(dtotallike$time[i], length(obs_vec)) , obs = obs_vec, model_val = rep( mu, length(obs_vec)) , likelihood = loglike) 
  }
  
  for (i in seq_len(nrow(matchedTime))) { 
    obs_cytoBcells <- pp38_Spleen %>% filter(time %in% matchedTime$time[i]) %>% pull(Bcell_no) # matching data to model time for B  
    obs_cytoTcells <- pp38_Spleen %>% filter(time %in% matchedTime$time[i]) %>% pull(Tcell_no)
    loglike_B_pp38 <- dbinom(obs_cytoBcells, as.integer(round(40000*matchedTime$Bprop[i])), prob = matchedTime$prob_Cb[i], log = TRUE) #binomial doesnt like when = 1?    
    loglike_T_pp38 <- dbinom(obs_cytoTcells, as.integer(round(40000*matchedTime$Tprop[i])), prob = matchedTime$prob_Ct[i], log = TRUE)
    LogLHoodStor=LogLHoodStor+sum(loglike_B_pp38,loglike_T_pp38, na.rm= TRUE)
    #result[[i]] <- data.frame(time = rep(dtotallike$time[i], length(obs_vec)) , obs = obs_vec, model_val = rep( mu, length(obs_vec)) , likelihood = loglike) 
    
  }  
  
  for (i in seq_len(nrow(matchedTime_Bu))) {
    obs_cytoBcell_Bu <- pp38_Bursa %>% filter(time %in% matchedTime_Bu$time[i]) %>% pull(Bcell_no) # matching data to model time for B
    obs_cytoTcell_Bu <- pp38_Bursa %>% filter(time %in% matchedTime_Bu$time[i]) %>% pull(Tcell_no)
    loglike_B_pp38_Bu <- dbinom(obs_cytoBcell_Bu, as.integer(round(40000*matchedTime_Bu$Bprop[i])), prob = matchedTime_Bu$prob_Cb_bu[i], log = TRUE) #binomial doesnt like when = 1?
    loglike_T_pp38_Bu <- dbinom(obs_cytoTcell_Bu, as.integer(round(40000*matchedTime_Bu$Tprop[i])), prob = matchedTime_Bu$prob_Ct_bu[i], log = TRUE)
    LogLHoodStor=LogLHoodStor+sum(loglike_B_pp38_Bu,loglike_T_pp38_Bu, na.rm= TRUE)
    #result[[i]] <- data.frame(time = rep(dtotallike$time[i], length(obs_vec)) , obs = obs_vec, model_val = rep( mu, length(obs_vec)) , likelihood = loglike)

  }
  for (i in seq_len(nrow(matchedTime_Th))) {
    obs_cytoBcell_Th <- pp38_Thymus %>% filter(time %in% matchedTime_Th$time[i]) %>% pull(Bcell_no) # matching data to model time for B
    obs_cytoTcell_Th <- pp38_Thymus %>% filter(time %in% matchedTime_Th$time[i]) %>% pull(Tcell_no)
    loglike_B_pp38_Th <- dbinom(obs_cytoBcell_Th, as.integer(round(40000*matchedTime_Th$Bprop[i])), prob = matchedTime_Th$prob_Cb_th[i], log = TRUE) #binomial doesnt like when = 1?
    loglike_T_pp38_Th <- dbinom(obs_cytoTcell_Th, as.integer(round(40000*matchedTime_Th$Tprop[i])), prob = matchedTime_Th$prob_Ct_th[i], log = TRUE)
    LogLHoodStor=LogLHoodStor+sum(loglike_B_pp38_Th,loglike_T_pp38_Th, na.rm= TRUE)
    #result[[i]] <- data.frame(time = rep(dtotallike$time[i], length(obs_vec)) , obs = obs_vec, model_val = rep( mu, length(obs_vec)) , likelihood = loglike)

  } 
  
  print(model_output)
  # 
  #sum_log <- result %>% map_dbl(~sum(.x[["likelihood"]])) 
  #likelihood <- -sum(sum_log) 
  # merged_data <- data %>%
  #   left_join(subset_I, by = "time")
  # 
  # MaxLike <- merged_data %>% mutate(logLike = log(dbinom(SampledInfected, SampleSize, I))) 
  # 
  # 
  # sumLog <- -sum(MaxLike$logLike)
  #print(c(-LogLHoodStor, pars["mu_o"], pars["beta"], pars["beta_2"], pars["mu_t"]))
  return(-LogLHoodStor) 
  
}

## DATA ##  
# PBL data for B and T cell infection in blood 
PBL_data <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/PayneRennie1976raw.xlsx")  

#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 
pp38_dat <- read_xlsx("baigent1998.xlsx", sheet = 4, na = "NA") 

#leaving NA for correct loop number, and filtering for Spleen organ only, Combining all T cells (CD4+CD8)
pp38_Spleen <- pp38_dat %>% filter(Organ == "spleen") %>% mutate(Tcell_no = CD4_no+CD8_no)
pp38_Bursa <- pp38_dat %>% filter(Organ == "bursa") %>% mutate(Tcell_no = CD4_no+CD8_no)
pp38_Thymus <- pp38_dat %>% filter(Organ == "thymus") %>% mutate(Tcell_no = CD4_no+CD8_no) 

#tried adding data from baigent paper..should I do mean or not? Also, thinking about the scaling of the model, but we did dbinom to the proportions so it should be fine 
long_pp38_Spleen <- pp38_Spleen %>% select(time, Organ, Bcell_no, Tcell_no) %>% melt(id.vars = c("time", "Organ")) %>% 
  group_by(time, variable, Organ) %>% summarise(mean = mean(value, na.rm = TRUE))

long_pp38_Bursa <- pp38_Bursa %>% select(time, Organ, Bcell_no, Tcell_no) %>% melt(id.vars = c("time", "Organ")) 
long_pp38_Thymus <- pp38_Thymus %>% select(time, Organ, Bcell_no, Tcell_no) %>% melt(id.vars = c("time", "Organ")) 

ggplot( data = long_pp38_Spleen, mapping = aes(x = time, y = mean, color = variable)) + geom_line()  

#blah blah adding someing
parameters_values <- c(
  M = 5.0e-14
  , beta = 1.347721e-100    #contact rate with B cells
  , beta_2 = 3.693628e-400  #contact rate with T cells
  , nu_A = 1/10            #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b =  1/48          #Activation rate of T cells by cytolytic T cells (hours)
  , nu_F =0.07              #Infection rate of follicular cells (hours)
  , mu_o =  1/72            #rate of circulation of any B cell lineage into the blood  
  , mu_p = 1/100             #rate of circulation to lymphoid organs
  , mu_t =  1/50            #rate of circulation of any T cell lineage into the blood 
  , alpha = 1/33            #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 1/48           #death rate of cytolytic T cells (every 48 hours)
  , alpha_B = 0.1
  , theta = 0.8            #population of activated T cells
  , g1 =10                #incoming B cells (every 15 hours)
  , g2 =1 
  , h1 = 10                #incoming T cells / determined no incoming T cells
  , h2 = 1
  , lambda = 0.005 #adding delay, how long latent cell 'exposed'  
  , epsilon = 0.01
)

initial_values <- c( 
  # Spleen
  B_cells = 521730000 ,
  Cb = 0,
  Ct = 0,
  T_cells = 454410000 ,
  At = 0,
  Lt = 0,
  Lt2 = 0,
  Lt3 = 0,
  Lt4 = 0,
  Lt5 = 0,
  
  # Bursa
  B_bu = 3819980000, 
  Cb_bu = 0,
  Ct_bu = 0, 
  T_bu = 30731000,
  At_bu = 0,
  Lt_bu = 0,
  Lt2_bu = 0,
  Lt3_bu = 0,
  Lt4_bu = 0,
  Lt5_bu = 0,
  
  # Thymus
  B_th = 12393000,
  Cb_th = 0, 
  Ct_th = 0, 
  T_th = 4193424000,
  At_th = 0,
  Lt_th = 0,
  Lt2_th = 0,
  Lt3_th = 0,
  Lt4_th = 0,
  Lt5_th = 0,
  
  # Blood pools
  Bb_cells = 0,
  Cbb_cells = 0,
  Tb_cells = 0,
  Atb_cells = 0,
  Ltb_cells = 0,
  Ctb_cells = 0,
  
  # Feather follicle & tumor
  f = 20,
  If = 0,
  dZ_sp = 0, 
  dZ_th = 0, 
  dZ_bu = 0 
) 

time_values <- seq(0, 1080, by = 1) # hours875875462 

# PBL TIME 
obs_hours <- as.numeric((unique(PBL_data$time)) * 24)   
#PP38 Time 
obs_hourspp38 <- c(72,96,120,144)


#data from paper 
PBL_B <- PBL_data %>% filter(variable_B == "infectBCells", variable_T == "InfectTcells") %>% mutate(time = time*24)

test_parms <- c(1.347721e-1, 3.693628e-1)


Likelihood(test_parms)

#### remember to transform any exp() numbers back, returned -loglikelihood so smaller number is larger 
answeroptim <- optim(par = test_parms, fn = Likelihood, method = "Nelder-Mead",control = list(maxit = 5000)) 

parameters <- print(answeroptim$par)
paste("beta=", parameters[1], "beta_2=", parameters[2])

