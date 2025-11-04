
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(dplyr) 
library(patchwork) 
library(writexl) 
library(readxl)  
library(purrr) 
PBL_data <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/PayneRennie1976raw.xlsx")

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

#blah blah adding someing
parameters_values <- c(
  M = 5.0e-14
  , beta = 1.802641e-10    #contact rate with B cells
  , beta_2 = 6.139670e-10  #contact rate with T cells
  , nu_A = 1/10            #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b =  1/100          #Activation rate of T cells by cytolytic T cells (hours)
  , nu_F =0.07              #Infection rate of follicular cells (hours)
  , mu_o =  1/80            #rate of circulation of any B cell lineage into the blood  
  , mu_p = 1/100            #rate of circulation to lymphoid organs
  , mu_t =  1/60            #rate of circulation of any T cell lineage into the blood 
  , alpha = 1/5            #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 1/30           #death rate of cytolytic T cells (every 48 hours)
  , alpha_B = 0
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
  B_cells = 130432500/0.125,
  Cb = 0,
  Ct = 0,
  T_cells = 113602500/0.125,
  At = 0,
  Lt = 0,
  Lt2 = 0,
  Lt3 = 0,
  Lt4 = 0,
  Lt5 = 0,
  
  # Bursa
  B_bu = 954995000/0.325, 
  Cb_bu = 0,
  Ct_bu = 0, 
  T_bu = 7682750/0.325,
  At_bu = 0,
  Lt_bu = 0,
  Lt2_bu = 0,
  Lt3_bu = 0,
  Lt4_bu = 0,
  Lt5_bu = 0,
  
  # Thymus
  B_th = 3098250/0.425,
  Cb_th = 0, 
  Ct_th = 0, 
  T_th = 1048356000/0.425,
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
  f = 400000,
  If = 0,
  Z_sp = 0, 
  Z_th = 0, 
  Z_bu = 0       #### I feel like this is wrong 
  
) 

time_values <- seq(0, 1080, by = 1) # hours
obs_hours <- as.numeric((unique(PBL_data$time)) * 24) 



sir_values_1 <- ode(
  y = initial_values,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)  




#added strings as factors because otherwise turns into factors for some reason 
sir_values_1 <- as.data.frame(sir_values_1, stringsAsFactors = FALSE) 

#turning into long dataframe 
df <- melt(sir_values_1, id.vars = "time") %>% filter(variable %in% c("B_cells","Cb","Ct","T_cells","At","Lt5", "Z_sp"))  
df2 <- melt(sir_values_1, id.vars = "time") %>% filter(variable %in% c("B_bu","Cb_bu","Ct_bu", "T_bu","At_bu","Lt5_bu","Z_bu"))   
df3 <-  melt(sir_values_1, id.vars = "time") %>% filter(variable %in% c("B_th","Cb_th", "Ct_th", "T_th","At_th","Lt5_th","Z_th" )) 
dtotal <-  sir_values_1 %>% mutate(B_total = rowSums(across(c(Bb_cells, Cbb_cells))/7500)) %>% mutate(T_total = rowSums(across(c(Tb_cells, Atb_cells, Ltb_cells, Ctb_cells))/7500)) %>% 
  melt(id.vars = "time") %>% filter(variable %in% c("B_total", "T_total"))  
#data from paper 
PBL_B <- PBL_data %>% filter(variable_B == "infectBCells", variable_T == "InfectTcells") %>% mutate(time = time*24)
PBL_C <- PBL_data %>% filter(variable_B == "BCellsControl", variable_T == "TCellsControl") %>% mutate(time = time*24)
#matching total model times with obs data times 
dtotallike <- dtotal %>% filter(time %in% obs_hours, variable == "B_total") 


p1 <-ggplot(data = df, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() + 
  scale_color_manual(values = c("B_cells" = "black", "T_cells" = "red", "Cb"="green", "Lt5" = "purple", 
                                "Ct" = "yellow", "Z_sp" = "lightblue", "At"="blue")) +
  labs(title = "WithinHost Delay - Spleen", color = "Cell Type") + theme(legend.position = "right") + theme_minimal() + 
  xlab(label = "Time (Days)") + ylab(label = "Cell Number")

p2 <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
  scale_color_manual(values = c("B_bu" = "black", "T_bu" = "red", "Cb_bu"="green", "Lt5_bu" = "purple",
                                "Ct_bu" = "yellow", "Z_bu" = "lightblue", "At_bu"="blue")) +
  labs(title = "WithinHost Delay - Bursa", color = "Cell Type") + theme(legend.position = "right") + theme_minimal() +
  xlab(label = "Time (Days)") + ylab(label = "Cell Number")

p3 <- ggplot(data = df3, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() +
  scale_color_manual(values = c("B_th" = "black", "T_th" = "red", "Cb_th"="green", "Lt_th5" = "purple",
                                "Ct_th" = "yellow", "Z_th" = "lightblue", "At_th"="blue")) +
  labs(title = "WithinHost Delay - Thymus", color = "Cell Type") + theme(legend.position = "right") + theme_minimal() +
  xlab(label = "Time (Days)") + ylab(label = "Cell Number")

p4 <- ggplot(data = PBL_B, mapping = aes(x = time/24, y = value_B, colour = "infectBCells")) + geom_point() +
  scale_color_manual(values = c("infectBCells" = "blue", "InfectTcells"="red", "B_total" = "blue", "T_total" = "red"))  + 
  geom_point(data =PBL_B, mapping = aes(x = time/24, y = value_T, colour = "InfectTcells")) +
  geom_line(data = dtotal %>% filter(variable == "B_total"), aes(x = time/24, y = value, colour = "B_total")) +
  geom_line(data = dtotal %>% filter(variable == "T_total"), aes(x = time/24, y = value, colour = "T_total")) + scale_y_log10(limits = c(1000, 100000)) + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),  legend.text = element_text(size = 12),
        legend.title = element_text(size = 12), axis.line = element_line(color = "black")) + labs( y = "Cell Number", x = "Time (Days)")



(p1 | p2| p3)/p4 + plot_layout(widths = c(2,2)) 

p4


