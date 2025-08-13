rm(list = ls())
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(dplyr) 
library(patchwork) 
library(writexl) 
library(readxl)
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

#blah blah adding someing 
parameters_values <- c( 
  M = 0.005
  , beta = 10.819e-4       #contact rate with B cells 
  , beta_2 = 5e-4           #contact rate with T cells 
  , nu_A = 0.05             #Activation rate of T cells by cytolytic B cells (hours)
  , nu_b = 0.001            #Activation rate of T cells by cytolytic T cells (hours)
  , nu_F =0.07              #Infection rate of follicular cells (hours)
  , mu_o = 0.005            #migration of blood cells out 
  , mu_p = 0.001            #rate of circulation to lymphoid organs 
  , mu_t = 0.008
  , alpha =  0.01           #death rate of cytolytic B cells (every 33 hours)
  , alpha_2 = 0.01          #death rate of cytolytic T cells (every 48 hours) 
  , alpha_B = 0.001
  , theta = 0.8            #population of activated T cells 
  , g1 =  0                #incoming B cells (every 15 hours)
  , g2 =0.001    
  , h1 = 0                     #incoming T cells / determined no incoming T cells 
  , h2 = 10  
  , lambda = 0.007          #adding delay, how long latent cell 'exposed' 
)

initial_values <- c( 
  B_cells = 56
  , Cb = 0 
  , T_cells =87
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

time_values <- seq(0, 1000) # hours

sir_values_1 <- ode(
  y = initial_values,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)  




#added strings as factors because otherwise turns into factors for some reason 
sir_values_1 <- as.data.frame(sir_values_1, stringsAsFactors = FALSE) 
#turning it longways 
df <- melt(sir_values_1, id.vars = "time") %>% filter(variable %in% c("B_cells", "T_cells", "Cb", "Lt5", "f", "If")) 
df2 <- melt(sir_values_1, id.vars = "time") %>% filter(variable %in% c("Bb_cells", "Cbb_cells", "Tb_cells", "Atb_cells", "Ltb_cells", "Ctb_cells"))  
dtotal <-  sir_values_1 %>% mutate(B_total = rowSums(across(c(Bb_cells, Cbb_cells)))) %>% mutate(T_total = rowSums(across(c(Tb_cells, Atb_cells, Ltb_cells, Ctb_cells)))) %>% 
  melt(id.vars = "time") %>% filter(variable %in% c("B_total", "T_total"))   


p1 <-ggplot(data = df, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() + 
  scale_color_manual(values = c("B_cells" = "black", "T_cells" = "red", "Cb"="green", "At" = "blue", "Lt5" = "purple", 
                                "Ct" = "yellow", "Z" = "lightblue", "f" = "aquamarine", "If" = "pink")) +
  labs(title = "WithinHost Delay - Spleen", color = "Cell Type") + theme(legend.position = "right") + theme_minimal() + xlab(label = "Time (Days)") + ylab(label = "Cell Number")

p2 <- ggplot(data = df2, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() + 
  scale_color_manual(values = c("Bb_cells" = "black", "Tb_cells" = "red", "Cbb_cells"="green", "Atb_cells" = "blue", "Ltb_cells" = "purple", 
                                "Ctb_cells" = "yellow")) +
  labs(title = "WithinHost Delay - Blood", color = "Cell Type") + theme(legend.position = "right") + theme_minimal() + xlab(label = "Time (Days)") + ylab(label = "Cell Number") + ylim(c(0,30)) #set ylim to see  


p3 <- ggplot(data = dtotal, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() + 
  scale_color_manual(values = c("B_total" = "black", "T_total" = "red")) +
  labs(title = "WithinHost Delay - Blood", color = "Cell Type") + theme(legend.position = "right") + theme_minimal() + xlab(label = "Time (Days)") + ylab(label = "Cell Number") + ylim(c(0,100))

(p1 | p2)/p3 + plot_layout(widths = c(2,2)) 


 