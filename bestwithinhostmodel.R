### NOTE: 9/14/2025 used this model for presentation for collaborators to illustrate adding delay to model improves fit ### 


rm(list = ls())
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(dplyr) 
library(readxl)
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
#get rid of macrophage and infect T cells and B cells population - could also have macrophages decay over time 
#blah blah adding someing  
# get rid of macrophages and add cB 
parameters_values <- c( 
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


#getting total B and T cell proportions 
sir_values_total <- sir_values_1 %>% mutate(total_lymphocyte = B_cells + T_cells )

#long for ggplot  
df <- melt(sir_values_1, id.vars = "time") %>% filter(!variable %in% c("Lt", "Lt2", "Lt4", "Lt3", "f", "If"))  

### DATA FOR PLOT ### 
baigent2016 <- read_xlsx("/Users/rayanhg/Downloads/baigent2016.xlsx", 2 ) 
baigent1998 <- read_xlsx("~/Desktop/WithinHostModel/WithinHostModel/baigent1998.xlsx", 4 ) 
baigent1998$mean.pp38 <- as.numeric(baigent1998$mean.pp38)
baigent1998_cyto <- baigent1998  %>% group_by(time) %>% summarise(total_pp38 = sum(mean.pp38, na.rm= TRUE))

df_feathers_only <- melt(sir_values_1, id.vars = "time") %>% filter(variable %in% c("If")) 


ggplot(data = df, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() + 
  scale_color_manual(values = c("B_cells" = "black", "T_cells" = "red", "Cb"="green", "At" = "blue", "Lt5" = "purple", 
                                "Ct" = "yellow", "Z" = "lightblue")) + 
  labs(title = "WithinHost Delay", color = "Cell Type") + theme(legend.position = "right") + 
  theme_minimal() + xlab(label = "Time (Days)") + ylab(label = "Cell Number") + 
  geom_point(data = baigent1998_cyto, aes(x = time/24, y = total_pp38), inherit.aes = FALSE, color = "red") + ylim(0,10000) 

 ggplot(data = df_feathers_only, aes(x = time/24, y = value, group = variable, colour = variable )) + geom_line() + 
  scale_color_manual(values = c( "If" = "pink")) + 
  geom_point(data = baigent2016, aes(x = time, y = mean_genomes), inherit.aes = FALSE, color = "red") + 
  labs(title = "WithinHost Delay", color = "Cell Type")  + xlab(label = "Time (Days)") + 
  ylab(label = "Cell Number")  +  scale_y_log10(limits = c(0.01, 10000000)) + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),   legend.text = element_text(size = 12),
        legend.title = element_text(size = 12), axis.line = element_line(color = "black")) + labs( y = "Cell Number", x = "Time (Days)", title = "WithinHost - Blood (B and T total)") 
 