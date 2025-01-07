library(deSolve)  

sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dB <- -M*B_cells - beta*Cb*B_cells - beta_2*Ct*B_cells + (g1*(Cb+Ct)/(g2+(Cb+Ct))) 
    dCb <- M*B_cells + beta*Cb*B_cells + beta_2*Ct*B_cells - alpha*Cb
    dT <- -nu_A*Cb*T_cells - nu_b*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct)))
    dAt <- nu_A*Cb*T_cells + nu_b*Ct*T_cells - beta_2*Ct*At - beta*Cb*At
    dLt <- theta*(beta_2*Ct*At + beta*Cb*At) - mu*Lt
    dCt <- (1-theta)*(beta_2*Ct*At + beta*Cb*At) - alpha_2*Ct
    dZ <- mu*Lt
    df <- -nu_F*Lt*f
    dIf <- nu_F*Lt*f
    return(list(c(dB, dCb, dT, dAt, dLt, dCt, dZ, df,dIf)))
  }) 
}

parameters_values <- c( 
  M = 0.1
  , beta =4.819e-4              #contact rate with B cells 
  , beta_2 =5e-4                #contact rate with T cells #0.5 
  , nu_A = 0.005               #Activation rate with T cells CD4+  
  , nu_b = 0.01                #Activation rate with B cells 
  , nu_F = 0.001                 #Infection rate of follicular cells 
  , mu = 0.1                 #Rate of Tumor Cells 
  , alpha = 0.015              #death rate of B cells 
  , alpha_2 = 0.010            #death rate of T cells 
  , theta = 0.5                #population of activated T cells 
  , g1 = 0.05                  #incoming B cells 
  , g2 =0.001    
  , h1 = 0                     #incoming T cells 
  , h2 = 10
)

initial_values <- c( 
  B_cells = 50  
  , Cb = 0 
  , T_cells =50
  , At = 0
  , Lt = 0
  , Ct = 0
  , Z = 0
  , f =5 
  , If =0
  
) 

time_values <- seq(0, 100) # hours

sir_values_1 <- ode(
  y = initial_values,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
) 

sir_values_1 <- as.data.frame(sir_values_1)


with(sir_values_1, {
  plot(x=time, y=B_cells, col="black", type="l", ylim= c(0,50),  xlab="Time (Hours)", ylab="Population density", main="Marek's Model", cex.lab = 1.5)          #Plot the data for S over time
  lines(time, T_cells, col="red")  
  lines(time, Cb, col= "green")                                  #Add a line for I over time
  lines(time, At, col="blue")
  lines(time,Lt, col="purple" ) 
  lines(time,Ct, col="yellow" )  
  lines(time, Z, col= "lightblue") 
  #lines(time, f, col="aquamarine") 
  #lines(time,If, col="pink" ) 
})

# adding a legend:
legend("topright", c("B_cells", "T_cells", "Cb", "At", "Lt", "Ct","Z"),
       col = c("black", "red", "green","blue","purple","yellow", "lightblue"), lty = 1, bty = "n")

#82 - 93% of B cells infected 3 dpi  
