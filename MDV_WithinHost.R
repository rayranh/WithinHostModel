library(deSolve)  

sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dB <- -M*B_cells - beta*Cb*B_cells - beta_2*Ct*B_cells + (g1*(Cb+Ct)/(g2+(Cb+Ct))) # beta = rate of cytolytically infected B cells by Cb and Ct 
    dCb <- M*B_cells + beta*Cb*B_cells + beta_2*Ct*B_cells - alpha*Cb # should alpha be larger here because cb keeps growing 
    dT <- -nu_A*Cb*T_cells - nu_b*Ct*T_cells + (h1*(Cb+Ct)/(h2+(Cb+Ct)))
    dAt <- nu_A*Cb*T_cells + nu_b*Ct*T_cells - beta_3*Ct*At - beta_4*Cb*At # beta = rate of activated T cells by Ct and Cb cells 
    dLt <- theta*(beta_3*Ct*At + beta_4*Cb*At) - mu*Lt
    dCt <- (1-theta)*(beta_3*Ct*At + beta_4*Cb*At) - alpha_2*Ct
    dZ <- mu*Lt
    df <- -nu_F*Lt*f
    dIf <- nu_F*Lt*f
    return(list(c(dB, dCb, dT, dAt, dLt, dCt, dZ, df,dIf)))
  }) 
}

parameters_values <- c( 
  M = 0.1
  , beta = 0.02898551           #contact rate with B cells (every 34 hours/ 1.4days) 
  , beta_2 = 0.00300            #contact rate with T cells (every 333 hour/ 13 days) 
   , beta_3 = beta_2              #contact rate with Ct cells for At 41 days activated T cells leaving 
   , beta_4 = beta              #contact rate with Cb cells for At 41 days activated T cells leaving
  , nu_A =  0.003               #Activation rate with T cells CD4+ (333 hours/ 13days)
  , nu_b = 0.006                #Activation rate with B cells (166 hours/ 41days)
  , nu_F = 0.006                #Infection rate of follicular cells (166 hours/ 41days)
  , mu = 0.01388889             #Rate of Tumor Cells (every 72 hours)
  , alpha = 0.03                #death rate of B cells (every 33 hours)
  , alpha_2 = 0.02083333        #death rate of T cells (every 48 hours)
  , theta = 0.7                 #population of activated T cells 
  , g1 = 0.06666667             #incoming B cells (every 15 hours)
  , g2 =0.001    
  , h1 = 0                      #incoming T cells / determined no incoming T cells 
  , h2 = 10
)

initial_values <- c( 
  B_cells = 100  
  , Cb = 0 
  , T_cells =100
  , At = 0
  , Lt = 0
  , Ct = 0
  , Z = 0
  , f =5 
  , If =0
  
) 

time_values <- seq(0, 200) # hours

sir_values_1 <- ode(
  y = initial_values,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
) 

sir_values_1 <- as.data.frame(sir_values_1)


with(sir_values_1, {
  plot(x=time, y=B_cells, col="black", type="l", ylim= c(0,100),  xlab="Time (Hours)", ylab="Population density", main="Marek's Model", cex.lab = 1.5, xlim = c(0,200))          #Plot the data for S over time
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
