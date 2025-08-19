rm(list = ls()) 

library(deSolve) 
library(dplyr) 


data <- as.data.frame(read.csv("~/Desktop/Modelling_Practice/Fitting-SIR-/SIR_FractionInfected_With_Noise_v2.csv"))  

data <- data %>% rename(time = Time) %>% mutate(fracInfected = SampledInfected/SampleSize)


SIR_odes <- function (t,x,params){ 
  S <- x[1] 
  I <- x[2]
  R <- x[3]  
  beta <- params[1] 
  gamma <- params[2] 
  dS <- -beta*S*I 
  dI <- beta*S*I - gamma*I 
  dR <- gamma*I 
  list(c(dS,dI,dR))
  
}

#Initial state parameters 
initial_state <- c(S=999/1000, I=1/1000, R=0)

#times 
t<- seq(0,150,by= 1) # time to solve ode 
timeValues <- seq(5,150,by=5) # time that matches data to subset 
#parameters 
params<- c(0.14,0.05)


### putting inputs to get model likelihood 
MximumLikelihood <- function(params){  
  beta <- params[1] 
  gamma <- params[2] 
  results <- as.data.frame(ode(y = initial_state, 
                               times = t, func = SIR_odes, parms =c(beta = beta, gamma = gamma))) 
  subset_I <- results %>%
    filter(time %in% timeValues) %>%
    select(time,I) 
  
  merged_data <- data %>%
    left_join(subset_I, by = "time")
  
  MaxLike <- merged_data %>% mutate(logLike = log(dbinom(SampledInfected, SampleSize, I))) 
  
  
  sumLog <- -sum(MaxLike$logLike)
  
  return(sumLog)
}

MximumLikelihood(params)


df_Model <- as.data.frame(ode(y = initial_state, 
                              times = t, func = SIR_odes, parms = params))


plot(x =data$time, y= data$fracInfected, type = 'p', col = 'red', ylim = c(0,1), main = 'sir', xlab = 'time', ylab = 'number of infected people')
lines(x = df_Model$S , type = 'l', col= 'black') 
lines(x = df_Model$I, type = 'l', col= 'green')
lines(x = df_Model$R,  type = 'l', col = 'blue')  
legend("topright", title="key", legend=c("susceptible", "infected", "recovered"), 
       col = c("black","green", "blue"), cex = 0.5,
       lty=c(1,1,1)) 


print(optim(par = params, fn = MximumLikelihood, method = "Nelder-Mead"))
