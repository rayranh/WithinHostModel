
######################
###### Events #######
####################
# x -> X + 1 , Rate = K 
# X -> X - 1 , Rate = gamma* X   

# dX/dt = k - gamma * x  

x <- 0 
t <- 0  

tend <- 1000 

k = 2 
gamma = 0.1  


while (tail(t,1) < tend) {

  current_x <- tail(x,1) 
  
  #gathering the rates, K is always K because it is incoming constantly, but because X (the population)
  #is changing and the rate at which they decay is gamma*X it is dynamic and can change because density dependent 
  rates = c(k, gamma*current_x)  
  
  rate_sum <- sum(rates) 
  
  #not 1/rate_sum because R uses rate not scale, this also creates a tau time that is between 0 and 1
  tau <- rexp(1,rate = rate_sum) 
  
  t<- c(t, tail(t,1) + tau)
  
  #generate random number between 0 & 1 
  rand <- runif(1,0,1)  
  
  # production event  
  #if the random number you generated * rate_sum which depends 
  # also on X is > 0 and if rand it is < rates[1] which is 2 (inflow rate) 
  # because dX/dt = k - gamma * x   when gamma*x is really small then K would be the max of what could be there 
  if (rand * rate_sum > 0 && rand * rate_sum < rates[1]) {  
    
    x<- c(x, current_x + 1)
    
  } 
  
  # decay event 
  # because dX/dt = k - gamma * x   when gamma*x is larger then K and gamma*X can be there 
  else if (rand * rate_sum > rates[1] && rand * rate_sum < rates[1] + rates[2]) { 
    x <-  c(x, current_x - 1)
    
    }
  
} 



ggplot() + geom_line(aes(x = t, y = x))

