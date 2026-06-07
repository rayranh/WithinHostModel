library(dplyr)
library(ggplot2)  
library(purrr)

set.seed(1)

u_diff <- c(6,5.9,5.8,5.7) 
sdv <- 0.5 
size <- 1000

observations <- map_dfr(u_diff, function(i){  
  map_df(1:10, function(rep){ 
    obs <- rnorm(size, i, sdv) # using i from outside function as putting in each new mean element of the vector 
    data.frame( 
      group = i,  
      obs = obs, 
      ID = rep
    )
    })

  }) 

#fit using optim to find mu 
#use LRT to compare the ML (fit mu) to other bersios one at a time  
# Cant use LRtest() get the difference in likelihood df = 1 , because in model fitting 1 and in other fitting 0  
# with  
#is value 6 about as good as fitted value, AND ANSWER WOULD BE YES OR NO DEPENDING ON BIGGER OR SMALLER THAN 0.05 



observationdf <- observations %>% group_by(group) %>% summarise(mean = mean(obs)) %>% ungroup()

observationdf_6 <- observationdf %>% filter(group == 6)

plot_obs <- ggplot(data = observationdf, aes(x = factor(group), y = mean, colour = factor(group) )) +
  geom_point() + coord_cartesian(ylim = c(2,7)) + 
  scale_y_continuous(breaks = seq(0,7,by=1)) + 
  scale_color_manual(values = c("6" = "orange", "5.7"= "red", "5.9" = "green", "5.8" = "blue" ))

plot_obs

#what is calculated likelihood using dnorm(u, sd)? 

likelihoodrepSD <- function(pars) { 
  
  likelihood <- dnorm(observationdf_6$mean, mean = pars[1], sd = sdv, log = TRUE)  
  sumLikelihood <- -sum(likelihood) 
  sumLikelihood 
  
} 

# this is my fitted model 
answeroptimSD <- optim(par = 6, fn = likelihoodrepSD,   lower = 4,
                     upper = 8, method = "Brent" ) 

# comparing fitted model to non-fitted model 
FixedUdiffLikelihood <- data.frame(fixedmu = u_diff) %>% 
  mutate(LikelihoodsFixedUdiff =  map_dbl(fixedmu, likelihoodrepSD), # take each element of vector and create a likelihood 
         LikelihoodOfFit = answeroptimSD$value, 
         LR = 2*(LikelihoodsFixedUdiff - answeroptimSD$value), # compare the likelihood to fitted optimized model 
         pvalue = pchisq(LR,df = 1, lower.tail = FALSE))
# if p < 0.05 then that means that the more complicated model is better   

FixedUdiffLikelihood 


#### now calculate using dnorm(u , sd/sqrt(1000)) , redo LRT  #### 
no_observations <- observations %>% group_by(group) %>% summarise(n = n()) %>% pull(n) %>% first()

likelihoodrepSE<- function(pars){ 
  
  likelihood <- dnorm(observationdf_6$mean, mean = pars[1], sd = sdv/sqrt(no_observations), log = TRUE)   
  sumLikelihood <- -sum(likelihood) 
  sumLikelihood 

}

answeroptimSE <- optim(par = 6, fn = likelihoodrepSE,   lower = 4,
                       upper = 8, method = "Brent" )  

answeroptimSE

FixedUdiffLikelihoodSE <- data.frame(fixedmu = u_diff) %>% 
  mutate(LikelihoodsFixedUdiffSE =  map_dbl(fixedmu, likelihoodrepSE), 
         LR = 2*(LikelihoodsFixedUdiffSE - answeroptimSE$value), 
         pvalue = pchisq(LR,df = 1, lower.tail = FALSE))

FixedUdiffLikelihoodSE 

# when n = 1000, you have more observations, and SE shows that the fitted model is better than not fitting 
# when n = 10 you have less observations and SD is still saying that simple model is okay and SE model is saying fitted is better 
# when n = 1 , you have even less observations and SE and SD are both saying that the simple model is okay  
