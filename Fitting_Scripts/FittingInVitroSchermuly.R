
#5.13.26 - Using this code to fit to in Vitro data using linear regression # 

rm(list = ls())
library(dplyr) 
library(tibble)
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(writexl) 
library(readxl)  
library(purrr)  
library(future) 
library(future.apply) 

## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 

BCells2015 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Schermuly2015/schermuly_Bcells_fin.csv")%>% 
  mutate(MDV_load = MDV_load/100, uppersd = uppersd/100, lowersd = lowersd/100)  

TCells2015 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Schermuly2015/schermuly_Tcells_fin.csv") %>%  
  mutate(MDV_load = MDV_load/100, uppersd = uppersd/100, lowersd = lowersd/100) %>% 
  filter(time != 72)


fit <- lm(MDV_load ~ time, data = BCells2015) 
summary(fit)

fit_Tcells<- lm(MDV_load ~ time, data = TCells2015)
summary(fit_Tcells)

ggplot(BCells2015,
       aes(x = time,
           y = MDV_load))+
  geom_point( color = "red") + scale_y_log10() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "black")

ggplot(TCells2015,
       aes(x = time,
           y = MDV_load))+
  geom_point( color = "red") + scale_y_log10() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "black")

