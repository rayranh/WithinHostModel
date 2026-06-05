
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
library(patchwork)

## DATA ##  
#cytolytic infection at a given time of B and T cells in Spleen, Thymus, Bursa 

BCells2015 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Schermuly2015/schermuly_Bcells_fin.csv")%>% 
  mutate(MDV_load = MDV_load/100, uppersd = uppersd/100, lowersd = lowersd/100, CellType = "B_cells")

TCells2015 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Schermuly2015/schermuly_Tcells_fin.csv") %>%  
  mutate(MDV_load = MDV_load/100, uppersd = uppersd/100, lowersd = lowersd/100, CellType = "T_cells") %>% 
  filter(time != 72) 

CellsCombined <- bind_rows(BCells2015, TCells2015)


fit <- lm(log10(MDV_load) ~ time, data = BCells2015) 
summary(fit)

fit_Tcells<- lm(log10(MDV_load) ~ time, data = TCells2015)
summary(fit_Tcells) 

#Model with same and different slope 
same <- lm(log10(MDV_load) ~ CellType + time, data = CellsCombined)
different <- lm(log10(MDV_load) ~ CellType * time, data = CellsCombined)





simple <- logLik(same)
complicated <- logLik(different) 

library(lmtest) 
lrtest(same, different)

LR <- 2*(complicated - simple ) 
complicated/simple

pchisq(3.3457, df = 1, lower.tail = FALSE)

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


x_hours <- seq(1,100, by = 1)
same_equation_B <- data.frame(Values = 10^(3.30244 + 0.03578*x_hours), Time = x_hours) # have to transform the values back because did log10(MDV_load) for lm and i scalelogy10 on raw values not logged  
same_equation_T <- data.frame(Values = 10^(2.91951 + 0.03578*x_hours), Time = x_hours)    
diff_equation_B <- data.frame(Values = 10^(3.229445 + 0.038660*x_hours), Time = x_hours) 
diff_equation_T <- data.frame(Values = 10^(2.99251 + 0.032897*x_hours), Time = x_hours) 

SimpleModelPlot <- ggplot(CellsCombined, aes(x = time/24, y = MDV_load)) + geom_point(aes(y = MDV_load, colour = CellType)) + 
  geom_line(aes(y = MDV_load, colour = CellType)) + 
  scale_y_log10(limits =  c(0.01,100000000),breaks = c(0.01,1,100,10000,1000000,100000000)) + 
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5, by = 1)) + 
  scale_color_manual(values = c(B_cells = "red", T_cells = "blue", Model = "black"))+  
  geom_line(data = same_equation_B, aes(x = Time/24, y = Values, colour = "Model"), inherit.aes = FALSE) + 
  labs(title = "in vitro B&T Cell Data - simplified model", x = "time(days)", y = "MDV Genome Copies/10^4 Cells") +  
  theme(plot.title = element_text(size = 12)) + 
  geom_line(data = same_equation_T, aes(x = Time/24, y = Values, colour = "Model"), inherit.aes = FALSE)   

DataPlot <- ggplot(CellsCombined, aes(x = time/24, y = MDV_load)) + geom_point(aes(y = MDV_load, colour = CellType)) + 
  geom_line(aes(y = MDV_load, colour = CellType)) + 
  scale_y_log10(limits =  c(0.01,100000000),breaks = c(0.01,1,100,10000,1000000,100000000)) + 
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5, by = 1)) + 
  scale_color_manual(values = c(B_cells = "red", T_cells = "blue", Model = "black"))+   
  theme(plot.title = element_text(size = 12)) + 
  labs(title = "in vitro B&T Cell Data", x = "time(days)", y = "MDV Genome Copies/10^4 Cells")  


ComplexModelPlot <- ggplot(CellsCombined, aes(x = time/24, y = MDV_load)) + geom_point(aes(y = MDV_load, colour = CellType)) + 
  geom_line(aes(y = MDV_load, colour = CellType)) + 
  scale_y_log10(limits =  c(0.01,100000000),breaks = c(0.01,1,100,10000,1000000,100000000)) + 
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5, by = 1)) + 
  scale_color_manual(values = c(B_cells = "red", T_cells = "blue", Model = "black"))+  
  geom_line(data = diff_equation_B, aes(x = Time/24, y = Values, colour = "Model"), inherit.aes = FALSE) + 
  labs(title = "in vitro B&T Cell Data - complicated model", x = "time(days)", y = "MDV Genome Copies/10^4 Cells") +  
  theme(plot.title = element_text(size = 12)) + 
  geom_line(data = diff_equation_T, aes(x = Time/24, y = Values, colour = "Model"), inherit.aes = FALSE) 

DataPlot + SimpleModelPlot + ComplexModelPlot

