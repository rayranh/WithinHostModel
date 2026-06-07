library(dplyr)
library(tidyr) 
library(patchwork) 
library(readxl) 
library(ggplot2)



### BAIGENT 2016 ###  #10^4 cells 
data_feathers2016 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Baigent2016/Unvax/feathers_noVax_Fin_2.csv", header = TRUE) %>%
  rename(time = V1, mean = V2) 

data_feathers2016$time <- round(data_feathers2016$time)

### BAIGENT 2011 ###   #10^4 cells 
feathers2011 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Baigent2011/Unvax/Feather_NoVax_RB1B_fin.csv") 
feathers2011$time <- round(feathers2011$time)

### Singh 2010 ### 
Feathers2010 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Singh2010/Unvax/singh2010_RB1B_feathers_only_Fin.csv") 
Feathers2010$time <- round(Feathers2010$time)

### Spatz 2007 ###  
Feathers2007 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Spatz2007/No vax/Spatz2007Feathers_Fin.csv") 


ggplot(data = data_feathers2016, aes(x = time, y = mean, color = "Baigent2016")) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = mean, ymax =  ConfInt)) +  
  scale_y_log10(limits = c(0.01, 10000000000)) + 
  scale_x_continuous(limits = c(0,90), breaks = seq(0,90,by=10)) + 
  geom_point(data = feathers2011, aes(x = time, y = mean, color = "Baigent2011"), inherit.aes = FALSE) +  
  geom_line(data = feathers2011, aes(x = time, y = mean, color = "Baigent2011", linetype = "RIR"), inherit.aes = FALSE) + 
  geom_errorbar(data = feathers2011, aes(x = time, ymin = mean, ymax = ConfInt, color = "Baigent2011"), inherit.aes = FALSE) + 
  geom_point(data = Feathers2010, aes(x = time, y = mean, color = "Singh2010"), inherit.aes = FALSE) +  
  geom_errorbar(data = Feathers2010, aes(x = time, ymin = mean, ymax = ConfInt, color = "Singh2010"), inherit.aes = FALSE) + 
  geom_line(data = Feathers2010, aes(x = time, y = mean, color = "Singh2010", linetype = "B19/B19"), inherit.aes = FALSE) + 
  geom_point(data = Feathers2007, aes(x = time, y = MDV_load, color = "Spatz2007"), inherit.aes = FALSE) + 
  geom_line(data = Feathers2007, aes(x = time, y = MDV_load, color = "Spatz2007"), inherit.aes = FALSE) + 
  geom_errorbar(data = Feathers2007, aes(x = time, ymin = lowersd, ymax = uppersd, color = "Spatz2007"), inherit.aes = FALSE, width = 4) + 
  scale_color_manual( values = c("Baigent2016" = "blue", "Baigent2011" = "purple", "Singh2010" = "seagreen", "Spatz2007" = "red"))  +  
  scale_linetype_manual(values = c("B19/B19" = "dotted", "RIR" = "solid"))  + 
  labs(title = "RB1B - Feather Follicle Infection Data", y = "RB1B Genomes per 10,000 Cells", color = "study", linetype = "Chicken Type ", x = "time (dpc)")
  
