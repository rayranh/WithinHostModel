library(ggplot2)
library(dplyr)
library(readxl)


dataL6_It <- read_xlsx("km_data.xlsx", sheet = "Intratracheal_L6") 
dataL7_It <- read_xlsx("km_data.xlsx", sheet = "Intratracheal_L7") Survival

dosedata <- read_xlsx("survivalButter20072.xlsx", sheet = 1) 
dosedata100000L6<- read_xlsx("survivalButter20072.xlsx", sheet = 2)
dosedata10000L7<- read_xlsx("survivalButter20072.xlsx", sheet = "L710000")
dosedata100000L7<- read_xlsx("survivalButter20072.xlsx", sheet = "L7100000")

IntratrachealPlot <- ggplot(data = dataL6_It, aes(x = time, y = proportion, colour = "L6")) +  
  geom_step() + geom_step(data = dataL7_It, aes(x = time, y = proportion, colour = "L7" )) +   
  scale_color_manual(values = c("L6" = "blue", "L7" = "red")) + 
  scale_x_continuous(breaks = seq(0,70,by=10))  + coord_cartesian(ylim = c(0,1), xlim = c(0,70)) + 
  labs(colour = "Chicken Line", title = "Butter - Intratracheal Challenge Susceptible vs Resistant", x = "Time (days)", y = "Proportion")


DoseInoculationPlots <- 
  ggplot(data = dosedata, aes(x = time, y = prop, colour = "L6_10000")) + 
  geom_step() + coord_cartesian(ylim = c(0,1)) + 
  geom_step(data = dosedata100000L6, aes(x = time_2, y = prop_2, colour = "L6_100000")) + 
  geom_step(data = dosedata10000L7, aes(x = time, y = prop, colour = "L7_10000")) + 
  geom_step(data = dosedata100000L7, aes(x = time, y = prop, colour = "L7100000")) + 
  scale_color_manual(values = c( "L6_10000" = "blue", "L6_100000" = "purple",  "L7_10000" = "green",  "L7100000" = "seagreen4" )) +
  labs(colour = "Line/Dose", title = "Butter - Dose Challenge Susceptible vs Resistant", x = "Time (days)", y = "Proportion")



#### Bertzbach2022 #### 
MHCdata <- read_xlsx("SurvivalCurvesBertzbach2022.xlsx", sheet = "15_15Line") 
MHC15_19 <-read_xlsx("SurvivalCurvesBertzbach2022.xlsx", sheet = "15_19Line")  
MHC15_21 <-read_xlsx("SurvivalCurvesBertzbach2022.xlsx", sheet = "15_21Line")  
MHC19_21 <-read_xlsx("SurvivalCurvesBertzbach2022.xlsx", sheet = "19_21Line")  
MHC21_21 <-read_xlsx("SurvivalCurvesBertzbach2022.xlsx", sheet = "21_21Line")  


MHCPlots <- ggplot(data = MHCdata, aes(x = time, y = prop, colour = "MHC15_15")) + geom_step() + coord_cartesian(ylim = c(0,1)) + 
  geom_step(data = MHC15_19, aes(x = time, y = prop, colour = "MHC15_19")) + geom_step(data = MHC15_21, aes(x = time, y = prop, colour = "MHC15_21")) + 
  geom_step(data = MHC19_21, aes(x = time, y = prop, colour = "MHC19_21")) + scale_x_continuous(limits = c(0,91), breaks = seq(0,91,by=7)) + 
  geom_step(data = MHC21_21, aes(x = time, y = prop, colour = "MHC21_21")) + 
  scale_color_manual(values = c( "MHC15_15" = "blue", "MHC15_21" = "green", "MHC15_19"= "purple","MHC19_21"= "orange","MHC21_21"="red"  )) + 
  labs(title = "Bertzbach - Disease Incidence and MHC", y = "Disease Incidence %", x = "dpi", colour = "MHC Type")

BertzTumorInc <- read_xlsx("SurvivalCurvesBertzbach2022.xlsx" , sheet = "TumorIncidence") 

BertzNumOrgans <- read_xlsx("SurvivalCurvesBertzbach2022.xlsx" , sheet = "NumOrgans") %>% 
  pivot_wider(names_from = type, values_from = NumOrgans)

BertzTumorIncPlots <- ggplot(data = BertzTumorInc, aes(x = as.factor(line), y = prop)) + geom_col(aes(fill = line)) + 
  coord_cartesian(ylim = c(0,1)) + labs(title = "Bertzbach - Tumor Incidence and MHC", x = "MHC Type", y = "Tumor Incidence Proportion", fill = "MHC Type")  

BertzNumorgansPlots <- ggplot(data = BertzNumOrgans, aes(x = as.factor(Line), y = mean)) + geom_col(aes(fill = Line)) + 
  geom_errorbar(aes(ymin = mean, ymax = upper, width = 0.2)) + 
  coord_cartesian(ylim = c(0,5)) + labs(title = "Bertzbach - Average # of Organs and MHC",x = "MHC Type", y = "Average # of Organs")   

#### Hulten 2020 ####  

# assuming they hatch at 21 days embryonation can approximate age 
HultenSurvdat <- read_xlsx("SurvivalCurveDataHulten2020.xlsx", sheet = "RB1B_Surv") #in ovo given at 18 days embryonation, by time of hatch its been ~3 days, and 6 days they are outside of egg
HultenSurv9 <- read_xlsx("SurvivalCurveDataHulten2020.xlsx", sheet = "Sheet5")
HultenSurv4 <- read_xlsx("SurvivalCurveDataHulten2020.xlsx", sheet = "Rb1braw4")


HultenSurvPlotCombo <- ggplot(data = HultenSurvdat, aes(x = time, y = prop, colour = "Ia ~6dph")) + geom_step() + 
  coord_cartesian(ylim = c(0,1)) + labs(title = "Hulten-RB1B days post vaccination", y = "Percent Survival", x = "Time (days)") + 
  geom_step(data = HultenSurv9, aes(x = time, y = prop, colour = "Ia ~10dph"), inherit.aes = FALSE) +  
  geom_step(data = HultenSurv4, aes(x = time, y = prop, colour= "Ia ~5dph"), inherit.aes = FALSE) + 
  scale_color_manual(values = c("Ia ~6dph" = "red", "Ia ~10dph" = "blue", "Ia ~5dph" = "purple")) 

HultenSurv9Plot <- ggplot(data = HultenSurv9, aes(x = time, y = prop )) + geom_step() + 
  coord_cartesian(xlim = c(0,80), ylim = c(0,1)) 

HultenSurv4Plot <- ggplot(data = HultenSurv4, aes(x = time, y = prop )) + geom_step() + 
  coord_cartesian(xlim = c(0,80), ylim = c(0,1)) 

#### Baaten 2009 #### 

dose5mg <- read_xlsx("SurvivalBaaten2009.xlsx", sheet = "5mg")  
dose2.5mg <- read_xlsx("SurvivalBaaten2009.xlsx", sheet = "2.5mg") 
dose1mg <-  read_xlsx("SurvivalBaaten2009.xlsx", sheet = "1mg")
dose0.1mg <-read_xlsx("SurvivalBaaten2009.xlsx", sheet = "0.1mg")  
dose0.001mg <- read_xlsx("SurvivalBaaten2009.xlsx", sheet = "0.001mg")    
dose0.002mg <- read_xlsx("SurvivalBaaten2009.xlsx", sheet = "0.002mg")

#combining using dplyr bind_rows 
CombinedBaaten <- bind_rows(dose5mg %>% mutate(group = "5mg"), dose2.5mg %>%
                              mutate(group = "2.5mg"), dose1mg %>% mutate(group = "1mg"), 
                            dose0.1mg %>% mutate(group = "0.01mg"), dose0.001mg %>% mutate(group = "0.001"), 
                            dose0.002mg %>% mutate(group = "0.002")) 

Baaten2009Plot <- ggplot(data = CombinedBaaten, aes(x = time, y = prop , colour = group)) + 
  geom_step() + coord_cartesian(xlim = c(0,70))  

#### Singh 2010 ####

nonvaxSingh <- read_excel("SurvivalCurvesSingh2010.xlsx", sheet = "nonvax")
sb1vaxSingh <- read_excel("SurvivalCurvesSingh2010.xlsx", sheet = "Sb")

combinedSingh <- bind_rows(nonvaxSingh %>% mutate(group = "nonvax"), sb1vaxSingh %>% mutate(group = "Sb1"))

ggplot(data = combinedSingh, aes(x = time, y = prop, colour = group)) +
  geom_step() + coord_cartesian(xlim = c(0,100), ylim = c(0,1)) +
  labs(title = "Singh - RB1B Challenged SB1 Vaccinated Chickens", x = "Time (days)", y = "Percentage Survival") + 
  annotate("text", x = 92, y = 0.90, label = "MD lesions 11%", color = "black", size = 3) + 
  annotate("text", x = 40, y = 0, label = "MD lesions 100%", color = "black", size = 3)

  
