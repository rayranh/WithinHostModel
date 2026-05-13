#### This code is plotting all RB1B within host data along with some Feather Follicle data. 
#### The Schermuly in vitro data is included in this code as well in addition to plaque assays  #### 



library(dplyr)
library(tidyr) 
library(patchwork) 
library(readxl) 
library(ggplot2)

### BAIGENT 2016 ### 
data_feathers2016 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Baigent2016/Unvax/feathers_noVax_Fin_2.csv", header = TRUE) 
data_pbl2016 <- read.csv( "~/Desktop/WithinHostModel/DataForProject/Baigent2016/Unvax/PBL_No_Vax_fin.csv")

data_feathers2016$V1 <- round(data_feathers2016$V1)
data_pbl2016$time <- round(data_pbl2016$time)


# Og scale - 10000000, x = (-5,40)
FeathersNoVax2016 <- ggplot(data = data_feathers2016, aes(x = V1, y = V2)) + 
  geom_point() + scale_y_log10() + 
  geom_linerange( aes(ymin = V2, ymax = ConfInt)) + 
  geom_errorbar(aes(ymin= ConfInt, ymax = ConfInt), width = 0.5) + 
  geom_line() + coord_cartesian(xlim = c(-5,90)) + scale_x_continuous(breaks = seq(-5,90, by = 5) ) + 
  labs(title = "Feathers RB1B - No Vaccine ", y = "RB1B Genomes per 10^4" , x = "time(days)") 

PBLNoVax2016 <- ggplot(data = data_pbl2016, aes(x = time, y = mean)) + 
  geom_point() + scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) + 
  geom_linerange( aes(ymin = mean, ymax = ConfInt)) + 
  geom_errorbar(aes(ymin= ConfInt, ymax = ConfInt), width = 0.5) + 
  geom_line() + coord_cartesian(xlim = c(-5,90), ylim = c(0.01,100000000)) + scale_x_continuous(breaks = seq(-5,90, by = 5) ) + 
  labs(title = "Baigent2016 - PBL RB1B - No Vaccine " , y = "RB1B Genomes per 10^4", x = "time(days)")  


### BAIGENT 2011 ### 

kidney2011 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Baigent2011/Unvax/Kidney_NoVax_RB1B_fin.csv") 
spleen2011 <-  read.csv("~/Desktop/WithinHostModel/DataForProject/Baigent2011/Unvax/Spleen_noVax_RB1B_fin.csv") 
liver2011 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Baigent2011/Unvax/Liver_NoVax_RB1B_fin.csv") 
feathers2011 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Baigent2011/Unvax/Feather_NoVax_RB1B_fin.csv")

kidney2011$time <- round(kidney2011$time)
spleen2011$time <- round(spleen2011$time)
liver2011$time <- round(liver2011$time)
feathers2011$time <- round(feathers2011$time)


# og y = (0.01,10000000), x = (-5,65)
KidneyNoVax2011 <- ggplot(data = kidney2011, aes(x = time, y = mean)) + 
  geom_point() + scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) +  
  geom_line(data = kidney2011 %>% slice(1:n()-1)) + 
  geom_linerange( aes(ymin = mean, ymax = ConfInt)) + 
  geom_errorbar(aes(ymin= ConfInt, ymax = ConfInt), width = 0.5) +
  coord_cartesian(xlim = c(-5,90), ylim = c(0.01,100000000)) + scale_x_continuous(breaks = seq(-5,90, by = 5)) + 
  labs(title = "Baigent2011 - Kidney RB1B - No Vaccine ", y = "RB1B Genomes per 10^4", x = "time(days)") 

SpleenNoVax2011 <- ggplot(data = spleen2011, aes(x = time, y = mean)) + 
  geom_point() + scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) +  
  geom_line(data = spleen2011 %>% slice(1:n()-1)) + 
  geom_linerange( aes(ymin = mean, ymax = ConfInt)) + 
  geom_errorbar(aes(ymin= ConfInt, ymax = ConfInt), width = 0.5) +
  coord_cartesian(xlim = c(-5,90), ylim = c(0.01,100000000) ) + scale_x_continuous(breaks = seq(-5,90, by = 5)) + 
  labs(title = "Baigent2011 - Spleen RB1B - No Vaccine ", y = "RB1B Genomes per 10^4", x = "time(days)") 

LiverNoVax2011 <- ggplot(data = liver2011, aes(x = time, y = mean)) + 
  geom_point() + scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) +  
  geom_line(data = liver2011 %>% slice(1:n()-1)) + 
  geom_linerange( aes(ymin = mean, ymax = ConfInt)) + 
  geom_errorbar(aes(ymin= ConfInt, ymax = ConfInt), width = 0.5) +
  coord_cartesian(xlim = c(-5,90),ylim = c(0.01,100000000)) + scale_x_continuous(breaks = seq(-5,90, by = 5)) + 
  labs(title = "Baigent 2011 - Liver RB1B - No Vaccine ", y = "RB1B Genomes per 10^4", x = "time(days)") 

FeathersNoVax2011 <- ggplot(data = feathers2011, aes(x = time, y = mean)) + 
  geom_point() + scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) +  
  geom_line(data = feathers2011 %>% slice(1:n()-1)) + 
  geom_linerange( aes(ymin = mean, ymax = ConfInt)) + 
  geom_errorbar(aes(ymin= ConfInt, ymax = ConfInt), width = 0.5) +
  coord_cartesian(xlim = c(-5,90), ylim = c(0.01,100000000)) + scale_x_continuous(breaks = seq(-5,90, by = 5)) + 
  labs(title = "Baigent 2011 - Feathers RB1B - No Vaccine ", y = "RB1B Genomes per 10^4", x = "time(days)") 


### Luo 2020 ###

CEF2020 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Luo2020/Novax/LuoCEF2_Fin.csv") %>% mutate(time = time/24)

CEF2020 <- CEF2020 %>%
  mutate(log10_load = log10(MDV_load),
         log10_lower = log10(lowersd),
         log10_upper = log10(uppersd))

# CEFNoVax2020 <- ggplot(data = CEF2020, aes(x = time, y = MDV_load)) + 
#   geom_point() + geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 1) + geom_line() + 
#   coord_cartesian(xlim = c(24,120)) + scale_x_continuous(breaks = seq(24,144, by = 24)) +
#   scale_y_continuous(breaks = seq(20000,120000, by = 20000)) + ylim(0,120000) + 
#   labs(title = "CEF - No Vaccine ", y = "RB1B Genomes per 10^4", x = "time(hours)")   

## CHANGED CEF TO LOG 

CEFNoVax2020 <- ggplot(data = CEF2020, aes(x = time, y = MDV_load)) + 
  geom_point() + geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.05) + geom_line() +  
  coord_cartesian(xlim = c(-5,90), ylim = c(0.01,100000000)) + 
  scale_x_continuous(breaks = seq(-5,90, by = 5)) +
  scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) +
  labs(title = "Luo2020 - CEF - No Vaccine ", y = "RB1B Genomes per 10^4", x = "time(days)")   


#### Conradie 2020 #### 
PBLConradie2020 <- read_excel("/Users/rayanhg/Desktop/WithinHostModel/DataForProject/Conradie2020/PBL_Genome_RB1B_Conradie2020.xlsx") %>% 
  pivot_wider(names_from = category, values_from = genome) %>% 
  mutate(mean = mean/100, upper = upper/100, lower = lower/100)



### Parcells 2001 ### 

#measured in pfu 
PBL2001 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Parcells2001/Parcells_PBL_fin.csv") %>% mutate(time = time * 7)
PBL2001_WL <- read.csv("~/Desktop/WithinHostModel/DataForProject/Parcells2001/Parcells_PBL_Pfu_whiteLeghorn_fin.csv") %>% mutate(time = time * 7)
Spleen2001 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Parcells2001/Parcells_Spleen_PFU_fin.csv") %>% mutate(time = time * 7)
Spleen2001_WL <- read.csv("~/Desktop/WithinHostModel/DataForProject/Parcells2001/Parcells_Spleen_PFU_whiteleghorn_fin.csv") %>% mutate(time_wk = time_wk * 7)
 

#white leghorns used a lower passage number RB1B ~ 12 vs B15 which used ~21 / theyre main objective was to compare to a vIL8 mutant not to compare to each other 

PBL2001NoVax <- ggplot(data = PBL2001, aes(x = time, y = pfuxmil)) + 
  geom_point() + geom_line() + 
  labs(title = "Parcells2001 - PBL Plaque Formation - B15 x 7 Chickens", y = "RB1B Plaques per 10^6 Plated Cells", x = "time(days)") + 
  scale_x_continuous(breaks = seq(0,28, by = 7)) + 
  scale_y_continuous(limits = c(0,1500), breaks = seq(0,1500, by = 250))

PBL2001WLNoVax <-  ggplot(data = PBL2001_WL, aes(x = time, y = pblxmil)) + geom_point() + 
  geom_line() + labs(title = "Parcells2001 - PBL Plaque Formation - White Leghorn Chickens", y = "RB1B Plaques per 10^6 Plated Cells", x = "time(days)") + 
  scale_y_continuous(limits = c(0, 1500), breaks = seq(0,1500, by = 250))  + 
  scale_x_continuous(limits = c(0, 28) , breaks = seq(0,28, by = 7)) 

Spleen2001NoVax <-  ggplot(data = Spleen2001, aes(x = time, y = PFUxmil)) + geom_point() + 
  geom_line() + labs(title = "Parcells2001 - Spleen Plaque Formation - B15x7 Chickens", y = "RB1B Plaques per 10^6 Plated Cells", x = "time(days)") + 
  scale_y_continuous(limits = c(0, 1500), breaks = seq(0,1500, by = 250)) +   scale_x_continuous(limits = c(0, 28) , breaks = seq(0,28, by = 7)) 


Spleen2001WLNoVax <-  ggplot(data = Spleen2001_WL, aes(x = time_wk, y = pfuxmil)) + geom_point() + 
  geom_line() + labs(title = "Parcells2001 - Spleen Plaque Formation - White Leghorn Chickens", y = "RB1B Plaques per 10^6 Plated Cells", x = "time(days)") + 
  scale_y_continuous(breaks = seq(0,1500, by = 250)) + 
  scale_x_continuous(limits = c(0, 28) , breaks = seq(0,28, by = 7)) 

### Sabsabi 2024 ### 

bursa2024 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Sabsabi2024/SabsabiWtbursa_Fin.csv") 
spleen2024 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Sabsabi2024/SabsabiWtSpleen_Fin.csv") 
thymus2024 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Sabsabi2024/SabsabiWtThymus_Fin.csv") 
PBL2024 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Sabsabi2024/sabsabi_PBL_WT_fin.csv")


PBL2024$time <- c(4,7,10,14,21,28)

bursa2024plot <- ggplot(data = bursa2024, aes(time, MDV_load)) + #factor turns time categorical 
  geom_point() + scale_y_log10(limits =  c(0.01,100000000) , breaks = c(0.01,1,100,10000,1000000,100000000)) + # not plotting because goes below 1 so changed to 0.95  
  scale_x_continuous(limits = c(-5,75), breaks = seq(-5,75, by = 5)) + 
  geom_linerange( aes(ymin = lowersd, ymax = uppersd)) + 
  geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.05) + 
  labs(title = "Sabsabi2024 - Bursa" , y = "RB1B Genomes per 10^6", x = "time(days)")  

spleen2024plot <-  ggplot(data = spleen2024, aes(time, MDV_load)) + #factor turns time categorical 
  geom_point() + scale_y_log10(limits =  c(0.01,100000000) , breaks = c(0.01,1,100,10000,1000000,100000000))+ 
  scale_x_continuous(limits = c(-5,75), breaks = seq(-5,75, by = 5)) + # not plotting because goes below 1 so changed to 0.95 
  geom_linerange( aes(ymin = lowersd, ymax = uppersd)) + 
  geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.05) + 
  labs(title = "Sabsabi2024 - Spleen" , y = "RB1B Genomes per 10^6", x = "time(days)")  

thymus2024plot <-  ggplot(data = thymus2024, aes(time, MDV_load)) + #factor turns time categorical 
  geom_point() + scale_y_log10(limits =  c(0.01,100000000) , breaks = c(0.01,1,100,10000,1000000,100000000)) + # not plotting because goes below 1 so changed to 0.95 
  scale_x_continuous(limits = c(-5,75), breaks = seq(-5,75, by = 5)) + 
  geom_linerange( aes(ymin = lowersd, ymax = uppersd)) + 
  geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.05) + 
  labs(title = "Sabsabi2024 - Thymus" , y = "RB1B Genomes per 10^6", x = "time(days)")  

PBL2024Plot <- ggplot(data = PBL2024, aes(time, MDV_load)) + #factor turns time categorical 
  geom_point() + geom_line() +  scale_y_log10(limits =  c(0.01,100000000) , breaks = c(0.01,1,100,10000,1000000,100000000))+ # not plotting because goes below 1 so changed to 0.95 
  scale_x_continuous(limits = c(-5,75), breaks = seq(-5,75, by = 5)) + 
  geom_linerange( aes(ymin = lowersd, ymax = uppersd)) + 
  geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.5) + 
  labs(title = "Sabsabi2024 - PBL" , y = "RB1B Genomes per 10^6", x = "time(days)")

### Schermuly 2015 ### 

BCells2015 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Schermuly2015/schermuly_Bcells_fin.csv") %>% mutate(time = time/24) %>% 
  mutate(MDV_load = MDV_load/100, uppersd = uppersd/100, lowersd = lowersd/100) #scaling per 10^4 cells instead of 10^6 
TCells2015 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Schermuly2015/schermuly_Tcells_fin.csv") %>% mutate(time = time/24) %>% 
  mutate(MDV_load = MDV_load/100, uppersd = uppersd/100, lowersd = lowersd/100)

#dont know if error bars are sd or something else 
BCells2015Plot <- ggplot(data = BCells2015, aes(time, MDV_load)) +
  geom_point() + geom_line() + scale_y_log10(limits =  c(0.01,100000000) , breaks = c(0.01,1,100,10000,1000000,100000000)) +  
  scale_x_continuous(limits = c(0,20), breaks = seq(0,20, by = 5)) + 
  geom_linerange( aes(ymin = lowersd, ymax = uppersd)) + 
  geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.5) + 
  labs(title = "B Cell In Vitro Infection" , y = "RB1B Genomes per 10^6", x = "time(days)") 

TCells2015Plot <- ggplot(data = TCells2015, aes(time, MDV_load)) +
  geom_point() + geom_line() + scale_y_log10(limits =  c(0.01,100000000) , breaks = c(0.01,1,100,10000,1000000,100000000)) + 
  scale_x_continuous(limits = c(0,20), breaks = seq(0,20, by = 5)) + 
  geom_linerange( aes(ymin = lowersd, ymax = uppersd)) + 
  geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.5) + 
  labs(title = "T Cell In Vitro Infection" , y = "RB1B Genomes per 10^6", x = "time(days)") 


### Singh 2010 ### 

PBL2010 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Singh2010/Unvax/Rb1b_only_PBL_Fin.csv") 
Feathers2010 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Singh2010/Unvax/singh2010_RB1B_feathers_only_Fin.csv") 


PBL2010$time <- round(PBL2010$time) 
Feathers2010$time <- round(Feathers2010$time)


PBL2010Plot <- ggplot(data = PBL2010, aes(time, mean)) +
  geom_point() + geom_line() + scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) +  
  coord_cartesian(xlim = c(-5,90), ylim = c(0.01,100000000)) +  
  scale_x_continuous(breaks = seq(-5,90, by = 5) )  + 
  geom_linerange( aes(ymin = mean, ymax = ConfInt)) + 
  geom_errorbar(aes(ymin= ConfInt, ymax = ConfInt), width = 0.5) + 
  labs(title = "Singh2010 - PBL RB1B" , y = "RB1B Genomes per 10^4", x = "time(days)") 

Feathers2010Plot <- ggplot(data = Feathers2010, aes(time, mean)) +
  geom_point() + geom_line() + scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) +  
  coord_cartesian(xlim = c(-5,90), ylim = c(0.01,100000000)) + 
  geom_linerange( aes(ymin = mean, ymax = ConfInt)) + 
  geom_errorbar(aes(ymin= ConfInt, ymax = ConfInt), width = 0.5) + 
  labs(title = "Singh2010 - Feathers RB1B" , y = "RB1B Genomes per 10^4", x = "time(days)") + 
  scale_x_continuous(limits = c(-5,90), breaks = seq(-5,90,by=5))  


### Spatz 2007 ###  

### Maybe CEF infected CEF cells make a lot of plaques, dont include 

InVitro2007 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Spatz2007/No vax/spatz_invitro_Rb1b_Fin.csv")
PBL2007 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Spatz2007/No vax/spatz_Pbl_Rb1b_Fin.csv") 
Feathers2007 <- read.csv("~/Desktop/WithinHostModel/DataForProject/Spatz2007/No vax/Spatz2007Feathers_Fin.csv") 



InVitro2007$time <- round(InVitro2007$time)  
InVitro2007 <- InVitro2007 %>% mutate(time = time/24)

PBL2007$time <- round(PBL2007$time) 

InVitro2007Plot <- ggplot(data = InVitro2007, aes(time, pfu)) +
  geom_point() + geom_line() + coord_cartesian(xlim = c(0, 28), ylim = c(0, 100000))  +
  scale_x_continuous(limit = c(0,28), breaks = seq(0,28,by = 7)) + 
  labs(title = "Spatz2007 - In Vitro CEF RB1B" , y = "PFU per mL", x = "time(days)")  


PBL2007Plot <- ggplot(data = PBL2007, aes(time, MDV_load)) +
  geom_point() + geom_line() + scale_y_log10( limit = c(0.01,100000000), breaks = c(0.01,1,100,10000,1000000,100000000))  + 
  geom_linerange( aes(ymin = LowerConfInf, ymax = UpperConfInt)) + 
  geom_errorbar(aes(ymin= LowerConfInf, ymax = UpperConfInt), width = 0.5) + 
  labs(title = "Spatz2007 - PBL RB1B" , y = "RB1B Genomes per 10^4", x = "time(Days)") + 
  scale_x_continuous(limits = c(-5,90), breaks = seq(-5,90,by=5)) 

Feathers2007Plot <- ggplot(data = Feathers2007, aes(time, MDV_load)) +
  geom_point() + geom_line() +  scale_y_log10( limit = c(0.01,100000000), breaks = c(0.01,1,100,10000,1000000,100000000)) + 
  geom_linerange( aes(ymin = lowersd, ymax = uppersd)) +                  ### REDO 
  geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.9) + 
  labs(title = "Spatz2007 - Feathers" , y = "RB1B Genomes per 10^4", x = "time(Days)") + 
  scale_x_continuous(limits = c(-5,90), breaks = seq(-5,90,by=5))   


#### Berthault2018 #### 
BerthaultThymus2018 <- read_excel("/Users/rayanhg/Desktop/WithinHostModel/DataForProject/Berthault2018/RB1B_Genomes_Organs_Berthault2018.xlsx", sheet = "Thymus") %>%  
  mutate(genomeper10k = genomepermil/100) %>% group_by(time) %>% mutate(median = median(genomeper10k))

BerthaultBursa2018 <- read_excel("/Users/rayanhg/Desktop/WithinHostModel/DataForProject/Berthault2018/RB1B_Genomes_Organs_Berthault2018.xlsx", sheet = "Bursa") %>%  
  mutate(genomeper10k = genomepermil/100) 

BerthaultPBMC2018 <- read_excel("/Users/rayanhg/Desktop/WithinHostModel/DataForProject/Berthault2018/RB1B_Genomes_Organs_Berthault2018.xlsx", sheet = "PBMC") %>%  
  mutate(genomeper10k = genomepermil/100) 

BerthaultSpleen2018 <- read_excel("/Users/rayanhg/Desktop/WithinHostModel/DataForProject/Berthault2018/RB1B_Genomes_Organs_Berthault2018.xlsx", sheet = "Spleen") %>%  
  mutate(genomeper10k = genomepermil/100)
  

##################### 
## Genomes per 10K ## 
####################
PBL2010Plot + PBLNoVax2016 +  PBL2007Plot +  KidneyNoVax2011 + SpleenNoVax2011 + LiverNoVax2011 + CEFNoVax2020 

GenomePer10KPlots <- ggplot(data = PBL2010, aes(time, mean, color = "PBL"))  + 
  geom_point() + geom_line() + scale_y_log10(breaks = c(0.01,1,100,10000,1000000,100000000)) +   
  geom_point(data = data_pbl2016, aes(x = time, y = mean, color = "PBL")) + 
  geom_linerange( data = data_pbl2016, aes(ymin = mean, ymax = ConfInt, color = "PBL")) +  
  geom_line(data = data_pbl2016, aes(x = time, y = mean, color = "PBL", linetype = "RIR")) + 
  geom_errorbar(data = data_pbl2016, aes(ymin= ConfInt, ymax = ConfInt, color = "PBL"), width = 0.5) +  
  geom_point(data = PBL2007, aes(time, MDV_load, color = "PBL"), inherit.aes = FALSE) +  
  geom_linerange(data = PBL2007, aes(x = time, ymin = LowerConfInf, ymax = UpperConfInt, color = "PBL"), inherit.aes = FALSE) + 
  geom_errorbar( data = PBL2007, aes(x = time, ymin= LowerConfInf, ymax = UpperConfInt, color = "PBL"), inherit.aes = FALSE, width = 0.5) +  
  geom_line(data = PBL2007, aes(time, MDV_load, color = "PBL", linetype = "Line P"), inherit.aes = FALSE) + 
  geom_point(data = kidney2011, aes(time, y = mean, color = "Kidney"), inherit.aes = FALSE) +  
  geom_line(data = kidney2011 %>% slice(1:n()-1), aes(time, y = mean, color = "Kidney"), inherit.aes = FALSE) + 
  geom_linerange( data = kidney2011, aes(time, ymin = mean, ymax = ConfInt, color = "Kidney"), inherit.aes = FALSE) + 
  geom_errorbar(data = kidney2011, aes(time, ymin= ConfInt, ymax = ConfInt, color = "Kidney"), width = 0.5, inherit.aes = FALSE) + 
  geom_point(data = spleen2011, aes(x = time, y = mean, color = "Spleen"), inherit.aes = FALSE) +  
  geom_line(data = spleen2011 %>% slice(1:n()-1), aes(time, mean, color = "Spleen"), inherit.aes = FALSE) + 
  geom_linerange( data = spleen2011, aes(time, ymin = mean, ymax = ConfInt, color = "Spleen"), inherit.aes = FALSE) + 
  geom_errorbar(data = spleen2011, aes(time, ymin= ConfInt, ymax = ConfInt, color = "Spleen"), width = 0.5, inherit.aes = FALSE) + 
  geom_point(data = liver2011, aes(x = time, y = mean, color = "Liver"), inherit.aes = FALSE) +   
  geom_line(data = liver2011 %>% slice(1:n()-1), aes(time, mean, color = "Liver"), inherit.aes = FALSE) + 
  geom_linerange( data = liver2011, aes(time, ymin = mean, ymax = ConfInt, color = "Liver"), inherit.aes= FALSE) + 
  geom_errorbar(data = liver2011, aes(time, ymin= ConfInt, ymax = ConfInt, color = "Liver"), inherit.aes = FALSE, width = 0.5) + 
  geom_point(data = CEF2020, aes(x = time, y = MDV_load, color = "CEF"), inherit.aes = FALSE) +  
  geom_errorbar(data = CEF2020, aes(time, ymin= lowersd, ymax = uppersd, color = "CEF"), width = 0.05, inherit.aes = FALSE) +
  geom_line(data = CEF2020, aes(x = time, y = MDV_load, color = "CEF"), inherit.aes = FALSE) + 
  geom_point(data = BCells2015, aes(time, MDV_load, color = "In Vitro B"), inherit.aes = FALSE) + 
  geom_line(data = BCells2015, aes(time, MDV_load, color = "In Vitro B"), inherit.aes = FALSE) + 
  geom_linerange(data = BCells2015, aes(x = time, ymin = lowersd, ymax = uppersd, color = "In Vitro B"), inherit.aes = FALSE) +
  geom_errorbar(data = BCells2015, aes(x = time, ymin = lowersd, ymax = uppersd, color = "In Vitro B"), inherit.aes = FALSE) +  
  geom_point(data = TCells2015, aes(time, MDV_load, color = "In Vitro T"), inherit.aes = FALSE) + 
  geom_line(data = TCells2015, aes(time, MDV_load, color = "In Vitro T"), inherit.aes = FALSE) + 
  geom_linerange(data = TCells2015, aes(x = time, ymin = lowersd, ymax = uppersd, color = "In Vitro T"), inherit.aes = FALSE) +
  geom_errorbar(data = TCells2015, aes(x = time, ymin = lowersd, ymax = uppersd, color = "In Vitro T"), inherit.aes = FALSE) + 
  coord_cartesian(ylim = c(0.01,100000000)) + geom_point(data = PBLConradie2020, aes(x = time, y = mean, color = "PBL"), inherit.aes = FALSE) + 
  geom_linerange(data = PBLConradie2020, aes(x = time, ymin = lower, ymax = upper, color = "PBL"), inherit.aes = FALSE) +  
  geom_line(data = PBLConradie2020, aes(x = time, y =mean, color = "PBL"), inherit.aes = FALSE) + 
  geom_errorbar(data = PBLConradie2020, aes(x = time, ymin= lower, ymax = upper, color = "PBL"), width = 0.5, inherit.aes = FALSE) + 
  scale_color_manual(name = "Organs", values = c("PBL" = "slateblue1", "Kidney" = "aquamarine", "Spleen" = "darkorchid2", "Liver" = "orange", 
                                                 "CEF" = "palegreen4", "In Vitro B" = "red", "In Vitro T" = "pink", "Bursa" = "green", "Thymus" = "blue" )) +  
  # geom_point(data = BerthaultThymus2018, aes(x = time, y = genomeper10k, color = "Thymus"), inherit.aes = FALSE) +   
  stat_summary(data = BerthaultThymus2018, aes(x = time, y = genomeper10k, colour = "Thymus", linetype = "WhiteLeghorn B19/19"), fun = "median", geom = "line", inherit.aes = FALSE) +  
  # geom_point(data = BerthaultSpleen2018, aes(x = time, y = genomeper10k, color = "Spleen"), inherit.aes = FALSE) + 
  stat_summary(data = BerthaultSpleen2018, aes(x = time, y = genomeper10k, colour = "Spleen", linetype = "WhiteLeghorn B19/19"), fun = "median", geom = "line", inherit.aes = FALSE) + 
  # geom_point(data = BerthaultPBMC2018, aes(x = time, y = genomeper10k, color = "PBL"), inherit.aes = FALSE) +  
  stat_summary(data = BerthaultPBMC2018, aes(x = time, y = genomeper10k, colour = "PBL", linetype = "WhiteLeghorn B19/19"), fun = "median", geom = "line", inherit.aes = FALSE) +
  # geom_point(data = BerthaultBursa2018, aes(x = time, y = genomeper10k, color = "Bursa"), inherit.aes = FALSE) +  
  stat_summary(data = BerthaultBursa2018, aes(x = time, y = genomeper10k, colour = "Bursa", linetype = "WhiteLeghorn B19/19"), fun = "median", geom = "line", inherit.aes = FALSE) +  
  scale_linetype_manual(name = "Chicken Type" , values = c("WhiteLeghorn B19/19" = "dashed", "RIR" = "solid", "Line P" = "dotted")) + 
  scale_x_continuous(limits = c(0,90), breaks = seq(0,90, by = 5) ) + labs(title = "RB1B WithinHost Data" , y = "RB1B Genomes per 10^4", x = "time(days)")  
  
  
   #dont know if error bars are sd or something else 
  BCells2015Plot <- ggplot(data = BCells2015, aes(time, MDV_load)) +
    geom_point() + geom_line() + scale_y_log10(limits =  c(0.01,100000000) , breaks = c(0.01,1,100,10000,1000000,100000000)) +  
    scale_x_continuous(limits = c(0,20), breaks = seq(0,20, by = 5)) + 
    geom_linerange( aes(ymin = lowersd, ymax = uppersd)) + 
    geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.5) + 
    labs(title = "B Cell In Vitro Infection" , y = "RB1B Genomes per 10^4", x = "time(days)") 
  
  TCells2015Plot <- ggplot(data = TCells2015, aes(time, MDV_load)) +
    geom_point() + geom_line() + scale_y_log10(limits =  c(0.01,100000000) , breaks = c(0.01,1,100,10000,1000000,100000000)) + 
    scale_x_continuous(limits = c(0,20), breaks = seq(0,20, by = 5)) + 
    geom_linerange( aes(ymin = lowersd, ymax = uppersd)) + 
    geom_errorbar(aes(ymin= lowersd, ymax = uppersd), width = 0.5) + 
    labs(title = "T Cell In Vitro Infection" , y = "RB1B Genomes per 10^4", x = "time(days)") 
  

##################### 
## Genomes per 10^6 ## 
#################### 

#Cannot add the B + T cell plots because they are per B and T cells 
#Sabsabi dont know if trust 
bursa2024plot + spleen2024plot + thymus2024plot + PBL2024Plot + BCells2015Plot + TCells2015Plot

#maybe not use spleen because comes up earlier  


##################### 
## Plaque Asssays  ## 
####################

InVitro2007Plot + PBL2001NoVax + PBL2001WLNoVax + Spleen2001NoVax + Spleen2001WLNoVax 
  
  
Plaque_Assays <-  ggplot(data = PBL2001, aes(x = time, y = pfuxmil)) + 
   geom_point(data = PBL2001, aes(x = time, y = pfuxmil, color = "PBL")) + geom_line(data = PBL2001, aes(x = time, y = pfuxmil, color = "PBL")) +
   geom_point(data = PBL2001_WL, aes(x = time, y = pblxmil, color = "PBL_WL"), inherit.aes = FALSE) + 
   geom_line(data = PBL2001_WL, aes(x = time, y = pblxmil, color = "PBL_WL"), inherit.aes = FALSE) + 
   geom_point(data = Spleen2001, aes(x = time, y = PFUxmil, color = "Spleen"), inherit.aes = FALSE) + 
   geom_line(data = Spleen2001, aes(x = time, y = PFUxmil, color = "Spleen"), inherit.aes = FALSE) + 
   geom_point(data = Spleen2001_WL, aes(x = time_wk, y = pfuxmil, color = "Spleen_WL"), inherit.aes = FALSE) + 
   geom_line(data = Spleen2001_WL, aes(x = time_wk, y = pfuxmil, color = "Spleen_WL"), inherit.aes = FALSE) +  
   # geom_point(data = InVitro2007, aes(x = time, y = pfu, color = "CEF"), inherit.aes = FALSE) + 
   # geom_line(data = InVitro2007, aes(x = time, y = pfu, color = "CEF"), inherit.aes = FALSE) not including because sucks 
   labs(title = "Plaque Assays - WL and B15x7", y = "RB1B Plaques per 10^6 Plated Cells", x = "time(days)") + 
   scale_x_continuous(breaks = seq(0,28, by = 7)) + 
   scale_y_continuous(limits = c(0,1500), breaks = seq(0,1500, by = 250)) 
 
 
    
