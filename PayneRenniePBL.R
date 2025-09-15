rm(list = ls())
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(dplyr) 
library(patchwork) 
library(writexl) 
library(readxl) 
library(scales)




Payne<- read_excel(path = "/Users/rayanhg/Downloads/payneRennie1976.xlsx", col_types = c("numeric", "numeric", "numeric","numeric", "numeric"), na = c("", "NA"))  

Payne_mean <- Payne %>% mutate( 
  logInfBcell = log10(infectBCells),
  logCtrlB = log10(BCellsControl), 
  logInfTcell = log10(InfectTcells),
  logCtrlTcell = log10(TCellsControl)) %>% 
  group_by(days) %>%
  summarise(
    meanInfBcell = mean(logInfBcell, na.rm = TRUE), 
    meanCtrlB = mean(logCtrlB), 
    meanInfTcell = mean(logInfTcell, na.rm = TRUE), 
    meanCtrlT = mean(logCtrlTcell), 
    antilog = 10^meanInfBcell, 
    antilogCtrl = 10^meanCtrlB, 
    antilogT = 10^meanInfTcell, 
    antilogCtrlT = 10^meanCtrlT) 
  
  
  

PayneMeanLongDf <- melt(Payne_mean, id.vars = c("days")) %>%
  filter(variable %in% c("antilog", "antilogCtrl"))

PayneMeanLongDfT <- melt(Payne_mean, id.vars = c("days")) %>% 
  filter( variable %in% c("antilogT", "antilogCtrlT")) 



PayneLongB <- melt(Payne, id.vars = c("days")) %>% filter(!is.na(value), variable %in% c("infectBCells", "BCellsControl"))
PayneLongT <- melt(Payne, id.vars = c("days")) %>% filter(!is.na(value), variable %in% c("InfectTcells", "TCellsControl"))  




p1 <- ggplot(data = PayneLongB) +geom_point(mapping =aes(x = days, y= value,group = variable, colour = variable), position = position_jitter(width = 0.55, height =0), size = 1.5) + 
  scale_y_log10(limits = c(100, 20000), breaks = c(100, 1000, 10000, 20000), minor_breaks = c(200, 300, 400, 500, 600, 700, 800, 900,2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000)) + 
  scale_x_continuous(limits = c(0,45), breaks = seq(0,45,by=5)) + geom_line(data = PayneMeanLongDf, mapping = aes(x = days, y = value, group = variable, colour = variable, linetype = variable), size = 1.2) + annotation_logticks(base = 10, sides = "l") +
  scale_linetype_manual(values = c("antilog" = "dashed", "antilogCtrl"= "solid")) + scale_color_manual(values = c("antilog"="blue", "antilogCtrl"="blue", "infectBCells"= "blue", "BCellsControl"= "blue"))+ 
  theme(panel.grid = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.line = element_line(color = "black"))


p2 <- ggplot(data = PayneLongT) +geom_point(mapping =aes(x = days, y= value,group = variable, colour = variable), position = position_jitter(width = 0.55, height =0), size = 1.5) + 
  scale_y_log10(limits = c(1000, 100000), labels = label_number(), breaks = c(1000, 10000, 100000), minor_breaks = c(2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000 ,20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000)) + 
  scale_x_continuous(limits = c(0,45), breaks = seq(0,45,by=5)) + geom_line(data = PayneMeanLongDfT, mapping = aes(x = days, y = value, group = variable, colour = variable, linetype = variable), size = 1.2) + annotation_logticks(base = 10, sides = "l") +
  scale_linetype_manual(values = c("antilogT" = "dashed", "antilogCtrlT"= "solid")) + scale_color_manual(values = c("antilogT"="hotpink", "antilogCtrlT"="hotpink", "InfectTcells"= "hotpink", "TCellsControl"= "hotpink"))+ 
theme(panel.grid = element_blank(), panel.background = element_blank(),  legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.line = element_line(color = "black"))

p1/p2

