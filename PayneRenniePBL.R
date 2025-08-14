rm(list = ls())
library(deSolve)  
library(reshape2) 
library(ggplot2)
library(dplyr) 
library(patchwork) 
library(writexl) 
library(readxl)




Payne<- read_excel(path = "/Users/rayanhg/Downloads/payneRennie1976.xlsx", col_types = c("numeric", "numeric", "numeric","numeric", "numeric"), na = c("", "NA"))  

Payne_mean <- Payne %>% mutate(logInfBcell = log10(infectBCells), logCtrlB = log10(BCellsControl)) %>% group_by(days) %>% 
  summarise(meanInfBcell = mean(logInfBcell), antilog = 10^meanInfBcell, meanCtrlB = mean(logCtrlB), antilogCtrl = 10^meanCtrlB ) %>%  
  filter(!is.na(antilog))  

PayneMeanLongDf <- melt(Payne_mean, id.vars = c("days", "meanInfBcell", "meanCtrlB")) %>% select(!c("meanInfBcell", "meanCtrlB"))

PayneLongB <- melt(Payne, id.vars = c("days")) %>% filter(!is.na(value), variable %in% c("infectBCells", "BCellsControl")) 


ggplot(data = PayneLongB) +geom_point(mapping =aes(x = days, y= value,group = variable, colour = variable), position = position_jitter(width = 0.55, height =0), size = 1.5) + 
  scale_y_log10(limits = c(100, 20000), breaks = c(100, 1000, 10000, 20000), minor_breaks = c(200, 300, 400, 500, 600, 700, 800, 900,2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000)) + 
  scale_x_continuous(limits = c(0,45), breaks = seq(0,45,by=5)) + geom_line(data = PayneMeanLongDf, mapping = aes(x = days, y = value, group = variable, colour = variable, linetype = variable)) + annotation_logticks(base = 10, sides = "l") +
  scale_linetype_manual(values = c("antilog" = "dashed", "antilogCtrl"= "solid")) + scale_color_manual(values = c("antilog"="indianred2", "antilogCtrl"="mediumseagreen", "infectBCells"= "indianred2", "BCellsControl"= "mediumseagreen"))+ theme()
