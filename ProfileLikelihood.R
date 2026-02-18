### Creating Pairwise Plots for Best Parameters ### 

#rm(list = ls())
library(deSolve)  
library(tidyr) 
library(ggplot2)
library(dplyr) 
library(readxl)  
library(purrr)  
library(GGally)

optim_data <- read.csv("/Users/rayanhg/Desktop/WithinHostModel/CodeOutputsRandNum/Feb.17.26.FittingDnbinom_mu_theta_ALL.csv") %>% 
  filter(Converged == 0, Likely <= 316.5815 ) %>% arrange(Likely) %>%  
  select(c(beta, beta_2, alpha, alpha_2,nu_a,nu_b,nu_f,mu,g1,g2,h1,h2,Pb)) %>% 
  mutate(across(c(beta, beta_2, alpha, alpha_2, nu_a, nu_b, nu_f, mu),log10, .names = "log_{.col}")) %>% 
  select(-c(beta, beta_2,alpha, alpha_2,nu_a,nu_b,nu_f,mu))



ggpairs(optim_data)

