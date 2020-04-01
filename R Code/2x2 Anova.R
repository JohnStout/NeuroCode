# Analyze behavioral data from terminal suppression paper
#
# To run line-by-line, highlight code, hit cntrl+enter
#
# This script requires an excel sheet to be saved as a .csv in long format

# load in tidyverse
library(tidyverse)

# set working directory - note that this needs to change. It is the directory
# where your data is stored.
setwd("X:/07. Manuscripts/In preparation/Terminal Suppression")

# Load in .csv data formatted long-ways
data = read.csv("Both Blocks 2x2 data_revised.csv")

####################################################################################################

###### 2(group) x 2(light) anova and consider each task phase as separate experiments

# split data according to Injected region
data_PFC = data %>% filter(Site != "Sub")
data_SUB = data %>% filter(Site != "Pfc")

# PFC->Re isolate each condition
data_PFC_sample = data_PFC %>% filter(Condition == "Sample")
data_PFC_delay  = data_PFC %>% filter(Condition == "Delay")
data_PFC_choice = data_PFC %>% filter(Condition == "Choice")

# HPC->Re 
data_SUB_sample = data_SUB %>% filter(Condition == "Sample")
data_SUB_delay  = data_SUB %>% filter(Condition == "Delay")
data_SUB_choice = data_SUB %>% filter(Condition == "Choice")

# ez anova
library(ez) # load library

## run anova pfc 
PFC_2x2_sample = ezANOVA(
  data       = data_PFC_sample
  , dv       = Accuracy           # dv
  , wid      = RatID              # subjects
  , within   = .(Light) # within subjects factors - condition is task phase
  , between  = Virus
  , detailed = TRUE
  , type     = 3 
)
print(PFC_2x2_sample)

PFC_2x2_delay = ezANOVA(
  data       = data_PFC_delay
  , dv       = Accuracy           # dv
  , wid      = RatID              # subjects
  , within   = .(Light) # within subjects factors - condition is task phase
  , between  = Virus
  , detailed = TRUE
  , type     = 3 
)
print(PFC_2x2_delay)

PFC_2x2_choice = ezANOVA(
  data       = data_PFC_choice
  , dv       = Accuracy           # dv
  , wid      = RatID              # subjects
  , within   = .(Light) # within subjects factors - condition is task phase
  , between  = Virus
  , detailed = TRUE
  , type     = 3 
)
print(PFC_2x2_choice)

## run anova hpc
SUB_2x2_sample = ezANOVA(
  data       = data_SUB_sample
  , dv       = Accuracy           # dv
  , wid      = RatID              # subjects
  , within   = .(Light) # within subjects factors - condition is task phase
  , between  = Virus
  , detailed = TRUE
  , type     = 3 
)
print(SUB_2x2_sample)

SUB_2x2_delay = ezANOVA(
  data       = data_SUB_delay
  , dv       = Accuracy           # dv
  , wid      = RatID              # subjects
  , within   = .(Light) # within subjects factors - condition is task phase
  , between  = Virus
  , detailed = TRUE
  , type     = 3 
)
print(SUB_2x2_delay)

SUB_2x2_choice = ezANOVA(
  data       = data_SUB_choice
  , dv       = Accuracy           # dv
  , wid      = RatID              # subjects
  , within   = .(Light) # within subjects factors - condition is task phase
  , between  = Virus
  , detailed = TRUE
  , type     = 3 
)
print(SUB_2x2_choice)

##################################################################################

# ~~~ pairwise t-tests ~~~ #

# for light on/off
#PFC_pwtt = pairwise.t.test(data_PFC_choice$Accuracy, ), 
#paired=TRUE, p.adjust.method="holm")
#print(PFC_pwtt)

