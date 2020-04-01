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

# Load in .csv data from terminal suppression
data = read.csv("2x2x3 data.csv")

###################################################################################################

## GOAL: 2x2x3 anova, 1 factor is Light condition (2 levels - within subject), 
# 1 factor is Virus type (2 levels - bw subject),
# 1 factor is task phase (3 levels of task phase - within subject)

# split data according to Injected region
data_PFC = data %>% filter(Site != "Sub")
data_SUB = data %>% filter(Site != "Pfc")

# ez anova
library(ez) # load library

# run anova using ezANOVA
PFC_2x2x3 = ezANOVA(
  data       = data_PFC
  , dv       = Accuracy           # dv
  , wid      = RatID              # subjects
  , within   = .(Light,Condition) # within subjects factors - condition is task phase
  , between  = Virus
  , detailed = TRUE
  , type     = 3 
)
print(PFC_2x2x3)

SUB_2x2x3 = ezANOVA(
  data       = data_SUB
  , dv       = Accuracy           # dv
  , wid      = RatID              # subjects
  , within   = .(Light,Condition) # within subjects factors - condition is task phase
  , between  = Virus
  , detailed = TRUE
  , type     = 3 
)
print(SUB_2x2x3)

