# Dissertation Chapter 2 - MINNESOTA - September 2026
# Gathering and cleaning all climate predictor data

# Load in packages

library(tidyverse)
library(lubridate)
library(zoo)
library(rLakeAnalyzer)
library(broom)
library(LAGOSNE)
library(ggplot2)
library(ggmap)
library(gridExtra)
library(mgcv)

# Calculate CLIMATE PREDICTORS -----------------------------------------------------------------------------------------

#PAR, surface wind speed, air temp, precip, ENSO, PDO, NAM, geopotential height 500

#PDO index ----
#PDO_index <- read.table("~/Desktop/Research/Data/Climate Index Data/PDO_index.txt", header = TRUE, sep = "", dec = ".")

#may-june average
PDO_index <- PDO_index %>%
  mutate(MJ_Avg = rowMeans(select(., May, Jun), na.rm = TRUE))

#JAS
PDO_index <- PDO_index %>%
  mutate(JAS_Avg = rowMeans(select(., Jul, Aug, Sep), na.rm = TRUE))

#add all to full predictor datasets
preseason_predictors <- preseason_waterq_predictors %>%
  left_join(PDO_index %>% select(Year, Jun, MJ_Avg), 
            by = c("year" = "Year")) %>%
  rename(PDO_June = Jun, PDO_MJ = MJ_Avg)

concurrent_predictors <- concurrent_waterq_predictors %>%
  left_join(PDO_index %>% select(Year, Jul, JAS_Avg), 
            by = c("year" = "Year")) %>%
  rename(PDO_Jul = Jul, PDO_JAS = JAS_Avg)

#NAM/AO index -------
#AO_index <- read.table("~/Desktop/Research/Data/Climate Index Data/AO_index.txt", header = TRUE, sep = "", dec = ".")

#may-june average
AO_index <- AO_index %>%
  mutate(MJ_Avg = rowMeans(select(., May, Jun), na.rm = TRUE))

#JAS
AO_index <- AO_index %>%
  mutate(JAS_Avg = rowMeans(select(., Jul, Aug, Sep), na.rm = TRUE))

#add all to full predictor datasets
preseason_predictors <- preseason_predictors %>%
  left_join(AO_index %>% select(Year, Jun, MJ_Avg), 
            by = c("year" = "Year")) %>%
  rename(AO_June = Jun, AO_MJ = MJ_Avg)

concurrent_predictors <- concurrent_predictors %>%
  left_join(AO_index %>% select(Year, Jul, JAS_Avg), 
            by = c("year" = "Year")) %>%
  rename(AO_Jul = Jul, AO_JAS = JAS_Avg)


#Nino 3.4 index -----

#Nino34_index <- read.table("~/Documents/R Repositories/Ch2_DataAnalysis/Nino3.4_Index.txt", header = TRUE, sep = "", dec = ".")

#may-june average
Nino34_index <- Nino34_index %>%
  mutate(MJ_Avg = rowMeans(select(., May, Jun), na.rm = TRUE))

#JAS
Nino34_index <- Nino34_index %>%
  mutate(JAS_Avg = rowMeans(select(., Jul, Aug, Sep), na.rm = TRUE))

#add all to full predictor datasets
preseason_predictors <- preseason_predictors %>%
  left_join(Nino34_index %>% select(Year, Jun, MJ_Avg), 
            by = c("year" = "Year")) %>%
  rename(ENSO_June = Jun, ENSO_MJ = MJ_Avg)

concurrent_predictors <- concurrent_predictors %>%
  left_join(Nino34_index %>% select(Year, Jul, JAS_Avg), 
            by = c("year" = "Year")) %>%
  rename(ENSO_Jul = Jul, ENSO_JAS = JAS_Avg)

#PNA index -------

#PNA_index <- read.table("~/Documents/R Repositories/Ch2_DataAnalysis/PNA_Index.txt", header = TRUE, sep = "", dec = ".")

#may-june average
PNA_index <- PNA_index %>%
  mutate(MJ_Avg = rowMeans(select(., May, Jun), na.rm = TRUE))

#JAS
PNA_index <- PNA_index %>%
  mutate(JAS_Avg = rowMeans(select(., Jul, Aug, Sep), na.rm = TRUE))

#add all to full predictor datasets
preseason_predictors <- preseason_predictors %>%
  left_join(PNA_index %>% select(Year, Jun, MJ_Avg), 
            by = c("year" = "Year")) %>%
  rename(PNA_June = Jun, PNA_MJ = MJ_Avg)

concurrent_predictors <- concurrent_predictors %>%
  left_join(PNA_index %>% select(Year, Jul, JAS_Avg), 
            by = c("year" = "Year")) %>%
  rename(PNA_Jul = Jul, PNA_JAS = JAS_Avg)

#NAM index ----

#NAM_index <- read.table("~/Documents/R Repositories/Ch2_DataAnalysis/NAM_Index.txt", header = TRUE, sep = "", dec = ".")

#may-june average
NAM_index <- NAM_index %>%
  mutate(MJ_Avg = rowMeans(select(., May, Jun), na.rm = TRUE))

#JAS
NAM_index <- NAM_index %>%
  mutate(JAS_Avg = rowMeans(select(., Jul, Aug, Sep), na.rm = TRUE))

#add all to full predictor datasets
preseason_predictors <- preseason_predictors %>%
  left_join(NAM_index %>% select(Year, Jun, MJ_Avg), 
            by = c("year" = "Year")) %>%
  rename(NAM_June = Jun, NAM_MJ = MJ_Avg)

concurrent_predictors <- concurrent_predictors %>%
  left_join(NAM_index %>% select(Year, Jul, JAS_Avg), 
            by = c("year" = "Year")) %>%
  rename(NAM_Jul = Jul, NAM_JAS = JAS_Avg)

#air temp and precip ----------

#MN_precip


#PAR, surface wind speed, geopotential height 500mb ------------

#this is sourced and processed in python worksheet " ", then exported and uploaded here in correct formats


