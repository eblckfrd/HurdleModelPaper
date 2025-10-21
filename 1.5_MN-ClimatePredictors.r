# Dissertation Chapter 2 - MINNESOTA - September 2025
# Gathering and cleaning climate predictor data
# Note that large scale climate data is processed in jupyter hub notebooks

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
library(sf)

# Calculate Climate Index Predictors -----------------------------------------------------------------------------------------

#PAR, surface wind speed, air temp, precip, ENSO, PDO, NAM, geopotential height 500

#PDO index ----
PDO_index <- read.table("~/Desktop/Research/Data/Climate Index Data/PDO_index.txt", header = TRUE, sep = "", dec = ".")

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
AO_index <- read.table("~/Desktop/Research/Data/Climate Index Data/AO_index.txt", header = TRUE, sep = "", dec = ".")

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

Nino34_index <- read.table("~/Documents/R Repositories/Ch2_DataAnalysis/Nino3.4_Index.txt", header = TRUE, sep = "", dec = ".")

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

PNA_index <- read.table("~/Documents/R Repositories/Ch2_DataAnalysis/PNA_Index.txt", header = TRUE, sep = "", dec = ".")

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

NAM_index <- read.table("~/Documents/R Repositories/Ch2_DataAnalysis/NAM_Index.txt", header = TRUE, sep = "", dec = ".")

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


# Calculate Seasonal Global Climate Data ----------------------------------

#load in csv file-- ERA5 Reanalysis data downloaded and processed in Python Jupiter Notebook

MN_seasonalclimate <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/MN_seasonalclimate.csv")

#remove number column

MN_seasonalclimate <- MN_seasonalclimate %>%
  select(-number)

#now, we need to match each lake into the closest climate data grid box

# Get unique grid coordinates
grid_coords <- MN_seasonalclimate %>% 
  distinct(latitude, longitude)

# Function to find nearest grid cell
find_nearest_grid <- function(lake_lat, lake_lon, grid_coords) {
  distances <- sqrt((grid_coords$latitude - lake_lat)^2 + 
                      (grid_coords$longitude - lake_lon)^2)
  nearest_idx <- which.min(distances)
  return(data.frame(
    grid_lat = grid_coords$latitude[nearest_idx],
    grid_lon = grid_coords$longitude[nearest_idx],
    distance = distances[nearest_idx]
  ))
}

#apply to lakes
unique_lakeIDs <- unique_lakeIDs %>%
  rowwise() %>%
  mutate(grid_match = list(find_nearest_grid(latitude, longitude, grid_coords))) %>%
  unnest_wider(grid_match) %>%
  ungroup()

# Calculate grid cell resolution (spacing between points)
lat_res <- min(diff(sort(unique(grid_coords$latitude))))
lon_res <- min(diff(sort(unique(grid_coords$longitude))))

# Create rectangles for each grid cell
grid_cells <- grid_coords %>%
  mutate(
    xmin = longitude - lon_res/2,
    xmax = longitude + lon_res/2,
    ymin = latitude - lat_res/2,
    ymax = latitude + lat_res/2
  )

# Plot with grid cells as rectangles
ggplot() +
  geom_rect(data = grid_cells,
            aes(xmin = xmin, xmax = xmax, 
                ymin = ymin, ymax = ymax),
            fill = "lightblue", color = "blue", 
            alpha = 0.3, linewidth = 0.3) +
  # Lakes
  geom_point(data = unique_lakeIDs, 
             aes(x = longitude, y = latitude), 
             color = "navy", size = 2, alpha = 0.8) +
  # Lines connecting lakes to their matched grid cell centers
  geom_segment(data = unique_lakeIDs,
               aes(x = longitude, y = latitude,
                   xend = grid_lon, yend = grid_lat),
               alpha = 0.4, color = "red", linewidth = 0.5) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  coord_fixed()  # Keep aspect ratio correct for geographic data

#check summary stats and histogram
summary(unique_lakeIDs$distance)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00943 0.07857 0.10262 0.10381 0.12822 0.23981

hist(unique_lakeIDs$distance, 
     main = "Distance from Lakes to Nearest Grid Cell",
     xlab = "Distance (degrees)")

#pull out preseason predictors from MN_seasonalclimate
preseason_climate <- MN_seasonalclimate %>%
  select(latitude, longitude, year, 
         ttr_may, ttr_may_june, 
         tsr_may, tsr_may_june, 
         windspeed_may, windspeed_may_june)

#merge with lake ids
preseason_climate_lakes <- unique_lakeIDs %>%
  select(lagoslakeid, grid_lat, grid_lon) %>%
  left_join(preseason_climate, 
            by = c("grid_lat" = "latitude", "grid_lon" = "longitude"))

#add to pre-season predictors
preseason_predictors <- preseason_predictors %>%
  left_join(preseason_climate_lakes %>% select(-grid_lat, -grid_lon),
            by = c("lagoslakeid" = "lagoslakeid", "year" = "year"))

#now for concurrent climate predictors

concurrent_climate <- MN_seasonalclimate %>%
  select(latitude, longitude, year,
         ttr_august, ttr_jas,
         tsr_august, tsr_jas,
         windspeed_august, windspeed_jas)

# Merge with lake IDs
concurrent_climate_lakes <- unique_lakeIDs %>%
  select(lagoslakeid, grid_lat, grid_lon) %>%
  left_join(concurrent_climate,
            by = c("grid_lat" = "latitude", "grid_lon" = "longitude"))

# Add to concurrent_predictors
concurrent_predictors <- concurrent_predictors %>%
  left_join(concurrent_climate_lakes %>% select(-grid_lat, -grid_lon),
            by = c("lagoslakeid" = "lagoslakeid", "year" = "year"))

#any naans?
#see which columns they are mostly in: 
sapply(preseason_predictors, function(x) sum(is.nan(x)))


#final datasets at this point:

#unique_lakeIDs: all lagoslakeid, latitude, longitude, grid locations, and category separators
#MN_selected_profiles: all DO and temp profiles for lakes in unique_lakeIDs
#concurrent_predictors: all processed concurrent predictors for August, and JAS
#preseason_predictors: all processed pre-season predictors for June, and MJ



#precipitation and temperature from MDNR according to watershed unit ------------------------------







# Continue on in R script 02-MN_VHT.R then switch to jupyter hub notebooks
#sessionInfo()
#R version 4.3.1 (2023-06-16)
#Platform: x86_64-apple-darwin20 (64-bit)
#Running under: macOS 15.4.1



