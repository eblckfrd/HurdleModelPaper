# Diss Chapter 2 - started June 2025
# Processing of Minnesota data from Gretchen and Mike


# Load in packages

library(tidyverse)
library(lubridate)
library(zoo)
library(rLakeAnalyzer)
library(broom)
library(LAGOSNE)
library(ggplot2)
library(ggmap)

# Load in MN_profiles data from shared google drive

#MN_profiles <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/MN_profiles.csv")

# Load in LAGOS data
#lagosne_get(overwrite = T,dest_folder = LAGOSNE:::lagos_path())


##Profile data processing ------------------------------------------------------

# Only keep columns that are needed
# Filter for only lakes that stratify-- filter for geom ratio <4

MN_profiles <- MN_profiles %>% #takes data down to 714410 and 15 variables
  dplyr::select(link_lagosid, year, date, doy, latitude, longitude, 
                monitoringlocationidentifier, monitoringlocationname,
                profile_ident, n_dep_inprofile,
                depth_m, temp_c, do_ppm, lake_geomratio, lake_maxdepth_m, lake_area_m2,
                flag_do, flag_temp, flag_dep_diff) %>%
  filter(lake_geomratio <= "4") 

#look at data
glimpse(MN_profiles)

#look at only lake IDs and # of years sampled
lake_years_sampled <- MN_profiles %>%
  group_by(link_lagosid) %>%
  summarise(years_sampled = n_distinct(year), .groups = "drop")

#look at histogram
ggplot(data = lake_years_sampled, aes(x = years_sampled)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

# Filter for lakes with >10 years of data (should be ~500)
ten_yrlakes <- lake_years_sampled %>%
  filter(years_sampled >= 10) 

#look at histogram again
ggplot(data = ten_yrlakes, aes(x = years_sampled)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

#filter MN_profiles dataset for only lakes with >10 yrs of sampling data

MN_profiles_10yrs <- MN_profiles %>%
  filter(link_lagosid %in% ten_yrlakes$link_lagosid)

#removing Mille Lacs Lake and Rainy Lake for being excessive size outliers (see commented out code later below for validation)
MN_profiles_10yrs <- MN_profiles_10yrs %>%
  filter(!link_lagosid %in% c(2498, 73))

#add new month column
MN_profiles_10yrs <- MN_profiles_10yrs %>%
  mutate(date = as.Date(date)) %>%
  mutate(month = month(date))

# Histogram of what months data was collected in
ggplot(data = MN_profiles_10yrs, aes(x = month)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

# See how many lakes have profile_ident values that dont match the lake's max depth: 482 out of 517
#chose 2m mismatch as arbitrary wiggle room..

MN_profiles_10yrs %>%
  mutate(depth_mismatch = abs(n_dep_inprofile - lake_maxdepth_m) > 2) %>%
  # Step 2: Group by lake ID
  group_by(link_lagosid) %>%
  # Step 3: Summarize if any mismatch exists for the lake
  summarise(has_mismatch = any(depth_mismatch), .groups = "drop") %>%
  # Step 4: Count how many lakes have a mismatch
  summarise(num_lakes_with_mismatch = sum(has_mismatch))

  #should i remove these lakes? so that im only looking at lakes that have full depth profiles at the deep hole?
  #or do we want to consider lakes with sampling at multiple basins?
  #i think we should have just one sampling site per lake.. at the deep point.

# Look at how many lakes have more than one monitoring location identifier: 387 lakes out of 517

MN_profiles_10yrs %>%
  group_by(link_lagosid) %>%
  summarise(num_sites = n_distinct(monitoringlocationidentifier), .groups = "drop") %>%
  filter(num_sites > 1) %>%
  summarise(num_lakes_with_multiple_sites = n())

#see if these are the same lakes.... intersection of 370 lakes! So kinda?

depth_mismatch_lakes <- MN_profiles_10yrs %>%
  mutate(depth_mismatch = abs(n_dep_inprofile - lake_maxdepth_m) > 2) %>%
  group_by(link_lagosid) %>%
  summarise(has_mismatch = any(depth_mismatch), .groups = "drop") %>%
  filter(has_mismatch) %>%
  pull(link_lagosid)

multi_site_lakes <- MN_profiles_10yrs %>%
  group_by(link_lagosid) %>%
  summarise(num_sites = n_distinct(monitoringlocationidentifier), .groups = "drop") %>%
  filter(num_sites > 1) %>%
  pull(link_lagosid)

intersect_lakes <- intersect(depth_mismatch_lakes, multi_site_lakes)
length(intersect_lakes)

# How to get down to just one central site per lake? need it like this to be able to bin lakes based on categories...


# Filter for any specific months? Mar-Sep?


# Histogram of how many sampling events per site


#Filter for any # of events?


# Just take most recent surface sample for each lake in new table, to easily match with LAGOS data
#look at only lake IDs, lat, lon, lake_maxdepth_m, lake_area_m2
unique_lakeIDs <- MN_profiles_10yrs %>%
  group_by(link_lagosid) %>%
  summarise(
    latitude = first(latitude),
    longitude = first(longitude),
    lake_maxdepth_m = first(lake_maxdepth_m),
    lake_area_m2 = first(lake_area_m2),
    .groups = "drop"
  )

# Histogram of depths
ggplot(data = unique_lakeIDs, aes(x = lake_maxdepth_m)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

# Change lake area to ha from m2
unique_lakeIDs$lake_area_ha <- unique_lakeIDs$lake_area_m2 / 10000

# Histogram of lake area
ggplot(data = unique_lakeIDs, aes(x = lake_area_ha)) +
  geom_histogram(binwidth = 100, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

# unique_lakeIDs %>%
#   filter(lake_area_ha > 50000) #ID: 2498
# 
# unique_lakeIDs %>%
#   filter(lake_area_ha >= 15000 & lake_area_ha <= 20000) #ID: 73
# 
# MN_profiles %>%
#   filter (link_lagosid == 2498) #MILLE LACS LAKE, MN
# 
# MN_profiles %>%
#   filter (link_lagosid == 73) #RAINY LK, MN

#these two lakes are just sooooo huge compared to the other ones, so im removing them as outliers up at 'MN_profiles_10yrs'
#and commenting out this section


# Plot lake lat/lon onto map to see north/south spread

ggplot(unique_lakeIDs, aes(x = longitude, y = latitude)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(color = "blue", size = 2, alpha = 0.7) +
  coord_quickmap(
    xlim = c(-97.5, -89),   # Longitude range for Minnesota
    ylim = c(43.5, 49.5)    # Latitude range for Minnesota
  ) +
  theme_minimal() +
  labs(
    title = "Lake Locations in Minnesota (≥10 Years of Data)",
    x = "Longitude",
    y = "Latitude"
  )

#now colored in by depth
ggplot(unique_lakeIDs, aes(x = longitude, y = latitude, color = lake_maxdepth_m)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_gradient(
    name = "Max Depth (m)",
    low = "skyblue2",
    high = "darkblue"
  ) +
  theme_minimal() +
  labs(
    title = "Minnesota Lakes Colored by Max Depth",
    x = "Longitude",
    y = "Latitude"
  )

##North and South investigation -------------------------------------------------

# Perform cluster analysis on rainfall and temp Mesonet data to define north and south regions?


##LAGOS merge -------------------------------------------------------------


# Filter for MN

# Merge with profile data using LAGOS_ID column

dataset$column <-merge(dataset,dataset2,by=c("lagoslakeid","year"),all.x=TRUE)

# First I need data so that i can categorize my lakes..
# Filter columns for: ID, most recent secchi depth, lake connectivity class, elevation, max lake depth

# Look at histogram of depths, then pick thresholds
  #

# From secchi depth, calculate TSI for each lake: 60 - (14.41 * ln(secci #))
# Look at histogram of TSI to see spread
# Then assign oligotrophic as <30, and eutrophic as >50 if this makes sense


# Assign new columns to table: depth (assign deep or shallow), productivity (oligotrophic or eutrophic), hydrology (seepage or drainage)



# Re-merge data with full profiles

