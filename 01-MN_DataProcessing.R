# Dissertation Chapter 2 - started June 2025
# Processing of Minnesota data from MGLP team (Mike, Gretchen, Heidi, Leslie)


# Load in packages

library(tidyverse)
library(lubridate)
library(zoo)
library(rLakeAnalyzer)
library(broom)
library(LAGOSNE)
library(ggplot2)
library(ggmap)



# Load in MN_profiles data from Gretchen shared google drive

#FULL_MN_profiles <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/MN_profiles.csv")

# Load in LAGOS data
#lagosne_get(overwrite = T,dest_folder = LAGOSNE:::lagos_path())


##Profile data processing ------------------------------------------------------


# Only keep columns that are needed
# Filter for only lakes that stratify-- filter for geom ratio <4

MN_profiles <- FULL_MN_profiles %>% #takes data down to 714,410 and 15 variables
  dplyr::select(link_lagosid, year, date, doy, latitude, longitude, 
                monitoringlocationidentifier, monitoringlocationname,
                profile_ident, n_dep_inprofile,
                depth_m, temp_c, do_ppm, lake_geomratio, lake_maxdepth_m, lake_area_m2,
                flag_do, flag_temp, flag_dep_diff) %>%
  filter(lake_geomratio <= "4") 

#just so i can run the borrowed code below
MN_profiles <- MN_profiles %>%
  rename(MonitoringLocationIdentifier = monitoringlocationidentifier)

# Remove Lake Superior Sites, Mines, and Pits ~ We only want Glacial Lakes; from MGLP
MN_profiles = MN_profiles %>%
  filter(MonitoringLocationIdentifier != "EPA_GLNPO-SU19") %>% #Lake Superior Site
  filter(MonitoringLocationIdentifier != "EPA_GLNPO-SU18") %>% #Lake Superior Site
  filter(MonitoringLocationIdentifier != "EPA_GLNPO-SU104550")%>% #Lake Superior Site
  filter(MonitoringLocationIdentifier != "21GPBCH-LS_PT_1") %>% #Lake Superior Site
  filter(MonitoringLocationIdentifier != "21GPBCH-LS_PT_9") %>% #Lake Superior Site
  filter(MonitoringLocationIdentifier != "GPC5_WQX-LS_PT_1") %>% #Lake Superior Site
  filter(MonitoringLocationIdentifier != "GPC5_WQX-LS_PT_9") %>% #Lake Superior Site
  filter(MonitoringLocationIdentifier != "GPC5_WQX-LS_PT_8") %>% #Lake Superior Site
  filter(MonitoringLocationIdentifier != "MNPCA-69-1297-00-202") %>% #Sherman Pit
  filter(MonitoringLocationIdentifier != "MNPCA-69-0429-00-201") %>% #Sabin Pit
  filter(MonitoringLocationIdentifier != "21GPBCH-LS_PT_8") %>%      #Lake Superior Site
  filter(MonitoringLocationIdentifier != "MNPCA-18-0440-01-100") %>% #Mahnomen Mine
  filter(MonitoringLocationIdentifier != "MNPCA-69-1297-00-201") %>% #Sherman Pit
  filter(MonitoringLocationIdentifier != "MNPCA-18-0439-00-100") %>% #Pennington Mine
  filter(MonitoringLocationIdentifier != "MNPCA-31-1325-04-201") %>% #Canisteo Mine Pit
  filter(MonitoringLocationIdentifier != "MNPCA-69-1298-00-203") %>% #South Twin City Pit
  filter(MonitoringLocationIdentifier != "MNPCA-18-0093-01-202") %>% #Rabbit Lake
  filter(MonitoringLocationIdentifier != "MNPCA-18-0041-01-201") %>% #Crooked Lake - Silver Bay
  filter(MonitoringLocationIdentifier != "MNPCA-69-0429-00-100") %>% #Embarrass Mine
  filter(MonitoringLocationIdentifier != "MNPCA-69-0428-00-201") #St. James Pit

  #now switch back name
MN_profiles <- MN_profiles %>%
    rename(monitoringlocationidentifier = MonitoringLocationIdentifier)

#filter for min depth: code from MGLP team
valid_profiles <- MN_profiles %>% #Profiles must have a minimum depth less than 3 meters; 
  group_by(profile_ident) %>% #otherwise the profiles shallowest read could be in 
  summarise(min.depth = min(depth_m, na.rm = TRUE)) %>% #the metalimnion and capture none of the epilimnion
  filter(min.depth < 3) %>% #Jane et al 2021
  pull(profile_ident)

MN_profiles <- MN_profiles %>%
  filter(profile_ident %in% valid_profiles)

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

# Histogram of # of depths in profile
ggplot(data = MN_profiles_10yrs, aes(x = n_dep_inprofile)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()


# See how many lakes have n_dep_inprofile values that dont match the lake's max depth: 482 out of 517
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

# from MGLP: we can filter for most sampled site at each lake and use that
mostsampledsite <- MN_profiles_10yrs %>%
  group_by(link_lagosid, monitoringlocationidentifier) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(link_lagosid) %>%
  mutate(MostSampled = max(n)) %>%
  ungroup() %>%
  filter(n == MostSampled) %>%
  distinct(link_lagosid, .keep_all = TRUE)

MN_profiles_10yrs <- MN_profiles_10yrs %>%
  filter(monitoringlocationidentifier %in% mostsampledsite$monitoringlocationidentifier)

#check that this worked (it did!)
MN_profiles_10yrs %>%
  group_by(link_lagosid) %>%
  summarise(n_sites = n_distinct(monitoringlocationidentifier)) %>%
  filter(n_sites > 1)

#another check
MN_profiles_10yrs %>%
  count(link_lagosid, monitoringlocationidentifier) %>%
  arrange(desc(n)) %>%
  head(10)

#now look at this number again; 464
MN_profiles_10yrs %>%
  mutate(depth_mismatch = abs(n_dep_inprofile - lake_maxdepth_m) > 2) %>%
  # Step 2: Group by lake ID
  group_by(link_lagosid) %>%
  # Step 3: Summarize if any mismatch exists for the lake
  summarise(has_mismatch = any(depth_mismatch), .groups = "drop") %>%
  # Step 4: Count how many lakes have a mismatch
  summarise(num_lakes_with_mismatch = sum(has_mismatch))

# Histogram of n of samples
ggplot(data = mostsampledsite, aes(x = n)) +
  geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

# Filter for any specific months? Mar-Oct, 405700 obs left

MN_profiles_10yrs <- MN_profiles_10yrs %>%
  filter(month >= 3 & month <= 10)

# Identify lakes with at least one sample during June–September (DOY 152–273)
#they all fit this criteria so dont need to filter!
lakes_with_summer_sample <- MN_profiles_10yrs %>%
  filter(doy >= 152 & doy <= 273) %>%
  distinct(link_lagosid)


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

