# Dissertation Chapter 2 - MINNESOTA - started June 2025
# Processing of Minnesota data from MGLP team (Mike, Gretchen Hansen, Heidi Rantala, Lesley Knoll)
# Processing LAGOS data
# Gathering and cleaning all other predictor data

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

#----------------------------------------------------------------------------------------------------------

# Load in MN_profiles data from Gretchen's shared google drive

#FULL_MN_profiles <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/MN_profiles.csv")

#-----------------------

# Load in LAGOS data files downloaded from EDI
#was originally using only LAGOS-NE data but it is less up to date than LAGOS-US
#so for water quality I can pull full country then filter to just my selected lakes, and later merge NE with US pre 2012 and post 2012

#LAGOSNE_state_105 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/LAGOSNE_state_105.csv")
        #state IDs
#LAGOSNE_lakeslocus101 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/LAGOSNE_lakeslocus101.csv")
         #lagos ID, name, lake area, lake perim, hu zones, elevation, stat and county zone id
#lake_characteristics <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/lake_characteristics.csv")
        #hydrology/lake connectivity

#this is from LAGOS-US, all water quality parameters
  #chemistry_limno <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/chemistry_limno.csv")

#LAGOS-NE water quality params
  #epi_waterquality10873 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/epi_waterquality10873.csv")
  

##Profile data processing ------------------------------------------------------------------------------------

# Only keep columns that are needed
# Filter for only lakes that stratify-- filter for geom ratio <4 (ratio already calculated in profiles dataset)

MN_profiles <- FULL_MN_profiles %>% #takes data down to 714,410 and 15 variables
  dplyr::select(link_lagosid, year, date, doy, latitude, longitude, 
                monitoringlocationidentifier, monitoringlocationname,
                profile_ident, n_dep_inprofile,
                depth_m, temp_c, do_ppm, lake_geomratio, lake_maxdepth_m, lake_area_m2,
                flag_do, flag_temp, flag_dep_diff) %>%
  filter(lake_geomratio <= "4") 

#rename just so i can run the borrowed code from Mike below
MN_profiles <- MN_profiles %>%
  rename(MonitoringLocationIdentifier = monitoringlocationidentifier)

# Remove Lake Superior Sites, Mines, and Pits ~ We only want Glacial Lakes; from Mike code
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

#filter for min depth: code from Mike, based on Jane et al 2021
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

#removing Mille Lacs Lake and Rainy Lake for being excessive size outliers (see commented out code later below for reasoning)
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

# Histogram of years included = most recent is 2020, goes back to ~1950 but mostly >1980
ggplot(data = MN_profiles_10yrs, aes(x = year)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

# trim data to only include 1980-2020, to standardize years of predictors
MN_profiles_10yrs <- MN_profiles_10yrs %>%
  filter(year >= 1980 & year <= 2020)


# See how many lakes have n_dep_inprofile values that dont match the lake's max depth: 406 
#chose 5m mismatch as arbitrary wiggle room..

MN_profiles_10yrs %>%
  mutate(depth_mismatch = abs(n_dep_inprofile - lake_maxdepth_m) > 5) %>%
  # Step 2: Group by lake ID
  group_by(link_lagosid) %>%
  # Step 3: Summarize if any mismatch exists for the lake
  summarise(has_mismatch = any(depth_mismatch), .groups = "drop") %>%
  # Step 4: Count how many lakes have a mismatch
  summarise(num_lakes_with_mismatch = sum(has_mismatch))

  #should i remove these lakes? so that im only looking at lakes that have full depth profiles at the deep hole?
  #or do we want to consider lakes with sampling at multiple basins?
  #i think we should have just one sampling site per lake..

# Look at how many lakes have more than one monitoring location identifier: 387 lakes out of 517

MN_profiles_10yrs %>%
  group_by(link_lagosid) %>%
  summarise(num_sites = n_distinct(monitoringlocationidentifier), .groups = "drop") %>%
  filter(num_sites > 1) %>%
  summarise(num_lakes_with_multiple_sites = n())

#see if these are the same lakes.... intersection of 364 lakes! So kinda?

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

# from Mike code: we can filter for most sampled site at each lake and use that; preserves the most amount of data
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

# Look at the distribution of depth mismatches
depth_distribution <- MN_profiles_10yrs %>%
  mutate(depth_diff = abs(n_dep_inprofile - lake_maxdepth_m)) %>%
  summarise(
    mean_diff = mean(depth_diff, na.rm = TRUE),
    median_diff = median(depth_diff, na.rm = TRUE),
    q75_diff = quantile(depth_diff, 0.75, na.rm = TRUE),
    q90_diff = quantile(depth_diff, 0.9, na.rm = TRUE),
    q95_diff = quantile(depth_diff, 0.95, na.rm = TRUE),
    max_diff = max(depth_diff, na.rm = TRUE)
  )

print(depth_distribution)

#mean_diff median_diff q75_diff  q90_diff  q95_diff  max_diff
#4.113012    1.811033  5.291329  10.61911  15.49476  62.35855

#First identified 72 lakes that dont have >70% coverage in any profiles
#FILTER AND REMOVE any profiles where the depth has >50% difference from lake max depth

problem_lakes <- MN_profiles_10yrs %>%
  mutate(
    relative_diff = abs(n_dep_inprofile - lake_maxdepth_m) / lake_maxdepth_m,
    depth_mismatch = relative_diff > 0.5
  ) %>%
  group_by(link_lagosid) %>%
  summarise(
    all_profiles_flagged = all(depth_mismatch),
    .groups = "drop"
  ) %>%
  filter(all_profiles_flagged) %>%
  pull(link_lagosid)

print(paste("Number of problem lakes to remove:", length(problem_lakes)))

# Remove the 72 problem lakes entirely
MN_profiles_final <- MN_profiles_10yrs %>%
  filter(!link_lagosid %in% problem_lakes) %>%
  # Apply the 50% relative threshold to remaining lakes
  mutate(
    relative_diff = abs(n_dep_inprofile - lake_maxdepth_m) / lake_maxdepth_m,
    depth_mismatch = relative_diff > 0.5
  ) %>%
  filter(!depth_mismatch) %>%
  select(-relative_diff, -depth_mismatch)

# Summary of final dataset
final_summary <- tibble(
  approach = c("Original", "Final (removed 72 lakes + 50% filter)"),
  lakes = c(
    n_distinct(MN_profiles_10yrs$link_lagosid),
    n_distinct(MN_profiles_final$link_lagosid)
  ),
  profiles = c(
    nrow(MN_profiles_10yrs),
    nrow(MN_profiles_final)
  )
) %>%
  mutate(
    lakes_removed = c(0, lakes[1] - lakes[2]),
    profiles_removed = c(0, profiles[1] - profiles[2]),
    percent_lakes_kept = round(lakes / lakes[1] * 100, 1),
    percent_profiles_kept = round(profiles / profiles[1] * 100, 1)
  )

print("\nFinal dataset summary:")
print(final_summary)


# Filter for any specific months? Mar-Oct, 371412 obs left

MN_profiles_final <- MN_profiles_final %>%
  filter(month >= 3 & month <= 10)

# Identify lakes with at least one sample during July–September (DOY 182–273)
#8 lakes dont fit this requirement.. only 421 lakes
lakes_with_summer_sample <- MN_profiles_final %>%
  filter(doy >= 182 & doy <= 273) %>%
  distinct(link_lagosid)

#do the same with spring sample April-June (91-181)
#only 416 lakes..
lakes_with_spring_sample <- MN_profiles_final %>%
  filter(doy >= 91 & doy <= 181) %>%
  distinct(link_lagosid)

#find which lakes have at least one spring and one summer sample
lakes_with_both <- intersect(
  lakes_with_spring_sample$link_lagosid,
  lakes_with_summer_sample$link_lagosid
)

# Confirm how many lakes remain: 410
length(lakes_with_both)

#view lakes
excluded_lakes <- setdiff(unique(MN_profiles_final$link_lagosid), lakes_with_both)

#filter dataset to remove those lakes
MN_profiles_final <- MN_profiles_final %>%
  filter(link_lagosid %in% lakes_with_both)

# Just take most recent surface sample for each lake in new table, to easily match with LAGOS data
#look at only lake IDs, lat, lon, lake_maxdepth_m, lake_area_m2
unique_lakeIDs <- MN_profiles_final %>%
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

## from Mike code, this selects the latest profile of the year; may be of use later?

EndofYear = MN_profiles_final %>%
  group_by(link_lagosid, year)%>%
  mutate(Max = max(doy))%>%
  filter(doy == Max)%>%
  ungroup()

#finally... create new profiles data with only the lake ideas narrowed down in the unique_lakeIDs dataset

MN_selected_profiles <- MN_profiles_final %>%
  filter(link_lagosid %in% unique_lakeIDs$link_lagosid)

#check that this worked; it did! 410 lakes included
n_distinct(MN_selected_profiles$link_lagosid)


####ignore this, ended up not being an issue with random forest models, they allow for non-consecutive samples as long as predictor years match

# #do I need only lakes that have consecutive sampling?
# lake_years <- MN_selected_profiles %>%
#   distinct(link_lagosid, year) %>%
#   arrange(link_lagosid, year)
# 
# #sort
# lakes_with_consec_10ish <- lake_years %>%
#   group_by(link_lagosid) %>%
#   summarise(sample_years = list(sort(unique(year))), .groups = "drop") %>%
#   rowwise() %>%
#   mutate(has_consec_block = {
#     years <- sample_years
#     found <- FALSE
#     if (length(years) >= 10) {
#       for (i in 1:(length(years) - 9)) {
#         window <- years[i:(i + 9)]
#         if (max(window) - min(window) <= 13) {  # Allow up to 1–2 year gaps
#           found <- TRUE
#           break
#         }
#       }
#     }
#     found
#   }) %>%
#   filter(has_consec_block) %>%
#   select(link_lagosid)
# 
# #see all streaks
# lake_streaks <- lake_years %>%
#   group_by(link_lagosid) %>%
#   mutate(streak_id = year - row_number()) %>%
#   group_by(link_lagosid, streak_id) %>%
#   summarise(
#     streak_length = n(),
#     min_year = min(year),
#     max_year = max(year),
#     .groups = "drop"
#   )
# 
# #identify the longest streak for each year
# longest_streaks <- lake_streaks %>%
#   group_by(link_lagosid) %>%
#   slice_max(order_by = streak_length, n = 1, with_ties = FALSE) %>%
#   ungroup()
# 
# #histogram
# ggplot(data = longest_streaks, aes(x = streak_length)) +
#   geom_histogram(binwidth = 2, fill = "blue", color = "black", alpha = 0.7) +
#   theme_minimal()


# Analysis on rainfall and temp Mesonet data to define north and south regions: SUMMER AVERAGE JJA ----------------------------------------------------
# To confirm ecoregion seperation

# NCEI 1991-2020 Climatological Average for Summer (JJA)
#https://mesonet.agron.iastate.edu/plotting/auto/?_wait=no&q=125&month=summer&sdate=0101&edate=1231&sector=state&src=ncei_climate91&state=MN&opt=both&var=avg_high&cmap=autumn&cmap_r=on&_r=t&dpi=100&_fmt=png

#Load in data
MN_precip <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/MN_precip.csv")

#select for needed columns
MN_precip <- MN_precip %>% 
  dplyr::select(station, lon, lat, total_precip, avg_high, avg_low, avg_temp, max_years)
                
#see where these stations are; filter based on PRECIP
ggplot(MN_precip, aes(x = lon, y = lat, color = total_precip)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_gradient(
    name = "Climatological average precipitation total JJA",
    low = "green",
    high = "darkblue"
  ) +
  theme_minimal() +
  labs(
    title = "Minnesota Mesonet Stations",
    x = "Longitude",
    y = "Latitude"
  )

#see where these stations are; filter based on avg high temp
ggplot(MN_precip, aes(x = lon, y = lat, color = avg_high)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_gradient(
    name = "Avg High (F)",
    low = "red",
    high = "darkblue"
  ) +
  theme_minimal() +
  labs(
    title = "Minnesota Mesonet Stations",
    x = "Longitude",
    y = "Latitude"
  )

#see where these stations are; filter based on avg low temp
ggplot(MN_precip, aes(x = lon, y = lat, color = avg_low)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_gradient(
    name = "Avg Low (F)",
    low = "darkblue",
    high = "lightblue2"
  ) +
  theme_minimal() +
  labs(
    title = "Minnesota Mesonet Stations",
    x = "Longitude",
    y = "Latitude"
  )

# now apply k-means clustering -----------------------------------------------------------
#to all three variables to see if there is agreeance on ranges

#assign 2 k's so we can split N/S

cluster_data <- MN_precip %>%
  select(avg_high, avg_low, total_precip) 

scaled_data <- scale(cluster_data)

set.seed(123)  #for reproducibility
kmeans_result <- kmeans(scaled_data, centers = 2)  

#add cluster assignments back to original dataset
MN_precip$cluster <- factor(kmeans_result$cluster)

#look at results; differs mostly by temp, but a bit by precip too
kmeans_result$centers

#Cluster  avg_high  avg_low  total_precip
#1         0.59      0.63      0.33
#2        -1.01     -1.07     -0.55


#plot
ggplot(MN_precip, aes(x = lon, y = lat)) +
  borders("state", regions = "minnesota", fill = "gray95", color = "grey70") +
  geom_point(aes(color = cluster), size = 2) +
  coord_fixed(1.3) +
  labs(title = "K-means Clustering of Mesonet Stations",
       color = "Cluster",
       x = "Longitude",
       y = "Latitude"
       ) +
  theme_minimal()

#elbow method for validation: elbow bend is at 2! validation for separating N/S

#create vector
wss <- vector()

for (k in 1:10) {
  km <- kmeans(scale(cluster_data), centers = k, nstart = 25)
  wss[k] <- km$tot.withinss
}

#plot elbow
plot(1:10, wss, type = "b",
     xlab = "Number of Clusters K",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method for Choosing K")


## LAGOS DATA CLEANING ------------------------------------------------------------------------------

#First filter all datasets for MN only = state 14

LAGOSNE_lakeslocus101 <- LAGOSNE_lakeslocus101 %>%
  filter(state_zoneid == "State_14")

#filter for only lagos lake IDs that we need and select columns needed

epi_waterquality10873 <- epi_waterquality10873 %>%
  filter(lagoslakeid %in% unique_lakeIDs$link_lagosid) %>% 
  dplyr::select(lagoslakeid, secchi, sampledate, eventida10873, programid,  chla, colora, doc, nh4, no2, no2no3, srp, tdn, tdp, tkn, tn, tp)

#only lake IDs that we need

chemistry_limno <- chemistry_limno %>%
  filter(lagoslakeid %in% unique_lakeIDs$link_lagosid)

##LAGOS data processing and category separation --------------------------------

#using unique_lakeIDs as my sorting dataset
# First I need data so that i can categorize my lakes..

#find most recent secchi measurement for each lake in lagos data
#first used LAGOS-NE data for this but it only went to 2012, went back and switched out for LAGOS-US data
#LAGOS-US most recent dates range generally from 2015-2020

recent_secchi <- chemistry_limno %>%
  filter(parameter_name == "secchi_m", !is.na(sample_date)) %>%
  group_by(lagoslakeid) %>%
  slice_max(order_by = sample_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(lagoslakeid, sample_date, parameter_name, parameter_value) %>%
  rename(secchi = parameter_value)

#rename link_lagosid to lagoslakeid so i can merge with LAGOS datasets
unique_lakeIDs <- unique_lakeIDs %>%
  rename(lagoslakeid = link_lagosid)

#now merge in secchi depth using lagoslakeid
unique_lakeIDs <- unique_lakeIDs %>%
  left_join(recent_secchi %>% select(lagoslakeid, secchi), by = "lagoslakeid")


# From secchi depth, calculate TSI for each lake: 60 - (14.41 * ln(secci #)) and put in new column

unique_lakeIDs <- unique_lakeIDs %>%
  mutate(TSI = 60 - (14.41 * log(secchi)))

#remove any nan or inf, takes us down to 407 lakes
unique_lakeIDs <- unique_lakeIDs %>%
  filter(is.finite(TSI))

# Look at histogram of TSI to see spread
ggplot(unique_lakeIDs, aes(x = TSI)) +
  geom_histogram(binwidth = 2, color = "black", fill = "skyblue")+
  theme_minimal()

# Then assign oligotrophic as <40, and eutrophic as >50 if this makes sense; based on carlson TSI classification

unique_lakeIDs <- unique_lakeIDs %>%
  mutate(trophic_state = case_when(
    TSI < 40 ~ "Oligotrophic",
    TSI > 50 ~ "Eutrophic",
    TRUE ~ "Mesotrophic" 
  ))

#see how many is in each column
table(unique_lakeIDs$trophic_state)

#Eutrophic  Mesotrophic Oligotrophic 
#169          154           84


#color histogram by trophic state
ggplot(unique_lakeIDs, aes(x = TSI, fill = trophic_state)) +
  geom_histogram(binwidth = 2, color = "black", position = "identity", alpha = 0.7) +
  labs(title = "Histogram of TSI by Trophic State",
       x = "Trophic State Index (TSI)",
       y = "Number of Lakes",
       fill = "Trophic State") +
  scale_fill_manual(values = c("Oligotrophic" = "skyblue", "Eutrophic" = "forestgreen")) +
  theme_minimal()

#color on map by trophic state; dont include mesotrophic lakes
ggplot(unique_lakeIDs %>% filter(trophic_state != "Mesotrophic"),
       aes(x = longitude, y = latitude, color = trophic_state)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_manual(
    name = "Trophic State",
    values = c("Oligotrophic" = "skyblue", "Eutrophic" = "forestgreen")
  ) +
  theme_minimal() +
  labs(
    title = "Minnesota Lakes Trophic State",
    x = "Longitude",
    y = "Latitude"
  )

#even though mean depth might be a better depth separator (according to Wetzel's Limno), but only have max depth available through LAGOS and profiles dataset

#use 33% quantiles as shallow and deep cutoffs
#this aligns with <10m or less being shallow for stratifying lakes

mean_depth <- mean(unique_lakeIDs$lake_maxdepth_m, na.rm = TRUE) #14.6
sd_depth <- sd(unique_lakeIDs$lake_maxdepth_m, na.rm = TRUE) #9.6

cutoffs <- quantile(unique_lakeIDs$lake_maxdepth_m, probs = c(0.33, 0.66), na.rm = TRUE)

shallow_cut <- cutoffs[1]
deep_cut <- cutoffs[2]

#classify lakes
unique_lakeIDs <- unique_lakeIDs %>%
  mutate(depth_category = case_when(
    lake_maxdepth_m < shallow_cut ~ "Shallow",
    lake_maxdepth_m > deep_cut ~ "Deep",
    TRUE ~ "Medium"  
  ))

#plot
ggplot(unique_lakeIDs, aes(x = lake_maxdepth_m)) +
  geom_histogram(binwidth = 2, color = "black", fill = "skyblue") +
  geom_vline(xintercept = shallow_cut, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = deep_cut, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Histogram of Lake Max Depth with Quantile-Based Cutoffs",
       x = "Max Depth (m)",
       y = "Number of Lakes") +
  theme_minimal()

#see how many is in each column
table(unique_lakeIDs$depth_category)

#Deep  Medium Shallow 
#139     134     134

#count categories; there is some crossover, but not entirely 
unique_lakeIDs %>%
  count(trophic_state, depth_category)

#Eutrophic     Deep              19
#Eutrophic     Medium            55
#Eutrophic     Shallow           95

#Oligotrophic  Deep              61
#Oligotrophic  Medium            18
#Oligotrophic  Shallow            5


#now for hydrology category, seepage (or "isolated" vs. drainage)

#from LAGOS dataset: 
#lake_connectivity_permanent = hydrologic connectivity class of the focal lake determined from the NHD network considering only permanent flow

#filted by only unique_lakeIDs lagos numbers
lake_characteristics <- lake_characteristics %>%
  filter(lagoslakeid %in% unique_lakeIDs$lagoslakeid)

#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_permanent)

#Drainage DrainageLk  Headwater   Isolated   Terminal TerminalLk 
#52        174         52        110         13          6


#I AM DECIDING TO ONLY USE PERMANENT FLOW.... NOT CONSIDERING BOTH PERMAMENT AND INTERMITTENT-EPHERMERAL FLOW -------

#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_fluctuates)

#N   Y 
#349  58

#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_class)

#Drainage DrainageLk  Headwater   Isolated   Terminal TerminalLk 
#62        202         41         91          8          3


#code definitions --
#drainage: lake with inlets and outlets; no upstream lakes greater than or equal to 10 ha
#drainagelk: lake with inlets and outlets; has upstream lakes greater than or equal to 10 ha
#isolated: lake with neither an inlet or an outlet

#so i am going to combine drainage and drainagelk into one category? to get seepage vs. drainage

#merge in connectivity column --------------

#now merge using lagoslakeid
unique_lakeIDs <- unique_lakeIDs %>%
  left_join(lake_characteristics %>% select(lagoslakeid, lake_connectivity_permanent), by = "lagoslakeid")

#assign categories
unique_lakeIDs <- unique_lakeIDs %>%
  mutate(hydro_category = case_when(
    lake_connectivity_permanent == "Isolated" ~ "Seepage",
    lake_connectivity_permanent == "Drainage" ~ "Drainage",
    lake_connectivity_permanent == "DrainageLk" ~ "Drainage",
    TRUE ~ "Other" 
  ))

#see how many is in each column
table(unique_lakeIDs$hydro_category)

#Drainage    Other  Seepage 
#226       71      110


#missing lat and lon in unique_lakeIDs for lake 79457 -- pull from LAGOSNE_lakeslocus 101

#find correct lat and lon: lat is 47.1, lon is -94.9
LAGOSNE_lakeslocus101 %>%
  filter (lagoslakeid == 79457)

#now fill that into unique_lakeIDs
unique_lakeIDs <- unique_lakeIDs %>%
  mutate(
    latitude = if_else(lagoslakeid == 79457, 47.07318, latitude),
    longitude = if_else(lagoslakeid == 79457, -94.91581, longitude)
  )

#check that it worked; it did
unique_lakeIDs %>%
  filter (lagoslakeid == 79457)

#now assign north and south to all lakes based on "K-means Clustering of Mesonet Stations" figure, 46 LAT looks like good divider

unique_lakeIDs <- unique_lakeIDs %>%
  mutate(geo_position = if_else(latitude > 46, "North", "South"))

#see how many is in each column
table(unique_lakeIDs$geo_position)

#North South 
#141   266 


##make a figure for prelim proposal --------------------------------------------------

p1 <- ggplot(unique_lakeIDs) +
  aes(x = longitude, y = latitude, color = geo_position) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_manual(
    name = "Location",
    values = c("North" = "darkorchid4", "South" = "orangered3")
  ) +
  theme_minimal() +
  labs(
    y = "Latitude",
    x = NULL
  )+ 
  theme(
    axis.title = element_text(size = 16),  # Axis title size
    legend.text = element_text(size = 12),  # Increases the size of the legend text
    legend.title = element_text(size = 12)# Facet label size (if used)
  )


p2 <- ggplot(unique_lakeIDs %>% filter(trophic_state != "Mesotrophic"),
       aes(x = longitude, y = latitude, color = trophic_state)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_manual(
    name = "Trophic State",
    values = c("Oligotrophic" = "royalblue3", "Eutrophic" = "forestgreen")
  ) +
  theme_minimal() +
  labs(
    y = NULL,
    x = NULL
  ) +
  theme(
    axis.title = element_text(size = 16),  # Axis title size
    legend.text = element_text(size = 12),  # Increases the size of the legend text
    legend.title = element_text(size = 12)# Facet label size (if used)
  )


p3 <- ggplot(unique_lakeIDs %>% filter(depth_category != "Medium"),
  aes(x = longitude, y = latitude, color = depth_category)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_manual(
    name = "Depth Category",
    values = c("Deep" = "darkblue", "Shallow" = "royalblue2")
  ) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) + 
    theme(
      axis.title = element_text(size = 16),  # Axis title size
      legend.text = element_text(size = 12),  # Increases the size of the legend text
      legend.title = element_text(size = 12)# Facet label size (if used)
    )

p4 <- ggplot(unique_lakeIDs %>% filter(hydro_category != "Other"), 
  aes(x = longitude, y = latitude, color = hydro_category)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_manual(
    name = "Hydrology",
    values = c("Seepage" = "forestgreen", "Drainage" = "turquoise3")
  ) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = NULL
  ) + 
  theme(
    axis.title = element_text(size = 16),  # Axis title size
    legend.text = element_text(size = 12),  # Increases the size of the legend text
    legend.title = element_text(size = 12)# Facet label size (if used)
  )

#combine
grid.arrange(p1, p2, p3, p4, ncol = 2)


# calculate WATER QUALITY PREDICTORS ------------------------------------------------------------

#filter LAGOS-NE data for only my unique lakes
epi_waterquality10873 <- epi_waterquality10873 %>%
  filter(lagoslakeid %in% unique_lakeIDs$lagoslakeid)

#filter for selected predictors
epi_waterquality10873 <- epi_waterquality10873 %>%
  dplyr::select (eventida10873, lagoslakeid, sampledate, chla, doc, nh4, no2no3, tp, secchi)

#add new month and year columns
epi_waterquality10873 <- epi_waterquality10873 %>%
  mutate(sampledate = as.Date(sampledate)) %>%
  mutate(month = month(sampledate))  %>%
  mutate(year = year(sampledate))

#see what years there are? 1938-2012
epi_waterquality10873 %>%
  distinct(year) %>%
  arrange(year)

#trim to 1980-2012
epi_waterquality10873 <- epi_waterquality10873 %>%
  filter(sampledate >= 1980)


#now trying to fill in 2012-2020 with LAGOS-US dataset...

#filter for my unique lakes
chemistry_limno <- chemistry_limno %>%
  filter(lagoslakeid %in% unique_lakeIDs$lagoslakeid)

#see what parameters are included in LAGOS-US dataset
chemistry_limno %>%
  distinct(parameter_name) %>%
  arrange(parameter_name)

#now filter for variables to match NE ones
chemistry_limno <- chemistry_limno %>%
  filter(parameter_name %in% c("chla_ugl", "doc_mgl", "secchi_m", "nh4n_ugl", "no2no3n_ugl", "tp_ugl"))


# trim data to only include 2012-2020
chemistry_limno <- chemistry_limno %>%
  filter(sample_date >= 2012 & sample_date <= 2020)

#average within months
chemistry_monthly_avg <- chemistry_limno %>%
  mutate(
    year = year(sample_date),
    month = month(sample_date)
  ) %>%
  group_by(lagoslakeid, year, month, parameter_name) %>%
  summarise(
    mean_value = mean(parameter_value, na.rm = TRUE),
    .groups = "drop"
  )

#pivot to wide
chemistry_wide <- chemistry_monthly_avg%>%
  pivot_wider(
    names_from = parameter_name,
    values_from = mean_value
  )

#nh4n_ugl", "no2no3n_ugl", doc_mgl are not consistently available for most lakes across the years, remove from predictor list
chemistry_wide <- chemistry_wide %>%
  select(-nh4n_ugl, -no2no3n_ugl, -doc_mgl)

#remove from lagos-NE data too
epi_waterquality10873 <- epi_waterquality10873 %>%
  select(-nh4, -no2no3, -doc)


#now average within months to match lagos-us
lagosne_monthlyave <- epi_waterquality10873 %>%
  group_by(lagoslakeid, year, month) %>%
  summarise(
    chla = mean(chla, na.rm = TRUE),
    tp = mean(tp, na.rm = TRUE),
    secchi = mean(secchi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(chla, tp, secchi), ~ ifelse(is.nan(.), NA, .)))


##need to combine datasets now; chemistry_wide with lagosne_monthlyave

#rename some columns in lagosne_monthlyave to make them match
lagosne_monthlyave <- lagosne_monthlyave %>%
  rename(chla_ugl = chla) %>%
  rename(tp_ugl = tp) %>%
  rename(secchi_m = secchi)

#combine into new dataset using lagoslakeid
LAGOSwaterquality <- bind_rows(lagosne_monthlyave, chemistry_wide)

#arrange
LAGOSwaterquality <- LAGOSwaterquality %>%
  arrange(lagoslakeid, year, month)

#still has values <1980... filter again
LAGOSwaterquality <- LAGOSwaterquality %>%
  filter(year >= 1980)


#pull surface water temp from profiles data (to include as a predictor) and average by month
surface_temps <- MN_selected_profiles %>%
  filter(depth_m <= 1) %>%
  group_by(link_lagosid, year, month) %>%
  summarise(surface_temp_c = mean(temp_c, na.rm = TRUE), .groups = "drop")

#merge into water quality predictor table
LAGOSwaterquality <- LAGOSwaterquality %>%
  left_join(surface_temps, by = c("lagoslakeid" = "link_lagosid", "year", "month"))


#now i need to average within months to create pre-season and concurrent predictors...


#### STARTING TO FORM PREDICTOR DATA TABLES ---------------------------------------------

##reminder: its okay if a lake doesnt have profile for a year, but predictors must be able to match all profile years available

#try may-june AND june for pre-season predictors
preseason_waterq_predictors <- LAGOSwaterquality %>%
  filter(month %in% c(5, 6)) %>%  # only keep May and June
  group_by(lagoslakeid, year) %>%
  summarise(
    chla_mayjune = mean(chla_ugl, na.rm = TRUE),
    chla_june    = mean(chla_ugl[month == 6], na.rm = TRUE),
    
    tp_mayjune   = mean(tp_ugl, na.rm = TRUE),
    tp_june      = mean(tp_ugl[month == 6], na.rm = TRUE),
    
    secchi_mayjune = mean(secchi_m, na.rm = TRUE),
    secchi_june    = mean(secchi_m[month == 6], na.rm = TRUE),
    
    wtemp_mayjune   = mean(surface_temp_c, na.rm = TRUE),
    wtemp_june      = mean(surface_temp_c[month == 6], na.rm = TRUE),
    
    .groups = "drop"
  )

#there are a lot of NaNs here... 36,550
sum(is.nan(as.matrix(preseason_waterq_predictors)))

#see which columns they are mostly in: 
sapply(preseason_waterq_predictors, function(x) sum(is.nan(x)))

#chla_mayjune  chla_june     tp_mayjune   tp_june  secchi_mayjune    secchi_june   wtemp_mayjune     wtemp_june 
#     3820       4136           6928        7159            129            388           6904           7086


## TRIM down to only year for each lake that have profiles in MN_selectedprofiles
preseason_waterq_predictors <- preseason_waterq_predictors %>%
  semi_join(
    MN_selected_profiles %>% distinct(lagoslakeid = link_lagosid, year),
    by = c("lagoslakeid", "year")
  )

#try july-aug-sep, and just aug
concurrent_waterq_predictors <- LAGOSwaterquality %>%
  filter(month %in% c(7, 8, 9)) %>%  # only keep May and June
  group_by(lagoslakeid, year) %>%
  summarise(
    chla_JAS = mean(chla_ugl, na.rm = TRUE),
    chla_aug    = mean(chla_ugl[month == 8], na.rm = TRUE),
    
    tp_JAS   = mean(tp_ugl, na.rm = TRUE),
    tp_aug      = mean(tp_ugl[month == 8], na.rm = TRUE),
    
    secchi_JAS = mean(secchi_m, na.rm = TRUE),
    secchi_aug    = mean(secchi_m[month == 8], na.rm = TRUE),
    
    wtemp_JAS   = mean(surface_temp_c, na.rm = TRUE),
    wtemp_aug      = mean(surface_temp_c[month == 8], na.rm = TRUE),
    
    .groups = "drop"
  )

#counts NaNs.. 37,830
sum(is.nan(as.matrix(concurrent_waterq_predictors))) 

#see which columns they are mostly in: 
sapply(concurrent_waterq_predictors, function(x) sum(is.nan(x)))

#chla_JAS    chla_aug      tp_JAS      tp_aug  secchi_JAS   secchi_aug   wtemp_JAS   wtemp_aug 
#3835           4372        7086        7434         118         579        7059        7347

## TRIM down to only year for each lake that have profiles in MN_selectedprofiles
concurrent_waterq_predictors <- concurrent_waterq_predictors %>%
  semi_join(
    MN_selected_profiles %>% distinct(lagoslakeid = link_lagosid, year),
    by = c("lagoslakeid", "year")
  )

#should i try to fill the NaNs? 


#predictors needed for all years: ice-on time, water level

#try calculating thermocline depth? to use as a predictor?


# Calculate CLIMATE PREDICTORS -----------------------------------------------------------------------------------------

#PAR, surface wind speed, air temp, precip, ENSO, PDO, NAM, geopotential height 500

#look back at mesonet data?
#PRISM data?

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


#ENSO index -----
#which one?

#PNA index





# Calculate VHT for 5 selected lake profiles ----------------------------------------------------------------------------------------

#calculate VHT function -------------

#Minnesota: based on TDO3… so 3mg/L and 21.5 C
#from MDNR water quality standard guide

#first take a subset of lakes to apply it to and make sure this code is working accurately

first_five_lakes <- MN_selected_profiles %>%
  distinct(link_lagosid) %>%
  slice(1:5)

firstfive_lakedata <- MN_selected_profiles %>%
  filter(link_lagosid %in% first_five_lakes$link_lagosid) %>%
  dplyr::select(link_lagosid, date, month, year, depth_m, temp_c, do_ppm, lake_maxdepth_m)


# Calculate VHT components for each lake profile
firstfive_lakedata <- firstfive_lakedata %>%
  # Group by lake and sample date to work with each profile separately
  group_by(link_lagosid, date, month, year) %>%
  arrange(link_lagosid, date, depth_m) %>%
  
  # Calculate temp_depth: first depth where temp_c <= 22
  mutate(
    temp_depth = {
      temp_indices <- which(temp_c <= 22)
      if (length(temp_indices) == 0) {
        0
      } else {
        depth_m[temp_indices[1]]
      }
    }
  ) %>%
  
  # Calculate o2_depth: first depth where do_ppm <= 4
  mutate(
    o2_depth = {
      o2_indices <- which(do_ppm <= 4)
      if (length(o2_indices) == 0) {
        first(lake_maxdepth_m)
      } else {
        depth_m[o2_indices[1]]
      }
    }
  ) %>%
  
  # Calculate VHT (assuming it's the difference between o2_depth and temp_depth)
  mutate(
    VHT = pmax(0, o2_depth - temp_depth)  # Ensure VHT is not negative
  ) %>%
  
  ungroup()

#check summary of results
summary(firstfive_lakedata[c("temp_depth", "o2_depth", "VHT")])

#temp_depth        o2_depth          VHT       
#Min.   : 0.000   Min.   : 3.96   Min.   : 0.00  
#1st Qu.: 0.170   1st Qu.:13.02   1st Qu.:10.92  
#Median : 0.450   Median :21.95   Median :18.02  
#Mean   : 2.537   Mean   :18.98   Mean   :16.49  
#3rd Qu.: 4.970   3rd Qu.:24.03   3rd Qu.:21.90  
#Max.   :15.940   Max.   :28.95   Max.   :28.92

#simplify
profile_summary <- firstfive_lakedata %>%
  group_by(link_lagosid, date, month, year) %>%
  summarise(
    temp_depth = first(temp_depth),
    o2_depth = first(o2_depth),
    VHT = first(VHT),
    max_depth_sampled = max(depth_m, na.rm = TRUE),
    lake_maxdepth = first(lake_maxdepth_m),
    .groups = "drop"
  )

#calculate monthly vht to plot
monthly_vht <- profile_summary %>%
  group_by(link_lagosid, month) %>%
  summarise(
    mean_VHT = mean(VHT, na.rm = TRUE),
    .groups = "drop"
  )

#plot seasonality

vht_seasonal_plot <- ggplot(monthly_vht, aes(x = month, y = mean_VHT, color = factor(link_lagosid))) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_x_continuous(breaks = 1:12, 
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  scale_color_viridis_d(name = "Lake ID") +  # Use viridis color palette for better distinction
  labs(
    title = "Seasonal Variation in Vertical Habitat Thickness (VHT)",
    x = "",
    y = "Mean VHT (m)",
    color = "LAGOS Lake ID"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.major.y = element_line(color = "gray90", size = 0.5)
  )

# Display the plot
print(vht_seasonal_plot)

# Calculate VHT for full lake profiles ----------------------------------------------------------------------------------------
#MN_selected_profiles 

#these lakes need their max depth's edited before running code again

# A tibble: 6 × 4
#link_lagosid date       temp_depth lake_maxdepth_m
#<int> <date>          <dbl>           <dbl>
#  1          197 1984-08-14       6               5.5 
#2          305 2017-06-07       5               2.13
#3          305 2017-08-29       5               2.13
#4          381 2006-07-24       4.5             4.27
#5          734 2009-07-14       7.77            7.62
#6         1171 2013-07-16      22.9            10.4 

MN_selected_profiles <- MN_selected_profiles %>%
  group_by(link_lagosid, date, month, year) %>%
  arrange(link_lagosid, date, depth_m) %>%
  
  # Calculate temp_depth: first depth where temp_c <= 22
  mutate(
    temp_depth = {
      valid_indices <- which(temp_c <= 22 & depth_m >= 0 & !is.na(temp_c))
      if (length(valid_indices) == 0) {
        # If temp never reaches 22°C, use the maximum depth (no suitable cool water)
        ifelse(!is.na(first(lake_maxdepth_m)), 
               first(lake_maxdepth_m), 
               max(depth_m[depth_m >= 0], na.rm = TRUE))
      } else {
        pmax(0, depth_m[valid_indices[1]])
      }
    }
  ) %>%
  
  # Calculate o2_depth: first depth where do_ppm <= 4
  mutate(
    o2_depth = {
      if (all(is.na(do_ppm))) {
        NA_real_  # Only return NA if we have no oxygen data at all
      } else {
        valid_indices <- which(do_ppm <= 4 & depth_m >= 0 & !is.na(do_ppm))
        if (length(valid_indices) == 0) {
          # Oxygen never drops below 4, use max depth (all water is oxygenated)
          ifelse(!is.na(first(lake_maxdepth_m)), 
                 first(lake_maxdepth_m), 
                 max(depth_m[depth_m >= 0], na.rm = TRUE))
        } else {
          pmax(0, depth_m[valid_indices[1]])
        }
      }
    }
  ) %>%
  
  # Calculate VHT: thickness of suitable habitat
  mutate(
    VHT = case_when(
      is.na(o2_depth) ~ NA_real_,  # Only NA if no oxygen data
      temp_depth >= o2_depth ~ 0,  # No suitable habitat
      TRUE ~ o2_depth - temp_depth  # Thickness of suitable habitat
    )
  ) %>%
  
  ungroup()

#check summary of results
summary(MN_selected_profiles[c("temp_depth", "o2_depth", "VHT")])

#temp_depth        o2_depth          VHT        
#Min.   : 0.000   Min.   : 0.00   Min.   : 0.000  
#1st Qu.: 0.000   1st Qu.: 5.00   1st Qu.: 1.000  
#Median : 0.230   Median : 7.00   Median : 6.000  
#Mean   : 2.085   Mean   :10.60   Mean   : 8.665  
#3rd Qu.: 4.020   3rd Qu.:11.00   3rd Qu.:10.333  
#Max.   :22.894   Max.   :63.36   Max.   :63.359

##QC CHECKS ------------------------------------------------------------------------------

# Check for any remaining issues
cat("Data integrity checks:\n")
cat("Any negative temp_depth:", any(MN_selected_profiles$temp_depth < 0, na.rm = TRUE), "\n")
cat("Any negative o2_depth:", any(MN_selected_profiles$o2_depth < 0, na.rm = TRUE), "\n")
cat("Any negative VHT:", any(MN_selected_profiles$VHT < 0, na.rm = TRUE), "\n")
cat("Any NA values in temp_depth:", sum(is.na(MN_selected_profiles$temp_depth)), "\n")
cat("Any NA values in o2_depth:", sum(is.na(MN_selected_profiles$o2_depth)), "\n")
cat("Any NA values in VHT:", sum(is.na(MN_selected_profiles$VHT)), "\n")

cat("Validation checks:\n")
cat("Profiles with NA oxygen data:", sum(is.na(MN_selected_profiles$o2_depth)), "\n")
cat("Profiles with NA VHT:", sum(is.na(MN_selected_profiles$VHT)), "\n")
cat("NAs match (should be TRUE):", sum(is.na(MN_selected_profiles$o2_depth)) == sum(is.na(MN_selected_profiles$VHT)), "\n\n")

logic_check <- MN_selected_profiles %>%
  filter(!is.na(VHT)) %>%
  mutate(
    # VHT should be 0 when temp_depth >= o2_depth
    should_be_zero = temp_depth >= o2_depth,
    actually_zero = VHT == 0,
    # VHT should equal o2_depth - temp_depth when temp_depth < o2_depth
    should_be_positive = temp_depth < o2_depth,
    correct_calculation = abs(VHT - (o2_depth - temp_depth)) < 0.001
  )

cat("Logic validation:\n")
cat("Cases where VHT should be 0:", sum(logic_check$should_be_zero), "\n")
cat("Cases where VHT actually is 0:", sum(logic_check$actually_zero), "\n")
cat("Cases where calculation is correct:", sum(logic_check$correct_calculation), "\n")
cat("Total non-NA cases:", nrow(logic_check), "\n\n")

#Look at the distribution of cases
case_summary <- MN_selected_profiles %>%
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  mutate(
    case_type = case_when(
      VHT == 0 ~ "No suitable habitat (VHT = 0)",
      VHT > 0 ~ "Suitable habitat exists (VHT > 0)"
    )
  ) %>%
  count(case_type)

cat("Distribution of cases:\n")
print(case_summary)

# Check that VHT = o2_depth - temp_depth (should be TRUE for all rows)
VHT_check <- MN_selected_profiles %>%
  mutate(
    calculated_VHT = pmax(0, o2_depth - temp_depth),
    VHT_matches = abs(VHT - calculated_VHT) < 0.001  # Allow for small floating point errors
  )

cat("VHT calculation consistency:", all(VHT_check$VHT_matches), "\n")
cat("Rows with VHT inconsistency:", sum(!VHT_check$VHT_matches), "\n")

# Look at a few random profiles to manually verify logic
set.seed(123)
sample_profiles <- MN_selected_profiles %>%
  select(link_lagosid, date) %>%
  distinct() %>%
  sample_n(5)

# Check these profiles in detail
for(i in 1:nrow(sample_profiles)) {
  profile_data <- MN_selected_profiles %>%
    filter(link_lagosid == sample_profiles$link_lagosid[i],
           date == sample_profiles$date[i]) %>%
    arrange(depth_m) %>%
    select(depth_m, temp_c, do_ppm, temp_depth, o2_depth, VHT)
  
  cat("\n=== Profile", i, "===\n")
  cat("Lake:", sample_profiles$link_lagosid[i], "Date:", sample_profiles$date[i], "\n")
  print(profile_data)
  
  # Check if temp_depth makes sense
  first_temp_22 <- profile_data %>%
    filter(temp_c <= 22, depth_m >= 0) %>%
    slice_head(n = 1) %>%
    pull(depth_m)
  
  if(length(first_temp_22) > 0) {
    cat("Expected temp_depth:", first_temp_22[1], "| Actual:", first(profile_data$temp_depth), "\n")
  }
  
  # Check if o2_depth makes sense
  first_o2_4 <- profile_data %>%
    filter(do_ppm <= 4, depth_m >= 0) %>%
    slice_head(n = 1) %>%
    pull(depth_m)
  
  if(length(first_o2_4) > 0) {
    cat("Expected o2_depth:", first_o2_4[1], "| Actual:", first(profile_data$o2_depth), "\n")
  }
}

# Check for unrealistic values
cat("\nDistribution checks:\n")
cat("Profiles where temp_depth > max lake depth:\n")
weird_temp <- MN_selected_profiles %>%
  filter(temp_depth > lake_maxdepth_m) %>%
  select(link_lagosid, date, temp_depth, lake_maxdepth_m) %>%
  distinct()
print(head(weird_temp))


cat("Profiles where o2_depth > max lake depth:\n")
weird_o2 <- MN_selected_profiles %>%
  filter(o2_depth > lake_maxdepth_m) %>%
  select(link_lagosid, date, o2_depth, lake_maxdepth_m) %>%
  distinct()
print(head(weird_o2))

# Check for extreme VHT values
cat("Profiles with VHT > 65m (potential outliers):\n")
extreme_vht <- MN_selected_profiles %>%
  filter(VHT > 65) %>%
  select(link_lagosid, date, temp_depth, o2_depth, VHT, lake_maxdepth_m) %>%
  distinct()
print(extreme_vht)

# Check for unrealistic values
cat("\nDistribution checks:\n")
cat("Profiles where temp_depth > max lake depth:\n")
weird_temp <- MN_selected_profiles %>%
  filter(temp_depth > lake_maxdepth_m) %>%
  select(link_lagosid, date, temp_depth, lake_maxdepth_m) %>%
  distinct()
print(head(weird_temp))

cat("Profiles where o2_depth > max lake depth:\n")
weird_o2 <- MN_selected_profiles %>%
  filter(o2_depth > lake_maxdepth_m) %>%
  select(link_lagosid, date, o2_depth, lake_maxdepth_m) %>%
  distinct()
print(head(weird_o2))

# Check for extreme VHT values
cat("Profiles with VHT > 100m (potential outliers):\n")
extreme_vht <- MN_selected_profiles %>%
  filter(VHT > 100) %>%
  select(link_lagosid, date, temp_depth, o2_depth, VHT, lake_maxdepth_m) %>%
  distinct()
print(extreme_vht)

# Check that each profile has the same VHT values across all depths
profile_consistency <- MN_selected_profiles %>%
  group_by(link_lagosid, date) %>%
  summarise(
    unique_temp_depths = n_distinct(temp_depth),
    unique_o2_depths = n_distinct(o2_depth),
    unique_VHT = n_distinct(VHT),
    .groups = "drop"
  ) %>%
  filter(unique_temp_depths > 1 | unique_o2_depths > 1 | unique_VHT > 1)

cat("Profiles with inconsistent VHT values:", nrow(profile_consistency), "\n")
if(nrow(profile_consistency) > 0) {
  print(head(profile_consistency))
}

# Distribution of VHT values
p1 <- ggplot(MN_selected_profiles %>% distinct(link_lagosid, date, .keep_all = TRUE), 
             aes(x = VHT)) +
  geom_histogram(bins = 50) +
  labs(title = "Distribution of VHT values")

print(p1)

plot_data <- MN_selected_profiles %>% 
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE)

ggplot(plot_data, aes(x = temp_depth, y = o2_depth, color = VHT)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_color_gradient(low = "blue", high = "yellow", name = "VHT") +
  labs(title = "temp_depth vs o2_depth", 
       subtitle = "Points below red line should have VHT = 0",
       x = "Temperature Depth (m)", 
       y = "Oxygen Depth (m)") +
  theme_minimal()

#assessment plots -----------------------------

# Monthly patterns
p1 <- MN_selected_profiles %>%
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  ggplot(aes(x = factor(month), y = VHT)) +
  geom_boxplot() +
  labs(title = "VHT by Month", 
       x = "Month", y = "VHT (m)") +
  theme_minimal()

print(p1)

# Year-over-year trends
p2 <- MN_selected_profiles %>%
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  group_by(year) %>%
  summarise(mean_VHT = mean(VHT, na.rm = TRUE),
            median_VHT = median(VHT, na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = mean_VHT, color = "Mean")) +
  geom_line(aes(y = median_VHT, color = "Median")) +
  labs(title = "VHT Trends Over Time", 
       x = "Year", y = "VHT (m)", color = "Statistic") +
  theme_minimal()

print(p2)

p3 <- MN_selected_profiles %>%
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  ggplot(aes(x = lake_maxdepth_m, y = VHT)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "VHT vs Lake Maximum Depth", 
       x = "Lake Maximum Depth (m)", y = "VHT (m)") +
  theme_minimal()

print(p3)

p8 <- MN_selected_profiles %>%
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  mutate(
    VHT_category = case_when(
      VHT == 0 ~ "No Habitat",
      VHT <= 2 ~ "Limited (0-2m)",
      VHT <= 5 ~ "Moderate (2-5m)",
      VHT <= 10 ~ "Good (5-10m)",
      TRUE ~ "Excellent (>10m)"
    ),
    VHT_category = factor(VHT_category, levels = c("No Habitat", "Limited (0-2m)", 
                                                   "Moderate (2-5m)", "Good (5-10m)", "Excellent (>10m)"))
  ) %>%
  ggplot(aes(x = VHT_category, fill = VHT_category)) +
  geom_bar() +
  labs(title = "Distribution of VHT Categories", 
       x = "VHT Category", y = "Number of Profiles") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p8)

library(maps)

# Get Minnesota map
minnesota <- map_data("state", region = "minnesota")


# Seasonal VHT patterns
seasonal_summary <- MN_selected_profiles %>%
  filter(!is.na(VHT)) %>%
  mutate(
    season = case_when(
      month %in% c(3, 4, 5) ~ "Spring", 
      month %in% c(7, 8, 9, 10) ~ "Summer"
    )
  ) %>%
  group_by(link_lagosid, latitude, longitude, season) %>%
  summarise(
    avg_VHT = mean(VHT, na.rm = TRUE),
    n_profiles = n(),
    .groups = "drop"
  ) %>%
  filter(n_profiles >= 2)  # Only lakes with multiple profiles per season

ggplot() +
  geom_polygon(data = minnesota, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = seasonal_summary, aes(x = longitude, y = latitude, color = avg_VHT), 
             size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red", name = "Average\nVHT (m)") +
  facet_wrap(~season, ncol = 2) +
  labs(title = "Seasonal VHT Patterns by Lake Location") +
  coord_fixed(1.3) +
  theme_void()

# Map showing proportion of profiles with suitable habitat
habitat_summary <- MN_selected_profiles %>%
  filter(!is.na(VHT)) %>%
  group_by(link_lagosid, latitude, longitude) %>%
  summarise(
    total_profiles = n(),
    profiles_with_habitat = sum(VHT > 0),
    prop_with_habitat = profiles_with_habitat / total_profiles,
    avg_VHT_when_available = mean(VHT[VHT > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_profiles >= 5)  # Only lakes with at least 5 profiles

ggplot() +
  geom_polygon(data = minnesota, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = habitat_summary, 
             aes(x = longitude, y = latitude, color = prop_with_habitat, size = total_profiles), 
             alpha = 0.8) +
  scale_color_gradient(low = "red", high = "green", name = "Proportion\nwith Habitat", 
                       labels = scales::percent) +
  scale_size_continuous(name = "Total\nProfiles", range = c(1, 4)) +
  labs(title = "Proportion of Profiles with Suitable Habitat (VHT > 0)") +
  coord_fixed(1.3) +
  theme_void()


# proceed -----------------

#simplify
MN_VHT_summary <- MN_selected_profiles %>%
  group_by(link_lagosid, date, month, year) %>%
  summarise(
    temp_depth = first(temp_depth),
    o2_depth = first(o2_depth),
    VHT = first(VHT),
    max_depth_sampled = max(depth_m, na.rm = TRUE),
    lake_maxdepth = first(lake_maxdepth_m),
    .groups = "drop"
  )

#calculate JAS and add to both predictor tables as VHT





#look at full seasonality averages of VHT, north vs south MN






