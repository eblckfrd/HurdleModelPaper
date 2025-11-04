# Dissertation Chapter 2 - MINNESOTA - started June 2025
# Processing of Minnesota data from MGLP team (Mike, Gretchen Hansen, Heidi Rantala, Lesley Knoll)
# Processing LAGOS data
# Gathering and cleaning all lake predictor data

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
library(data.table)

#----------------------------------------------------------------------------------------------------------

# Load in MN_profiles data from Gretchen's shared google drive

FULL_MN_profiles <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/MN_profiles.csv")

#-----------------------

# Load in LAGOS data files downloaded from EDI
#was originally using only LAGOS-NE data but it is less up to date than LAGOS-US
#so for water quality I can pull full country then filter to just my selected lakes, and later merge NE with US pre 2012 and post 2012

LAGOSNE_state_105 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/LAGOSNE_state_105.csv")
        #state IDs
LAGOSNE_lakeslocus101 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/LAGOSNE_lakeslocus101.csv")
         #lagos ID, name, lake area, lake perim, hu zones, elevation, stat and county zone id
lake_characteristics <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/lake_characteristics.csv")
        #hydrology/lake connectivity
epi_waterquality10873 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/epi_waterquality10873.csv")
        #LAGOS-NE water quality params
chemistry_limno <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/chemistry_limno.csv")
      #this is from LAGOS-US, all water quality parameters

#-----------------------

#ice duration from Minnesota Lake Ice Phenology dataset: https://conservancy.umn.edu/items/940c6197-486b-437c-96ca-cca9b4534dfa

ice_data <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/ice_duration.csv")

#DOW : Department of Waterbody Lake ID
#winter.year : If the ice came on in fall of 2000 and melted in spring of 2001 the winter.year variable would be 2001
    #using this one
#fall.year : The fall year of the observation. If the ice came on in fall of 2000 and melted in the spring of 2001 the fall.year variable would be 2000


##Profile data processing ------------------------------------------------------------------------------------

# Only keep columns that are needed
# Filter for only lakes that stratify-- filter for geom ratio <4 (ratio already calculated in profiles dataset)

MN_profiles <- FULL_MN_profiles %>% #takes data down to 714,410 and 19 variables
  dplyr::select(link_lagosid, year, date, doy, latitude, longitude, 
                monitoringlocationidentifier, monitoringlocationname,
                profile_ident, n_dep_inprofile,
                depth_m, temp_c, do_ppm, lake_geomratio, lake_maxdepth_m, lake_area_m2,
                flag_do, flag_temp, flag_dep_diff) %>%
  filter(lake_geomratio <= "4") 

uniqueN(FULL_MN_profiles$link_lagosid) #3397 lakes

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

#expirimenting with removing the below block...

# #filter for min depth: code from Mike, based on Jane et al 2021
# valid_profiles <- MN_profiles %>% #Profiles must have a minimum depth less than 3 meters; 
#   group_by(profile_ident) %>% #otherwise the profiles shallowest read could be in 
#   summarise(min.depth = min(depth_m, na.rm = TRUE)) %>% #the metalimnion and capture none of the epilimnion
#   filter(min.depth < 3) %>% #Jane et al 2021
#   pull(profile_ident)

# MN_profiles <- MN_profiles %>%
#   filter(profile_ident %in% valid_profiles)

uniqueN(MN_profiles$link_lagosid) #1815 lakes

#replace incorrect max depths (identified in code lines ~700 and MN_VHT script and moved up in code order)

#East Moore Lake in Minneapolis, MN; DNR lists max depth as 6.7m, ID = 197
MN_profiles$lake_max_depth_m[MN_profiles$link_lagosid == 197] <- 6.7

#MN DNR says 8.5 m max depth for ID = 734
MN_profiles$lake_max_depth_m[MN_profiles$link_lagosid == 734] <- 8.5

#NEW ISSUE: remove lakes that have a max depth of 3m or less
MN_profiles <- MN_profiles %>% 
  filter(lake_maxdepth_m > 3)

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

# Filter for lakes with >10 years of data 
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

# Look at how many lakes have more than one monitoring location identifier: 364 lakes

MN_profiles_10yrs %>%
  group_by(link_lagosid) %>%
  summarise(num_sites = n_distinct(monitoringlocationidentifier), .groups = "drop") %>%
  filter(num_sites > 1) %>%
  summarise(num_lakes_with_multiple_sites = n())

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

uniqueN(MN_profiles_10yrs$link_lagosid)

# Look at the distribution of depth mismatches between # of depths in profile and lake max depth
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


#find lake profiles that cover >=50% of max depth OR is >=10m deep
good_lakes_any <- MN_profiles_10yrs %>%
  mutate(
    # Calculate what percentage of max depth this profile covers
    depth_coverage = n_dep_inprofile / lake_maxdepth_m,
    # Profile is adequate if it covers >=50% of max depth OR is >=10m deep
    adequate_profile = (depth_coverage >= 0.5) | (n_dep_inprofile >= 10)
  ) %>%
  group_by(link_lagosid) %>%
  summarise(
    has_adequate_profiles = any(adequate_profile),
    max_coverage = max(depth_coverage),
    max_profile_depth = max(n_dep_inprofile),
    n_profiles = n(),
    .groups = "drop"
  ) %>%
  filter(has_adequate_profiles) %>%
  pull(link_lagosid)

#filter
MN_profiles_filtered <- MN_profiles_10yrs %>%
  filter(link_lagosid %in% good_lakes_any)

# Now analyze temp and DO patterns in deep profiles (>=10m)
deep_profile_analysis <- MN_profiles_filtered %>%
  filter(n_dep_inprofile >= 10) %>%  # Only profiles that are 10m or deeper
  group_by(link_lagosid, profile_ident) %>%
  arrange(depth_m) %>%
  mutate(
    # Flag depths >= 10m within each profile
    deep_zone = depth_m >= 10
  ) %>%
  # First check if we have any data below 10m before calculating stats
  filter(any(deep_zone)) %>%
  summarise(
    max_depth = max(depth_m),
    n_points_below_10m = sum(deep_zone, na.rm = TRUE),
    # Get vectors of temp and DO below 10m (removing NAs)
    temp_vals_below_10m = list(temp_c[deep_zone & !is.na(temp_c)]),
    do_vals_below_10m = list(do_ppm[deep_zone & !is.na(do_ppm)]),
    .groups = "drop"
  ) %>%
  # Calculate stats only if we have enough valid measurements
  mutate(
    n_temp_below_10m = map_int(temp_vals_below_10m, length),
    n_do_below_10m = map_int(do_vals_below_10m, length),
    # Temperature stats below 10m
    temp_min_below_10m = map_dbl(temp_vals_below_10m, ~ ifelse(length(.x) >= 2, min(.x), NA_real_)),
    temp_max_below_10m = map_dbl(temp_vals_below_10m, ~ ifelse(length(.x) >= 2, max(.x), NA_real_)),
    temp_range_below_10m = temp_max_below_10m - temp_min_below_10m,
    temp_cv_below_10m = map_dbl(temp_vals_below_10m, ~ ifelse(length(.x) >= 2 & mean(.x) != 0, sd(.x) / mean(.x), NA_real_)),
    # DO stats below 10m
    do_min_below_10m = map_dbl(do_vals_below_10m, ~ ifelse(length(.x) >= 2, min(.x), NA_real_)),
    do_max_below_10m = map_dbl(do_vals_below_10m, ~ ifelse(length(.x) >= 2, max(.x), NA_real_)),
    do_range_below_10m = do_max_below_10m - do_min_below_10m,
    do_cv_below_10m = map_dbl(do_vals_below_10m, ~ ifelse(length(.x) >= 2 & mean(.x) != 0, sd(.x) / mean(.x), NA_real_))
  ) %>%
  # Remove the list columns we don't need anymore
  select(-temp_vals_below_10m, -do_vals_below_10m) %>%
  # Only keep profiles with sufficient data below 10m
  filter(n_temp_below_10m >= 2 | n_do_below_10m >= 2) %>%
  # Flag profiles where temp/DO are relatively stable below 10m
  mutate(
    temp_stable = ifelse(!is.na(temp_range_below_10m), 
                         temp_range_below_10m <= 2.0 | (!is.na(temp_cv_below_10m) & temp_cv_below_10m <= 0.1), 
                         NA),
    do_stable = ifelse(!is.na(do_range_below_10m), 
                       do_range_below_10m <= 2.0 | (!is.na(do_cv_below_10m) & do_cv_below_10m <= 0.2), 
                       NA),
    both_stable = ifelse(!is.na(temp_stable) & !is.na(do_stable), 
                         temp_stable & do_stable, 
                         NA)
  )

stable_lakes_summary <- deep_profile_analysis %>%
  group_by(link_lagosid) %>%
  summarise(
    n_deep_profiles = n(),
    n_temp_stable = sum(temp_stable, na.rm = TRUE),
    n_do_stable = sum(do_stable, na.rm = TRUE),
    n_both_stable = sum(both_stable, na.rm = TRUE),
    prop_temp_stable = n_temp_stable / n_deep_profiles,
    prop_do_stable = n_do_stable / n_deep_profiles,
    prop_both_stable = n_both_stable / n_deep_profiles,
    avg_temp_range = mean(temp_range_below_10m, na.rm = TRUE),
    avg_do_range = mean(do_range_below_10m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Flag lakes where majority of deep profiles show stable conditions
  mutate(
    lake_has_stable_temp = prop_temp_stable >= 0.5,
    lake_has_stable_do = prop_do_stable >= 0.5,
    lake_has_stable_both = prop_both_stable >= 0.5
  )

print("Summary of deep profile stability:")
print(paste("Lakes with deep profiles (>=10m):", nrow(stable_lakes_summary)))
print(paste("Lakes with stable temperature below 10m:", sum(stable_lakes_summary$lake_has_stable_temp)))
print(paste("Lakes with stable DO below 10m:", sum(stable_lakes_summary$lake_has_stable_do)))
print(paste("Lakes with both temp and DO stable below 10m:", sum(stable_lakes_summary$lake_has_stable_both)))

#lakes with deep profiles = 268
#lakes with stable temp below 10m = 178
#lakes with stable DO below 10m = 205
#lakes with both temp and DO stable below 10m = 145

# Get list of lakes that are stable below 10m (we'll keep ALL profiles from these lakes); 145 lakes
stable_lakes <- stable_lakes_summary %>%
  filter(lake_has_stable_both) %>%  # Change this criteria as needed
  pull(link_lagosid)

# Now apply the two-tier filtering approach
MN_profiles_final <- MN_profiles_filtered %>%
  mutate(
    depth_coverage = n_dep_inprofile / lake_maxdepth_m,
    is_stable_lake = link_lagosid %in% stable_lakes
  ) %>%
  filter(
    # Keep ALL profiles from stable lakes
    is_stable_lake | 
      # OR keep profiles that cover >80% of max depth from non-stable lakes
      (!is_stable_lake & depth_coverage > 0.8)
  )

print("Filtering Summary:")
print(paste("Total lakes after initial filtering:", length(unique(MN_profiles_filtered$link_lagosid)))) #392
print(paste("Lakes classified as stable below 10m:", length(stable_lakes))) #146
print(paste("Lakes in final dataset:", length(unique(MN_profiles_final$link_lagosid)))) #349
print(paste("Total profiles after initial filtering:", nrow(MN_profiles_filtered))) #399394
print(paste("Total profiles in final dataset:", nrow(MN_profiles_final))) #327732


# Check what depth coverage we're getting in non-stable lakes
non_stable_coverage <- MN_profiles_final %>%
  filter(!link_lagosid %in% stable_lakes) %>%
  mutate(depth_coverage = n_dep_inprofile / lake_maxdepth_m) %>%
  pull(depth_coverage)

if(length(non_stable_coverage) > 0) {
  print(paste("\nDepth coverage in non-stable lakes (kept profiles):"))
  print(paste("Mean coverage:", round(mean(non_stable_coverage, na.rm = TRUE), 3))) #1.116
  print(paste("Min coverage:", round(min(non_stable_coverage, na.rm = TRUE), 3))) #0.8
  print(paste("Profiles with >90% coverage:", sum(non_stable_coverage > 0.9, na.rm = TRUE))) #96546
}

print(paste("Lakes in final dataset:", length(unique(MN_profiles_final$link_lagosid))))
#349

#check what years are included here
ggplot(data = MN_profiles_final, aes(x = year)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

#---------------------------------------------------------------------

#FIGURING OUT THIS SECTION, DO NOT MOVE FORWARD UNTIL DECIDING WHAT LENGTH YEARS TO USE
  #testing 10+ year then looking at spread, can go down to 7+ years if need to

# Filter for any specific months? Mar-Oct, 371412 obs left
MN_profiles_final <- MN_profiles_final %>%
  filter(month >= 3 & month <= 10)

# Windows by month
is_summer_month <- function(m) m %in% 7:9   # Jul–Sep
is_spring_month <- function(m) m %in% 5:6   # May–Jun

# All (lake, year) combos observed (baseline to test "every year" against)
years_by_lake <- MN_profiles_final %>%
  distinct(link_lagosid, year)

# Count how many distinct years each lake has ≥1 sample in Jul–Sep
summer_years_per_lake <- MN_profiles_final %>%
  mutate(month = month(date)) %>%
  filter(month %in% 7:9) %>%                      # July–September
  distinct(link_lagosid, year) %>%               # 1 row per lake-year with a summer sample
  count(link_lagosid, name = "n_summer_years")    # tally years per lake

ggplot(data = summer_years_per_lake, aes(x = n_summer_years)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal()

# Keep lakes with at least 10 years of summer samples
lakes_with_10plus_summer_years <- summer_years_per_lake %>%
  filter(n_summer_years >= 10) %>%
  arrange(desc(n_summer_years))

#156 with +10 years, 228 with 7+ years, +259 with 5 years

hist(lakes_with_10plus_summer_years$n_summer_years)

#filter profiles dataset with only those 156 lakes

MN_profiles_final <- MN_profiles_final %>%
  filter(link_lagosid %in% lakes_with_10plus_summer_years$link_lagosid)

#----------------


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

#check that this worked; it did! 
n_distinct(MN_selected_profiles$link_lagosid)


####ignore this, ended up not being an issue with random forest models, they allow for non-consecutive samples as long as predictor years match --------

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

#try using chlorophyll-a first

# recent_chla <- chemistry_limno %>%
#   filter(parameter_name == "chla_ugl", !is.na(sample_date)) %>%
#   # Extract month and year from sample_date
#   mutate(
#     month = month(sample_date),
#     year = year(sample_date)
#   ) %>%
#   # Filter for summer months (7-9)
#   filter(month >= 7 & month <= 9) %>%
#   group_by(lagoslakeid) %>%
#   # Get the three most recent years for each lake
#   filter(year %in% tail(sort(unique(year)), 3)) %>%
#   # Calculate average chla across those years
#   summarise(
#     chla = mean(parameter_value, na.rm = TRUE),
#     n_observations = n(),
#     years_used = paste(sort(unique(year)), collapse = ", "),
#     .groups = "drop"
#   )
# 
# hist(recent_chla$chla)
# 
# #rename link_lagosid to lagoslakeid so i can merge with LAGOS datasets
# unique_lakeIDs <- unique_lakeIDs %>%
#   rename(lagoslakeid = link_lagosid)
# 
# #now merge in 
# unique_lakeIDs <- unique_lakeIDs %>%
#   left_join(recent_chla %>% select(lagoslakeid, chla), by = "lagoslakeid")
# 
# unique_lakeIDs <- unique_lakeIDs %>%
#   mutate(tsi_chla = 9.81 * log(chla) + 30.6)
# 
# # Look at histogram of TSI to see spread
# ggplot(unique_lakeIDs, aes(x = tsi_chla)) +
#   geom_histogram(binwidth = 2, color = "black", fill = "skyblue")+
#   theme_minimal()


#now try using secchi
#follow rules for using average 3 most recent years of measurements from between 7-9

recent_secchi <- chemistry_limno %>%
  filter(parameter_name == "secchi_m", !is.na(sample_date)) %>%
  # Extract month and year from sample_date
  mutate(
    month = month(sample_date),
    year = year(sample_date)
  ) %>%
  # Filter for summer months (7-9)
  filter(month >= 7 & month <= 9) %>%
  group_by(lagoslakeid) %>%
  # Get the three most recent years for each lake
  filter(year %in% tail(sort(unique(year)), 3)) %>%
  # Calculate average chla across those years
  summarise(
    secchi = mean(parameter_value, na.rm = TRUE),
    n_observations = n(),
    years_used = paste(sort(unique(year)), collapse = ", "),
    .groups = "drop"
  )

# recent_secchi <- chemistry_limno %>%
#   filter(parameter_name == "secchi_m", !is.na(sample_date)) %>%
#   # Extract month from sample_date and filter for summer months (6-9)
#   mutate(month = month(sample_date)) %>%
#   filter(month >= 7 & month <= 9) %>%
#   group_by(lagoslakeid) %>%
#   slice_max(order_by = sample_date, n = 1, with_ties = FALSE) %>%
#   ungroup() %>%
#   select(lagoslakeid, sample_date, parameter_name, parameter_value) %>%
#   rename(secchi = parameter_value)

#rename link_lagosid to lagoslakeid so i can merge with LAGOS datasets
unique_lakeIDs <- unique_lakeIDs %>%
  rename(lagoslakeid = link_lagosid)

#now merge in secchi depth using lagoslakeid
unique_lakeIDs <- unique_lakeIDs %>%
  left_join(recent_secchi %>% select(lagoslakeid, secchi), by = "lagoslakeid")

# ----


# From secchi depth, calculate TSI for each lake: 60 - (14.41 * ln(secci #)) and put in new column

unique_lakeIDs <- unique_lakeIDs %>%
  mutate(TSI = 60 - (14.41 * log(secchi)))

#remove any secchi depth of 0 (max depth of 1m) cause results in inf
unique_lakeIDs <- unique_lakeIDs %>%
  filter(is.finite(TSI))

# Look at histogram of TSI to see spread
ggplot(unique_lakeIDs, aes(x = TSI)) +
  geom_histogram(binwidth = 2, color = "black", fill = "skyblue")+
  theme_minimal()

# Look at histogram of secchi depth to see spread
ggplot(unique_lakeIDs, aes(x = secchi)) +
  geom_histogram(binwidth = 2, color = "black", fill = "skyblue")+
  theme_minimal()

# Then assign oligotrophic as <40, and eutrophic as >50; based on carlson TSI classification

unique_lakeIDs <- unique_lakeIDs %>%
  mutate(trophic_state = case_when(
    TSI < 40 ~ "Oligotrophic",
    TSI > 50 ~ "Eutrophic",
    TRUE ~ "Mesotrophic" 
  ))

#see how many is in each column
table(unique_lakeIDs$trophic_state)

#Eutrophic  Mesotrophic Oligotrophic 
#67          61         28


#color histogram by trophic state
ggplot(unique_lakeIDs, aes(x = TSI, fill = trophic_state)) +
  geom_histogram(binwidth = 2, color = "black", position = "identity", alpha = 0.7) +
  labs(title = "Histogram of TSI by Trophic State",
       x = "Trophic State Index (TSI)",
       y = "Number of Lakes",
       fill = "Trophic State") +
  scale_fill_manual(values = c("Oligotrophic" = "skyblue", "Eutrophic" = "forestgreen")) +
  theme_minimal()

#color on map by trophic state
ggplot(unique_lakeIDs, aes(x = longitude, y = latitude, color = trophic_state)) +
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

#repair max depth issue found while calculating VHT; see mismatch below

# A tibble: 6 × 4
#link_lagosid date       temp_depth lake_maxdepth_m
#<int> <date>          <dbl>           <dbl>
#1          197 1984-08-14       6               5.5 
#2          305 2017-06-07       5               2.13
#3          305 2017-08-29       5               2.13
#4          381 2006-07-24       4.5             4.27
#5          734 2009-07-14       7.77            7.62
#6         1171 2013-07-16      22.9            10.4 

#this depth data is the same as in LAGOS; find lakes on internet by their coordinates and name

#fix these DEPTHS earlier up: lakes 197, 734, 375

#investigate profiles

lake_ids <- c(305, 305, 381, 1171)
dates <- c("2017-06-07", "2017-08-29", "2006-07-24", "2013-07-16")

# Filter to all four specific profiles
profiles_to_investigate <- MN_profiles[
  (MN_profiles$link_lagosid == 305 & MN_profiles$date == "2017-06-07") |
    (MN_profiles$link_lagosid == 305 & MN_profiles$date == "2017-08-29") |
    (MN_profiles$link_lagosid == 381 & MN_profiles$date == "2006-07-24") |
    (MN_profiles$link_lagosid == 1171 & MN_profiles$date == "2013-07-16"), ]

#In 381 2006-07-24, error is in depth_m, recorded as 4.5 but lets change it to 4.27
#in 1171 i think they recorded measurements in feet and not in meters

#lets just remove these 4 profile dates
MN_profiles <- MN_profiles[
  !((MN_profiles$link_lagosid == 305 & MN_profiles$date == "2017-06-07") |
      (MN_profiles$link_lagosid == 305 & MN_profiles$date == "2017-08-29") |
      (MN_profiles$link_lagosid == 381 & MN_profiles$date == "2006-07-24") |
      (MN_profiles$link_lagosid == 1171 & MN_profiles$date == "2013-07-16")), ]

#and in final dataset

MN_selected_profiles <- MN_selected_profiles[
  !((MN_selected_profiles$link_lagosid == 305 & MN_selected_profiles$date == "2017-06-07") |
      (MN_selected_profiles$link_lagosid == 305 & MN_selected_profiles$date == "2017-08-29") |
      (MN_selected_profiles$link_lagosid == 381 & MN_selected_profiles$date == "2006-07-24") |
      (MN_selected_profiles$link_lagosid == 1171 & MN_selected_profiles$date == "2013-07-16")), ]


#depth category -----------------------------

#use 33% quantiles as shallow and deep cutoffs
#this aligns with <10m or less being shallow for stratifying lakes

mean_depth <- mean(unique_lakeIDs$lake_maxdepth_m, na.rm = TRUE) #14.1
sd_depth <- sd(unique_lakeIDs$lake_maxdepth_m, na.rm = TRUE) #6.68

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
#53     51     52

#count categories; there is some crossover, but not entirely 
unique_lakeIDs %>%
  count(trophic_state, depth_category)

#Eutrophic     Deep              20
#Eutrophic     Medium            21
#Eutrophic     Shallow           26
#Mesotrophic   Deep              18
#Mesotrophic   Medium            20
#Mesotrophic   Shallow           23
#Oligotrophic  Deep              15
#Oligotrophic  Medium            10
#Oligotrophic  Shallow            3


#now for hydrology category, seepage (or "isolated" vs. drainage)

#from LAGOS dataset: 
#lake_connectivity_permanent = hydrologic connectivity class of the focal lake determined from the NHD network considering only permanent flow

#filted by only unique_lakeIDs lagos numbers
lake_characteristics <- lake_characteristics %>%
  filter(lagoslakeid %in% unique_lakeIDs$lagoslakeid)

#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_permanent)

#Drainage DrainageLk  Headwater   Isolated   Terminal TerminalLk 
#22        67          17           41        8          1


#I AM DECIDING TO ONLY USE PERMANENT FLOW.... NOT CONSIDERING BOTH PERMAMENT AND INTERMITTENT-EPHERMERAL FLOW -------

#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_fluctuates)

#N   Y 
#138 18

#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_class)

#Drainage DrainageLk  Headwater   Isolated   Terminal 
#21         77        16.          36         6 


#code definitions --
#drainage: lake with inlets and outlets; no upstream lakes greater than or equal to 10 ha
#drainagelk: lake with inlets and outlets; has upstream lakes greater than or equal to 10 ha
#headwater: lake with only outlet
#terminallk: lake with only inlet
#isolated: lake with neither an inlet or an outlet

#so i am going to combine drainage, drainagelk, headwater, and terminal into one category? to get seepage vs. drainage

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
    lake_connectivity_permanent == "Headwater" ~ "Drainage",
    lake_connectivity_permanent == "TerminalLk" ~ "Drainage",
    lake_connectivity_permanent == "Terminal" ~ "Drainage"
  ))

#see how many is in each column
table(unique_lakeIDs$hydro_category)

#Drainage  Seepage 
#115        41

ggplot(unique_lakeIDs, aes(x = hydro_category)) +
  geom_histogram(binwidth = 2, color = "black", fill = "skyblue", stat = "count") +
  geom_vline(xintercept = shallow_cut, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = deep_cut, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Hydrologic Connectivity",
       x = "Connectivity Type",
       y = "Number of Lakes") +
  theme_minimal()

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

#now assign north and south to all lakes, following advice from Gretchen
#follows ecoregion boundaries-- other Minnesota lit uses a squiggly line, not a line at 46 lat

#define points along drawn boundary
boundary_data <- data.frame(
  lon = c(-97, -96.5, -96, -95, -94.5, -94, -93.5, -93, -92, -91, -90),
  lat = c(47.2, 47.0, 46.8, 46.5, 46.3, 46.6, 46.4, 46.2, 45.0, 45.0, 45.0)
)

#fit a smooth function
boundary_model <- gam(lat ~ s(lon), data = boundary_data)

#apply to lakes
unique_lakeIDs <- unique_lakeIDs %>%
  mutate(
    boundary_lat = predict(boundary_model, newdata = data.frame(lon = longitude)),
    geo_position = if_else(latitude > boundary_lat, "North", "South")
  )

#see how many is in each column
table(unique_lakeIDs$geo_position)

#North South 
#36    120


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

print(p1)


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


#there are a lot of NaNs here... 12311
sum(is.nan(as.matrix(preseason_waterq_predictors)))

#see which columns they are mostly in: 
sapply(preseason_waterq_predictors, function(x) sum(is.nan(x)))

#chla_mayjune      chla_june     tp_mayjune        tp_june    secchi_mayjune    secchi_june   wtemp_mayjune  wtemp_june
#1377                1506           2854           2958             57            174           1639            1746


## TRIM down to only years for each lake that have profiles in MN_selected_profiles
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

#counts NaNs.. 12929
sum(is.nan(as.matrix(concurrent_waterq_predictors))) 

#see which columns they are mostly in: 
sapply(concurrent_waterq_predictors, function(x) sum(is.nan(x)))

#chla_JAS    chla_aug      tp_JAS      tp_aug   secchi_JAS   secchi_aug   wtemp_JAS   wtemp_aug 
#1389          1608        2924        3090          52         282         1697         1887

## TRIM down to only year for each lake that have profiles in MN_selectedprofiles
concurrent_waterq_predictors <- concurrent_waterq_predictors %>%
  semi_join(
    MN_selected_profiles %>% distinct(lagoslakeid = link_lagosid, year),
    by = c("lagoslakeid", "year")
  )


#ice data ----------------------------------------------------------------------------------------

#predictors needed for all years: ice-on time, water level
#use winter.year to filter years
ice_data <- ice_data %>%
  filter(winter.year >= 1980 & winter.year <= 2020)

#how many lakes are there with ice data: 1188 
print(paste("Unique lake IDs:", length(unique(ice_data$DOW))))

#ice duration processing; first trim down columns
#keep winter.year, DOW, min_duration
ice_data <- ice_data %>%
  select(DOW, winter.year, min_duration)

#pull data from full_MN_profiles and lagosne to try to find lagosid's for all DOW ice lakes
testing_bridge <- LAGOSNE_lakeslocus101 %>%
  inner_join(
    FULL_MN_profiles %>% distinct(nhdid = link_nhd_id, lagoslakeid = link_lagosid, state_ident, latitude, longitude),
    by = c("nhdid", "lagoslakeid")
  ) %>%
  select(lagoslakeid, nhdid, state_ident, latitude, longitude)

#now join using bridge
ice_data <- ice_data %>%
  left_join(
    testing_bridge %>% select(state_ident, lagoslakeid),
    by = c("DOW" = "state_ident")
  )

#how many DOWs are missing a lagos_id? #222
ice_data %>%
  filter(is.na(lagoslakeid)) %>%
  summarise(unique_DOWs_without_lagoslakeid = n_distinct(DOW))


#how many lakes match?
#there are 69 lakes in unique_lakeIDs that are missing from ice_data
missing_lakes <- setdiff(unique_lakeIDs$lagoslakeid, ice_data$lagoslakeid)
cat("Number of lakes missing from ice_data:", length(missing_lakes), "\n")

missing_lakes <- data.frame(
  lagos_lakeid = missing_lakes
)

#plot to see on graph
unique_lakeIDs$has_ice_data <- unique_lakeIDs$lagoslakeid %in% ice_data$lagoslakeid

ggplot(unique_lakeIDs, aes(x = longitude, y = latitude, color = has_ice_data)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"),
                     labels = c("TRUE" = "Has ice data", "FALSE" = "Missing ice data")) +
  coord_quickmap(
    xlim = c(-97.5, -89),   # Longitude range for Minnesota
    ylim = c(43.5, 49.5)    # Latitude range for Minnesota
  ) +
  theme_minimal() +
  labs(
    title = "Minnesota Lakes: Ice Data Coverage",
    x = "Longitude",
    y = "Latitude",
    color = "Status"
  )


# #try to match ids through lat and lon..
#this didnt work well

# missing_dow_coords <- FULL_MN_profiles %>%
#   filter(state_ident %in% missing_lagos_DOWs$DOW) %>%
#   select(state_ident, latitude, longitude) %>%  # adjust column names if needed
#   distinct(state_ident, .keep_all = TRUE)
# 
# dow_check <- FULL_MN_profiles %>%
#   filter(state_ident %in% missing_lagos_DOWs$DOW) %>%
#   distinct(state_ident)
# 
# cat("Number of missing DOWs found in FULL_MN_profiles:", nrow(dow_check), "\n")
# 
# # Check the coordinates for these DOWs
# missing_dow_coords <- FULL_MN_profiles %>%
#   filter(state_ident %in% missing_lagos_DOWs$DOW) %>%
#   select(state_ident, latitude, longitude) %>%
#   distinct(state_ident, .keep_all = TRUE)
# 
# cat("Number of missing DOWs with coordinates:", nrow(missing_dow_coords), "\n")
# 
# cat("\nNA values in FULL_MN_profiles coordinates:", 
#     sum(is.na(FULL_MN_profiles$latitude)) + sum(is.na(FULL_MN_profiles$longitude)), "\n")
# 
# cat("NA values in unique_lakeIDs coordinates:", 
#     sum(is.na(unique_lakeIDs$latitude)) + sum(is.na(unique_lakeIDs$longitude)), "\n")
# 
# missing_dow_coords <- FULL_MN_profiles %>%
#   filter(state_ident %in% missing_lagos_DOWs$DOW) %>%
#   filter(!is.na(latitude) & !is.na(longitude)) %>%  # Remove NAs
#   select(state_ident, latitude, longitude) %>%
#   distinct(state_ident, .keep_all = TRUE)
# 
# missing_dow_coords <- missing_dow_coords %>%
#   left_join(unique_lakeIDs %>% 
#               select(latitude, longitude, link_lagosid) %>% 
#               distinct(link_lagosid, .keep_all = TRUE), 
#             by = c("latitude" = "latitude", "longitude" = "longitude"))
# 
# # Round coordinates to reduce floating-point precision issues
# missing_dow_coords <- missing_dow_coords %>%
#   mutate(
#     lat_rounded = round(latitude, 2),    # 6 decimal places ≈ 0.1 meter precision
#     lon_rounded = round(longitude, 2)
#   )
# 
# unique_lakeIDs <- unique_lakeIDs %>%
#   mutate(
#     lat_rounded = round(latitude, 2),
#     lon_rounded = round(longitude, 2)
#   )
# 
# # Join on rounded coordinates
# result <- missing_dow_coords %>%
#   left_join(
#     unique_lakeIDs %>% 
#       select(lat_rounded, lon_rounded, link_lagosid) %>%
#       distinct(link_lagosid, .keep_all = TRUE),
#     by = c("lat_rounded" = "lat_rounded", "lon_rounded" = "lon_rounded")
#   )



# Continue on in R script 1.5-MN_ClimatePredictors, then 02-MN_VHT.R, then switch to jupyter hub notebooks

#sessionInfo()
#R version 4.3.1 (2023-06-16)
#Platform: x86_64-apple-darwin20 (64-bit)
#Running under: macOS 15.4.1
