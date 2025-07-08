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

----------------------------------------------

# Load in MN_profiles data from Gretchen shared google drive

#FULL_MN_profiles <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/MN_profiles.csv")

-------------------------------------------------------

# Load in LAGOS data files downloaded from edi; dont end up using all of these

#lakeslimno10873 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/lakeslimno10873.csv")
        #lake lat, lon, meandepth, maxdepth
#LAGOSNE_state_105 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/LAGOSNE_state_105.csv")
        #state ID
#LAGOSNE_lakeslocus101 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/LAGOSNE_lakeslocus101.csv")
         #lagos ID, name, lake area, lake perim, hu zones, elevation, stat and county zone id
#epi_waterquality10873 <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/epi_waterquality10873.csv")
          #secchi depth, other water quality params
#lake_characteristics <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/lake_characteristics.csv")
          #hydrology/lake connectivity

    
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

# Histogram of years included = most recent is 2020, goes back to ~1950 but mostly >1980
ggplot(data = MN_profiles_10yrs, aes(x = year)) +
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

## from MGLP code, this selects the latest profile of the year 

#EndofYear = MN_profiles_10yrs %>%
  #group_by(link_lagosid, year)%>%
  #mutate(Max = max(doy))%>%
  #filter(DOY == Max)%>%
  #ungroup()

##North and South investigation -------------------------------------------------



# Analysis on rainfall and temp Mesonet data to define north and south regions --------------------------
# To confirm ecoregion seperation

# NCEI 1991-2020 Climatological Average for Summer (JJA)
#https://mesonet.agron.iastate.edu/plotting/auto/?_wait=no&q=125&month=summer&sdate=0101&edate=1231&sector=state&src=ncei_climate91&state=MN&opt=both&var=avg_high&cmap=autumn&cmap_r=on&_r=t&dpi=100&_fmt=png

#Load in data
#MN_precip <- read.csv("~/Documents/R Repositories/Ch2_DataAnalysis/MN_precip.csv")

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

#see where these stations are; filter based on avg high
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

#see where these stations are; filter based on avg low
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

# now apply k-means clustering to all three variables to see if there is agreeance on ranges

#assign 2 ks so we can split N/S

cluster_data <- MN_precip %>%
  select(avg_high, avg_low, total_precip) 

scaled_data <- scale(cluster_data)

set.seed(123)  # for reproducibility
kmeans_result <- kmeans(scaled_data, centers = 2)  # adjust number of clusters as needed

# Add cluster assignments back to your original data
MN_precip$cluster <- factor(kmeans_result$cluster)

#look at results
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

#elbow method for validation: elbow bend is at 2! validation for seperating N?S

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


#am I using this data for my predictors? or just to cluster? 


## LAGOS DATA CLEANING -----------------------------------------------------

#First filter all datasets for MN only = state 14

#might not even need this dataset....
LAGOSNE_lakeslocus101 <- LAGOSNE_lakeslocus101 %>%
  filter(state_zoneid == "State_14")

#filter for only lagos lake IDs that we need and select columns needed

epi_waterquality10873 <- epi_waterquality10873 %>%
  filter(lagoslakeid %in% unique_lakeIDs$link_lagosid) %>% 
  dplyr::select(lagoslakeid, secchi, sampledate, eventida10873, programid,  chla, colora, doc, nh4, no2, no2no3, srp, tdn, tdp, tkn, tn, tp)

lakeslimno10873 <- lakeslimno10873 %>%
  filter(lagoslakeid %in% unique_lakeIDs$lagoslakeid)

##LAGOS data processing and category separation -------------------------------------------------------------

#using unique_lakeIDs as my sorting dataset
# First I need data so that i can categorize my lakes..

#find most recent secchi measurement for each lake in lagos data; 2012 is most recent.............

recent_secchi <- epi_waterquality10873 %>%
  filter(!is.na(secchi)) %>%
  group_by(lagoslakeid) %>%
  slice_max(order_by = sampledate, n = 1) %>%
  ungroup() %>%
  select(lagoslakeid, sampledate, secchi)

#rename link_lagosid to lagoslakeid so i can merge with LAGOS datasets

unique_lakeIDs <- unique_lakeIDs %>%
  rename(lagoslakeid = link_lagosid)

#now merge in secchi depth using lagoslakeid
unique_lakeIDs <- unique_lakeIDs %>%
  left_join(recent_secchi %>% select(lagoslakeid, secchi), by = "lagoslakeid")

# From secchi depth, calculate TSI for each lake: 60 - (14.41 * ln(secci #)) and put in new column

unique_lakeIDs <- unique_lakeIDs %>%
  mutate(TSI = 60 - (14.41 * log(secchi)))

#remove any nan or inf, takes us down to 500 lakes
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
#202          194          104


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

#even though mean depth might be a better depth separator, but only have max depth available through LAGOS and profiles dataset

#use 33% quantiles as shallow and deep cutoffs
#this aligns with <10m or less being shallow for stratifying lakes

mean_depth <- mean(unique_lakeIDs$lake_maxdepth_m, na.rm = TRUE) #14.6
sd_depth <- sd(unique_lakeIDs$lake_maxdepth_m, na.rm = TRUE) #9.6

cutoffs <- quantile(unique_lakeIDs$lake_maxdepth_m, probs = c(0.33, 0.66), na.rm = TRUE)

shallow_cut <- cutoffs[1]
deep_cut <- cutoffs[2]

# Classify lakes
unique_lakeIDs <- unique_lakeIDs %>%
  mutate(depth_category = case_when(
    lake_maxdepth_m < shallow_cut ~ "Shallow",
    lake_maxdepth_m > deep_cut ~ "Deep",
    TRUE ~ "Medium"  
  ))

#plot
ggplot(unique_lakeIDs, aes(x = lake_maxdepth_m)) +
  geom_histogram(binwidth = 2, color = "black", fill = "skyblue") +
  geom_vline(xintercept = shallow_cut, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = deep_cut, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Histogram of Lake Max Depth with Quantile-Based Cutoffs",
       x = "Max Depth (m)",
       y = "Number of Lakes") +
  theme_minimal()

#see how many is in each column
table(unique_lakeIDs$depth_category)

#  Deep  Medium Shallow 
#170     165     165 

#count categories; there is some crossover, but not entirely split up
unique_lakeIDs %>%
  count(trophic_state, depth_category)

#1 Eutrophic     Deep             21
# Eutrophic     Medium            69
# Eutrophic     Shallow          112

# Oligotrophic  Deep              62
# Oligotrophic  Medium            29
# Oligotrophic  Shallow           13


#now for hydrology category, seepage (or "isolated" vs. drainage)

#from LAGOS dataset: 
#lake_connectivity_permanent = hydrologic connectivity class of the focal lake determined from the NHD network considering only permanent flow

#filted by only unique_lakeIDs lagos numbers
lake_characteristics <- lake_characteristics %>%
  filter(lagoslakeid %in% unique_lakeIDs$lagoslakeid)

#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_permanent)

#Drainage DrainageLk  Headwater   Isolated   Terminal TerminalLk 
#62        213         62        144         13          6


#I AM DECIDING TO ONLY USE PERMANENT FLOW.... NOT CONSIDERING BOTH PERMAMENT AND INTERMITTENT-EPHERMERAL FLOW --------

#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_fluctuates)

#N   Y 
#430  70


#see how many is in each hydro column
table(lake_characteristics$lake_connectivity_class)

#Drainage DrainageLk  Headwater   Isolated   Terminal TerminalLk 
#74        245         47        121         10          3 


#code definitions ------
#drainage: lake with inlets and outlets; no upstream lakes greater than or equal to 10 ha
#drainagelk: lake with inlets and outlets; has upstream lakes greater than or equal to 10 ha
#isolated: lake with neither an inlet or an outlet

#so i am going to combine drainage and drainagelk into one category? to get seepage vs. drainage

#merge in connectivity column

#now merge  using lagoslakeid
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
#275       81      144 

#now assign north and south based on "K-means Clustering of Mesonet Stations" figure, 46 LAT looks like good divider

unique_lakeIDs <- unique_lakeIDs %>%
  mutate(geo_position = if_else(latitude > 46, "North", "South"))

#see how many is in each column
table(unique_lakeIDs$geo_position)

#North South 
#186   313

#map this to see if it seperated well

ggplot(unique_lakeIDs) +
  aes(x = longitude, y = latitude, color = geo_position) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(size = 2, alpha = 0.8) +
  coord_quickmap(
    xlim = c(-97.5, -89),
    ylim = c(43.5, 49.5)
  ) +
  scale_color_manual(
    name = "Trophic State",
    values = c("North" = "skyblue", "South" = "red3")
  ) +
  theme_minimal() +
  labs(
    title = "Minnesota Lakes Split North/South",
    x = "Longitude",
    y = "Latitude"
  )



