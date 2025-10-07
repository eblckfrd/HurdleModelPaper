# Dissertation Chapter 2 - MINNESOTA - July 2025
# Calculating and visualizing VHT

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


# Calculate VHT for 5 randomly selected lake profiles ----------------------------------------------------------------------------------------

#calculate VHT function 

#Minnesota: based on TDO3… so 3mg/L and 21.5 C
#from MDNR water quality standard guide

#first take a subset of lakes to apply it to and make sure this code is working accurately

first_five_lakes <- MN_selected_profiles %>%
  distinct(link_lagosid) %>%
  slice(200:205)

firstfive_lakedata <- MN_selected_profiles %>%
  filter(link_lagosid %in% first_five_lakes$link_lagosid) %>%
  dplyr::select(link_lagosid, date, month, year, depth_m, temp_c, do_ppm, lake_maxdepth_m, longitude, latitude)

#see where they are located

ggplot(firstfive_lakedata, aes(x = longitude, y = latitude)) +
  borders("state", regions = "minnesota", fill = "gray95", colour = "gray70") +
  geom_point(color = "blue", size = 1, alpha = 0.7) +
  coord_quickmap(
    xlim = c(-97.5, -89),   # Longitude range for Minnesota
    ylim = c(43.5, 49.5)    # Latitude range for Minnesota
  ) +
  theme_minimal() +
  labs(
    title = "Selected Lakes for VHT Test",
    x = "Longitude",
    y = "Latitude"
  )

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
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(linewidth = 2.5, alpha = 0.9) +
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
#Use MN_selected_profiles dataset

#check for any missing O2 or temp data

missing_data_summary <- MN_selected_profiles %>%
  group_by(link_lagosid, date) %>%
  summarise(
    has_temp = any(!is.na(temp_c)),
    has_o2 = any(!is.na(do_ppm)),
    .groups = "drop"
  ) %>%
  mutate(
    missing_category = case_when(
      !has_temp & !has_o2 ~ "Missing both",
      !has_temp ~ "Missing temp only",
      !has_o2 ~ "Missing O2 only",
      TRUE ~ "Has both"
    )
  )

# Check the overlap
table(missing_data_summary$missing_category)

#Has both   Missing O2 only Missing temp only 
#25908              4462               990


#Switching here to using MN_VHT_calculator going forward

#filter to only keep only profiles with BOTH temp and O2
MN_VHT_calculator <- MN_VHT_calculator %>%
  group_by(link_lagosid, date) %>%
  filter(
    any(!is.na(temp_c)) & any(!is.na(do_ppm))
  ) %>%
  ungroup()

#now calculate VHT
MN_VHT_calculator <- MN_VHT_calculator %>%
  group_by(link_lagosid, date, month, year) %>%
  arrange(link_lagosid, date, depth_m) %>%
  
  # Calculate temp_depth: first depth where temp_c <= 22
  mutate(
    temp_depth = {
      if (all(is.na(temp_c))) {
        NA_real_
      } else {
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
summary(MN_VHT_calculator[c("temp_depth", "o2_depth", "VHT")])

#temp_depth        o2_depth           VHT        
#Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
#1st Qu.: 0.000   1st Qu.: 5.000   1st Qu.: 0.000  
#Median : 0.500   Median : 7.000   Median : 5.000  
#Mean   : 2.217   Mean   : 8.672   Mean   : 6.648  
#3rd Qu.: 4.064   3rd Qu.:10.000   3rd Qu.: 9.000  
#Max.   :30.000   Max.   :56.038   Max.   :56.038 

#plot
MN_VHT_calculator %>%
  group_by(link_lagosid, date) %>%
  slice(1) %>%
  ungroup() %>%
  select(temp_depth, o2_depth) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "depth") %>%
  ggplot(aes(x = depth, fill = metric)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  theme_minimal() +
  labs(title = "Distribution of temp_depth vs o2_depth",
       x = "Depth (m)", y = "Count")

##QC CHECKS ------------------------------------------------------------------------------

# Check for any remaining issues
cat("Data integrity checks:\n")
cat("Any negative temp_depth:", any(MN_VHT_calculator$temp_depth < 0, na.rm = TRUE), "\n")
cat("Any negative o2_depth:", any(MN_VHT_calculator$o2_depth < 0, na.rm = TRUE), "\n")
cat("Any negative VHT:", any(MN_VHT_calculator$VHT < 0, na.rm = TRUE), "\n")
cat("Any NA values in temp_depth:", sum(is.na(MN_VHT_calculator$temp_depth)), "\n")
cat("Any NA values in o2_depth:", sum(is.na(MN_VHT_calculator$o2_depth)), "\n")
cat("Any NA values in VHT:", sum(is.na(MN_VHT_calculator$VHT)), "\n")

cat("Validation checks:\n")
cat("Profiles with NA oxygen data:", sum(is.na(MN_VHT_calculator$o2_depth)), "\n")
cat("Profiles with NA VHT:", sum(is.na(MN_VHT_calculator$VHT)), "\n")
cat("NAs match (should be TRUE):", sum(is.na(MN_VHT_calculator$o2_depth)) == sum(is.na(MN_VHT_calculator$VHT)), "\n\n")

logic_check <- MN_VHT_calculator %>%
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
case_summary <- MN_VHT_calculator %>%
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

#No suitable habitat (VHT = 0)      7951
#Suitable habitat exists (VHT > 0) 17957

# Check that VHT = o2_depth - temp_depth (should be TRUE for all rows)
VHT_check <- MN_VHT_calculator %>%
  mutate(
    calculated_VHT = pmax(0, o2_depth - temp_depth),
    VHT_matches = abs(VHT - calculated_VHT) < 0.001  # Allow for small floating point errors
  )

cat("VHT calculation consistency:", all(VHT_check$VHT_matches), "\n")
cat("Rows with VHT inconsistency:", sum(!VHT_check$VHT_matches), "\n")

# Look at a few random profiles to manually verify logic
set.seed(123)
sample_profiles <- MN_VHT_calculator %>%
  select(link_lagosid, date) %>%
  distinct() %>%
  sample_n(5)

# Check these profiles in detail
for(i in 1:nrow(sample_profiles)) {
  profile_data <- MN_VHT_calculator %>%
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
weird_temp <- MN_VHT_calculator %>%
  filter(temp_depth > lake_maxdepth_m) %>%
  select(link_lagosid, date, temp_depth, lake_maxdepth_m) %>%
  distinct()
print(head(weird_temp))


cat("Profiles where o2_depth > max lake depth:\n")
weird_o2 <- MN_VHT_calculator %>%
  filter(o2_depth > lake_maxdepth_m) %>%
  select(link_lagosid, date, o2_depth, lake_maxdepth_m) %>%
  distinct()
print(head(weird_o2))

# Check for extreme VHT values
cat("Profiles with VHT > 65m (potential outliers):\n")
extreme_vht <- MN_VHT_calculator %>%
  filter(VHT > 65) %>%
  select(link_lagosid, date, temp_depth, o2_depth, VHT, lake_maxdepth_m) %>%
  distinct()
print(extreme_vht)

# Check for unrealistic values
cat("\nDistribution checks:\n")
cat("Profiles where temp_depth > max lake depth:\n")
weird_temp <- MN_VHT_calculator %>%
  filter(temp_depth > lake_maxdepth_m) %>%
  select(link_lagosid, date, temp_depth, lake_maxdepth_m) %>%
  distinct()
print(head(weird_temp))

cat("Profiles where o2_depth > max lake depth:\n")
weird_o2 <- MN_VHT_calculator %>%
  filter(o2_depth > lake_maxdepth_m) %>%
  select(link_lagosid, date, o2_depth, lake_maxdepth_m) %>%
  distinct()
print(head(weird_o2))

# Check that each profile has the same VHT values across all depths
profile_consistency <- MN_VHT_calculator %>%
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
p1 <- ggplot(MN_VHT_calculator %>% distinct(link_lagosid, date, .keep_all = TRUE), 
             aes(x = VHT)) +
  geom_histogram(bins = 50) +
  labs(title = "Distribution of VHT values")

print(p1)

plot_data <- MN_VHT_calculator %>% 
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
p1 <- MN_VHT_calculator %>%
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  ggplot(aes(x = factor(month), y = VHT)) +
  geom_boxplot() +
  labs(title = "VHT by Month", 
       x = "Month", y = "VHT (m)") +
  theme_minimal()

print(p1)

# Year-over-year trends

# Join the datasets to get geographic position for each lake
VHT_with_geo <- MN_VHT_calculator %>%
  left_join(unique_lakeIDs, by = c("link_lagosid" = "lagoslakeid"))


# Year-over-year trends by geographic position
p2 <- VHT_with_geo %>%
  filter(!is.na(VHT), !is.na(geo_position)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  group_by(year, geo_position) %>%
  summarise(mean_VHT = mean(VHT, na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = year, y = mean_VHT, color = geo_position)) +
  geom_line(size = 1) +
  labs(title = "VHT Trends: Northern vs Southern Lakes", 
       x = "Year", 
       y = "Mean VHT (m)", 
       color = "Geographic Position") +
  theme_minimal() +
  scale_color_manual(values = c("North" = "#2E86C1", "South" = "#E74C3C"))

print(p2)

p3 <- MN_VHT_calculator %>%
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  ggplot(aes(x = lake_maxdepth_m, y = VHT)) +
  geom_point(alpha = 0.5) +
#  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "VHT vs Lake Maximum Depth", 
       x = "Lake Maximum Depth (m)", y = "VHT (m)") +
  theme_minimal()

print(p3)

p8 <- MN_VHT_calculator %>%
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
seasonal_summary <- MN_VHT_calculator %>%
  filter(!is.na(VHT)) %>%
  mutate(
    season = case_when(
      month %in% c(3, 4, 5, 6) ~ "Spring", 
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

# Map 2: Average VHT at each location with capped scale
habitat_summary_capped <- habitat_summary %>%
  mutate(avg_VHT_capped = pmin(avg_VHT_all, 30))  # Cap values at 30m

ggplot() +
  geom_polygon(data = minnesota, aes(x = long, y = lat, group = group), 
               fill = "grey90", color = "white") +  # Lighter grey
  geom_point(data = habitat_summary_capped, 
             aes(x = longitude, y = latitude, color = avg_VHT_capped), 
             alpha = 0.8, size = 2) +
  scale_color_viridis_c(name = "Average\nVHT (m)", 
                        limits = c(0, 30),
                        breaks = c(0, 10, 20, 30),
                        labels = c("0", "10", "20", "30+")) +
  labs(title = "Average VHT at Each Lake Location (Capped at 30m)") +
  coord_fixed(1.3) +
  theme_void()


# proceed -----------------

#simplify
MN_VHT_summary <- MN_VHT_calculator %>%
  group_by(link_lagosid, date, month, year) %>%
  summarise(
    temp_depth = first(temp_depth),
    o2_depth = first(o2_depth),
    VHT = first(VHT),
    max_depth_sampled = max(depth_m, na.rm = TRUE),
    lake_maxdepth = first(lake_maxdepth_m),
    .groups = "drop"
  )

#calculate August, and JAS

VHT_JAS <- MN_VHT_summary %>%
  # Filter for July, August, September only
  filter(month %in% c(7, 8, 9)) %>%
  group_by(link_lagosid, year) %>%
  summarise(
    # July/August/September average
    JAS_VHT = mean(VHT, na.rm = TRUE),
    
    .groups = "drop"
  )

summary(VHT_JAS)

#histogram
ggplot(concurrent_VHT, aes(x = JAS_VHT)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal()

#look at map plot

lake_VHT_map <- MN_VHT_calculator %>%
  select(link_lagosid, latitude, longitude) %>%
  distinct() %>%
  left_join(
    VHT_JAS %>%
      group_by(link_lagosid) %>%
      summarise(mean_JAS_VHT = mean(JAS_VHT, na.rm = TRUE), .groups = "drop"),
    by = "link_lagosid"
  )

ggplot() +
  geom_polygon(data = minnesota, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  # Lake points colored by VHT
  geom_point(data = lake_VHT_map, 
             aes(x = longitude, y = latitude, color = mean_JAS_VHT),
             size = 2, alpha = 0.7) +
  scale_color_viridis_c(name = "JAS VHT (m)", option = "plasma") +
  theme_minimal() +
  coord_fixed(1.3) +
  labs(
    title = "July-August-September VHT across Minnesota Lakes",
    x = "Longitude",
    y = "Latitude"
  )

mean_VHT_by_year <- concurrent_VHT %>%
  group_by(year) %>%
  summarise(mean_VHT = mean(JAS_VHT, na.rm = TRUE), .groups = "drop")

# Create the plot
ggplot() +
  # Individual lake lines (thin and transparent)
  geom_line(data = concurrent_VHT, 
            aes(x = year, y = JAS_VHT, group = link_lagosid),
            alpha = 0.2, color = "gray50") +
  # Mean line (thick and prominent)
  geom_line(data = mean_VHT_by_year, 
            aes(x = year, y = mean_VHT),
            color = "red", size = 1.5) +
  geom_point(data = mean_VHT_by_year, 
             aes(x = year, y = mean_VHT),
             color = "red", size = 3) +
  theme_minimal() +
  labs(
    title = "July-August-September VHT Over Time",
    subtitle = "Gray lines: individual lakes | Red line: mean across all lakes",
    x = "Year",
    y = "JAS VHT (m)"
  )

#need to get 25 random lakes and put into csv to transfer over to python

#use VHT_JAS

#new dataset called lake_VHT, with lagoslakeid, VHT, latitude, longitude, and year

lake_VHT <- VHT_JAS %>%
  left_join(unique_lakeIDs, by = c("link_lagosid" = "lagoslakeid")) %>%
  select(
    lagoslakeid = link_lagosid,
    VHT = JAS_VHT,
    latitude = latitude,
    longitude = longitude,
    year = year
  )

write.csv(lake_VHT, file = "~/Desktop/lake_VHT.csv", row.names = FALSE)
