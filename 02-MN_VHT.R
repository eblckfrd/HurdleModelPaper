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

#calculate VHT function -------------

#Minnesota: based on TDO3… so 3mg/L and 21.5 C
#from MDNR water quality standard guide

#first take a subset of lakes to apply it to and make sure this code is working accurately

first_five_lakes <- MN_selected_profiles %>%
  distinct(link_lagosid) %>%
  slice(1:5)

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
#MN_selected_profiles 

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

# Join the datasets to get geographic position for each lake
VHT_with_geo <- MN_selected_profiles %>%
  left_join(unique_lakeIDs, by = c("link_lagosid" = "lagoslakeid"))

## FIND OUT WHY 2019 IS BEING WEIRD - LOTS OF 0'S

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

# Get unique profiles with their VHT values from 2019
profiles_with_vht_2019 <- MN_selected_profiles %>%
  filter(year == 2019) %>%
  distinct(profile_ident, .keep_all = TRUE) %>%
  select(profile_ident, VHT, date, link_lagosid) %>%
  arrange(profile_ident)

print(profiles_with_vht_2019)

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
    avg_VHT_all = mean(VHT, na.rm = TRUE),  # Average VHT including zeros
    .groups = "drop"
  ) %>%
  filter(total_profiles >= 5)

ggplot() +
  geom_polygon(data = minnesota, aes(x = long, y = lat, group = group), 
               fill = "grey90", color = "white") +  # Lighter grey
  geom_point(data = habitat_summary, 
             aes(x = longitude, y = latitude, color = prop_with_habitat, size = total_profiles), 
             alpha = 0.8) +
  scale_color_gradient(low = "tomato", high = "blue2", name = "Proportion\nwith Habitat", 
                       labels = scales::percent) +
  scale_size_continuous(name = "Total\nProfiles", range = c(1, 4)) +
  labs(title = "Proportion of Profiles with Suitable Habitat (VHT > 0)") +
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

#calculate JAS


