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
library(hexbin)
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

# visualizations --------------------------------------------------------------------

#seasonality in bar plots, truncuated at 25m

MN_VHT_calculator %>%
  filter(!is.na(VHT)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  ggplot(aes(x = factor(month), y = VHT, fill = factor(month))) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
  scale_fill_viridis_d(option = "viridis") +
  scale_x_discrete(labels = month.abb) +
  coord_cartesian(ylim = c(0, 25)) +  # Focus on 0-25m range
  labs(title = "Seasonal Patterns in VHT", 
       subtitle = "Y-axis truncated at 25m for clarity",
       x = "Month", 
       y = "VHT (m)") +
  theme_minimal() +
  theme(legend.position = "none")

#same plot but in range of blues with jitter and transparency
MN_VHT_calculator %>%
  filter(!is.na(VHT), VHT <= 25) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  ggplot(aes(x = factor(month), y = VHT)) +
  geom_jitter(aes(color = factor(month)), alpha = 0.3, width = 0.2, size = 0.8) +
  geom_boxplot(aes(fill = factor(month)), alpha = 0.6, outlier.shape = NA, width = 0.5) +
  scale_fill_manual(values = colorRampPalette(c("lightblue", "steelblue", "darkblue"))(12)) +
  scale_color_manual(values = colorRampPalette(c("lightblue", "steelblue", "darkblue"))(12)) +
  scale_x_discrete(labels = month.abb) +
  labs(title = "Seasonal Patterns in VHT", 
       x = "Month", 
       y = "VHT (m)") +
  theme_minimal() +
  theme(legend.position = "none")


#add all categories into dataset I am using to plot
MN_VHT_calculator <- MN_VHT_calculator %>%
  left_join(unique_lakeIDs %>% 
            select(lagoslakeid, trophic_state, depth_category, hydro_category, geo_position),
            by = c("link_lagosid" = "lagoslakeid"))


#Now can plot seasonality differences based on categories ------------------------------------------------

#Create loops to make plots for all the categories

# Define categories to loop through
categories <- list(
  list(var = "hydro_category", values = c("Drainage", "Seepage")),
  list(var = "trophic_state", values = c("Eutrophic", "Oligotrophic")),
  list(var = "depth_category", values = c("Deep", "Shallow")),
  list(var = "geo_position", values = c("South", "North"))
)

# Define categories to loop through
categories <- list(
  list(var = "hydro_category", values = c("Drainage", "Seepage"), title = "Hydrology"),
  list(var = "trophic_state", values = c("Eutrophic", "Oligotrophic"), title = "Trophic State"),
  list(var = "depth_category", values = c("Deep", "Shallow"), title = "Mean Depth"),
  list(var = "geo_position", values = c("South", "North"), title = "Geographic Position")
)

# Loop through each category for line plots by year
for (cat in categories) {
  var_name <- cat$var
  values <- cat$values
  custom_title <- cat$title
  
  # Prepare the data first
  plot_data <- MN_VHT_calculator %>%
    filter(!is.na(VHT), VHT <= 25) %>%
    filter(!!sym(var_name) %in% values) %>%
    distinct(link_lagosid, date, .keep_all = TRUE) %>%
    group_by(!!sym(var_name), year, month) %>%
    summarise(mean_VHT = mean(VHT, na.rm = TRUE),
              se_VHT = sd(VHT, na.rm = TRUE) / sqrt(n()), .groups = "drop")
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = month, y = mean_VHT, color = factor(year), group = year)) +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 2) +
    scale_color_viridis_d(option = "plasma") +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Keep month abbreviations
    facet_wrap(as.formula(paste("~", var_name))) +
    labs(title = paste("Seasonal Patterns in VHT Across Years by", custom_title),
         x = NULL,  # Remove "Month" label
         y = "Mean VHT (m)") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      legend.position = "none",  # Remove year legend
      strip.text = element_text(size = 14, face = "bold"),  # Increase facet label size
      axis.text = element_text(size = 14),  # Increase axis text size (includes months)
      axis.text.x = element_text(size = 14),  # Specifically increase month labels
      axis.title.y = element_text(size = 14, face = "bold")  # Make y-axis title bigger and bold
    )
  
  print(p)
  
  # Save to desktop
  #ggsave(filename = paste0("~/Desktop/", custom_title, "_Seasonal_Patterns.png"), 
  #       plot = p, width = 10, height = 6, dpi = 300)
}

# Loop through each category for ribbon plots
for (cat in categories) {
  var_name <- cat$var
  values <- cat$values
  custom_title <- cat$title
  
  # Prepare the data first
  plot_data <- MN_VHT_calculator %>%
    filter(!is.na(VHT), VHT <= 25) %>%
    filter(!!sym(var_name) %in% values) %>%
    distinct(link_lagosid, date, .keep_all = TRUE) %>%
    group_by(!!sym(var_name), month) %>%
    summarise(mean_VHT = mean(VHT, na.rm = TRUE),
              sd_VHT = sd(VHT, na.rm = TRUE), .groups = "drop")
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = month)) +
    geom_ribbon(aes(ymin = mean_VHT - sd_VHT, ymax = mean_VHT + sd_VHT, 
                    fill = !!sym(var_name)), 
                alpha = 0.3) +
    geom_line(aes(y = mean_VHT, color = !!sym(var_name)), linewidth = 1.5) +
    scale_color_manual(values = c("steelblue", "darkblue")) +
    scale_fill_manual(values = c("steelblue", "darkblue")) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Keep month abbreviations
    labs(title = paste("Seasonal VHT Patterns by", custom_title),
         subtitle = "Line shows mean, ribbon shows ±1 SD",
         x = NULL,  # Remove "Month" label
         y = "Mean VHT (m)", 
         color = custom_title, 
         fill = custom_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      strip.text = element_text(size = 14, face = "bold"),  # Increase facet label size
      axis.text = element_text(size = 14),  # Increase axis text size
      axis.text.x = element_text(size = 14),  # Specifically increase month labels
      axis.title.y = element_text(size = 14, face = "bold"),  # Make y-axis title bigger and bold
      legend.key.size = unit(1.5, "cm"),  # Adjusts the size of the legend key (box)
      legend.text = element_text(size = 14),  # Increase legend text
      legend.title = element_text(size = 12, face = "bold")  # Increase legend title
    )
  
  print(p)
  
  # Save to desktop
  #ggsave(filename = paste0("~/Desktop/", custom_title, "_Ribbon_Plot.png"), 
  #       plot = p, width = 10, height = 6, dpi = 300)
}

# Year-over-year trends by category




# Year-over-year trend by geographic position ------------------------------------------------------------------

# Join the datasets to get geographic position for each lake
VHT_with_geo <- MN_VHT_calculator %>%
  left_join(unique_lakeIDs, by = c("link_lagosid" = "lagoslakeid"))


# Year-over-year trends by geographic position
VHT_with_geo %>%
  filter(!is.na(VHT), !is.na(geo_position)) %>%
  distinct(link_lagosid, date, .keep_all = TRUE) %>%
  group_by(year, geo_position) %>%
  summarise(mean_VHT = mean(VHT, na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = year, y = mean_VHT, color = geo_position)) +
  geom_line(size = 1) +
  labs(title = "VHT Trends: Northern vs Southern Lakes", 
       x = NULL,
       y = "Mean VHT (m)", 
       color = "Geographic Position") +
  theme_minimal() +
  scale_color_manual(values = c("North" = "#2E86C1", "South" = "#E74C3C")) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label size
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.text.x = element_text(size = 14),  # Specifically increase month labels
    axis.title.y = element_text(size = 14, face = "bold"),  # Make y-axis title bigger and bold
    legend.key.size = unit(1.5, "cm"),  # Adjusts the size of the legend key (box)
    legend.text = element_text(size = 14),  # Increase legend text
    legend.title = element_text(size = 14, face = "bold")  # Increase legend title
  )

#plots on map ---------------------------------------------------------

library(maps)

# Get Minnesota map
minnesota <- map_data("state", region = "minnesota")
mn_map <- map_data("state", region = "minnesota")

# Seasonal VHT patterns
seasonal_summary <- MN_VHT_calculator %>%
  filter(!is.na(VHT)) %>%
  mutate(
    season = case_when(
      month %in% c(5, 6) ~ "Preseason", 
      month %in% c(7, 8, 9) ~ "Concurrent"
    )
  ) %>%
  group_by(link_lagosid, latitude, longitude, season) %>%
  summarise(
    avg_VHT = mean(VHT, na.rm = TRUE),
    n_profiles = n(),
    .groups = "drop"
  ) %>%
  filter(n_profiles >= 4)  # Only lakes with multiple profiles per season

#create table showing seasonal changes
seasonal_change <- seasonal_summary %>%
  filter(!is.na(season)) %>%  # Remove rows where season is NA
  pivot_wider(
    id_cols = c(link_lagosid, latitude, longitude),
    names_from = season, 
    values_from = avg_VHT,
    values_fn = mean
  ) %>%
  filter(!is.na(Preseason) & !is.na(Concurrent)) %>%
  mutate(VHT_change = Concurrent - Preseason)

#summarize
summary(seasonal_change$VHT_change)

#plot changes on map
ggplot() +
  geom_polygon(data = mn_map, aes(x = long, y = lat, group = group), 
               fill = "gray95", color = "gray50") +
  geom_point(data = seasonal_change, 
             aes(x = longitude, y = latitude, color = VHT_change),
             size = 3, alpha = 0.8) +
  scale_color_viridis_c(
    option = "plasma",  # Try "plasma", "inferno", "magma", or "turbo"
    limits = c(-5, 5),
    oob = scales::squish,
    name = "VHT Change (m)"
  ) +
  labs(title = "Seasonal Change in VHT Across Minnesota Lakes") +
  coord_map() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )


#show using hexbin to avoid dots overlapping
ggplot() +
  geom_polygon(data = mn_map, aes(x = long, y = lat, group = group), 
               fill = "gray95", color = "gray50") +
  stat_summary_hex(data = seasonal_change,
                   aes(x = longitude, y = latitude, z = VHT_change),
                   bins = 30, alpha = 0.8) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = 1,
    limits = c(-5, 5),  # Symmetric around 0
    oob = scales::squish,
    name = "VHT Change (m)"
  ) +
  labs(title = "Regional Patterns in Seasonal VHT Change") +
  coord_map() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank()
  )

#investigate lakes that are showing increases...
#how did these lakes make it this far?! they do not have 10 different years of data

# Lakes with increases (VHT_change >= 0)
increasing_lakes <- seasonal_change %>%
  filter(VHT_change >= 0) %>%
  arrange(desc(VHT_change))

# View them
View(increasing_lakes)

increasing_vht <- increasing_lakes$link_lagosid

MN_VHT_calculator %>%
  filter(link_lagosid %in% increasing_vht,
         month %in% c(5, 6, 7, 8, 9)) %>%
  arrange(link_lagosid, date) %>%
  select(link_lagosid, date, month, VHT, year) %>%
  View()

MN_VHT_calculator %>%
  filter(link_lagosid %in% increasing_vht,
         month %in% c(5, 6, 7, 8, 9)) %>%
  mutate(season = case_when(
    month %in% c(5, 6) ~ "Preseason",
    month %in% c(7, 8, 9) ~ "Concurrent"
  )) %>%
  ggplot(aes(x = date, y = VHT, color = season)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(alpha = 0.3) +
  scale_color_manual(values = c("Preseason" = "lightblue", "Concurrent" = "darkblue")) +
  facet_wrap(~link_lagosid, scales = "free") +
  labs(title = "VHT Over Time for 'Increasing' Lakes",
       x = "Date", y = "VHT (m)", color = "Season") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

MN_VHT_calculator %>%
  filter(link_lagosid %in% increasing_vht,
         month %in% c(5, 6, 7, 8, 9)) %>%
  mutate(season = case_when(
    month %in% c(5, 6) ~ "Preseason",
    month %in% c(7, 8, 9) ~ "Concurrent"
  )) %>%
  group_by(link_lagosid, season) %>%
  summarise(
    n_measurements = n(),
    min_VHT = min(VHT),
    median_VHT = median(VHT),
    mean_VHT = mean(VHT),
    max_VHT = max(VHT),
    years_sampled = paste(unique(year), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(link_lagosid, season) %>%
  View()

data_balance <- seasonal_summary %>%
  filter(!is.na(season)) %>%
  select(link_lagosid, season, n_profiles)

balance_comparison <- seasonal_summary %>%
  filter(!is.na(season)) %>%
  select(link_lagosid, season, n_profiles) %>%
  pivot_wider(names_from = season, values_from = n_profiles, values_fill = 0)

ggplot(balance_comparison, aes(x = Preseason, y = Concurrent)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 1:1 line
  labs(title = "Sampling Balance: Preseason vs Concurrent Profiles",
       subtitle = "Points on red line = equal sampling",
       x = "Number of Preseason Profiles",
       y = "Number of Concurrent Profiles") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

balance_comparison <- balance_comparison %>%
  mutate(ratio = pmax(Preseason, Concurrent) / pmin(Preseason, Concurrent))

ggplot(balance_comparison, aes(x = ratio)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Sampling Balance Ratio Distribution",
       subtitle = "Red line at 2:1 ratio threshold",
       x = "Ratio (Max profiles / Min profiles)",
       y = "Number of Lakes") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

# Check the increasing lakes specifically
increasing_lake_ids <- increasing_lakes$link_lagosid

seasonal_summary %>%
  filter(link_lagosid %in% increasing_lake_ids) %>%
  select(link_lagosid, season, n_profiles, avg_VHT) %>%
  arrange(link_lagosid, season)


balance_comparison <- balance_comparison %>%
  mutate(increasing = link_lagosid %in% increasing_lakes$link_lagosid)

ggplot(balance_comparison, aes(x = Preseason, y = Concurrent)) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
  geom_point(data = balance_comparison %>% filter(!increasing), 
             aes(color = increasing), alpha = 0.6, size = 3) +
  geom_point(data = balance_comparison %>% filter(increasing), 
             aes(color = increasing), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "steelblue"),
                     labels = c("Decreasing/Stable", "Increasing")) +
  labs(title = "Are Increasing Lakes Due to Unbalanced Sampling?",
       x = "Number of Preseason Profiles",
       y = "Number of Concurrent Profiles",
       color = "VHT Trend") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12))

##as a result of this, i am going back up in code to only use lakes in plotting that have at least 4 samples from each season


# simplify VHT and prepare for model --------------------------------------------------------------

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

#write.csv(lake_VHT, file = "~/Desktop/lake_VHT.csv", row.names = FALSE)

#dataset of merged VHT and pre-season predictors

# Make sure join keys have the same type
lake_VHT <- lake_VHT %>%
  mutate(
    lagoslakeid = as.character(lagoslakeid),
    year        = as.integer(year)
  )

preseason_predictors <- preseason_predictors %>%
  mutate(
    lagoslakeid = as.character(lagoslakeid),
    year        = as.integer(year)
  )

Model_insert <- inner_join(
  lake_VHT,
  preseason_predictors,
  by = c("lagoslakeid", "year")
)

#checks
setdiff(unique_lakeIDs$lagoslakeid, Model_insert$lagoslakeid)

length(unique(Model_insert$lagoslakeid))

#save
#write.csv(Model_insert, file = "~/Desktop/Model_insert.csv", row.names = FALSE)

