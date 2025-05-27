---
title: "Data exploration practical"
excerpt: "Metagenome-associated environmental data exploration"
layout: single
mathjax: true
author: "Paida Mataranyika"
---

By Paida Mataranyika {% include toc %}

# Download the data:
Click here to download the data: [Data File](../../assets/environmental_metadata.csv)

# Exploring Environmental Data: Relationships and Visualisation

Take a look at the environmental data collected. Could we derive relationships or interactions based on the data? To do this, we will begin by loading the necessary packages.

### Load Required Packages

``` r

# Install Bioconductor if not already installed and use it to install genomic annotation packages
install.packages("Bioconductor")

# Install and load 'tidyverse', a collection of R packages for data manipulation and visualization
if (!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

# Install and load 'devtools', used for installing packages from GitHub
install.packages("devtools")
library("devtools")

# Install and load 'conflicted' to manage function name conflicts
devtools::install_github("r-lib/conflicted")
install.packages("gtable")  # For manipulating 'ggplot2' graphic objects

# Load additional libraries
library("gtable")
library(dplyr)
library(conflicted)

# Install and load clustering tools
if (!require(cluster)) install.packages("cluster")
library(cluster)
```

### Load and Inspect Environmental Data

Loading and inspecting the data is a critical first step. It allows us to ensure data is loaded correctly, to filter unnecessary and missing data, to understand the structure and data types, and to summarize the data.

``` r
# Read data from CSV and filter out rows with empty site names
data <- read.csv("../data/environmental_metadata.csv")
refined_metadata <- data %>% filter(site_name != "")

# Display structure of the dataset
str(refined_metadata)

# Count total missing values in the dataset
sum(is.na(refined_metadata))

# Provide summary statistics for each variable
summary(refined_metadata)

# Summarise categorical data into frequency tables
refined_metadata %>%
  summarise(
    crop_rotation_primary_y0_Count = table("crop_rotation_primary_y0"),
    takeall_seen_Count = table("takeall_seen"),
    fertiliser_use_Count = table("fertiliser_use"),
    carbon_Count = table("carbon")
  )
```

### Visualise Crop Rotation Distribution

We visualise crop rotation distribution in the script to gain a clear and immediate understanding of how crop rotation practices are spread across the data set. This section predominantly sets out to understand or paint a picture (literally) to understand category prevalence.

``` r
# Create a bar plot showing the distribution of 'crop_rotation_primary_y0'
refined_metadata %>%
  count(crop_rotation_primary_y0) %>%
  ggplot(aes(x = crop_rotation_primary_y0, y = n, fill = crop_rotation_primary_y0)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Distribution of crop rotation primary y0", x = "crop rotation primary y0", y = "Count")
```

### Wheat Type by Soil Texture

``` r
# Group data by crop rotation and soil texture and count observations
# Create a grouped bar plot to visualize the relationship
refined_metadata %>%
  group_by(crop_rotation_primary_y0, soil_texture) %>%
  summarise(Count = n()) %>%
  ggplot(aes(x = crop_rotation_primary_y0, y = Count, fill = soil_texture)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Wheat variety by Soil texture", x = "Wheat type", y = "Count")
```

### Chi-Square Test for Associations

The Chi-square test helps validate observed patterns by determining whether there is a statistically significant association or relationship between two categorical variables.

``` r
# Create a contingency table between crop rotation and take-all disease visibility
contingency_table <- table(refined_metadata$crop_rotation_primary_y0, refined_metadata$takeall_seen)

# Perform a Chi-square test to see if variables are independent
chi_square <- chisq.test(contingency_table)
chi_square

# Run Chi-square test with simulation for more robust p-values
chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000)
```

### Convert Categorical Variables to Factors

This step sets the data in a format that R can compute- factors. This prevents misinterpretation and enables us to actually plot the data.

``` r
# Convert several categorical variables to factors for analysis
refined_metadata <- refined_metadata %>%
  mutate(
    takeall_seen = factor(takeall_seen, levels = c("Y", "N")),
    fertiliser_use = factor(fertiliser_use),
    carbon = factor(carbon),
    crop_rotation_primary_y0 = factor(crop_rotation_primary_y0)
  )
```

### Clustering Analysis

By performing clustering analysis we can explore natural groupings or patterns in the environmental data that might not be obvious through basic summaries or visualisations. Essentially it reveals natural groupings that can guide further investigation and improve understanding. Other cluster options are hierarchical and visualisations can be done via PCA.

``` r


# Convert categorical variables into a one-hot encoded matrix
refined_metadata_encoded <- model.matrix(~ crop_rotation_primary_y0 + 
                                         soil_texture + tillage_method, 
                                         data = refined_metadata)

# Perform k-means clustering with 3 clusters
kmeans_result <- kmeans(refined_metadata_encoded, centers = 3)

# Add the cluster assignments back to the dataset
refined_metadata$Cluster <- as.factor(kmeans_result$cluster)

# Visualize clusters using jitter plot
refined_metadata %>%
  ggplot(aes(x = crop_rotation_primary_y0, y = soil_texture, color = Cluster)) +
  geom_jitter() +
  theme_minimal() +
  labs(title = "Clustering environmental data", x = "Wheat type", y = "Soil texture")
```
