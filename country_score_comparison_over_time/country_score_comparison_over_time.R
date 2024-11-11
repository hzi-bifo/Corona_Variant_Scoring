#!/usr/bin/env Rscript

#### Country Score Comparison Over Time
#Author: Katrina Norwood
#Last Updated: 01/11/24

# Script to plot the change in country antigenic scores over the three year period 
# of 01-2020 to 12-2023. The cumulative scores file contains countries and their 
# given antigenic score that eiher have 1% of that month's submitted isolates or 
# 500 sequences in that month. This was done so that the in the earlier months in 
# the pandemic, when there were fewer sequences, countries that were still sequencing
# could be represented, despite not having 500 sequences. Please see the 
# country_frequency_threshold.py script to see how this was done.

# To run:
# 1. start in the /Corona_Variant_Scoring/ home directory
# 2. in command line run: ```Rscript country_score_comparison_over_time/country_score_comparison_over_time.R```

library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
#output <- args[1]
#cumulative_scores <- read.csv(args[2], sep = "\t")

output <- "country_score_comparison_over_time/"
cumulative_scores <- read.csv("country_score_comparison_over_time/country_antigenic_score_with_threshold.tsv", sep = "\t")
print(head(cumulative_scores))

## Data Cleaning and Set up:

# Calculating the zscore of each country per month (to identify significant countries)
monthly_df_list <- list()

for (month in unique(cumulative_scores$date)) {
  print(month)
  # Separate by month
  month_df <- cumulative_scores[cumulative_scores$date == month, ]
  # Calculate zscores
  month_df$zscore <- (month_df$country_score - mean(month_df$country_score)) / sd(month_df$country_score)
  print(head(month_df))
  # Check output
  print(paste("Number of rows: ", nrow(month_df)))
  print(paste("Monthly SD: ", sd(month_df$country_score)))
  print(paste("Monthly Mean: ", mean(month_df$country_score)))
  # Add to list of DFs
  monthly_df_list <- append(monthly_df_list, list(month_df))
}

cumulative_scores_zscores <- do.call(rbind, monthly_df_list)
print(head(cumulative_scores_zscores))
print(paste("Length of cumulative_scores_zscores: ", nrow(cumulative_scores_zscores)))
print(paste("Length of cumulative_scores: ", nrow(cumulative_scores)))

# Creating a subplot with just a select number of countries
selected_countries <- c("USA", "Germany", "India", "South Africa", "Brazil", "United Kingdom")
cumulative_scores_countries_df <- cumulative_scores_zscores[cumulative_scores_zscores$Country %in% 
                                                      selected_countries, ]
cumulative_scores_countries_df <- cumulative_scores_countries_df %>%
  mutate(date = as.Date(paste0("01-", date), format="%d-%m-%Y")) %>%
  arrange(date)
print(paste("cumulative_scores_countries_df nrows: ", nrow(cumulative_scores_countries_df)))
      
# Creating a subplot with all other countries in it (to be shown as grey in the plot)
cumulative_scores_others_df <- cumulative_scores_zscores[!(cumulative_scores_zscores$Country  %in% 
                                                   selected_countries), ]
# Filtering out countries that have high sequencing for less than 30 months (to make plot cleaner)
cumulative_scores_others_df_subset <- cumulative_scores_others_df %>%
  group_by(Country) %>%
  filter(n() >= 30) %>%
  mutate(date = as.Date(paste0("01-", date), format="%d-%m-%Y")) %>%
  arrange(date)
print(paste("cumulative_scores_others_df_subset nrows: ", nrow(cumulative_scores_others_df_subset)))

sapply(cumulative_scores_countries_df, class)
sapply(cumulative_scores_others_df_subset, class)

## Plotting the results

plt1 <- ggplot() + 
  geom_line(data = cumulative_scores_others_df_subset, aes(x = date, y = country_score, 
                                                           group = Country, color = "Other")) +
  geom_line(data = cumulative_scores_countries_df, aes(x = date, y = country_score, 
                                                       group = Country, color = Country)) + 
  scale_x_date(date_breaks = "2 months", date_labels = "%m-%Y", expand = c(0, 0)) + 
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "grey", "#F0E442",
                                "#0072B2", "#009E73")) +
  labs(x = "Month", y = "Country Antigenic Alteration Score", color = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 7, vjust = 0.5, hjust = 1),
        axis.text.x = element_text(angle = 45), 
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))
plt1
ggsave(plt1, 
       filename = paste(output, "country_dynamics_ant_alt_scores_2020-01-2023-12_80x180.pdf", sep = ""),
       device = "pdf",
       height = 60, width = 180, units = "mm")

plt2 <- ggplot() + 
  geom_line(data = cumulative_scores_others_df_subset, aes(x = date, y = zscore, 
                                                           group = Country, color = "Other")) +
  geom_line(data = cumulative_scores_countries_df, aes(x = date, y = zscore, 
                                                       group = Country, color = Country)) + 
  scale_x_date(date_breaks = "2 months", date_labels = "%m-%Y", expand = c(0, 0)) + 
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "grey", "#F0E442",
                                "#0072B2", "#009E73")) +
  labs(x = "Month", y = "Country Standardized Z-Score", color = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 7, vjust = 0.5, hjust = 1),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 8), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))
plt2
ggsave(plt2, 
       filename = paste(output, "country_dynamics_ant_alt_zscores_2020-01-2023-12_80x180.pdf", sep = ""),
       device = "pdf",
       height = 50, width = 180, units = "mm")

knitr::knit_exit()