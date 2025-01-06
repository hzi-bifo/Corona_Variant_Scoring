#!/usr/bin/env Rscript

#### Country Score Comparison Over Time
#Author: Katrina Norwood
#Last Updated: 19/11/24

# Script to plot the change in country antigenic scores from 01-2020 to the most 
# recent month. Built to run in the CoVerage server pipeline. 

library(dplyr)
library(ggplot2)
library(gridExtra)
#install.packages("plotly")
library(plotly)

print <- base::print

args = commandArgs(trailingOnly=TRUE)
output <- args[1]
cumulative_scores <- read.csv(args[2], sep = ",")
filtered_countries <- read.csv(args[3], sep = "\t")

print("Output: ")
print(output)
print("cumulative_scores: ")
head(cumulative_scores)
print("filtered_countries")
head(filtered_countries)

# The filtered_countries df are countries used for each month (these are countries that pass the
# predetermined threshold, ie. they have either 1% of the global sequences or 500 sequences
# for the selected month, this is why there may be some gaps in the plot for different countries)

## Data Cleaning and Set up:

# Calculating the zscore of each country per month (to identify significant countries)
monthly_df_list <- list()

print("Creating DF with the zscores for the selected countries in each month: ")
for (month in unique(cumulative_scores$date)) {
  print(month)
  # Getting Countries that pass threshold test for that month
  month_countries_df <- filtered_countries[filtered_countries$date == month, ]
  month_countries_list <- as.list(month_countries_df$Country)
  # Separate by month
  month_df_unfiltered <- cumulative_scores[cumulative_scores$date == month, ]
  month_df <- month_df_unfiltered[month_df_unfiltered$Country %in% month_countries_list, ]
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

# Combining with the lineage information
print("Combining with lineage information: ")
cumulative_scores_zscores <- merge(cumulative_scores_zscores, filtered_countries,
                                                       by.x = c("Country", "date"),
                                                       by.y = c("Country", "date"))
print(head(cumulative_scores_zscores))
print(paste("Length of cumulative_scores_zscores: ", nrow(cumulative_scores_zscores)))

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

# Recombining the plots so that we have the desired countries (even if they don't
# meet the desired threshold) as well as the countries meeting the 30 month sequencing 
# threshold.
combined_df <- rbind(cumulative_scores_others_df_subset, cumulative_scores_countries_df) # check here, its changing the column type
print("Combined Df: ")
print(nrow(combined_df))
print(head(combined_df))

# Adding two columns for plotting countries with zscore above and below 1
print("Adding columns with thresholds: ")
combined_df$zscore_above <-  NA_real_
combined_df <- combined_df %>% mutate(zscore_above = if_else(zscore < 1, zscore_above, zscore))
print(head(combined_df))

# Formatting lineage information column, for easier reading
combined_df$lineages_information <- gsub("},", "},\n", combined_df$lineages_information)
combined_df$lineages_information <- gsub("'Pango lineage'", "Pango lineage", combined_df$lineages_information)
combined_df$lineages_information <- gsub("'lineage_frequency'", "Frequency (countrywise)", combined_df$lineages_information)
combined_df$lineages_information <- gsub("'antigenic_score'", "Antigenic Score", combined_df$lineages_information)
combined_df$lineages_information <- gsub("'zscore'", "Z-Score", combined_df$lineages_information)
print(head(combined_df))

## Plotting the results

# Country Score plot
plt1 <- ggplot() + 
  geom_line(data = cumulative_scores_others_df_subset, aes(x = date, y = country_score, 
                                                           group = Country, color = "Other",
                                                           text = paste("<span style='font-size:10px;'><b>Date:</b>", date, 
                                                                        "<br><b>Country:</b>", Country, 
                                                                        "<br><b>Country Score:</b>", round(country_score, 2), "</span>"))) +
  #aes(text = paste("Date:", date, "<br>Country:", Country, "<br>Country Score:", round(country_score, 2), "<br>Z-Score:", round(zscore, 2))) +
  geom_line(data = cumulative_scores_countries_df, aes(x = date, y = country_score, 
                                                       group = Country, color = Country,
                                                       text = paste("<span style='font-size:10px;'><b>Date:</b>", date, 
                                                                    "<br><b>Country:</b>", Country, 
                                                                    "<br><b>Country Score:</b>", round(country_score, 2), "</span>"))) + 
  scale_x_date(date_breaks = "2 months", date_labels = "%m-%Y", expand = c(0, 0)) + 
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "grey", "#F0E442",
                                "#0072B2", "#009E73")) +
  labs(x = "Month", y = "Country Antigenic Alteration Score", color = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, vjust = 0.5, hjust = 1),
        axis.text.x = element_text(angle = 45), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
#plt1

# Z-Scores plot (with given countries colored)
plt2 <- ggplot() + 
  geom_line(data = combined_df, 
            aes(x = date, y = zscore, group = Country,
                text = paste("<span style='font-size:10px;'><b>Date</b>:", date, 
                             "<br><b>Country:</b>", Country,  
                             "<br><b>Country Z-Score:</b>", round(zscore, 2), 
                             "<br><b>Lineage Information:</b>", lineages_information, "</span>")), 
            color = "grey") +
  geom_line(data = combined_df, 
            aes(x = date, y = zscore_above, group = Country, color = Country,
                text = paste("<span style='font-size:10px;'><b>Date:</b>", date, 
                             "<br><b>Country:</b>", Country, 
                             "<br><b>Country Z-Score:</b>", round(zscore, 2), 
                             "<br><b>Lineage Information:</b>", lineages_information, "</span>"))) +
  scale_x_date(date_breaks = "2 months", date_labels = "%m-%Y", expand = c(0, 0)) +
  labs(x = "Month", y = "Country Standardized Z-Score", color = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
#plt2

# Convert both ggplot objects to plotly .html objects
final_plt1 <- ggplotly(plt1, tooltip = "text")
final_plt2 <- ggplotly(plt2, tooltip = "text")
#final_plt2

final_subplot <- subplot(final_plt2, final_plt1, nrows = 2, shareX = TRUE, shareY = TRUE)
final_subplot <- final_subplot %>% layout(showlegend = FALSE)

# Saving final plot
htmlwidgets::saveWidget(final_subplot, file = paste(output, "country_dynamics_ant_alt_zscores.html", sep = ""))
