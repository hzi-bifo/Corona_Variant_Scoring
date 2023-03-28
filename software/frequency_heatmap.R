#!/usr/bin/env Rscript

#install.packages("devtools")
#library(devtools)
#install_version("d3heatmap", version = "0.6.1.2")

library(d3heatmap)
library(htmlwidgets)
library(countrycode)
library(DT)

args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
input <- read.csv(args[2], sep = "\t")
months_file <- args[3]
cutoff <- as.numeric(args[4])
output <- args[5]

print("Folder")
print(folder)
#print(input)
print("Months_file")
print(months_file)
print("Cutoff")
print(cutoff)
print("Output")
print(output)

get_heatmap <- function(folder, countries, current_month, output_name, rename, cutoff, lineages, output){
  
  # get significant lineages from all countries
  sign_lineages_df <- lineages[lineages[['antigenic_score']] >= 3.85, ]
  # remove duplicates
  sign_lineages_df[!duplicated(sign_lineages_df), ]
  sign_lineages <- sign_lineages_df[['Pango.lineage']]
  
  # get data
  data <- matrix(data = 0, nrow = length(sign_lineages), ncol = length(countries))
  colnames(data) <- countries
  rownames(data) <- sign_lineages
  
  for (country in countries){
    tmp <- read.csv(paste(folder, "/", country, "/", country, "_frequencies.txt", sep = ""), stringsAsFactors = FALSE, sep = "\t", header = TRUE)
    for (lineage in sign_lineages){
      freq <- tmp[which(tmp$lineage == lineage), paste("X", gsub("-", ".", current_month, fixed = TRUE), sep = "")]
      if (length(freq) > 0 && ! is.na(freq)){
        data[lineage, country] <- freq
      } else {
        data[lineage, country] <- 0
      }
    }
  }
  
  # remove lineages with frequency 0 in all countries
  data <- data[rowSums(data) > 0,,drop=FALSE]
  
  # write table
  write.csv(data, paste(output, "/", output_name, "_statistics.csv", sep = ""), quote = FALSE)
  
  # further filtering using the cutoff for plotting
  data <- data[rowSums(data) >= cutoff,,drop=FALSE]
  # keeping only the top 20 lineages for visualization
  data <- data[1:20,]
  
  # saving dataframe with identified pVOIs and their antigenic scores for later visualization
  pvois <- row.names(data)
  pvoi_df <- sign_lineages_df[sign_lineages_df$Pango.lineage %in% pvois, ]
  # adding frequency data 
  freq_df <- data
  global_freq <- rowSums(freq_df)
  freq_df <- cbind(freq_df, global_freq)
  freq_df_subset = subset(freq_df, select = c(global_freq))
  Pango.lineage <- row.names(freq_df_subset)
  freq_df_subset <- cbind(freq_df_subset, Pango.lineage)
  # binding to significant lineages df
  pvoi_df_final <- pvoi_df[c("Pango.lineage", "antigenic_score")]
  pvoi_final <- merge(pvoi_df_final, freq_df_subset, by=c("Pango.lineage"))
  #ordering by antigenic score
  pvoi_final <- pvoi_final[order(pvoi_final$antigenic_score, decreasing = TRUE),]  
  write.csv(pvoi_final, paste(output, "/", output_name, "_pVOI_table.csv", sep = ""), quote = FALSE, row.names = FALSE)

  # plot heatmap and save as html
  if (rename == TRUE){
    country_names <- countrycode(colnames(data), origin = 'iso2c', destination = 'country.name')
    country_names[country_names == "Hong Kong SAR China"] <- "Hong Kong" # rename Hong Kong to make it shorter
    colnames(data) <- country_names
  }
  #data <- data[order(rowSums(data), decreasing = TRUE),]
  heatmap <- d3heatmap(data, Rowv = NULL, Colv = NULL, colors = "Reds", xaxis_font_size = "12pt", yaxis_font_size = "12pt", xaxis_height = 120)
  saveWidget(heatmap, paste(output, "/", output_name, "_heatmap.html", sep = ""), selfcontained = TRUE)
}

# get all countries to be analyzed based on the existing folders
#countries <- dir(pattern = "^[A-Z][A-Z]") # countries and states (all folders with two upper case letters)
countries <- dir(path = folder, pattern = "^[A-Z][A-Z]$") # only countries, not states
states <- dir(path = folder, pattern = "^DE_") # only German states

months <- as.vector(t(read.table(months_file, stringsAsFactors = FALSE)))
current_month <- months[length(months)-1] # look at frequencies of the previous month to get representative values

get_heatmap(folder, countries, current_month, "antigenic_scoring_summary", TRUE, cutoff, input, output)
#get_heatmap(folder, countries, current_month, "antigenic_scoring_summary", FALSE, cutoff, input)
get_heatmap(folder, states, current_month, "antigenic_scoring_summary_states", FALSE, cutoff, input, output)

# saving dataframe with identified pVOIs and their antigenic scores
#data <- read.csv(paste(output, "/antigenic_scoring_summary_statistics.csv", sep = ""), row.names = 1, sep = ",")
#data <- data[rowSums(data) >= cutoff,,drop=FALSE]
#data <- data[1:20,]
#pvois <- row.names(data)

#sign_lineages_df <- input[input[['antigenic_score']] >= 3.85, ]
#sign_lineages_df[!duplicated(sign_lineages_df), ]

#pvoi_df <- sign_lineages_df[sign_lineages_df$Pango.lineage %in% pvois, ]
#pvoi_df_final <- pvoi_df[c("Pango.lineage", "antigenic_score")]
#print(pvoi_df_final)

#pvoi_df_final %>% datatable() # , options = list(pageLength = 1), rownames = FALSE)
#saveWidget(pvoi_table, paste(output, "/ntigenic_scoring_pVOItable.html", sep = ""), selfcontained = TRUE)
