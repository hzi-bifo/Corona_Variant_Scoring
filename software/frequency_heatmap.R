#!/usr/bin/env Rscript

#library(devtools)
#install_version("d3heatmap", version = "0.6.1.2")

library(d3heatmap)
library(htmlwidgets)
library(countrycode)

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
  sign_lineages_df <- lineages[lineages[['antigenic_score']] >= 2.67, ]
  sign_lineages <- sign_lineages_df[['Pango.lineage']]
  # remove duplicates
  sign_lineages <- sort(unique(sign_lineages))
  
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
  
  # plot heatmap and save as html
  if (rename == TRUE){
    country_names <- countrycode(colnames(data), origin = 'iso2c', destination = 'country.name')
    country_names[country_names == "Hong Kong SAR China"] <- "Hong Kong" # rename Hong Kong to make it shorter
    colnames(data) <- country_names
  }
  data <- data[order(rowSums(data), decreasing = TRUE),]
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
