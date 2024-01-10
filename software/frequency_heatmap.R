#!/usr/bin/env Rscript

#install.packages("devtools")
#library(devtools)
#install_version("d3heatmap", version = "0.6.1.2")

library(d3heatmap)
library(htmlwidgets)
library(countrycode)
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
input <- read.csv(args[2], sep = "\t")
months_file <- args[3]
cutoff <- as.numeric(args[4])
output <- args[5]
variantDf <- read.csv(args[6], sep = "\t", colClasses = c("NULL","NULL","NULL", NA, NA, "NULL", "NULL"))

print("Folder")
print(folder)
print("Months_file")
print(months_file)
print("Cutoff")
print(cutoff)
print("Output")
print(output)

get_heatmap <- function(folder, countries, current_month, output_name, rename, cutoff, lineages, output, variantDf){
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
    current_file <- paste(folder, "/", country, "/", country, "_frequencies.txt", sep = "")
    if (file.exists(current_file) && file.info(current_file)$size > 0) {
    	tmp <- read.csv(paste(folder, "/", country, "/", country, "_frequencies.txt", sep = ""), stringsAsFactors = FALSE, sep = "\t", header = TRUE)}
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
  if (nrow(data) > 2) {
	  stop <- FALSE
	  print(data)
  }
  else {
  	stop <- TRUE}
  print("Stop:")
  print(stop)

  if (stop == FALSE) { 

  	# write table
  	write.csv(data, paste(output, "/", output_name, "_statistics.csv", sep = ""), quote = FALSE)
  
  	# further filtering using the cutoff for plotting
  	data <- data[rowSums(data) >= cutoff,,drop=FALSE]
  	# keeping only the top 20 lineages with highest frequency for visualization
  	freq_order <- order(rowSums(data), decreasing = TRUE)
 	 data <- data[freq_order, ]
  	if (nrow(data) > 21) {
	  	data_subset <- data[1:20,]
  	} else {
  		data_subset <- data
  	}
  	# dropping countries with low frequency, to reduce heatmap noise
  	data_subset <- subset(data_subset, select = (colSums(data_subset) > 0))
  
  	# saving dataframe with identified antigenically altering lineages and their antigenic scores for later visualization
  	top_variants <- row.names(data_subset)
  	input[ , c('rank', 'WHO_label')] <- list(NULL)
  	# Calculating global frequency of lineages
  	variantFreq <- as.data.frame(table(variantDf$Pango.lineage))
  	variantFreq$frequent <- variantFreq$Freq / nrow(variantDf)
  	variantFreq[ , c('Freq')] <- list(NULL)
  	colnames(variantFreq)[which(names(variantFreq) == "Var1")] <- "Pango.lineage"
  	#input$frequent <- sapply(input$Pango.lineage, function(x) x %in% top_variants)
  	# creating list of lineages and mutations (taking the first mutation list occurrence for each lineage)
  	mutations <- variantDf[match(unique(variantDf$Pango.lineage), variantDf$Pango.lineage),]
  	tempDF <- merge(x = input, y = mutations, by = "Pango.lineage", all.x = TRUE)
  	finalDF <- merge(x = tempDF, y = variantFreq, by = "Pango.lineage", all.x = TRUE)
  	finalDF$antigenic_score <- round(finalDF$antigenic_score, 2)
  	finalDF <- unique(finalDF)
  	# Formatting dataframe for JSON file
  	colnames(finalDF)[which(names(finalDF) == "Pango.lineage")] <- "lineage"
  	colnames(finalDF)[which(names(finalDF) == "AA.Substitutions")] <- "substitutions"
  	# Adding spaces to string list of substitutions
  	finalDF$substitutions <- gsub(",([A-Za-z])", ", \\1", finalDF$substitutions)
  	# Saving as a json file
  	jsondata <- toJSON(finalDF, pretty = TRUE)
  	write(jsondata, file=paste(output, "/", output_name, "_lineages_table.JSON", sep = ""))

  	# plot heatmap and save as html
  	if (rename == TRUE){
    		country_names <- countrycode(colnames(data_subset), origin = 'iso2c', destination = 'country.name')
    		country_names[country_names == "Hong Kong SAR China"] <- "Hong Kong" # rename Hong Kong to make it shorter
    		colnames(data_subset) <- country_names
  	}
  	#data <- data[order(rowSums(data), decreasing = TRUE),]
  	heatmap <- d3heatmap(data_subset, Rowv = NULL, Colv = NULL, colors = "Reds", xaxis_font_size = "12pt", yaxis_font_size = "12pt", xaxis_height = 120, main = current_month)
  	saveWidget(heatmap, paste(output, "/", output_name, "_heatmap.html", sep = ""), selfcontained = TRUE)
  }
  return(stop)
}

# get all countries to be analyzed based on the existing folders
#countries <- dir(pattern = "^[A-Z][A-Z]") # countries and states (all folders with two upper case letters)
countries <- dir(path = folder, pattern = "^[A-Z][A-Z]$") # only countries, not states
states <- dir(path = folder, pattern = "^DE_") # only German states

months <- as.vector(t(read.table(months_file, stringsAsFactors = FALSE)))
month_len <- 1
current_month <- months[length(months)-month_len]
#current_month <- months[length(months)-1] # look at frequencies of the previous month to get representative values

print("Countries:")
x <- get_heatmap(folder, countries, current_month, "antigenic_scoring_summary", TRUE, cutoff, input, output, variantDf)
if (x == TRUE) {
	while (x == TRUE){
		print("Running Again")
		month_len <- month_len + 1
    		current_month <- months[length(months)-month_len]
    		print(current_month)
		x <- get_heatmap(folder, countries, current_month, "antigenic_scoring_summary", TRUE, cutoff, input, output, variantDf)
	}
}
#get_heatmap(folder, countries, current_month, "antigenic_scoring_summary", FALSE, cutoff, input)

print("States:")
y <- get_heatmap(folder, states, current_month, "antigenic_scoring_summary_states", FALSE, cutoff, input, output, variantDf)
if (y == TRUE) {
	while (y == TRUE){
		print("Running again")
        	month_len <- month_len + 1
                current_month <- months[length(months)-month_len]
                print(current_month)
                y <- get_heatmap(folder, states, current_month, "antigenic_scoring_summary_states", FALSE, cutoff, input, output, variantDf)
	}
}
