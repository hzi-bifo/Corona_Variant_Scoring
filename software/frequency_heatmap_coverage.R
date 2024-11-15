#!/usr/bin/env Rscript

#install.packages("devtools")
#library(devtools)
#install_version("d3heatmap", version = "0.6.1.2")

library(htmltools)
library(dplyr)
library(d3heatmap)
library(plotly)
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

print <- base::print

print("Folder")
print(folder)
print("Months_file")
print(months_file)
print("Cutoff")
print(cutoff)
print("Output")
print(output)

get_heatmap <- function(folder, countries, current_month, output_name, rename, cutoff, lineages, output, variantDf){
  # Function to build either a global heatmap with zscores or a German regional 
  # heatmap of frequency of selected lineages per region
  
  ## Data Cleaning
  
  # get significant lineages from all countries
  sign_lineages_df <- lineages[lineages["significant"] == "yes", ]
  
  # remove duplicates and find lineages of interest
  sign_lineages_df[!duplicated(sign_lineages_df), ]
  sign_lineages_list <- sign_lineages_df[["Pango.lineage"]]
  number_sign_lineages <- length(sign_lineages_list)
  lineages_filtered <- lineages[!lineages$Pango.lineage %in% sign_lineages_list, ]
  top_lineages_df <- lineages_filtered[lineages[['antigenic_score']] >= 0, ] # should there be no significant lineages that month, use the top scoring lineages
  top_lineages_list <- top_lineages_df[["Pango.lineage"]]
  sign_lineages <- unique(c(sign_lineages_list, top_lineages_list)) # Combines any top scoring lineages with significant lineages 
  sign_lineages <- sign_lineages[!is.na(sign_lineages)]
  
  # get data
  print("Getting Data: ")
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
  print("Filtering using the cutoff for plotting")
  data <- data[rowSums(data) >= cutoff,,drop=FALSE]
  # keeping only the top 20 lineages with highest frequency for visualization
  freq_order <- order(rowSums(data), decreasing = TRUE)
  data <- data[freq_order, ]
  if (nrow(data) > 20){
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
  
  # creating list of lineages and mutations (taking the first mutation list occurrence for each lineage)
  mutations <- variantDf[match(unique(variantDf$Pango.lineage), variantDf$Pango.lineage),]
  tempDF <- merge(x = input, y = mutations, by = "Pango.lineage", all.x = TRUE)
  finalDF <- merge(x = tempDF, y = variantFreq, by = "Pango.lineage", all.x = TRUE)
  finalDF$antigenic_score <- round(finalDF$antigenic_score, 2)
  finalDF <- unique(finalDF)
  
  # Adding value for missing z-scores
  finalDF$zscore[is.na(finalDF$zscore)] <- "NA"
  
  # Formatting dataframe for JSON file
  colnames(finalDF)[which(names(finalDF) == "Pango.lineage")] <- "lineage"
  colnames(finalDF)[which(names(finalDF) == "AA.Substitutions")] <- "substitutions"
  
  # Adding spaces to string list of substitutions
  finalDF$substitutions <- gsub(",([A-Za-z])", ", \\1", finalDF$substitutions)
  
  # Saving as a json file to be viewed on the CoVerage website
  jsondata <- toJSON(finalDF, pretty = TRUE)
  write(jsondata, file=paste(output, "/", output_name, "_lineages_table.json", sep = ""))

  ## Data Visualization 
  
  # plot heatmap and save as html
  if (rename == TRUE){
    country_names <- countrycode(colnames(data_subset), origin = 'iso2c', destination = 'country.name')
    country_names[country_names == "Hong Kong SAR China"] <- "Hong Kong" # rename Hong Kong to make it shorter
    colnames(data_subset) <- country_names
  }
  
  # Adding a asterisk to significant lineages
  #data <- data[order(rowSums(data), decreasing = TRUE),]
  # Make a copy of the data_subset (which will be used later with the zscores)
  data_subset_no_asterisk <- data.frame(data_subset)
  if (length(sign_lineages_list) > 0){
	  print("Adjusting row names: ")
  	  for (i in sign_lineages_list){
		  print(i)
		  row.names(data_subset)[row.names(data_subset) == i] <- paste0(i, "*")}
  }
  print("final data:")
  print(data_subset)
  
  # Saving heatmap data
  heatmap <- d3heatmap(data_subset, Rowv = NULL, Colv = NULL, colors = "Reds", xaxis_font_size = "12pt", yaxis_font_size = "12pt", xaxis_height = 120)
  saveWidget(heatmap, paste(output, "/", output_name, "_heatmap.html", sep = ""), selfcontained = TRUE)
  
  # For global lineages we add a heatmap with the zscores as well as frequencies
  if (rename == TRUE) {
    # Plotting Z-Scores as a heatmap
    print("Plotting z-scores as a heatmap: ")
    finalDF_filtered <- finalDF %>% filter(lineage %in% top_variants)
    zscores <- subset(finalDF_filtered, select = c(lineage, zscore))
    zscores$zscore <- as.numeric(zscores$zscore)  
    #zscores$zscore <- ifelse(is.na(zscores$zscore), NA, as.numeric(zscores$zscore))
    
    # Reordering rows by frequency (like that of the data_subset)
    print("Reordering rows by frequency: ")
    lineages_order <- rownames(data_subset_no_asterisk)
    zscores <- zscores[match(lineages_order, zscores$lineage),]
    
    # Setting rownames as lineages
    print("Setting rownames as lineages: ")
    print("rownames: ")
    print(rownames(zscores))
    print("zscores$lineage: ")
    print(zscores$lineage)
    rownames(zscores) <- zscores$lineage
    
    # Renaming column
    zscores_df <- subset(zscores, select = -c(1))
    colnames(zscores_df) <- "Zscore"
    print("zscores_df")
    print(zscores_df)
    
    # Combining zscore dataframe with the data_subset
    #combined_data <- cbind(data_subset, zscores_df)
    
    # Renaming row if they are significantly antigenically altered
    print("Renaming significant rows: ")
    if (length(sign_lineages_list) > 0){
      print("Adjusting row names: ")
      for (i in sign_lineages_list){
        print(i)
        row.names(zscores_df)[row.names(zscores_df) == i] <- paste0(i, "*")}
    }

    # Calculate the minimum and maximum values for the z-score column to set the range
    zscore_min <- min(zscores_df$Zscore, na.rm = TRUE)
    zscore_max <- max(zscores_df$Zscore, na.rm = TRUE)
    
    # Defining the two colorscales
    frequency_colorscale <- list(
      list(0, "#fff5f0"),
      list(0.5, "#fcbba1"), 
      list(1, "#a50f15"))
    zscore_colorscale <- list(
      list(0, "blue"),     # For minimum values (negative)
      list(0.5, "#cbc9e2"),  # Neutral point (zero)
      list(1, "#54278f")       # For maximum values (positive)
    )
    
    # Reverse row order for heatmap
    print("reversing rows for order in heatmap: ")
    reversed_rows_data_subset <- rev(rownames(data_subset))
    reversed_rows_zscores_df <- rev(rownames(zscores_df))
    print(reversed_rows_data_subset)
    print(reversed_rows_zscores_df)
    
    # Create separate heatmaps for frequency and z-score
    frequency_heatmap <- plot_ly(
      #z = as.matrix(data_subset),
      z = as.matrix(data_subset[reversed_rows_data_subset, , drop = FALSE]), 
      x = colnames(data_subset),
      #y = rownames(data_subset),
      y = reversed_rows_data_subset,
      type = "heatmap",
      colorscale = frequency_colorscale,
      colorbar = list(title = list(text = "Frequency", font = list(size = 12))),
      hovertemplate = paste(
        "Country: %{x}<br>",
        "Lineage: %{y}<br>",
        "Frequency: %{z}<br>",
        "<extra></extra>")
    )
    
    zscore_heatmap <- plot_ly(
      #z = matrix(zscores_df$Zscore, ncol = 1),  
      z = matrix(zscores_df[reversed_rows_zscores_df, "Zscore"], ncol = 1),
      x = list("Zscore"),                      
      #y = rownames(zscores_df),
      y = reversed_rows_zscores_df,
      type = "heatmap",
      colorscale = zscore_colorscale,
      zmin = zscore_min,
      zmax = zscore_max,
      colorbar = list(title = list(text = "Z-score", font = list(size = 12)), 
                      tickvals = c(zscore_min, 0, zscore_max)),
      hovertemplate = paste(
        "Lineage: %{y}<br>",
        "Z-score: %{z}<br>",
        "<extra></extra>")
    )
    
    # Use subplot to combine them closely
    heatmap <- subplot(
      frequency_heatmap,
      zscore_heatmap,
      nrows = 1,           
      shareY = TRUE,     
      widths = c(0.95, 0.05)  # Adjust size of z-score column
    ) %>%
      layout(
        title = list(text = current_month,
          font = list(size = 16, weight = "bold"), x = 0.06), 
        xaxis = list(
          tickfont = list(size = 12),
          tickangle = 45
        ),
        xaxis2 = list(
          tickfont = list(size = 12),
          tickangle = 45
          ),
        yaxis = list(tickfont = list(size = 12))#,
        #width = 800,  
        #height = 600 
      )
    
    htmlwidgets::saveWidget(heatmap, paste(output, "/", output_name, "_combined_heatmap_plotly.html", sep = ""), selfcontained = TRUE)
  }
}

# get all countries to be analyzed based on the existing folders
countries <- dir(path = folder, pattern = "^[A-Z][A-Z]$") # only countries, not states
states <- dir(path = folder, pattern = "^DE_") # only German states

months <- as.vector(t(read.table(months_file, stringsAsFactors = FALSE)))
current_month <- months[length(months)-1] # look at frequencies of the previous month to get representative values

get_heatmap(folder, countries, current_month, "antigenic_scoring_summary", TRUE, cutoff, input, output, variantDf)
get_heatmap(folder, states, current_month, "antigenic_scoring_summary_states", FALSE, cutoff, input, output, variantDf)
