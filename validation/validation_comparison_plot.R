#!/usr/bin/env Rscript

# Script to validate the different proposed methodologies for scoring SARS-CoV-2 
# circulating lineages based on their potential antigenicity. Compares these scores
# to both antigenic cartography distances as well as averaged mFRN values from
# previous literature.

library(dplyr)
library(ggplot2)
library(psych)

## Importing Data Inputs

args = commandArgs(trailingOnly=TRUE)
output <- args[1]
mfrn <- read.csv(args[2], sep = "\t")
antigenic_cartography_distances <- args[3]
## Without weights at antigenic sites (method1) 
method1_antigenicScores_dir <- args[4]
## With weights at all sites (method2)
method2_antigenicScores_dir <- args[5]
## Weights at antigenic sites (method3)
method3_antigenicScores_dir <- args[6]
## Without weights at all sites - Baseline (method4)
method4_antigenicScores_dir <- args[7]
## Weights without directionality at all sites
method5_antigenicScores_dir <- args[8]
## With new weights at all sites (using amino acid changes that occurred at least 3 times throughout the tree)
method6_antigenicScores_dir <- args[9]

## Data Cleaning and Normalization for the mAb neutralization scores

# mFRN values method 1 - taking means without normalization

# Setting mAb names as the row names for clearer organization
rownames(mfrn) <- mfrn$mAb
mfrn$mAb <- NULL
print(mfrn)
# Getting means of mfrn values per variant
mfrn_means <- c()
for (column in colnames(mfrn)){
  print(column)
  mean_val <- mean(mfrn[[column]], na.rm = TRUE)
  mfrn_means <- append(mfrn_means, mean_val)
}
# Making dataframe of mfrn mean values
variant <- colnames(mfrn)
mfrn_df <- data.frame(variant, mfrn_means)
mfrn_df

# mFRN values method 2 - taking the median value without normalization
#rownames(mfrn) <- mfrn$mAb
#mfrn$mAb <- NULL
#print(mfrn)
## Getting means of mfrn values per variant
#mfrn_medians <- c()
#for (column in colnames(mfrn)){
#  print(column)
#  median_val <- median(mfrn[[column]], na.rm = TRUE)
#  mfrn_medians <- append(mfrn_medians, median_val)
#}
# Making dataframe of mfrn mean values
#variant <- colnames(mfrn)
#mfrn_df_median <- data.frame(variant, mfrn_medians)
#mfrn_df_median
#mfrn_df

# mFRN values method 3 - taking mean value of normalized data

# Normalizing mFRN data prior to taking the mean
#min_max_normalization <- function(x) {
#  (x - min(x)) / (max(x) - min(x))}
#mfrn_no_missing_values <- na.omit(mfrn)
#mfrn_norm <- as.data.frame(lapply(mfrn_no_missing_values[1:9], min_max_normalization))
#mfrn_norm_means <- c()
#for (column in colnames(mfrn_norm)){
#  print(column)
#  mean_val <- mean(mfrn_norm[[column]], na.rm = TRUE)
#  mfrn_norm_means <- append(mfrn_norm_means, mean_val)
#}
#variant <- colnames(mfrn_norm)
#mfrn_df_norm <- data.frame(variant, mfrn_norm_means)
#mfrn_df_norm

# Raw mFRN values
#variant_mFRN_values_conversion <- function(mfrn_df, column_name) {
#  variant_mfrn_values <- mfrn[c(column_name)]
#  variant_mfrn_values$variant <- column_name
#  variant_mfrn_values <- variant_mfrn_values %>% rename_at(column_name, ~'mfrn_value')
#  return(variant_mfrn_values)
#}

#alpha_mfrn_df <- variant_mFRN_values_conversion(mfrn, "Alpha")
#beta_mfrn_df <- variant_mFRN_values_conversion(mfrn, "Beta")
#delta_mfrn_df <- variant_mFRN_values_conversion(mfrn, "Delta")
#gamma_mfrn_df <- variant_mFRN_values_conversion(mfrn, "Gamma")
#ba1_mfrn_df <- variant_mFRN_values_conversion(mfrn, "BA.1")
#ba2_mfrn_df <- variant_mFRN_values_conversion(mfrn, "BA.2")
#ba11_mfrn_df <- variant_mFRN_values_conversion(mfrn, "BA.1.1")
#ba45_mfrn_df <- variant_mFRN_values_conversion(mfrn, "BA.4_BA.5")
#ba2121_mfrn_df <- variant_mFRN_values_conversion(mfrn, "BA.2.12.1")

#mfrn_raw_df <- rbind(alpha_mfrn_df, beta_mfrn_df, delta_mfrn_df, gamma_mfrn_df,
#                     ba1_mfrn_df, ba2_mfrn_df, ba11_mfrn_df, ba45_mfrn_df, ba2121_mfrn_df)
#na.omit(mfrn_raw_df)
#rownames(mfrn_raw_df) <- NULL
#mfrn_raw_df

## Data Cleaning for antigenic cartography distances

# Calculating the euclidean distances from the D614G variant
antigenic_cartography_distances

# Changing the names of lineages to WHO format
antigenic_cartography_distances$variant[antigenic_cartography_distances$variant == 'B.1.1.7'] <- "Alpha"
antigenic_cartography_distances$variant[antigenic_cartography_distances$variant == 'B.1.351'] <- "Beta"
antigenic_cartography_distances$variant[antigenic_cartography_distances$variant == 'P.1'] <- "Gamma"
antigenic_cartography_distances$variant[antigenic_cartography_distances$variant == 'B.1.617.2'] <- "Delta"
antigenic_cartography_distances$variant[antigenic_cartography_distances$variant == 'C.37'] <- "Lambda"
antigenic_cartography_distances$variant[antigenic_cartography_distances$variant == 'B.1.617.1'] <- "Kappa"
antigenic_cartography_distances$variant[antigenic_cartography_distances$variant == 'B.1.429'] <- "Epsilon"
antigenic_cartography_distances$variant[antigenic_cartography_distances$variant == 'B.1.621'] <- "Mu"
antigenic_cartography_distances

## Data Cleaning and Normalization for VOC Antigenic Scores

create_dataframe <- function(indir){
  # Function for data import and cleaning for the antigenic scoring results
  
  antigenicScores_df <- data.frame(matrix(ncol = 4, nrow = 0))
  variants_pangoLineage <- c("B.1.1.7", "B.1.351", "P.1", "B.1.429", "B.1.617.2", "C.37", "B.1.621", "BA.1", "BA.2", "BA.1.1", "BA.2.12.1", "BA.4", "BA.5", "XBB", "XBB.1.16", "XBB.1.5", "XBB.2.3", "EG.5", "EG.5.1", "JN.1")
  for (file in list.files(path = indir, pattern = ".csv", all.files = TRUE, full.names = TRUE)) {
    print(file)
    df <- read.csv(file, sep = "\t")
    print(nrow(df))
    variants_df <- df[df$Pango.lineage %in% variants_pangoLineage, ]
    print(nrow(variants_df))
    antigenicScores_df <- rbind(antigenicScores_df, variants_df)}
  
  print(antigenicScores_df)
  
  # Taking median antigenic score per pango lineage
  antigenicScores_df <- select(antigenicScores_df, c("Pango.lineage", "antigenic_score"))
  #antigenicScores_df_mean <- aggregate(.~Pango.lineage, antigenicScores_df, mean) # Here we take the mean score of the lineages
  antigenicScores_df_mean <- aggregate(.~Pango.lineage, antigenicScores_df, median) # Takes the median score of the lineages as there may have been high variance in each lineage based on the boxplot comparison
  #antigenicScores_df_mean %>% mutate_at(c("antigenic_score"), as.numeric)
  print("antigenicScores_df_mean")
  print(antigenicScores_df_mean)
  
  # Taking the mean of lineages BA.4 and BA.5 to compare with mfrn values
  ba4_ba5 <- c("BA.4_BA.5")
  value1_df <- antigenicScores_df_mean[antigenicScores_df_mean$Pango.lineage == 'BA.4', ]
  value2_df <- antigenicScores_df_mean[antigenicScores_df_mean$Pango.lineage == 'BA.5', ]
  ba45_mean <- mean(c(value1_df$antigenic_score[1], value2_df$antigenic_score[1]))
  print("BA4_5 mean:")
  print(ba45_mean)
  #ba45_mean <- mean(c(antigenicScores_df_mean$antigenic_score[8], antigenicScores_df_mean$antigenic_score[9]))
  ba4_ba5 <- c(ba4_ba5, ba45_mean)
  antigenicScores_df_mean <- rbind(antigenicScores_df_mean, ba4_ba5)
  #antigenicScores_df_mean <- antigenicScores_df_mean[-c(8,9), ]
  print(antigenicScores_df_mean)
  print(colnames(antigenicScores_df_mean))
  
  # Renaming Pango lineages
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'B.1.1.7'] <- "Alpha"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'B.1.351'] <- "Beta"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'P.1'] <- "Gamma"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'B.1.617.2'] <- "Delta"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'C.37'] <- "Lambda"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'B.1.621'] <- "Mu"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'B.1.429'] <- "Epsilon"
  antigenicScores_df_renamed <- antigenicScores_df_mean
  antigenicScores_df_renamed <- antigenicScores_df_renamed %>% rename_at('Pango.lineage', ~'variant')
  antigenicScores_df_renamed <- antigenicScores_df_renamed[-c(10,11), ]
  print(antigenicScores_df_renamed)
  
  return(antigenicScores_df_renamed)
}

method1_antigenicScores <- create_dataframe(method1_antigenicScores_dir)
method2_antigenicScores <- create_dataframe(method2_antigenicScores_dir)
method3_antigenicScores <- create_dataframe(method3_antigenicScores_dir)
method4_antigenicScores <- create_dataframe(method4_antigenicScores_dir)
method5_antigenicScores <- create_dataframe(method5_antigenicScores_dir)
method6_antigenicScores <- create_dataframe(method6_antigenicScores_dir)

## Data Visualization

data_visualization <- function(antigenic_scores_df, antigenic_distances_df, mfrn_df) {
  # Function for merging mFRN values, antigenic distances, and antigenic scores for a final DF for visualization
  
  # Creating one data frame with all comparisons
  antigenicScores_df_filtered <- antigenic_scores_df # Copying input data
  antigenicScores_df_filtered$antigenic_score <- as.numeric(antigenicScores_df_filtered$antigenic_score)
  antigenic_distances_df$comparison <- 'Antigenic Distances'
  antigenic_distances_df <- antigenic_distances_df %>% rename_at('distance_from_D614G', ~'value')
  merged1 <- merge(x = antigenic_distances_df, y = antigenicScores_df_filtered, by = c("variant")) 
  mfrn_df$comparison <- 'mFRN Values'
  
  if ('mfrn_means' %in% names(mfrn_df)) {
    mfrn_df <- mfrn_df %>% rename_at('mfrn_means', ~'value')} else {
    mfrn_df <- mfrn_df %>% rename_at('mfrn_medians', ~'value')}
  merged2 <- merge(x = mfrn_df, y = antigenicScores_df_filtered, by = c("variant")) 
  
  print("merged1: ")
  print(merged1)
  print("merged2: ")
  print(merged2)
  
  df4 <- rbind(merged1, merged2)
  df4$value <- round(df4$value, digits = 2)
  df4$antigenic_score <- round(df4$antigenic_score, digits = 2)
  df4
  
  return(df4)
}

# Mean mFRN values used
method1_df_vis <- data_visualization(method1_antigenicScores, antigenic_cartography_distances, mfrn_df)
method2_df_vis <- data_visualization(method2_antigenicScores, antigenic_cartography_distances, mfrn_df)
method3_df_vis <- data_visualization(method3_antigenicScores, antigenic_cartography_distances, mfrn_df)
method4_df_vis <- data_visualization(method4_antigenicScores, antigenic_cartography_distances, mfrn_df)
method5_df_vis <- data_visualization(method5_antigenicScores, antigenic_cartography_distances, mfrn_df)
method6_df_vis <- data_visualization(method6_antigenicScores, antigenic_cartography_distances, mfrn_df)
# Raw mFRN values used
#method2_df_vis_mfrn_median <- data_visualization(method2_antigenicScores, antigenic_cartography_distances, mfrn_df_median)

## Visualization

scatterplot <- function (method_df_vis, ylim, file_name) {
  method_plt <- ggplot(method_df_vis, aes(x = value, y = antigenic_score, group = 1, colour = variant)) +
    theme_bw() +
    geom_point(size = 2.5) +
    ylim(0, ylim) +
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#332288", "#b66dff", "#999933", "slategray1", "violetred1", "springgreen", "slateblue", "tomato", "tan3", "lightpink4")) +   
    facet_wrap(~comparison, scales = "free_x") + #, strip.position = "bottom"
    labs(x = NULL, y = "Antigenic Alterations Score") +
    theme(text = element_text(family = "Helvetica"), strip.placement = "outside", strip.text.x = element_text(size = 10), axis.text = element_text(size=10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.title = element_text(size=10), legend.text = element_text(size=8)) #strip.background = element_blank(),
  method_plt
  # Saving plot
  ggsave(method_plt, 
         filename = paste(output, file_name, sep = ""),
         device = "pdf",
         height = 45, width = 180, units = "mm")
  return(method_plt)
}

# Scatterplots using averaged mFRN values
method1_plt <- scatterplot(method1_df_vis, 14, "antigenicScores_without_weights_at_antigenic_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf")
method2_plt <- scatterplot(method2_df_vis, 25, "antigenicScores_weights_at_all_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf")
method3_plt <- scatterplot(method3_df_vis, 7, "antigenicScores_weights_at_antigenic_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf")
method4_plt <- scatterplot(method4_df_vis, 59, "antigenicScores_noweights_at_all_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf")
method5_plt <- scatterplot(method5_df_vis, 27, "antigenicScores_weights_at_all_sites_nodirectionality_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf")
method6_plt <- scatterplot(method6_df_vis, 17, "antigenicScores_new_weights_at_all_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf")

# Scatterplots using "raw" mFRN values
#method2_plt_raw_mfrn <- scatterplot(method2_df_vis_mfrn_median, 25, "antigenicScores_without_weights_at_antigenic_sites_comparison_raw-mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf")
#method2_plt_raw_mfrn
#method1_plt

## Statistical Testing Between Antigenic Distances and Scores

# Data Set Up

split_df_cartography <- function(input_df) {
  df_split_cartography <- input_df[input_df$comparison == 'Antigenic Distances', ]
  df_split_cartography$antigenic_score <- as.numeric(df_split_cartography$antigenic_score)
  return(df_split_cartography)
}

split_df_neutralization <- function(input_df) {
  split_df_neutralization <- input_df[input_df$comparison == 'mFRN Values', ]
  split_df_neutralization$antigenic_score <- as.numeric(split_df_neutralization$antigenic_score)
  return(split_df_neutralization)
}

# Method 1
method1_split_cartography <- split_df_cartography(method1_df_vis)
method1_split_neutralization <- split_df_neutralization(method1_df_vis)
cartography_distances <- as.matrix(method1_split_cartography['value']) # Only need to do this once as its the same across all methods
neutralization_values <- as.matrix(method1_split_neutralization['value']) # Only need to do this once as its the same across all methods
method1_antigenic_scores <- as.matrix(method1_split_cartography['antigenic_score']) #as.vector
method1_antigenic_scores_mfrn <- as.matrix(method1_split_neutralization['antigenic_score']) #as.vector
cartography_distances
neutralization_values
method1_antigenic_scores
method1_antigenic_scores_mfrn
# Method 2
method2_split_cartography <- split_df_cartography(method2_df_vis)
method2_split_neutralization <- split_df_neutralization(method2_df_vis)
method2_antigenic_scores <- as.matrix(method2_split_cartography['antigenic_score']) #as.vector
method2_antigenic_scores_mfrn <- as.matrix(method2_split_neutralization['antigenic_score']) #as.vector
method2_antigenic_scores
method2_antigenic_scores_mfrn
# Method 3
method3_split_cartography <- split_df_cartography(method3_df_vis)
method3_split_neutralization <- split_df_neutralization(method3_df_vis)
method3_antigenic_scores <- as.matrix(method3_split_cartography['antigenic_score']) #as.vector
method3_antigenic_scores_mfrn <- as.matrix(method3_split_neutralization['antigenic_score']) #as.vector
method3_antigenic_scores
method3_antigenic_scores_mfrn
# Method 4
method4_split_cartography <- split_df_cartography(method4_df_vis)
method4_split_neutralization <- split_df_neutralization(method4_df_vis)
method4_antigenic_scores <- as.matrix(method4_split_cartography['antigenic_score']) #as.vector
method4_antigenic_scores_mfrn <- as.matrix(method4_split_neutralization['antigenic_score']) #as.vector
method4_antigenic_scores
method4_antigenic_scores_mfrn
# Method 5
method5_split_cartography <- split_df_cartography(method5_df_vis)
method5_split_neutralization <- split_df_neutralization(method5_df_vis)
method5_antigenic_scores <- as.matrix(method5_split_cartography['antigenic_score']) #as.vector
method5_antigenic_scores_mfrn <- as.matrix(method5_split_neutralization['antigenic_score']) #as.vector
method5_antigenic_scores
method5_antigenic_scores_mfrn
# Method 6
method6_split_cartography <- split_df_cartography(method6_df_vis)
method6_split_neutralization <- split_df_neutralization(method6_df_vis)
method6_antigenic_scores <- as.matrix(method6_split_cartography['antigenic_score']) #as.vector
method6_antigenic_scores_mfrn <- as.matrix(method6_split_neutralization['antigenic_score']) #as.vector
method6_antigenic_scores
method6_antigenic_scores_mfrn

# Test with median mFRN values
#method2_split_neutralization_median <- split_df_neutralization(method2_df_vis_mfrn_median)
#neutralization_values_median <- as.matrix(method2_split_neutralization_median['value'])
#method2_antigenic_scores_mfrn_median <- as.matrix(method2_split_neutralization_median['antigenic_score']) #as.vector
#method2_antigenic_scores_mfrn_median

# Testing for outliers within the data
boxplot_mfrn_values <- ggplot(method1_split_neutralization, aes(x = comparison, y = value)) + geom_boxplot()
boxplot_mfrn_values

boxplot_cartography_values <- ggplot(method1_split_cartography, aes(x = comparison, y = value)) + geom_boxplot()
boxplot_cartography_values

boxplot_method1_antigenic_scores <- ggplot(method1_df_vis, aes(x = comparison, y = antigenic_score)) + geom_boxplot()
boxplot_method1_antigenic_scores

boxplot_method2_antigenic_scores <- ggplot(method2_df_vis, aes(x = comparison, y = antigenic_score)) + geom_boxplot()
boxplot_method2_antigenic_scores

boxplot_method3_antigenic_scores <- ggplot(method3_df_vis, aes(x = comparison, y = antigenic_score)) + geom_boxplot()
boxplot_method3_antigenic_scores

boxplot_method4_antigenic_scores <- ggplot(method4_df_vis, aes(x = comparison, y = antigenic_score)) + geom_boxplot()
boxplot_method4_antigenic_scores

boxplot_method5_antigenic_scores <- ggplot(method5_df_vis, aes(x = comparison, y = antigenic_score)) + geom_boxplot()
boxplot_method5_antigenic_scores

boxplot_method6_antigenic_scores <- ggplot(method6_df_vis, aes(x = comparison, y = antigenic_score)) + geom_boxplot()
boxplot_method6_antigenic_scores

# Perform Correlation Test
# Ultimately the Spearman's Correlation was used due to potential outliers 

# Perform the Pearson's Correlation
# Against Antigenic Distances
method1_pearsons <- cor.test(cartography_distances, method1_antigenic_scores)
method2_pearsons <- cor.test(cartography_distances, method2_antigenic_scores)
method3_pearsons <- cor.test(cartography_distances, method3_antigenic_scores)
method4_pearsons <- cor.test(cartography_distances, method4_antigenic_scores)
method5_pearsons <- cor.test(cartography_distances, method5_antigenic_scores)
method6_pearsons <- cor.test(cartography_distances, method6_antigenic_scores)
print("Method 1 Pearson's Results: ")
method1_pearsons
print("Method 2 Pearson's Results: ")
method2_pearsons
print("Method 3 Pearson's Results: ")
method3_pearsons
print("Method 4 Pearson's Results: ")
method4_pearsons
print("Method 5 Pearson's Results: ")
method5_pearsons
print("Method 6 Pearson's Results: ")
method6_pearsons
# Against mFRN Values
method1_pearsons_mfrn <- cor.test(neutralization_values, method1_antigenic_scores_mfrn)
method2_pearsons_mfrn <- cor.test(neutralization_values, method2_antigenic_scores_mfrn)
method3_pearsons_mfrn <- cor.test(neutralization_values, method3_antigenic_scores_mfrn)
method4_pearsons_mfrn <- cor.test(neutralization_values, method4_antigenic_scores_mfrn)
method5_pearsons_mfrn <- cor.test(neutralization_values, method5_antigenic_scores_mfrn)
method6_pearsons_mfrn <- cor.test(neutralization_values, method6_antigenic_scores_mfrn)
print("Method 1 Pearson's Results on mFRN values: ")
method1_pearsons_mfrn
print("Method 2 Pearson's Results on mFRN values: ")
method2_pearsons_mfrn
print("Method 3 Pearson's Results on mFRN values: ")
method3_pearsons_mfrn
print("Method 4 Pearson's Results on mFRN values: ")
method4_pearsons_mfrn
print("Method 5 Pearson's Results on mFRN values: ")
method5_pearsons_mfrn
print("Method 6 Pearson's Results on mFRN values: ")
method6_pearsons_mfrn

# Perform Spearman's Correlation Test 
# Against Antigenic Distance
method1_spearmans <- cor.test(cartography_distances, method1_antigenic_scores, method = 'spearman', exact = FALSE)
method2_spearmans <- cor.test(cartography_distances, method2_antigenic_scores, method = 'spearman', exact = FALSE)
method3_spearmans <- cor.test(cartography_distances, method3_antigenic_scores, method = 'spearman', exact = FALSE)
method4_spearmans <- cor.test(cartography_distances, method4_antigenic_scores, method = 'spearman', exact = FALSE)
method5_spearmans <- cor.test(cartography_distances, method5_antigenic_scores, method = 'spearman', exact = FALSE)
method6_spearmans <- cor.test(cartography_distances, method6_antigenic_scores, method = 'spearman', exact = FALSE)
print("Method 1 Spearman's Results: ")
method1_spearmans
print("Method 2 Spearman's Results: ")
method2_spearmans
print("Method 3 Spearman's Results: ")
method3_spearmans
print("Method 4 Spearman's Results: ")
method4_spearmans
print("Method 5 Spearman's Results: ")
method5_spearmans
print("Method 6 Spearman's Results: ")
method6_spearmans
# Against mFRN Values
method1_spearmans_mfrn <- cor.test(neutralization_values, method1_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
method2_spearmans_mfrn <- cor.test(neutralization_values, method2_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
method3_spearmans_mfrn <- cor.test(neutralization_values, method3_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
method4_spearmans_mfrn <- cor.test(neutralization_values, method4_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
method5_spearmans_mfrn <- cor.test(neutralization_values, method5_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
method6_spearmans_mfrn <- cor.test(neutralization_values, method6_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
print("Method 1 Spearman's Results on mFRN values: ")
method1_spearmans_mfrn
print("Method 2 Spearman's Results on mFRN values: ")
method2_spearmans_mfrn
print("Method 3 Spearman's Results on mFRN values: ")
method3_spearmans_mfrn
print("Method 4 Spearman's Results on mFRN values: ")
method4_spearmans_mfrn
print("Method 5 Spearman's Results on mFRN values: ")
method5_spearmans_mfrn
print("Method 6 Spearman's Results on mFRN values: ")
method6_spearmans_mfrn

# r.test to compare correlation values (the difference between two independent correlations, baseline vs method)
method1_corr_cartography <- method1_pearsons$estimate
method1_corr_mfrn <- method1_pearsons_mfrn$estimate
method2_corr_cartography <- method2_pearsons$estimate
method2_corr_mfrn <- method2_pearsons_mfrn$estimate
method3_corr_cartography <- method3_pearsons$estimate
method3_corr_mfrn <- method3_pearsons_mfrn$estimate
method4_corr_cartography <- method4_pearsons$estimate
method4_corr_mfrn <- method4_pearsons_mfrn$estimate

method1_rtest <- r.test(length(cartography_distances), method1_corr_cartography, 
                        r34 = method4_corr_cartography, pooled=TRUE, twotailed = TRUE)
method1_rtest_mfrn <- r.test(length(neutralization_values), method1_corr_mfrn, 
                        r34 = method4_corr_mfrn, pooled=TRUE, twotailed = TRUE)
method2_rtest <- r.test(length(cartography_distances), method2_corr_cartography, 
                        r34 = method4_corr_cartography, pooled=TRUE, twotailed = TRUE)
method2_rtest_mfrn <- r.test(length(neutralization_values), method2_corr_mfrn, 
                             r34 = method4_corr_mfrn, pooled=TRUE, twotailed = TRUE)
method3_rtest <- r.test(length(cartography_distances), method3_corr_cartography, 
                        r34 = method4_corr_cartography, pooled=TRUE, twotailed = TRUE)
method3_rtest_mfrn <- r.test(length(neutralization_values), method3_corr_mfrn, 
                             r34 = method4_corr_mfrn, pooled=TRUE, twotailed = TRUE)

print("Method 1 versus Method 4 (Baseline) R.Test for Cartography: ")
method1_rtest
print("Method 1 versus Method 4 (Baseline) R.Test for Neutralization Values: ")
method1_rtest_mfrn
print("Method 2 versus Method 4 (Baseline) R.Test for Cartography: ")
method2_rtest
print("Method 2 versus Method 4 (Baseline) R.Test for Neutralization Values: ")
method2_rtest_mfrn
print("Method 3 versus Method 4 (Baseline) R.Test for Cartography: ")
method3_rtest
print("Method 3 versus Method 4 (Baseline) R.Test for Neutralization Values: ")
method3_rtest_mfrn

# r.test to compare correlation values (the difference between two independent correlations, two best methods)
method2_corr_cartography <- method2_spearmans$estimate
method2_corr_mfrn <- method2_spearmans_mfrn$estimate
method5_corr_cartography <- method5_spearmans$estimate
method5_corr_mfrn <- method5_spearmans_mfrn$estimate

method_rtest <- r.test(length(cartography_distances), method2_corr_cartography, 
                        r34 = method5_corr_cartography, pooled=TRUE, twotailed = TRUE)
method_rtest_mfrn <- r.test(length(neutralization_values), method2_corr_mfrn, 
                             r34 = method5_corr_mfrn, pooled=TRUE, twotailed = TRUE)

print("Method 2 versus Method 5 R.Test for Cartography: ")
method_rtest
print("Method 2 versus Method 5 R.Test for Neutralization Values: ")
method_rtest_mfrn
