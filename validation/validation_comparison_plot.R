#!/usr/bin/env Rscript

#### Methods Validation Script
#Author: Katrina Norwood
#Last Updated: 23/01/25

# Script to validate the different proposed methodologies for scoring SARS-CoV-2 
# circulating lineages based on their potential antigenicity. Compares these scores
# to both antigenic cartography distances as well as averaged mFRN values from
# previous literature.

# To run:
# 1. start in the /Corona_Variant_Scoring/ home directory
# 2. Confirm that the library writexl is installed otherwise, comment out line 245
# 3. in command line run: ```Rscript validation/validation_comparison.R```

library(dplyr)
library(ggplot2)
library(scales)
library(stats)
#install.packages('writexl')     
library(writexl)

## Importing Data Inputs - to run with args comment out the below lines:

#args = commandArgs(trailingOnly=TRUE)
#output <- args[1]
output <- "/Users/katrina/Corona_Variant_Scoring/validation/output/"
#mfrn <- read.csv(args[2], sep = "\t")
mfrn <- read.csv("validation/mfrn_mabs_vocs.csv", sep = "\t")
#antigenic_cartography_distances <- read.csv(args[3], sep = "\t")
antigenic_cartography_distances <- read.csv("validation/antigenic_cartography_distances.csv", sep = ",")
## Without weights at antigenic sites (method1) 
#method1_antigenicScores_dir <- args[4]
method1_antigenicScores_dir <- "validation/methods_validation_data/WHO_ranked_without_weights_at_antigenic_sites"
## With weights at all sites (method2)
#method2_antigenicScores_dir <- args[5]
method2_antigenicScores_dir <- "validation/methods_validation_data/WHO_ranked_weights_at_all_sites"
## Weights at antigenic sites (method3)
#method3_antigenicScores_dir <- args[6]
method3_antigenicScores_dir <- "validation/methods_validation_data/WHO_ranked_weights_at_antigenic_sites"
## Without weights at all sites - Baseline (method4)
#method4_antigenicScores_dir <- args[7]
method4_antigenicScores_dir <- "validation/methods_validation_data/WHO_ranked_noweights_all_sites"
## Weights without directionality at all sites
#method5_antigenicScores_dir <- args[8]
method5_antigenicScores_dir <- "validation/methods_validation_data/WHO_ranked_all_sites_reversible_weights"
## With new weights at all sites (using amino acid changes that occurred at least 3 times throughout the tree)
#method6_antigenicScores_dir <- args[9]
method6_antigenicScores_dir <- "validation/methods_validation_data/WHO_ranked_all_sites_new_weights"

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
  variants_pangoLineage <- c("B.1.1.7", "B.1.351", "P.1", "B.1.429", "B.1.617.2", 
                             "C.37", "B.1.621", "BA.1", "BA.2", "BA.1.1", "BA.2.12.1", 
                             "BA.4", "BA.5", "XBB", "XBB.1.16", "XBB.1.5", "XBB.2.3", 
                             "EG.5", "EG.5.1", "JN.1")
  for (file in list.files(path = indir, pattern = ".csv", all.files = TRUE, full.names = TRUE)) {
    print(file)
    df <- read.csv(file, sep = "\t")
    #print(nrow(df))
    variants_df <- df[df$Pango.lineage %in% variants_pangoLineage, ]
    #print(nrow(variants_df))
    antigenicScores_df <- rbind(antigenicScores_df, variants_df)}
  
  #print(antigenicScores_df)
  
  # Taking median antigenic score per pango lineage
  antigenicScores_df <- select(antigenicScores_df, c("Pango.lineage", "antigenic_score"))
  #antigenicScores_df_mean <- aggregate(.~Pango.lineage, antigenicScores_df, mean) # Here we take the mean score of the lineages
  antigenicScores_df_mean <- aggregate(.~Pango.lineage, antigenicScores_df, median) # Takes the median score of the lineages as there may have been high variance in each lineage based on the boxplot comparison
  #antigenicScores_df_mean %>% mutate_at(c("antigenic_score"), as.numeric)
  #print("antigenicScores_df_mean")
  #print(antigenicScores_df_mean)
  
  # Taking the mean of lineages BA.4 and BA.5 to compare with mfrn values
  ba4_ba5 <- c("BA.4_BA.5")
  value1_df <- antigenicScores_df_mean[antigenicScores_df_mean$Pango.lineage == 'BA.4', ]
  value2_df <- antigenicScores_df_mean[antigenicScores_df_mean$Pango.lineage == 'BA.5', ]
  ba45_mean <- mean(c(value1_df$antigenic_score[1], value2_df$antigenic_score[1]))
  #print("BA4_5 mean:")
  #print(ba45_mean)
  #ba45_mean <- mean(c(antigenicScores_df_mean$antigenic_score[8], antigenicScores_df_mean$antigenic_score[9]))
  ba4_ba5 <- c(ba4_ba5, ba45_mean)
  antigenicScores_df_mean <- rbind(antigenicScores_df_mean, ba4_ba5)
  #antigenicScores_df_mean <- antigenicScores_df_mean[-c(8,9), ]
  #print(antigenicScores_df_mean)
  #print(colnames(antigenicScores_df_mean))
  
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
  #print(antigenicScores_df_renamed)
  
  return(antigenicScores_df_renamed)
}

method1_antigenicScores <- create_dataframe(method1_antigenicScores_dir)
method2_antigenicScores <- create_dataframe(method2_antigenicScores_dir)
method3_antigenicScores <- create_dataframe(method3_antigenicScores_dir)
method4_antigenicScores <- create_dataframe(method4_antigenicScores_dir)
method5_antigenicScores <- create_dataframe(method5_antigenicScores_dir)
method6_antigenicScores <- create_dataframe(method6_antigenicScores_dir)

# Compiling antigenic scores and antigenic distances into one dataframe and saving
scoring_comparison <- antigenic_cartography_distances
method1_copy <- method1_antigenicScores
method2_copy <- method2_antigenicScores
method3_copy <- method3_antigenicScores
method4_copy <- method4_antigenicScores
method5_copy <- method5_antigenicScores
method6_copy <- method6_antigenicScores

colnames(scoring_comparison)[colnames(scoring_comparison) == "distance_from_D614G"] <- "cartography_distances"
colnames(method1_copy)[colnames(method1_copy) == "antigenic_score"] <- "method1_scores"
colnames(method2_copy)[colnames(method2_copy) == "antigenic_score"] <- "method2_scores"
colnames(method3_copy)[colnames(method3_copy) == "antigenic_score"] <- "method3_scores"
colnames(method4_copy)[colnames(method4_copy) == "antigenic_score"] <- "baseline_scores"
colnames(method5_copy)[colnames(method5_copy) == "antigenic_score"] <- "method5_scores"
colnames(method6_copy)[colnames(method6_copy) == "antigenic_score"] <- "method6_scores"

# Merging antigenic scores with antigenic cartography
dfs_to_join <- list(method1_copy, method2_copy, method3_copy, method5_copy, method6_copy, method4_copy)

for (df in dfs_to_join) {
  scoring_comparison <- left_join(scoring_comparison, df, by = "variant")
}

scoring_comparison_numeric <- scoring_comparison %>%
  mutate(across(-variant, ~as.numeric(as.character(.))))

scoring_comparison_rounded <- scoring_comparison_numeric %>%
  mutate(across(where(is.numeric), ~round(. , 2)))

print("Scoring Comparison - Antigenic Distances")
print(scoring_comparison_rounded)
write.csv(scoring_comparison_rounded, paste(output, "validation_scores.csv", sep = ""), row.names=FALSE)

# Merging antigenic scores with mFRN values
scoring_comparison_mfrn <- mfrn_df

for (df in dfs_to_join) {
  scoring_comparison_mfrn <- left_join(scoring_comparison_mfrn, df, by = "variant")
}

scoring_comparison_mfrn_numeric <- scoring_comparison_mfrn %>%
  mutate(across(-variant, ~as.numeric(as.character(.))))

scoring_comparison_mfrn_rounded <- scoring_comparison_mfrn_numeric %>%
  mutate(across(where(is.numeric), ~round(. , 2)))

print("Scoring Comparison - mFRN Values")
print(scoring_comparison_mfrn_rounded)
write.csv(scoring_comparison_mfrn_rounded, paste(output, "validation_scores_mFRN.csv", sep = ""), row.names=FALSE)

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
  
  #print("merged1: ")
  #print(merged1)
  #print("merged2: ")
  #print(merged2)
  
  df4 <- rbind(merged1, merged2)
  df4$value <- round(df4$value, digits = 2)
  df4$antigenic_score <- round(df4$antigenic_score, digits = 2)

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

write_xlsx(list(method1_validation = method1_df_vis, method2_validation = method2_df_vis,
                method3_validation = method3_df_vis, method4_validation = method4_df_vis,
                method5_validation = method5_df_vis, method6_validation = method6_df_vis),
           path = paste(output, "method_validation.xlsx", sep = ""), col_names = TRUE)

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

## Statistical Testing Between Antigenic Distances and Scores

# Testing for outliers within the data
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

# Perform Spearman's Correlation Test
# Against Antigenic Distance
method1_spearmans <- cor.test(scoring_comparison_rounded$cartography_distances, scoring_comparison_rounded$method1_scores, method = 'spearman', exact = FALSE)
method2_spearmans <- cor.test(scoring_comparison_rounded$cartography_distances, scoring_comparison_rounded$method2_scores, method = 'spearman', exact = FALSE)
method3_spearmans <- cor.test(scoring_comparison_rounded$cartography_distances, scoring_comparison_rounded$method3_scores, method = 'spearman', exact = FALSE)
method4_spearmans <- cor.test(scoring_comparison_rounded$cartography_distances, scoring_comparison_rounded$baseline_scores, method = 'spearman', exact = FALSE)
method5_spearmans <- cor.test(scoring_comparison_rounded$cartography_distances, scoring_comparison_rounded$method5_scores, method = 'spearman', exact = FALSE)
method6_spearmans <- cor.test(scoring_comparison_rounded$cartography_distances, scoring_comparison_rounded$method6_scores, method = 'spearman', exact = FALSE)
print("Method 1 Spearman's Results: ")
print(method1_spearmans)
print("Method 2 Spearman's Results: ")
print(method2_spearmans)
print("Method 3 Spearman's Results: ")
print(method3_spearmans)
print("Method 4 Spearman's Results: ")
print(method4_spearmans)
print("Method 5 Spearman's Results: ")
print(method5_spearmans)
print("Method 6 Spearman's Results: ")
print(method6_spearmans)

# p-value correction (Benjamini-Hochberg)
methods <- c("method1", "method2", "method3", "method4", "method5", "method6")
pvalues_cartography <- c(method1_spearmans$p.value, method2_spearmans$p.value, method3_spearmans$p.value,
                         method4_spearmans$p.value, method5_spearmans$p.value, method6_spearmans$p.value)
pvalues_cartography_adjusted <- data.frame(methods, pvalues_cartography)
pvalues_cartography_adjusted$pvalues_cartography_adjusted <- p.adjust(pvalues_cartography_adjusted$pvalues_cartography, method = "BH")
pvalues_cartography_adjusted

# Against mFRN Values
method1_spearmans_mfrn <- cor.test(scoring_comparison_mfrn_rounded$mfrn_means, scoring_comparison_mfrn_rounded$method1_scores, method = 'spearman', exact = FALSE)
method2_spearmans_mfrn <- cor.test(scoring_comparison_mfrn_rounded$mfrn_means, scoring_comparison_mfrn_rounded$method2_scores, method = 'spearman', exact = FALSE)
method3_spearmans_mfrn <- cor.test(scoring_comparison_mfrn_rounded$mfrn_means, scoring_comparison_mfrn_rounded$method3_scores, method = 'spearman', exact = FALSE)
method4_spearmans_mfrn <- cor.test(scoring_comparison_mfrn_rounded$mfrn_means, scoring_comparison_mfrn_rounded$baseline_scores, method = 'spearman', exact = FALSE)
method5_spearmans_mfrn <- cor.test(scoring_comparison_mfrn_rounded$mfrn_means, scoring_comparison_mfrn_rounded$method5_scores, method = 'spearman', exact = FALSE)
method6_spearmans_mfrn <- cor.test(scoring_comparison_mfrn_rounded$mfrn_means, scoring_comparison_mfrn_rounded$method6_scores, method = 'spearman', exact = FALSE)
print("Method 1 Spearman's Results on mFRN values: ")
print(method1_spearmans_mfrn)
print("Method 2 Spearman's Results on mFRN values: ")
print(method2_spearmans_mfrn)
print("Method 3 Spearman's Results on mFRN values: ")
print(method3_spearmans_mfrn)
print("Method 4 Spearman's Results on mFRN values: ")
print(method4_spearmans_mfrn)
print("Method 5 Spearman's Results on mFRN values: ")
print(method5_spearmans_mfrn)
print("Method 6 Spearman's Results on mFRN values: ")
print(method6_spearmans_mfrn)

# p-value correction (Benjamini-Hochberg)
pvalues_mfrn <- c(method1_spearmans_mfrn$p.value, method2_spearmans_mfrn$p.value, method3_spearmans_mfrn$p.value,
                         method4_spearmans_mfrn$p.value, method5_spearmans_mfrn$p.value, method6_spearmans_mfrn$p.value)
pvalues_mfrn_adjusted <- data.frame(methods, pvalues_mfrn)
pvalues_mfrn_adjusted$pvalues_mfrn_adjusted <- p.adjust(pvalues_mfrn_adjusted$pvalues_mfrn, method = "BH")
pvalues_mfrn_adjusted

# Wilcoxon test to compare correlation coefficients of selected method and baseline (method 4)

# Apply Min-Max Scaling to all numeric columns (except 'variant')
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))}

scaled_data <- scoring_comparison_rounded %>%
  mutate(across(where(is.numeric) & !variant, min_max_scale))

cartography_distances_scaled <- scaled_data$cartography_distances
method2_scores_scaled <- scaled_data$method2_scores
baseline_scores_scaled <- scaled_data$baseline_scores

# Calculate absolute deviations (how much the scores deviate from the cartography distances)
deviation_method2_scaled <- abs(cartography_distances_scaled - method2_scores_scaled)
deviation_baseline_scaled <- abs(cartography_distances_scaled - baseline_scores_scaled)

# Wilcoxon signed-rank test
wilcoxon_test <- wilcox.test(deviation_method2_scaled, deviation_baseline_scaled, alternative = "less", exact = FALSE, correct = TRUE, paired = FALSE) # Added paired = FALSE, though this should be default

# Summarize results
results <- list(
  "Test Statistic (Wilcoxon rank sum statistic)" = wilcoxon_test$statistic, # Wilcoxon rank sum statistic
  "P-Value" = wilcoxon_test$p.value,
  "Conclusion" = ifelse(wilcoxon_test$p.value < 0.05, 
                        "Method 2 is significantly better than the baseline (p < 0.05).", 
                        "No significant difference between Method 2 and the baseline.")
)
print(results)

knitr::knit_exit()
