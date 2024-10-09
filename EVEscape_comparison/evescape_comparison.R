#### EVEscape Methods Comparison
#Author: Katrina Norwood
#Last Updated: 20/08/24

#Script to compare the evescape immune escpae prediction algorithm to our methods 
#weighing amino acid changes across the spike protein. 

library(dplyr)
library(ggplot2)
library(psych)

# Data Importation

args = commandArgs(trailingOnly=TRUE)
output <- args[1]
mfrn <- read.csv(args[2], sep = "\t")
antigenic_cartography_distances <- args[3]
## With weights at all spike protein sites
method_antigenicScores_dir <- args[4]
## EVEscape scores 
method_evescores <- args[5]

## Data Cleaning and Normalization for the mAb neutralization scores

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

## Data Cleaning for EVEscape Scores

# Getting lineages for comparison
variants_pangoLineage <- c("B.1.1.7", "B.1.351", "P.1", "B.1.429", "B.1.617.2", 
                           "C.37", "B.1.621", "BA.1", "BA.2", "BA.1.1", "BA.2.12.1", 
                           "BA.4", "BA.5", "XBB", "XBB.1.16", "XBB.1.5", "XBB.2.3", 
                           "EG.5", "EG.5.1", "JN.1")
method_evescores %>% mutate_at(c("EVEscape.score_sigmoid.avg", "EVEscape.score_sigmoid.max",
                                 "EVEscape.score_pos.avg", "EVEscape.score_pos.max"), as.numeric)
#method_evescores <- subset(method_evescores, select = c(pango_lineages, EVEscape.score_pos.avg))
#method_evescores <- method_evescores %>% rename(variant = pango_lineages, antigenic_score = EVEscape.score_pos.avg)
method_evescores <- subset(method_evescores, select = c(pango_lineages, EVEscape.score_sigmoid.avg))
method_evescores <- method_evescores %>% rename(variant = pango_lineages, antigenic_score = EVEscape.score_sigmoid.avg)
#method_evescores <- subset(method_evescores, select = c(pango_lineages, EVEscape.score_sigmoid.max)) 
#method_evescores <- method_evescores %>% rename(variant = pango_lineages, antigenic_score = EVEscape.score_sigmoid.max)
method_evescores <- method_evescores[method_evescores$variant %in% variants_pangoLineage, ]

# Taking the mean of lineages BA.4 and BA.5 to compare with mfrn values
ba4_ba5_evescape <- c("BA.4_BA.5")
value1_df <- method_evescores[method_evescores$variant == 'BA.4', ]
value2_df <- method_evescores[method_evescores$variant == 'BA.5', ]
ba45_mean_evescape <- mean(c(value1_df$antigenic_score[1], value2_df$antigenic_score[1]))
print(ba45_mean_evescape) # should be 98.31754 for averaged column, 
ba4_ba5_evescape <- c(ba4_ba5_evescape, ba45_mean_evescape)
method_evescores <- rbind(method_evescores, ba4_ba5_evescape)
print(method_evescores)
print(colnames(method_evescores))

# Renaming lineages
method_evescores$variant[method_evescores$variant == 'B.1.1.7'] <- "Alpha"
method_evescores$variant[method_evescores$variant == 'B.1.351'] <- "Beta"
method_evescores$variant[method_evescores$variant == 'P.1'] <- "Gamma"
method_evescores$variant[method_evescores$variant == 'B.1.617.2'] <- "Delta"
method_evescores$variant[method_evescores$variant == 'C.37'] <- "Lambda"
method_evescores$variant[method_evescores$variant == 'B.1.621'] <- "Mu"
method_evescores$variant[method_evescores$variant == 'B.1.429'] <- "Epsilon"

## Data Cleaning for VOC Antigenic Scores

create_dataframe <- function(indir){
  # Function for data import and cleaning for the antigenic scoring results
  
  antigenicScores_df <- data.frame(matrix(ncol = 4, nrow = 0))
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
  
  # Taking the mean of lineages BA.4 and BA.5 to compare with mfrn values
  ba4_ba5 <- c("BA.4_BA.5")
  ba45_mean <- mean(c(antigenicScores_df_mean$antigenic_score[8], antigenicScores_df_mean$antigenic_score[9]))
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

method_antigenicScores <- create_dataframe(method_antigenicScores_dir)

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
method_df_vis <- data_visualization(method_antigenicScores, antigenic_cartography_distances, mfrn_df)
evescape_df_vis <- data_visualization(method_evescores, antigenic_cartography_distances, mfrn_df)

## Visualization

scatterplot <- function (method_df_vis, ylim, file_name, ylabel) {
  method_plt <- ggplot(method_df_vis, aes(x = value, y = antigenic_score, group = 1, colour = variant)) +
    theme_bw() +
    geom_point(size = 2.5) +
    ylim(0, ylim) +
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#332288", "#b66dff", "#999933", "slategray1", "violetred1", "springgreen", "slateblue", "tomato", "tan3", "lightpink4")) +   
    facet_wrap(~comparison, scales = "free_x") + #, strip.position = "bottom"
    labs(x = NULL, y = ylabel) +
    theme(text = element_text(family = "Helvetica"), strip.placement = "outside", strip.text.x = element_text(size = 10), axis.text = element_text(size=10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.title = element_text(size=10), legend.text = element_text(size=8)) #strip.background = element_blank(),
  method_plt
  # Saving plot
  ggsave(method_plt, 
         filename = paste(output, file_name, sep = ""),
         device = "pdf",
         height = 90, width = 180, units = "mm")
  return(method_plt)
}

method_plt <- scatterplot(method_df_vis, 25, "antigenicScores_wit_weights_at_all_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_90x180_test.pdf", "Antigenic Alterations Score")
evescape_plt <- scatterplot(evescape_df_vis, 11, "evescape_scores_sigmoid_avg_comparison_mFRNA_antigenicCartography_90x180.pdf", "EVEscape Score Averaged")
method_plt
evescape_plt

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

# Method 
method_split_cartography <- split_df_cartography(method_df_vis)
method_split_neutralization <- split_df_neutralization(method_df_vis)
cartography_distances <- as.matrix(method_split_cartography['value']) # Only need to do this once as its the same across all methods
neutralization_values <- as.matrix(method_split_neutralization['value']) # Only need to do this once as its the same across all methods
method_antigenic_scores <- as.matrix(method_split_cartography['antigenic_score']) #as.vector
method_antigenic_scores_mfrn <- as.matrix(method_split_neutralization['antigenic_score']) #as.vector
cartography_distances
neutralization_values
method_antigenic_scores
method_antigenic_scores_mfrn

evescape_split_cartography <- split_df_cartography(evescape_df_vis)
evescape_split_neutralization <- split_df_neutralization(evescape_df_vis)
evescape_cartography_distances <- as.matrix(evescape_split_cartography['value']) # Only need to do this once as its the same across all methods
evescape_neutralization_values <- as.matrix(evescape_split_neutralization['value']) # Only need to do this once as its the same across all methods
evescape_antigenic_scores <- as.matrix(evescape_split_cartography['antigenic_score']) #as.vector
evescape_antigenic_scores_mfrn <- as.matrix(evescape_split_neutralization['antigenic_score']) #as.vector
evescape_cartography_distances
evescape_neutralization_values
evescape_antigenic_scores
evescape_antigenic_scores_mfrn

# Perform Spearman's Correlation Test 
# Against Antigenic Distance
method_spearmans <- cor.test(cartography_distances, method_antigenic_scores, method = 'spearman', exact = FALSE)
evescape_spearmans <- cor.test(evescape_cartography_distances, evescape_antigenic_scores, method = 'spearman', exact = FALSE)
print("Method Spearman's Results: ")
method_spearmans
print("Evescape Spearman's Results: ")
evescape_spearmans
# Against mFRN Values
method1_spearmans_mfrn <- cor.test(neutralization_values, method_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
evescape_spearmans_mfrn <- cor.test(evescape_neutralization_values, evescape_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
print("Method 1 Spearman's Results on mFRN values: ")
method1_spearmans_mfrn
print("Evescape Spearman's Results on mFRN values: ")
evescape_spearmans_mfrn

