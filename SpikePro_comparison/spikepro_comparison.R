#### SpikePro Methods Comparison
#Author: Katrina Norwood
#Last Updated: 28/10/24

#Script to compare the SpikePro immune escape prediction algorithm to our methods 
#weighing amino acid changes across the spike protein.

# To run:
# 1. start in the /Corona_Variant_Scoring/ home directory
# 2. in command line run: ```Rscript SpikePro_comparison/spikepro_comparison.R```

library(dplyr)
library(ggplot2)
library(psych)

# Data Importation

output <- "SpikePro_comparison/"
mfrn <- read.csv("validation/mfrn_mabs_vocs.csv", sep = "\t")
antigenic_cartography_distances <- read.csv("validation/antigenic_cartography_distances.csv", sep = ",")
method_antigenicScores_dir <- "validation/methods_validation_data/WHO_ranked_weights_at_all_sites"
method_spikepro_df <- read.csv("SpikePro_comparison/spikepro_scores.csv", sep = ",")
  
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

## Data Cleaning for SpikePro Scores

# Getting lineages for comparison
# Here we don't take the BA.4 lineage since in the mFRN values its combined with BA.5 and we don't have that score for SpikePro
variants_pangoLineage <- c("B.1.1.7", "B.1.351", "P.1", "B.1.617.2", "BA.1", "BA.2", "BA.2.12.1", "XBB.1.16", "XBB.1.5", "EG.5", "EG.5.1")
method_spikepro <- subset(method_spikepro_df, select = c(Pango.lineage, score))
method_spikepro <- method_spikepro %>% rename(variant = Pango.lineage, antigenic_score = score)
method_spikepro <- method_spikepro[method_spikepro$variant %in% variants_pangoLineage, ]
method_spikepro

# Renaming lineages
method_spikepro$variant[method_spikepro$variant == 'B.1.1.7'] <- "Alpha"
method_spikepro$variant[method_spikepro$variant == 'B.1.351'] <- "Beta"
method_spikepro$variant[method_spikepro$variant == 'P.1'] <- "Gamma"
method_spikepro$variant[method_spikepro$variant == 'B.1.617.2'] <- "Delta"
method_spikepro

## Data Cleaning for VOC Antigenic Scores

create_dataframe <- function(indir){
  # Function for data import and cleaning for the antigenic scoring results
  
  antigenicScores_df <- data.frame(matrix(ncol = 4, nrow = 0))
  print(indir)
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
  print("antigenicScores_df_mean: ")
  print(antigenicScores_df_mean)
  
  # Renaming Pango lineages
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'B.1.1.7'] <- "Alpha"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'B.1.351'] <- "Beta"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'P.1'] <- "Gamma"
  antigenicScores_df_mean$Pango.lineage[antigenicScores_df_mean$Pango.lineage == 'B.1.617.2'] <- "Delta"
  antigenicScores_df_renamed <- antigenicScores_df_mean
  print("antigenicScores_df_mean: ")
  print(antigenicScores_df_mean)
  antigenicScores_df_renamed <- antigenicScores_df_renamed %>% rename_at('Pango.lineage', ~'variant')
  #antigenicScores_df_renamed <- antigenicScores_df_renamed[-c(10,11), ] 
  print(antigenicScores_df_renamed)
  
  return(antigenicScores_df_renamed)
}

method_antigenicScores <- create_dataframe(method_antigenicScores_dir)
method_antigenicScores

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
  print("Final df: ")
  print(df4)
  
  return(df4)
}

# Mean mFRN values used
method_df_vis <- data_visualization(method_antigenicScores, antigenic_cartography_distances, mfrn_df)
spikepro_df_vis <- data_visualization(method_spikepro, antigenic_cartography_distances, mfrn_df)

## Visualization

scatterplot <- function (method_df_vis, ylim, file_name, ylabel) {
  method_plt <- ggplot(method_df_vis, aes(x = value, y = antigenic_score, group = 1, colour = variant)) +
    theme_bw() +
    geom_point(size = 2.5) +
    ylim(0, ylim) +
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#332288", "#b66dff", "#999933", "slategray1", "violetred1", "springgreen", "slateblue", "tomato", "tan3", "lightpink4")) +   
    facet_wrap(~comparison, scales = "free_x") + #, strip.position = "bottom"
    labs(x = NULL, y = ylabel) +
    theme(text = element_text(family = "Helvetica"), strip.placement = "outside", 
          strip.text.x = element_text(size = 10), axis.text = element_text(size=10), 
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), 
          legend.title = element_text(size=8), legend.text = element_text(size=6)) +
    guides(col = guide_legend(ncol = 3))
  method_plt
  # Saving plot
  ggsave(method_plt, 
         filename = paste(output, file_name, sep = ""),
         device = "pdf",
         height = 45, width = 180, units = "mm")
  return(method_plt)
}

method_plt <- scatterplot(method_df_vis, 22, "antigenicScores_with_weights_at_all_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf", "Antigenic Alteration Scores")
spikepro_plt <- scatterplot(spikepro_df_vis, 18, "spikepro_scores_comparison_mFRNA_antigenicCartography_45x180.pdf", "SpikePro Scores")
method_plt
spikepro_plt

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

spikepro_split_cartography <- split_df_cartography(spikepro_df_vis)
spikepro_split_neutralization <- split_df_neutralization(spikepro_df_vis)
spikepro_cartography_distances <- as.matrix(spikepro_split_cartography['value']) # Only need to do this once as its the same across all methods
spikepro_neutralization_values <- as.matrix(spikepro_split_neutralization['value']) # Only need to do this once as its the same across all methods
spikepro_antigenic_scores <- as.matrix(spikepro_split_cartography['antigenic_score']) #as.vector
spikepro_antigenic_scores_mfrn <- as.matrix(spikepro_split_neutralization['antigenic_score']) #as.vector
spikepro_cartography_distances
spikepro_neutralization_values
spikepro_antigenic_scores
spikepro_antigenic_scores_mfrn

# Perform Spearman's Correlation Test 
# Against Antigenic Distance
method_spearmans <- cor.test(cartography_distances, method_antigenic_scores, method = 'spearman', exact = FALSE)
spikepro_spearmans <- cor.test(spikepro_cartography_distances, spikepro_antigenic_scores, method = 'spearman', exact = FALSE)
print("Method Spearman's Results on Antigenic Cartography Distances: ")
method_spearmans
print("SpikePro Spearman's Results on Antigenic Cartography Distances: ")
spikepro_spearmans
# Against mFRN Values
method_spearmans_mfrn <- cor.test(neutralization_values, method_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
spikepro_spearmans_mfrn <- cor.test(spikepro_neutralization_values, spikepro_antigenic_scores_mfrn, method = 'spearman', exact = FALSE)
print("Method Spearman's Results on mFRN Values: ")
method_spearmans_mfrn
print("SpikePro Spearman's Results on mFRN Values: ")
spikepro_spearmans_mfrn

knitr::knit_exit()