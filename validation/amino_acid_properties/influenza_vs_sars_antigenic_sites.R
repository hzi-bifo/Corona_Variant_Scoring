#!/usr/bin/env Rscript

#### Influenza vs. SARS-CoV-2 antigenic amio acid changes
#Author: Katrina Norwood
#Last Updated: 11/02/25

# Script to compare (using a Fishers Exact Test and ROC Curve) amino acid changes
# in influenza that are antigenic, to those in SARS-CoV-2 that are found to be 
# antigenic in literature, to see any similarities between the two.

# To run:
# 1. start in the /Corona_Variant_Scoring/ home directory
# 2. in command line run: ```Rscript validation/amino_acid_properties/influenza_vs_sars_antigenic_sites.R```

library(dplyr)
library(stats)
library(ggplot2)
#install.packages("pROC")
library(pROC)
#install.packages("plotROC")
library(plotROC)

## Importing Data Inputs 

output <- "validation/amino_acid_properties/" 
aa_changes_master <- read.csv("validation/amino_acid_properties/viral_amino_acid_properties.csv", sep = ",")

## Data Cleaning & Contingency Table

# Making the comparison columns factors
aa_changes_master <- aa_changes_master %>% mutate_at(c('Influenza', 'SARS.CoV.2'), as.factor)

# Defining a contingency table for the fisher's exact test
contingency_table <- table(aa_changes_master$Influenza, aa_changes_master$SARS.CoV.2)
contingency_table

# Creating a dataframe for the ROC curve analysis
aa_changes_comparison <- aa_changes_master[c("Weight", "Influenza", "SARS.CoV.2")]

# Creating a dataframe for the F1-score
# Here we turn the NA values in the weight column to 0, as the aa changes are not
# considered antigenic
aa_changes_f1score <- aa_changes_comparison %>% replace(is.na(.), 0)

## Fishers Exact Test

fisher.test(contingency_table)

## ROC Comparison

roc_curve_sars <- roc(as.numeric(aa_changes_f1score$SARS.CoV.2), aa_changes_f1score$Weight)

auc_value <- auc(roc_curve_sars)
print(paste("AUC:", round(auc_value, 2)))

## Fishers Exact Test and F1-Score at Different Thresholds (Without Shuffling)

threshold_fisher <- function(weights, df) { #actual_labels
  # Calculates the Fishers Exact Test at different weight thresholds
  weights <- weights[is.finite(weights)]  
  thresholds <- seq(min(weights, na.rm = TRUE), max(weights, na.rm = TRUE), length.out = 100)
  #print(thresholds)
  p_values <- sapply(thresholds, function(t) {
    #df_threshold <- df %>% filter(df$Weight > t)
    df_temp <- df
    df_temp$Influenza <- as.logical(df_temp$Influenza)
    df_temp$Influenza <- ifelse(df_temp$Weight <= t, FALSE, df_temp$Influenza)
    #print(head(df_temp))
    table_values <- table(df_temp$Influenza, df_temp$SARS.CoV.2)
    
    if (nrow(table_values) == 2 && ncol(table_values) == 2) {
      return(fisher.test(table_values)$p.value)
    } else {
      return(NA)}})
  
  f1_scores <- sapply(thresholds, function (t) {
    df_temp2 <- df
    df_temp2$Influenza <- as.logical(df_temp2$Influenza)
    df_temp2$Influenza <- ifelse(df_temp2$Weight <= t, FALSE, df_temp2$Influenza)
    #print(head(df_temp2))
    #print(t)
    pred <- df_temp2$SARS.CoV.2
    actual <- df_temp2$Influenza
    
    tp <- sum(pred == TRUE & actual == TRUE)
    fp <- sum(pred == TRUE & actual == FALSE)
    fn <- sum(pred == FALSE & actual == TRUE)
    
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    
    if (is.nan(precision) || is.nan(recall)) {
      f1_score <- 0} else {
        f1_score <- 2 * (precision * recall) / (precision + recall)
      } 
    return(f1_score)})
    #fisher.test(table_values)$p.value})
  return(list(thresholds, p_values, f1_scores))}

fishers_threshold_results <- threshold_fisher(aa_changes_f1score$Weight, aa_changes_f1score)
fishers_threshold_df <- data.frame(fishers_threshold_results[[1]], fishers_threshold_results[[2]], 
                                   fishers_threshold_results[[3]])

## Visualizations

# Fishers Exact Test and F1 Scores at Different Thresholds Without Shuffling
fishers_threshold_df_plot <- na.omit(fishers_threshold_df)

pvalue_threshold_plt <- ggplot(fishers_threshold_df_plot, aes(x = fishers_threshold_results..1..)) +
  geom_line(aes(y = fishers_threshold_results..2.., color = "P-Value"), linewidth = 1.2) +  
  geom_line(aes(y = fishers_threshold_results..3.., color = "F1 Score"), 
            linewidth = 1.2, linetype = "solid", show.legend = TRUE) +  
  scale_y_continuous(
    name = "P-Value",
    sec.axis = sec_axis(~ . * max(fishers_threshold_df_plot$fishers_threshold_results..3..) / max(fishers_threshold_df_plot$fishers_threshold_results..2..), 
                        name = "F1 Score")
  ) +
  scale_color_manual(values = c("P-Value" = "#1C86EE", "F1 Score" = "orange")) +
  labs(x = "Antigenic Weight Threshold", color = "Metric") +  # Removed title
  theme_minimal(base_size = 10) +  
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 10))
pvalue_threshold_plt

#ggsave(paste(output, "sarscov2_v_influenza_fishers_and_f1score_nonshuffled_90x90.pdf", sep = ""), plot = fishers_plt, 
#       width = 90, height = 90, units = c("mm"), dpi = 300)

# ROC Curve
# Extract ROC curve data (fpr, tpr, thresholds)
plot.roc(roc_curve_sars, col = "blue")
text(x = 0.6, y = 0.2, labels = paste("AUC = ", round(auc_value, 3)), 
     col = "red", cex = 1.2)

knitr::knit_exit()
