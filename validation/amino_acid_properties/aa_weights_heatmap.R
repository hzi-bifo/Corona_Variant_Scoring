#!/usr/bin/env Rscript

#### Amino Acid Weight Changes Script
#Author: Katrina Norwood
#Last Updated: 23/01/25

# Script to compare the weight changes of different amino acid changes that we
# have weights for.

# To run:
# 1. start in the /Corona_Variant_Scoring/ home directory
# 2. in command line run: ```Rscript validation/amino_acid_properties/aa_weights_heatmap.R```

library(dplyr)
library(ggplot2)
#install.packages('hrbrthemes')
library(hrbrthemes)

args = commandArgs(trailingOnly=TRUE)
#output <- args[1]
#mfrn <- read.csv(args[2], sep = "\t")
#antigenicScores_dir <- read.csv(args[3], sep = "\t")

output <- "validation/amino_acid_properties/"
aaweights <- read.csv("reference/antigenic_weights.csv", sep = "\t")
print(aaweights)

# Give extreme colors:
heatmap <- ggplot(aaweights, aes(Original_aa, Changed_aa, fill = Weight)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(
    x = "Original Amino Acid",   
    y = "Changed Amino Acid") +
  theme_ipsum(base_family = "Helvetica", base_size = 8) +
  theme(
    axis.title.x = element_text(hjust = 0.5),  
    axis.title.y = element_text(hjust = 0.5))
ggsave("aa_weight_changes_heatmap.pdf", heatmap, path = output, 
       width = 90, height = 90, units = c("mm"), dpi = 300) 
heatmap
