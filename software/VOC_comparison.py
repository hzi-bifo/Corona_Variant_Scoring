#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

columns = ['Pango lineage', 'antigenic_score', 'WHO_label']
# voc_mutationScore_df = pd.read_csv(sys.argv[0], sep = '\t', usecols = columns)
# voc_df = pd.read_csv(sys.argv[1], sep = '\t')
# output = sys.argv[2]

voc_mutationScore_df = pd.read_csv('/home/knorwood/software/Corona_Variant_Scoring/test/output/variant_scoring_without_weights/antigenic_scores_ranked_without_weights_with_WHO.csv', sep = '\t', usecols = columns)
# voc_df = pd.read_csv('/home/knorwood/software/Corona_Variant_Scoring/reference/known_variants_of_concern.csv', sep = '\t')
output = '/home/knorwood/software/Corona_Variant_Scoring/test/output/variant_scoring_without_weights/'

# Getting mutation score for the VOCs
# lineages_df = df[["Pango lineage", "antigenic_score"]]
# lineages_df = lineages_df.drop_duplicates()
# voc_mutationScore_df = pd.merge(voc_df, lineages_df, how = "inner", on = "Pango lineage") # Want to keep only the VOCs that have mutation scores
# voc_mutationScore_df = voc_mutationScore_df.drop_duplicates()
# voc_mutationScore_df

# Calculating the threshold for finding other variants of interest, based on the average mutation score of the known VOCs
voc_threshold = voc_mutationScore_df['antigenic_score'].mean()
# Calculating mean of Omicron, Beta, and Gamma lineages (to serve as a more stringent threshold)
omicron_voc_mutationScore = voc_mutationScore_df[voc_mutationScore_df["WHO_label"] == "Omicron"]
beta_voc_mutationScore = voc_mutationScore_df[voc_mutationScore_df["WHO_label"] == "Beta"]
gamma_voc_mutationScore = voc_mutationScore_df[voc_mutationScore_df["WHO_label"] == "Gamma"]
top_voc_threshold = pd.concat([omicron_voc_mutationScore, beta_voc_mutationScore, gamma_voc_mutationScore])
top_voc_threshold_score = top_voc_threshold["antigenic_score"].mean()
# Saving thresholds to txt file
# output_txt = open(output + "voc_threshold.txt", "a")
# output_txt.write("\n")
# output_txt.write(str(datetime.now()))
# output_txt.write("\nVariant of Concern Antigenic Score Threshold: \n")
# output_txt.write(str(voc_threshold))
# output_txt.write("\nOmicron, Gamma, and Beta Average Mutation Score Threshold: \n")
# output_txt.write(str(top_voc_threshold_score))
# output_txt.write("\n")
# output_txt.close()

# plt_data = pd.merge(lineages_df, voc_df, how = "left", on = "Pango lineage")
plt_data = pd.DataFrame.copy(voc_mutationScore_df)
# plt_data['WHO_label'] = plt_data['WHO_label'].fillna("Non Variant of Concern") # Just in case but not necessary

# Density plot
# fig = sns.kdeplot(data = plt_data, x = "antigenic_score", hue = "WHO_label")
# fig.axvline(voc_threshold, color = 'r')
# fig.set(xlabel = "Antigenic Score", ylabel = "Density",
#         title = "Antigenic Score Densities for Known Variants of Concern")
# plt.savefig(output + "antigenic_score_VOC_density_plot.jpeg")

# Boxplot
fig = sns.boxplot(x = plt_data["WHO_label"], y = plt_data["antigenic_score"], palette = 'hls')
fig.axhline(voc_threshold, color = 'b')
fig.axhline(top_voc_threshold_score, color = 'r')
plot_xticks = np.arange(0, 7, 1)
#plot_xlabels = ('B.1.617.2\n(Delta)', 'B.1.1.7\n(Alpha)', 'Not\nVariant of Concern', 'B.1.1.529\n(Omicron)', 'B.1.351\n(Beta)', 'P.1\n(Gamma)', 'B.1.427/B.1.429\n(Epsilon)') # antigenic scores with weights based on TP Sites
plot_xlabels = ('B.1.1.529\n(Omicron)', 'Not\nVariant of Concern', 'B.1.617.2\n(Delta)', 'P.1\n(Gamma)', 'B.1.351\n(Beta)', 'B.1.1.7\n(Alpha)', 'B.1.427/B.1.429\n(Epsilon)') # analysis without weights
fig.set_xticks(plot_xticks, labels = plot_xlabels, fontsize = 15);
sns.set(rc = {"figure.figsize":(16, 16)});
plt.xlabel("Lineage (Including All Sublineages)", fontsize = 16, fontweight = 'bold');
plt.ylabel("Averaged Antigenic Score", fontsize = 16, fontweight = 'bold');
plt.yticks(fontsize = 15)
plt.title("Antigenic Score Averaged Across All Locations", fontsize = 18)
plt.savefig(output + "mutation_score_VOC_box_plot.jpeg")
