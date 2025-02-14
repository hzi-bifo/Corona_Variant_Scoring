#!/usr/bin/env python
# coding: utf-8

#### Amino Acid Properties Script
#Author: Katrina Norwood
#Last Updated: 14/02/25

# Script to create a boxplot comparing the antigenic weights and properties of given amino acid changes and to 
# explore the various property differences among the amino acid changes

import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from statsmodels.stats import multitest

output = "/validation/amino_acid_properties/"
aa_weights = pd.read_csv("/reference/antigenic_weights.csv", sep = "\t", header = 0)
aa_properties = pd.read_csv("/reference/amino_acid_properties.csv", sep = "\t", header = 0)
aa_weights_all = pd.read_csv("/reference/antigenic_weights_original.csv", sep = "\t", header = 0)

# ### Data Clean Up
# We used both the averaged antigenic weights (which were ultimately used in CoVerage) and also all the antigenic weights for ALL changes occurring throughout the influenza A antigenic tree. Here there are multiple amino acid changes occurring more than once but at different positions on the HA protein. For instance D2N and D31N, each with different weights.

# #### Averaged Antigenic Weights

# Creating a copy of the properties dataframes to be merged with the weights df
original_aa_properties = aa_properties.copy()
original_aa_properties.rename(columns={'amino_acid': 'original_amino_acid', 'code': 'Original_aa', 'hydropathy':'original_hydropathy',
                                      'charge':'original_charge', 'pKa_NH2':'original_pKa_NH2', 'pKa_COOH':'original_pKa_COOH',
                                      'pKR':'original_pKR', 'solubility':'original_solubility', 'molecular_weight':'original_molecular_weight'}, inplace=True)
changed_aa_properties = aa_properties.copy()
changed_aa_properties.rename(columns={'amino_acid': 'changed_amino_acid', 'code': 'Changed_aa', 'hydropathy':'changed_hydropathy',
                                      'charge':'changed_charge', 'pKa_NH2':'changed_pKa_NH2', 'pKa_COOH':'changed_pKa_COOH',
                                      'pKR':'changed_pKR', 'solubility':'changed_solubility','molecular_weight':'changed_molecular_weight'}, inplace=True)
#print("original_aa_properties: ")
#print(pd.DataFrame.head(original_aa_properties))
#print("\nchanged_aa_properties: ")
#print(pd.DataFrame.head(changed_aa_properties))

# Merging weights and properties dataframes together
aa_properties_weights = pd.merge(aa_weights, original_aa_properties, how="left", on=["Original_aa"])
aa_properties_weights = pd.merge(aa_properties_weights, changed_aa_properties, how="left", on=["Changed_aa"])
#print("aa_properties_weights: ")
#print(pd.DataFrame.head(aa_properties_weights))

# Creating categories columns for the boxplot visualization
def hydropathy(row):
    if row["original_hydropathy"] == row["changed_hydropathy"]:
        return "unchanged" 
    return f"{row['original_hydropathy']} - {row['changed_hydropathy']}"

def charge(row):
    charge_mapping = {
        "N": "Neutral",
        "+": "Positive",
        "-": "Negative"}
    
    if row["original_charge"] == row["changed_charge"]:
        return "unchanged" 
    
    original_charge = charge_mapping.get(row["original_charge"], row["original_charge"])
    changed_charge = charge_mapping.get(row["changed_charge"], row["changed_charge"])
    
    return f"{original_charge} - {changed_charge}"

def size(row):
    if row["original_molecular_weight"] == row["changed_molecular_weight"]:
        return "unchanged" 
    return row['original_molecular_weight'] - row['changed_molecular_weight']

aa_properties_weights["Hydropathy"] = aa_properties_weights.apply(hydropathy, axis=1)
aa_properties_weights["Charge"] = aa_properties_weights.apply(charge, axis=1)
aa_properties_weights["Molecular Weight Change"] = aa_properties_weights.apply(size, axis=1)

#print(pd.DataFrame.head(aa_properties_weights))
aa_properties_weights

# Creating a dataframe for boxplot visualization
df_hydropathy = aa_properties_weights[["Original_aa","Changed_aa", "Weight", "Hydropathy"]]
df_hydropathy["Molecular_Property"] = "Hydropathy"
df_hydropathy = df_hydropathy.rename(columns = {"Hydropathy": "Change"})

df_charge = aa_properties_weights[["Original_aa","Changed_aa", "Weight", "Charge"]]
df_charge["Molecular_Property"] = "Charge"
df_charge = df_charge.rename(columns = {"Charge": "Change"})

df_visualization = pd.concat([df_hydropathy, df_charge])

# Molecular weight dataframe for line plot
df_molecular_weight = aa_properties_weights[["Original_aa","Changed_aa", "Weight", "Molecular Weight Change"]]
df_molecular_weight.loc[df_molecular_weight["Molecular Weight Change"] == "unchanged", "Molecular Weight Change"] = 0

print(pd.DataFrame.tail(df_visualization))

# Creating a sub dataframe for boxplot visualization (To show large charge changes or neutral to changes)
df_charge["Change_Group"] = "unchanged"  # default value
df_charge.loc[df_charge["Change"].isin(["Neutral - Positive", "Positive - Neutral"]), "Change_Group"] = "Neutral / Positive"
df_charge.loc[df_charge["Change"].isin(["Neutral - Negative", "Negative - Neutral"]), "Change_Group"] = "Neutral / Negative"
df_charge.loc[df_charge["Change"].isin(["Positive - Negative", "Negative - Positive"]), "Change_Group"] = "Positive / Negative"

print(df_charge.tail())

# Saving the dataframes
aa_properties_weights.to_csv(output+"aa_properties_weights.tsv", sep = "\t", header = True)


# #### All Antigenic Weights
# Includes weights that occur across the whole influenza antigenic tree, so weights of the changes at different positions, rather than the averages of the specific amion acid change. 

# Merging weights and properties dataframes together
aa_properties_weights_all = pd.merge(aa_weights_all, original_aa_properties, how="left", on=["Original_aa"])
aa_properties_weights_all = pd.merge(aa_properties_weights_all, changed_aa_properties, how="left", on=["Changed_aa"])
#print("aa_properties_weights: ")
# print(pd.DataFrame.head(aa_properties_weights))
# print(len(aa_weights_all))
# print(len(aa_properties_weights_all))
# print(aa_properties_weights_all.columns)

aa_properties_weights_all["Hydropathy"] = aa_properties_weights_all.apply(hydropathy, axis=1)
aa_properties_weights_all["Charge"] = aa_properties_weights_all.apply(charge, axis=1)
aa_properties_weights_all["Molecular Weight Change"] = aa_properties_weights_all.apply(size, axis=1)
# print(pd.DataFrame.head(aa_properties_weights))
# print(len(aa_weights_all))
# print(len(aa_properties_weights_all))
# print(aa_properties_weights_all.columns)

# Creating a dataframe for boxplot visualization
df_hydropathy_all = aa_properties_weights_all[["Original_aa","Changed_aa", "Position", "Weight", "Hydropathy"]]
df_hydropathy_all["Molecular_Property"] = "Hydropathy"
df_hydropathy_all = df_hydropathy_all.rename(columns = {"Hydropathy": "Change"})

df_charge_all = aa_properties_weights_all[["Original_aa","Changed_aa", "Position", "Weight", "Charge"]]
df_charge_all["Molecular_Property"] = "Charge"
df_charge_all = df_charge_all.rename(columns = {"Charge": "Change"})

df_visualization_all = pd.concat([df_hydropathy_all, df_charge_all])

# Molecular weight dataframe for line plot
df_molecular_weight_all = aa_properties_weights_all[["Original_aa","Changed_aa", "Position", "Weight", "Molecular Weight Change"]]
df_molecular_weight_all.loc[df_molecular_weight_all["Molecular Weight Change"] == "unchanged", "Molecular Weight Change"] = 0

print(pd.DataFrame.tail(df_visualization_all))

# Creating a sub dataframe for boxplot visualization (To show large charge changes or neutral to changes)
df_charge_all["Change_Group"] = "unchanged"  # default value
df_charge_all.loc[df_charge_all["Change"].isin(["Neutral - Positive", "Positive - Neutral"]), "Change_Group"] = "Neutral / Positive"
df_charge_all.loc[df_charge_all["Change"].isin(["Neutral - Negative", "Negative - Neutral"]), "Change_Group"] = "Neutral / Negative"
df_charge_all.loc[df_charge_all["Change"].isin(["Positive - Negative", "Negative - Positive"]), "Change_Group"] = "Positive / Negative"

# Display the result
print(df_charge_all.tail())

aa_properties_weights_all.to_csv(output+"aa_properties_weights_all.tsv", sep = "\t", header = True)

# ### Data Exploration
# Confirming how many data point are in each section, making sure everything is put together correctly

# #### Averaged Antigenic Weights

# aa_properties_weights - used for analysis
print("Length of properties dataframe (", len(aa_weights),") matches the original list of weights: ", "Y" if len(aa_properties_weights) == len(aa_weights) else "N")

# HYDROPATHY
print("\nDataFrame Unique Properties - Hydropathy: \n")
print("Unique Hydropathic Properties:\n", aa_properties_weights["Hydropathy"].unique())
aa_properties_weights_hydropathy_unchanged = aa_properties_weights[aa_properties_weights["Hydropathy"] == "unchanged"]
aa_properties_weights_hydrophobic_hydrophilic = aa_properties_weights[aa_properties_weights["Hydropathy"] == "hydrophobic - hydrophilic"]
aa_properties_weights_moderate_hydrophobic = aa_properties_weights[aa_properties_weights["Hydropathy"] == "moderate - hydrophobic"]
aa_properties_weights_hydrophilic_hydrophobic = aa_properties_weights[aa_properties_weights["Hydropathy"] == "hydrophilic - hydrophobic"]
aa_properties_weights_moderate_hydrophilic = aa_properties_weights[aa_properties_weights["Hydropathy"] == "moderate - hydrophilic"]
aa_properties_weights_hydrophobic_moderate = aa_properties_weights[aa_properties_weights["Hydropathy"] == "hydrophobic - moderate"]
aa_properties_weights_hydrophilic_moderate = aa_properties_weights[aa_properties_weights["Hydropathy"] == "hydrophilic - moderate"]
print("\nNumber of AA Changes - Unchanged in Hydropathy (Total): ", len(aa_properties_weights_hydropathy_unchanged))
print("Number of AA Changes - Unchanged in Hydropathy (Hydrophobic - Hydrophobic): ", len(aa_properties_weights_hydropathy_unchanged[aa_properties_weights_hydropathy_unchanged["changed_hydropathy"] == "hydrophobic"]))
print("Number of AA Changes - Unchanged in Hydropathy (Hydrophilic - Hydrophilic): ", len(aa_properties_weights_hydropathy_unchanged[aa_properties_weights_hydropathy_unchanged["changed_hydropathy"] == "hydrophilic"]))
print("Number of AA Changes - Unchanged in Hydropathy (Moderate - Moderate): ", len(aa_properties_weights_hydropathy_unchanged[aa_properties_weights_hydropathy_unchanged["changed_hydropathy"] == "moderate"]), "\n")
print("\nNumber of AA Changes - hydrophobic - hydrophilic: ", len(aa_properties_weights_hydrophobic_hydrophilic))
print("Number of AA Changes - moderate - hydrophobic: ", len(aa_properties_weights_moderate_hydrophobic))
print("Number of AA Changes - hydrophilic - hydrophobic: ", len(aa_properties_weights_hydrophilic_hydrophobic))
print("Number of AA Changes - moderate - hydrophilic: ", len(aa_properties_weights_moderate_hydrophilic))
print("Number of AA Changes - hydrophobic - moderate: ", len(aa_properties_weights_hydrophobic_moderate))
print("Number of AA Changes - hydrophilic - moderate: ", len(aa_properties_weights_hydrophilic_moderate))
print("\nUnchanged Hydropathy + Changed Hydropathy = ", len(aa_properties_weights_hydropathy_unchanged) + (len(aa_properties_weights_hydrophobic_hydrophilic) +
                                                                                                          len(aa_properties_weights_moderate_hydrophobic) + len(aa_properties_weights_hydrophilic_hydrophobic) +
                                                                                                          len(aa_properties_weights_moderate_hydrophilic) + len(aa_properties_weights_hydrophobic_moderate) +
                                                                                                          len(aa_properties_weights_hydrophilic_moderate)))
# CHARGE
print("\nDataFrame Unique Properties - Charge: \n")
print("Unique Charge Properties:\n", aa_properties_weights["Charge"].unique())
aa_properties_weights_charge_unchanged = aa_properties_weights[aa_properties_weights["Charge"] == "unchanged"]
aa_properties_weights_neutral_negative = aa_properties_weights[aa_properties_weights["Charge"] == "Neutral - Negative"]
aa_properties_weights_negative_neutral = aa_properties_weights[aa_properties_weights["Charge"] == "Negative - Neutral"]
aa_properties_weights_negative_positive = aa_properties_weights[aa_properties_weights["Charge"] == "Negative - Positive"]
aa_properties_weights_neutral_positive = aa_properties_weights[aa_properties_weights["Charge"] == "Neutral - Positive"]
aa_properties_weights_positive_neutral = aa_properties_weights[aa_properties_weights["Charge"] == "Positive - Neutral"]
aa_properties_weights_positive_negative = aa_properties_weights[aa_properties_weights["Charge"] == "Positive - Negative"]
print("\nNumber of AA Changes - Unchanged in Charge (Total): ", len(aa_properties_weights_charge_unchanged))
print("Number of AA Changes - Unchanged in Charge (Neutral - Neutral): ", len(aa_properties_weights_charge_unchanged[aa_properties_weights_charge_unchanged["changed_charge"] == "N"]))
print("Number of AA Changes - Unchanged in Charge (Positive - Positive): ", len(aa_properties_weights_charge_unchanged[aa_properties_weights_charge_unchanged["changed_charge"] == "+"]))
print("Number of AA Changes - Unchanged in Charge (Negative - Negative): ", len(aa_properties_weights_charge_unchanged[aa_properties_weights_charge_unchanged["changed_charge"] == "-"]), "\n")
print("\nNumber of AA Changes - neutral - negative: ", len(aa_properties_weights_neutral_negative))
print("Number of AA Changes - negative - neutral: ", len(aa_properties_weights_negative_neutral))
print("Number of AA Changes - negative - positive: ", len(aa_properties_weights_negative_positive))
print("Number of AA Changes - neutral - positive: ", len(aa_properties_weights_neutral_positive))
print("Number of AA Changes - positive - neutral: ", len(aa_properties_weights_positive_neutral))
print("Number of AA Changes - positive - negative: ", len(aa_properties_weights_positive_negative))
print("\nUnchanged Charge + Changed Charge = ", len(aa_properties_weights_charge_unchanged) + (len(aa_properties_weights_neutral_negative) +
                                                                                                          len(aa_properties_weights_negative_neutral) + len(aa_properties_weights_negative_positive) +
                                                                                                          len(aa_properties_weights_neutral_positive) + len(aa_properties_weights_positive_neutral) +
                                                                                                          len(aa_properties_weights_positive_negative)))


# df_visualization - used for visualization
print("Length of df_visualization dataframe (", len(aa_weights),"x 2 ) matches the original list of weights: ", "Y" if len(df_visualization) == len(aa_weights)*2 else "N")

print("\nUnique Changes:\n", df_visualization["Change"].unique())
print("Negative to Positive Charge Changes: ")
print(len(df_visualization[df_visualization["Change"] == "Negative - Positive"]))

aa_properties_weights_positive_positive = aa_properties_weights.loc[(aa_properties_weights['original_charge'] == "+") & (aa_properties_weights['changed_charge'] == "+")]
aa_properties_weights_negative_negative = aa_properties_weights.loc[(aa_properties_weights['original_charge'] == "-") & (aa_properties_weights['changed_charge'] == "-")]
aa_properties_weights_neutral_neutral = aa_properties_weights.loc[(aa_properties_weights['original_charge'] == "N") & (aa_properties_weights['changed_charge'] == "N")]

# Finding Averaged Weights for the different Properties
print("Averaged Antigenic Weight (Negative - Neutral): ", aa_properties_weights_negative_neutral["Weight"].mean())
print("Averaged Antigenic Weight (Neutral - Negative): ", aa_properties_weights_neutral_negative["Weight"].mean())
print("Averaged Antigenic Weight (Negative - Positive): ", aa_properties_weights_negative_positive["Weight"].mean())
print("Averaged Antigenic Weight (Neutral - Positive): ", aa_properties_weights_neutral_positive["Weight"].mean())
print("Averaged Antigenic Weight (Positive - Neutral): ", aa_properties_weights_positive_neutral["Weight"].mean())
print("Averaged Antigenic Weight (Positive - Negative): ", aa_properties_weights_positive_negative["Weight"].mean())

print("Averaged Antigenic Weight (Positive - Positive): ", aa_properties_weights_positive_positive["Weight"].mean())
print("Averaged Antigenic Weight (Negative - Negative): ", aa_properties_weights_negative_negative["Weight"].mean())
print("Averaged Antigenic Weight (Neutral - Neutral): ", aa_properties_weights_neutral_neutral["Weight"].mean())


# ### Statistical Analysis

# #### Averaged Antigenic Weights

# Comparing Hydropathy to Antigenic Weight
p_values = []
properties = []

hydropathy_changed_data = aa_properties_weights[aa_properties_weights["Hydropathy"] != "unchanged"]
hydropathy_unchanged_data = aa_properties_weights[aa_properties_weights["Hydropathy"] == "unchanged"]

# Wilcoxon signed rank test for Hydropathy
for property in hydropathy_changed_data["Hydropathy"].unique():
    print(property)
    changed_data = hydropathy_changed_data[hydropathy_changed_data['Hydropathy'] == property]
    
    w_stat, p_value = stats.mannwhitneyu(hydropathy_unchanged_data["Weight"], changed_data["Weight"], alternative = 'two-sided')
    p_values.append(p_value)
    properties.append(property)
    
    if p_value < 0.05:
        print(f"Statistically significant difference between {property} and unchanged.\n")
    else:
        print(f"No significant difference between {property} between and unchanged.\n")

# Adding Bonferroni Correction to the p-values
hydropathy_stats = pd.DataFrame(list(zip(properties, p_values)), columns = ["hydropathy", "p-value"])
rejected, hydropathy_stats["p_value_adjusted"], alphacSidak, alphacBonf = multitest.multipletests(p_values, alpha = 0.05, method = "fdr_bh")
print(hydropathy_stats)
print("Adjusted alpha: ", alphacBonf)

# Comparing Charge to Antigenic Weight
p_values_charge = []
properties_charge = []

charge_changed_data = aa_properties_weights[aa_properties_weights["Charge"] != "unchanged"]
charge_unchanged_data = aa_properties_weights[aa_properties_weights["Charge"] == "unchanged"]

# Wilcoxon signed rank test for Hydropathy
for property in charge_changed_data["Charge"].unique():
    print(property)
    changed_data = charge_changed_data[charge_changed_data['Charge'] == property]
    
    w_stat, p_value = stats.mannwhitneyu(charge_unchanged_data["Weight"], changed_data["Weight"], alternative = 'two-sided')
    p_values_charge.append(p_value)
    properties_charge.append(property)
    
    if p_value < 0.05:
        print(f"Statistically significant difference between {property} and unchanged.\n")
    else:
        print(f"No significant difference between {property} between and unchanged.\n")

# Adding Bonferroni Correction to the p-values
charge_stats = pd.DataFrame(list(zip(properties_charge, p_values_charge)), columns = ["charge", "p-value"])
rejected, charge_stats["p_value_adjusted"], alphacSidak, alphacBonf = multitest.multipletests(p_values_charge, alpha = 0.05, method = "fdr_bh")
print(charge_stats)
print("Adjusted alpha: ", alphacBonf)

# Comparing Charge GROUPS to Antigenic Weight
p_values_charge_group = []
properties_charge_group = []

charge_changed_data = df_charge[df_charge["Change_Group"] != "unchanged"]
charge_unchanged_data = df_charge[df_charge["Change_Group"] == "unchanged"]

# Wilcoxon signed rank test for Hydropathy
for property_charge_group in charge_changed_data["Change_Group"].unique():
    print(property_charge_group)
    changed_data = charge_changed_data[charge_changed_data['Change_Group'] == property_charge_group]
    
    w_stat, p_value = stats.mannwhitneyu(charge_unchanged_data["Weight"], changed_data["Weight"], alternative = 'two-sided')
    p_values_charge_group.append(p_value)
    properties_charge_group.append(property_charge_group)
    
    if p_value < 0.05:
        print(f"Statistically significant difference between {property_charge_group} and unchanged.\n")
    else:
        print(f"No significant difference between {property_charge_group} between and unchanged.\n")

# Adding Bonferroni Correction to the p-values
charge_group_stats = pd.DataFrame(list(zip(properties_charge_group, p_values_charge_group)), columns = ["Change_Group", "p-value"])
rejected, charge_group_stats["p_value_adjusted"], alphacSidak, alphacBonf = multitest.multipletests(p_values_charge_group, alpha = 0.05, method = "fdr_bh")
print(charge_group_stats)
print("Adjusted alpha: ", alphacBonf)

# Antigenic Weight to Molecular Weight Correlation
stats.spearmanr(df_molecular_weight["Weight"], df_molecular_weight["Molecular Weight Change"], alternative='two-sided')

# #### All Antigenic Weights

# Comparing Hydropathy to Antigenic Weight
p_values_all = []
properties_all = []

hydropathy_changed_data_all = aa_properties_weights_all[aa_properties_weights_all["Hydropathy"] != "unchanged"]
hydropathy_unchanged_data_all = aa_properties_weights_all[aa_properties_weights_all["Hydropathy"] == "unchanged"]

# Wilcoxon signed rank test for Hydropathy
for property_all in hydropathy_changed_data_all["Hydropathy"].unique():
    print(property_all)
    changed_data_all = hydropathy_changed_data_all[hydropathy_changed_data_all['Hydropathy'] == property_all]
    
    w_stat_all, p_value_all = stats.mannwhitneyu(hydropathy_unchanged_data_all["Weight"], changed_data_all["Weight"], alternative = 'two-sided')
    p_values_all.append(p_value_all)
    properties_all.append(property_all)
    
    if p_value_all < 0.05:
        print(f"Statistically significant difference between {property_all} and unchanged.\n")
    else:
        print(f"No significant difference between {property_all} between and unchanged.\n")

# Adding Bonferroni Correction to the p-values
hydropathy_stats_all = pd.DataFrame(list(zip(properties_all, p_values_all)), columns = ["hydropathy", "p-value"])
rejected_all, hydropathy_stats_all["p_value_adjusted"], alphacSidak_all, alphacBonf_all = multitest.multipletests(p_values_all, alpha = 0.05, method = "fdr_bh")
print(hydropathy_stats_all)
print("Adjusted alpha: ", alphacBonf_all)

# Comparing Charge to Antigenic Weight
p_values_charge_all = []
properties_charge_all = []

charge_changed_data_all = aa_properties_weights_all[aa_properties_weights_all["Charge"] != "unchanged"]
charge_unchanged_data_all = aa_properties_weights_all[aa_properties_weights_all["Charge"] == "unchanged"]

# Wilcoxon signed rank test for Hydropathy
for property_all in charge_changed_data_all["Charge"].unique():
    print(property_all)
    changed_data_all = charge_changed_data_all[charge_changed_data_all['Charge'] == property_all]
    
    w_stat_all, p_value_all = stats.mannwhitneyu(charge_unchanged_data_all["Weight"], changed_data_all["Weight"], alternative = 'two-sided')
    p_values_charge_all.append(p_value_all)
    properties_charge_all.append(property_all)
    
    if p_value_all < 0.05:
        print(f"Statistically significant difference between {property_all} and unchanged.\n")
    else:
        print(f"No significant difference between {property_all} between and unchanged.\n")

# Adding Bonferroni Correction to the p-values
charge_stats_all = pd.DataFrame(list(zip(properties_charge_all, p_values_charge_all)), columns = ["charge", "p-value"])
rejected_all, charge_stats_all["p_value_adjusted"], alphacSidak_all, alphacBonf_all = multitest.multipletests(p_values_charge_all, alpha = 0.05, method = "fdr_bh")
print(charge_stats_all)
print("Adjusted alpha: ", alphacBonf_all)

# Comparing Charge GROUPS to Antigenic Weight
p_values_charge_group_all = []
properties_charge_group_all = []

charge_changed_data_all = df_charge_all[df_charge_all["Change_Group"] != "unchanged"]
charge_unchanged_data_all = df_charge_all[df_charge_all["Change_Group"] == "unchanged"]

# Wilcoxon signed rank test for Charge
for property_charge_group in charge_changed_data_all["Change_Group"].unique():
    print(property_charge_group)
    changed_data_all = charge_changed_data_all[charge_changed_data_all['Change_Group'] == property_charge_group]
    
    w_stat_all, p_value_all = stats.mannwhitneyu(charge_unchanged_data_all["Weight"], changed_data_all["Weight"], alternative = 'two-sided')
    p_values_charge_group_all.append(p_value_all)
    properties_charge_group_all.append(property_charge_group)
    
    if p_value_all < 0.05:
        print(f"Statistically significant difference between {property_charge_group} and unchanged.\n")
    else:
        print(f"No significant difference between {property_charge_group} between and unchanged.\n")

# Adding Bonferroni Correction to the p-values
charge_group_stats_all = pd.DataFrame(list(zip(properties_charge_group_all, p_values_charge_group_all)), columns = ["Change_Group", "p-value"])
rejected, charge_group_stats_all["p_value_adjusted"], alphacSidak_all, alphacBonf_all = multitest.multipletests(p_values_charge_group_all, alpha = 0.05, method = "fdr_bh")
print(charge_group_stats_all)
print("Adjusted alpha: ", alphacBonf_all)

# Antigenic Weight to Molecular Weight Correlation
stats.spearmanr(df_molecular_weight_all["Weight"], df_molecular_weight_all["Molecular Weight Change"], alternative='two-sided')

# ### Data Visualization
# Once again we did data visuals for both the average antigenic weights per amino acid change and the weights for ALL amino acid changes across the influenza antigenic tree. 

# #### Averaged Antigenic Weights

color_blind_palette = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 
                       "#999999", "#000000", "#F1A4D1", "#A9D18E", "#FF6A13", "#D72631"]

sns.set_theme(style = "whitegrid")
sns.boxplot(x = 'Molecular_Property', y = 'Weight', data = df_visualization, hue = 'Change', 
            showfliers = False, palette = color_blind_palette, dodge = True, boxprops = dict(alpha = 0.5))

ax = sns.swarmplot(x = 'Molecular_Property', y = 'Weight', data = df_visualization, hue = 'Change', 
              dodge = True, palette = color_blind_palette, size = 5, edgecolor = 'black')
handles, labels = ax.get_legend_handles_labels()

plt.xticks(ticks=range(len(df_visualization['Molecular_Property'].unique())), 
           labels=df_visualization['Molecular_Property'].unique(), rotation=45)

plt.legend(handles[0:13], labels[0:13], title = "Moleculare Change", bbox_to_anchor = (1.05, 1), loc = 'upper left')

plt.show()

# Boxplot looking at only Charge changes
sns.set_theme(style = "whitegrid")
sns.boxplot(x = 'Molecular_Property', y = 'Weight', data = df_charge, hue = 'Change_Group', 
            showfliers = False, palette = color_blind_palette, dodge = True, boxprops = dict(alpha = 0.5))

ax = sns.swarmplot(x = 'Molecular_Property', y = 'Weight', data = df_charge, hue = 'Change_Group', 
              dodge = True, palette = color_blind_palette, size = 5, edgecolor = 'black')
handles, labels = ax.get_legend_handles_labels()

plt.xticks(ticks=range(len(df_charge['Molecular_Property'].unique())), 
           labels=df_charge['Molecular_Property'].unique(), rotation=45)

plt.legend(handles[0:4], labels[0:4], title = "Moleculare Change", bbox_to_anchor = (1.05, 1), loc = 'upper left')

plt.show()

# Lineplot comparing antigenic weight to molecular weight
sns.set_theme(style = "whitegrid")
sns.lineplot(data = df_molecular_weight, x = 'Weight', y = 'Molecular Weight Change')

# #### All Antigenic Weights

sns.set_theme(style = "whitegrid")
sns.boxplot(x = 'Molecular_Property', y = 'Weight', data = df_visualization_all, hue = 'Change', 
            showfliers = False, palette = color_blind_palette, dodge = True, boxprops = dict(alpha = 0.5))

ax = sns.swarmplot(x = 'Molecular_Property', y = 'Weight', data = df_visualization_all, hue = 'Change', 
              dodge = True, palette = color_blind_palette, size = 5, edgecolor = 'black')
handles, labels = ax.get_legend_handles_labels()

plt.xticks(ticks=range(len(df_visualization_all['Molecular_Property'].unique())), 
           labels=df_visualization_all['Molecular_Property'].unique(), rotation=45)

plt.legend(handles[0:13], labels[0:13], title = "Moleculare Change", bbox_to_anchor = (1.05, 1), loc = 'upper left')

plt.show()

import matplotlib.patches as mpatches

# Boxplot looking at only Charge changes
cm = 1/2.54 
plt.figure(figsize = (9*cm, 9*cm))
order = ['unchanged', 'Neutral / Negative', 'Neutral / Positive', 'Positive / Negative']

sns.set_theme(style = "whitegrid")
sns.boxplot(x = 'Change_Group', y = 'Weight', data = df_charge_all, 
            showfliers = False, palette = color_blind_palette, dodge = True, boxprops = dict(alpha = 0.5), 
            order = order)

ax = sns.stripplot(x = 'Change_Group', y = 'Weight', data = df_charge_all, # Can also do a swarmplot() here
             dodge = True, palette = color_blind_palette, size = 5, edgecolor = 'black')

plt.xticks(rotation=45)
ax.set_xticks([])
ax.set_xlabel("")

legend_patches = [mpatches.Patch(color=color_blind_palette[i], label=order[i]) for i in range(len(order))]
plt.legend(handles=legend_patches, title="Charge Change", bbox_to_anchor=(1.05, 1))

# Plotting adjust p-values over the boxes
p_values = charge_group_stats_all["p_value_adjusted"].values  # Ensure it has 3 values

x_positions = [1, 2, 3]  # 'Neutral / Negative', 'Neutral / Positive', 'Positive / Negative'
y_max = df_charge_all["Weight"].max() * 1.1  

for i in range(3):
    ax.text(x_positions[i], y_max, f"p = {p_values[i]:.3g}", ha='center', fontsize=10, color='black')

# Saving plot
plt.savefig(output+'all_amino_acid_weights_grouped_molecular_change_90x90.pdf', dpi=300, bbox_inches='tight')
plt.show()

print(len(df_charge_all[df_charge_all["Change_Group"] == "unchanged"]))
print(len(df_charge_all[df_charge_all["Change_Group"] == "Neutral / Negative"]))
print(len(df_charge_all[df_charge_all["Change_Group"] == "Neutral / Positive"]))
print(len(df_charge_all[df_charge_all["Change_Group"] == "Positive / Negative"]))

# Lineplot comparing antigenic weight to molecular weight
sns.set_theme(style = "whitegrid")
sns.lineplot(data = df_molecular_weight_all, x = 'Weight', y = 'Molecular Weight Change')