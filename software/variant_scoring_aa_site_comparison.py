#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import re
from matplotlib import pyplot as plt

columns = ['Location', 'Collection date', 'Accession ID', 'Pango lineage', 'AA Substitutions']

metadata = pd.read_csv(sys.argv[0], sep = '\t', usecols = columns)
tpSites = pd.read_csv(sys.arg[1], sep = ',')
tpSites_list = tpSites['tp_sites'].tolist()
tpSites_list = map(str, tpSites_list) 
output = sys.arg[2]
weights = pd.read_csv(sys.arg[3], sep = '\t')

def comparing_lists(list1):
# Function to find the intersection of two lists and returns the length of the intersection
    common = list(set(list1).intersection(set(top_pos_list)))
    return(len(common))

def conversion_position(x):
# Function to convert aa list in the metadata file to their position and original / changed aa
    subsList = x.split(',')
    mutationsSubsList = []
    spike_mutations = []
    for i in subsList:
        if "Spike" in i:
            spike_mutations.append(i)
    for i in spike_mutations:
        if i.__contains__("_"):
            pos = i[i.index("_"):]
            x = pos[1:]
            pos = "".join([i for i in pos if i.isdigit()])
            mutationsSubsList.append(pos)
        else:
            mutationsSubsList.append(None)
    return(mutationsSubsList)

def conversion_originalAA(x):
# Function to convert aa list in the metadata file to a list of the original aa 
    subsList = x.split(',')
    original_aaList = []
    spike_mutations = []
    for i in subsList:
        if "Spike" in i:
            spike_mutations.append(i)
    for i in spike_mutations:
        if i.__contains__("_"):
            pos = i[i.index("_"):]
            x = pos[1:]
            ogAAList = re.findall('([a-zA-Z]*)\d*.*', x)
            original_aaList.append(ogAAList[0])
    return(original_aaList)

def conversion_changedAA(x):
# Function to convert aa list in the metadata file to a list of the original aa 
    subsList = x.split(',')
    changed_aaList = []
    spike_mutations = []
    for i in subsList:
        if "Spike" in i:
            spike_mutations.append(i)
    for i in spike_mutations:
        if i.__contains__("_"):
            pos = i[i.index("_"):]
            x = pos[1:]
            cAAList = re.findall('\d([a-zA-Z]*)', x)
            cAA = cAAList[-1:]
            q = ''.join([str(i) for i in cAA])
            changed_aaList.append(q)
    return(changed_aaList)

def position_scores(input_df): # , output_df
# Function to add mutation scores to the metadata file based on amino acid changes
    input_df['AA Substitutions'] = input_df['AA Substitutions'].astype('string') 
    # Getting a list of all substitutions and amino acid changes per sequence
    substitutions = []
    original_AAs = []
    changed_AAs = []
    for i in input_df['AA Substitutions']:
        if len(i) == 0 or i == np.NaN:
            x = np.NaN
            substitutions.append(x)
            original_AAs.append(x)
            changed_AAs.append(x)
        else:
            substitution_pos = conversion_position(i)
            substitutions.append(substitution_pos)
            ogAA = conversion_originalAA(i)
            original_AAs.append(ogAA)
            cAA = conversion_changedAA(i)
            changed_AAs.append(cAA)
    if len(substitutions) == len(original_AAs) & len(original_AAs) == len(changed_AAs):
        input_df["Substitution_Positions"] = substitutions
        input_df["Original_aa"] = original_AAs
        input_df["Changed_aa"] = changed_AAs
    else:
        print("Length of spike substitution positions does not match the length of amino acid changes at those positions!")
        sys.exit()
    # Calculating the weight of spike amino acid changes
    cols = ['Substitution_Positions','Original_aa','Changed_aa']
    metadata_expanded = input_df.explode(cols)
    metadata_weights_expanded = pd.merge(metadata_expanded, weights, how = 'left', on = ['Original_aa', 'Changed_aa'])
    metadata_weights_expanded.drop('AA Substitutions', axis = 1, inplace = True)
    metadata_weights_expanded.drop_duplicates(inplace = True)
    # Summing the weight across each position per lineage
    # metadata_weights_expanded['weight_sum_per_position'] = metadata_weights_expanded.groupby(['Pango lineage','Substitution_Positions'])['Weight'].transform('sum')
    # metadata_weights_expanded['antigenic_score'] = metadata_weights_expanded.groupby('Substitution_Positions')['weight_sum_per_position'].transform('mean')
    # Averaging the antigenic weight per position across all lineages
    metadata_weights_expanded['Weight'] = metadata_weights_expanded['Weight'].fillna(0)
    metadata_weights_expanded['antigenic_score'] = metadata_weights_expanded.groupby('Substitution_Positions')['Weight'].transform('mean')
    metadata_weights_expanded['antigenic_score'] = metadata_weights_expanded['antigenic_score'].fillna(0)
    metadata_weights_expanded.drop(['Pango lineage', 'Original_aa', 'Changed_aa', 'Weight'], axis = 1, inplace = True) #, 'weight_sum_per_position'
    metadata_weights_expanded.drop_duplicates(inplace = True)
    metadata_weights_expanded['rank'] = metadata_weights_expanded['antigenic_score'].rank(ascending = False)
    metadata_weights_expanded = metadata_weights_expanded.sort_values('rank', ascending = True)
    # Adding column for TP Sites
    metadata_weights_expanded['tp_site'] = np.where(metadata_weights_expanded['Substitution_Positions'].isin(tpSites_list), "True", "False")
    return(metadata_weights_expanded)

# ### Selecting the Lineages of Comparison

# Number of sequences per lineage
freqs = metadata['Pango lineage'].value_counts()
seq_freq = pd.DataFrame(freqs)
seq_freq.reset_index(inplace = True)
seq_freq.columns = ['Pango lineage', 'sequence_frequencies']
seq_freq

# Selecting all pango lineages above the mean sequence frequency for manual curation
# freq_mean = seq_freq['sequence_frequencies'].mean()
freq_upper_quartile = seq_freq['sequence_frequencies'].quantile(q=0.95)
print(freq_upper_quartile)
seq_freq_upper_quartile = seq_freq[seq_freq['sequence_frequencies'] > freq_upper_quartile]
print(seq_freq_upper_quartile)
seq_freq_upper_quartile.to_csv(output + 'upper_quartile_all.csv')
# seq_freq_aboveMean

# Manually removing the sub lineages from the dataFrame:
remove = ['BA.1.1', 'None', 'AY.25.1', 'B.1.2', 'B.1', 'AY.4.2', 'B.1.1', 'AY.9.2',
         'AY.98.1', 'AY.4.2.2', 'AY.39.1', 'B.1.1.519', 'P.1.14', 'AY.4.2.1',
         'B.1.1.214', 'AY.4.6', 'AY.4.5', 'B.1.177.21', 'AY.119.2', 'B', 'AY.4', 'AY.43', 'AY.103', 'AY.44',
          'AY.122', 'AY.3', 'AY.25', 'AY.100', 'AY.29', 'AY.5', 'AY.39', 'AY.26', 'AY.98', 'AY.126', 'AY.47',
          'AY.121', 'AY.119', 'AY.20', 'AY.120', 'AY.125', 'AY.6', 'AY.42', 'AY.75', 'AY.127', 'AY.46.6', 'AY.99.2',
          'AY.27', 'AY.33', 'AY.46', 'AY.118', 'AY.23', 'AY.102', 'AY.112', 'AY.117', 'AY.9', 'AY.3.1', 'AY.36',
          'AY.129', 'AY.14', 'AY.54', 'AY.113', 'AY.74', 'AY.34']
seq_freq_filtered = seq_freq_upper_quartile[seq_freq_upper_quartile['Pango lineage'].isin(remove) == False]
seq_freq_filtered

# Establishing a list of lineages to look at average amino acid change antigenic score per site change
VOC_list = seq_freq_filtered['Pango lineage']
VOC_list

# #### Visualization of Selected Lineages and Frequency of Sequences

# Box plot of distribution
# fig = plt.figure(figsize = (8,8))
# plt.boxplot(seq_freq['sequence_frequencies'])
# plt.title('Distribution of the Number of Sequences \nper Pango Lineage', fontsize = 18, fontweight = 'bold')
# plt.xlabel('All Pango Lineages', fontsize = 15)
# plt.ylabel('Number of Sequences (1e6)', fontsize = 15)
# plt.yticks(fontsize = 15)
# plt.axhline(y = freq_upper_quartile, linestyle = '--', label = '95% Quartile \n(total of 80 pango lineages)')
# plt.legend(fontsize = 11)
# plt.savefig(output + 'sequence_frequency_box_plot.png')
# plt.close(fig = None)

# #### Creation of Final Dataframe

# Creating final dataframe to assign antigenic weights based on the selected lineages
VOC_metadata = metadata[metadata["Pango lineage"].isin(VOC_list)]
VOC_metadata.drop(['Accession ID', 'Collection date', 'Location'], axis = 1, inplace = True)
print(len(VOC_metadata))
VOC_metadata = VOC_metadata.dropna(subset=['AA Substitutions'])
print(len(VOC_metadata))
VOC_metadata.drop_duplicates(inplace = True)
print(len(VOC_metadata))
VOC_metadata

# ### Assigning Weights per Position

df = position_scores(VOC_metadata)
pd.DataFrame.head(df)
print(df.tp_site)
# df.to_csv(output + 'amino_acid_antigenic_weights.csv')
df

# #### Visualization

topDF = df[0:49]
tpDF = topDF[topDF["tp_site"] == "True"]
nontpDF = topDF[topDF["tp_site"] == "False"]

fig = plt.figure(figsize = (15,8))
#colors = ['silver' if (x == 'False') else 'darkorange' for x in df.tp_site]
# plt.bar(df.Substitution_Positions[0:49], df.antigenic_score[0:49], color = colors)
plt.bar(tpDF.Substitution_Positions, tpDF.antigenic_score, color = 'darkorange', label = 'TP Site')
plt.bar(nontpDF.Substitution_Positions, nontpDF.antigenic_score, color = 'steelblue', label = 'Non TP Site')
plt.title('Top 50 Antigenic Weights per Amino Acid Position in the Spike Protein \n (Averaged Across Known Variants of Concern)', fontsize = 18)
plt.xlabel('Amino Acid Position', fontsize = 16)
plt.xticks(fontsize = 14, rotation = 90)
plt.ylabel('Averaged Antigenic Weight', fontsize = 16)
plt.yticks(fontsize = 14)
# tpSites = ['False', 'True'], 
plt.legend(loc="best", fontsize = 14)
# plt.savefig(output + 'aaPosition_antigenic_weight_bar_plot.png')
plt.close(fig = None)

# #### Matching top changes to Pango lineages

# Matching the top 50 amino acid changes identified to the pango lineages that contain them
# spike_positions = []
# for x in metadata['AA Substitutions']:
#     positions = conversion_position(x)
#     spike_positions.append(positions)
# metadata['AA_Spike_Positions'] = spike_positions
metadata.drop(['Accession ID', 'Collection date', 'Location'], axis = 1, inplace = True)
metadata = metadata.dropna(subset=['AA Substitutions'])
metadata.drop_duplicates(inplace = True)
metadata['Substitution_Positions'] = metadata['AA Substitutions'].apply(conversion_position)
print(pd.DataFrame.head(metadata))

top_pos_list = list(topDF['Substitution_Positions'])
metadata['number_of_top_sites'] = metadata["Substitution_Positions"].apply(comparing_lists)
top_lineages = metadata[metadata['number_of_top_sites'] > 0]
top_lineages['averaged_number_of_top_sites'] = top_lineages.groupby('Pango lineage')['number_of_top_sites'].transform('mean')
top_lineages.drop(['Substitution_Positions', 'number_of_top_sites'], axis = 1, inplace = True)
top_lineages.drop_duplicates(inplace = True)
top_lineages['number_of_sites_rank'] = top_lineages['averaged_number_of_top_sites'].rank(ascending = False)
top_lineages = top_lineages.sort_values('number_of_sites_rank', ascending = True)
print(pd.DataFrame.head(top_lineages))
top_lineages.to_csv(output + 'lineages_number_of_top_sites.csv', sep = '\t', index = False)
