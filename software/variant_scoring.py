#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import numpy as np
import re
import plotly.express as px
from itertools import chain
from time import process_time

columns = ['Location', 'Collection date', 'Accession ID', 'Pango lineage', 'AA Substitutions']
metadata = pd.read_csv(sys.argv[0], sep = '\t', usecols = columns)
tpSites = pd.read_csv(sys.arg[1], sep = ',')
output = sys.arg[2]
weights = pd.read_csv(sys.arg[3], sep = '\t')
voc_df = pd.read_csv(sys.arg[4], sep = '\t')

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

def mutation_scores(input_df, tpSites_list): # , output_df
# Function to add mutation scores to the metadata file based on amino acid changes
    print("Dropping NaN rows ...")
    input_df = input_df.dropna(subset = ['AA Substitutions'])

    # Converting AA Substitutions a list of aa positions for easier comprehension
    print("Converting AA Substitutions to a string...")
    input_df['AA Substitutions'] = input_df['AA Substitutions'].astype('string')
    print("Converting AA Substitutions to a list of aa positions...")
    print("Substitution_Positions")
    substitutions = input_df['AA Substitutions'].apply(conversion_position)
    input_df["Substitution_Positions"] = substitutions
    print("Original Amino Acids")
    original_AAs = input_df['AA Substitutions'].apply(conversion_originalAA)
    input_df["Original_aa"] = original_AAs
    print("Changed Amino Acids")
    changed_AAs = input_df['AA Substitutions'].apply(conversion_changedAA)
    input_df["Changed_aa"] = changed_AAs

    # Expanding list of original and changed amino acids 
    try:
        cols = ['Substitution_Positions','Original_aa','Changed_aa']
        metadata_expanded = input_df.explode(cols)
        # Calculating the antigenic weight scores dependent on the amino acid changes (from a reference file)
        metadata_filtered = metadata_expanded[metadata_expanded['Substitution_Positions'].isin(tpSites_list)]
        metadata_filtered_weights = pd.merge(metadata_filtered, weights, how='left', on=['Original_aa', 'Changed_aa'])

    except: 
#         input_df.drop()
        metadata_filtered_weights = input_df
        metadata_filtered_weights['Substitution_Positions'] = np.NAN
        metadata_filtered_weights['Weight'] = np.NAN
        
    print(pd.DataFrame.head(metadata_filtered_weights))
    return(metadata_filtered_weights)

# Creating list of true positive antigenic sites for comparison
print("Creating list of known antigenic sites...")
tpSites_list = tpSites['tp_sites'].values.tolist()
tpSites_list = [str(item) for item in tpSites_list]
print("Done")

# Mutation Scores for metadata parts
print("Calculating Mutation Scores....")
metadata_filtered_weights = pd.DataFrame()
metadata_filtered_combinedWeights = pd.DataFrame()
n = 1500
t = 0
df_chunks = np.array_split(metadata, n)
tstart = process_time()
for chunk in df_chunks:
    print(chunk.shape[0])
    t = t+1
    print(t)
    mutation_df = mutation_scores(chunk, tpSites_list) 
    metadata_filtered_weights = pd.concat([mutation_df, metadata_filtered_weights])
    print("Output DataFrame Shape: ")
    print(metadata_filtered_weights.shape)
    print(pd.DataFrame.head(metadata_filtered_weights))
tstop = process_time()
print("Process Time: ", tstop - tstart)
print("FINAL METADATA_FILTERED_WEIGHTS DF: ")
pd.DataFrame.head(metadata_filtered_weights)

# Summing Weights by accession ID
metadata_filtered_combinedWeights = metadata_filtered_weights.groupby(metadata_filtered_weights['Accession ID']).aggregate({'Weight':'sum'})
df = pd.merge(metadata, metadata_filtered_combinedWeights, how = 'left', left_on = 'Accession ID', right_index = True)
df['Weight'].fillna(0, inplace = True)

print("Creating final dataframe...")
# Calculating final mutation score by taking the mean of the lineages that occur per location
# df['mutation_score'] = df.groupby(['Location','Pango lineage'])['Weight'].transform('mean')

# Calculating final antigenic score by averaging the scores across the lineages
df['antigenic_score'] = df.groupby('Pango lineage')['Weight'].transform('mean')
df_final = df[["Accession ID","Collection date", "Location", "Pango lineage", "AA Substitutions", "Weight", "antigenic_score"]]

print("Saving dataframe...")
df_final.to_csv(output + "antigenic_scores_all.csv", sep = '\t', index = False, header = True)
# Merging lineages per location and saving as a separate file
# df_filtered = df[['Location', 'Pango lineage', 'mutation_score']]
# df_filtered = df_filtered.drop_duplicates()
# df_filtered.to_csv(output + "mutation_scores.csv", sep = '\t', index = False, header = True)
print("Antigenic Scores COMPLETE")

# Ranking pango lineages based on average mutation score across all the sequences
df_ranked = df_final[["Pango lineage", "antigenic_score"]]
df_ranked = df_ranked.drop_duplicates()
df_ranked.dropna(axis = 0, inplace = True)
df_ranked.sort_values(by = ['antigenic_score'], ascending = False, inplace = True)
df_ranked.reset_index(inplace = True)
df_ranked["rank"] = df_ranked['antigenic_score'].rank(ascending = False)
df_ranked = df_ranked[["Pango lineage", "antigenic_score", "rank"]]
df_merged_ranked = df_ranked.merge(voc_df, how = "left", on = "Pango lineage")
df_merged_ranked['WHO_label'] = df_merged_ranked['WHO_label'].fillna("Non Variant of Concern")
df_merged_ranked.to_csv(output + "antigenic_scores_ranked_with_WHO.csv", sep = '\t', index = False, header = True)
print("Ranked Antigenic Scores COMPLETE")


# ### Visualization

# # BASED ON MUTATION SCORE THRESHOLD
# # Using the averaged mutation score per country for each lineage
# df_vis = df_final[['Location', 'Pango lineage', 'antigenic_score']]
# df_vis["Continent"] = df_vis["Location"].apply(lambda x: x.split("/")[0])
# df_vis["Country"] = df_vis["Location"].apply(lambda x: x.split("/")[1])
# country = []
# continent = []
# for item in df_vis["Country"]:
#     x = item.lstrip()
#     y = x.rstrip()
#     country.append(y)
# for item in df_vis["Continent"]:
#     x = item.lstrip()
#     y = x.rstrip()
#     continent.append(y)
# df_vis['Country'] = country
# df_vis['Continent'] = continent
# df_vis['averaged_antigenic_score'] = df_vis.groupby(['Country', 'Pango lineage'])['mutation_score'].transform('mean')
# # Dropping redundant columns and duplicates
# # df_vis = df_vis['Pango lineage', 'Country', ' averaged_mutation_score']
# df_vis.drop(["Location"], axis = 1, inplace = True)
# df_vis = df_vis.drop_duplicates()
#
# # Keeping variants that have an average mutation score greater than the assigned VOC threshold
# df_vis_threshold = df_vis[df_vis["averaged_antigenic_score"] > 1.025]
# df_vis_threshold['averaged_antigenic_score_per_country'] = df_vis_threshold.groupby('Country')['averaged_mutation_score'].transform('mean')
# df_vis_threshold_averaged = df_vis_threshold[['Continent','Country', 'averaged_antigenic_score_per_country']]
# df_vis_threshold_averaged = df_vis_threshold_averaged.drop_duplicates()
# df_vis_threshold_averaged
#
# # European only df for focused visualization:
# df_vis_threshold_averaged_eu = df_vis_threshold_averaged[df_vis_threshold_averaged['Continent'] == "Europe"]
# df_vis_threshold_averaged_eu
#
# # Creating global and european map of averaged antigenic scores
# labels_dict = dict(zip(df_vis_threshold_averaged.Country, df_vis_threshold_averaged.averaged_mutation_score_per_country))
# fig = px.choropleth(df_vis_threshold_averaged, locations = "Country",
#                            locationmode = "country names",
#                            color = "averaged_mutation_score_per_country",
#                            color_continuous_scale = 'spectral_r',
#                            labels = labels_dict,
#                            title = 'Mutation Scores per Country')
# fig.write_html(output + 'mutation_score_map.html')
# fig_eu = px.choropleth(df_vis_threshold_averaged_eu, locations = "Country",
#                            locationmode = "country names",
#                            color = "averaged_mutation_score_per_country",
#                            color_continuous_scale = 'spectral_r',
#                            labels = labels_dict,
#                            title = 'Mutation Scores per Country',
#                            scope = 'europe')
# fig_eu.write_html(output + 'mutation_score_map_europe.html')

