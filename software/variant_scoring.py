#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import numpy as np
import re
from itertools import chain
from time import process_time
from datetime import date
#from datetime import datetime
import datetime
start_time = datetime.datetime.now()

print(sys.argv[0])
print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])
print(sys.argv[5])
print(sys.argv[6])
print(len(sys.argv), "\n")

columns = ['Location', 'Collection date', 'Accession ID', 'Pango lineage', 'AA Substitutions']
metadata = pd.read_csv(sys.argv[1], sep = '\t', usecols = columns)
tpSites = pd.read_csv(sys.argv[2], sep = ',')
output = sys.argv[3]
weights = pd.read_csv(sys.argv[4], sep = '\t')
voc_df = pd.read_csv(sys.argv[5], sep = ',')
seqsUI = pd.read_csv(sys.argv[6], sep = '\t')

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
    return (mutationsSubsList)


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
    return (original_aaList)


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
    return (changed_aaList)

def seqsUI_filtration(seqsUI, max_month, max_year):
    # Function to return a list of sequences under review for removal prior to analysis
    print("Identifying Sequences Under Review")
    seqsUI = seqsUI[seqsUI["Collection date"].str.contains("-")]
    seqsUI['year'] = seqsUI["Collection date"].apply(lambda x: x.split("-")[0])
    seqsUI['month'] = seqsUI["Collection date"].apply(lambda x: x.split("-")[1])
    seqsUI_filtered = seqsUI[seqsUI.year == max_year]
    seqsUI_filtered = seqsUI_filtered[seqsUI_filtered.month == max_month]
    seqsUI_list = seqsUI_filtered["Accession ID"].to_list()
    return (seqsUI_list)

def month_filtration(metadata, seqsUI_list, max_month, max_year):
    # Function to filter metadata file by most recent month by chunking dataframe
    metadata = metadata[metadata["Collection date"].str.contains("-")]
    metadata['year'] = metadata["Collection date"].apply(lambda x: x.split("-")[0])
    metadata['month'] = metadata["Collection date"].apply(lambda x: x.split("-")[1])
    dftemp = metadata[metadata.year == max_year]
    dftemp = dftemp[dftemp.month == max_month]
    dftemp = dftemp[dftemp['Accession ID'].isin(seqsUI_list) == False]
    return (dftemp)

def mutation_scores(input_df, tpSites_list):  # , output_df
    # Function to add mutation scores to the metadata file based on amino acid changes
    input_df = input_df.dropna(subset=['AA Substitutions'])

    # Converting AA Substitutions a list of aa positions for easier comprehension
    try:
        print("Converting AA Substitutions and Positions")
        input_df['AA Substitutions'] = input_df['AA Substitutions'].astype('string')
        substitutions = input_df['AA Substitutions'].apply(conversion_position)
        input_df["Substitution_Positions"] = substitutions
        original_AAs = input_df['AA Substitutions'].apply(conversion_originalAA)
        input_df["Original_aa"] = original_AAs
        changed_AAs = input_df['AA Substitutions'].apply(conversion_changedAA)
        input_df["Changed_aa"] = changed_AAs
        if len(substitutions) != len(original_AAs):
            print("Conversion of amino acid list was incorrect and the lengths of the substitution positions defers from the amino acids!")
            print("Substitution Positions: ", len(substitutions))
            print("Number of amino acids: ", len(original_AAs))
            exit()
    except:
        print("Unable to parse the AA Substitutions column of the metadata file, please investigate.")
        exit()

    try:
        print("Expanding Substitution Positions")
        metadata_expanded = input_df.explode(['Substitution_Positions','Original_aa','Changed_aa'])
        print(pd.DateFrame.head(metadata_expanded))
        # Calculating the antigenic weight scores dependent on the amino acid changes (from a reference file)
        print("Filtering based on tp list")
        metadata_filtered = metadata_expanded[metadata_expanded['Substitution_Positions'].isin(tpSites_list)]
        print("Merging with the weights file")
        metadata_filtered_weights = pd.merge(metadata_filtered, weights, how='left', on=['Original_aa', 'Changed_aa'])
        print("mutation score weight complete")
    except:
        print("Unable to expand Substitution_Positions, chunking instead...")
        metadata_filtered = pd.DataFrame()
        n = 100
        t = 0
        df_chunks = np.array_split(input_df, n) 
        for chunk in df_chunks:
            t = t + 1
            metadata_expanded = chunk.explode(['Substitution_Positions','Original_aa','Changed_aa'])
            metadata_filtered_chunk = metadata_expanded[metadata_expanded['Substitution_Positions'].isin(tpSites_list)]
            metadata_filtered = pd.concat([metadata_filtered, metadata_filtered_chunk])
        metadata_filtered_weights = pd.merge(metadata_filtered, weights, how='left', on=['Original_aa', 'Changed_aa'])
        print("mutation score weight complete")

    return (metadata_filtered_weights)

# Creating list of true positive antigenic sites which will be used to define antigenic weights
print("Creating list of known antigenic sites...")
tpSites_list = tpSites['tp_sites'].values.tolist()
tpSites_list = [str(item) for item in tpSites_list]
print("Done\n")

# Calculating Antigenic Scores by chunking Metadata
print("Calculating Mutation Scores....")
print("Getting most recent month")
if len(sys.argv) == 9:  
    print("Using user input month and year settings")
    # Using a predetermined month and year for the analysis that the user inputs
    max_month = str(sys.argv[6])  # add one back
    max_year = str(sys.argv[7])  # add one back

elif len(sys.argv) == 7: 
    # Calculating the most recent whole month and year based on today's date
    today = date.today()
    firstOfMonth = today.replace(day=1)
    lastMonth = firstOfMonth - datetime.timedelta(days=1)
    max_month = str(lastMonth.month)
    # Calculating the year, needs to change if the previous month was December
    if max_month == "12":
        lastYear = firstOfMonth - datetime.timedelta(days=366)
        max_year = str(lastYear) 
    else:
        max_year = str(today.year)
    if len(max_month) < 2:
        max_month = "0" + str(max_month)

    print("Previous Month and Year:")
    print(max_year)
    print(max_month, "\n")

else:
    print("Wrong number of arguments, please use -h for more information")
    sys.exit(1)
    
max_year = str(max_year) #New line
max_month = str(max_month) # New line

print("Filtering Metadata by most recent month")
n = len(metadata) // 50  # 1500
print("n: ", n)

monthly_metadata = pd.DataFrame()
seqsUI_list = seqsUI_filtration(seqsUI, max_month, max_year)
df_chunks = np.array_split(metadata, n)
for chunk in df_chunks:
    dfMonth_chunk = month_filtration(chunk, seqsUI_list, max_month, max_year)
    monthly_metadata = pd.concat([monthly_metadata, dfMonth_chunk], axis = 0)
print("Monthly Metadata 1st Attempt: ")
print("Length of monthly_metadata: ", len(monthly_metadata))
print("Column names: ", monthly_metadata.columns.values)
print(monthly_metadata.head())

if len(monthly_metadata) == 0:
    print("No isolate data for: ", max_month, "-", max_year)
    max_month = int(max_month) - 1
    if max_month == 0:
        max_month = "12"
    else:
        max_month = str(max_month)
        if len(max_month) < 2:
            max_month = "0" + str(max_month)
    print("Now using data from: ", max_month, "-", max_year)

    monthly_metadata = pd.DataFrame()
    seqsUI_list = seqsUI_filtration(seqsUI, max_month, max_year)
    for chunk in df_chunks:
        dfMonth_chunk = month_filtration(chunk, seqsUI_list, max_month, max_year)
        monthly_metadata = pd.concat([monthly_metadata, dfMonth_chunk], axis = 0)
if len(monthly_metadata) == 0:
    print("ERROR - despite using data from two months ago, the monthly_metadata file is still empty, please check!")
    exit()
print("Monthly Filtration Complete\n")

print("Writing Months Text File")
month_file = open(output + "month_vis.txt", "w")
month_file.writelines(str(max_month) + "-" + str(max_year))
month_file.close()

print("\nRunning Analysis - Calculating Antigenic Scores for Pango Lineages")
metadata_filtered_weights = pd.DataFrame()
print("Length of Monthly Metadata File to be chunked: ", len(monthly_metadata))
n = len(monthly_metadata) // 50 # 1500
print("Number of chunks for analysis: ", n, "\n")
t = 0
df_chunks = np.array_split(monthly_metadata, n)  # (metadata, n)
tstart = process_time()
for chunk in df_chunks:
    print(chunk.shape[0])
    t = t + 1
    print(t)

    if len(chunk) == 0:
        print("Empty df after filtering")
        continue
    else:
        # Calculating Weights
        mutation_df = mutation_scores(chunk, tpSites_list)
        # Agreggating by accession id
        mutation_df = mutation_df.groupby(mutation_df['Accession ID']).aggregate({'Weight': 'sum'})
        metadata_filtered_weights = pd.concat([mutation_df, metadata_filtered_weights])

    print("Output DataFrame Shape: ")
    print(metadata_filtered_weights.shape)
    print("Chunk complete")
tstop = process_time()
print("\nProcess Time: ", tstop - tstart)

# Summing Weights by accession ID
print("Merging Dataframes")
df = pd.merge(monthly_metadata, metadata_filtered_weights, how='left', left_on='Accession ID', right_index=True)
df['Weight'].fillna(0, inplace=True)

print("Creating final dataframe...")
# Calculating final antigenic score by averaging the scores across the lineages
df['antigenic_score'] = df.groupby('Pango lineage')['Weight'].transform('mean')
df_final = df[["Accession ID", "Collection date", "Location", "Pango lineage", "AA Substitutions", "Weight", "antigenic_score"]]

print("Saving dataframe...")
df_final.to_csv(output + "antigenic_scores_all.csv", sep='\t', index=False, header=True)

# Ranking pango lineages based on average mutation score across all the sequences
df_ranked = df_final[["Pango lineage", "antigenic_score"]]
df_ranked = df_ranked.drop_duplicates()
df_ranked.dropna(axis=0, inplace=True)
df_ranked.sort_values(by=['antigenic_score'], ascending=False, inplace=True)
df_ranked.reset_index(inplace=True)
df_ranked["rank"] = df_ranked['antigenic_score'].rank(ascending=False)
df_ranked = df_ranked[["Pango lineage", "antigenic_score", "rank"]]
df_merged_ranked = df_ranked.merge(voc_df, how="left", on="Pango lineage")
df_merged_ranked['WHO_label'] = df_merged_ranked['WHO_label'].fillna("Non Variant of Concern")
df_merged_ranked.to_csv(output + "antigenic_scores_ranked_with_WHO.csv", sep='\t', index=False, header=True)
print("Ranked Antigenic Scores COMPLETE")

### Creation of Dataframes for Visualizations

# Creating Dataframe for global map visualization
df_vis = df_final[['Location', 'Pango lineage', 'antigenic_score', 'Collection date']]
df_vis["Continent"] = df_vis["Location"].apply(lambda x: x.split("/")[0])
df_vis["Country"] = df_vis["Location"].apply(lambda x: x.split("/")[1])
country = []
continent = []
for item in df_vis["Country"]:
    x = item.lstrip()
    y = x.rstrip()
    country.append(y)
for item in df_vis["Continent"]:
    x = item.lstrip()
    y = x.rstrip()
    continent.append(y)
df_vis['Country'] = country
df_vis['Continent'] = continent
df_vis.drop(["Location"], axis=1, inplace=True)

# Calculating the frequency for each lineage per country:
df_vis = df_vis[df_vis['Pango lineage'] != 'None']
n_lineages = pd.Series(df_vis.groupby(['Country', 'Pango lineage'])['Pango lineage'].count(), name="n_lineages")
n_lineages_df = n_lineages.to_frame().reset_index()
n_lineages_df.ffill(axis=0, inplace=True)
df_vis = df_vis.merge(n_lineages_df, how='left', on=['Country', 'Pango lineage'])
total_lineages = pd.Series(df_vis.groupby('Country')['Country'].count(), name='total_lineages')
df_vis = df_vis.merge(total_lineages.to_frame(), how='left', on='Country')
df_vis['frequency'] = df_vis['n_lineages'].div(df_vis['total_lineages'])

# Calculating average antigenic score for each lineage in the selected month:
df_vis.drop(['Collection date', 'n_lineages', 'total_lineages'], axis=1, inplace=True)  # 'collection_date_list'
df_vis_threshold_averaged = df_vis.drop_duplicates()

# Calculating antigenic score per country:
df_vis_threshold_averaged['score'] = df_vis_threshold_averaged['antigenic_score'] * df_vis_threshold_averaged[
    'frequency']
df_vis_threshold_averaged['country_score'] = df_vis_threshold_averaged.groupby('Country')['score'].transform('sum')
df = df_vis_threshold_averaged.drop(['Pango lineage', 'frequency', 'antigenic_score', 'score'], axis=1, )
df.drop_duplicates(inplace=True)

# Saving visualization dataframe:
df.to_csv(output + "antigenic_scores_map_visualization.csv", sep='\t', index=False, header=True)

end_time = datetime.datetime.now()
print('Duration: {}'.format(end_time - start_time))