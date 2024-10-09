#!/usr/bin/env python
# coding: utf-8

# The purpose of this script is to test the analysis without using assigned weights but rather a weight of 1 at every true positive site

import os.path
import sys
import pandas as pd
import numpy as np
import re
from itertools import chain
from time import process_time
from datetime import date
import datetime
from dateutil.relativedelta import relativedelta

start_time = datetime.datetime.now()

print(sys.argv[0])
print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])
print(sys.argv[5])
print(sys.argv[6])
print(len(sys.argv), "\n")

print("Reading in the metadata file: ")
columns = ['Location', 'Collection date', 'Accession ID', 'Pango lineage', 'AA Substitutions', 'Host']
print("Chunking metadata file")
chunked_df = pd.read_csv(sys.argv[1], sep = '\t', usecols = columns, dtype = {"Host": "category", "Location": "category", "Pango lineage": "category"}, iterator = True, chunksize = 1000)
metadata = pd.concat(chunked_df, ignore_index = True)
    #metadata = pd.DataFrame()
    #for chunk in pd.read_csv(sys.argv[1], sep = '\t', usecols = columns, chunksize = 10):
    #    metadata = pd.concat([metadata, chunk], axis = 0)
#metadata = pd.read_csv(sys.argv[1], sep = '\t', usecols = columns)
tpSites = pd.read_csv(sys.argv[2], sep = ',')
output = sys.argv[3]
weights = pd.read_csv(sys.argv[4], sep = '\t')
voc_df = pd.read_csv(sys.argv[5], sep = ',')
try:
    seqsUI = pd.read_csv(sys.argv[6], sep = '\t')
except pd.errors.EmptyDataError:
    print("Sequences under review file was empty, skipping.")
    seqsUI = pd.DataFrame(columns = ['Accession ID','Collection date','Submission date','Location'])

metadata_time = datetime.datetime.now()
print('import data Duration: {}'.format(metadata_time - start_time))

def aa_substitution_filter(x, *tp_list):
    # Function to filter out the mutations occuring on the spike protein at known antigenic sites and return separate lists of the positions and mutations
    remove = ['180','181','182','183','184','185','186','187','188','189','118','218','318','418','518','618','718','818','918','1018','1118','1218']
    subslist = x.split(',')
    mutations_list = [i.split("_")[1]for i in subslist if '_' in i and 'Spike' in i] 
    mutations_list = [i for i in mutations_list if any(x in i for x in tp_list)]
    #mutations_list = [i for x in tp_list for i in mutations_list if x in i]
    mutations_list = set(mutations_list)-{i for x in remove for i in mutations_list if x in i} # removing incorrect matches, ie. that contain "18" but are not TP sites
    subs_positions = [re.findall(r'\d+', i)[0] for i in mutations_list] 
    originalAA = [re.findall('([a-zA-Z]*)\d*.*', i)[0] for i in mutations_list]
    changedAA = [re.findall('\d([a-zA-Z]*)', i)[-1] for i in mutations_list]
    return(subs_positions, originalAA, changedAA)

def seqsUI_filtration(seqsUI, max_month, max_year):
    # Function to return a list of sequences under review for removal prior to analysis
    print("Identifying Sequences Under Review")
    if len(seqsUI) == 0:
        print("Returning empty sequences under reivew list as df is empty.")
        seqsUI_list = []
    else:
        #seqsUI = seqsUI[seqsUI["Collection date"].str.contains("-")]
        #seqsUI['year'] = seqsUI["Collection date"].apply(lambda x: x.split("-")[0])
        #seqsUI['month'] = seqsUI["Collection date"].apply(lambda x: x.split("-")[1])
        #seqsUI_filtered = seqsUI[seqsUI.year == max_year]
        #seqsUI_filtered = seqsUI_filtered[seqsUI_filtered.month == max_month]
        seqsUI_filtered = seqsUI.loc[seqsUI['Collection date'].str.contains(max_year+'-'+max_month)]
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

def mutation_scores(input_df, tpSites_list):
    # Function to add mutation scores to the metadata file based on amino acid changes
    input_df = input_df.dropna(subset=['AA Substitutions'])

    # Converting AA Substitutions a list of aa positions for easier comprehension
    print("Converting AA Substitutions and Positions")
    tp_isolate_mutations = input_df['AA Substitutions'].apply(aa_substitution_filter, args=tpSites_list)
    input_df['Substitution_Positions'] = [i[0] for i in tp_isolate_mutations]
    input_df['Original_aa'] = [i[1] for i in tp_isolate_mutations]
    input_df['Changed_aa'] = [i[2] for i in tp_isolate_mutations]

    print("Expanding Substitution Positions")
    # metadata_expanded = aa_positions_unnest(input_df, ['Original_aa','Changed_aa'])
    #metadata_expanded = input_df.explode(['Substitution_Positions', 'Original_aa', 'Changed_aa'])
    metadata_expanded = input_df.set_index(['Accession ID','Collection date','Location','Host','Pango lineage','AA Substitutions']).apply(pd.Series.explode).reset_index()
    # Calculating the antigenic weight scores dependent on the amino acid changes (from a reference file)
    #print("Merging with the weights file")
    #metadata_filtered_weights = pd.merge(metadata_expanded, weights, how='left', on=['Original_aa', 'Changed_aa'])
    print("Calculating weights")
    metadata_filtered_weights = metadata_expanded
    metadata_filtered_weights["Weight"] = np.nan
    for value in metadata_filtered_weights['Substitution_Positions']:
        print(value)
        if value not in tpSites_list:
           metadata_filtered_weights["Weight"] = 0
        else:
           metadata_filtered_weights["Weight"] = 1
    print("mutation score weight complete")

    return (metadata_filtered_weights)

# Creating list of true positive antigenic sites which will be used to define antigenic weights
print("Creating list of known antigenic sites...")
tpSites_list = tpSites['tp_sites'].values.tolist()
tpSites_list = [str(item) for item in tpSites_list]
print("Done\n")

tplist_time = datetime.datetime.now()
print('tp site list Duration: {}'.format(tplist_time - metadata_time))

# Filtering by most recent month by chunking the Metadata file to reduce computational load
print("Calculating Mutation Scores....")
print("Getting most recent month")
if len(sys.argv) == 9:  
    print("Using user input month and year settings")
    # Using a predetermined month and year for the analysis that the user inputs
    max_month = str(sys.argv[7])
    max_year = str(sys.argv[8])
    print("MONTH: ", max_month)
    print("YEAR: ", max_year)

elif len(sys.argv) == 7: 
    # Calculating the most recent whole month and year based on today's date
    today = date.today()
    lastMonth = (today - relativedelta(months = 1)).month
    year = (today - relativedelta(months = 1)).year
    max_month = str(lastMonth)
    max_year = str(year)
    if len(max_month) < 2:
        max_month = "0" + str(max_month)
    print("Previous Month and Year:")
    print(max_year)
    print(max_month, "\n")

else:
    print("Wrong number of arguments, please use -h for more information")
    sys.exit(1)
    
max_year = str(max_year)
max_month = str(max_month)

month_calc_time = datetime.datetime.now()
print('month calculation Duration: {}'.format(month_calc_time - tplist_time))

print("Filtering Metadata by most recent month")
n = len(metadata) // 50  # 1500
t=0
print("n: ", n)

monthly_metadata = pd.DataFrame()
seqsUI_list = seqsUI_filtration(seqsUI, max_month, max_year)
df_chunks = np.array_split(metadata, n)
for chunk in df_chunks:
    t = t + 1
    dfMonth_chunk = chunk.loc[chunk['Collection date'].str.contains(max_year+'-'+max_month)]
    dfMonth_chunk = dfMonth_chunk[dfMonth_chunk['Accession ID'].isin(seqsUI_list) == False]
    #dfMonth_chunk = month_filtration(chunk, seqsUI_list, max_month, max_year)
    monthly_metadata = pd.concat([monthly_metadata, dfMonth_chunk], axis = 0)
print("Monthly Metadata 1st Attempt: ")
print("Length of monthly_metadata: ", len(monthly_metadata))
print(monthly_metadata.head())

# Fail safe to go back one month if metadata file empty for current month
if len(monthly_metadata) == 0:
    print("No isolate data for: ", max_month, "-", max_year)
    # Making file to remove months from month file, used later for frequency heatmap
    currentMonth = str(today.month)
    if len(currentMonth) < 2:
        currentMonth = "0" + str(currentMonth)
    currentYear = str(today.year)
    month_remove_file = open(output + "month_remove.txt", "w")
    month_remove_file.writelines(str(currentYear) + "-" + str(currentMonth))
    month_remove_file.writelines("\n" + str(max_year) + "-" + str(max_month))
    month_remove_file.close()

    twoMonths_ago = (today - relativedelta(months = 2)).month
    twoMonths_year = (today - relativedelta(months = 2)).year
    max_month = str(twoMonths_ago)
    max_year = str(twoMonths_year)
    if len(max_month) < 2:
        max_month = "0" + str(max_month)
    print("Now trying data from: ", max_month, "-", max_year)

    monthly_metadata = pd.DataFrame()
    seqsUI_list = seqsUI_filtration(seqsUI, max_month, max_year)
    for chunk in df_chunks:
        dfMonth_chunk = chunk.loc[chunk['Collection date'].str.contains(max_year + '-' + max_month)]
        #dfMonth_chunk = month_filtration(chunk, seqsUI_list, max_month, max_year)
        monthly_metadata = pd.concat([monthly_metadata, dfMonth_chunk], axis = 0)

# Going back two months if data still empty
if len(monthly_metadata) == 0:
    print("No isolate data for: ", max_month, "-", max_year)
    month_remove_file = open(output + "month_remove.txt", "a")
    month_remove_file.writelines("\n" + str(max_year) + "-" + str(max_month))
    month_remove_file.close()

    threeMonths_ago = (today - relativedelta(months = 3)).month
    threeMonths_year = (today - relativedelta(months = 3)).year
    max_month = str(threeMonths_ago)
    max_year = str(threeMonths_year)
    if len(max_month) < 2:
        max_month = "0" + str(max_month)
    print("Now trying data from: ", max_month, "-", max_year)

    monthly_metadata = pd.DataFrame()
    seqsUI_list = seqsUI_filtration(seqsUI, max_month, max_year)
    for chunk in df_chunks:
        dfMonth_chunk = chunk.loc[chunk['Collection date'].str.contains(max_year + '-' + max_month)]
        #dfMonth_chunk = month_filtration(chunk, seqsUI_list, max_month, max_year)
        monthly_metadata = pd.concat([monthly_metadata, dfMonth_chunk], axis = 0)

# Throwing error if data from two months ago is still empty
if len(monthly_metadata) == 0:
    print("ERROR - despite using data from two months ago, the monthly_metadata file is still empty, please check!")
    exit()
print("Monthly Filtration Complete\n")
month_filter_time = datetime.datetime.now()
print('monthly filtration Duration: {}'.format(month_filter_time - month_calc_time))

print("Writing Months Text File")
month_file = open(output + "month_vis.txt", "w")
month_file.writelines(str(max_month) + "-" + str(max_year))
month_file.close()

month_file_time = datetime.datetime.now()
print('month file Duration: {}'.format(month_file_time - month_filter_time))

# Changing dtype of columns for analysis
monthly_metadata['Pango lineage'] = monthly_metadata['Pango lineage'].astype(str)
monthly_metadata['Location'] = monthly_metadata['Location'].astype(str)
# Filtering by "Human" host data only
monthly_metadata = monthly_metadata.loc[monthly_metadata['Host'] == 'Human']

print("\nRunning Analysis - Calculating Antigenic Scores for Pango Lineages")
#metadata_filtered_weights = pd.DataFrame()
metadata_filtered_weights_list = []
print("Length of Monthly Metadata File to be chunked: ", len(monthly_metadata))
n = len(monthly_metadata) // 10 # 1500
print("Number of chunks for analysis: ", n, "\n")
t = 0
df_chunks = np.array_split(monthly_metadata, n)  # (metadata, n)
tstart = process_time()
for chunk in df_chunks:
    t = t + 1
    print(t)
    # Calculating Weights
    print("Calculating Weights")
    mutation_df = mutation_scores(chunk, tpSites_list)
    # Agreggating by accession id
    print("Aggregating by Accession ID")
    mutation_df = mutation_df.groupby(mutation_df['Accession ID']).aggregate({'Weight': 'sum'})
    metadata_filtered_weights_list.append(mutation_df)
    #metadata_filtered_weights = pd.concat([mutation_df, metadata_filtered_weights])

    print("Output DataFrame Shape: ")
    #print(metadata_filtered_weights.shape)
    print(mutation_df.shape)
    print("Chunk complete")
metadata_filtered_weights = pd.concat(metadata_filtered_weights_list)
tstop = process_time()
print("\nProcess Time: ", tstop - tstart)

analysis_time = datetime.datetime.now()
print('analysis time Duration: {}'.format(analysis_time - month_file_time))

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

df_final_time = datetime.datetime.now()
print('df_final Duration: {}'.format(df_final_time - analysis_time))

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
