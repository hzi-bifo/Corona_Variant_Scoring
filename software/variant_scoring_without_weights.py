#!/usr/bin/env python
# coding: utf-8

import sys
from os.path import exists
import pandas as pd
import numpy as np
from time import process_time

columns = ['Location', 'Collection date', 'Accession ID', 'Pango lineage', 'AA Substitutions']
metadata = pd.read_csv(sys.argv[0], sep = '\t', usecols = columns)
output = sys.arg[1]
voc_df = pd.read_csv(sys.arg[2], sep = '\t')
tpSites = pd.read_csv(sys.arg[3], sep = ',')

### Defining analysis functions

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
            # x = pos[1:]
            pos = "".join([i for i in pos if i.isdigit()])
            mutationsSubsList.append(pos)
        else:
            mutationsSubsList.append(None)
    return mutationsSubsList

def comparing_lists(list1):
# Function to find the intersection of two lists and returns the length of the intersection
#     common = list(set(list1).intersection(set(tpSites_list)))
    common = len(list(set(list1) & set(tpSites_list)))
    if pd.isna(common) == True:
        common = 0
    return common

def mutation_scores(input_df): # , output_df
# Function to add mutation scores to based on the number of known antigenic sites with aa changes
#     print("Dropping NaN rows ...")
    # input_df = input_df.dropna(subset = ['AA Substitutions'])
    # Converting AA Substitutions a list of aa positions for easier comprehension
    input_df['AA Substitutions'] = input_df['AA Substitutions'].astype('string')
    try:
        print("\nSubstitutions")
        # substitutions = input_df['AA Substitutions'].apply(conversion_position)
        substitutions = []
        antigenic_scores = []
        for i in input_df['AA Substitutions']:
            if len(i) == 0:
                score = 0
                substitution_pos = np.NaN
                substitutions.append(substitution_pos)
                antigenic_scores.append(score)
            else:
                substitution_pos = conversion_position(i)
                substitutions.append(substitution_pos)
                score = comparing_lists(substitution_pos)
                if pd.isna(score):
                    score = 0
                antigenic_scores.append(score)
        input_df["Substitution_Positions"] = substitutions
        input_df["Weight"] = antigenic_scores
        # antigenic_scores = input_df['Substitution_Positions'].apply(comparing_lists)
        # input_df['Weight'] = antigenic_scores
        # print("Scoring")
    except:
        print("Unable to run normally - filling with NA")
        print("Substitutions")
        input_df["Substitution_Positions"] = np.NaN
        print("Scoring")
        input_df['Weight'] = 0
    return input_df

def analysis(input_df):
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    df3 = pd.DataFrame()
    df4 = pd.DataFrame()
    df5 = pd.DataFrame()
    df6 = pd.DataFrame()
    df7 = pd.DataFrame()
    df8 = pd.DataFrame()
    len_df = input_df.shape[0]
    print("Number of Sequences: ", len_df)
    n = np.ceil(len_df / 100)
    print("Number of Chunks: ", n)
    l1 = int(n / 8)
    l2 = l1 + l1
    l3 = l2 + l1
    l4 = int(n / 2)
    l5 = l4 + l1
    l6 = l5 + l1
    l7 = l6 + l1
    print("Dataframe partitions: ", l1, l2, l3, l4, l5, l6, l7, n)
    t = 0
    df_chunks = np.array_split(input_df, n)
    tstart = process_time()
    for chunk in df_chunks:
        print("Chunk Shape: ", chunk.shape[0])
        t = t+1
        print("Iteration: ", t)
        if (t < l1):
            print("DF1: 0 - ", l1)
            mutation_df = mutation_scores(chunk)
            if chunk.shape[0] == mutation_df.shape[0]:
                df1 = pd.concat([mutation_df, df1])
            else:
                df_diff = chunk.merge(mutation_df, indicator=True, how='left').loc[lambda x: x['_merge'] != 'both']
                print(df_diff)
                print("Size of Chunk: ", chunk.shape)
                print(pd.DataFrame.tail(chunk))
                print("Size of Output DF: ", mutation_df.shape)
                print(pd.DataFrame.tail(mutation_df))
                print("THE ACCESSION IDS DO NOT MATCH BETWEEN THE CHUNK AND OUTPUT DF")
                sys.exit()
            print("Output DataFrame Shape: ")
            print(df1.shape)
            print("\n")
        if (l2 > t >= l1):
            print("DF2: ", l1, " - ", l2)
            mutation_df = mutation_scores(chunk)
            if chunk.shape[0] == mutation_df.shape[0]:
                df2 = pd.concat([mutation_df, df2])
            else:
                df_diff = chunk.merge(mutation_df, indicator=True, how='left').loc[lambda x: x['_merge'] != 'both']
                print(df_diff)
                print("Size of Chunk: ", chunk.shape)
                print(pd.DataFrame.tail(chunk))
                print("Size of Output DF: ", mutation_df.shape)
                print(pd.DataFrame.tail(mutation_df))
                print("THE ACCESSION IDS DO NOT MATCH BETWEEN THE CHUNK AND OUTPUT DF")
                sys.exit()
            print("Output DataFrame Shape: ")
            print(df2.shape)
            print("\n")
        if (l3 > t >= l2):
            print("DF3: ", l2, " - ", l3)
            mutation_df = mutation_scores(chunk)
            if chunk.shape[0] == mutation_df.shape[0]:
                df3 = pd.concat([mutation_df, df3])
            else:
                df_diff = chunk.merge(mutation_df, indicator=True, how='left').loc[lambda x: x['_merge'] != 'both']
                print(df_diff)
                print("Size of Chunk: ", chunk.shape)
                print(pd.DataFrame.tail(chunk))
                print("Size of Output DF: ", mutation_df.shape)
                print(pd.DataFrame.tail(mutation_df))
                print("THE ACCESSION IDS DO NOT MATCH BETWEEN THE CHUNK AND OUTPUT DF")
                sys.exit()
            print("Output DataFrame Shape: ")
            print(df3.shape)
            print("\n")
        if (l4 > t >= l3):
            print("DF4: ", l3, " - ", l4)
            mutation_df = mutation_scores(chunk)
            if chunk.shape[0] == mutation_df.shape[0]:
                df4 = pd.concat([mutation_df, df4])
            else:
                df_diff = chunk.merge(mutation_df, indicator=True, how='left').loc[lambda x: x['_merge'] != 'both']
                print(df_diff)
                print("Size of Chunk: ", chunk.shape)
                print(pd.DataFrame.tail(chunk))
                print("Size of Output DF: ", mutation_df.shape)
                print(pd.DataFrame.tail(mutation_df))
                print("THE ACCESSION IDS DO NOT MATCH BETWEEN THE CHUNK AND OUTPUT DF")
                sys.exit()
            print("Output DataFrame Shape: ")
            print(df4.shape)
            print("\n")
        if (l5 > t >= l4):
            print("DF5: ", l4, " - ", l5)
            mutation_df = mutation_scores(chunk)
            if chunk.shape[0] == mutation_df.shape[0]:
                df5 = pd.concat([mutation_df, df5])
            else:
                df_diff = chunk.merge(mutation_df, indicator=True, how='left').loc[lambda x: x['_merge'] != 'both']
                print(df_diff)
                print("Size of Chunk: ", chunk.shape)
                print(pd.DataFrame.tail(chunk))
                print("Size of Output DF: ", mutation_df.shape)
                print(pd.DataFrame.tail(mutation_df))
                print("THE ACCESSION IDS DO NOT MATCH BETWEEN THE CHUNK AND OUTPUT DF")
                sys.exit()
            print("Output DataFrame Shape: ")
            print(df5.shape)
            print("\n")
        if (l6 > t >= l5):
            print("DF6: ", l5, " - ", l6)
            mutation_df = mutation_scores(chunk)
            if chunk.shape[0] == mutation_df.shape[0]:
                df6 = pd.concat([mutation_df, df6])
            else:
                df_diff = chunk.merge(mutation_df, indicator=True, how='left').loc[lambda x: x['_merge'] != 'both']
                print(df_diff)
                print("Size of Chunk: ", chunk.shape)
                print(pd.DataFrame.tail(chunk))
                print("Size of Output DF: ", mutation_df.shape)
                print(pd.DataFrame.tail(mutation_df))
                print("THE ACCESSION IDS DO NOT MATCH BETWEEN THE CHUNK AND OUTPUT DF")
                sys.exit()
            print("Output DataFrame Shape: ")
            print(df6.shape)
            print("\n")
        if (l7 > t >= l6):
            print("DF7: ", l6, " - ", l7)
            mutation_df = mutation_scores(chunk)
            if chunk.shape[0] == mutation_df.shape[0]:
                df7 = pd.concat([mutation_df, df7])
            else:
                df_diff = chunk.merge(mutation_df, indicator=True, how='left').loc[lambda x: x['_merge'] != 'both']
                print(df_diff)
                print("Size of Chunk: ", chunk.shape)
                print(pd.DataFrame.tail(chunk))
                print("Size of Output DF: ", mutation_df.shape)
                print(pd.DataFrame.tail(mutation_df))
                print("THE ACCESSION IDS DO NOT MATCH BETWEEN THE CHUNK AND OUTPUT DF")
                sys.exit()
            print("Output DataFrame Shape: ")
            print(df7.shape)
            print("\n")
        if (n >= t >= l7):
            print("DF8: ", l7, " - ", n)
            mutation_df = mutation_scores(chunk)
            if chunk.shape[0] == mutation_df.shape[0]:
                df8 = pd.concat([mutation_df, df8])
            else:
                df_diff = chunk.merge(mutation_df, indicator=True, how='left').loc[lambda x: x['_merge'] != 'both']
                print(df_diff)
                print("Size of Chunk: ", chunk.shape)
                print(pd.DataFrame.tail(chunk))
                print("Size of Output DF: ", mutation_df.shape)
                print(pd.DataFrame.tail(mutation_df))
                print("THE ACCESSION IDS DO NOT MATCH BETWEEN THE CHUNK AND OUTPUT DF")
                sys.exit()
            print("Output DataFrame Shape: ")
            print(df8.shape)
            print("\n")
    tstop = process_time()
    print("Process Time: ", tstop - tstart)
    output_df = pd.concat([df1, df2, df3, df4, df5, df6, df7, df8])
    return output_df

### Applying analyses to metadata file

# Creating list of known antigenic sites
tpSites_list = tpSites['tp_sites'].values.tolist()
tpSites_list = [str(item) for item in tpSites_list]

# Splitting Metadata file into 4 parts:
pt1 = int(len(metadata) / 4)
pt2 = int(len(metadata) / 2)
pt3 = pt2 + pt1
metadata1 = metadata.iloc[:pt1]
metadata2 = metadata.iloc[(pt1 + 1):pt2]
metadata3 = metadata.iloc[(pt2 + 1):pt3]
metadata4 = metadata.iloc[(pt3 + 1):(len(metadata))]

# Calculating mutation scores
if exists(output + "antigenic_scores_pt1_without_weights_all.csv") == True:
    output_df1 = pd.read_csv(output + "antigenic_scores_pt1_without_weights_all.csv", sep = '\t')
else:
    print("Part 1")
    output_df1 = analysis(metadata1)
    output_df1.to_csv(output + "antigenic_scores_pt1_without_weights_all.csv", sep = '\t', index = False, header = True)
if exists(output + "antigenic_scores_pt2_without_weights_all.csv") == True:
    output_df2 = pd.read_csv(output + "antigenic_scores_pt2_without_weights_all.csv", sep='\t')
else:
    print("Part 2")
    output_df2 = analysis(metadata2)
    output_df2.to_csv(output + "antigenic_scores_pt2_without_weights_all.csv", sep = '\t', index = False, header = True)
if exists(output + "antigenic_scores_pt3_without_weights_all.csv") == True:
    output_df3 = pd.read_csv(output + "antigenic_scores_pt3_without_weights_all.csv", sep='\t')
else:
    print("Part 3")
    output_df3 = analysis(metadata3)
    output_df3.to_csv(output + "antigenic_scores_pt3_without_weights_all.csv", sep = '\t', index = False, header = True)
if exists(output + "antigenic_scores_pt4_without_weights_all.csv") == True:
    output_df4 = pd.read_csv(output + "antigenic_scores_pt4_without_weights_all.csv", sep='\t')
else:
    print("Part 4")
    output_df4 = analysis(metadata4)
    output_df4.to_csv(output + "antigenic_scores_pt4_without_weights_all.csv", sep = '\t', index = False, header = True)

# Combining to one final dataframe
print("FINAL METADATA FILTERED DF: ")
metadata_filtered = pd.concat([output_df1, output_df2, output_df3, output_df4])

# if metadata_filtered.shape[0] == metadata.shape[0]:
#     metadata_filtered['antigenic_score'] = metadata_filtered.groupby('Pango lineage')['Weight'].transform('mean')
#     metadata_filtered.to_csv(output + "antigenic_scores_without_weights_all.csv", sep = '\t', index = False, header = True)
# else:
#     print("Metadata: ", metadata.shape[0])
#     print("DF1: ", output_df1.shape[0])
#     print("DF2: ", output_df2.shape[0])
#     print("DF3: ", output_df3.shape[0])
#     print("DF4: ", output_df4.shape[0])
#     print("Metadata Filtered: ", metadata_filtered.shape[0])
#     raise Exception("The length of the two dataframes do not match, please check!")

metadata_filtered['antigenic_score'] = metadata_filtered.groupby('Pango lineage')['Weight'].transform('mean')
metadata_filtered.to_csv(output + "antigenic_scores_without_weights_all.csv", sep = '\t', index = False, header = True)

# Ranking pango lineages based on average mutation score across all the sequences
df_ranked = metadata_filtered[["Pango lineage", "antigenic_score"]]
df_ranked = df_ranked.drop_duplicates()
df_ranked.dropna(axis = 0, inplace = True)
df_ranked.sort_values(by = ['antigenic_score'], ascending = False, inplace = True)
df_ranked.reset_index(inplace = True)
df_ranked["rank"] = df_ranked['antigenic_score'].rank(ascending = False)
df_ranked = df_ranked[["Pango lineage", "antigenic_score", "rank"]]
df_merged_ranked = df_ranked.merge(voc_df, how = "left", on = "Pango lineage")
df_merged_ranked['WHO_label'] = df_merged_ranked['WHO_label'].fillna("Non Variant of Concern")
df_merged_ranked.to_csv(output + "antigenic_scores_ranked_without_weights_with_WHO.csv", sep = '\t', index = False, header = True)
print("Ranked Antigenic Scores COMPLETE")
