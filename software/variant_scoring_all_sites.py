#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import numpy as np
import re
import plotly.express as px
from itertools import chain
from os.path import exists
from time import process_time

columns = ['Location', 'Collection date', 'Accession ID', 'Pango lineage', 'AA Substitutions']
metadata = pd.read_csv(sys.argv[0], sep = '\t', usecols = columns)
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

def combining_accessionID(input_df):
# Function to sum antigenic weights per assession ID, requires a column called 'Weight'
    metadata_filtered_combinedWeights = input_df.groupby(input_df['Accession ID']).aggregate({'Weight': 'sum'})
    x = pd.DataFrame.copy(input_df)
    x.drop(['Substitution_Positions', 'Original_aa', 'Changed_aa', 'Weight'], axis=1, inplace=True)
    x.drop_duplicates(inplace=True)
    df = pd.merge(x, metadata_filtered_combinedWeights, how='left', left_on='Accession ID', right_index=True)
    df.reset_index(drop=True, inplace=True)
    df['Weight'].fillna(0, inplace=True)
    return df

def mutation_scores(input_df): # , output_df
# Function to add mutation scores to the metadata file based on amino acid changes
    print("Dropping NaN rows ...")
    # input_df = input_df.dropna(subset = ['AA Substitutions'])
    input_df['AA Substitutions'] = input_df['AA Substitutions'].astype('string')

    # print("Converting AA Substitutions to a list of aa positions...")
    # print("Substitution_Positions")
    # substitutions = input_df['AA Substitutions'].apply(conversion_position)
    # input_df["Substitution_Positions"] = substitutions
    # print("Original Amino Acids")
    # original_AAs = input_df['AA Substitutions'].apply(conversion_originalAA)
    # input_df["Original_aa"] = original_AAs
    # print("Changed Amino Acids")
    # changed_AAs = input_df['AA Substitutions'].apply(conversion_changedAA)
    # input_df["Changed_aa"] = changed_AAs

    # Expanding list of original and changed amino acids 
    try:
        # substitutions = input_df['AA Substitutions'].apply(conversion_position)
        # input_df["Substitution_Positions"] = substitutions
        # original_AAs = input_df['AA Substitutions'].apply(conversion_originalAA)
        # input_df["Original_aa"] = original_AAs
        # changed_AAs = input_df['AA Substitutions'].apply(conversion_changedAA)
        # input_df["Changed_aa"] = changed_AAs
        substitutions = []
        original_AAs = []
        changed_AAs = []
        for i in input_df['AA Substitutions']:
            if len(i) == 0:
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
        if len(substitutions) == len(original_AAs) & len(substitutions) == len(changed_AAs):
            input_df["Substitution_Positions"] = substitutions
            input_df["Original_aa"] = original_AAs
            input_df["Changed_aa"] = changed_AAs
        else:
            print(
                "Length of spike substitution positions does not match the length of amino acid changes at those positions!")
            sys.exit()
        # Calculating the weight of spike amino acid changes
        cols = ['Substitution_Positions','Original_aa','Changed_aa']
        metadata_expanded = input_df.explode(cols)
        metadata_weights_expanded = pd.merge(metadata_expanded, weights, how = 'left', on = ['Original_aa', 'Changed_aa'])
        metadata_weights = combining_accessionID(metadata_weights_expanded)
    except:
#         input_df.drop()
        metadata_weights = input_df
        # metadata_weights['Substitution_Positions'] = np.NAN
        metadata_weights['Weight'] = 0
    # print(pd.DataFrame.head(metadata_weights))
    return(metadata_weights)

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

    # # Summing Weights by accession ID
    # metadata_filtered_combinedWeights = metadata_filtered_weights.groupby(metadata_filtered_weights['Accession ID']).aggregate({'Weight': 'sum'})
    # x = pd.DataFrame.copy(metadata_filtered_weights)
    # x.drop(['Substitution_Positions', 'Original_aa', 'Changed_aa', 'Weight'], axis=1, inplace=True)
    # x.drop_duplicates(inplace=True)
    # df = pd.merge(x, metadata_filtered_combinedWeights, how='left', left_on='Accession ID', right_index=True)
    # df.reset_index(drop=True, inplace=True)
    # df['Weight'].fillna(0, inplace=True)

    return output_df

### Applying analyses to metadata file

# Splitting Metadata file into 4 parts:
pt1 = int(len(metadata) / 4)
pt2 = int(len(metadata) / 2)
pt3 = pt2 + pt1
metadata1 = metadata.iloc[:pt1]
metadata2 = metadata.iloc[(pt1 + 1):pt2]
metadata3 = metadata.iloc[(pt2 + 1):pt3]
metadata4 = metadata.iloc[(pt3 + 1):(len(metadata))]

# Calculating mutation scores
if exists(output + "antigenic_scores_pt1_without_TPsites_all.csv") == True:
    output_df1 = pd.read_csv(output + "antigenic_scores_pt1_without_TPsites_all.csv", sep = '\t')
else:
    print("Part 1")
    output_df1 = analysis(metadata1)
    output_df1.to_csv(output + "antigenic_scores_pt1_without_TPsites_all.csv", sep = '\t', index = False, header = True)
if exists(output + "antigenic_scores_pt2_without_TPsites_all.csv") == True:
    output_df2 = pd.read_csv(output + "antigenic_scores_pt2_without_TPsites_all.csv", sep='\t')
else:
    print("Part 2")
    output_df2 = analysis(metadata2)
    output_df2.to_csv(output + "antigenic_scores_pt2_without_TPsites_all.csv", sep = '\t', index = False, header = True)
if exists(output + "antigenic_scores_pt3_without_TPsites_all.csv") == True:
    output_df3 = pd.read_csv(output + "antigenic_scores_pt3_without_TPsites_all.csv", sep='\t')
else:
    print("Part 3")
    output_df3 = analysis(metadata3)
    output_df3.to_csv(output + "antigenic_scores_pt3_without_TPsites_all.csv", sep = '\t', index = False, header = True)
if exists(output + "antigenic_scores_pt4_without_TPsites_all.csv") == True:
    output_df4 = pd.read_csv(output + "antigenic_scores_pt4_without_TPsites_all.csv", sep='\t')
else:
    print("Part 4")
    output_df4 = analysis(metadata4)
    output_df4.to_csv(output + "antigenic_scores_pt4_without_TPsites_all.csv", sep = '\t', index = False, header = True)

# Combining to one final dataframe
print("FINAL METADATA FILTERED DF: ")
metadata_filtered = pd.concat([output_df1, output_df2, output_df3, output_df4])
# Summing weights another time by accession ID, as the output dataframes have been combined:
# metadata_filtered_combinedWeights = metadata_filtered.groupby(metadata_filtered['Accession ID']).aggregate({'Weight':'sum'})
# x = pd.DataFrame.copy(metadata_filtered)
# x.drop(['Weight'], axis = 1, inplace = True)
# x.drop_duplicates(inplace = True)
# df = pd.merge(x, metadata_filtered_combinedWeights, how = 'left', left_on = 'Accession ID', right_index = True)
# df.reset_index(drop=True, inplace = True)
# df['Weight'].fillna(0, inplace = True)
# Now adding an antigenic score
metadata_filtered['antigenic_score'] = metadata_filtered.groupby('Pango lineage')['Weight'].transform('mean')
metadata_filtered.to_csv(output + "antigenic_scores_without_TPsites_all.csv", sep = '\t', index = False, header = True)

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
df_merged_ranked.to_csv(output + "antigenic_scores_ranked_without_TPsites_with_WHO.csv", sep = '\t', index = False, header = True)
print("Ranked Antigenic Scores COMPLETE")
