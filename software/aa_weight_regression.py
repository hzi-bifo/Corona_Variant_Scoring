#!/usr/bin/env python
# coding: utf-8

# Purpose of script is to run two regression predicting first the impact of each amino acid on antigenicity and then
# both the amino acid properties on antigenicity.

# To Do List:
# - import necessary functions from previous analysis (for cleaner code)

import sys
import glob
import os.path
from os.path import exists
import re
import pandas as pd
import numpy as np
import plotly as plt
import sklearn.linear_model

indir = sys.argv[1]
who_indir = sys.argv[2]
outdir = sys.argv[3]
ref_dir = sys.argv[4]

columns = ['Collection date', 'Accession ID', 'Pango lineage', 'AA Substitutions', 'antigenic_score']
# Getting list of known sites of antigenicity:
tpSites = pd.read_csv(ref_dir + "tp_sites.csv", sep = ",")
tpSites_list = tpSites['tp_sites'].values.tolist()
tpSites_list = [str(item) for item in tpSites_list]
# Importing reference files for weight and amino acid properties
weights = pd.read_csv(ref_dir + "antigenic_weights.csv", sep = "\t")
aa_prop_df = pd.read_csv(ref_dir + "amino_acid_properties.csv", sep = "\t")

# Defining Necessary Functions

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

# Data Preparation

## Variant Dataframe - dataframe with variants, amino acid changes and their averaged antigenic score from January 2020 to August 2022

lineages = set() # Creating a unique set of all lineages from antigenic scoring results within given time frame (Jan2020 - Aug2022)
AAs = tuple(aa_prop_df["code"].tolist())
df = pd.DataFrame()
df_variant_antigenicity = pd.DataFrame()

    # Importing antigenic score results from Jan2020 to Aug2022 and getting amino acid changes per variant with assigned weight
for file in glob.glob(indir + "*.csv"):
    print(file)
    tmp = pd.read_csv(file, sep = "\t", usecols = columns)
    #####
    # Creating a dataframe with variants and their averaged antigenic scores used for regression later
    tmp_antigenic = tmp.drop(["Accession ID", "Collection date", "AA Substitutions"], axis = 1)
    tmp_antigenic = tmp_antigenic.drop_duplicates()
    df_variant_antigenicity = df_variant_antigenicity.append(tmp_antigenic)
    #print(pd.DataFrame.head(tmp_antigenic))
    #print(pd.DataFrame(tmp))
    #####
    tmp = tmp.drop(["antigenic_score"], axis = 1)
    tmp = tmp[tmp["AA Substitutions"].notna()]
    #for i in tmp["Pango lineage"].unique(): lineages.add(i) 
    tmp['AA Substitutions'] = tmp['AA Substitutions'].astype('string')
    substitutions = tmp['AA Substitutions'].apply(conversion_position)
    tmp["Substitution_Positions"] = substitutions
    original_AAs = tmp['AA Substitutions'].apply(conversion_originalAA)
    tmp["Original_aa"] = original_AAs
    changed_AAs = tmp['AA Substitutions'].apply(conversion_changedAA)
    tmp["Changed_aa"] = changed_AAs
    tmp = tmp[tmp.Substitution_Positions.map(len) > 0] #dropping all non-spike protein changes
    tmp = tmp.explode(['Substitution_Positions', 'Original_aa', 'Changed_aa']) # AA changes as individual rows
    # adding known antigenic weights to the aa positions:
    tmp_filtered = tmp[tmp["Substitution_Positions"].isin(tpSites_list)]
    tmp_filtered_weights = pd.merge(tmp_filtered, weights, how = "left", on = ["Original_aa", "Changed_aa"])
    tmp_final = pd.merge(tmp, tmp_filtered_weights, how = "left", on = ["Accession ID", "Substitution_Positions", "Original_aa", "Changed_aa"])
    print("tmp_final before: ")
    print(pd.DataFrame.head(tmp_final))
    tmp_final = tmp_final.drop(["Accession ID","Collection date_x", "AA Substitutions_x", "Collection date_y", "Pango lineage_y", "AA Substitutions_y"], axis = 1)
    tmp_final.rename(columns = {"Pango lineage_x": "Pango lineage"}, inplace = True)
    # Dropping duplicate amino acid changes (per variant) and adding each variant to lineages set
    tmp_final = tmp_final.drop_duplicates()
    for i in tmp_final["Pango lineage"].unique(): lineages.add(i)
    # adding tmp_final to complete dataframe
    df = df.append(tmp_final)
    print("tmp_final length: ", len(tmp_final))
    print(pd.DataFrame.head(tmp_final))

# Taking the average antigenic scores per variant across the designated time frame (Jan 2020 to Aug 2022)
df_variant_antigenicity['antigenic_score'] = df_variant_antigenicity.groupby('Pango lineage')['antigenic_score'].transform('mean')
variant_antigenicity = df_variant_antigenicity.to_records(index = False)
print("VARIANT ANTIGENICITY")
print(variant_antigenicity)
print(repr(variant_antigenicity))
# Final data cleaning of amino acid changes df
print("df length: ", len(df))
df = df.drop_duplicates()
print("filtered df length: ", len(df))
print(pd.DataFrame.head(df))
# Creating an ordered set of lineages that will be used as the baseline for the structured array with variants and aa changes
lineages = tuple(lineages)
print(lineages)
print("length of lineages: ", len(lineages))
print(AAs)
print("length of AAs: ", len(AAs))

## Matrix A - variants as rows and amino acids as columns

matrixA = np.zeros((len(lineages), len(AAs))) # Array in which the lineages will be the rows and amino acids the columns
print(matrixA)
print(matrixA.shape)

# Counting the number of amino acid changes in each lineage in comparison to reference strain (Wuhan)
for variant in lineages:
    print("Variant: ", variant)
    variant_df = df[df["Pango lineage"] == variant]
    variant_coord = lineages.index(variant)
    print(variant_coord)
    ##variant_coord = matrixA[np.isin(matrixA['lineage'], variant)]
    #variant_coord = matrixA[np.isin(matrixA['lineage'], variant)]
    #print("Variant_coord: ", variant_coord)
    #print(pd.DataFrame.head(variant_df))
    for i in variant_df["Original_aa"]:
        #print(i)
        if i in AAs:
            aa_coord = AAs.index(i)
            print(aa_coord)
            matrixA[variant_coord, aa_coord] = matrixA[variant_coord, aa_coord] - 1
            
            #aa_coord = AAs.index(i) + 1 # Added 1 to the index to take into account that lineage is index 0 in the array
            #print(aa_coord)
            #matrixA[variant_coord][aa_coord] = matrixA[variant_coord][aa_coord] - 1
            ##test_aa_coord = np.where(matrixA[variant_coord] == i)
            ##print("test_aa_coord: ", test_aa_coord)
            ##variant_coord[i] - 1
            ##print(variant_coord)
            ##matrixA[i][variant_coord] - 1
            ##aa_coord = matrixA["lineage"][variant]
            ##print(aa_coord)
            ##matrixA[variant, i] = matrixA[variant, i] - 1
        else:
            continue
    for i in variant_df["Changed_aa"]:
        if i in AAs:
            aa_coord = AAs.index(i)
            matrixA[variant_coord, aa_coord] = matrixA[variant_coord, aa_coord] + 1
            #matrixA[variant, i] = matrixA[variant, i] + 1
        else:
            continue
print(matrixA)
print(matrixA.shape)
#x = (matrixA > 0).sum()
#print(x)

#import numpy.lib.recfunctions as rf
#rf.unstructured_to_structured(matrixA, dtype = matrixA_def)
#print(matrixA)

x = len(lineages)
y = len(AAs)
matrixA = np.array(matrixA, dtype = str).reshape(x,y)
print(matrixA)
print(matrixA.shape)

tmp = np.array([])

for variant in lineages:
    variant_coord = lineages.index(variant)
    #matrixA[variant_coord][0] = variant
    a = np.insert(matrixA[variant_coord], 0, variant)
    print(a)
    tmp = np.append(tmp, a)
    #matrixA = np.insert(matrixA, variant_coord, variant)
tmp = tmp.reshape(len(lineages),(len(AAs)+1))
print(tmp)
print(tmp.shape)

#matrixA = np.array(matrixA, dtype = matrixA_def)
matrixA = np.core.records.fromarrays(tmp.transpose(), names = 'lineage, R, N, D, E, Q, K, S, T, C, H, M, A, V, G, I, L, F, P, W, Y', 
        formats = '<U10, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4, f4')
print(matrixA)
print(matrixA.shape)

#### To Do:
# Set up MatrixA as a two dimensional array rather than a one dimensional
# linear regression with scikit learn (https://datagy.io/python-sklearn-linear-regression/)

## Matrix B - amino acids as rows and properties as columns

#matrixB = pd.read_csv(ref_dir + "amino_acid_properties.csv", sep = "\t").to_records(index=False)
matrixB = aa_prop_df.to_records(index = False)
print(matrixB)

# Regression Models

## First regression (A(x) = antigenicity - vector) - impact of each amino acid on antigenicity

## Second regression