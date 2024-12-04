#!/usr/bin/env python
# coding: utf-8

# ## Country Frequency Threshold
# #### Author: Katrina Norwood
# #### Last Updated: 19/11/2024
# 
# Script to identify which countries have at least 1% of the monthly sequences (mostly for the earlier months ie. Jan 2020) or have 500 sequences. These countries will then be used for the country_score_comparison_over_time analysis. 

import os
import sys
import pandas as pd
import numpy as np

output = sys.argv[1]
monthlyComparison_dir = sys.argv[2]
reference_dir = sys.argv[3]
month_file = open(sys.argv[4])

# Reading in current month
month_file.seek(0)
month = month_file.readline()
month_file.close()

print("output: ", output)
print("monthly comparison dir: ", monthlyComparison_dir)
print("month: ", month)

# Reading in cumulative file and if no data exists quitting script
if os.path.exists(reference_dir + "country_list_with_threshold.tsv") == True:
    df_cumulative = pd.read_csv(reference_dir + "country_list_with_threshold.tsv", sep='\t')
    print("Cumulative DF File: ")
    print(pd.DataFrame.head(df_cumulative))
else:
    print("No cumulative country-wise file, please check! Now exiting...")
    quit()

# Parsing through the month directories to find countries that had at least 1% of sequences
print("Reading In Input DFs: ")
antigenic_score_df = pd.read_csv(monthlyComparison_dir+"/antigenic_scores_all.csv", sep = '\t', usecols = ["Location", "Pango lineage"])
country_score_df = pd.read_csv(monthlyComparison_dir+"/antigenic_scores_map_visualization.csv", sep = "\t", usecols = ["Country", "country_score"])
zscore_df = pd.read_csv(monthlyComparison_dir+"/antigenic_scores_ranked_with_WHO.csv", sep = "\t", usecols = ["Pango lineage", "antigenic_score", "zscore"])

# Checking to see which countries have a frequency of greater than 1% of the global submitted sequences
print("Adding number of seqeunces and parsing out country information: ")
num_sequences = len(antigenic_score_df)
antigenic_score_df['Country'] = antigenic_score_df['Location'].apply(lambda x: x.split("/")[1])
        
# Some countries have spaces following them, here we need to remove that to get a proper count
print("Formatting of country names: ")
country = []
for item in antigenic_score_df["Country"]:
    x = item.lstrip()
    y = x.rstrip()
    country.append(y)
antigenic_score_df["Country"] = country
print(pd.DataFrame.head(antigenic_score_df))

# Getting number of isolates and their frequency per country (to identify countries that match the threshold)
print("Getting number of isolates per country: ")
antigenic_score_df["num_isolates"] = antigenic_score_df.groupby("Country")["Country"].transform("count")
antigenic_score_df["seq_frequency"] = antigenic_score_df["num_isolates"] / num_sequences

# Calculating frequency per pango lineage 
print("Calculating lineage frequency: ")
antigenic_score_df["num_lineage_isolates"] = antigenic_score_df.groupby(['Country', 'Pango lineage'])['Pango lineage'].transform("count")
antigenic_score_df["lineage_frequency"] = antigenic_score_df["num_lineage_isolates"] / antigenic_score_df["num_isolates"]
print("After adding lineage and lineage frequency: ")
print(pd.DataFrame.head(antigenic_score_df))

# Merging antigenic score and zscore
print("Merging antigenic score and zscore: ")
antigenic_score_countries_df = antigenic_score_df.merge(zscore_df, how='left', on='Pango lineage')
print(pd.DataFrame.head(antigenic_score_countries_df))

# Creating a list of the countries with greater than 1% representation
print("Creating list of countries that meet threshold: ")
freq_countries_monthly_df = antigenic_score_countries_df.drop(['Location'], axis = 1)
freq_countries_monthly_df = freq_countries_monthly_df[(freq_countries_monthly_df['seq_frequency'] >= 0.01) | (freq_countries_monthly_df['num_isolates'] >= 500)]

# Dropping duplicates
print("Dropping duplicates: ")
freq_countries_monthly_df = freq_countries_monthly_df.drop_duplicates()
print(pd.DataFrame.head(freq_countries_monthly_df, n = 25))

# Keeping lineages with the highest frequency per country (country score = antigenic score x frequency)
print("Keeping most frequent lineages per country: ")
freq_countries_monthly_df = freq_countries_monthly_df.groupby('Country', group_keys=False).apply(
    lambda x: x.nlargest(5, 'lineage_frequency')
).reset_index(drop=True)
# Cleaning up decimal places
freq_countries_monthly_df["lineage_frequency"] = freq_countries_monthly_df["lineage_frequency"].round(3)
freq_countries_monthly_df["antigenic_score"] = freq_countries_monthly_df["antigenic_score"].round(3)
freq_countries_monthly_df["zscore"] = freq_countries_monthly_df["zscore"].round(3)
print(pd.DataFrame.head(freq_countries_monthly_df, n = 25))

# Putting together the lineage information
print("Putting together lineage information per country: ")
freq_countries_monthly_df['lineages_information'] = freq_countries_monthly_df.apply(lambda row: {
    'Pango lineage': row['Pango lineage'],
    'lineage_frequency': row['lineage_frequency'],
    'antigenic_score': row['antigenic_score'],
    'zscore': row['zscore']}, axis=1)
print(pd.DataFrame.head(freq_countries_monthly_df))
lineage_information = freq_countries_monthly_df.groupby('Country')['lineages_information'].apply(list).reset_index()
print(lineage_information)
# Removing redundant columns
freq_countries_monthly_df = freq_countries_monthly_df.drop(['Pango lineage', 'lineage_frequency', 'zscore', 'num_lineage_isolates', 'antigenic_score', 'lineages_information'], axis = 1)
# Dropping duplicates (so theres a dataframe of country and country info only)
freq_countries_monthly_df = freq_countries_monthly_df.drop_duplicates()
# Adding dictionary of information to df
freq_countries_monthly_df = freq_countries_monthly_df.merge(lineage_information, how = "left", on = "Country")
print(pd.DataFrame.head(freq_countries_monthly_df))

# Creating a df with the high freq countries and their country antigenic score for that month
print("Adding date column: ")
# Adding date column
print(pd.DataFrame.head(freq_countries_monthly_df))
freq_countries_monthly_df['date'] = str(month)

# Checking to see if the month has already been added to the cumulative file, and if not, adding it
print("Adding new countries to the cumulative file: ")
if not (df_cumulative['date'].eq(month)).any():
    df_cumulative = pd.concat([df_cumulative, freq_countries_monthly_df])
    df_cumulative.to_csv(reference_dir + "country_list_with_threshold.tsv", sep='\t', index=False, header=True)
else:
    print("Already added month, passing...")
    pass

freq_countries_monthly_df.to_csv(reference_dir + "country_list_with_threshold.tsv", sep='\t', index=False, header=True)
