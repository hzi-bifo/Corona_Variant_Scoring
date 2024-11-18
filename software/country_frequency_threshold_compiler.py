#!/usr/bin/env python
# coding: utf-8

# ## Country Frequency Threshold
# #### Author: Katrina Norwood
# #### Last Updated: 18/11/2024
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
antigenic_score_df = pd.read_csv(monthlyComparison_dir+"/antigenic_scores_all.csv", sep = '\t', usecols = ["Location"])
country_score_df = pd.read_csv(monthlyComparison_dir+"/antigenic_scores_map_visualization.csv", sep = "\t", usecols = ["Country", "country_score"])

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
        
# Getting number of isolates and their frequency per country
antigenic_score_countries_df = antigenic_score_df.value_counts(['Country']).reset_index().rename(columns={0:'num_isolates'})
antigenic_score_countries_df["seq_frequency"] = antigenic_score_countries_df["num_isolates"] / num_sequences

# Creating a list of the countries with greater than 1% representation
print("Creating list of countries that meet threshold: ")
freq_countries_monthly_df = antigenic_score_countries_df[(antigenic_score_countries_df['seq_frequency'] >= 0.01) | (antigenic_score_countries_df['num_isolates'] >= 500)]

# Creating a df with the high freq countries and their country antigenic score for that month
print("Adding date column and dropping duplicates: ")
# Adding date column
print(pd.DataFrame.head(freq_countries_monthly_df))
freq_countries_monthly_df['date'] = str(month)
# Dropping duplicates
print("dropping duplicates")
freq_countries_monthly_df = freq_countries_monthly_df.drop_duplicates()

# Checking to see if the month has already been added to the cumulative file, and if not, adding it
print("Adding new countries to the cumulative file: ")
if not (df_cumulative['date'].eq(month)).any():
    df_cumulative = pd.concat([df_cumulative, freq_countries_monthly_df])
    df_cumulative.to_csv(reference_dir + "country_list_with_threshold.tsv", sep='\t', index=False, header=True)
else:
    print("Already added month, passing...")
    pass
