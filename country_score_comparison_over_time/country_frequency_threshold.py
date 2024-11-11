#!/usr/bin/env python
# coding: utf-8

# ## Country Frequency Threshold
# #### Author: Katrina Norwood
# #### Last Updated: 30/10/2024
# 
# Script to identify which countries have at least 1% of the monthly sequences over the course of the three year analysis period (01-2020 to 12-2023). These countries will then be used for the country_score_comparison_over_time analysis. 

import os
import sys
import pandas as pd
import numpy as np

output = sys.argv[1]
monthlyComparison_dir = sys.argv[2]

print("output: ", output)
print("monthly comparison dir: ", monthlyComparison_dir)

print(output)

final_df = pd.DataFrame()

for directory in os.listdir(monthlyComparison_dir):
    # Parsing through the month directories to find countries that had at least 1% of sequences
    print("Directory: ", directory)
    if "-20" in str(directory):
        antigenic_score_df = pd.read_csv(monthlyComparison_dir+directory+"/antigenic_scores_all.csv", sep = '\t', usecols = ["Location"])
        country_score_df = pd.read_csv(monthlyComparison_dir+directory+"/antigenic_scores_map_visualization.csv", sep = "\t", usecols = ["Country", "country_score"])
        # Checking to see which countries have a frequency of greater than 1% of the global submitted sequences
        num_sequences = len(antigenic_score_df)
        antigenic_score_df['Country'] = antigenic_score_df['Location'].apply(lambda x: x.split("/")[1])
        # Some countries have spaces following them, here we need to remove that to get a proper count
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
        freq_countries_monthly_df = antigenic_score_countries_df[(antigenic_score_countries_df['seq_frequency'] >= 0.01) | (antigenic_score_countries_df['num_isolates'] >= 500)]
        # Creating a df with the high freq countries and their country antigenic score for that month
        freq_countries_score_df = freq_countries_monthly_df.merge(country_score_df, how = "left", on = "Country")
        # Adding date column
        freq_countries_score_df["date"] = str(directory)
        # Adding to final dataframe
        final_df = final_df.append(freq_countries_score_df, ignore_index=True)
    else:
        continue

print(final_df)

final_df.to_csv(output + "country_antigenic_score_with_threshold.tsv", sep = "\t", index = False, header = True)

