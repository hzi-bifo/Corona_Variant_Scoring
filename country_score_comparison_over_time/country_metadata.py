#!/usr/bin/env python
# coding: utf-8

#### Country Score Comparison Over Time
#Author: Katrina Norwood
#Last Updated: 29/11/24

# Script to put together the metadata (including lineages circulating and their frequencies) for the country-wise
# analysis. This script requires a directory with the results of the antigenic alterations
# pipeline saved in indivudal months directories (format MM-YYYY). It then goes through
# those directories for the antigenic_scores_all.csv to find and calculate circulating
# lineages and their frequencies (based on the sequences in the country - and ignoring
# sequences under review since these are removed in the analysis). This file is then saved
# to be used as a reference for the country-wise analysis done for CoVerage.

import sys
import os
import pandas as pd

# Importing the necessary directory
results_directory = sys.argv[1]
countries_df = pd.read_csv(sys.argv[2], sep = "\t")
output = sys.argv[2]

print(pd.DataFrame.head(countries_df))

final_df = pd.DataFrame()

for folder in os.listdir(results_directory):
    print(folder)
    if "-20" in str(folder):
        print("Month: ", folder)
        
        # Importing necessary DFs
        antigenic_df = pd.read_csv(results_directory+folder+"/antigenic_scores_all.csv", 
                                   sep = "\t", usecols = ["Location", "Pango lineage", "antigenic_score"])
        countries_month_df = countries_df[countries_df["date"] == folder]
        
        # Adding a country column to results df
        antigenic_df['Country'] = antigenic_df['Location'].apply(lambda x: x.split("/")[1])
        # Some countries have spaces following them, here we need to remove that to get a proper count
        country = []
        for item in antigenic_df["Country"]:
            x = item.lstrip()
            y = x.rstrip()
            country.append(y)
        antigenic_df["Country"] = country
        
        # Calculating lineage frequency for countries in the countries_df file
        countries = list(countries_month_df["Country"].unique())
        filtered_antigenic_df = antigenic_df[antigenic_df["Country"].isin(countries)]
        
        #country_antigenic_df = pd.DataFrame()
        for country in countries:
            print(country)
            country_antigenic_df_temp = filtered_antigenic_df[filtered_antigenic_df["Country"] == country]
            num_country_sequences = len(country_antigenic_df_temp)
            #print("Total number of sequences: ", num_country_sequences)
            country_antigenic_df_temp = country_antigenic_df_temp.value_counts(['Pango lineage']).reset_index().rename(columns={0:'lineage_num_isolates'})
            country_antigenic_df_temp["lineage_frequency"] = country_antigenic_df_temp["lineage_num_isolates"] / num_country_sequences
            country_antigenic_df_temp["Country"] = country
            country_antigenic_df_temp["date"] = folder
            #print(pd.DataFrame.head(country_antigenic_df_temp))
            final_df = final_df.append(country_antigenic_df_temp, ignore_index = True)
    else:
        print("Not a results directory")

print(pd.DataFrame.head(final_df))

final_df.to_csv(output + "country_score_metadata.tsv", sep = "\t", index = False, header = True)