#!/usr/bin/env python
# coding: utf-8

#### EVEscape Methods Comparison
#Author: Katrina Norwood
#Last Updated: 06/11/24

# Script to copy accession IDs from the GISAID data as a .csv file for 
# reproducibility.

# To Run:
# In the script directory run '''python accession_ids </path to results dir> <output dir> <filename>'''
# The results dir is the directory with the Corona_Variant_Scoring monthly results (keep this in mind) and the
# output_dir is where you want the zip file saved while the filename is the name of the file you want saved.

import sys
import os
import pandas as pd

results_dir = sys.argv[1]
output = sys.argv[2]
filename = sys.argv[3]

# Concatenating all the IDs into one file

ids = pd.DataFrame()

for month_dir in os.listdir(results_dir):
    # Parsing through the results directory to pull accession ids from each month
    print(month_dir)
    if '-202' in month_dir:
        print("Parsing through: ", month_dir)
        df = pd.read_csv(results_dir+month_dir+"/antigenic_scores_all.csv", sep = "\t", usecols = ["Accession ID"])
        print("Number of IDs: ", len(df))
        ids = pd.concat([ids, df], ignore_index = True)
        print("Total Number of IDs: ", len(ids))

print(pd.DataFrame.head(ids))

# Saving as a zipped file
ids.to_csv(output+filename+'.csv.gz', index=False, compression='gzip')

