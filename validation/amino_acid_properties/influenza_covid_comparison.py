#!/usr/bin/env python
# coding: utf-8

#### Amino Acid Properties Script
#Author: Katrina Norwood
#Last Updated: 31/01/25

# Script to create a dataframe of influenza true positive sites and SARS-CoV-2 true positive sites from literature

import os
import itertools
import pandas as pd

columns = ["Original_aa","Changed_aa","Hydropathy","Charge","Molecular Weight Change"]
output = "/validation/amino_acid_properties/"
aa_weights = pd.read_csv("/reference/antigenic_weights.csv", sep = "\t", header = 0)
aa_properties = pd.read_csv("/reference/amino_acid_properties.csv", sep = "\t", header = 0)

# Creating a dataframe of all amino acid combinations (but minus the changes to themselves ie an A to A change)
amino_acids = aa_properties["code"] 

# Finding all aa combinations
aa_combinations = list(itertools.product(amino_acids, repeat=2))
aa_combinations_filtered = [pair for pair in aa_combinations if pair[0] != pair[1]]
aa_changes_df = pd.DataFrame(aa_combinations_filtered, columns=['Original_aa', 'Changed_aa'])

# Adding Influenza and SARS-CoV-2 properties to the dataframe

aa_weights["Influenza"] = "True"
merged_df = pd.merge(aa_changes_df, aa_weights, on=['Original_aa', 'Changed_aa'], how='left')
merged_df

# Adding molecular properties

# Creating a dataframe with the aa changes and their properties
original_aa_properties = aa_properties.copy()
original_aa_properties.rename(columns={'amino_acid': 'original_amino_acid', 'code': 'Original_aa', 'hydropathy':'original_hydropathy',
                                      'charge':'original_charge', 'pKa_NH2':'original_pKa_NH2', 'pKa_COOH':'original_pKa_COOH',
                                      'pKR':'original_pKR', 'solubility':'original_solubility', 'molecular_weight':'original_molecular_weight'}, inplace=True)
changed_aa_properties = aa_properties.copy()
changed_aa_properties.rename(columns={'amino_acid': 'changed_amino_acid', 'code': 'Changed_aa', 'hydropathy':'changed_hydropathy',
                                      'charge':'changed_charge', 'pKa_NH2':'changed_pKa_NH2', 'pKa_COOH':'changed_pKa_COOH',
                                      'pKR':'changed_pKR', 'solubility':'changed_solubility','molecular_weight':'changed_molecular_weight'}, inplace=True)
aa_properties_weights = pd.merge(merged_df, original_aa_properties, how="left", on=["Original_aa"])
aa_properties_weights = pd.merge(aa_properties_weights, changed_aa_properties, how="left", on=["Changed_aa"])

# Creating Changed category columns
def hydropathy(row):
    if row["original_hydropathy"] == row["changed_hydropathy"]:
        return "unchanged" 
    return f"{row['original_hydropathy']} - {row['changed_hydropathy']}"

def charge(row):
    charge_mapping = {
        "N": "Neutral",
        "+": "Positive",
        "-": "Negative"}
    
    if row["original_charge"] == row["changed_charge"]:
        return "unchanged" 
    
    original_charge = charge_mapping.get(row["original_charge"], row["original_charge"])
    changed_charge = charge_mapping.get(row["changed_charge"], row["changed_charge"])
    
    return f"{original_charge} - {changed_charge}"

def size(row):
    if row["original_molecular_weight"] == row["changed_molecular_weight"]:
        return "unchanged" 
    return row['original_molecular_weight'] - row['changed_molecular_weight']

aa_properties_weights["Hydropathy"] = aa_properties_weights.apply(hydropathy, axis=1)
aa_properties_weights["Charge"] = aa_properties_weights.apply(charge, axis=1)
aa_properties_weights["Molecular Weight Change"] = aa_properties_weights.apply(size, axis=1)

# Creating final dataframe
aa_properties_final = aa_properties_weights[["Original_aa","Changed_aa","Weight","Influenza","Hydropathy","Charge","Molecular Weight Change"]]
aa_properties_final

# Saving Dataframe
aa_properties_final.to_csv(output+"influenza_covid_comparision_aa_properties_final.tsv", sep = "\t", header = True, index = False)