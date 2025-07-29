#!/usr/bin/env python
# coding: utf-8

# Script to add new antigenically altered lineages to a compiled list of antigenic lineages (from 2025-10-02 to now)
# Created with ChatGPT
# Last Updated: 29/07/2025

import sys
import csv
import json
from datetime import datetime

# Input files
compiled_tsv = sys.argv[1]
new_json = sys.argv[2]
output_tsv = sys.argv[3]

# Helper to parse dates
def parse_date(d):
    try:
        return datetime.strptime(d, "%Y-%m-%d")
    except:
        return None

# Today's date for new JSON entries
run_date = datetime.today().strftime("%Y-%m-%d")

# Load existing TSV into a dict
lineage_records = {}
with open(compiled_tsv, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        lineage = row["lineage"]
        lineage_records[lineage] = row

# Load new JSON
data = []
with open(new_json, "r") as f:
    try:
        data = json.load(f)
    except ValueError:
        print("Invalid JSON file: %s" % new_json)

# Add new entries if lineage not already in records
for entry in data:
    if entry.get("significant", "").lower() == "yes":
        lineage = entry.get("lineage", "NA")
        if lineage not in lineage_records:
            lineage_records[lineage] = {
                "lineage": lineage,
                "antigenic_score": entry.get("antigenic_score", "NA"),
                "normalized_antigenic_score": entry.get("zscore", "NA"),
                "significantly_antigenically_altered": entry.get("significant", "NA"),
                "run_date": run_date
            }

# Sort by run_date ascending
compiled_data = sorted(lineage_records.values(), key=lambda x: parse_date(x["run_date"]))

# Write back to TSV
with open(output_tsv, "w") as out_f:
    fieldnames = ["lineage", "antigenic_score", "normalized_antigenic_score", "significantly_antigenically_altered", "run_date"]
    writer = csv.DictWriter(out_f, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for row in compiled_data:
        writer.writerow(row)

print("Appended new unique lineages from %s to %s with today's date %s" % (new_json, output_tsv, run_date))