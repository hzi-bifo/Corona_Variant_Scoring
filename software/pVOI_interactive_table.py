#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import plotly.graph_objects as go

# Importing pVOI antigenic data
data = pd.read_csv(sys.argv[1], sep = ',')
output = sys.argv[2]

# Creating Figure
data['antigenic_score'] = data['antigenic_score'].round(3)
lineages =  list(data["Pango.lineage"])
scores = list(data["antigenic_score"])
fig = go.Figure(data=[go.Table(header=dict(values=['pVOI (Pango Lineage)', 'Antigenic Score']),
                 cells=dict(values=[lineages, scores]))])
# Saving Figure
fig.write_html(output + 'pVOI_interactive_table.html')

