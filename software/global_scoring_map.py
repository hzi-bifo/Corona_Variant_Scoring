#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import plotly.express as px

df = pd.read_csv(sys.argv[0], sep = "\t")
output = sys.argv[1]

# Keeping variants that have an average mutation score greater than the assigned VOC threshold
df_vis_threshold = df_vis[df_vis["antigenic_score"] > 2.67]
df_vis_threshold['averaged_antigenic_score_per_country'] = df_vis_threshold.groupby('Country')['antigenic_score'].transform('mean')
df_vis_threshold_averaged = df_vis_threshold[['Continent','Country', 'averaged_antigenic_score_per_country']]
df_vis_threshold_averaged = df_vis_threshold_averaged.drop_duplicates()
df_vis_threshold_averaged

# European only df for focused visualization:
df_vis_threshold_averaged_eu = df_vis_threshold_averaged[df_vis_threshold_averaged['Continent'] == "Europe"]
df_vis_threshold_averaged_eu

# Creating global and european map of averaged antigenic scores
labels_dict = dict(zip(df_vis_threshold_averaged.Country, df_vis_threshold_averaged.averaged_antigenic_score_per_country))
fig = px.choropleth(df_vis_threshold_averaged, locations = "Country",
                           locationmode = "country names",
                           color = "averaged_antigenic_score_per_country",
                           color_continuous_scale = 'spectral_r',
                           labels = labels_dict,
                           title = 'Antigenic Scores per Country')
fig.write_html(output + 'antigenic_score_map.html')
fig_eu = px.choropleth(df_vis_threshold_averaged_eu, locations = "Country",
                           locationmode = "country names",
                           color = "averaged_antigenic_score_per_country",
                           color_continuous_scale = 'spectral_r',
                           labels = labels_dict,
                           title = 'Antigenic Scores per Country',
                           scope = 'europe')
fig_eu.write_html(output + 'antigenic_score_map_europe.html')
