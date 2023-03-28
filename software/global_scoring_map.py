#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import plotly.express as px

df = pd.read_csv(sys.argv[1], sep = "\t")
output = sys.argv[2]
month_file = open(sys.argv[3])

# Reading in current month
month_file.seek(0)
month = month_file.readline()
month_file.close()

# European only df for focused visualization:
df_eu = df[df['Continent'] == "Europe"]

# Creating global and european map of averaged antigenic scores
labels_dict = dict(zip(df.Country, df.country_score))
fig = px.choropleth(df, locations = "Country",
                           locationmode = "country names",
                           color = "country_score",
                           color_continuous_scale = 'spectral_r',
                           labels = labels_dict,
			               range_color = [0, 10],
                           title = '<b>Regional Antigenic Scores for %s</b>' %month)
fig['layout']['title']['font'] = dict(size = 20)
fig.update_layout(coloraxis_colorbar = dict(title = "Regional Antigenic Score", dtick = 1),
                  font = dict(size = 18)) #orientation = "h"
fig.write_html(output + 'antigenic_score_map.html')

fig_eu = px.choropleth(df_eu, locations = "Country",
                           locationmode = "country names",
                           color = "country_score",
                           color_continuous_scale = 'spectral_r',
                           labels = labels_dict,
                           title = '<b>Regional Antigenic Scores for %s</b>' %month,
                           scope = 'europe')
fig['layout']['title']['font'] = dict(size = 20)
fig.update_layout(coloraxis_colorbar = dict(title = "Regional Antigenic Score", dtick = 1))
fig_eu.write_html(output + 'antigenic_score_map_europe.html')
