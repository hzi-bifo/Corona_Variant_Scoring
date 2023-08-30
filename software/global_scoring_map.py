#!/usr/bin/env python
# coding: utf-8

import sys
import os
import pandas as pd
import plotly.express as px

df = pd.read_csv(sys.argv[1], sep = "\t")
output = sys.argv[2]
month_file = open(sys.argv[3])
reference_dir = sys.argv[4]

# Reading in current month
month_file.seek(0)
month = month_file.readline()
month_file.close()

# European only df for focused visualization:
df_eu = df[df['Continent'] == "Europe"]

# Creating global and european map of averaged antigenic scores
## Global Map (scaled 0-10)
labels_dict = dict(zip(df.Country, df.country_score))
fig = px.choropleth(df, locations = "Country",
                           locationmode = "country names",
                           color = "country_score",
                           color_continuous_scale = 'spectral_r',
                           labels = labels_dict,
			               range_color = [0, 10]) #,
#                           title = '<b>Country Antigenic Scores for %s</b>' %month)
#fig['layout']['title']['font'] = dict(size = 20)
fig.update_layout(geo = dict(showframe = False), coloraxis_colorbar = dict(title = "Country<br>Antigenic Score", dtick = 1, len = 0.75), font = dict(size = 14))
fig.add_annotation(x = 0, y = 0.97, text = '<b>Country Antigenic Scores for %s</b>' %month, showarrow = False)
fig.update_annotations(font = dict(size = 20))
fig.write_html(output + 'antigenic_score_map.html')
## European Map (unscaled)
fig_eu = px.choropleth(df_eu, locations = "Country",
                           width = 850,
                           locationmode = "country names",
                           color = "country_score",
                           color_continuous_scale = 'spectral_r',
                           labels = labels_dict,
                           scope = 'europe')
fig_eu.update_layout(geo = dict(showframe = False), coloraxis_colorbar = dict(title = "Country<br>Antigenic Score", dtick = 0.20, len = 0.75), 
        font = dict(size = 14))
fig_eu.add_annotation(x = 0, y = 0.87, text = '<b>Country Antigenic Scores for %s</b>' %month, showarrow = False)
fig_eu.update_annotations(font = dict(size = 20))
fig_eu.write_html(output + 'antigenic_score_map_europe.html')

# Creating global map with slider through time frame
## Checking to see if antigenic_scores_map_visualization_cumulative.csv exists, and if not creating one
df['date'] = month
if os.path.exists(reference_dir + "antigenic_scores_map_visualization_cumulative.csv") == False:
    df.to_csv(reference_dir + "antigenic_scores_map_visualization_cumulative.csv", sep=',', index=False, header=True)
    df_cumulative = df
else:
    df_cumulative = pd.read_csv(reference_dir + "antigenic_scores_map_visualization_cumulative.csv", sep = ',')

if not (df_cumulative['date'].eq(month)).any():
    df_cumulative = pd.concat([df_cumulative, df])
    df_cumulative.to_csv(reference_dir + "antigenic_scores_map_visualization_cumulative.csv", sep=',', index=False, header=True)
else:
    pass

labels_dict = dict(zip(df_cumulative.Country, df_cumulative.country_score))
fig = px.choropleth(df_cumulative, locations = "Country",
                           locationmode = "country names",
                           scope = "world",
                           animation_frame= "date",
                           color = "country_score",
                           color_continuous_scale = 'spectral_r',
                           labels = labels_dict,
			               range_color = [0, 10],
                           title = '<b>Country Antigenic Scores from 01-2020 to %s</b>' %month)
fig['layout']['title']['font'] = dict(size = 20)
fig.update_layout(coloraxis_colorbar = dict(title = "Country<br>Antigenic Score", dtick = 1, len = 0.75),
                  font = dict(size = 12)) #orientation = "h"
fig.write_html(output + 'antigenic_score_map_cumulative.html')
