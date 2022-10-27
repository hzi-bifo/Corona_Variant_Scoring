#!/usr/bin/env python
# coding: utf-8
import sys
import pandas as pd
import plotly.express as px

df = pd.read_csv(sys.argv[1], sep = "\t")
output = sys.argv[2]

#def visuals (df, month, year, output):
#    df_eu = df[df['Continent'] == "Europe"]
#    months = str(month) + "-" + str(max)
#    labels_dict = dict(zip(df.Country, df.country_score))
#    fig = px.choropleth(df, locations = "Country",
#                               locationmode = "country names",
#                               color = "country_score",
#                               color_continuous_scale = 'spectral_r',
#                               labels = labels_dict,
#                               title = 'Antigenic Scores per Country for %s' %months)
#    fig.write_html(output + 'antigenic_score_map.html')
#    fig_eu = px.choropleth(df_eu, locations = "Country",
#                               locationmode = "country names",
#                               color = "country_score",
#                               color_continuous_scale = 'spectral_r',
#                               labels = labels_dict,
#                               title = 'Antigenic Scores per Country for %s' %months,
#                               scope = 'europe')
#    fig_eu.write_html(output + 'antigenic_score_map_europe.html')
    
def df_selection (df, month, year):
# Function to select given months from sequencing data
    df_vis = df[df.month == month]
    df_vis = df_vis[df_vis.year == year]
    
    # Calculating the frequency for each lineage per country:
    df_vis = df_vis[df_vis['Pango lineage'] != 'None']
    n_lineages = pd.Series(df_vis.groupby(['Country', 'Pango lineage'])['Pango lineage'].count(), name = "n_lineages")
    n_lineages_df = n_lineages.to_frame().reset_index()
    n_lineages_df.ffill(axis = 0, inplace = True)
    df_vis = df_vis.merge(n_lineages_df, how = 'left', on = ['Country','Pango lineage'])
    total_lineages = pd.Series(df_vis.groupby('Country')['Country'].count(), name = 'total_lineages')
    df_vis = df_vis.merge(total_lineages.to_frame(), how = 'left', on = 'Country')
    df_vis['frequency'] = df_vis['n_lineages'].div(df_vis['total_lineages'])
    
    # Filtering duplicates & assigning a country score:
    df_vis.drop(['Collection date', 'collection_date_list', 'n_lineages', 'total_lineages'], axis = 1, inplace = True)
    df_vis_threshold_averaged = df_vis.drop_duplicates()
    df_vis_threshold_averaged['score'] = df_vis_threshold_averaged['antigenic_score']*df_vis_threshold_averaged['frequency']
    df_vis_threshold_averaged['country_score'] = df_vis_threshold_averaged.groupby('Country')['score'].transform('sum')

    # Saving dataframe
    df_vis_threshold_averaged.to_csv(output + "country_scores" + month + "_" + year + ".csv", sep = '\t', 
                                     index = False, header = True)

    return df_vis_threshold_averaged

# Creating country and continent column for later analysis
df_vis = df[['Location', 'Pango lineage', 'antigenic_score', 'Collection date']]
df_vis["Continent"] = df_vis["Location"].apply(lambda x: x.split("/")[0])
df_vis["Country"] = df_vis["Location"].apply(lambda x: x.split("/")[1])
country = []
continent = []
for item in df_vis["Country"]:
    x = item.lstrip()
    y = x.rstrip()
    country.append(y)
for item in df_vis["Continent"]:
    x = item.lstrip()
    y = x.rstrip()
    continent.append(y)
df_vis['Country'] = country
df_vis['Continent'] = continent
df_vis.drop(["Location"], axis = 1, inplace = True)

# Adding a year and month column to filter data by most current month
df_vis['year'] = df_vis["Collection date"].apply(lambda x: x.split("-")[0])
df_vis['collection_date_list'] = df_vis["Collection date"].apply(lambda x: x.split("-"))
keep_list = []
for item in df_vis["collection_date_list"]:
    if len(item) <= 1:
        keep_list.append("N")
    else:
        keep_list.append("Y")
df_vis['keep'] = keep_list
df_vis = df_vis[df_vis.keep != "N"]
df_vis['month'] = df_vis["Collection date"].apply(lambda x: x.split("-")[1])

# Running for each of the past months
print("CREATING DATAFRAMES:")
aug2022_df = df_selection(df_vis, "08", "2022")
print("AUGUST 2022 DF COMPLETE")
jul2022_df = df_selection(df_vis, "07", "2022")
print("JULY 2022 DF COMPLETE")
jun2022_df = df_selection(df_vis, "06", "2022")
print("JUNE 2022 DF COMPLETE")
may2022_df = df_selection(df_vis, "05", "2022")
print("MAY DF COMPLETE")
apr2022_df = df_selection(df_vis, "04", "2022")
print("APRIL 2022 DF COMPLETE")
mar2022_df = df_selection(df_vis, "03", "2022")
print("MARCH 2022 DF COMPLETE")
feb2022_df = df_selection(df_vis, "02", "2022")
print("FEB 2022 DF COMPLETE")
jan2022_df = df_selection(df_vis, "01", "2022")
print("JAN 2022 DF COMPLETE")
dec2021_df = df_selection(df_vis, "12", "2021")
print("DEC 2021 DF COMPLETE")
nov2021_df = df_selection(df_vis, "11", "2021")
print("NOV 2021 DF COMPLETE")
oct2021_df = df_selection(df_vis, "10", "2021")
print("OCT 2021 DF COMPLETE")
sept2021_df = df_selection(df_vis, "09", "2021")
print("SEPT 2021 DF COMPLETE")
aug2021_df = df_selection(df_vis, "08", "2021")
print("AUG 2021 DF COMPLETE")
jul2021_df = df_selection(df_vis, "07", "2021")
print("JUL 2021 DF COMPLETE")
jun2021_df = df_selection(df_vis, "06", "2021")
print("JUN 2021 DF COMPLETE")
may2021_df = df_selection(df_vis, "05", "2021")
print("MAY 2021 DF COMPLETE")
apr2021_df = df_selection(df_vis, "04", "2021")
print("APRIL 2021 DF COMPLETE")
mar2021_df = df_selection(df_vis, "03", "2021")
print("MARCH 2021 DF COMPLETE")
feb2021_df = df_selection(df_vis, "02", "2021")
print("FEB 2021 DF COMPLETE")
jan2021_df = df_selection(df_vis, "01", "2021")
print("JAN 2021 DF COMPLETE")
dec2020_df = df_selection(df_vis, "12", "2020")
print("DEC 2020 DF COMPLETE")
nov2020_df = df_selection(df_vis, "11", "2020")
print("NOV 2020 DF COMPLETE")
oct2020_df = df_selection(df_vis, "10", "2020")
print("OCT 2020 DF COMPLETE")
sept2020_df = df_selection(df_vis, "09", "2020")
print("SEPT 2020 DF COMPLETE")
aug2020_df = df_selection(df_vis, "08", "2020")
print("AUG 2020 DF COMPLETE")
jul2020_df = df_selection(df_vis, "07", "2020")
print("JUL 2020 DF COMPLETE")
jun2020_df = df_selection(df_vis, "06", "2020")
print("JUN 2020 DF COMPLETE")
may2020_df = df_selection(df_vis, "05", "2020")
print("MAY 2020 DF COMPLETE")
apr2020_df = df_selection(df_vis, "04", "2020")
print("APRIL 2020 DF COMPLETE")
mar2020_df = df_selection(df_vis, "03", "2020")
print("MAR 2020 DF COMPLETE")
feb2020_df = df_selection(df_vis, "02", "2020")
print("FEB 2020 DF COMPLETE")
jan2020_df = df_selection(df_vis, "01", "2020")
print("JAN 2020 DF COMPLETE")
