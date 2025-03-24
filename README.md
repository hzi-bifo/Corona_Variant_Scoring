# Corona Variant Prediction

Pipeline to assign an antigenic score to pangolin lineages based on the amino acid changes across the SARS-CoV-2 spike glycoprotein. The score is based on the summation of antigenic weights from the influenza dataset (available here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3330098/) that are then averaged across the different positions that they occured. The repo includes data and results from the current preprint (https://www.biorxiv.org/content/10.1101/2024.03.07.583829v2) under the validation, country_score_comparison_over_time, SpikePro_comparison, and EVEscape_comparison directories. 

## Installation

The analysis pipeline is designed to run on BIFO servers, but can be run individually as well using the provided environment, sarscoverage.yml. The sarscoverage.yml can be used on linux and the sarscoverage_without_builds.yml file can be used on MacOS or Windows as it doesn't include package builds. 

To install the environment run the following comand:
```console
conda env create -f env_local.yml
```
Then the environment can be activated with:
```console
conda activate aa_scoring
```

## Required Pipeline Inputs

The Corona_Variant_Scoring pipeline requires the following inputs:

- GISAID metadata file: labeled as metadata.tsv, you will need to provide the path to this file
- GISAID metadata file of isolates under review: this is the metadata file for all sequences under review by GISAID, these will be removed from the analysis
- SD Plots frequencies directory: this is a specific output run with a concurrent analysis for CoVerage, but it contains the frequencies of all circulating lineages in each country and their relative frequencies for all months
- Months text file: file which list all months from the onset of SARS-CoV-2 to the desired month of analysis

The pipeline will require the path to the frequencies directory and months file, both of which are outputs from the SD Plots analysis. There are examples of what these files look like in the test_run directory. The months file contains a list of all months up to the desired month of analysis, while the frequencies file contains frequencies for each of the lineages circulating in individual countries for all months. These two files are necessary for the heatmap visual and as such the pipeline can be run without them if the heatmap visual is not run. Simply comment out the heatmap script from the variant_scoring.sh as seen below:

```
# Rscript "$SOFTWAREPATH""frequency_heatmap_coverage.R" "$FREQUENCY" "$OUTDIR""output/antigenic_scores_ranked_with_WHO.csv" "$MONTHS" "0.1" "$OUTDIR""output/" "$OUTDIR""output/antigenic_scores_all.csv" >> "$OUTDIR""STDOUT.txt"
```

## Running the pipeline

To run the pipeline you will need to specify the output and input directory, path to the Corona_Variant_Scoring repo, path to frequencies data, and the paths to a months file. 

Run the following command to start the pipeline:

```console
$bash variant_scoring_local.sh \
  -o /path to output directory file/ \
  -i /path to input directory/ 
  -v /path to Corona_Variant_Scoring repo/ \
  -f /path to frequencies data/ \
  -m /path to months text file/ 
  -u /path to sequences under review file/
```

The user also has the option of using -q and -w to add a specified month and year respectively for analysis rather than running it on the previous month as it is set to do for the CoVerage website. Each input file is outlined below:

- **(o) directory file:** this is the desired directory where all the output (as discussed below) for the antigenic scoring analysis will be saved.
- **(i) input directory:** this is the directory that contains all the input data, which is the metadata.tsv file from GISAID
  (a metadata file containing SARS-CoV-2 isolates and their respective spike protein changes)
- **(v) Corona_Variant_Scoring repo:** this repo, the analysis script will use different scripts in the software directory
- **(f) frequencies data:** an output of the SD Plots pipeline that contains the frequencies of different circulating lineages,
used in the heatmap visualization
- **(m) months text file:** also an output of the SD Plots pipeline that contains a list of the previous months of analysis 
as well as the current month
- **(u) sequences under review file:** a .csv file containing a list of the isolate IDs that are under review by GISAID for
that current month, these will need to be removed in the analysis

**Use bash variant_scoring_local.sh --help for more details.**

## Analysis Outputs

The pipeline will output two .csv files:

- antigenic_scores_all.csv : file containing all sequences per country and their associated antigenic scores, raw file in which scores have not been averaged per pango lineage\
- antigenic_scores_ranked_with_WHO.csv : file containing the averaged antigenic scores across all sequences per Pango lineage, these scores are than ranked to show the lineages with the highest antigenic scores
- antigenic_scores_map_visualization.csv : file containg the country antigenic score for all unique countries in the GISAID metadata file for the desired analysis month
- antigenic_scoring_summary_lineages_table.json : json file with all lineages, their zscores (whether significantly altered or not), and mutations
- antigenic_scoring_summary_states_lineages_table.json : same as above but for the German states only
- antigenic_scoring_summary_states_statistics.csv : circulating lineages and their frequencies for the German states
- antigenic_scoring_summary_statistics.csv : global circulating lineages and their frequencies
- month_vis.txt : text file with the month the analysis looked at

Visualization outputs:

- antigenic_score_map.png : Global map showing country scores, these have been calculated based on a weighted antigenic score based on each lineages frequency. For instance, in CountryA there is lineageA with an averaged score of 2.0 / frequency of 0.5 and lineageB with an averaged score of 1.5 / frequency of 0.5. The country score would then be equal to 2.0(0.5) + 1.5(0.5).\
- antigenic_score_map_cumulative.html : plotly output with a cumulative collection of global maps from Jan. 2020 to the current month of analysis
- antigenic_score_map_europe : Same concept as the global map but with a focus on European countries only
- antigenic_scoring_summary_heatmap.html : Heatmap with the top lineages, their frequencies and calculated z-scores
- antigenic_scoring_summary_states_heatmap.html : d3heatmap of the german states and their circulating lineages and respective antigenic scores


## Run pipeline on test data

There is provided test data under the test_run directory to test the pipeline. It contains all the necessary input for dummy data. Should you wish to recreate the results from the manuscript please download the GISAID metadata for the given sequence IDs (in the data/ dir of this repo) and then run the pipeline with the commands shown above. 

Otherwise to run the test please follow the following:

1. Activate the working environment

```console
conda activate aa_scoring
```

2. Navigate to the Corona_Variant_Scoring head directory (/Corona_Variant_Scoring/) - wherever this was installed on your machine
   
```console
cd /Corona_Variant_Scoring/
mkdir -p ../aa_scoring_out
```

3. Run the Corona_Variant_Pipeline with the following command. Keep in mind that if the output directory is same as the test_run directory it will overwrite the current output directory which contains the expected results. The pipeline is set to create an output/ directory with the results wherever the user points to as the output dir.

```console
bash variant_scoring.sh -o ../aa_scoring_out \
  -i ./test_run/ \
  -v ./ \
  -f ./test_run/SDplots_frequencies/ \
  -m ./test_run/09-2024_month.txt \
  -u ./test_run/metadata_under_review.tsv 09 2024
```
**Please note that here the -v argument is the current directory, if your current working directory in not in the Corona_Variant_Scoring please change this -v argument**

4. Check results - they should contain the same as what is currently in the test_run/output/ directory with all of the files outlined above in the analysis output section

Running this test pipeline took approximately 1 minute 10 seconds on a Macbook Air M2 16gb of ram. Times may vary depending on computer. 

### Table of Contents:

| Head Directory                  | Sub Directories                       | About                                                                                                                                                                                                                                                                                                                                                                  |
|---------------------------------|---------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Corona_Variant_Scoring/software/ |                                       |                                                                                                                                                                                                                                                                                                                                                                        |
|                                 | variant_scoring.py                    | python script to assign mutation scores to the different pango lineages occuring in each country. These mutation scores are based on the antigenic weights from the amino acid changes assigned for influenza. Outputs mutation score csv with european and global visualization (via plotly) of the scores                                                            |
| 				                            | variant_scoring_all_sites.py          | assigns antigenic scores to the different sequences dependent on ALL amino acid changes, not just amino acid changes occuring at known antigenic sites. Used to compare results with variant scoring analysis                                                                                                                                                          |
| 				                            | variant_scoring_without_weights.py    | assigns antigenic scores to lineages dependent on the number of mutations present at known antigenic sites, does not weight the score based on the up/down weight of the change. Used to compare results with the weighted and all sites results                                                                                                                       |
|                                 | VOC_comparison.py                     | outputs a visual comparison of known VOC's and their mutation scores in comparison to other pango lineages, used as a control comparison to make sure that we were seeing known VOCs with larger mutation scores. Also used to calculate the threshold for identifying variants of interest, based on the average mutation score of all the known variants of concern. |
|                                 | variant_scoring_aa_site_comparison.py | script to output a bar graph comparing amino acid site antigenic scores of VOC's (not including their sublineages)                                                                                                                                                                                                                                                     |
 |                                 | months_comparison_loop.sh             | used to run the frequency heatmap visualization on each month since the beginning of the pandemic (01/2020), creates its own month file and then requires the output of the antigenic scoring analysis (the antigenic_scores_all.csv files for each month from variant_scoring.py) as the input. Input .csv files must be in their own separate directory.             | 
|  								| time_comparison.py                    | runs the antigenic scoring analysis monthly for a given time frame and calculates each lineages' frequency for that month. Output for this was also used for the monthly visualizations from 01/2020 to 10/2022                                                                                                                                                        |
|                                   | LICENSE.md          | License for use of this software |
|                                   | sarscoverage.yml  | environment file for the pipeline (specifically for linux users) |
|                                   | variant_scoring.sh | shell script that runs all of the seperate analyses and visualization scripts    |
|                                   | sarscoverage_without_builds.yml  | environment file for the pipeline without builds, works with other operating systems |
|                                   | time_comparison_month_calc.py         | calculates the month based on the input files basename, used for the above loop                                                                                                                                                                                                                                                                                        |
| Corona_Variant_Scoring/reference/ |                                      |                                                                                                                                                                                                                                                                                                                                                                        |
|                                    | antigenic_weights.csv                 | csv file that contains the summation of the antigenic weights averaged across the different positions that they occured, from the influenza data (published here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3330098/), methodology used for this process can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4989879/                                                                                                                                             |
|                                    | known_variants_of_concern.csv | file that contains all of the variants of concern listed on the ecdc website (https://www.ecdc.europa.eu/en/covid-19/variants-concern) as well as descalated variants of concern as of 24 May 2022. They are listed as the pango lineage name on the ecdc website, and for this list the subvariants were also used (as listed on the cov-linead.org site) as the ecdc website states: "All sub-lineages of the listed lineages are also included in the variant, e.g., BA.2 is included in Omicron as it is a sub-lineage of B.1.1.529." (Ultimately included - Epsilon, Alpha (de-escalated variants), Beta, Gamma, Delta, & Omicron)|
|                                    | tp_sites.csv | file that contains curated list of known antigenic sites in the S1 subunit of the spike protein, from literature sources up through June 30, 2021 |
|                                    | amino_acid_properties.csv | contains the physical properties of each amino acid, copied from https://www.thermofisher.com/de/de/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/amino-acid-physical-properties.html, used for the prediction of antigenic weights, the molecular weights of the amino acids were retrieved from (https://www.thermofisher.com/de/de/home/references/ambion-tech-support/rna-tools-and-calculators/proteins-and-amino-acids.html)|
|                                    | antigenic_scores_map_visualization_cumulative.csv | cumulative file with all of the country antigenic scores from 2020 to October 2024  |
|                                    | antigenic_weights_higher_threshold.csv | contains a different set of weights, used when testing different methods, these weights were calculated from changes that occurred at least 3 times throughout the influenza A antigenic tree |
|                                    | country_list_with_threshold.tsv        | cumulative file of countries that have either 500 sequences in the month or represent at least 1% of the months global's sequences, this is used for the country_score_comparison_over_time analysis |
|                                    | months.txt          | example months text for input |
|                                    | weights_reversible.csv | reversible weights used in one of the tested methods, if weights not present for the reverse change then the antigenic weights for the original change would be used, for instance, A-D = 3 then D-A = 3 |
| Corona_Variant_Scoring/EVEscape_comparison/|                              |                                |
|                                    | EVEscape-01-2020_12-2023.csv      | contains the EVEscape results calculated using data from 01-2020 to 12-2023 |
|                                    | antigenicScores_wit_weights_at_all_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf | resultant scatterplot for the antigenic scoring method against both mFRN values and antigenic distances of select VOCs |
|                                    | evescape_comparison.R | script used to compare the antigenic scoring method and the EVEscape method |
|                                    | evescape_comparison.html | results of the script in html format |
|                                    | /evescape_scores_sigmoid_avg_comparison_mFRNA_antigenicCartography_45x180.pdf | resultant scatterplot for the EVEscape method against both mFRN values and antigenic distances of select VOCs |
| Corona_Variant_Scoring/SpikePro_comparison/ |                             |                                        |
|                                    | antigenicScores_with_weights_at_all_sites_comparison_mFRNA_antigenicCartography_2020-01-2023-12_45x180.pdf | resultant scatterplot for the antigenic scores vs the mFRN values and antigenic distances |
|                                    | spikepro_comparison.R  | analysis script for comparing antigenic scores and the SpikePro scores to mFRN values and antigenic distancess |


