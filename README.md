## Corona Variant Prediction

Pipeline to assign an antigenic score to pangolin lineages based on the amino acid changes at known antigenic sites on the SARS-CoV-2 spike glycoprotein. The score is based on the summation of antigenic weights from the influenza dataset (available here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3330098/) that are then averaged across the different positions that they occured.

### Required Environment

The analysis pipeline is designed to run on BIFO servers, specifically along with the Sarscoverage environment that can be found here: https://github.com/hzi-bifo/coverage/tree/main/generate_web_data

It will require the path to the frequencies directory and months file, both of which are outputs from the SD Plots analysis.

#### Required Software & Version:

- python
- R
- pandas
- numpy
- re
- itertools
- plotly

### Running the pipeline:

To run the pipeline you will need to specify the output and input directory, path to the Corona_Variant_Scoring repo, path to frequencies data, and the paths to a months file (an example of this can be found in the reference/ directory). 

Run the following command to start the pipeline:

```
$bash variant_scoring.sh -o /path to output directory file/ -i /path to input directory/ 
-v /path to Corona_Variant_Scoring repo/ -f /path to frequencies data/ -m /path to months text file/ 
-u /path to sequences under review file/
```

The user also has the option of using -q and -w to add a specified month and year respectively for analysis. Each input
file is outlined below:
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

**Use bash variant_scoring.sh --help for more details.**

### Analysis Output:

The pipeline will output two .csv files:
- antigenic_scores_all.csv : file containing all sequences per country and their associated antigenic scores, raw file in which scores have not been averaged per pango lineage\
- antigenic_scores_ranked_with_WHO.csv : file containing the averaged antigenic scores across all sequences per Pango lineage, these scores are than ranked to show the lineages with the highest antigenic scores

Visualization outputs:
- antigenic_score_map.png : Global map showing country scores, these have been calculated based on a weighted antigenic score based on each lineages frequency. For instance, in CountryA there is lineageA with an averaged score of 2.0 / frequency of 0.5 and lineageB with an averaged score of 1.5 / frequency of 0.5. The country score would then be equal to 2.0(0.5) + 1.5(0.5).\
- igenic_score_map_europe : Same concept as the global map but with a focus on European countries only.

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
|                                   | time_comparison_month_calc.py         | calculates the month based on the input files basename, used for the above loop                                                                                                                                                                                                                                                                                        |
| Corona_Variant_Scoring/reference/ |                                      |                                                                                                                                                                                                                                                                                                                                                                        |
|                                    | antigenic_weights.csv                 | csv file that contains the summation of the antigenic weights averaged across the different positions that they occured, from the influenza data (published here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3330098/), methodology used for this process can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4989879/                                                                                                                                             |
|                                    | known_variants_of_concern.csv | file that contains all of the variants of concern listed on the ecdc website (https://www.ecdc.europa.eu/en/covid-19/variants-concern) as well as descalated variants of concern as of 24 May 2022. They are listed as the pango lineage name on the ecdc website, and for this list the subvariants were also used (as listed on the cov-linead.org site) as the ecdc website states: "All sub-lineages of the listed lineages are also included in the variant, e.g., BA.2 is included in Omicron as it is a sub-lineage of B.1.1.529." (Ultimately included - Epsilon, Alpha (de-escalated variants), Beta, Gamma, Delta, & Omicron)|
|                                    | tp_sites.csv | ile that contains curated list of known antigenic sites in the S1 subunit of the spike protein |
|                                    | amino_acid_properties.csv | contains the physical properties of each amino acid, copied from https://www.thermofisher.com/de/de/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/amino-acid-physical-properties.html, used for the prediction of antigenic weights|

