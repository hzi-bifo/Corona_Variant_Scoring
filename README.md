## Corona Variant Prediction

Pipeline to assign an antigenic score to pangolin lineages based on the amino acid changes at known antigenic sites on the SARS-CoV-2 spike glycoprotein. The score is based on the summation of antigenic weights from the influenza dataset (available here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3330098/) that are then averaged across the different positions that they occured.

### Required Environment

The analysis pipeline is designed to run on BIFO servers, specifically along with the Sarscoverage environment that can be found here: https://github.com/hzi-bifo/coverage/tree/main/generate_web_data

It will require the path to the frequencies directory and months file, both of which are outputs from the SD Plots analysis.

### Running the pipeline:

To run the pipeline you will need to specify the output and input directory, path to the Corona_Variant_Scoring repo, path to frequencies data, and the paths to a months file (an example of this can be found in the reference/ directory). 

Run the following command to start the pipeline:

$bash variant_scoring.sh -o <path to output directory file> -i /path to input directory/ -v /path to Corona_Variant_Scoring repo/ -f /path to frequencies data/ -m /path to months text file/

Use bash variant_scoring.sh --help for more details.

### Analysis Output:

The pipeline will output two .csv files:
	- antigenic_scores_all.csv : file containing all sequences per country and their associated antigenic scores, raw file in which scores have not been averaged per pango lineage
	- antigenic_scores_ranked_with_WHO.csv : file containing the averaged antigenic scores across all sequences per Pango lineage, these scores are than ranked to show the lineages with the highest antigenic scores

Visualization outputs:
	- antigenic_score_map.png : Global map showing country scores, these have been calculated based on a weighted antigenic score based on each lineages frequency. For instance, in CountryA there is lineageA with an averaged score of 2.0 / frequency of 0.5 and lineageB with an averaged score of 1.5 / frequency of 0.5. The country score would then be equal to 2.0(0.5) + 1.5(0.5). 
	- antigenic_score_map_europe : Same concept as the global map but with a focus on European countries only.

For more information about each of the individual files please see the table_of_contents.md file 
