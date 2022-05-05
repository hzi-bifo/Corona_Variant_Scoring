### Table of Contents

#### Software
Corona_Variant_Scoring/software/ - contains all the necessary analysis software used in the pipeline
	- ../software/antigenic_weights_reference.ipynb : converts the antigenic_weights_reference_original.csv to the antigenic_weights_reference.csv file, it averages the weights across the different changes and outputs the finalized antigenic weight reference file
	- ../software/variant_scoring.py : python script to assign mutation scores to the different pango lineages occuring in each country. These mutation scores are based on the antigenic weights from the amino acid changes assigned for influenza. Outputs mutation score csv with european and global visualization (via plotly) of the scores
	- ../software/variant_scoring_all_sites.py : assigns antigenic scores to the different sequences dependent on ALL amino acid changes, not just amino acid changes occuring at known antigenic sites. Used to compare results with variant scoring analysis
	- ../software/variant_scoring_without_weights.py : assigns antigenic scores to lineages dependent on the number of mutations present at known antigenic sites, does not weight the score based on the up/down weight of the change. Used to compare results with the weighted and all sites results
	- ../software/VOC_comparison.py : outputs a visual comparison of known VOC's and their mutation scores in comparison to other pango lineages, used as a control comparison to make sure that we were seeing known VOCs with larger mutation scores. Also used to calculate the threshold for identifying variants of interest, based on the average mutation score of all the known variants of concern. 
	- ../software/variant_scoring_aa_site_comparison.py : script to output a bar graph comparing amino acid site antigenic scores of VOC's (not including their sublineages) 
	- ../software/voc_box_plot_visualization.py : outputs a box plot comparing the antigenic scores of the variants of concern (omicron, gamma, beta, alpha, delta, episilon) to other non-variant of concern lineages

#### Reference files
Corona_Variant_Scoring/reference - contains reference files for assigning mutation scores (based on the antigenic weights for amino acid changes)
	- ../reference/antigenic_weights.csv : .csv file that contains the summation of the antigenic weights averaged across the different positions that they occured, from the influenza data (published here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3330098/)
		- methodology used for this process can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4989879/
	- ../reference/known_variants_of_concern.csv : file that contains all of the variants of concern listed on the ecdc website (https://www.ecdc.europa.eu/en/covid-19/variants-concern) as well as descalated variants of concern as of 24 May 2022. They are listed as the pango lineage name on the ecdc website, and for this list the subvariants were also used (as listed on the cov-linead.org site) as the ecdc website states: "All sub-lineages of the listed lineages are also included in the variant, e.g., BA.2 is included in Omicron as it is a sub-lineage of B.1.1.529." (Ultimately included - Epsilon, Alpha (de-escalated variants), Beta, Gamma, Delta, & Omicron)
	- ../reference/tp_sites.csv : file that contains curated list of known antigenic sites in the S1 subunit of the spike protein

#### Output
Corona_Variant_Scoring/test/ - contains an initial test run of the pipeline as an example of the required inputs and outputs
	- ../test/metadata.tsv : initial input for pipeline, obtained from GISAID EpiCoV metadata on 11.02.2022
	- ../test/voc_threshold.txt : text file containing the averaged mutation score of the current variants of concern (the threshold for finding other variants of interest) and the date to which it was run
Corona_Variant_Scoring/test/output - output file for the test run on the given metadata.tsv file, contains the different output for each of the methodology tests (later used to compare and identify the optimal method), base file contains the test results for the methodology weighing known TP sites based on the influenza phylogenetic tree and then assigning antigenic scores based on those values
	- ../output/antigenic_scores_all.csv : antigenic scores PER SEQUENCE, not averaged per location or lineage
	- ../output/antigenic_scores_ranked_with_WHO.csv : averaged antigenic scores per sequence across pango lineages and then ranked based on antigenic score. The WHO label assigned is based on the WHO label given to each known variant of concern and its sublineages
	- ../output/mutation_score_VOC_box_plot.png : box plot that compares the distribution of antigenic scores between the known variants of concern (omicron, gamma, beta, alpha, delta, and episilon) and the other pango lineages, used to compare methodology for antigenic scoring and confirm results
	- ../output/variant_scoring_without_weights/ :
	
	
	- ../test/mutation_score_VOC_box_plot.jpeg : boxplot of the averaged mutation scores for each pango lineage. Broken down into known variants of concern and their sublineages (according to the reference file) and then other non-variants of concern lineages to compare the antigenic scores between the two
	- ../test/mutation_score_VOC_density_plot.jpeg : density plot of the averaged mutation scores for each pango lineage. Broken down into known variants of concern and their sublineages (according to the reference file) and then other non-variants of concern lineages to compare the antigenic scores between the two
	
Corona_Variant_Scoring/test/output - output for the test analysis
	- ../output/brazil_mutation_scores.csv : mutation scores for only Brazil, which was noted to have a higher average mutation score across all lineages. 
	- ../output/geramany_mutation_scores.csv : mutation scores for only Germany, for a closer look at regions across the country
	- ../output/mutation_scores.csv : averaged mutation scores across each lineage occurring in each unique location
	- ../output/mutation_scores_all.csv : mutation score per SEQUENCE, not averaged per location. Mutation score assigned is based on the antigenic weight of the amino acid changes occurring at known antigenic sites (TP sites)
	- ../output/mutation_scores_ranked_with_WHO.csv : averaged mutation scores per sequence across pango lineages and then ranked based on mutation score. The WHO label assigned is based on the WHO label given to each known variant of concern and its sublineages
	- ../output/mutation_score_map.html/png : global map of averaged mutation score across all lineages occurring per country
	- ../output/mutation_score_map_europe.html/png : european map of averaged mutation score across all lineages occurring per country
