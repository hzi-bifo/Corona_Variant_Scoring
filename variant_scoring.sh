#!/bin/bash 

#----------
# Set Working Directory and Input
#----------

function usage () {
	echo
	echo 'Antigenic Scoring Pipeline'
	echo
	echo 'to run:'
        echo '$bash variant_scoring.sh -o <path to output directory file> -i <path to input directory> -v <path to Corona_Variant_Scoring repo> -f <path to frequencies data> -m <path to months text file -u <path to the sequences under review file>'
        echo
	echo '-o / --o : desired directory for the output files'
	echo '-i / --i : desired input file, this is the GISAID metadata file that includes global patient sequences and amino acid changes in the spike protein for those sequences'
	echo '-v / --v : path to Corona_Variant_Scoring directory (../Corona_Variant_Scoring/)'
	echo '-f / --f : path to frequencies data (output of SD plot analysis)'
	echo '-m / --m : path to month text file which contains a list of months used to determine current month for the heatmap visualization'
	echo '-u / --u : path to the under investigation sequences csv, this is a list of accession ids that are currently under review from GISAID, they will be removed prior to analysis'
	echo
}

if [ "$#" == 0 ]; then
	echo
	echo "Please provide the path to a valid In/Out Directory. For additional help use -h / --h" 
	echo
	exit
else
	while getopts 'o:i:v:f:m:u:hqw' OPTION; do
		case "${OPTION}" in
			h) usage; exit ;;
			o)
				if [ -d "${OPTARG}" ]; then
					OUTDIR=${OPTARG}
				else
					echo "Please provide a valid output directory, see -h for additional information"; exit
				fi
				;;
			i)
				if [ -d "${OPTARG}" ]; then
					INDIR=${OPTARG}
				else
					echo "Please provide a valid input directory, see -h for additional information"; exit
				fi
				;;
			v)
				if [ -d "${OPTARG}" ]; then
					AntigenicScoring=${OPTARG}
					SOFTWAREPATH=$AntigenicScoring"software/"
				else
					echo "Please provide a path to the Corona Variant Scoring directory, see -h for additional information"; exit
				fi
				;;
			f)
				if [ -d "${OPTARG}" ]; then
					FREQUENCY=${OPTARG}
				else
					echo "Please provide a path to the frequency data, see -h for additional information"; exit
				fi
				;;
			m)
				if [ -f "${OPTARG}" ]; then
					MONTHS=${OPTARG}
				else
					echo "Please provide a path to the months data, see -h for additional information"; exit
				fi
				;;
			u)
				if [ -f "${OPTARG}" ]; then
					SEQUI=${OPTARG}
				else
					echo "Please provide a path to the sequences under review file, see -h for additional information"; exit
				fi
				;;
			q) # Month input required for the variant_scoring.py script should, used for selecting specific months for analysis (ie months comparison)
				MONTH=${OPTARG}
				;;
			w) # Year input require for the variant_scoring.py script
				YEAR=${OPTARG}
				;;
			?)
				echo "$1 $2 is not an appropriate argument, please see the usage instructions: "; usage; exit
		esac
	#shift
	done
fi
cd "$OUTDIR"

#----------
# Analysis
#----------

# Creating an ouput directory if one has not been created already
if [ ! -d "$OUTDIR"'output/' ]; then mkdir "$OUTDIR"'output/'; fi

eval "$(conda shell.bash hook)"
conda activate sarscoverage

# Running Variant Scoring Analysis and Visualization
echo "Running Variant Scoring Analysis and Visualization"
# Antigenic Scoring Analysis
python "$SOFTWAREPATH""variant_scoring_all_sites.py" "$INDIR""metadata.tsv" "$AntigenicScoring""reference/tp_sites.csv" "$OUTDIR""output/" "$AntigenicScoring""reference/antigenic_weights.csv" "$AntigenicScoring""reference/known_variants_of_concern.csv" "$SEQUI" > "$OUTDIR""STDOUT.txt" # month (including 0 before value if < 10 (ex. 07 for july) # year

# Frequency Heatmap
echo "Creating Frequency Heatmap"
if [ -f "$OUTDIR"'output/month_remove.txt' ]; then grep -Fvxf "$OUTDIR"'output/month_remove.txt' "$MONTHS" > "$OUTDIR"'output/month_corrected.txt'; MONTHS="$OUTDIR"'output/month_corrected.txt'; fi
#Rscript "$SOFTWAREPATH""frequency_heatmap.R" "$FREQUENCY" "$OUTDIR""output/antigenic_scores_ranked_with_WHO.csv" "$MONTHS" "0.1" "$OUTDIR""output/" "$OUTDIR""output/antigenic_scores_all.csv"
Rscript "$SOFTWAREPATH""frequency_heatmap_coverage.R" "$FREQUENCY" "$OUTDIR""output/antigenic_scores_ranked_with_WHO.csv" "$MONTHS" "0.1" "$OUTDIR""output/" "$OUTDIR""output/antigenic_scores_all.csv" >> "$OUTDIR""STDOUT.txt"

# Global Map of Antigenic Scores
echo "Creating Global Map"
python "$SOFTWAREPATH""global_scoring_map.py" "$OUTDIR""output/antigenic_scores_map_visualization.csv" "$OUTDIR""output/" "$OUTDIR""output/month_vis.txt" "$AntigenicScoring""reference/" >> "$OUTDIR""STDOUT.txt"

# Country-wise Plot
#echo "Creating Country-Wise Plot"
python "$SOFTWAREPATH""country_frequency_threshold_compiler.py" "$OUTDIR""output/" "$OUTDIR""output/" "$AntigenicScoring""reference/" "$OUTDIR""output/month_vis.txt" >> "$OUTDIR""STDOUT.txt"
Rscript "$SOFTWAREPATH""country_score_over_time_coverage.R" "$OUTDIR""output/" "$AntigenicScoring""reference/antigenic_scores_map_visualization_cumulative.csv" "$AntigenicScoring""reference/country_list_with_threshold.tsv" >> "$OUTDIR""STDOUT.txt"

# Selected pVOI table
#python "$SOFTWAREPATH""pVOI_interactive_table.py" "$OUTDIR""output/antigenic_scoring_summary_pVOI_table.csv" "$OUTDIR""output/"
echo "COMPLETE"
