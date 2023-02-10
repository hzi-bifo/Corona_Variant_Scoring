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
	while getopts 'o:i:v:f:m:u:hqw' opt; do
		case "$opt" in
			h) usage; exit ;;
			o)
				if [ -d "$OPTARG" ]; then
					OUTDIR=$OPTARG
				else
					echo "Please provide a valid output directory, see -h for additional information"; exit
				fi
				;;
			i)
				if [ -d "$OPTARG" ]; then
					INDIR=$OPTARG
				else
					echo "Please provide a valid input directory, see -h for additional information"; exit
				fi
				;;
			v)
				if [ -d "$OPTARG" ]; then
					AntigenicScoring=$OPTARG
					SOFTWAREPATH=$AntigenicScoring"software/"
				else
					echo "Please provide a path to the Corona Variant Scoring directory, see -h for additional information"; exit
				fi
				;;
			f)
				if [ -d "$OPTARG" ]; then
					FREQUENCY=$OPTARG
				else
					echo "Please provide a path to the frequency data, see -h for additional information"; exit
				fi
				;;
			m)
				if [ -f "$OPTARG" ]; then
					MONTHS=$OPTARG
				else
					echo "Please provide a path to the months data, see -h for additional information"; exit
				fi
				;;
			u)
				if [ -f "$OPTARG" ]; then
					SEQUI=$OPTARG
				else
					echo "Please provide a path to the sequences under review file, see -h for additional information"; exit
				fi
				;;
			q) # Month input required for the variant_scoring.py script should, used for selecting specific months for analysis (ie months comparison)
				MONTH=$OPTARG
				;;
			w) # Year input require for the variant_scoring.py script
				YEAR=$OPTARG
				;;
			*)
				echo "$1 $2 is not an appropriate argument, please see the usage instructions: "; usage; exit
		esac
	#shift
	done
fi
cd $OUTDIR

#----------
# Analysis
#----------

# Creating an ouput directory if one has not been created already
if [ ! -d $OUTDIR'output/' ]; then mkdir $OUTDIR'output/'; fi

eval "$(conda shell.bash hook)"
conda activate sarscoverage

# Running Variant Scoring Analysis and Visualization
python "$SOFTWAREPATH""variant_scoring.py" "$INDIR""metadata.tsv" "$AntigenicScoring""reference/tp_sites.csv" "$OUTDIR""output/" "$AntigenicScoring""reference/antigenic_weights.csv" "$AntigenicScoring""reference/known_variants_of_concern.csv" "$SEQUI" > "$OUTDIR""STDOUT.txt" # month (including 0 before value if < 10 (ex. 07 for july) # year
#python $SOFTWAREPATH"variant_scoring.py" $AntigenicScoring"reference/tp_sites.csv" $OUTDIR"output/" $AntigenicScoring"reference/antigenic_weights.csv" $AntigenicScoring"reference/known_variants_of_concern.csv" $SEQUI $MONTH $YEAR
Rscript "$SOFTWAREPATH""frequency_heatmap.R" "$FREQUENCY" "$OUTDIR""output/antigenic_scores_ranked_with_WHO.csv" "$MONTHS" "0.1" "$OUTDIR""output/" 
python "$SOFTWAREPATH""global_scoring_map.py" "$OUTDIR""output/antigenic_scores_map_visualization.csv" "$OUTDIR""output/" "$OUTDIR""output/month_vis.txt"
echo "COMPLETE"
