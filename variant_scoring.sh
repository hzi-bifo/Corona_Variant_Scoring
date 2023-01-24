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
	while [ "$#" -gt 0 ]; do
	case $1 in
		-h | --h) 
			usage; exit
			;;	
		-o | --o)
			if [ -d "$2" ]; then
				OUTDIR=$2
			else
				echo "Please provide a valid output directory, see -h for additional information"; exit
			fi
			shift
			;;
#		-i | --i)
#			if [ -d "$2" ]; then
#				INDIR=$2
#			else
#				echo "Please provide a valid input directory, see -h for additional information"
#			fi
#			shift
#			;;
		-v | --v)
			if [ -d "$2" ]; then
				AntigenicScoring=$2
				SOFTWAREPATH=$AntigenicScoring"software/"
			else
				echo "Please provide a path to the Corona Variant Scoring directory, see -h for additional information"
			fi
			shift
			;;
		-f | --f)
			if [ -d "$2" ]; then
				FREQUENCY=$2
			else
				echo "Please provide a path to the frequency data, see -h for additional information"
			fi
			shift
			;;
		-m | --m)
			if [ -f "$2" ]; then
				MONTHS=$2
			else
				echo "Please provide a path to the months data, see -h for additional information"
			fi
			shift
			;;
		-u | --u)
			if [ -f "$2" ]; then
				SEQUI=$2
			else
				echo "Please provide a path to the sequences under review file, see -h for additional information"
			fi
			shift
			;;
		-q | --q) # Month input required for the variant_scoring.py script should, used for selecting specific months for analysis (ie months comparison)
			#if [ -f "$2" ]; then
			MONTH=$2
			#else
			#	echo "Please provide a month for the months comparison to work (in MM format), see -h for additional information"
			#fi
			shift
			;;
		-w | --w) # Year input require for the variant_scoring.py script
			#if [ -f "$2" ]; then
			YEAR=$2
			#else
			#	echo "Please provide a year for the months comparison to work (in YYYY format), see -h for additional information"
			#fi
			shift
			;;
		*)
			echo "$1 $2 is not an appropriate argument, please see the usage instructions: "; usage; exit
	esac
	shift
	done
fi
cd $OUTDIR

echo $OUTDIR
#echo $INDIR
echo $SOFTWAREPATH
echo $FREQUENCY
echo $MONTHS
echo $MONTH
echo $YEAR

#----------
# Analysis
#----------

# Creating an ouput directory if one has not been created already
if [ ! -d $OUTDIR'output/' ]; then mkdir $OUTDIR'output/'; fi

# Running Variant Scoring Analysis and Visualization
python $SOFTWAREPATH"variant_scoring.py" $INDIR'metadata.tsv' $AntigenicScoring"reference/tp_sites.csv" $OUTDIR'output/' $AntigenicScoring"reference/antigenic_weights.csv" $AntigenicScoring"reference/known_variants_of_concern.csv" $SEQUI # month (including 0 before value if < 10 (ex. 07 for july) # year
#python $SOFTWAREPATH"variant_scoring.py" $AntigenicScoring"reference/tp_sites.csv" $OUTDIR"output/" $AntigenicScoring"reference/antigenic_weights.csv" $AntigenicScoring"reference/known_variants_of_concern.csv" $SEQUI $MONTH $YEAR
#Rscript $SOFTWAREPATH"frequency_heatmap.R" $FREQUENCY $OUTDIR"output/antigenic_scores_ranked_with_WHO.csv" $MONTHS "0.1" $OUTDIR"output" 
Rscript $SOFTWAREPATH"frequency_heatmap.R" $OUTDIR"output/antigenic_scores_ranked_with_WHO.csv" $MONTHS "0.1" $OUTDIR"output/"
python $SOFTWAREPATH"global_scoring_map.py" $OUTDIR"output/antigenic_scores_map_visualization.csv" $OUTDIR"output/" $OUTDIR"output/month_vis.txt"
echo "COMPLETE"
