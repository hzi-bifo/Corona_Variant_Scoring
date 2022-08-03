#!bin/bash

#----------
# Set Working Directory and Input
#----------

function usage () {
	echo
	echo 'Antigenic Scoring Pipeline'
	echo
	echo 'to run:'
        echo '$bash variant_scoring.sh -o <path to output directory file> -i <path to input directory> -v <path to Corona_Variant_Scoring repo> -f <path to frequencies data> -m <path to months text file>'
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
		-i | --i)
			if [ -d "$2" ]; then
				INDIR=$2
			else
				echo "Please provide a valid input directory, see -h for additional information"
			fi
			shift
			;;
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
			if [ -d "$2" ]; then
				MONTHS=$2
			else
				echo "Please provide a path to the months data, see -h for additional information"
			fi
			shift
			;;
		*)
			echo "$1 $2 is not an appropriate argument, please see the usage instructions: "; usage; exit
	esac
	shift
	done
fi
cd $OUTDIR

#----------
# Analysis
#----------

# Creating an ouput directory if one has not been created already
if [ ! -d $OUTDIR'output/' ]; then mkdir $OUTDIR'output/'; fi

# Running Variant Scoring Analysis and Visualization
python $SOFTWAREPATH"variant_scoring.py" $INDIR  $AntigenicScoring"reference/tp_sites.csv" $OUTDIR $AntigenicScoring"reference/antigenic_weights.csv" $AntigenicScoring"reference/known_variants_of_concern.csv"
Rscript $SOFTWAREPATH"frequency_heatmap.R" $FREQUENCY $OUTDIR"output/antigenic_scores_ranked_with_WHO.csv" $MONTHS "0.1"
python $SOFTWAREPATH"global_scoring_map.py" $OUTDIR"output/antigenic_scores_map_visualization.csv" $OUTDIR"output" $OUTDIR"month_vis.txt"
#echo "COMPLETE"
