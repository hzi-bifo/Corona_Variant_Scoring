#!/bin/bash

OUTPUT=/net/viral_genomics/corona/katrina/Corona_Variant_Scoring/08_2022_results/month_comparison/
SOFTWAREPATH=/net/viral_genomics/corona/katrina/Corona_Variant_Scoring/software/
FREQUENCY=

cd $SOFTWAREPATH

# Creating a global map visual for all input files in the month_comparison directory (which are the output of the time_comparison.py)
echo Creating global visualizations time comparisons
#for file in $OUTPUT*.csv;
#do
#	python time_comparison_month_calc.py $file > $OUTPUT'/month.txt'
#	python global_scoring_map.py $file $OUTPUT $OUTPUT'/month.txt'
#done
echo COMPLETE

# Creating heatmap visuals for all input files in the month_comparison directory
echo Creating frequency heatmaps

cd $OUTPUT
echo "2020-01" > $OUTPUT"months.txt" # This should be the starting month of the analysis, or timeframe, you want to look at

ls *.csv | sort -n -t _ -k 4 | while read file;
do
	# Creating months file and output name
	echo $file
	name=$(basename $file | cut -c 31-37)
	echo $name
	#out_name=${string//[_]/-}
	month=${name:0:2}
	year=${name:3:7}

	if [ "$month" == "12" ]; then
		month="01"
		#year=year+1
		declare -i year
		year+=1
	else
		#month=month+1
		monthInt=${month#0}
		declare -i monthInt
		monthInt+=1
		if ((monthInt<10)); then
			month="0"$monthInt
		else
			month=$monthInt
		fi	
	fi
	echo "YEAR: "$year
	echo "MONTH: "$month
	echo $year"-"$month >> $OUTPUT"months.txt"
	
	# Making output directory if it does not already exist
#	if [[ ! -d $OUTPUT$name ]]; then 
#		mkdir $OUTPUT$name
#	fi
#	# Running heatmap visualization
#	Rscript $SOFTWAREPATH"frequency_heatmap.R" $FREQUENCY $file $OUTPUT"months.txt" "0.1" $OUTPUT$name 
done
echo COMPLETE
