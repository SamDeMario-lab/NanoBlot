#!/bin/bash

#Default Values 
PLOTS="./user_input_files/plot_data.tsv" 
PROBES="./user_input_files/probes.bed"
META_DATA="./user_input_files/data_metadata.tsv"
NANO_BLOT_RSCRIPT="./scripts/nano_blot_generation.R"
ANNOTATION_FILE="./user_input_files/Saccharomyces_cerevisiae.R64-1-1.107.gtf"
PRINT_HELP=FALSE
SUBSET_BAMS=TRUE
CDNA=FALSE
NORM=TRUE
RT_PCR=FALSE
CLEAN_ALL=FALSE
OVERWRITE_COUNTS=FALSE
NORM_FACTOR=0 
# Values of NORM_FACTOR, 0=differential, 1=size, 2=skip 

echo -e "Starting Nanoblot" '\u2622' 
echo "======="

while getopts ":HFCWOR:M:B:T:Y:A:N:" opt; do
	case $opt in
		H ) 
		PRINT_HELP=TRUE
		echo "Help printed"
			;;
		F ) 
		SUBSET_BAMS=FALSE
		echo "Subsetting BAM files skipped"
			;;
		C ) 
		CDNA=TRUE
		echo "Treating reads as cDNA"
			;;
		R ) 
		NANO_BLOT_RSCRIPT=$OPTARG
		echo "Using custom R script $NANO_BLOT_RSCRIPT"
			;;
		M ) 
		META_DATA=$OPTARG
		echo "Using location of metadata file $METADATA"
			;;
		B ) 
		PLOTS=$OPTARG
		echo "Using blots metadata file $PLOTS"
			;;
		T ) 
		PROBES=$OPTARG
		echo "Using probes bed file $PROBES"
			;;
		N )
		if [ $OPTARG = "differential" ]; then
			NORM_FACTOR=0
			echo "Normalization method: Differential"
		elif [ $OPTARG = "size" ]; then
			NORM_FACTOR=1
			echo "Normalization method: Size"
		elif [ $OPTARG = "skip" ]; then
			NORM=FALSE
			NORM_FACTOR=2
			echo "Skipping data normalization"
		else
			echo "Invalid option for -N: please choose from {differential, size, skip}"
			exit 1
		fi
			;;
		A ) 
		ANNOTATION_FILE=$OPTARG
		echo "Using annotation file $ANNOTATION_FILE"
			;;
		O ) 
		OVERWRITE_COUNTS=TRUE
		echo "Overwriting counts generation"
		COUNT_FOLDER="./temp/count_tables"
		if [ ! -d "$COUNT_FOLDER" ]
		then
			echo "Invalid option for -O, counts tables folder needs to exist first"
			exit 1
		fi
		rm -r $COUNT_FOLDER
			;;
		Y )
		RT_PCR=TRUE
		echo "Running Nano RT-PCR"
		PLOTS=$OPTARG
			;;
		W )
		CLEAN_ALL=TRUE
		echo "Clear all files from temp after plot generation"
			;;
		\? ) echo "Invalid option: -$OPTARG"
		PRINT_HELP=TRUE
			;;
		: )
		if [ $OPTARG = "Y" ]
		then
			RT_PCR=TRUE
			echo "Running default Nano RT-PCR script"
			PLOTS="./user_input_files/plot_data_RTPCR.tsv"
		else
    	echo "Invalid Option: -$OPTARG requires an argument" 1>&2
    	exit 1
    fi
      ;;
	esac
done

#Check if user asked for the help text and echo the help text.
if [ $PRINT_HELP = TRUE ]
then
	echo "
Nanoblot (Version 1.0)

O       o O       o O       o O       o
| O   o | | O   o | | O   o | | O   o |
| | O | | | | O | | | | O | | | | O | |
| o   O | | o   O | | o   O | | o   O |
o       O o       O o       O o       O
             xxxxxxx
         x xxxxxxxxxxxxx x
      x     xxxxxxxxxxx     x
            xxxxxxxxx
  x          xxxxxxx          x
              xxxxx
 x             xxx             x
                x
xxxxxxxxxxxxxxx   xxxxxxxxxxxxxxx
 xxxxxxxxxxxxx     xxxxxxxxxxxxx
  xxxxxxxxxxx       xxxxxxxxxxx
   xxxxxxxxx         xxxxxxxxx
     xxxxxx           xxxxxx
       xxx             xxx
            x         x
                 x
O       o O       o O       o O       o
| O   o | | O   o | | O   o | | O   o |
| | O | | | | O | | | | O | | | | O | |
| o   O | | o   O | | o   O | | o   O |
o       O o       O o       O o       O

Current Dependencies
==========================================
Package Name   | Version 
Bedtools       | 2.30.0
Samtools       | 1.15.1
HTSeq          | 2.0.2
R              | 4.2.1
R (ggplot2)    | 3.3.6
R (samtools)   | 3.15
R (ggridges)   | 0.5.3
R (DeSeq2)     | 3.15
R (dplyr)      | 1.0.10

For an explanation of the required input files see the README.md
==========================================
-H  |  Print help menu
-T  |  Probes bed file 
-B  |  Plotting file
-M  |  Location of metadata file
-A  |  Annotation file 
-R  |  Use custom R script
-Y  |  RT-PCR mode, supply with own metadata file 
-N  |  Normalization function {differential (default), size, skip}
-O  |  Overwrite count tables and recalculate
-C  |  Treat reads as cDNA (disregard strand) 
-F  |  Skip subsetting BAM files for plot generation
-W  |  Clear all files from ./temp/ after plot generation
==========================================
"
exit
fi

#Checks to see if there are duplicate probe names
if [ $(awk '{print $4}' $PROBES | sort | uniq -d | wc -l) -ne 0 ]
then
	echo "Duplicate probes found in ${PROBES}, please fix and rerun"
	exit
fi

#Checks to see if there are duplicate sample names
if [ $(awk '{print $1}' $META_DATA | sort | uniq -d | wc -l) -ne 0 ]
then
	echo "Duplicate metadata samples found in ${META_DATA}, please fix and rerun"
	exit
fi

sleep 1 #Pauses for one second to actually let the user know that the program is being started
echo "R Script: $NANO_BLOT_RSCRIPT";
echo "Meta Data File: $META_DATA";

declare -i END_META=$(grep -c '.' $META_DATA) 

echo "Probes Bed File: $PROBES";
echo "Plots File: $PLOTS";

#This fixes files that have bene generated in windows computers 
TEMP_PLOT="./temp/temp_plot.tsv"
TEMP_META="./temp/temp_metadata.tsv"
TEMP_PROBES="./temp/temp_probes.bed"
cp $PLOTS $TEMP_PLOT
cp $META_DATA $TEMP_META 
cp $PROBES $TEMP_PROBES
awk 'sub("\r$", "")+1' $TEMP_PLOT > $PLOTS
awk 'sub("\r$", "")+1' $TEMP_META > $META_DATA
awk 'sub("\r$", "")+1' $TEMP_PROBES > $PROBES
rm $TEMP_PLOT $TEMP_META $TEMP_PROBES 

declare -i END_PLOT=$(grep -c '.' $PLOTS) 

# Checking to make sure all input files in the correct format
ALLFILES=($PLOTS $PROBES $META_DATA)
for file in $ALLFILES
do
	NUM_SPACES=$(grep -o ' ' $file | wc -l)
	if [ $NUM_SPACES -gt 0 ]
	then
		echo "Incorrect file format for $file. There are spaces. Please fix and then rerun"
		exit 1
	fi
done

if [ $NORM = TRUE ] 
then
	for (( c=2; c<=$END_PLOT; c++ ))
	do
		P_LINE=$(head -n $c $PLOTS | tail -n -1) #This line gets the individual row for each row of the plot_data
		IFS=$'\t'; read -a fields <<<"$P_LINE" 
		
		# Check to see if first character in the first column of the plot_csv file is a #, which if it is, 
		# will skip that line's normalization in addition to plot generation
		if [ ${fields[0]::1} = "#" ]
		then
			echo "Skipping ${fields[0]} normalization"
			continue
		fi
		BAMS=${fields[1]} 
		
		if ! [[ $BAMS =~ , ]]
		then
			NORM=FALSE
			NORM_FACTOR=2
			echo "Skipping data normalization for ${fields[0]} due to single sample"
			continue
		fi
		
		# This if is for the DeSeq2 normalization method which first generates count tables
		if [ $NORM_FACTOR = 0 ]; then
			echo -e "=======\nNormalization: Generating Count Tables for ${fields[0]}"
			echo Annotation File Location: $ANNOTATION_FILE
			
			COUNT_FOLDER="./temp/count_tables"
			if [ ! -d "$COUNT_FOLDER" ] # checks if the norm folder exists and is a directory 
			then
				mkdir -p $COUNT_FOLDER
			fi
			
			IFS=',';
			for sample in $BAMS
				do
					# Creates count tables for each of the samples based off of their metadata location file
					# Stores those count tables into a NORM folder
					# This norm folder will then be later accessed in the nano_blot_generation which will call the normalization.R script
					COUNT_FILE_NAME="${COUNT_FOLDER}/${sample}-htseq_counts.tsv"
					# First check if htseq-count is necessary, if not, then print that it was already counted
					if [ ! -f $COUNT_FILE_NAME ] 
					then
						echo "Running htseq-count for $sample"
						DATA_LINE_T=$(awk -v var="$sample" '$1==var {print $0}' $META_DATA)
						IFS=$'\t' read -a EELS <<<"$DATA_LINE_T" 
						DATA_LOCATION=${EELS[1]}
						python3 -m HTSeq.scripts.count -m union $DATA_LOCATION $ANNOTATION_FILE> $COUNT_FILE_NAME
					fi
				done
		# This if is for the size normalization method, essentially CPM 
		elif [ $NORM_FACTOR = 1 ]; then
			echo -e "=======\nNormalization: Counts Per Million for ${fields[0]}"
			IFS=',';
			for sample in $BAMS
				do
				DATA_LINE_T=$(awk -v var="$sample" '$1==var {print $0}' $META_DATA)
				IFS=$'\t' read -a EELS <<<"$DATA_LINE_T" 
				DATA_LOCATION=${EELS[1]}
				echo $sample: $(samtools view -c -F 260 $DATA_LOCATION) reads 
			done
		
		# This only runs if there is an error with setting the norm factor, this should never run 
		else
			echo Error with NORM_FACTOR, value is $NORM_FACTOR
			exit 1
		fi
	done # finishes the first loop
fi #finishes the if to check if NORM is true 

declare -i PLOT_NUM=1
# Probe data + configuring plots 
for (( j=2; j<=$END_PLOT; j++ )) # For loop that goes through each row of plot_data
do
	P_LINE=$(head -n $j $PLOTS | tail -n -1) #Gets each individual row 
	echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="
	IFS=$'\t'; read -a fields <<<"$P_LINE"
	
	# Check to see if first character in the first column of the plot_csv file is a #, which if it is, 
		# will skip that line's plot generation
		if [ ${fields[0]::1} = "#" ]
		then
			echo "Skipping ${fields[0]} plotting"
			continue
		fi
		
	echo "Generating plot" $PLOT_NUM
	((PLOT_NUM++))
	echo "======="
	
	if [ $RT_PCR = FALSE ]
	then
		PROBE_FIELD=${fields[2]}
		ANTIPROBE_FIELD=${fields[3]}
	else 
		VIEWING_WINDOW=${fields[2]}
		PROBE_FIELD=${fields[3]}
		ANTIPROBE_FIELD=${fields[4]}
	fi
	
	echo 'Probe(s):'$PROBE_FIELD
	if [ -z "$ANTIPROBE_FIELD" ] #checks to see if there is an antiprobe, -z checks for empty string 
	then
		echo "No Negative Probe Used" #Print statement
	else
		echo 'Negative Probe(s):'$ANTIPROBE_FIELD #Print statement
	fi
	echo "======="

	# Setting variables 
	# TARGET is an array variable 
	IFS=',' read -a TARGET <<< "$PROBE_FIELD"
	
	# TARG_NEGS is an array variable 
	IFS=',' read -a TARG_NEGS <<< "$ANTIPROBE_FIELD"
	
	# BAMS is a 
	BAMS=${fields[1]}
	echo $BAMS

	PREVIOUS_PROBE="" #Needs a tracker to see if the previous probe was already subsetted

	# Go through each target probe
	for probe in "${TARGET[@]}"
		do
			PR_LINE=$(awk -v var="$probe" '$4==var {print $0}' $PROBES)
			IFS=$'\t'; read -a feels <<< "$PR_LINE"

			if [[ -z "${feels[@]}" ]]
			then
				echo "Probe: $probe not found. Check bed file and blots metadata file. Exiting script"
				exit 1
			fi
			echo 'Probe'
			echo $probe": "${feels[0]}:${feels[1]}-${feels[2]} 
			echo "$PR_LINE" > "./temp/temp_bed.bed" 

			IFS=$',';
			for sample in $BAMS;
			do
				DATA_LINE_T=$(awk -v var="$sample" '$1==var {print $0}' $META_DATA)
				IFS=$'\t' read -a EELS <<<"$DATA_LINE_T" 
				SAMPLE_NAME=${EELS[0]} #Need to filter individual data line to get the sample name
				
				if [[ -z "DATA_LINE_T" ]]
				then
					echo "Metadata: $sample not found in $META_DATA. Check meta_data file. Exiting script"
					exit 1
				fi
				
				if ! [ -z $PREVIOUS_PROBE ] #Checks if its not an empty string
				then
					DATA_LOCATION="./temp/"$SAMPLE_NAME"_"$PREVIOUS_PROBE".bam"
					TEMP_NAME=${SAMPLE_NAME}_${PREVIOUS_PROBE}_${probe}.bam
				else
					DATA_LOCATION=${EELS[1]}
					TEMP_NAME=$SAMPLE_NAME"_"$probe".bam"
				fi
				
				if [[ "$SUBSET_BAMS"  == TRUE ]]
				then
				echo -e "Subsetting: "$DATA_LOCATION"\nNaming Subset:  "$TEMP_NAME
					
					if [[ "$CDNA"  == TRUE ]]
					then
						bedtools intersect -a $DATA_LOCATION -b "./temp/temp_bed.bed" -wa -split -nonamecheck > "./temp/$TEMP_NAME"
						samtools index "./temp/$TEMP_NAME"
					else
						bedtools intersect -a $DATA_LOCATION -b "./temp/temp_bed.bed" -wa -split -s -nonamecheck > "./temp/$TEMP_NAME"
						samtools index "./temp/$TEMP_NAME"
						# the -s forces overlap of B and A on the same strand, which is what we want most of the time
					fi
				fi
			done
			if [[ "$SUBSET_BAMS" == FALSE ]]
			then
				echo "Skipping subsetting target probe $probe BAM files. If subsetting is desired remove -F flag."
			fi
			
			if ! [ -z $PREVIOUS_PROBE ]
			then
				PREVIOUS_PROBE=${PREVIOUS_PROBE}_$probe #update previous probe right before the for loop for each probe updates 
			else
				PREVIOUS_PROBE=$probe
			fi
	done
	
	#Go through each anti-target probe 
	PREVIOUS_ANTI_PROBE=$PREVIOUS_PROBE

	if [ -z "$ANTIPROBE_FIELD" ]
	then
		echo "No negative probe used"
	else
		echo 'Negative Probe used'
	
		for antiprobe in "${TARG_NEGS[@]}"
			do
				APR_LINE=$(awk -v var="$antiprobe" '$4==var {print $0}' $PROBES)
				IFS=$'\t'; read -a deels <<< "$APR_LINE"
	
				if [[ -z "${deels[@]}" ]]
				then
					echo "Antiprobe: $antiprobe not found. Check bed file and blots metadata file. Exiting script"
					exit #Not sure if I want this to exit just yet
				fi
				echo 'Antiprobe'
				echo $antiprobe": "${deels[0]}:${deels[1]}-${deels[2]} 
				echo "$APR_LINE" > "./temp/temp_anti_bed.bed"
				
				IFS=$',';
				for sample in $BAMS; ##You need to only create subsets of the samples you are using 
				do
					DL_ANTI=$(awk -v var="$sample" '$1==var {print $0}' $META_DATA)
					IFS=$'\t' read -a bells <<<"$DL_ANTI" 
					SAMPLE_NAME=${bells[0]} #Neeed to filter individual data line to get the sample name
					DATA_LOCATION="./temp/"$SAMPLE_NAME"_"$PREVIOUS_ANTI_PROBE".bam"
					TEMP_NAME=${SAMPLE_NAME}_${PREVIOUS_ANTI_PROBE}_anti_${antiprobe}.bam
					
					if [[ "$SUBSET_BAMS"  == TRUE ]]
					then
						echo -e "Subsetting: "$DATA_LOCATION"\nNaming Subset:  "$TEMP_NAME
						if [[ "$CDNA"  == TRUE ]]
							then
								bedtools intersect -a $DATA_LOCATION -b "./temp/temp_anti_bed.bed" -wa -split -v -nonamecheck > "./temp/$TEMP_NAME"
								samtools index "./temp/$TEMP_NAME"
							else
								bedtools intersect -a $DATA_LOCATION -b "./temp/temp_anti_bed.bed" -wa -split -v -s -nonamecheck > "./temp/$TEMP_NAME"
								samtools index "./temp/$TEMP_NAME"
						fi
					fi
				done
				
				if [[ "$SUBSET_BAMS" == FALSE ]]
				then
					echo "Skipping subsetting anti probe $antiprobe BAM files. If subsetting is desired remove -F flag."
				fi
				
				IFS=$'\t';
				PREVIOUS_ANTI_PROBE=${PREVIOUS_ANTI_PROBE}_anti_$antiprobe #update previous antiprobe right before the for loop for each antiprobe updates
		done
	fi
	
	if [ $RT_PCR = TRUE ] 
	then
		BUFFER_SIZE=5
		echo "======="
		echo "Running viewing window $VIEWING_WINDOW subset now for RT-PCR mode"
		VW_LINE=$(awk -v var="$VIEWING_WINDOW" '$4==var {print $0}' $PROBES)
		IFS=$'\t'; read -a veels <<< "$VW_LINE"
		if [[ -z "${veels[@]}" ]]
		then
			echo "Viewing window: $VIEWING_WINDOW not found. Check probes file. Exiting script"
			exit 
		fi
		# Psuedo code is esentially it takes the start position, puts that in a bed, does a -wa intersect
		# Then takes the end position, puts that in the bed, does another -wa intersect
		# Lastly, takes the full length view window, then does a hard intersect without -wa, which will then only print the 
		# regions that intersect in that view window 
		WINDOW_START=${veels[1]}
		WINDOW_END=${veels[2]}
		if [[ "$SUBSET_BAMS"  == TRUE ]]
		then
			IFS=$',';
			for sample in $BAMS; ##You need to only create subsets of the samples you are using 
			do
				DATA_LOCATION="./temp/${sample}_${PREVIOUS_ANTI_PROBE}.bam"
				TEMP_DATA_LOCATION="./temp/RTPCR_temp.bam"
				cp $DATA_LOCATION $TEMP_DATA_LOCATION
				echo -e "${veels[0]}\t$(($WINDOW_START-$BUFFER_SIZE))\t$WINDOW_START" > "./temp/temp_start.bed"
				echo -e "${veels[0]}\t$WINDOW_END\t$(($WINDOW_END+$BUFFER_SIZE))" > "./temp/temp_end.bed"
				echo "$VW_LINE" > "./temp/temp.bed"
				bedtools intersect -a $TEMP_DATA_LOCATION -b "./temp/temp_start.bed" -wa -split -nonamecheck > $DATA_LOCATION
				bedtools intersect -a $DATA_LOCATION -b "./temp/temp_end.bed" -wa -split -nonamecheck > $TEMP_DATA_LOCATION
				bedtools intersect -a $TEMP_DATA_LOCATION -b "./temp/temp.bed" -split -nonamecheck > $DATA_LOCATION
				rm $TEMP_DATA_LOCATION "./temp/temp_start.bed" "./temp/temp_end.bed"
			done
		fi
	fi
	
	echo "======="
	echo "======="
	IFS=$'\t';
	echo "Running R scripts"
	BAMS=${fields[1]} #I dont know why I need this but I do
	Rscript $NANO_BLOT_RSCRIPT $BAMS $PROBE_FIELD $NORM_FACTOR ${PREVIOUS_ANTI_PROBE} $META_DATA $ANTIPROBE_FIELD
	# this order has to be this way because if there is no antiprobe, then it collapses to an empty
	# string and the number of arguments passed to the script decreases by one, that is why the antiprobe
	# argument has to be the last one
	echo "======="
	
done 
echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="

# This code essentially just removes all the temp files that were created by the program
if [ $CLEAN_ALL = TRUE ]
then
	echo "Clearing ./temp/"
	rm -r ./temp/*
fi