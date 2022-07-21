#!/bin/bash

#Default Values 
PLOTS="./user_input_files/plot_data.csv"
PROBES="./user_input_files/probes.bed"
META_DATA="./user_input_files/data_metadata.csv"
NANO_BLOT_RSCRIPT="./scripts/nano_blot_generation.R"
PRINT_HELP=FALSE
SUBSET_BAMS=TRUE
MAKE_PLOT=TRUE
CDNA=FALSE
NORM=TRUE
CLEAN_ALL=FALSE

while getopts ":HFPCNWR:M:B:T:" opt; do
	case ${opt} in
		H ) 
		PRINT_HELP=TRUE
			;;
		F ) 
		SUBSET_BAMS=FALSE
			;;
		P ) 
		MAKE_PLOT=FALSE
			;;
		C ) 
		CDNA=TRUE
			;;
		R ) 
		NANO_BLOT_RSCRIPT=${OPTARG}
			;;
		M ) 
		META_DATA=${OPTARG}
			;;
		B ) 
		PLOTS=${OPTARG}
			;;
		T ) 
		PROBES=${OPTARG}
			;;
		N ) 
		NORM=FALSE
			;;
		W )
		CLEAN_ALL=TRUE
			;;
		\? ) echo "Usage: "
			;;
	esac
done

echo -e "Starting Nanoblot "'\u2622' 
echo "======="
echo "R Script: $NANO_BLOT_RSCRIPT";
echo "Meta Data File: $META_DATA";
declare -i END_META=$(wc -l < $META_DATA)
echo "Probes Bed File: $PROBES";
echo "Plots File: $PLOTS";

#Check if user asked for the help text and echo the help text.
if [[ "$PRINT_HELP" == TRUE ]]
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

For an explanation of the required input files see the README.md

==========================================
-H  |  Print help menu
-T  |  Probes bed file 
-B  |  Blots metadata file
-M  |  Location of metadata file
-R  |  Use custem R script
-N  |  Skip data normalization (Work in progress)
-C  |  Treat reads as cDNA (disregard strand) 
-F  |  Skip subsetting BAM files for plot generation
-P  |  Skip nanoblots generation
-W  |  Clear all files from ./temp/ after plot generation
==========================================
"
exit
fi

declare -i END_PLOT=$(wc -l < $PLOTS)

if [[ $NORM == TRUE ]]
then
	echo "======="
	echo "Normalization"
	for (( c=2; c<=$END_PLOT; c++ ))
	do
		P_LINE=$(head -n $c $PLOTS | tail -n -1)
		IFS=$'\t'; read -a fields <<<"$P_LINE"
		BAMS=${fields[1]}
		NORM_FOLDER="./temp/"$(echo "$BAMS" | sed -e 's/,/_/g')"_NORM"
		if [ ! -d "$NORM_FOLDER" ]
		then
			echo $NORM_FOLDER " not found"
			mkdir -p $NORM_FOLDER;
			NORM_METADATA_FILE=$NORM_FOLDER"/data_metadata.csv"
			NORM_METADATA_HEADER="Sample_name (This must be unique for each sample)  Type (FAST5 or BAM)  Location (For BAM inputs the path to the bam file should be given.)"
			echo "Starting Normalization"
			echo "Samples: "$BAMS
			IFS=','
			echo ${NORM_METADATA_HEADER} >$NORM_METADATA_FILE
			COMPARE=999999999999999 #This is probobly not the best way to do this
			for f in $BAMS;
			do
				for (( o=2; o<=$END_META; o++ ))
				do
					DATA_LINE_o=$(head -n $o $META_DATA | tail -n -1)
					IFS=$'\t'; read -a EELS <<<"$DATA_LINE_o"
					DATA_LOCATION=${EELS[1]}
					SAMP_NAME=${EELS[0]}
					if [[ $f == $SAMP_NAME ]]
					then
						echo "Found "$f
						READ_COUNT=$(samtools view -c -F 260 $DATA_LOCATION)
						if [[ $READ_COUNT -lt $COMPARE ]]
						then
							COMPARE=$READ_COUNT
						fi
					fi
				done
			done
			echo $COMPARE
			IFS=','
			for f in $BAMS;
			do
				for (( e=2; e<=$END_META; e++ ))
				do
					DATA_LINE=$(head -n $e $META_DATA | tail -n -1)
					IFS=$'\t'; read -a EELS <<<"$DATA_LINE"
					DATA_LOCATION=${EELS[1]}
					SAMP_NAME=${EELS[0]}
					if [[ $f == $SAMP_NAME ]]
					then
						echo "Subsetting "$SAMP_NAME" to "$COMPARE" reads."
						OUTPUT_NORM_SAM=$NORM_FOLDER"/"$f"_NORM.sam"
						OUTPUT_NORM_BAM=$NORM_FOLDER"/"$f"_NORM.bam"
						cat <(samtools view -H $DATA_LOCATION) <(samtools view ${DATA_LOCATION} | shuf -n $COMPARE) > $OUTPUT_NORM_SAM #THIS TAKES A LOT OF RAM
						samtools view -S -b $OUTPUT_NORM_SAM > $OUTPUT_NORM_BAM
						samtools sort $OUTPUT_NORM_BAM -o $OUTPUT_NORM_BAM
						samtools index $OUTPUT_NORM_BAM
						echo $SAMP_NAME"		"$OUTPUT_NORM_BAM >> $NORM_METADATA_FILE
					fi
				done
			done
		else 
			echo ""
			echo $NORM_FOLDER" already exists."
			echo "If data normalization should be re-run delete "$NORM_FOLDER
			echo ""
		fi
	done
fi

for (( j=2; j<=$END_PLOT; j++ ))
do
	P_LINE=$(head -n $j $PLOTS | tail -n -1)
	echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="
	PLOT_NUM=$((j-1))
	echo "Generating plot" $PLOT_NUM
	echo "======="
	
	IFS=$'\t'; read -a fields <<<"$P_LINE"
	
	echo 'Probe:'${fields[2]}
	if [ -z "${fields[4]}" ]
	then
		echo "No Negative Probe Used"
	else
		echo 'Negative Probe:'${fields[4]}
	fi
	echo 'Duplication Factor:'${fields[3]}
	echo "======="
	
	DUP_FACTOR=${fields[3]}
	TARGET=${fields[2]}
	TARG_NEGS=${fields[4]}
	BAMS=${fields[1]}
	
	if [[ $NORM == "TRUE" ]]
	then
		NORM_FOLDER="./temp/"$(echo "$BAMS" | sed -e 's/,/_/g')"_NORM"
		NORM_METADATA_FILE=$NORM_FOLDER"/data_metadata.csv"
		META_DATA=$NORM_METADATA_FILE
	fi

#Apply the first filter
	declare -i END_PROBE=$(wc -l < $PROBES)
	END_PROBE=$((END_PROBE+1))

	for (( d=1; d<=$END_PROBE; d++ ))
	do
		if [[ $d == $END_PROBE ]]
		then
			echo "Probe not found. Check bed file and blots metadata file."
		fi

		IFS=$','
		PR_LINE=$(head -n $d $PROBES | tail -n -1)

		IFS=$'\t'; read -a feels <<<"$PR_LINE"
		TARGET_PROBE=${feels[3]}

		if [[ "$TARGET_PROBE" == "$TARGET" ]]
		then
			echo 'Probe'
			echo ${fields[2]}": "${feels[0]}:${feels[1]}-${feels[2]}

			IFS=$','; echo $PR_LINE > "./temp/temp_bed.bed"

			if [[ "$SUBSET_BAMS"  == TRUE ]]
			then
				for (( t=2; t<=$END_META; t++ ))
				do
					DATA_LINE_T=$(head -n $t $META_DATA | tail -n -1)
					IFS=$'\t'; read -a EELS <<<"$DATA_LINE_T"
					DATA_LOCATION=${EELS[1]}
					TEMP_NAME=${EELS[0]}_$TARGET_PROBE.bam
					echo "Subsetting: "$DATA_LOCATION"
Naming Subset: " $TEMP_NAME
					if [[ "$CDNA"  == TRUE ]]
					then
						bedtools intersect -a $DATA_LOCATION -b "./temp/temp_bed.bed" -wa -split -nonamecheck > "./temp/$TEMP_NAME"
						samtools index "./temp/$TEMP_NAME"
					else
						bedtools intersect -a $DATA_LOCATION -b "./temp/temp_bed.bed" -wa -split -s -nonamecheck > "./temp/$TEMP_NAME"
						samtools index "./temp/$TEMP_NAME"
					fi
				done
			else 
				echo "Skipping filtering BAM files. If filtering is desired remove -F flag."
			fi
				break
		fi
	done
	echo "======="
	if [ -z "${fields[4]}" ]
	then
		echo "No negative probe used"
	else
		echo 'Negative Probe'

		declare -i END_PROBE=$(wc -l < $PROBES)
		END_PROBE=$((END_PROBE+1))

		for (( h=1; h<=$END_PROBE; h++ ))
		do
			if [[ $h == $END_PROBE ]]
			then
				echo "Probe not found. Check bed file and blots metadata file."
				break
			fi

			IFS=$','
			APR_LINE=$(head -n $h $PROBES | tail -n -1)
			echo $APR_LINE > "./temp/temp_anti_bed.bed"

			IFS=$'\t'; read -a deels <<<"$APR_LINE"
			TARG_ANTI=${deels[3]}

			if [[ $TARG_NEGS == $TARG_ANTI ]]
			then
				echo $TARG_ANTI": "${deels[0]}:${deels[1]}-${deels[2]}
				if [[ "$SUBSET_BAMS"  == TRUE ]]
				then
					for (( p=2; p<=$END_META; p++ ))
					do
						DL_ANTI=$(head -n $p $META_DATA | tail -n -1)
						IFS=$'\t'; read -a bells <<<"$DL_ANTI"
						FIRST_NAME="./temp/"${bells[0]}"_"$TARGET".bam"
						echo "Subsetting:" $FIRST_NAME
						NEG_NAME="./temp/"${bells[0]}"_"$TARGET"_anti_"$TARG_ANTI".bam"
						echo "Naming Subset: "$NEG_NAME
						if [[ "$CDNA"  == TRUE ]]
						then
							bedtools intersect -a $FIRST_NAME -b "./temp/temp_anti_bed.bed" -wa -split -v -nonamecheck > $NEG_NAME
							samtools index $NEG_NAME
						else
							bedtools intersect -a $FIRST_NAME -b "./temp/temp_anti_bed.bed" -wa -split -v -s -nonamecheck > $NEG_NAME
							samtools index $NEG_NAME
						fi
					done
				fi
			break
			fi
		done
	fi
	echo "======="
	if [[ "$MAKE_PLOT" == TRUE ]]
	then
		echo "======="
		echo "Running R script"
		BAMS=${fields[1]} #I dont know why I need this but I do
		Rscript $NANO_BLOT_RSCRIPT $BAMS $TARGET $DUP_FACTOR $TARG_NEGS
		echo "======="
	else
		echo "Skipping plot generation. If plot generation is desired remove -P flag."
	fi
done 
echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="

if [[ "$CLEAN_ALL" == TRUE ]]
then
	echo "Clearing ./temp/"
	rm -r ./temp/*
fi