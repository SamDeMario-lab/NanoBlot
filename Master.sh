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


while getopts ":HFPCNR:M:B:T:" opt; do
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
    \? ) echo "Usage: "
      ;;
  esac
done

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
==========================================
"

exit
fi

declare -i END_PLOT=$(wc -l < $PLOTS)

if [[ $NORM == TRUE ]]
then
	for (( c=2; c<=$END_PLOT; c++ ))
	do
		P_LINE=$(head -n $c $PLOTS | tail -n -1)
		IFS=$'\t'; read -a fields <<<"$P_LINE"
		BAMS=${fields[1]}
		echo "Starting Normalization"
		declare -A count_array
		echo "Samples: "$BAMS
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
					echo "Shit worked BB. Found "$f
					READ_COUNT=$DATA_LOCATION"	""TEST_THIS_NEEDS_TO_BE_A_NUMBER"
					count_array[$DATA_LOCATION]=$READ_COUNT
				fi
			done
		done
		echo ${count_array[*]}
	done
fi

for (( c=2; c<=$END_PLOT; c++ ))
do
  P_LINE=$(head -n $c $PLOTS | tail -n -1)
  echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="
  PLOT_NUM=$((c-1))
  echo "Generating plot" $PLOT_NUM
  echo "======="
  
  IFS=$'\t'; read -a fields <<<"$P_LINE"
  
  echo 'Getting Probe:'${fields[2]}
  echo 'Duplication Factor:'${fields[3]}
  
  DUP_FACTOR=${fields[3]}
  TARGET=${fields[2]}
  BAMS=${fields[1]}
  
  declare -i END_PROBE=$(wc -l < $PROBES)
  END_PROBE=$((END_PROBE+1))
  
  for (( d=1; d<=$END_PROBE; d++ ))
  do
     if [[ $d == $END_PROBE ]]
     then
       echo "Probe not found. Check bed file and blots metadata file."
       break
     fi
     
     PR_LINE=$(head -n $d $PROBES | tail -n -1)
     
     IFS=$'\t'; read -a feels <<<"$PR_LINE"
     TARGET_PROBE=${feels[3]}
     
     if [[ "$TARGET_PROBE" == "$TARGET" ]]
     then
       echo ${fields[2]}": "${feels[0]}:${feels[1]}-${feels[2]}
       IFS=','; read -a samples <<<"$BAMS" #I'm like 80% sure I don't need this line
      
       if [[ "$SUBSET_BAMS"  == TRUE ]]
       then
	       for (( e=2; e<=$END_META; e++ ))
	       do
	         DATA_LINE=$(head -n $e $META_DATA | tail -n -1)
	         IFS=$'\t'; read -a EELS <<<"$DATA_LINE"
	         DATA_LOCATION=${EELS[1]}
	         TEMP_NAME=${EELS[0]}_$TARGET_PROBE.bam
	         echo "Subsetting: "$DATA_LOCATION"
Naming Subset: " $TEMP_NAME
	
					 if [[ "$CDNA"  == TRUE ]]
					 then
			         samtools view -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
					 else
							 echo "Strand: "${feels[5]}
							 if [[ "${feels[5]}"  == "+" ]]
							 then
							   samtools view -F 20 -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
							 elif [[ "${feels[5]}"  == "-" ]]
							 then
							   samtools view -f 16 -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
							 else
							   echo "Warning: No strand found for probe treating region disregarding strand information"
							   samtools view -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
							 fi
					 fi
	       done
       else 
         echo "Skipping filtering BAM files. If filtering is desired remove -F flag."
       fi
       break
     fi
  done
  
  #Normalize Reads if desired
  
  if [[ "$MAKE_PLOT" == TRUE ]]
  then
  	echo ""
  	echo "======="
  	echo "Running R script"
  	echo "======="
  	Rscript $NANO_BLOT_RSCRIPT $BAMS $TARGET $DUP_FACTOR
  	echo "======="
  	echo "======="
  else
  	echo "Skipping plot generation. If plot generation is desired remove -P flag."
  fi
echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="
echo ""
done 

