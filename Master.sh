#!/bin/bash

PLOTS="./user_input_files/plot_data.csv"
PROBES="./user_input_files/probes.bed"
META_DATA="./user_input_files/data_metadata.csv"
DUP_FACTOR="1"

NANO_BLOT_RSCRIPT="./scripts/nano_blot_generation.R"

PRINT_HELP=TRUE
SUBSET_BAMS=TRUE
MAKE_PLOT=TRUE

if [[ "$PRINT_HELP" == TRUE ]]
then
	echo "
Nanoblot (Version 1.0)
==========================================
-H  |  Print help menu
-F  |  Filter BAM files for plot generation
-P  |  Generate nanoblots
=========================================="
fi

declare -i END_PLOT=$(wc -l < $PLOTS)

for (( c=2; c<=$END_PLOT; c++ ))
do
  P_LINE=$(head -n $c $PLOTS | tail -n -1)
#  echo $P_LINE
  echo "======="  
  echo "Generating plot" $c
  
  IFS=$'\t'; read -a fields <<<"$P_LINE"
  echo 'Getting Probe:'${fields[2]}
  
  DUP_FACTOR=${fields[3]}
  TARGET=${fields[2]}
  BAMS=${fields[1]}
#  echo $TARGET
  
  declare -i END_PROBE=$(wc -l < $PROBES)
  
  for (( d=1; d<=$END_PROBE; d++ ))
  do
     PR_LINE=$(head -n $d $PROBES | tail -n -1)
     
     IFS=$'\t'; read -a feels <<<"$PR_LINE"
     TARGET_PROBE=${feels[3]}
     
     if [[ "$TARGET_PROBE" == "$TARGET" ]]
     then
      echo ${feels[0]}:${feels[1]}-${feels[2]}

      declare -i END_META=$(wc -l < $META_DATA)
      IFS=','; read -a samples <<<"$BAMS"
      
      for (( e=2; e<=$END_META; e++ ))
      do
        DATA_LINE=$(head -n $e $META_DATA | tail -n -1)
        IFS=$'\t'; read -a EELS <<<"$DATA_LINE"
        DATA_LOCATION=${EELS[1]}
        echo $DATA_LOCATION
        
        TEMP_NAME=${EELS[0]}_$TARGET_PROBE.bam
        echo $TEMP_NAME
        
        if [[ "$SUBSET_BAMS" == TRUE ]]
				then
        	samtools view -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
        else 
        	echo "Skipping filtering BAM files. If filtering is desired use -F."
        fi
        
      done
     fi
  done 
  if [[ "$MAKE_PLOT" == TRUE ]]
  then
  	echo "Making plots --- WIP"
  	#Rscript $NANO_BLOT_RSCRIPT $BAMS $TARGET $DUP_FACTOR
  else
  	echo "Skipping plot generation. If plot generation is desired use -P."
  fi

done 

