#!/bin/bash

PLOTS="./user_input_files/plots_example.csv"
PROBES="./user_input_files/example.bed"
META_DATA="./user_input_files/Data_metadata.csv"

declare -i END_PLOT=$(wc -l < $PLOTS)

for (( c=2; c<=$END_PLOT; c++ ))
do
  P_LINE=$(head -n $c $PLOTS | tail -n -1)
#  echo $P_LINE
  echo "======="  
  echo "Generating plot" $c
  
  IFS=$'\t'; read -a fields <<<"$P_LINE"
  echo 'Getting Probe:'${fields[2]}
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
        
        samtools view -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
      done
     fi
  done 
   $BAMS $TARGET
done 

