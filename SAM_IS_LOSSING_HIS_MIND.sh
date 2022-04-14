#!/bin/bash

PLOTS="./plots_example.csv"
PROBES="./probes/example.bed"
META_DATA="./data/Data_metadata.csv"

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
#     echo "CHECKING TARGET:"$TARGET
#     echo "CHECKING PROBE:"$TARGET_PROBE
     
     if [[ "$TARGET_PROBE" == "$TARGET" ]]
     then
      echo "FOUND THAT SHIT"
     fi
  done 
  
  declare -i END_META=$(wc -l < $META_DATA)
  IFS=','; read -a samples <<<"$BAMS"
  
  for (( e=1; e<=$END_META; e++ ))
  do
    DATA_LINE=$(head -n $e $META_DATA | tail -n -1)
    echo $DATA_LINE
  done
done 

