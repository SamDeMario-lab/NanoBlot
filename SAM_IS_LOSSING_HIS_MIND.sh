#!/bin/bash

PLOTS="./plots_example.csv"
PROBES="./probes/example.bed"

declare -i END=$(wc -l < $PLOTS)

for (( c=2; c<=$END; c++ ))
do
  P_LINE=$(head -n $c $PLOTS | tail -n -1)
  echo $P_LINE
  echo "Generating plot" $c
  
  IFS=$'\t'; read -a fields <<<"$P_LINE"
  echo 'Getting Probe:'${fields[2]}
  TARGET=${fields[2]}
  
done

