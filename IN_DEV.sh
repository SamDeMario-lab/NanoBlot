#!/bin/bash

PLOTS="./plots_example.csv"
META_DATA="./data/Data_metadata.csv"

while read plot_name loading_order probe_black probe_red probe_blue probe_yellow
do
  echo "---"
  echo "Plot Name : $plot_name"
  echo "Loading order : $loading_order"
  echo "Probe : $probe_black"
  echo "---"
done < $PLOTS

while read sample_name BAM_vs_FAST5 BAM_location
do 
  echo "---"
  echo "Sample Name : $sample_name"
  echo "BAM vs FAST5 : $BAM_vs_FAST5"
  echo "File location: $BAM_location" 
done < $META_DATA


  #check formatting of data
#read in probes
  #check formatting of probes
#read in blots
  #check formatting of blots
  #check if data is present for each blot
  #check if probe is present for each blot
#For loop through blots 
  #Subset data based on loading order and probes save as temp data
  #Make an array of file paths?
#Pass file paths to an R script
#generate plots
