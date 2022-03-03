#!/bin/bash

PLOTS="./plots_example.csv"

while read plot_name loading_order probe_black probe_red probe_blue probe_yellow
do
  echo "Plot Name : $plot_name"
done < $PLOTS

#read in data
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