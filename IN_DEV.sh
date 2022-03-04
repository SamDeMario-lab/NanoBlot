#!/bin/bash

PLOTS="./plots_example.csv"
<<<<<<< HEAD
META_DATA= "./data/Data_metadata.csv"

while read plot_name loading_order probe_black probe_red probe_blue probe_yellow
do
  echo "---"
  echo "Plot Name : $plot_name"
  echo "Loading order : $loading_order"
  echo "Probe : $probe_black"
  echo "---"
=======

while read plot_name loading_order probe_black probe_red probe_blue probe_yellow
do
  echo "Plot Name : $plot_name"
>>>>>>> b89364d5b7a8057d2bffa1b4ed78e3b3f4b305bf
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