# To Do

## Code

### RACE flag
Add a RACE flag that essentially doesn't do a hard check for 3' or 5' overlap when doing the RT-PCR mode
A clear example of this is for snR37 between rrp6 and trf4, look at the manuscript for more information --> 
1) Make edits to script
2) Reflect changes in Readme

### RTpcr
Works but only with a newer patched version of samtools, for now it works in yeast but not in human NMD data
An additional thought I had, when we are running RT mode, do we want to add a check to ask if the user would like to run in cDNA mode? 

### I want to change the RDS output. I think it should be after the normailization but before the oversampling

## Logistic

### Add example datasets I'm gunna do this last so I can use it to test stuff.
ideally they would be small and approximately the same sequencing depth.
I'm regenrating data from WT and rrp6 directRNA sequencing. 
Example data sets should be generated from merged subset files of the samples that are being visualized, this way it retains all the data needed for example visualizations while keeping the size of the data small 
