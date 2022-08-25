# To Do

## Big Projects
### RTPCR
Add a flag that accepts a separate metadata file and cuts the reads down to a user specified length before plotting.

### Remove duplication factor and have the Rscript automatically generate an appopriate sampling rate
First, find a way using Rsamtools to count the correct number of bam reads for the size normalization method's size factors. 

## Small Projects 
### Skipping subsetting currently doesn't work because the R script doesn't know where those bam files are 

### Add user_input_files to the gitignore

## Final steps
### Fix names in R script
I want it to add the negative probe to the plot name if a negative probe is used

## Write a manuscript + ReadMe
Describe limitations

### Add example datasets I'm gunna do this last so I can use it to test stuff.
ideally they would be small and approximately the same sequencing depth.
I'm regenrating data from WT and rrp6 directRNA sequencing. 
