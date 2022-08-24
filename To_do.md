# To Do

### RTPCR
Add a flag that accepts a separate metadata file and cuts the reads down to a user specificed length before plotting.

### Re-add simple library size normalization function, when adding this feature, it will be like a mode that the user can select for normalization, so like -N size or -N differential 
#### We can make the default option -N differential, so the user will need to specify the annotation file

### Add a feature to have count tables generated again, this could be due to Htseq not running properly, or the user happened to use a new metadata file, etc. 

### Skipping subsetting currently doesn't work because the R script doesn't know where those bam files are 

### Add user_input_files to the gitignore

### Remove duplication factor and have the Rscript automatically generate an appopriate sampling rate

### Fix names in R script
I want it to add the negative probe to the plot name if a negative probe is used

## Write a manuscript + ReadMe

Describe limitations

### Add example datasets I'm gunna do this last so I can use it to test stuff.

ideally they would be small and approximately the same sequencing depth.
I'm regenrating data from WT and rrp6 directRNA sequencing. 
