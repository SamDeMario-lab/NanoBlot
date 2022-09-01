# To Do

## Big Projects
### RTPCR
Add a flag that accepts a separate metadata file and cuts the reads down to a user specified length before plotting.

## Small Projects 
### Skipping subsetting currently doesn't work because the R script doesn't know where those bam files are 

### Update R version to R (>4.2.1)
Bioconductor was updated and now needs >4.2

### Add user_input_files to the gitignore
I'm pretty sure this works for some of the branches already.

### \? ) echo "Invalid option: -$OPTARG" should output help text and exit script
Current it just says you have an Invalid option and then just runs the script anyway

### Changed input file name to "tsv" form "csv"
Don't forget to also change the default user inputs in the Nanoblot.sh script

## Final steps
### Fix names in R script
I want it to add the negative probe to the plot name if a negative probe is used

### Add dependences to the -H help text

## Write a manuscript + ReadMe
Describe limitations

### Add example datasets I'm gunna do this last so I can use it to test stuff.
ideally they would be small and approximately the same sequencing depth.
I'm regenrating data from WT and rrp6 directRNA sequencing. 
