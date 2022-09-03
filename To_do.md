# To Do

## Big Projects
### RTPCR
RT-PCR is also being really janky and finnicky on Sam's linux build, need to see what shell he is running it in

### Nanoplot plotting
If a sample has no reads, it still needs an empty column to simulate a real Northern probe --> do a set x-axis label and window size based on number of samples 
It would also be nice if the plots could be labelled with sample names? 
Show Sam what I have so far, I have empty columns for nanoblot, but not nanoridge and nanoviolin, I mean, that makes sense right? 

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
