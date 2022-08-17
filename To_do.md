# To Do

### Change R script to have it write out the data used to generate the plots as an RDS. 
Should be a one line addition. Maybe a small chunk to make a name. The Sam branch has an alternative R script which does this.

### Find an alternative to normalization of reads between samples
Since our visualization method keeps each read as a discrete data point, there's no way for us to change how each read is represented. However, we can change how abundant that read appears to the user through its darkness, essentially a psuedo for duplication factor. My idea is to take DeSeq2's scaling factor for each sample and automatically applying it as a duplication factor to that sample's plot. 

### Fix names in R script
I want it to add the negative probe to the plot name if a negative probe is used

## Write a manuscript + ReadMe

Describe limitations

### Add example datasets I'm gunna do this last so I can use it to test stuff.

ideally they would be small and approximately the same sequencing depth.
I'm regenrating data from WT and rrp6 directRNA sequencing. 
