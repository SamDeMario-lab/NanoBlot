# To Do

### Change R script to have it write out the data used to generate the plots as an RDS. 
Should be a one line addition. Maybe a small chunk to make a name. The Sam branch has an alternative R script which does this.

### Allow user to comment out plots in the plot_data.csv file
I want to add a check to see if the leading character in the csv line is a # and if it is skip that plot. Should be easy just add and if stament to see if the first character of the plot line is a # and then break the for loop.

### Fix names in R script
I want it to add the negative probe to the plot name if a negative probe is used

## Write a manuscript + ReadMe

Describe limitations

### Add example datasets I'm gunna do this last so I can use it to test stuff.

ideally they would be small and approximately the same sequencing depth.
I'm regenrating data from WT and rrp6 directRNA sequencing. 
