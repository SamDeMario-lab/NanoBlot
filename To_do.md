# To Do

## Code

### Where do we want the master script to run
Do we want it to be only run in the directory we create and pull from
Or do we want it to be run absolute path, which will inevitable create a bunch of sub directories in whatever directory the user is in --> we can do a simple check to prevent this

### RTpcr
Works but only with a newer patched version of samtools, for now it works in yeast but not in human NMD data

### Check how density plots work
I noticed that we are using the same dataset for making the Nanoblot as the density plots (ridge/violin) I don't think we should do this.

### I want to change the RDS output. I think it should be after the normailization but before the oversampling

## Logistic

### Add example datasets I'm gunna do this last so I can use it to test stuff.
ideally they would be small and approximately the same sequencing depth.
I'm regenrating data from WT and rrp6 directRNA sequencing. 
