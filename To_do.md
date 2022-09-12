# To Do

## Code

### RTpcr
The error in the RTpcr not calculating read length correctly is because even if there is correct viewing window clipping, the actual CIGAR string of the read never changes, which causes the qwdith calculation in scanToBam to always stay the same. Don't know how to work around this?

### Check how density plots work
I noticed that we are using the same dataset for making the Nanoblot as the density plots (ridge/violin) I don't think we should do this.

## Logistic

### Add example datasets I'm gunna do this last so I can use it to test stuff.
ideally they would be small and approximately the same sequencing depth.
I'm regenrating data from WT and rrp6 directRNA sequencing. 

### Human NMD plots

### Update Read_me
Package versions
Basic usage 
Extended Methods section
