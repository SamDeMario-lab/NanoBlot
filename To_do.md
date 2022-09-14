# To Do

## Code

### Fix Folder Names
I'm not sure how we want to do this. I think we should add an input to the R Script that is the plot name form the user input file.

### RTpcr
The error in the RTpcr not calculating read length correctly is because even if there is correct viewing window clipping, the actual CIGAR string of the read never changes, which causes the qwdith calculation in scanToBam to always stay the same. Don't know how to work around this?

Found a way to solve this
Psuedocode
Essentially this revolves around the method samtools ampliconclip
First need to generate a complement of the bed regions we want to keep, since ampliconclip will clip these regions
Then run ampliconclip with the complement bed file, and use hard clipping (this is important since it won't factor in the size calculation), also need the option of both ends, this clip is regardless of strand since reads are treated as cDNA at this point
Think about whether we still need the prefiltering step of inclusive end and start 
Make sure read lengths are correct, something about using samtools fixmate but not sure what that does
Make sure that after generating the clipped file, we re-index it

A lot more work but robust way of generating bedtools complement files is to use the bedtools complement command but we would need a separate genome text file that would take a lot of work

Check this post for clarification on why certain regions are not clipped: https://bioinformatics.stackexchange.com/questions/16112/precisely-clipping-bam-file-to-bed-coordinates

### Add check for input bams being sorted and indexed and if they're not then sort and index

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
