# Nano_blot

Nanoblot is a simple command line based tool for the creation of digital northern blots.

Advancements in long read sequencing have yielded unprecedented information about isoform usage. However, 
due to the high information density, visualization of long-read sequencing data remains challenging. 
Northern blots have long been used to study isoform usage. Here we present NanoBlot, a simple, open-source, 
command line tool to produce northern blot-like images from Oxford Nanopore sequencing data. One advantage 
of NanoBlots is that probes can be designed to visualize isoforms which would normally be difficult to 
observe on traditional northern blots. For example, transcripts can be excluded from a blot based on the 
presence of a specified region. Additionally, multiple colors can be used to highlight specific isoforms. 
NanoBlot can accept either raw Nanopore data or processed bam files. It is based around ggplot and is 
easily customizable. 

## Dependencies 

**Bedtools (> 2.30.0)**<br/>
**Samtools (> v1.15.1)**<br/>
**HTSeq (v> 2.0.2)**<br/>
**R (> v4.1.2)**<br/>
  ggplot2<br/>
  Rsamtools<br/>
  ggridges<br/>
  Deseq2 (installed using Bioconductor)<br/>
  dplyr<br/>
  
The Nanoblot.sh script will require the different dependices to run from the bash terminal. When running using the command line interface, I recommend attaching it to $PATH. 
A sample function would look like<br/>
```
export PATH="$PATH:/Library/Frameworks/R.framework/Resources"
source ~/.bash_profile
```

## Usage:

Nanoblot can be run in its entirety via the included bash script "Nanoblot.sh"

Nanoblot (Version 1.0)

| Flag | Description |
| ---  | --- |
| -H   |  Print help menu |
| -T   |  Probes bed file |
| -B   |  Blots metadata file |
| -M   |  Location of metadata file |
| -A   |  Annotation file |
| -R   |  Use custom R script |
| -Y   |  RT-PCR mode, supply with own metadata file |
| -N   |  Normalization function {differential (default), size, skip} |
| -C   |  Treat reads as cDNA (disregard strand) |
| -F   |  Skip subsetting BAM files for plot generation |
| -W   |  Clear all files from ./temp/ after plot generation |


It requires 3 inputs:

```./Master.sh -FP -R 'Test1' -M 'Test2' -B 'Test3' -T 'Test4'```

##### 1) A set of probes to be used in standard bed format "example.bed" 
```
  chrIV	1359922	1359969	RPS18A_Exon1	.	+
  chrVI	54686	54696	ACT1_Exon1	.	-
  chrIV	1236558	1236842	YRA1_Exon1	.	+
  chrXI	431906	432034	RPL14A_Exon1	.	+
```

##### 2) A csv file listing the plots to be produced 
```
  plot_name	loading_order	probe	antiprobe
  ACT1_5exon	WT,RRP6,SLU7,RRP6SLU7	ACT1_Exon1	
```
Special note: Adding a # in the first character of a plot line will cause the script to skip plotting of that line. I.e. in the above example, adding #ACT1_5exon will cause the script to skip normalization and plotting of the ACT1_5exon plot
In addition, it is not necessary to have an antiprobe, leaving it blank is fine

##### 3) A csv listing the names of input data file and their locations
```
  Sample_name	Location
  WT	/home/guillaume-chanfreau/Sequencing_Data/slu7_rrp6/pass/barcode01/sorted_merged.bam
```
Make sure each line is tab deliminated and is exactly one tab apart
