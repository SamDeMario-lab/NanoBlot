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

Nanoblot is in active development. 

# Usage:

Nanoblot can be run in its entirety via the included bash script "Master.sh"


It requires 3 inputs:

## 1) A set of probes to be used in standard bed format "example.bed" 

  chrIV	1359922	1359969	RPS18A_Exon1	.	+
  chrVI	54686	54696	ACT1_Exon1	.	-
  chrIV	1236558	1236842	YRA1_Exon1	.	+
  chrXI	431906	432034	RPL14A_Exon1	.	+

## 2) A csv file listing the plots to be produced 

  plot_name	loading_order	probe_black	probe_red	probe_blue	probe_yellow
  ACT1_5exon	WT,RRP6,SLU7,RRP6SLU7	ACT1_Exon1	

## 3) A csv listing the names of input data file and their locations

  Sample_name (This must be unique for each sample)	Type (FAST5 or BAM)	Location (For FAST5 inputs this should be a directory ending in a /, For BAM inputs the path to the bam file should be given.)
  WT		/home/guillaume-chanfreau/Sequencing_Data/slu7_rrp6/pass/barcode01/sorted_merged.bam
