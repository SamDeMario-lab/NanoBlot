<h1 align="center">
  <img src="https://i.pinimg.com/736x/86/0f/b9/860fb969002c23702a9db5908e335d8f--healthy-lifestyle-science-safety.jpg" width=200 height=200/><br>
  Nano_blot
</h1>

<h4 align="center">A Command Line Tool for Visualization of Isoform Usage From Oxford Nanopore RNA-seq</h4>

<div align="center">
  <a href="https://bedtools.readthedocs.io/en/latest/" target="_blank">
    <img src="https://img.shields.io/badge/Dependencies-Bedtools-informational" />
  </a>
  <a href="http://www.htslib.org" target="_blank">
    <img src="https://img.shields.io/badge/Dependencies-Samtools-informational" />
  </a>
  <a href="https://htseq.readthedocs.io/en/master/" target="_blank">
    <img src="https://img.shields.io/badge/Dependencies-HTSeq-informational" />
  </a>
  <a href="http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html" target="_blank">
      <img src="https://img.shields.io/badge/Bioconductor-DESeq2-important">
  </a>
</div>
<br/>

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
  Rsamtools (installed using Bioconductor)<br/>
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
| -O   |  Overwrite count tables and recalculate
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

## Hoffman2 Cluster Notes
The most important thing when running in Hoffman2 is loading all the required dependencies

**Bedtools (> 2.30.0)** --> to load this, run ```module load bedtools``` <br/>
**Samtools (> v1.15.1)** --> to load this, run ```module load samtools``` <br/>
**HTSeq (v> 2.0.2)** --> to load this, run ```module av python``` ---> ```pip install HTSeq``` <br/>
**R (> v4.1.2)** --> to load this, run ```module load intel/2022.1.1``` --> ```module load R/4.2.1``` <br/>
  Each dependency can then be installed after running R in the terminal and installing the packages for the first time <br/>
  ggplot2<br/>
  Rsamtools (installed using Bioconductor)<br/>
  ggridges<br/>
  Deseq2 (installed using Bioconductor)<br/>
  dplyr<br/>
