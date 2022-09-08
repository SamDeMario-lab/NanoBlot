<h1 align="center">
  <img src="https://i.pinimg.com/736x/86/0f/b9/860fb969002c23702a9db5908e335d8f--healthy-lifestyle-science-safety.jpg" width=200 height=200/><br>
  Nanoblot
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

```./Nanoblot.sh -M 'Test2' -B 'Test3' -T 'Test4'```

#### 1) A set of probes to be used in standard bed format "example.bed" 
```
  chrIV	1359922	1359969	RPS18A_Exon1	.	+
  chrVI	54686	54696	ACT1_Exon1	.	-
  chrIV	1236558	1236842	YRA1_Exon1	.	+
  chrXI	431906	432034	RPL14A_Exon1	.	+
```

#### 2) A tsv file listing the plots to be produced 
```
  plot_name	loading_order	probe	antiprobe
  ACT1_5exon	WT,RRP6,SLU7,RRP6SLU7	ACT1_Exon1	
  #RPL18A WT,RRP6 RPL18A_EXon1
```
Special note: Adding a # in the first character of a plot line will cause the script to skip plotting of that line. I.e. in the above example, adding #RPL18A will cause the script to skip normalization and plotting of the RPL18A plot
In addition, it is not necessary to have an antiprobe, leaving it blank is fine

#### 3) A tsv listing the names of input data file and their locations
```
  Sample_name	Location
  WT	/home/guillaume-chanfreau/Sequencing_Data/slu7_rrp6/pass/barcode01/sorted_merged.bam
```

## Methods

### Normalization
The purpose of normalization is to adjust for factors that prevent for direct comparison of expression of genes between two biological samples. Common normalization procedures include counts per million (CPM), transcripts per kilobase million (TPM), RPKM/FPKM, DeSeq2’s median of ratios, and EdgeR’s trimmed mean of M values (TMM). Because we are not concerned with gene comparisons within samples, the only other normalization method that allows for differential expression analysis of genes between samples is DESeq2’s median of ratios. This normalization method accounts for not only sequencing depth but also RNA composition. This will be Nanoblot’s default normalization method.  

**Common Normalization Methods** [^1]
| **Normalization method** | **Description** | **Accounted factors** | **Recommendations for use** | 
| ---  | --- | --- | --- |
| **CPM** (counts per million) | counts scaled by total number of reads | sequencing depth | gene count comparisons between replicates of the same samplegroup; **NOT for within sample comparisons or DE analysis** | 
| **TPM** (transcripts per kilobase million) | counts per length of transcript (kb) per million reads mapped | sequencing depth and gene length | gene count comparisons within a sample or between samples of the same sample group; **NOT for DE analysis** | 
| **RPKM/FPKM** (reads/fragments per kilobase of exon per million reads/fragments mapped) | similar to TPM | sequencing depth and gene length | gene count comparisons between genes within a sample; **NOT for between sample comparisons or DE analysis** | 
| **DESeq2’s median of ratios** | counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene | sequencing depth and RNA composition | gene count comparisons between samples and for **DE analysis; NOT for within sample comparisons** | 
| EdgeR’s **trimmed mean of M values (TMM)** | uses a weighted trimmed mean of the log expression ratios between samples | sequencing depth, RNA composition, and gene length | gene count comparisons **between and within samples and for DE analysis** | 

[^1]: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
 
For the purposes of Nanoblot, because we are only accepting technical replicates of n=1, we will not be using the complete DeSeq2 standard workflow that includes estimating size factors, estimating dispersions, and then conducting a negative binomial Wald test (CITATION NEEDED). We will only be estimating size factors which produces an appropriate normalization factor for each sequencing sample.  
 
There are many ways to run DeSeq2. For Nanoblot, we will be first generating count tables using the htseq-count method, which is part of the high-throughput sequence analysis in Python (HTSeq). This method will require an annotation file as input, which the user can specify using the –A flag. HTSeq performs best when the annotation file is supplied using the Ensembl database. Nanoblot will then run HTSeq for each sample required using the sequence bam files and store the count tables .tsv file in the temp/count-tables folder.  
 
These count tables will then be used by DeSeq2’s estimateSizeFactors() method, with the size factors printed to the console. Because Nanoblot deals with discrete read representations instead of continuous count numbers, the best way to visually represent the normalization is to duplicate the existing number of reads by a certain number that we call the duplication factor. This number is first calculated by taking the inverse of the size factor, since the size factor is divided during DeSeq2’s count normalization, and then multiplying the inversed number by 10 and rounding to the nearest digit to account for minimal data loss. Although this duplication factor would scale each sample respectively to their normalized counts, certain samples that inherently have lower reads than others would be visually hard to see. We then assigned an arbitrary DUPLICATION_CONSTANT a value of 2000 in order to scale all samples up or down respectively to ensure equal plotting exposure. This is essentially like automatically determining the exposure for Northern blots. The duplication factor from the previous step is then multiplied by taking the DUPLICATION_CONSTANT divided by the max number of reads across all samples plotted and rounded to the nearest integer. The plotting script then takes each sample’s reads and upscales it n times according to the duplication factor. 
 
In the event the user is interested in plotting reads that are not accurately annotated or do not have annotation regions, the HTSeq counts table might misrepresent the true sequencing depth of the samples. In this specific case example, we have provided a normalization method based off of counts per million that will calculate a scaling factor based off of sequencing depth alone. We understand this is an imperfect method as CPM should not be used for differential expression analysis and should be proceeded with caution. We recommend users who are interested in using this normalization method to cross-validate this approach with the DeSeq normalization as well. 
 
Should the user choose to normalize their sequenced bam files themselves, we also provided a way to skip normalization that sets each samples size factor to 1 and would only calculate duplication factors according to the maximum number of reads across all samples. 

### RT-PCR


## Computing Cluster
The most computationally intensive part of Nanoblot is generating count tables and subsetting bam files for target probes. Should you wish to submit this script to a computing cluster, here are some tips to help you do so. UCLA provides the Hoffman2 cluster for free to use with computing nodes available. The most important thing when running computing clusters is loading all the required dependencies. Here is an example of what it's like on Hoffman 2. 

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
