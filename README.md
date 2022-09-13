<h1 align="center">
  <img src="https://i.pinimg.com/736x/86/0f/b9/860fb969002c23702a9db5908e335d8f--healthy-lifestyle-science-safety.jpg" width=200 height=200/><br>
  Nanoblot
</h1>

<h3 align="center">A Command Line Tool for Visualization of Isoform Usage From Oxford Nanopore RNA-seq</h3>

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
 
For the purposes of Nanoblot, because we are only accepting technical replicates of n=1, we will not be using the complete DeSeq2 standard workflow that includes estimating size factors, estimating dispersions, and then conducting a negative binomial Wald test [^2]. We will only be estimating size factors which produces an appropriate normalization factor for each sequencing sample.  

[^2]: Anders, S., Huber, W. Differential expression analysis for sequence count data. Genome Biol 11, R106 (2010). https://doi.org/10.1186/gb-2010-11-10-r106
 
There are many ways to run DeSeq2. For Nanoblot, we will be first generating count tables using the **htseq-count method**, which is part of the high-throughput sequence analysis in Python (HTSeq) [^3]. This method will require an annotation file as input, which the user can specify using the –A flag. HTSeq performs best when the annotation file is supplied using the **Ensembl** [^4] database. Nanoblot will then run HTSeq for each sample required using the sequence bam files and store the count tables .tsv file in the ``` temp/count-tables ``` folder.  

[^3]: G Putri, S Anders, PT Pyl, JE Pimanda, F Zanini Analysing high-throughput sequencing data in Python with HTSeq 2.0 https://doi.org/10.1093/bioinformatics/btac166 (2022)
[^4]: Ensembl 2022. Nucleic Acids Res. 2022, vol. 50(1):D988-D995 PubMed PMID: 34791404. doi:10.1093/nar/gkab1049
 
These count tables will then be used by DeSeq2’s **estimateSizeFactors()** method [^2], with the size factors printed to the console. In the event the user is interested in plotting reads that are not accurately annotated or do not have annotation regions, the HTSeq counts table might misrepresent the true sequencing depth of the samples. In this specific case example, we have provided a normalization method based off of counts per million that will calculate a scaling factor based off of sequencing depth alone. We understand this is an imperfect method as CPM should not be used for differential expression analysis and should be proceeded with caution. We recommend users who are interested in using this normalization method to cross-validate this approach with the DeSeq normalization as well. 
 
Should the user choose to normalize their sequenced bam files themselves, we also provided a way to skip normalization that sets each samples size factor to 1. 

### Subsetting
<img src=https://bedtools.readthedocs.io/en/latest/_images/intersect-glyph.png>[^5]

[^5]: https://bedtools.readthedocs.io/en/latest/

After the normalization step is completed off of the sequenced bam files, the next step is to create subsets of each sequenced bam files to the corresponding probes and antiprobes that the user has specified. If skipping subsetting is desired, usually in the case where the user wants to control downstream analysis such as plotting and save computation time, the flag option ```./Nanoblot.sh -F``` can be used

The subsetting code goes through each target probe first, then uses those targetted bam files and subsequently goes through each antiprobe. The heart of the subsetting command relies on the bedtools **intersect** function, which you can read about at their wiki page [^6]. Based on whether the user has specified the input sequence files to be treated as cDNA strands or not with the ```./Nanoblot.sh -C```, the bedtools command will be run with different flag options. 

While Nanoblot allows for multiple target probes and antiprobes, the bedtools intersect function must be run one at a time for each intersection performed. If you look at the intersect diagram shown above, providing mulitple regions at once to the intersect function will render an effective **OR** intersection when we are looking for **AND** intersection. In order to subset the bam files appropriately, an intersect function must be run for each probe, one at a time. 

It is important to note that formatting of the plotting file is extremely important, as extra spaces and commas will lead to errors in processing column data. For the basic non RT-PCR mode, the probe field is found in column 2 (starting with index = 1), and the antiprobe field which is optional is found in column 3. 

[^6]: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

### RT-PCR
An additional feature of Nanoblot is the ability to process reads in RT-PCR mode, which the user can specify with the ```./Nanoblot.sh -Y ``` command. This flag does take positional arguments, namely the location of the RT-PCR plotting file. The plotting file for RT-PCR is almost identical to the original plotting format, with the exception of one extra added column. The RT-PCR plotting file takes arguments with the following column information: plot_name, loading_order, viewing window, probe, antiprobe 

An example of this format is shown below 
```
plot_name loading_order	viewing_window	probe	antiprobe		
RPS18A_full	WT,RRP6,SLU7,RRP6SLU7	RPS18A_vw	RPS18A_Exon1
RPS18A_spliced	WT,RRP6,SLU7,RRP6SLU7	RPS18A_vw	RPS18A_Exon1	RPS18A_intron		
RPS18A_unspliced	WT,RRP6,SLU7,RRP6SLU7	RPS18A_vw	RPS18A_Exon1,RPS18A_intron
```

The viewing window is a genomic window with the same features as a probe, and for that reason, is included in the probes metadata file. As a result, if the genomic window is not included, the code will exit with a corresponding error that the viewing window was not found. 

The RT-PCR mode essentially takes everything the base Nanoblot does, with normalization and subsetting of probes, and adds an additional step that does a hard cut at sites that only overlap with the viewing window. This is based off of how wet-bench RT-PCR reactions are performed, where the flanking sites represent extension primers that amplify a certain region of interest. First, we filter the subsetted bam files to only reads that include both the start and the end of the viewing window, as primers in a wet bench RT-pCR reaction would not amplify sequences where the flanking sites are not present. We use a bedtools intersect start and end window with a **BUFFER_SIZE** variable of 5 nucleotides to allow for indels at the flanking regions. After selecting for reads that include the start and end sites, we then perform a hard genomic cut that only includes the regions inside the genomic window. To do this, we use the function ```ampliconclip``` from the ```samtools``` package to clip the end of read alignments if they intersect with regions defined in a BED file. [^7] 

Since ```ampliconclip``` will clip reads that overlap with the BED files, we will use the viewing window and find the complement bed entries to intersect, keeping effectively the reads that are only found within the viewing window. The options we use for ampliconclip include ```--hard-clip``` to ensure that the read width does not calculate clipped bases and ```--both-ends``` to ensure that both the 5' and the 3' ends where the regions match are cut. This genomic window is irrespective of strandedness since the reads are effectively cDNA at this point. 

[^7]: http://www.htslib.org/doc/samtools-ampliconclip.html

An example of an input BED file for ampliconclip is as follows. Let's say our viewing window of interest is this bed file, which is the entire RPL18A gene of the Saccharomyces cerevisiae. 
<img width="946" alt="Screen Shot 2022-09-12 at 9 51 12 PM" src="https://user-images.githubusercontent.com/26608622/189811599-bd48c80a-38ca-45c1-9158-910a9f3eb722.png">

```
chrXV	93395	94402	RPL18A_vw	.	-
```
We then create a bed tool that represents the complement genomic positions of that viewing window, which this ampliconclip will hard clip all reads that intersect this region, effectively, keeping reads that only are contained within the viewing window. The max length found in the second row represents the maximum genome size of any organism, which is found in the Japanese flower, Paris japonica, with 149 billion base pairs. 
```
chrXV 0 93395 
chrXV 94402 149000000000
```

Here is an example of what the result of ampliconclip will look like in IGV. 
<img width="736" alt="Screen Shot 2022-09-12 at 9 55 45 PM" src="https://user-images.githubusercontent.com/26608622/189812111-76b54d01-4018-426d-b4a9-5c215d339e60.png">


### R Plotting

Because Nanoblot deals with discrete read representations instead of continuous count numbers, the best way to visually represent the normalization is to duplicate the existing number of reads by a certain number that we call the duplication factor. This number is first calculated by taking the inverse of the size factor, since the size factor is divided during DeSeq2’s count normalization, and then multiplying the inversed number by 10 and rounding to the nearest digit to account for minimal data loss. Although this duplication factor would scale each sample respectively to their normalized counts, certain samples that inherently have lower reads than others would be visually hard to see. We then assigned an arbitrary **DUPLICATION_CONSTANT** a value of 2000 in order to scale all samples up or down respectively to ensure equal plotting exposure. This is essentially like automatically determining the exposure for Northern blots. The duplication factor from the previous step is then multiplied by taking the **DUPLICATION_CONSTANT** divided by the max number of reads across all samples plotted and rounded to the nearest integer. The plotting script then takes each sample’s reads and upscales it n times according to the duplication factor. 

For anyone interested in writing their own R script, you can do so by using the ```./Nanoblot.sh -R``` flag. The inputs that Nanoblot.sh calls to the Rscript are as follows ```Rscript $NANO_BLOT_R_SCRIPT $BAMS $PROBE_FIELD $NORM_FACTOR $PREVIOUS_ANTI_PROBE $META_DATA $ANTIPROBE_FIELD```
Since the first argument is the R script itself, bash is essentially passing 6 arguments to the R script
- Argument 1 **$BAMS**: These are all the samples that are being plotted separated by a comma, e.g. "WT,RRP6,SLU7,RRP6SLU7"
- Argument 2 **$PROBE_FIELD**: These are all the target probes, separated by commas if there is more than 1, e.g. "RPL18A_Exon1,RPL18A_Intron"
- Argument 3 **$NORM_FACTOR**: This is a key that tells the R script what type of normalization was done, 0 is for the default differential (DESeq2), 1 is for size (CPM), and 2 is for skipping normalization
- Argument 4 **$PREVIOUS_ANTI_PROBE**: This is the string that contains the naming convention to find the subset bam files, e.g. "RPL18A_Exon1_anti_RPL18A_Intron"
- Argument 5 **$META_DATA**: This is the metadata file that will be passed to the R script in case it needs to access it (only needs it for CPM normalization)
- Argument 6 **$ANTIPROBE_FIELD**: This is the only argument that could be an empty string; these are all the target antiprobes separated by commas if there is more than 1 or empty if there are none, e.g. "RPL18A_Intron"

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
