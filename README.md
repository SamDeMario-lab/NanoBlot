<h1 align="center">
  <img src="https://i.pinimg.com/736x/86/0f/b9/860fb969002c23702a9db5908e335d8f--healthy-lifestyle-science-safety.jpg" width=200 height=200/><br>
  Nanoblot
</h1>

<h3 align="center">NanoBlot: An R Package for Visualization of RNA Isoforms from Long Read RNA-sequencing Data</h3>

<div align="center">
  <a href="https://bedtools.readthedocs.io/en/latest/" target="_blank">
    <img src="https://img.shields.io/badge/Dependencies-Bedtools-informational" />
  </a>
  <a href="http://www.htslib.org" target="_blank">
    <img src="https://img.shields.io/badge/Dependencies-Samtools-informational" />
  </a>
  <a href="http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html" target="_blank">
      <img src="https://img.shields.io/badge/Bioconductor-DESeq2-important">
  </a>
  <a href="https://zenodo.org/badge/latestdoi/459793151"><img src="https://zenodo.org/badge/459793151.svg" alt="DOI"></a>
  </a>
</div>
<br/>

RT-PCR and Northern blots have long been used to study RNA isoforms usage for single genes. Recently, advancements in long read sequencing have yielded unprecedented information about the usage and abundance of these RNA isoforms. However, visualization of long-read sequencing data remains challenging due to the high information density. To alleviate these issues we have developed NanoBlot, a simple, open-source, command line tool, which generates Northern blot and RT-PCR-like images from third generation sequencing data. NanoBlot accepts processed bam files. Plotting is based around ggplot2 and is easily customizable. Advantages of NanoBlots include: designing probes to visualize isoforms which would be impossible with traditional RT-PCR or Northern blots, excluding reads from the Nanoblots based on the presence or absence of a specified region and, multiplexing plots with multiple colors. We present examples of Nanoblots compared to actual northern blot data. In addition to traditional gel-like images, NanoBlot also outputs other visualizations such as violin plots. The use of Nanoblot should provide a simple answer to the challenge of visualization of long-read RNA sequencing data. 

## Basic Usage 

The NanoBlot R package can be loaded entirely using the devtools package built into R. First, make sure that you have the devtools package in R. <br>
``` install.packages("devtools") ```

Then, make sure you are in the right directory for the R NanoBlot package folder, which should end in something like "/NanoBlotPackage". Then, run this code to compile NanoBlot <br>
```devtools::install()``` 
After installing, you can load the package using a simple library command <br>
```library("NanoBlotPackage")``` 

NanoBlot is now ready for use! 

## Setting up conda environment 

Some core Nanoblot functions will rely on system() commands that are expected to be built into the user's console. These commands include packages like ```bedtools``` and ```samtools``` as there is no R equivalent that works as efficiently. Thus, if the user does not have previous path commands in their console environment, the user can utilize the prewritten conda environment for easy installation. The provided conda .yml environment is listed under the scripts folder and titled ```nanoblotenv.yml``` 

The steps to install the conda environment are as followed <br>
Step 1: Install anaconda if not already on computer. Then find the system path to conda using the command ```which conda``` in any terminal of your choice. <br>
Step 2: Call the function ```conda create -f {your path}/scripts/nanoblotenv.yml``` in either the terminal, or you can call it in R using the function <br>
```system2("{your path to conda}, args = c("env", "create", "-f", "{your path}/scripts/nanoblotenv.yml"))```

Step 3: Add the newly created environment into R's system path so that R will automatically search for it using this code block <br>
```
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "{path to conda environment}", sep = ":"))
```


## Extended Methods

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

[^1]: Meeta Mistry, Mary Piper, Jihe Liu, & Radhika Khetani. (2021, May 24). hbctraining/DGE_workshop_salmon_online: Differential Gene Expression Workshop Lessons from HCBC (first release). Zenodo. https://doi.org/10.5281/zenodo.4783481
 
Depending on the user's intended usage of NanoBlot, the user can either run normalization using one technical replicate or multiple replicates for normalization. The standard workflow includes running one technical replicate with a n=1, and skips the estimation of dispersion and negative bionomial wald test, and instead only uses an estimation of size factors which we call in the script as ```size factors```[^2]. This size factor produces an appropriate normalization factor for each sequencing sample.  

[^2]: Anders, S., Huber, W. Differential expression analysis for sequence count data. Genome Biol 11, R106 (2010). https://doi.org/10.1186/gb-2010-11-10-r106
 
There are many ways to run DeSeq2. For Nanoblot, we will be first generating count tables through the Rsubread R package, which is used for mapping, quantification, and variant analysis of sequencing data. [^3] In the backend, the function ```normalizeNanoblot``` will call the helper function ```calculateDESeqSizeFactors``` which takes as input an annotation file as well as the unnormalized files to calculate the normalization based on the entire sequencing library as opposed to a subset of it. The annotation file should be in the gtf format supplied using the **Ensembl** [^4] database. 

[^3]: Liao Y, Smyth GK, Shi W (2019). “The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads.” Nucleic Acids Research, 47, e47. doi: 10.1093/nar/gkz114.
[^4]: Ensembl 2022. Nucleic Acids Res. 2022, vol. 50(1):D988-D995 PubMed PMID: 34791404. doi:10.1093/nar/gkab1049
 
These count tables will then be used by DeSeq2’s **estimateSizeFactors()** method [^2], with the size factors printed to the console. Should the user choose to normalize technical replicates of more than 1 using DeSeq2, there is an optional argument to ```normalizeNanoblotData``` which takes a coldata argument, which will be passed into DeSeq2's command ```DeSeq2::DESeqDataSetFromMatrix```. The user should read the documentation for DeSeq2 for more information on how to specify the coldata argument. 

In the event the user is interested in plotting reads that are not accurately annotated or do not have annotation regions, the Rsubread counts table might misrepresent the true sequencing depth of the samples. In these instances we have provided a normalization method based off of counts per million that will calculate a scaling factor based off of sequencing depth alone. We understand this is an imperfect method as CPM should not be used for differential expression analysis and should be proceeded with caution. We recommend users who are interested in using this normalization method to cross-validate this approach with the DeSeq normalization as well. 

### Subsetting
<img src=https://bedtools.readthedocs.io/en/latest/_images/intersect-glyph.png>[^5]

[^5]: Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842

After the normalization step is completed off of the sequenced bam files, the next step is to create subsets of each sequenced bam files to the corresponding probes and antiprobes that the user has specified. If skipping subsetting is desired, usually in the case where the user wants to control downstream analysis such as plotting and save computation time, the flag option ```./Nanoblot.sh -F``` can be used

The subsetting code goes through each target probe first, then uses those targetted bam files and subsequently goes through each antiprobe. The heart of the subsetting command relies on the bedtools **intersect** function, which you can read about at their wiki page [^5]. Based on whether the user has specified the input sequence files to be treated as cDNA strands or not with the ```./Nanoblot.sh -C```, the bedtools command will be run with different flag options. 

While Nanoblot allows for multiple target probes and antiprobes, the bedtools intersect function must be run one at a time for each intersection performed. If you look at the intersect diagram shown above, providing mulitple regions at once to the intersect function will render an effective **OR** intersection when we are looking for **AND** intersection. In order to subset the bam files appropriately, an intersect function must be run for each probe, one at a time. 

It is important to note that formatting of the plotting file is extremely important, as extra spaces and commas will lead to errors in processing column data. For the basic non RT-PCR mode, the probe field is found in column 2 (starting with index = 1), and the antiprobe field which is optional is found in column 3. 

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

The RT-PCR mode essentially takes everything the base Nanoblot does, with normalization and subsetting of probes, and adds an additional step that does a hard cut at sites that only overlap with the viewing window. This is based off of how wet-bench RT-PCR reactions are performed, where the flanking sites represent extension primers that amplify a certain region of interest. First, we filter the subsetted bam files to only reads that include both the start and the end of the viewing window, as primers in a wet bench RT-pCR reaction would not amplify sequences where the flanking sites are not present. We use a bedtools intersect start and end window with a **BUFFER_SIZE** variable of 5 nucleotides to allow for indels at the flanking regions. After selecting for reads that include the start and end sites, we then perform a hard genomic cut that only includes the regions inside the genomic window. To do this, we use the function ```ampliconclip``` from the ```samtools``` package to clip the end of read alignments if they intersect with regions defined in a BED file. [^6] 

Since ```ampliconclip``` will clip reads that overlap with the BED files, we will use the viewing window and find the complement bed entries to intersect, keeping effectively the reads that are only found within the viewing window. The options we use for ampliconclip include ```--hard-clip``` to ensure that the read width does not calculate clipped bases and ```--both-ends``` to ensure that both the 5' and the 3' ends where the regions match are cut. This genomic window is irrespective of strandedness since the reads are effectively cDNA at this point. 

[^6]: Whitwham, Andrew, and Rob Davies. “Samtools Ampliconclip.” Samtools-Ampliconclip(1) Manual Page, Sanger Institute , https://www.htslib.org/doc/samtools-ampliconclip.html

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

### RACE: Rapid Amplification of CDNA Ends

Alongside the RT-PCR mode, there is an added RACE (Rapid Amplification of CDNA Ends) option with the ```-P``` flag. It uses the same inputs as the RT-PCR with a separate plotting file template, the only difference is that it does not check for inclusive ends as it is effectively a one-sided PCR. Conceptually, this makes sense as the RACE protocol is used to determine regions of unknown sequences. For example, a 3' RACE experiment to determine unknown 3' mRNA sequences that lie between the exon and the poly(A) tail uses a gene-specific primer that anneals to a region of known exon sequences and an adapter primer that targets the poly(A) tail. [^7]  For more information on usage, read 5' or 3' RACE protocols to see what you need. Reference manuscript for example usage. 


[^7]: “3´ RACE System for Rapid Amplification of CDNA Ends.” Thermo Fisher Scientific - US, https://www.thermofisher.com/us/en/home/references/protocols/nucleic-acid-amplification-and-expression-profiling/cdna-protocol/3-race-system-for-rapid-amplification-of-cdna-ends.html. 

### R Plotting

Because Nanoblot deals with discrete read representations instead of continuous count numbers, the best way to visually represent the normalization is to duplicate the existing number of reads by a certain number that we call the duplication factor. This number is first calculated by taking the inverse of the size factor, since the size factor is divided during DeSeq2’s count normalization, and then multiplying the inversed number by 10 and rounding to the nearest digit to account for minimal data loss. Although this duplication factor would scale each sample respectively to their normalized counts, certain samples that inherently have lower reads than others would be visually hard to see. We then assigned an arbitrary **DUPLICATION_CONSTANT** a value of 2000 in order to scale all samples up or down respectively to ensure equal plotting exposure. This is essentially like automatically determining the exposure for Northern blots. The duplication factor from the previous step is then multiplied by taking the **DUPLICATION_CONSTANT** divided by the max number of reads across all samples plotted and rounded to the nearest integer. The plotting script then takes each sample’s reads and upscales it n times according to the duplication factor. 

For anyone interested in writing their own R script, you can do so by using the ```./Nanoblot.sh -R``` flag. The inputs that Nanoblot.sh calls to the Rscript are as follows ```Rscript $NANO_BLOT_R_SCRIPT $BAMS $PROBE_FIELD $NORM_FACTOR $PREVIOUS_ANTI_PROBE $META_DATA $ANTIPROBE_FIELD```
Since the first argument is the R script itself, bash is essentially passing 6 arguments to the R script
- Argument 1 **$BAMS**: These are all the samples that are being plotted separated by a comma, e.g. "WT,RRP6,SLU7,RRP6SLU7"
- Argument 2 **$PROBE_FIELD**: These are all the target probes, separated by commas if there is more than 1, e.g. "RPL18A_Exon1,RPL18A_Intron"
- Argument 3 **$NORM_FACTOR**: This is a key that tells the R script what type of normalization was done, 0 is for the default differential (DESeq2), 1 is for size (CPM), and 2 is for skipping normalization
- Argument 4 **$PREVIOUS_ANTI_PROBE**: This is the string that contains the naming convention to find the subset bam files, e.g. "RPL18A_Exon1_anti_RPL18A_Intron"
- Argument 5 **$META_DATA**: This is the metadata file that will be passed to the R script in case it needs to access it (only needs it for CPM normalization)
- Argument 6 **PLOT_TITLE**: This is the title of the plot which will be used to create the folder that the plots are stored in
- Argument 7 **$ANTIPROBE_FIELD**: This is the only argument that could be an empty string; these are all the target antiprobes separated by commas if there is more than 1 or empty if there are none, e.g. "RPL18A_Intron"

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
