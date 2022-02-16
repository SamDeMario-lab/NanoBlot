
library(ggplot2)
library(Rsamtools)
library(GenomicRanges)

Make_probe <- function(chr, start_pos, stop_pos, strand) {
  probe <- GRanges(seqnames = chr, ranges = IRanges(start = start_pos, end = stop_pos), strand = strand)
  return(probe)
}

RPL18A <- Make_probe("chrXV", 93395, 94402, "-")
RPL18A_3Exon <- Make_probe("chrXV", 93395, 93843, "-")

bamPath_wt <- "/var/lib/minknow/data/11-17-2021_WT_direct_RNA_for_SAM/11-17-2021_WT/20211117_1637_MN36788_FAQ93339_102dd49a/fastq/pass/sorted_merged.bam"
bamPath_rrp6 <- "/var/lib/minknow/data/11-18-2021_rrp6d_RNA_seq_for_SAM/11-18-2021_rrp6d/20211118_1508_MN36788_FAQ91220_398a4946/fastq/pass/sorted_merged.bam"

wt_bam <- BamFile(bamPath_wt)
rrp6_bam <- BamFile(bamPath_rrp6)

params <- ScanBamParam(which = RPL18A, what = scanBamWhat())

aln_wt <- scanBam(wt_bam, param = params)
aln_rrp6 <- scanBam(rrp6_bam, param = params)

blot_data_wt <- data.frame("qname"=c(aln_wt$`chrXV:93395-94402`$qname,aln_rrp6$`chrXV:93395-94402`$qname),
                           "qwidth"=c(aln_wt$`chrXV:93395-94402`$qwidth,aln_rrp6$`chrXV:93395-94402`$qwidth))
blot_data_wt[1:length(aln_wt$`chrXV:93395-94402`$qname),"sample_ID"] <- "wt"
blot_data_wt[1:length(aln_wt$`chrXV:93395-94402`$qname),"sample_pos"] <- 1
blot_data_wt[((length(aln_wt$`chrXV:93395-94402`$qname)+1):(length(aln_rrp6$`chrXV:93395-94402`$qname)+length(aln_wt$`chrXV:93395-94402`$qname))),"sample_ID"] <- "rrp6"
blot_data_wt[((length(aln_wt$`chrXV:93395-94402`$qname)+1):(length(aln_rrp6$`chrXV:93395-94402`$qname)+length(aln_wt$`chrXV:93395-94402`$qname))),"sample_pos"] <- 2

coding <-ggplot(data = blot_data_wt, aes(x = sample_pos, y = qwidth))+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  scale_fill_distiller(palette=1, direction=1) +
  scale_x_continuous(limits = c(0,3))+
  scale_y_continuous(limits = c(0,800))+
  ylab(label = "Size in nts")+
  labs(title = "RPL18A - \"Nanoblot\"")


params <- ScanBamParam(which = RPL18A_3Exon, what = scanBamWhat())

aln_wt <- scanBam(wt_bam, param = params)
aln_rrp6 <- scanBam(rrp6_bam, param = params)

blot_data <- data.frame("qname"=c(aln_wt$`chrXV:93395-93843`$qname,aln_rrp6$`chrXV:93395-93843`$qname),
                        "qwidth"=c(aln_wt$`chrXV:93395-93843`$qwidth,aln_rrp6$`chrXV:93395-93843`$qwidth))
blot_data[1:length(aln_wt$`chrXV:93395-93843`$qname),"sample_ID"] <- "wt"
blot_data[1:length(aln_wt$`chrXV:93395-93843`$qname),"sample_pos"] <- 1
blot_data[((length(aln_wt$`chrXV:93395-93843`$qname)+1):(length(aln_rrp6$`chrXV:93395-93843`$qname)+length(aln_wt$`chrXV:93395-93843`$qname))),"sample_ID"] <- "rrp6"
blot_data[((length(aln_wt$`chrXV:93395-93843`$qname)+1):(length(aln_rrp6$`chrXV:93395-93843`$qname)+length(aln_wt$`chrXV:93395-93843`$qname))),"sample_pos"] <- 2

exon3exon <- ggplot(data = blot_data, aes(x = sample_pos, y = qwidth))+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  scale_fill_distiller(palette=1, direction=1) +
  scale_x_continuous(limits = c(0,3))+
  scale_y_continuous(limits = c(0,800))+
  ylab(label = "Size in nts")+
  labs(title = "RPL18A 3'Exon Probe - \"Nanoblot\"")

print(exon3exon)
ggsave(plot = exon3exon,filename = "exon2_probe.png", device = "png")
print(coding)
ggsave(plot = coding,filename = "wholegene_probe.png", device = "png")
