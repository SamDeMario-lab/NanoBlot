
library(ggplot2)
library(Rsamtools)
library(GenomicRanges)
library(ggblur)

Make_probe <- function(chr, start_pos, stop_pos, strand) {
  probe <- GRanges(seqnames = chr, ranges = IRanges(start = start_pos, end = stop_pos), strand = strand)
  return(probe)
}

#Make_probe is a basic function that takes 4 inputs and returns a GRange object based on the inputs. 

RPL18A <- Make_probe("chrXV", 93395, 94402, "-")
RPL18A_3Exon <- Make_probe("chrXV", 93395, 93843, "-")

IMD4_5exon <- Make_probe("chrXIII", 163717, 164176, "-")
snR54 <- Make_probe("chrXIII", 163535, 163620, "-")

bamPath_wt <- "/var/lib/minknow/data/11-17-2021_WT_direct_RNA_for_SAM/11-17-2021_WT/20211117_1637_MN36788_FAQ93339_102dd49a/fastq/pass/sorted_merged.bam"
bamPath_rrp6 <- "/var/lib/minknow/data/11-18-2021_rrp6d_RNA_seq_for_SAM/11-18-2021_rrp6d/20211118_1508_MN36788_FAQ91220_398a4946/fastq/pass/sorted_merged.bam"
bamPath_wt_IVP <- "/home/guillaume-chanfreau/Sequencing_Data/slu7_rrp6/pass/barcode01/sorted_merged.bam"
bamPath_slu7_IVP <- "/home/guillaume-chanfreau/Sequencing_Data/slu7_rrp6/pass/barcode04/sorted_merged.bam"

wt_bam <- BamFile(bamPath_wt)
rrp6_bam <- BamFile(bamPath_rrp6)
wt_IVP_bam <- BamFile(bamPath_wt_IVP)
slu7_IVP_bam <- BamFile(bamPath_slu7_IVP)

params_IMD4 <- ScanBamParam(which = IMD4_5exon, what = scanBamWhat())
params_snR54 <- ScanBamParam(which = snR54, what = scanBamWhat())
params <- ScanBamParam(which = RPL18A, what = scanBamWhat())

aln_wt <- scanBam(wt_bam, param = params)
aln_rrp6 <- scanBam(rrp6_bam, param = params)
filterBam(file = wt_IVP_bam, destination = "./wt_IVP_IMD4.bam", param = params_IMD4)
filterBam(file = slu7_IVP_bam, destination = "./slu7_IVP_IMD4.bam", param = params_IMD4)

wt_IVP_IMD4 <- BamFile("./wt_IVP_IMD4.bam")
slu7_IVP_IMD4 <- BamFile("./slu7_IVP_IMD4.bam")

wt_IVP_IMD4_snR54 <- scanBam(wt_IVP_IMD4, param = params_snR54)
slu7_IVP_IMD4_snR54 <- scanBam(slu7_IVP_IMD4, param = params_snR54)

blot_data_wt <-   data.frame("qname"=c(aln_wt$`chrXV:93395-94402`$qname,aln_rrp6$`chrXV:93395-94402`$qname),
                             "qwidth"=c(aln_wt$`chrXV:93395-94402`$qwidth,aln_rrp6$`chrXV:93395-94402`$qwidth))

blot_data_wt[1:length(aln_wt$`chrXV:93395-94402`$qname),"sample_ID"] <- "wt"
blot_data_wt[1:length(aln_wt$`chrXV:93395-94402`$qname),"sample_pos"] <- 1
blot_data_wt[((length(aln_wt$`chrXV:93395-94402`$qname)+1):(length(aln_rrp6$`chrXV:93395-94402`$qname)+length(aln_wt$`chrXV:93395-94402`$qname))),"sample_ID"] <- "rrp6"
blot_data_wt[((length(aln_wt$`chrXV:93395-94402`$qname)+1):(length(aln_rrp6$`chrXV:93395-94402`$qname)+length(aln_wt$`chrXV:93395-94402`$qname))),"sample_pos"] <- 2

blot_data_imd4 <- data.frame("qname"=c(wt_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname,slu7_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname),
                             "qwidth"=c(wt_IVP_IMD4_snR54$`chrXIII:163535-163620`$qwidth,slu7_IVP_IMD4_snR54$`chrXIII:163535-163620`$qwidth))
blot_data_imd4[1:length(wt_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname),"sample_ID"] <- "wt"
blot_data_imd4[1:length(wt_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname),"sample_pos"] <- 1
blot_data_imd4[((length(wt_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname)+1):(length(slu7_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname)+length(wt_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname))),"sample_ID"] <- "slu7"
blot_data_imd4[((length(wt_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname)+1):(length(slu7_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname)+length(wt_IVP_IMD4_snR54$`chrXIII:163535-163620`$qname))),"sample_pos"] <- 2

saveRDS(object = blot_data_wt, file = "./RPL18A_data")
saveRDS(object = blot_data_imd4, file = "./IMD4_data")

test <- rbind(blot_data_imd4, blot_data_imd4, blot_data_imd4, blot_data_imd4, blot_data_imd4, blot_data_imd4, blot_data_imd4, blot_data_imd4)
blot_data_imd4_spag <- rbind(test, test, test, test, test, test, test, test)

IMD4 <- ggplot(data = blot_data_imd4)+
  scale_y_continuous(limits = c(0,800))+
  geom_tile(aes(x = sample_pos, y = qwidth))+
  ylab(label = "Size in nts")+
  labs(title = "IMD4 5'Exon + snR54 - \"Nanoblot\"")

IMD4 <- ggplot(data = blot_data_imd4)+
  geom_point(aes(x = sample_pos, y = qwidth))+
  ylab(label = "Size in nts")+
  labs(title = "IMD4 5'Exon + snR54 - \"Nanoblot\"")

print(IMD4)

coding <-ggplot(data = blot_data_wt, aes(x = sample_pos, y = qwidth))+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  scale_fill_distiller(palette=2, direction=1) +
  scale_x_continuous(limits = c(0,3))+
  scale_y_continuous(limits = c(0,800))+
  ylab(label = "Size in nts")+
  labs(title = "RPL18A - \"Nanoblot\"")

ggplot(mtcars) +
  geom_point_blur(aes(mpg, wt, blur_size = disp), blur_steps = 6) +
  scale_blur_size_continuous(range = c(1, 15)) +
  theme_bw() +
  labs(title = "Larger blur indicates larger engine displacement")


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
