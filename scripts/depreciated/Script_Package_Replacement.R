###To Do
#Normalization
#plot title
#Make all types of plots
#Save plots to a reasonable place. 

### Notes/todo
### For loop for creating filenames should probobly be changed. 
# $BAMS $PROBE_FIELD $NORM_FACTOR ${PREVIOUS_ANTI_PROBE} $META_DATA $PLOT_TITLE $ANTIPROBE_FIELD
# [1] "WT,RRP6"
# [2] "ACT1"
# [3] "0"
# [4] "ACT1"
# [5] "./user_input_files/data_metadata.tsv"
# [6] "ACT1"

library(NanoBlotPackage)

args = commandArgs(trailingOnly=TRUE)
test_args <-
	c("WT,RRP6",
		"ACT1",
		"0",
		"ACT1",
		"./user_input_files/data_metadata.tsv",
		"ACT1")
args <- test_args
###
bio_samples <- strsplit(x = args[1], split = ",") #Splits the sample loading order into a vector
# Arg 2 is no longer being used because args[4] is sufficient 
probe_name <- args[2]
# Arg 3 is no longer being used because normalization occurs in R.
# Arg 3 needs to be replaced with a normalization style
NORM_FACTOR <- args[3]
SampleProbes <- paste("_", args[4], sep="") 
meta_data_file <- args[5]
plot_title <- args[6]
neg_probe_msg <- paste("Negative Probe(s):", args[7])
if (length(args)==6) {
	NegProbe <- FALSE
	neg_probe_msg <- "No Negative Probes"
} else if (length(args)==7) {
	SampleProbes <- paste(SampleProbes, args[7], sep="_")
}

 

FileNames <- as.list("ERROR")



#Make locations of bamfiles
<<<<<<< HEAD
print(args)

# print(file_path)
# for (i in 1:length(bio_samples[[1]])) {
# 	filenames[i] <- paste(getwd(), "/temp/", bio_samples[[1]][i], "_", file_path, ".bam", sep="")
# 	bam_samples[[i]] <- scanBam(filenames[i])
# 	}
# ###
# 
# subset_bam_files <- 
# plot_data <- 
# 
# 
# NanoBlotPackage::makeNanoblot(plotTitle = plot_title)
=======
for (i in seq_along(bio_samples[[1]])) {
	FileNames[i] <- paste(getwd(), "/temp/", bio_samples[[1]][i], SampleProbes, ".bam", sep="")
}

SampleNames <- sapply(bio_samples[[1]], paste, SampleProbes, sep = "")

PlotInformation <- data.frame(
	SampleID = SampleNames,
	SampleLanes = seq_along(bio_samples[[1]]),
	SampleColors = rep_len(x = 'black', length.out = length(bio_samples[[1]]))
)
###

NanoBlotData <- scanBamFiles(SampleID = SampleNames, BamFileLocations = FileNames)

makeNanoblot(nanoblotData = NanoBlotData,
						plotInfo = PlotInformation,
						blotType = 'blot',
						plotTitle = )
#subset_bam_files <- 
#plot_data <- 


#NanoBlotPackage::makeNanoblot(plotTitle = plot_title)
>>>>>>> 565d4522a7330f92c366453751617f3576e45eab
