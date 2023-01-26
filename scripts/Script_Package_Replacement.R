### Notes/todo
### For loop for creating filenames should probobly be changed. 

library(NanoBlotPackage)

args = commandArgs(trailingOnly=TRUE)
###
sample_names <- args[1]
probe_name <- args[2]
NORM_FACTOR <- args[3]
file_path <- args[4]
meta_data_file <- args[5]
plot_title <- args[6]
neg_probe_msg <- paste("Negative Probe(s):", args[7])
if (length(args)==6) {
	neg_probe_msg <- "No Negative Probes"
}
bio_samples <- strsplit(x = args[1], split = ",") #Splits the sample loading order into a vector 
filenames <- as.vector("ERROR") 

#Make locations of bamfiles

print(file_path)
for (i in 1:length(bio_samples[[1]])) {
	filenames[i] <- paste(getwd(), "/temp/", bio_samples[[1]][i], "_", file_path, ".bam", sep="")
	bam_samples[[i]] <- scanBam(filenames[i])
	}
###

subset_bam_files <- 
plot_data <- 


NanoBlotPackage::makeNanoblot(plotTitle = plot_title)
