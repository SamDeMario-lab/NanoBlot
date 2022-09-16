#Load libraries 
if (!require("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

suppressPackageStartupMessages(library("BiocManager", quietly = TRUE))
suppressPackageStartupMessages(library("ggplot2", quietly = TRUE))
suppressPackageStartupMessages(library("Rsamtools", quietly = TRUE))
suppressPackageStartupMessages(library("ggridges", quietly = TRUE))
suppressPackageStartupMessages(library("DESeq2", quietly = TRUE))
suppressPackageStartupMessages(library("dplyr", quietly = TRUE))

#The 3 args are currently 1) loading order 2) probe 3) duplication factor
args = commandArgs(trailingOnly=TRUE)
# This retrieves the commands that are passed to the R script, with the arguments stored inside
# a character vector

sample_msg <- paste("Sample loading order:", args[1])
probe_msg <- paste("Probe(s):", args[2])
NORM_FACTOR <- args[3]
plot_title <- args[6]
neg_probe_msg <- paste("Negative Probe(s):", args[7])
print("Starting plot generation.")
print(sample_msg)
print(probe_msg)
if (length(args)==6) {
	neg_probe_msg <- "No Negative Probes"
}
print(neg_probe_msg)

bio_samples <- strsplit(x = args[1], split = ",") #Splits the sample loading order into a vector 
filenames <- as.vector("ERROR") 

#Make locations of bamfiles
file_path <- args[4]
print(file_path)
for (i in 1:length(bio_samples[[1]])) {
	filenames[i] <- paste(getwd(), "/temp/", bio_samples[[1]][i], "_", file_path, ".bam", sep="")
}

#Read in BAM files
bam_samples <- as.list("ERROR") #The bam_samples directory is created as an empty list
for (i in 1:length(filenames)) {
  bam_samples[[i]] <- scanBam(filenames[i])
}

size_factors <- c() #Creating variable name
if (NORM_FACTOR == 0) {
	# Differential normalization type
	count_filenames <- c()
	for (i in 1:length(bio_samples[[1]])) {
		count_filenames <- append(count_filenames, paste(bio_samples[[1]][i], "-htseq_counts.tsv", sep=""))
	}
	table <- data.frame(sampleName = count_filenames, fileName = count_filenames)
	path=paste(getwd(), "/temp/count_tables", sep="")
	#Create Deseq2 data set
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table, directory = path, design= ~ 1)
	#Remove low counts
	keep <- rowSums(counts(ddsHTSeq)) >= 5
	ddsHTSeq <- ddsHTSeq[keep,]
	#Run DESeq2, which is split into these different functions
	dds <- estimateSizeFactors(ddsHTSeq)
	size_factors <- sizeFactors(dds)
	cat("DESeq2 Size Factors\n-------\n")
	print(size_factors) 
	
} else if (NORM_FACTOR == 1) {
	# Size normalization type
	meta_data_file <- args[5]
	meta_data_table <- read.table(meta_data_file, sep='\t', header = TRUE)
	
	#setting param so it does not count any unmapped reads or secondary read alignments
	param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE))
	for (i in 1:length(bio_samples[[1]])) {
		data_location <- (meta_data_table %>% filter(Sample_name == bio_samples[[1]][i]))[["Location"]]
		#Gets the number of reads in the raw bam file, then converts to CPM 
		raw_read_number <- countBam(file = c(data_location), param = param)[[6]]
		size_factors[i] <- strtoi(raw_read_number) / 1000000
	}
	cat("Counts Per Million Size Factor\n-------\n")
	names(size_factors) <- bio_samples[[1]]
	print(size_factors) 
	
} else if (NORM_FACTOR == 2) {
	# Skipping normalization 
	size_factors <- rep(c(1), times= length(bio_samples[[1]]))
	cat("Skipping Normalization Factors\n-------\n")
}

duplication_factors <- c()
raw_reads <- c()
for (i in 1:length(size_factors)) {
	param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE))
	raw_read_number <- countBam(file = c(filenames[i]), param = param)[[6]]
	# This calculation is taking the "normalization_factor" and finding the inverse, then multiplying by
	# 10 to have meaningful effect, and then rounding to the nearest digit
	duplication_factors[i] <- round((1/size_factors[[i]]) * 10, digits = 0)
	if (raw_read_number == 0){
		duplication_factors[i] <- NaN #This allows it so that there are no reads mapped to this sample
	}
	raw_reads[i] <- raw_read_number
}

DUPLICATION_CONSTANT = 2000
max_raw_read <- max(raw_reads)
cat("Max raw read count:", max_raw_read, "\n")
duplication_factors <- round(duplication_factors * (DUPLICATION_CONSTANT / max_raw_read), digits = 0)
names(duplication_factors) <- bio_samples[[1]]
cat("Duplication Factors\n-------\n")
print(duplication_factors) 

#Extract relevant data from BAM files into dataframe
blot_data <- data.frame("qname"=character(),
												"qwidth"=integer(),
												"sample_name"=character(),
												"row_number"=integer())

for (i in 1:length(bio_samples[[1]])) {
	blot_data <- blot_data %>% add_row("qname"=c(bam_samples[[i]][[1]][["qname"]]),
																		 "qwidth"=c(bam_samples[[i]][[1]][["qwidth"]]),
																		 "sample_name"=bio_samples[[1]][[i]],
																		 "row_number"=i)
}

for (i in 1:length(duplication_factors)) {
	added_data <- blot_data %>% filter(row_number == i)
	if (duplication_factors[[i]] <= 1 | is.nan(duplication_factors[[i]])) {next;}
	for (j in 1:duplication_factors[[i]]) {
		blot_data <- rbind(blot_data, added_data)
	}
}

#Make folder names + plot names 

pre_plot_name <- paste("_" ,c(strsplit(args[1], split = ",")[[1]]), sep = "", collapse = "")

# Makes the plots folder if it doesn't exist already
output_dir <- file.path(getwd(), "plots/")
if (!dir.exists(output_dir)){
	dir.create(output_dir)
} 
	
folder_name <-
	paste("./plots/",
				plot_title,
				"/",
				sep = "")

dir.create(folder_name)

plot_name_nano <-
	paste(folder_name,
				"nanoblot__",
				file_path,
				"_",
				pre_plot_name,
				".png",
				sep = "")

plot_name_ridge <-
	paste(folder_name,
				"nanoridge__",
				file_path,
				"_",
				pre_plot_name,
				".png",
				sep = "")

plot_name_violin <-
	paste(folder_name,
				"nanoviolin__",
				file_path,
				"_",
				pre_plot_name,
				".png",
				sep = "")
				
blot_data_RDS <-
	paste(folder_name,
				"blot_data__",
				file_path,
				"_",
				pre_plot_name,
				"_",
				sep = "")


# Default, -0.45, 0.45 
COLUMN_WIDTH <- 0.25
blot_data$row_number_fuzz <- blot_data$row_number + runif(nrow(blot_data), min = -COLUMN_WIDTH, max = COLUMN_WIDTH)

#Plotting graphs
plot_fuzzed <- ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
  geom_point(alpha = 0.01, color = "black", size = 1)+
  theme(axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
	scale_x_continuous(breaks = c(1:length(bio_samples[[1]])), 
									 labels = c(bio_samples[[1]]),
									 limits = c(1-COLUMN_WIDTH, length(bio_samples[[1]]) + COLUMN_WIDTH)) +
  ylab(label = "Size in nts")+
  xlab(label = "")

plot_ridge <- ggplot(data = blot_data)+
	geom_density_ridges2(aes(x = qwidth, y = row_number, group = row_number, fill = row_number), show.legend = FALSE)+
	theme(axis.line = element_line(colour = "white"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.ticks.y = element_blank())+
	xlab(label = "Size in nts")+
	ylab(label = "") +
	scale_y_continuous(breaks = c(1:length(bio_samples[[1]])), 
										 labels = c(bio_samples[[1]]))

plot_violin <- ggplot(data = blot_data)+
	geom_violin(aes(x = row_number, y = qwidth, group = row_number, fill = row_number), show.legend = FALSE)+
	theme(axis.line = element_line(colour = "white"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.ticks.x = element_blank())+
	ylab(label = "Size in nts")+
	xlab(label = "") + 
	scale_x_continuous(breaks = c(1:length(bio_samples[[1]])), 
										 labels = c(bio_samples[[1]]))

saveRDS(object = blot_data, file = blot_data_RDS)

ggsave(filename = plot_name_nano ,plot = plot_fuzzed)
ggsave(filename = plot_name_ridge ,plot = plot_ridge)
ggsave(filename = plot_name_violin ,plot = plot_violin)
