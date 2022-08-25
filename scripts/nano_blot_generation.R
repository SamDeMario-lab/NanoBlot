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
neg_probe_msg <- paste("Negative Probe(s):", args[6])
print("Starting plot generation.")
print(sample_msg)
print(probe_msg)
if (length(args)==5) {
	neg_probe_msg <- "No Negative Probes"
}
print(neg_probe_msg)

bio_samples <- strsplit(x = args[1], split = ",") #Splits the sample loading order into a vector 
filenames <- as.vector("ERROR") #declaring the variable so it exists outside the loop scopes, so just assigning it any random value 

#Make locations of bamfiles
file_path <- args[4]
print(file_path)
for (i in 1:length(bio_samples[[1]])) {
	filenames[i] <- paste(getwd(), "/temp/", bio_samples[[1]][i], "_", file_path, ".bam", sep="")
}

#Read in BAM files
bam_samples <- as.list("ERROR") #The bam_samples directory is created as an empty list
# Reminder, lists can contain any object type, which makes it easier compared to vectors
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
for (i in 1:length(size_factors)) {
	# This calculation is taking the "normalization_factor" and finding the inverse, then multiplying by
	# 10 to have meaningful effect, and then rounding to the nearest digit
	duplication_factors[i] <- round((1/size_factors[[i]]) * 10, digits = 0)
}
names(duplication_factors) <- bio_samples[[1]]
cat("Duplication Factors\n-------\n")
print(duplication_factors) 

# So the variable bam_samples is a list of scanBam, which is in itself
# a list, so bam_samples is a list of lists 

#Extract relevant data from BAM files into dataframe

# To understand this --> need to look at documentation of scanBam
# qname is the name of the read
# qwidth is the length of the read
# sample_name is the name of this sample, taken from an earlier character vector
# row_number is just 1
# This code below basically creates the first column of sample_name with all the 
# qnames and qwidths corresponding to that sample_name 
# Walking through the indices, first set gets the first sample in bam_samples, second
# set gets the first set of alignments, here we only get one set
# the third set of brackets get's the specific list element using the name of the element
blot_data <- data.frame("qname"=c(bam_samples[[1]][[1]][["qname"]]),
                        "qwidth"=c(bam_samples[[1]][[1]][["qwidth"]]),
                        "sample_name"=bio_samples[[1]][[1]],
                        "row_number"=1)

# Now this is putting all the rest of the sample names into data frames
for (i in 2:length(bio_samples[[1]])) {
  #print(i) #Don't know why we need this print statement?
  temp_blot_data <- data.frame("qname"=c(bam_samples[[i]][[1]][["qname"]]),
                          "qwidth"=c(bam_samples[[i]][[1]][["qwidth"]]),
                          "sample_name"=bio_samples[[1]][[i]],
                          "row_number"=i)
  blot_data <- rbind(blot_data,temp_blot_data)
  # This joins the other data frame together, creating a master data frame called blot_data
}

for (i in 1:length(duplication_factors)) {
	added_data <- blot_data %>% filter(row_number == i)
	for (j in 1:duplication_factors[[i]]) {
		if (duplication_factors[[i]] <= 1) {break;}
		blot_data <- rbind(blot_data, added_data)
	}
}

#Make folder names + plot names 

pre_plot_name <- paste("_" ,c(strsplit(args[1], split = ",")[[1]]), sep = "", collapse = "")

folder_name <-
	paste("./plots/",
				file_path,
				"_",
				pre_plot_name,
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

#Add fuzz, this is adding it so that the duplication factor does not produce random points
# This is important to note because of the row_number here, this is essentially creating the width of 
# the lane, similar to what a northern blot would produce
# Default, -0.45, 0.45 
blot_data$row_number_fuzz <- blot_data$row_number + runif(nrow(blot_data), min = -0.25, max = 0.25)

#Plotting graphs
plot_fuzzed <- ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
  geom_point(alpha = 0.01, color = "black", size = 1)+
  theme(axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
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
	ylab(label = "")

plot_violin <- ggplot(data = blot_data)+
	geom_violin(aes(x = row_number, y = qwidth, group = row_number, fill = row_number), show.legend = FALSE)+
	theme(axis.line = element_line(colour = "white"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.ticks.x = element_blank())+
	ylab(label = "Size in nts")+
	xlab(label = "")

ggsave(filename = plot_name_nano ,plot = plot_fuzzed)
ggsave(filename = plot_name_ridge ,plot = plot_ridge)
ggsave(filename = plot_name_violin ,plot = plot_violin)
# ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
#   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
#   theme(axis.line = element_line(colour = "white"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   ylab(label = "Size in nts")+
#   xlab(label = "")
