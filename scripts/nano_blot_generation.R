#Load libraries 

suppressPackageStartupMessages(library("BiocManager", quietly = TRUE))
suppressPackageStartupMessages(library("ggplot2", quietly = TRUE))
suppressPackageStartupMessages(library("Rsamtools", quietly = TRUE))
suppressPackageStartupMessages(library("ggridges", quietly = TRUE))

#The 3 args are currently 1) loading order 2) probe 3) duplication factor
args = commandArgs(trailingOnly=TRUE)
# This retrieves the commands that are passed to the R script, with the arguments stored inside
# a character vector

# This assumes that there is no antiprobe
if (length(args)==3) {
	sample_msg <- paste("Sample loading order:", args[1]) #ARGS 1 is the loading order, separated by commas
	probe_msg <- paste("Probe:", args[2]) # ARGS 2 is the target probe sequence
	dup_msg <- paste("Duplication Factor:", args[3]) # ARGS 3 is the duplication factor
	print("Starting plot generation.")
	print(sample_msg)
	print(probe_msg)
	print(dup_msg)
} else if (length(args)==4) { # This does assume that there is an antiprobe. This needs to be changed 
	# if we want to have multiple probes, since we can't be having a set number of arguments
  sample_msg <- paste("Sample loading order:", args[1])
  probe_msg <- paste("Probe:", args[2])
  dup_msg <- paste("Duplication Factor:", args[3])
  neg_probe_msg <- paste("Negative Probe:", args[4]) #ARGS 4 if it exists, is the antiprobe
  print("Starting plot generation.")
  print(sample_msg)
  print(probe_msg)
  print(dup_msg)
  print(neg_probe_msg)
}

bio_samples <- strsplit(x = args[1], split = ",") #Splits the sample loading order into a vector 
probe <- args[2] # Probe sequence, currently only accepts 1 inclusive probe and one exclusive probe
filenames <- as.vector("ERROR") #declaring the variable so it exists outside the loop scopes, so just assigning it any random value 

#Make locations of bamfiles

if (length(args)==3) {
	for (i in 1:length(bio_samples[[1]])) { #For loops through each line of bio samples, for each sample
		# REMINDER THAT R STARTS WITH INDEX 1!
		# These filenames already exist, its just finding those generated subsetted .bam files in the folders now
		filenames[i] <- paste(getwd(),
													"/temp/",
													bio_samples[[1]][i],
													"_",
													probe,
													".bam",
													sep = "")
	}
} else if (length(args)==4) {
	for (i in 1:length(bio_samples[[1]])) {
		filenames[i] <- paste(getwd(),
													"/temp/",
													bio_samples[[1]][i],
													"_",
													probe,
													"_anti_",
													args[4],
													".bam",
													sep = "")
	}
}


#Read in BAM files

bam_samples <- as.list("ERROR") #The bam_samples directory is created as an empty list
# Reminder, lists can contain any object type, which makes it easier compared to vectors
for (i in 1:length(filenames)) {
  bam_samples[[i]] <- scanBam(filenames[i])
}

#Extract relevant data from BAM files into dataframe

# To understand this --> need to look at documentation of scanBam
blot_data <- data.frame("qname"=c(bam_samples[[1]][[1]][["qname"]]),
                        "qwidth"=c(bam_samples[[1]][[1]][["qwidth"]]),
                        "sample_name"=bio_samples[[1]][[1]],
                        "row_number"=1)
                        
for (i in 2:length(bio_samples[[1]])) {
  print(i)
  temp_blot_data <- data.frame("qname"=c(bam_samples[[i]][[1]][["qname"]]),
                          "qwidth"=c(bam_samples[[i]][[1]][["qwidth"]]),
                          "sample_name"=bio_samples[[1]][[i]],
                          "row_number"=i)
  blot_data <- rbind(blot_data,temp_blot_data)
}

#Increase depth
for (i in 1:args[3]) {
  blot_data <- rbind(blot_data,blot_data)
}

#Make folder names + plot names 

pre_plot_name <- paste("_" ,c(strsplit(args[1], split = ",")[[1]]), sep = "", collapse = "")

folder_name <-
	paste("./plots/",
				args[2],
				"_",
				pre_plot_name,
				"_",
				args[3],
				"/",
				sep = "")

dir.create(folder_name)

plot_name_nano <-
	paste(folder_name,
				"nanoblot__",
				args[2],
				"_",
				pre_plot_name,
				"_",
				args[3],
				".png",
				sep = "")

plot_name_ridge <-
	paste(folder_name,
				"nanoridge__",
				args[2],
				"_",
				pre_plot_name,
				"_",
				args[3],
				".png",
				sep = "")

plot_name_violin <-
	paste(folder_name,
				"nanoviolin__",
				args[2],
				"_",
				pre_plot_name,
				"_",
				args[3],
				".png",
				sep = "")

#Add fuzz
blot_data$row_number_fuzz <- blot_data$row_number + runif(nrow(blot_data), min = -0.45, max = 0.45)

#Plotting graphs
plot_fuzzed <- ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
  geom_point(alpha = 0.05, color = "black", size = 2)+
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
