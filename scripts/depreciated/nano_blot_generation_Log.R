#Load libraries 

suppressPackageStartupMessages(library("BiocManager", quietly = TRUE))
suppressPackageStartupMessages(library("ggplot2", quietly = TRUE))
suppressPackageStartupMessages(library("Rsamtools", quietly = TRUE))
suppressPackageStartupMessages(library("ggridges", quietly = TRUE))

#The 3 args are currently 1) loading order 2) probe 3) duplication factor
args = commandArgs(trailingOnly=TRUE)
# This retrieves the commands that are passed to the R script, with the arguments stored inside
# a character vector

sample_msg <- paste("Sample loading order:", args[1])
probe_msg <- paste("Probe(s):", args[2])
dup_msg <- paste("Duplication Factor:", args[3])
neg_probe_msg <- paste("Negative Probe(s):", args[5])
print("Starting plot generation.")
print(sample_msg)
print(probe_msg)
print(dup_msg)
if (length(args)==4) {
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

#Increase depth, this is based off of the duplication factor
# There has to be some better way to do this compared to just duplicating the # of plot points
for (i in 1:args[3]) {
  blot_data <- rbind(blot_data,blot_data)
}

#Make folder names + plot names 

pre_plot_name <- paste("_" ,c(strsplit(args[1], split = ",")[[1]]), sep = "", collapse = "")

folder_name <-
	paste("./plots/",
				file_path,
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
				file_path,
				"_",
				pre_plot_name,
				"_",
				args[3],
				".png",
				sep = "")

plot_name_ridge <-
	paste(folder_name,
				"nanoridge__",
				file_path,
				"_",
				pre_plot_name,
				"_",
				args[3],
				".png",
				sep = "")

plot_name_violin <-
	paste(folder_name,
				"nanoviolin__",
				file_path,
				"_",
				pre_plot_name,
				"_",
				args[3],
				".png",
				sep = "")

blot_data_RDS <-
	paste(folder_name,
				"blot_data",
				file_path,
				"_",
				pre_plot_name,
				"_",
				args[3],
				sep = "")

#Add fuzz, this is adding it so that the duplication factor does not produce random points
# This is important to note because of the row_number here, this is essentially creating the width of 
# the lane, similar to what a northern blot would produce
# Default, -0.45, 0.45 
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
  xlab(label = "")+
	scale_y_continuous(trans = "log2")

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

saveRDS(object = blot_data, file = blot_data_RDS)
# ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
#   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
#   theme(axis.line = element_line(colour = "white"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   ylab(label = "Size in nts")+
#   xlab(label = "")
