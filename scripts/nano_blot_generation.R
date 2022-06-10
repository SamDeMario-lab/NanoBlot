library(ggplot2)
library(Rsamtools)

#The 3 args are currently 1) loading order 2) probe 3) duplication factor
#args = commandArgs(trailingOnly=TRUE)
args <- "WT,RRP6,SLU7,RRP6SLU7"
args[2] <- "RPL14A_Exon1"
args[3] <- 1

if (length(args)==0) {
  stop("This script requires 2 inputs. None detected.", call.=FALSE)
} else if (length(args)==1) {
  stop("This script requires 2 inputs. Only one detected.", call.=FALSE)
} else if (length(args)==2) {
  sample_msg <- paste("Sample loading order:", args[1])
  probe_msg <- paste("Probe:", args[2])
  print("Starting plot generation.")
  print(sample_msg)
  print(probe_msg)
}

bio_samples <- strsplit(x = args[1], split = ",")
probe <- args[2]
filenames <- as.vector("ERROR")

#Make locations of bamfiles

for (i in 1:length(bio_samples[[1]])) {
  filenames[i] <- paste(getwd(),
        "/temp/",
        bio_samples[[1]][i],
        "_",
        probe,
        ".bam",
        sep = "")
}

#Read in BAM files

bam_samples <- as.list("ERROR")
for (i in 1:length(filenames)) {
  bam_samples[[i]] <- scanBam(filenames[i])
}

#Extract relevant data from BAM files into dataframe

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

#Add fuzz
blot_data$row_number_fuzz <- blot_data$row_number + runif(nrow(blot_data), min = -0.45, max = 0.45)

#Plot graph
plot_fuzzed <- ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
  geom_point(alpha = 0.05, color = "black", size = 2)+
  theme(axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ylab(label = "Size in nts")+
  xlab(label = "")

pre_plot_name <- paste("_" ,c(strsplit(args[1], split = ",")[[1]]), sep = "", collapse = "")

plot_name <- paste("./plots/nanoblot__", args[2] ,"_",pre_plot_name, ".png", sep = "")

ggsave(filename = plot_name ,plot = plot_fuzzed)
# ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
#   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
#   theme(axis.line = element_line(colour = "white"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+
#   ylab(label = "Size in nts")+
#   xlab(label = "")
