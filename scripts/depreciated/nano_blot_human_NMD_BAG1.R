
suppressPackageStartupMessages(library("BiocManager", quietly = TRUE))
suppressPackageStartupMessages(library("ggplot2", quietly = TRUE))
suppressPackageStartupMessages(library("Rsamtools", quietly = TRUE))
suppressPackageStartupMessages(library("ggridges", quietly = TRUE))
suppressPackageStartupMessages(library("DESeq2", quietly = TRUE))
suppressPackageStartupMessages(library("dplyr", quietly = TRUE))

COLUMN_WIDTH <- 0.25

blot_data <- readRDS(file = "~/github/Nano_blot/plots/BAG1_RT/blot_data__BAG1_Exon__scr1_scr2_upf1_smg6_smg7_smg67_")

plot_fuzzed <- ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
	geom_point(alpha = 0.01, color = "black", size = 1)+
	theme(axis.line = element_line(colour = "white"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank())+
	scale_x_continuous()+
	scale_y_continuous(limits = c(60,280))+
	ylab(label = "Size in nts")+
	xlab(label = "")

print(plot_fuzzed)

ggsave(filename = "./plots/BAG1_RT/tweaked_BAG1.png", device = "png", units = "in", width = 6, height = 4)
