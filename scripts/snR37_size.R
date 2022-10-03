if (!require("BiocManager", quietly = TRUE))
	install.packages("BiocManager")

suppressPackageStartupMessages(library("BiocManager", quietly = TRUE))
suppressPackageStartupMessages(library("ggplot2", quietly = TRUE))
suppressPackageStartupMessages(library("Rsamtools", quietly = TRUE))
suppressPackageStartupMessages(library("ggridges", quietly = TRUE))
suppressPackageStartupMessages(library("DESeq2", quietly = TRUE))
suppressPackageStartupMessages(library("dplyr", quietly = TRUE))

blot_data <- readRDS(file = "/home/guillaume-chanfreau/github/Nano_blot/plots/snR37/blot_data__snR37__slu7_term_")

size <- data.frame(qwidth = 386, row_number = 2, row_number_fuzz = 1.5 + runif(1000,max = 0.5))

plot_fuzzed <- ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
	geom_point(alpha = 0.01, color = "black", size = 1)+
	theme(axis.line = element_line(colour = "white"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank())+
	scale_x_continuous(limits = c(0,3)) +
	ylab(label = "Size in nts")+
	xlab(label = "")+
	geom_point(data = size ,alpha = 0.04, color = "red", size = 1)

print(plot_fuzzed)

ggsave(plot = plot_fuzzed, filename = "/home/guillaume-chanfreau/github/Nano_blot/plots/snR37/sizeestimate.png", device = "png", units = "in", height = 9, width = 4)
