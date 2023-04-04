library(ggplot2)

blot_data <- readRDS(file = "C:/Users/samde/Downloads/blot_data__snR4__WT_RRP6_")

plot_fuzzed <- ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth))+
	geom_point(alpha = .006, size = 2, shape = 16, colour = "black")+
	theme(axis.line = element_line(colour = "white"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank())+
	ylab(label = "Size in nts")+
	xlab(label = "")+
	ylim(c(200,1000))+
	xlim(c(0,8))

print(plot_fuzzed)

ggsave(filename = "./plots/snR4_size.png", plot = plot_fuzzed, device = "png",units = "in",width = 4, height = 8)
