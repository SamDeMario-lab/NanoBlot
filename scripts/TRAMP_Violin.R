library(ggplot2)
library(RColorBrewer)

blot_data <- readRDS(file = "plots/snR4__WT_RRP6_AIR1RRP6_AIR2RRP6_AIR1AIR2RRP6_1/blot_datasnR4__WT_RRP6_AIR1RRP6_AIR2RRP6_AIR1AIR2RRP6_1")

blot_d_2 <- data.frame("row_number" = 6 , "qwidth" = 187)

# blot_data[5755:5845,"qwidth"] <- rnorm(90 , mean=187, sd=1)
# blot_data[5846:5855,"qwidth"] <- runif(11 , min = 150, max = 700)
# blot_data[5755:5855,"row_number"] <- 6
# blot_data[5755:5855,"sample_name"] <- "Control"

plot_violin <- ggplot(data = blot_data)+
	geom_violin(aes(x = row_number, y = qwidth, group = row_number, fill = sample_name), show.legend = FALSE)+
	theme(axis.line = element_line(colour = "white"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.ticks.x = element_blank())+
	ylab(label = "Size in nts")+
	xlab(label = "")+
	scale_y_continuous(limits = c(150,700))+
	scale_fill_brewer(palette="Accent")+
	geom_point(data = blot_d_2, aes(x = row_number, y = qwidth), shape = "-", size = 50, color = 'red4')

print(plot_violin)

ggsave(filename = "plots/snR4__WT_RRP6_AIR1RRP6_AIR2RRP6_AIR1AIR2RRP6_1/custem_violin_16AUG2022.png" ,plot = plot_violin)
