library(ggplot2)

df_1 <- readRDS(file = "plots/snR4__WT_RRP6_1/blot_datasnR4__WT_RRP6_1")
df_2 <- readRDS(file = "plots/RPS7B__WT_RRP6_1/blot_dataRPS7B__WT_RRP6_1")

df_1$Probe <- "snR4"
df_2$Probe <- "RPS7B"

blot_data <- rbind(df_1, df_2)

plot_fuzzed <- ggplot(data = blot_data, aes(x = row_number_fuzz, y = qwidth, colour = Probe))+
	geom_point(alpha = .15, size = 2, shape = 16)+
	theme(axis.line = element_line(colour = "white"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank())+
	ylab(label = "Size in nts")+
	xlab(label = "")+
	ylim(c(200,1000))

print(plot_fuzzed)
