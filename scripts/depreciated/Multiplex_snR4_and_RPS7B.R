library(ggplot2)

df_1 <- readRDS(file = "C:/Users/samde/Downloads/blot_data__RPS7B__WT_RRP6_")
df_2 <- readRDS(file = "C:/Users/samde/Downloads/blot_data__snR4__WT_RRP6_")

df_2$Probe <- "snR4"

df_1$Probe <- "RPS7B"

for (i in 1:nrow(df_1)) {
	if (df_1[i,"row_number"]==2) {
		df_1[i,"row_number_fuzz"] <- df_1[i,"row_number_fuzz"] + 1 
	}
}

for (i in 1:nrow(df_2)) {
	if (df_2[i,"row_number"]==1) {
		df_2[i,"row_number_fuzz"] <- df_2[i,"row_number_fuzz"] + 1 
	}
	if (df_2[i,"row_number"]==2) {
		df_2[i,"row_number_fuzz"] <- df_2[i,"row_number_fuzz"] + 2 
	}
}

blot_data <- rbind(df_1, df_2)

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
