## acceptable plotType 'blot', 'violin', 'ridge', 'dataRDS'

isUnique <-
	function(vector) {
		return(!any(duplicated(vector)))
	}

extractNanoblotData <-
	function(SampleList) {
		NBdf <- data.frame("qname"=SampleList[[1]][['qname']],
							 "qwidth"=SampleList[[1]][['qwidth']])
		return(NBdf)
	}

scanBamFiles <-
	function(SampleID,
					 BamFileLocations) {
		print("Noice")
## check that sampleIDs and locations are unique
		if (!isUnique(SampleID)) {
			stop("SampleID contains non-unique names. All IDs must be unique.")
		}
		if (!isUnique(BamFileLocations)) {
			stop("BamFileLocations contains non-unique files. Check filenames.")
		}
		BamFileList <- lapply(BamFileLocations, Rsamtools::scanBam)
		names(BamFileList) <- SampleID
		return(BamFileList)
	}

bamFilesToNanoblotData <-
	function(BamFileList) {
		DataframesExtract <- lapply(BamFileList, extractNanoblotData)
		nanoblotData <- do.call("rbind", DataframesExtract)
		nanoblotData$'SampleID' <-
			as.factor(vapply(strsplit(row.names(nanoblotData), "\\."), `[`, 1, FUN.VALUE =
											 	character(1)))
		return(nanoblotData)
	}

makeNanoblot <-
	function(nanoblotData,
					 plotInfo,
					 blotType = 'blot',
					 plotTitle = "") {
## Start with the logical checks
		infoCol = c("SampleID",
								"SampleLanes", 
								"SampleColors")
		if (!identical(names(plotInfo),
									infoCol)) {
			stop("names(plotInfo) does not match expected values. Check formatting.")
		}
## Add check for if all of the SampleIDs exist in the nanoblotData		
## Add sample lanes to nanoblotData
		nanoblotData <- merge(nanoblotData, plotInfo, by="SampleID")
## Add fuzz to lanes
		columnWidth <- 0.25
		nanoblotData$SampleLanesFuzz <-
			nanoblotData$SampleLanes + runif(nrow(nanoblotData),
																			 min = -columnWidth,
																			 max = columnWidth)
## Create Nanoblot
		if (blotType == 'blot') {
			NanoPlot <- ggplot2::ggplot(data = nanoblotData, ggplot2::aes(x = SampleLanesFuzz, y = qwidth))+
				ggplot2::theme(axis.line = ggplot2::element_line(colour = "white"),
							panel.grid.major = ggplot2::element_blank(),
							panel.grid.minor = ggplot2::element_blank(),
							panel.border = ggplot2::element_blank(),
							panel.background = ggplot2::element_blank())+
				ggplot2::scale_x_continuous(breaks = c(1:max(plotInfo$"SampleLanes")),
			limits = c(
				min(plotInfo$"SampleLanes") - columnWidth,
				max(plotInfo$"SampleLanes") + columnWidth
			))+
				ggplot2::ylab(label = "Size in nts")+
				ggplot2::xlab(label = "")
			for (Sample in levels(nanoblotData$SampleID)) {
				NanoPlot <- NanoPlot + ggplot2::geom_point(
					data = nanoblotData[nanoblotData$SampleID == Sample,],
					alpha = 0.01,
					color = nanoblotData[nanoblotData$SampleID == Sample,]$SampleColors[[1]],
					size = 1
				)
			}
		}
		
		print("At least the function ran...")
		print(NanoPlot)
	}

TestDataNames <- c("WT_RPL18A", 
									 "RRP6_RPL18A")
TestDataFiles <- c("./temp/WT_RPL18A_Exon1.bam",
										"./temp/RRP6_RPL18A_Exon1.bam")

WT_RRP6Test <- data.frame(
	SampleID = c("WT_RPL18A","RRP6_RPL18A"),
	SampleLanes = c(1,2),
	SampleColors = c('red','blue')
)

BamsWTRRP6 <- scanBamFiles(TestDataNames, TestDataFiles)

NanoBlotDataWTRRP6 <- bamFilesToNanoblotData(BamsWTRRP6)

makeNanoblot(nanoblotData = NanoBlotDataWTRRP6, plotInfo = WT_RRP6Test)
