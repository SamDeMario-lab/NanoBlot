#' Generate Nanoblot from NanoblotData
#'
#' This function takes a NanoblotData object, a plotInfo dataframe, and
#' returns a ggplot object
#' @param nanoblotData The result of running bamFilesToNanoblotData
#' @param plotInfo A data frame.with 3 columns of equal lengths called SampleID, SampleLanes and, SampleColors
#' @param blotType The type of Nanoblot to produce. There are 3 acceptable options: 'blot', 'violin', and 'ridge'. Default 'blot'
#' @param plotTitle An optional title for the resulting plot. Default blank
#' @export
#' @examples
#'
#' ExampleDataset <- bamFilesToNanoblotData()
#'
#' ExamplePlotInfo <- data.frame(
#'   SampleID = c(WT,Test),
#'   SampleLanes = c(1,2),
#'   SampleColors = c('Red','Blue')
#'     )
#'
#' makeNanoblot(nanoblotData = ExampleDataset,
#' plotInfo = ExamplePlotInfo,
#' plotTitle = "Test plot",
#' )
#'

makeNanoblot <-
	function(nanoblotData,
					 plotInfo,
					 blotType = 'blot',
					 plotTitle = "",
					 size_factors = rep(1,length(levels(nanoblotData$SampleID)))) {
		## Start with the logical checks
		infoCol = c("SampleID",
								"SampleLanes",
								"SampleColors")
		if (!identical(names(plotInfo),
									 infoCol)) {
			stop("names(plotInfo) does not match expected values. Check formatting.")
		}
		## Duplicates data if blotType "blot" was selected
		if (blotType == 'blot') {
			duplication_factors <- calculateDuplicationFactors(nanoblotData, size_factors)
			nanoblotData <- duplicateNanoblotData(nanoblotData, duplication_factors)
		}

		# Add check for if all of the SampleIDs in infoCol exist in the nanoblotData
		if (sum(!(plotInfo$SampleID %in% as.vector(levels(nanoblotData$SampleID)))) != 0)    {
			stop("not all samples ID's in plotInfo param match with sample names in nanoblotData")
		}

		## Add sample lanes to nanoblotData
		nanoblotData <- merge(nanoblotData, plotInfo, by="SampleID")
		## Add fuzz to lanes
		columnWidth <- 0.25
		nanoblotData$SampleLanesFuzz <-
			nanoblotData$SampleLanes + runif(nrow(nanoblotData),
																			 min = -columnWidth,
																			 max = columnWidth)
		nanoblotData$SampleColors <- as.factor(nanoblotData$SampleColors)
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
			for (Sample in unique(plotInfo$SampleID)) {
				NanoPlot <- NanoPlot + ggplot2::geom_point(
					data = nanoblotData[nanoblotData$SampleID == Sample,],
					alpha = 0.01,
					color = nanoblotData[nanoblotData$SampleID == Sample,]$SampleColors[[1]],
					size = 1
				)
			}
		}

		if (blotType == 'violin') {
			##Check for multiplexing
			##Check for unnormalized data
			NanoPlot <- ggplot2::ggplot(data = nanoblotData, ggplot2::aes(x = SampleLanes, y = qwidth, group = SampleLanes, fill = SampleColors))+
				ggplot2::theme(axis.line = ggplot2::element_line(colour = "white"),
											 panel.grid.major = ggplot2::element_blank(),
											 panel.grid.minor = ggplot2::element_blank(),
											 panel.border = ggplot2::element_blank(),
											 panel.background = ggplot2::element_blank())+
				ggplot2::ylab(label = "Size in nts")+
				ggplot2::xlab(label = "")+
				ggplot2::geom_violin(
					data = nanoblotData,
					show.legend = FALSE
				)+
				ggplot2::scale_fill_manual(values=levels(nanoblotData$SampleColors))
		}
		if (blotType == 'ridge') {
			##Check for multiplexing
			##Check for unnormalized data
			NanoPlot <- ggplot2::ggplot(data = nanoblotData, ggplot2::aes(x = qwidth, y = SampleLanes, group = SampleLanes, fill = SampleColors))+
				ggplot2::theme(axis.line = ggplot2::element_line(colour = "white"),
											 panel.grid.major = ggplot2::element_blank(),
											 panel.grid.minor = ggplot2::element_blank(),
											 panel.border = ggplot2::element_blank(),
											 panel.background = ggplot2::element_blank())+
				ggplot2::ylab(label = "")+
				ggplot2::xlab(label = "Size in nts")+
				ggridges::geom_density_ridges2(
					data = nanoblotData,
					show.legend = FALSE
				)+
				ggplot2::scale_fill_manual(values=levels(nanoblotData$SampleColors))+
				ggplot2::ggtitle(label = plotTitle)
		}

		return(NanoPlot)
	}
