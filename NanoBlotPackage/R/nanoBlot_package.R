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
		nanoblotData$'SampleID' <-as.factor(vapply(strsplit(row.names(nanoblotData), "\\."), `[`, 1, FUN.VALUE = character(1)))
		#tried adding mutate to the extractNanoBlotData command instead using lapply, seq_along, and names to no luck, 
		#can try again if need be to fix the SampleID 
		return(nanoblotData)
	}

# Auxillary function to calculate DESeq size factors
calculateDESeqSizeFactors <- function(nanoblotData, unnormalizedFiles, annotationFile) {
	sampleNames <- as.vector(levels(nanoblotData$SampleID))
	fc_SE <- Rsubread::featureCounts(files = unnormalizedFiles,
																	 annot.ext = annotationFile, 
																	 isGTFAnnotationFile = TRUE,
																	 isLongRead = TRUE)
	coldata <- data.frame(condition = rep("treated", length(levels(nanoblotData$SampleID))))
	rownames(coldata) <- colnames(fc_SE$counts)
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = fc_SE$counts,
																				colData = coldata,
																				design = ~ 1)
	#Remove low counts
	keep <- rowSums(DESeq2::counts(dds)) >= 5
	dds <- dds[keep,]
	#Run DESeq2, which is split into these different functions
	dds <- DESeq2::estimateSizeFactors(dds)
	size_factors <- DESeq2::sizeFactors(dds)
	names(size_factors) <- sampleNames
	cat("DESeq2 Size Factors\n-------\n")
	print(size_factors) 
	return(size_factors)
}

# Auxillary function to calculate CPM size factors
calculateCPMSizeFactors <- function(nanoblotData, unnormalizedFiles) {
	
	cat("=======\nNormalization: Counts Per Million\n")
	sampleNames <- as.vector(levels(nanoblotData$SampleID))
	size_factors <- c()
	param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE))
	for (i in 1:length(levels(nanoblotData$SampleID))) {
		raw_read_number <- Rsamtools::countBam(file = unnormalizedFiles[i], param = param)[[6]]
		cat(paste(sampleNames[i], raw_read_number, "reads\n", sep = " "))
		size_factors[i] <- strtoi(raw_read_number) / 1000000
	}
	cat("Counts Per Million Size Factor\n-------\n")
	names(size_factors) <- sampleNames
	print(size_factors) 
	return(size_factors)
}

# Auxillary function to calculate duplication factors from size factors
calculateDuplicationFactors <- function(nanoblotData, size_factors, DUPLICATION_CONSTANT = 2000) {
	
	sampleNames <- as.vector(levels(nanoblotData$SampleID))
	duplication_factors <- c()
	raw_reads <- c()
	for (i in 1:length(size_factors)) {
		raw_read_number <- nrow(dplyr::filter(nanoblotData, SampleID == sampleNames[i]))
		# This calculation is taking the "normalization_factor" and finding the inverse, then multiplying by
		# 10 to have meaningful effect, and then rounding to the nearest digit
		duplication_factors[i] <- round((1/size_factors[[i]]) * 10, digits = 0)
		if (raw_read_number == 0){
			duplication_factors[i] <- NaN #This allows it so that there are no reads mapped to this sample
		}
		raw_reads[i] <- raw_read_number
	}
	
	max_raw_read <- max(raw_reads)
	cat("Max raw read count:", max_raw_read, "\n")
	duplication_factors <- round(duplication_factors * (DUPLICATION_CONSTANT / max_raw_read), digits = 0)
	names(duplication_factors) <- sampleNames
	cat("Duplication Factors\n-------\n")
	print(duplication_factors) 
}

#Auxillary function that duplicates Nanoblot data frame based on vector of duplication factors 
duplicateNanoblotData <- function(nanoblotData, duplication_factors) {
	sampleNames <- as.vector(levels(nanoblotData$SampleID))
	for (i in 1:length(duplication_factors)) {
		added_data <- dplyr::filter(nanoblotData, SampleID == sampleNames[i])
		if (duplication_factors[[i]] <= 1 | is.nan(duplication_factors[[i]])) {next;}
		for (j in 1:duplication_factors[[i]]) {
			nanoblotData <- rbind(nanoblotData, added_data)
		}
	}
	return(nanoblotData)
}

normalizeNanoblotData <-
	function(nanoblotData,
					 normalizationType = "differential", 
					 unnormalizedFiles, 
					 annotationFile = NA) {
		#Different types of checks
		# vector length of unnormalizedFiles has to be the same as number of samples
		# all files have to exist
		# annotationFile has to exist if type is differential
		# input of nanoblotData has to be correct
		
		if (normalizationType == "differential") { size_factors <- calculateDESeqSizeFactors(nanoblotData, unnormalizedFiles, annotationFile) }
		else if (normalizationType == "size") { size_factors <- calculateCPMSizeFactors(nanoblotData, unnormalizedFiles)}
		else { stop("Normalization type can only be differential or size. Please check spelling and try again") }
		
		duplication_factors <- calculateDuplicationFactors(nanoblotData, size_factors)
		return(duplicateNanoblotData(nanoblotData, duplication_factors))
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
		print(NanoPlot)
	}

TestDataNames <- c("WT_RPL18A", 
									 "RRP6_RPL18A")
TestDataFiles <- c("./temp/WT_RPL18A_Exon1.bam",
										"./temp/RRP6_RPL18A_Exon1.bam")
annotation <- "./user_input_files/Saccharomyces_cerevisiae.R64-1-1.107.gtf"
unnormalizedLocations <- c("../TRAMP_Direct/WT/sorted_merged.bam", "../TRAMP_Direct/rrp6/sorted_merged.bam")

WT_RRP6Test <- data.frame(
	SampleID = c("WT_RPL18A","RRP6_RPL18A"),
	SampleLanes = c(1,2),
	SampleColors = c('red','blue')
)

BamsWTRRP6 <- scanBamFiles(TestDataNames, TestDataFiles)
NanoBlotDataWTRRP6 <- bamFilesToNanoblotData(BamsWTRRP6)
NanoBlotDataWTRRP6 <- normalizeNanoblotData(NanoBlotDataWTRRP6, "differential", unnormalizedLocations, annotation)
makeNanoblot(nanoblotData = NanoBlotDataWTRRP6, plotInfo = WT_RRP6Test)