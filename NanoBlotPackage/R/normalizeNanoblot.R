#' @name calculateDESeqSizeFactors
#' @title Calculate DESeq Size Factors using DESeq2
#' @description This function is required to normalize Nanoblots using the DeSeq normalization method.
#' It is not meant to be used by the end user.
#' It uses a GTF annotation file and
#' creates a count table using unnormalized files which is then passed to DeSeq2's sizeFactors() function
#' @param nanoblotData This is a data frame that is generated from bamFilesToNanoblotData
#' @param unnormalizedFiles This is a vector that contains the string file locations of all the samples in nanoblotData
#' @param annotationFile This is a gtf file obtained from a database that contains the organism of interest's genomic annotations
#' @param coldata This is a DeSeq2 coldata argument that is called with the DeSeq2 function DESeqDataSetFromMatrix
#' @export
#' @examples
#' TestDataNames <- c("WT_RPL18A", "RRP6_RPL18A")
#' TestDataFiles <- c("./temp/WT_RPL18A_Exon1.bam", "./temp/RRP6_RPL18A_Exon1.bam")
#' annotation <- "./user_input_files/Saccharomyces_cerevisiae.R64-1-1.107.gtf"
#' unnormalizedLocations <- c("./data/example/WT_sorted_merged.bam", "./data/example/RRP6_sorted_merged.bam")
#'
#' BamsWTRRP6 <- scanBamFiles(TestDataNames, TestDataFiles)
#' NanoblotDataWTRRP6 <- bamFilesToNanoblotData(BamsWTRRP6)
#' calculateDESeqSizeFactors(NanoblotDataWTRRP6, unnormalizedLocations, annotation)
#'
#' returns a vector with 2 elements
calculateDESeqSizeFactors <- function(nanoblotData, unnormalizedFiles, annotationFile, coldata) {
	sampleNames <- as.vector(levels(nanoblotData$SampleID))
	fc_SE <- Rsubread::featureCounts(files = unnormalizedFiles,
																	 annot.ext = annotationFile,
																	 isGTFAnnotationFile = TRUE,
																	 isLongRead = TRUE)
	if (is.null(coldata)) {
	  coldata <- data.frame(condition = rep("treated", length(levels(nanoblotData$SampleID))))
	  rownames(coldata) <- colnames(fc_SE$counts)
	  dds <- DESeq2::DESeqDataSetFromMatrix(countData = fc_SE$counts,
	                                        colData = coldata,
	                                        design = ~ 1)
	  #Remove low counts
	  keep <- rowSums(DESeq2::counts(dds)) >= 5
	  dds <- dds[keep,]
	}
	else {
	  dds <- DESeq2::DESeqDataSetFromMatrix(countData = fc_SE$counts,
	                                        colData = coldata,
	                                        design = ~ condition)
	  #Remove low counts
	  keep <- rowSums(DESeq2::counts(dds)) >= 5
	  dds <- dds[keep,]
	  dds <- DESeq2::DESeq(dds)
	  res <- DESeq2::results(dds)
	  res #Prints out the DESeq2 result
	  summary(res) #Prints out a summary of the DESeq2 result

	}

	#Run DESeq2, which is split into these different functions
	dds <- DESeq2::estimateSizeFactors(dds)
	size_factors <- DESeq2::sizeFactors(dds)
	names(size_factors) <- sampleNames
	cat("DESeq2 Size Factors\n-------\n")
	print(size_factors)
	return(size_factors)
}

#' @name calculateLibrarySizeFactors
#' @title Calculate Normalization Size Factors for NanoBlot without annotation file
#' @description This function is required to normalize Nanoblots using a counts per greatest common factor normalization method.
#' It is not meant to be used by the end user.
#' It does not use a GTF annotation file.
#' @param nanoblotData This is a data frame that is generated from bamFilesToNanoblotData
#' @param unnormalizedFiles This is a vector that contains the string file locations of all the samples in nanoblotData
#' @export
#' @examples
#' TestDataNames <- c("WT_RPL18A", "RRP6_RPL18A")
#' TestDataFiles <- c("./temp/WT_RPL18A_Exon1.bam", "./temp/RRP6_RPL18A_Exon1.bam")
#' unnormalizedLocations <- c("./data/example/WT_sorted_merged.bam", "./data/example/RRP6_sorted_merged.bam")
#'
#' BamsWTRRP6 <- scanBamFiles(TestDataNames, TestDataFiles)
#' NanoblotDataWTRRP6 <- bamFilesToNanoblotData(BamsWTRRP6)
#' calculateDESeqSizeFactors(NanoblotDataWTRRP6, unnormalizedLocations)
#'
#' returns a vector with 2 elements
calculateLibrarySizeFactors <- function(nanoblotData, unnormalizedFiles) {

	cat("=======\nNormalization: Counts Per Library Depth\n")
	sampleNames <- as.vector(levels(nanoblotData$SampleID))
	size_factors <- c()
	param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE))
	for (i in 1:length(levels(nanoblotData$SampleID))) {
		raw_read_number <- Rsamtools::countBam(file = unnormalizedFiles[i], param = param)[[6]]
		cat(paste(sampleNames[i], raw_read_number, "reads\n", sep = " "))
		size_factors[i] <- strtoi(raw_read_number)
	}
	size_factors <- max(size_factors) / size_factors
	cat("Counts Per Library Depth Size Factor\n-------\n")
	names(size_factors) <- sampleNames
	print(size_factors)
	return(size_factors)
}

#' @name calculateDuplicationFactors
#' @title Calculate Duplication Factors for NanoBlot based off of size factors for plot type "blot"
#' @description This function is required to plot Nanoblots of type "blot" by creating duplication factors to standardize
#' visualization opacity. It is not meant to be used by the end user.
#' This function is called in the makeNanoblot function when the plot type of "blot" is specified
#' @param nanoblotData This is a data frame that is generated from bamFilesToNanoblotData
#' @param size_factors This is a vector that contains the size factors for each sample in the nanoblotData
#' @param DUPLICATION_CONSTANT This is essentiall a value that determines how many "points" a standard Nanoblot should have.
#' The higher this number, the greater the opacity of the visualization. Default value is set to 2000
#' @export
#' @examples
#' TestDataNames <- c("WT_RPL18A", "RRP6_RPL18A")
#' TestDataFiles <- c("./temp/WT_RPL18A_Exon1.bam", "./temp/RRP6_RPL18A_Exon1.bam")
#' unnormalizedLocations <- c("./data/example/WT_sorted_merged.bam", "./data/example/RRP6_sorted_merged.bam")
#'
#' BamsWTRRP6 <- scanBamFiles(TestDataNames, TestDataFiles)
#' NanoblotDataWTRRP6 <- bamFilesToNanoblotData(BamsWTRRP6)
#' size_factors <- calculateDESeqSizeFactors(NanoblotDataWTRRP6, unnormalizedLocations)
#' calculateDuplicationFactors(NanoblotDataWTRRP6, size_factors)
#'
#' returns a vector with 2 elements of data type containing only whole numbers
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

#' @name duplicateNanoblotData
#' @title Duplicate NanoblotData for the purposes of "blot" type visualization
#' @description This function is required to plot Nanoblots of type "blot" by duplicating existing Nanoblot data for
#' visualization opacity. It is not meant to be used by the end user.
#' This function is called in the makeNanoblot function when the plot type of "blot" is specified
#' @param nanoblotData This is a data frame that is generated from bamFilesToNanoblotData
#' @param duplication_factors This is a vector that contains the duplication factors for each sample in the nanoblotData and generated
#' from the function calculateDuplicationFactors
#' @export
#' @examples
#' TestDataNames <- c("WT_RPL18A", "RRP6_RPL18A")
#' TestDataFiles <- c("./temp/WT_RPL18A_Exon1.bam", "./temp/RRP6_RPL18A_Exon1.bam")
#' unnormalizedLocations <- c("./data/example/WT_sorted_merged.bam", "./data/example/RRP6_sorted_merged.bam")
#'
#' BamsWTRRP6 <- scanBamFiles(TestDataNames, TestDataFiles)
#' NanoblotDataWTRRP6 <- bamFilesToNanoblotData(BamsWTRRP6)
#' size_factors <- calculateDESeqSizeFactors(NanoblotDataWTRRP6, unnormalizedLocations)
#' duplication_factors <- calculateDuplicationFactors(NanoblotDataWTRRP6, size_factors)
#' duplicateNanoblotData(NanoblotDataWTRRP6, duplication_factors)
#'
#' returns the original NanoblotDataWTRRP6 data frame with duplicated data for each sample
duplicateNanoblotData <- function(nanoblotData, duplication_factors) {
	sampleNames <- as.vector(levels(nanoblotData$SampleID))
	progressBar <- txtProgressBar(min = 0, max = sum(duplication_factors),
																style = 3, width = 50, char = "=")
	for (i in 1:length(duplication_factors)) {
		added_data <- dplyr::filter(nanoblotData, SampleID == sampleNames[i])
		if (duplication_factors[[i]] <= 1 | is.nan(duplication_factors[[i]])) {next;}
		for (j in 1:duplication_factors[[i]]) {
		  nanoblotData <- rbind(nanoblotData, added_data)
		  if (i == 1) {
				setTxtProgressBar(progressBar, j)
			} else {
				setTxtProgressBar(progressBar, sum(duplication_factors[1:(i-1)]) + j)
			}
		}
	}
	close(progressBar)
	return(nanoblotData)
}

#' @name normalizeNanoblotData
#' @title Normalizes Nanoblot data
#' @description This function is used to generate size factors for Nanoblot data normalization. The user can specify two
#' types of normalization, either 'differential' which uses DeSeq2 normlization or 'size' which uses a library size normalization
#' method for sequencing data that is not based off of previously annotated genomic regions
#' @param nanoblotData This is a data frame that is generated from bamFilesToNanoblotData
#' @param normalizationType This is the type of normalization, either 'differential' or 'size'. DEFAULT is 'differential'
#' @param unnormalizedFiles This is a vector that contains the string file locations of all the samples in nanoblotData.
#' @param annotationFile This is a gtf file obtained from a database that contains the organism of interest's genomic annotations.
#' This is only needed for the 'differential' normalization method
#' @export
#' @examples
#' TestDataNames <- c("WT_RPL18A", "RRP6_RPL18A")
#' TestDataFiles <- c("./temp/WT_RPL18A_Exon1.bam", "./temp/RRP6_RPL18A_Exon1.bam")
#' unnormalizedLocations <- c("./data/example/WT_sorted_merged.bam", "./data/example/RRP6_sorted_merged.bam")
#'
#' BamsWTRRP6 <- scanBamFiles(TestDataNames, TestDataFiles)
#' NanoblotDataWTRRP6 <- bamFilesToNanoblotData(BamsWTRRP6)
#' size_factors <- calculateDESeqSizeFactors(NanoblotDataWTRRP6, normalizationType = 'size', unnormalizedLocations)
#'
#' returns a vector with 2 elements of data type containing only whole numbers
normalizeNanoblotData <-
	function(nanoblotData,
					 normalizationType = "differential",
					 unnormalizedFiles,
					 annotationFile = NA,
					 coldata = NULL) {
		#Different types of checks
		# vector length of unnormalizedFiles has to be the same as number of samples
		# all files have to exist
		# annotationFile has to exist if type is differential
		# input of nanoblotData has to be correct

		if (normalizationType == "differential") { size_factors <- calculateDESeqSizeFactors(nanoblotData, unnormalizedFiles, annotationFile,coldata) }
		else if (normalizationType == "size") { size_factors <- calculateLibrarySizeFactors(nanoblotData, unnormalizedFiles)}
		else { stop("Normalization type can only be differential or size. Please check spelling and try again") }
		return(size_factors)
	}
