#' @name scanBamFiles
#' @title Reading in bam files with scanBam
#' @description A wrap around function which runs Rsamtools::scanBam on a series
#' of .bam files. 
#' 
#' Note: All .bam files are loaded into RAM. This can result in issues if the 
#' bam files are large. 
#' @param SampleID A vector of unique sample names.
#' @param BamFileLocations A vector containing the paths to the bam files to be
#' read in. .bam files are imported to R via Rsamtools::scanBam. 
#' @export
#' @examples
#' sampleNames <- c('WT', 'test')
#' sampleLocations <- c('/path/to/WT.bam', '/path/to/test.bam')
#' 
#' SampleBamFiles <- scanBamFiles(SampleID=sampleNames, BamFileLocations=sampleLocations)
#' 
#' SampleNanoblotData <- bamFilesToNanoblotData(SampleBamFiles)
#' 
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
		DataframesExtract <- lapply(seq_along(BamFileList), extractNanoblotData, 
																SampleNames = names(BamFileList),
																SampleList = BamFileList)
		nanoblotData <- do.call("rbind", DataframesExtract)
		#nanoblotData$'SampleID' <-as.factor(vapply(strsplit(row.names(nanoblotData), "\\."), `[`, 1, FUN.VALUE = character(1)))
		return(nanoblotData)
	}