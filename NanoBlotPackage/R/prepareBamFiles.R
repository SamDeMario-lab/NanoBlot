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