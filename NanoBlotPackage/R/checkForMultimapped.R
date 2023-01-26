checkForMultimapped <-
	function(BamFileLocation) {
TotalReads<- Rsamtools::countBam(BamFileLocation)$records
SecondaryAlignments <- Rsamtools::countBam(
	BamFileLocation,
	param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSecondaryAlignment = TRUE))
)$records
SecondaryPercent <- (SecondaryAlignments/TotalReads)*100
if (SecondaryPercent>10) {
	warning(paste("More then 10% of reads in ",BamFileLocation," reported as secondary alignment."))
}
}