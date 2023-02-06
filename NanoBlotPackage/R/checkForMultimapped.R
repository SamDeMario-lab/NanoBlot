checkForMultimapped <-
	function(BamFileLocation, WarnPercent = 10) {
TotalReads<- Rsamtools::countBam(BamFileLocation)$records
SecondaryAlignments <- Rsamtools::countBam(
	BamFileLocation,
	param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSecondaryAlignment = TRUE))
)$records
SecondaryPercent <- (SecondaryAlignments/TotalReads)*100
if (SecondaryPercent>WarnPercent) {
	warning(paste("More then 10% of reads in ",BamFileLocation," reported as secondary alignment."))
}
}