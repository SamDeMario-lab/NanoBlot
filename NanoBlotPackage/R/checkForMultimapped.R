#' Check for Multimapping
#' 
#' @param BamFileLocation The path to the BAM file to be
#' read in. The BAM files is imported to R via Rsamtools::scanBam.  
#' @param WarnPercent A number between 0 and 100. If more then the specified 
#' percent are reads in a BAM file are specified as secondary alignments then
#' a warning is issued.  
#' @export
#' @examples
#' checkForMultimapped("Path/to/bam/file.bam")
#' 
#' 
#' 
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