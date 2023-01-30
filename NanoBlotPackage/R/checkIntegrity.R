#' Calculate Integrity 
#' 
#' @param GeneTargets A GRanges object with the target genomic regions. 
#' These genomic regions should theoretically have no 5' transcripts ends. 
#' @param BamFiles A BamFileList() from Rsamtools 
#' @export
#' @examples
#' makeNanoblot()
#' 
#' 

calculateIntegrity <- function(GeneTargets, BamFiles) {
GRanFilter <- ScanBamParam(which = GeneTargets)
totalCounts <- countBam(file = BamFiles, param = GRanFilter)

GRanFilterEnds <- ScanBamParam(which = GeneTargets )

return(totalCounts)

	#Make filter for loading in bam files
#Load in the strand positon and qwidth
#get number of reads store in a variable
#Check strand of reads
	#Check if read has a 5' end
	#Count numbers of 5' ends
#Calculate a number of incomplete reads
#Calculate  of complete reads/incomplete reads
#return a dataframe with the number of incomplete reads, complete reads and ratio
}
