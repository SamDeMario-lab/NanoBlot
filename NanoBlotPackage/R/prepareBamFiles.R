#' @name bamFileListToNanoblotData
#' @title Converting BamFileList object to Nanoblot dataframe suitable for Nanoblot analysis
#' @description Extracts relevant fields from BamFile within BamFileList and converts into a dataframe format
#'
#' Note: All .bam files are loaded into RAM. This can result in issues if the
#' bam files are large. Make sure that bam files are subset already through subsetNanoblot
#'
#' All paths and names of BamFileList object must be unique. Names of each element of BamFileList entry will be automatically
#' determined based on file path. User can override by using names(BamFileList) <- c()
#' @param BamFileList A BamFileList object from the Rsamtools class
#' @export
#' @examples
#' TestBamFileList <- Rsamtools::BamFileList(c("./temp/WT_sorted_merged_RPL18A_Exon1.bam", "./temp/RRP6_sorted_merged_RPL18A_Exon1.bam"))
#' names(TestBamFileList) <- c("WT","RRP6")
#' NanoblotDataWTRRP6 <- bamFileListToNanoblotData(TestBamFileList)
#'
#'
bamFileListToNanoblotData <-
  function(BamFileList) {
    ## check that sampleIDs and locations are unique
    BamFileListNames <- names(BiocGenerics::path(BamFileList))
    if (!isUnique(BiocGenerics::path(BamFileList))) {
      stop("BamFileList paths contain non-unique names. All file paths must be unique.")
    }

    if (!isUnique(BamFileListNames)) {
      stop("BamFileList names are non unique. All names must be unique.")
    }

    ListOfBams <- lapply(BamFileList, Rsamtools::scanBam, param = Rsamtools::ScanBamParam(what = c("qname","qwidth")))
    lapply(as.vector(BiocGenerics::path(BamFileList)), checkForMultimapped)
    names(ListOfBams) <- BamFileListNames
    DataframesExtract <- lapply(seq_along(ListOfBams), extractNanoblotData,
                                SampleNames = names(ListOfBams),
                                SampleList = ListOfBams)
    nanoblotData <- do.call("rbind", DataframesExtract)
    #nanoblotData$'SampleID' <-as.factor(vapply(strsplit(row.names(nanoblotData), "\\."), `[`, 1, FUN.VALUE = character(1)))
    return(nanoblotData)
  }



















