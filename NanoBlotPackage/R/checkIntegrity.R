#' @name checkIntegrity
#' @title Check Integrity of Bam files for Nanoblot
#' @description Checks the percent of intact reads for all gene targets and then outputs a scatter plot of integrity
#' along with the cumulative distribution plot.
#'
#' @param GeneTargets A GRanges object with the target genomic regions.
#' These genomic regions should theoretically have no 5' transcripts ends.
#' @param BamFiles A BamFileList() from Rsamtools
#' @export
#' @examples
#' makeNanoblot()
#'
#' YeastCDSGFF <- rtracklayer::import.gff("~/Downloads/saccharomyces_cerevisiae.20210411.CDS.no_mito.gff")
#' TestBamFileList <- Rsamtools::BamFileList(c("./../data/example/WT_sorted_merged.bam", "./../data/example/RRP6_sorted_merged.bam"))
#'
#' counts <- calculateIntegrity(GeneTargets = YeastCDSGFF, BamFiles = TestBamFileList)

checkIntegrity <- function(GeneTargets, BamFiles) {

  GRanFilter <- Rsamtools::ScanBamParam(which = GeneTargets)
  totalBamCounts <- Rsamtools::countBam(file = BamFiles)
  totalRecords <- sum(totalBamCounts$records)
  uniqueBamNames <- unique(totalBamCounts$file)
  totalCounts <- c()
  progressBar <- txtProgressBar(min = 0, max = totalRecords,
  															style = 3, width = 50, char = "=")
  currentRecordCount <- 0
  setTxtProgressBar(progressBar, currentRecordCount)
  for (SampleNum in seq_along(BamFiles)) {
    countsPerBam <- Rsamtools::countBam(file = BamFiles, param = GRanFilter)
    countsPerBam <- dplyr::select(countsPerBam, -records)
    countsPerBam <- dplyr::mutate(countsPerBam, records = c(0), FivePrimeEnds = c(0))

    currentBamFile <- BamFiles[[SampleNum]]
    Rsamtools::yieldSize(currentBamFile) <- 100000
    bf <- Rsamtools::open.BamFile(currentBamFile)
    while(length(chunk <- GenomicAlignments::readGAlignments(bf))) {
      currentRecords <- countsPerBam$records
      currentFivePrimeEnds <- countsPerBam$FivePrimeEnds
      newRecords <- GenomicRanges::countOverlaps(GeneTargets, chunk)
      newFivePrimeEnds <- GenomicRanges::countOverlaps(GeneTargets, chunk, type = "within")
      countsPerBam <- dplyr::mutate(countsPerBam,
                                    records = currentRecords + newRecords,
                                    FivePrimeEnds = currentFivePrimeEnds + newFivePrimeEnds)
      currentRecordCount <- currentRecordCount + length(chunk)
      if (currentRecordCount > totalRecords) {currentRecordCount = totalRecords}
      setTxtProgressBar(progressBar, currentRecordCount)
    }
    Rsamtools::close.BamFile(bf)

    countsPerBam <- dplyr::mutate(countsPerBam, IntactReads = (FivePrimeEnds/records) * 100)
    # Combine the countsPerBam into the totalCounts data frame
    totalCounts <- rbind(totalCounts, countsPerBam)
  }
  totalCounts$file <- as.factor(totalCounts$file)

  scatterOutput <- ggplot2::ggplot(data = totalCounts)+
    ggplot2::geom_point(ggplot2::aes(x=width,y=IntactReads, colour=file)) +
    ggplot2::ylab(label = "% Integrity")
  print(scatterOutput)

  cdfOutput <-ggplot2::ggplot(data = totalCounts)+
    ggplot2::stat_ecdf(ggplot2::aes(x = IntactReads, color = file))+
    ggplot2::ylab(label = "Probability")+
    ggplot2::xlab(label = "% Partial Sequence")+
    ggplot2::theme_bw()+
    ggplot2::ggtitle(label = "Integrity comparison of RNA Samples")
  print(cdfOutput)

  return(totalCounts)
}
