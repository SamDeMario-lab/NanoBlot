#kevin's attempt to make this function run faster

YeastCDSGFF <- rtracklayer::import.gff("~/Downloads/saccharomyces_cerevisiae.20210411.CDS.no_mito.gff")
TestBamFileList <- Rsamtools::BamFileList(c("./../data/example/WT_sorted_merged.bam",
                                            "./../data/example/RRP6_sorted_merged.bam"))

GeneTargets = YeastCDSGFF
BamFiles = TestBamFileList

counts1 <- calculateIntegrity(GeneTargets, BamFiles)
counts2 <- calculateIntegrity2(GeneTargets, BamFiles)

calculateIntegrity2 <- function(GeneTargets, BamFiles) {

  GRanFilter <- Rsamtools::ScanBamParam(which = GeneTargets)
  totalCounts <- c()
  for (SampleNum in seq_along(BamFiles)) {
    countsPerBam <- Rsamtools::countBam(file = BamFiles[[SampleNum]], param = GRanFilter)
    countsPerBam <- dplyr::select(countsPerBam, -records)
    countsPerBam <- dplyr::mutate(countsPerBam, records = c(0), FivePrimeEnds = c(0))

    currentBamFile <- BamFiles[[SampleNum]]
    Rsamtools::yieldSize(currentBamFile) <- 50000
    bf <- Rsamtools::open.BamFile(currentBamFile)
    while(length(chunk <- GenomicAlignments::readGAlignments(bf))) {
      currentRecords <- countsPerBam$records
      currentFivePrimeEnds <- countsPerBam$FivePrimeEnds
      newRecords <- GenomicAlignments::countOverlaps(GeneTargets, chunk)
      newFivePrimeEnds <- GenomicAlignments::countOverlaps(GeneTargets, chunk, type = "within")
      countsPerBam <- dplyr::mutate(countsPerBam,
                                    records = currentRecords + newRecords,
                                    FivePrimeEnds = currentFivePrimeEnds + newFivePrimeEnds)
    }
    Rsamtools::close.BamFile(bf)

    countsPerBam <- dplyr::mutate(countsPerBam, ReadEndsInRegion = (FivePrimeEnds/records) * 100)
    # Combine the countsPerBam into the totalCounts data frame
    totalCounts <- rbind(totalCounts, countsPerBam)
  }
  totalCounts$file <- as.factor(totalCounts$file)
  plotoutput <- ggplot2::ggplot(data = totalCounts)+
    ggplot2::geom_point(ggplot2::aes(x=width,y=ReadEndsInRegion, colour=file))
  print(plotoutput)

  return(totalCounts)
}
